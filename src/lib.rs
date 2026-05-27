//! Core logic for rsomics-vcf-fixref.
//!
//! Ports the FASTA-based check and fix operations from bcftools +fixref (MIT,
//! Genome Research Ltd). The dbSNP-ID modes (`id`, `stats`) and Illumina TOP
//! conversion (`top`) require a dbSNP VCF or Illumina manifest and are out of
//! scope; they are documented but not implemented.
//!
//! # Algorithm (fixref.c, MIT)
//!
//! For each biallelic SNP record the reference base `ir` is looked up by a
//! random-access fetch from a FASTA+FAI. The VCF REF and ALT single bases are
//! encoded as `A=0 C=1 G=2 T=3`; reverse complement is `revint(x) = 3 - x`
//! (A↔T, C↔G). Four cases:
//!
//! | ref FASTA | action |
//! |-----------|--------|
//! | ir == REF | no change (`FIX_NONE`) |
//! | ir == ALT | swap REF/ALT + swap GT 0↔1 (`FIX_SWAP`) |
//! | ir == rc(REF) | flip both alleles, keep GT (`FIX_FLIP`) |
//! | ir == rc(ALT) | flip + swap REF/ALT + swap GT (`FIX_FLIP_SWAP`) |
//! | none | error (`FIX_ERR`) |
//!
//! In `flip` mode (not `flip-all`) ambiguous pairs — A/T (`pair=0x9`) and C/G
//! (`pair=0x6`) — are skipped as unresolvable without a dbSNP anchor.
//!
//! Non-SNP records (REF or ALT longer than 1 bp, or multi-allelic) and
//! non-ACGT alleles are tallied but passed through unmodified.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::items_after_statements,
    clippy::too_many_lines,
    clippy::too_many_arguments
)]

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;

use rsomics_common::{Result, RsomicsError};

// ── Nucleotide encoding ───────────────────────────────────────────────────────

/// Sentinel for non-ACGT bases.
const NON_ACGT: u8 = u8::MAX;

/// A=0 C=1 G=2 T=3; `NON_ACGT` for anything else.
#[inline]
fn nt2int(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => NON_ACGT,
    }
}

/// Integer → nucleotide byte: only valid for 0..=3.
#[inline]
fn int2nt(x: u8) -> u8 {
    b"ACGT"[x as usize]
}

/// Reverse complement in integer space: A↔T (0↔3), C↔G (1↔2).
#[inline]
fn revint(x: u8) -> u8 {
    3 - x
}

// ── FASTA index ───────────────────────────────────────────────────────────────

/// One sequence entry in a .fai file.
struct FaiEntry {
    length: u64,
    /// Byte offset of the first base in the FASTA file.
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// Load a .fai index into a name → entry map (keys as `Vec<u8>` for zero-copy
/// lookup against the VCF CHROM field without UTF-8 conversion).
fn load_fai(fai_path: &Path) -> Result<HashMap<Vec<u8>, FaiEntry>> {
    let reader = BufReader::new(
        File::open(fai_path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fai_path.display())))?,
    );
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.splitn(5, '\t').collect();
        if cols.len() < 5 {
            return Err(RsomicsError::InvalidInput(format!(
                "malformed .fai line: {line}"
            )));
        }
        let parse_u64_str = |s: &str| {
            s.parse::<u64>()
                .map_err(|e| RsomicsError::InvalidInput(format!("bad .fai field '{s}': {e}")))
        };
        map.insert(
            cols[0].as_bytes().to_vec(),
            FaiEntry {
                length: parse_u64_str(cols[1])?,
                offset: parse_u64_str(cols[2])?,
                line_bases: parse_u64_str(cols[3])?,
                line_width: parse_u64_str(cols[4])?,
            },
        );
    }
    Ok(map)
}

/// Load an entire contig into a contiguous byte buffer, stripping newlines.
///
/// Bases are stored as raw FASTA ASCII (uppercase); callers index directly with
/// the 0-based position. Loading once per contig eliminates per-record seek
/// syscalls — the dominant cost for chromosome-sorted VCF inputs.
fn load_contig(fasta: &mut File, entry: &FaiEntry) -> Option<Arc<Vec<u8>>> {
    let newline_bytes = entry.line_width - entry.line_bases;
    let full_lines = entry.length / entry.line_bases;
    let remainder = entry.length % entry.line_bases;
    let raw_len = full_lines * entry.line_width
        + if remainder > 0 {
            remainder + newline_bytes
        } else {
            0
        };

    fasta.seek(SeekFrom::Start(entry.offset)).ok()?;
    let mut raw = vec![0u8; raw_len as usize];
    fasta.read_exact(&mut raw).ok()?;

    let mut bases = Vec::with_capacity(entry.length as usize);
    for b in raw {
        if b != b'\n' && b != b'\r' {
            bases.push(b);
        }
    }
    bases.truncate(entry.length as usize);
    Some(Arc::new(bases))
}

/// Look up a single base (0-based position) from a cached contig buffer.
#[inline]
fn base_from_cache(contig: &[u8], pos0: u64) -> Option<u8> {
    let b = *contig.get(pos0 as usize)?;
    let v = nt2int(b);
    if v == NON_ACGT { None } else { Some(v) }
}

// ── Fix mode ─────────────────────────────────────────────────────────────────

/// The operation to apply for this record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixMode {
    /// Check-only: pass through, emit stats to stderr.
    Check,
    /// Swap/flip based on single-base REF/ALT matching; skip A/T and C/G pairs.
    Flip,
    /// Same as flip but also processes ambiguous A/T and C/G pairs.
    FlipAll,
}

impl FixMode {
    #[must_use]
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "check" | "stats" => Some(Self::Check),
            "flip" => Some(Self::Flip),
            "flip-all" => Some(Self::FlipAll),
            _ => None,
        }
    }
}

// ── Per-record counters ────────────────────────────────────────────────────────

#[derive(Default, Debug)]
struct Stats {
    nsite: u64,
    nok: u64,
    nflip: u64,
    nswap: u64,
    nflip_swap: u64,
    nunresolved: u64,
    nerr: u64,
    nskip: u64,
    non_acgt: u64,
    non_snp: u64,
    non_biallelic: u64,
    /// Substitution type count indexed by `[ref_int][alt_int]`.
    count: [[u64; 4]; 4],
}

// ── Tab-split helpers ─────────────────────────────────────────────────────────

/// Find the byte offset of the `n`-th tab (0-indexed) in `line`, or `None`.
#[inline]
fn nth_tab(line: &[u8], n: usize) -> Option<usize> {
    let mut count = 0usize;
    for (i, &b) in line.iter().enumerate() {
        if b == b'\t' {
            if count == n {
                return Some(i);
            }
            count += 1;
        }
    }
    None
}

/// Field `n` (0-indexed) of a tab-split line, as a byte slice.
#[inline]
fn field(line: &[u8], n: usize) -> &[u8] {
    let start = if n == 0 {
        0
    } else {
        match nth_tab(line, n - 1) {
            Some(i) => i + 1,
            None => return b"",
        }
    };
    let end = match nth_tab(line, n) {
        Some(i) => i,
        None => {
            // Last field: strip trailing newline.
            let e = line.len();
            if e > 0 && (line[e - 1] == b'\n' || line[e - 1] == b'\r') {
                let e = if e >= 2 && line[e - 2] == b'\r' {
                    e - 2
                } else {
                    e - 1
                };
                return &line[start..e];
            }
            return &line[start..e];
        }
    };
    &line[start..end]
}

/// Write GT-swapped FORMAT+samples (from the 9th column onward) to `w`.
///
/// Swaps allele 0↔1 in every sample's GT field. Operates on byte slices to
/// avoid allocations.
fn write_swapped_gt_column(fmt_samples: &[u8], w: &mut dyn Write) -> std::io::Result<()> {
    // Find the GT field index within the FORMAT column.
    let tab_pos = fmt_samples.iter().position(|&b| b == b'\t');
    let format_col = match tab_pos {
        Some(p) => &fmt_samples[..p],
        None => fmt_samples,
    };

    let gt_idx = format_col.split(|&b| b == b':').position(|f| f == b"GT");

    let Some(gt_idx) = gt_idx else {
        // No GT field — pass through unchanged.
        return w.write_all(fmt_samples);
    };

    let samples_start = match tab_pos {
        Some(p) => p + 1,
        None => return w.write_all(fmt_samples),
    };

    w.write_all(format_col)?;

    let samples_bytes = &fmt_samples[samples_start..];
    // Iterate samples (tab-separated).
    for (si, sample) in samples_bytes.split(|&b| b == b'\t').enumerate() {
        if si == 0 {
            w.write_all(b"\t")?;
        } else {
            w.write_all(b"\t")?;
        }
        // Iterate FORMAT fields (colon-separated).
        for (fi, fld) in sample.split(|&b| b == b':').enumerate() {
            if fi > 0 {
                w.write_all(b":")?;
            }
            if fi == gt_idx {
                write_swapped_gt(fld, w)?;
            } else {
                w.write_all(fld)?;
            }
        }
    }
    Ok(())
}

/// Write a GT token with alleles 0 and 1 swapped, preserving phase separators.
#[inline]
fn write_swapped_gt(gt: &[u8], w: &mut dyn Write) -> std::io::Result<()> {
    let sep = if let Some(p) = gt.iter().position(|&b| b == b'/') {
        Some((p, b'/'))
    } else if let Some(p) = gt.iter().position(|&b| b == b'|') {
        Some((p, b'|'))
    } else {
        None
    };
    if let Some((pos, sep_byte)) = sep {
        write_swapped_allele(&gt[..pos], w)?;
        w.write_all(&[sep_byte])?;
        write_swapped_allele(&gt[pos + 1..], w)?;
    } else {
        write_swapped_allele(gt, w)?;
    }
    Ok(())
}

#[inline]
fn write_swapped_allele(tok: &[u8], w: &mut dyn Write) -> std::io::Result<()> {
    match tok {
        b"0" => w.write_all(b"1"),
        b"1" => w.write_all(b"0"),
        other => w.write_all(other),
    }
}

// ── Core per-line processor (zero-allocation hot path) ────────────────────────

/// Process one VCF data line (as raw bytes, with newline stripped by caller).
///
/// Writes the (possibly modified) line to `w` followed by `\n`, updating
/// `stats`. Returns the action string for FIXREF annotation (or `None` in
/// check mode / error / malformed input).
fn process_line_bytes(
    line: &[u8],
    contig_cache: Option<&[u8]>,
    mode: FixMode,
    stats: &mut Stats,
    skipped_chroms: &mut HashMap<Vec<u8>, bool>,
    annotate: bool,
    w: &mut dyn Write,
) -> std::io::Result<()> {
    stats.nsite += 1;

    // Quick field count check: need at least 8 tab-separated fields.
    let tab3 = nth_tab(line, 2); // after col2 (ID)
    if tab3.is_none() {
        stats.nerr += 1;
        w.write_all(line)?;
        w.write_all(b"\n")?;
        return Ok(());
    }

    let ref_col = field(line, 3);
    let alt_col = field(line, 4);
    let info_col_start = nth_tab(line, 6).map(|i| i + 1);

    if info_col_start.is_none() {
        stats.nerr += 1;
        w.write_all(line)?;
        w.write_all(b"\n")?;
        return Ok(());
    }

    // Multi-allelic: pass through.
    if memchr(b',', alt_col) {
        stats.non_biallelic += 1;
        stats.nskip += 1;
        return write_line_annotated(line, b"skip", annotate, w);
    }

    // Non-SNP: pass through.
    if ref_col.len() != 1 || alt_col.len() != 1 {
        stats.non_snp += 1;
        stats.nskip += 1;
        return write_line_annotated(line, b"skip", annotate, w);
    }

    let ia = nt2int(ref_col[0]);
    let ib = nt2int(alt_col[0]);

    if ia == NON_ACGT || ib == NON_ACGT {
        stats.non_acgt += 1;
        stats.nskip += 1;
        return write_line_annotated(line, b"skip", annotate, w);
    }

    // Parse POS (field 1).
    let pos_bytes = field(line, 1);
    let pos1: u64 = match parse_u64(pos_bytes) {
        Some(v) => v,
        None => {
            stats.nerr += 1;
            w.write_all(line)?;
            w.write_all(b"\n")?;
            return Ok(());
        }
    };
    let pos0 = pos1.saturating_sub(1);

    // Look up reference base.
    let Some(ir) = contig_cache.and_then(|seq| base_from_cache(seq, pos0)) else {
        let chrom = field(line, 0);
        if !skipped_chroms.contains_key(chrom) {
            let chrom_str = String::from_utf8_lossy(chrom);
            eprintln!("Ignoring sequence \"{chrom_str}\"");
            skipped_chroms.insert(chrom.to_vec(), true);
        }
        stats.nskip += 1;
        return write_line_annotated(line, b"skip", annotate, w);
    };

    stats.count[ia as usize][ib as usize] += 1;

    // In flip (not flip-all) mode skip ambiguous A/T and C/G pairs.
    if mode == FixMode::Flip {
        let pair = (1u8 << ia) | (1u8 << ib);
        if pair == 0x9 || pair == 0x6 {
            stats.nunresolved += 1;
            return write_line_annotated(line, b"skip", annotate, w);
        }
    }

    if mode == FixMode::Check {
        if ir == ia {
            stats.nok += 1;
        } else if ir == ib {
            stats.nswap += 1;
        } else if ir == revint(ia) {
            stats.nflip += 1;
        } else if ir == revint(ib) {
            stats.nflip_swap += 1;
        } else {
            stats.nerr += 1;
        }
        w.write_all(line)?;
        w.write_all(b"\n")?;
        return Ok(());
    }

    if ir == ia {
        // No change needed.
        stats.nok += 1;
        return write_line_annotated(line, b"none", annotate, w);
    }

    if ir == ib {
        // FIX_SWAP: swap REF and ALT, swap GT 0↔1.
        stats.nswap += 1;
        return write_modified_line(line, alt_col[0], ref_col[0], true, b"swap,GT", annotate, w);
    }

    if ir == revint(ia) {
        // FIX_FLIP: complement both alleles, keep GT.
        stats.nflip += 1;
        let new_ref = int2nt(revint(ia));
        let new_alt = int2nt(revint(ib));
        return write_modified_line(line, new_ref, new_alt, false, b"flip", annotate, w);
    }

    if ir == revint(ib) {
        // FIX_FLIP|FIX_SWAP|FIX_GT: flip+swap REF/ALT, swap GT.
        stats.nflip_swap += 1;
        let new_ref = int2nt(revint(ib));
        let new_alt = int2nt(revint(ia));
        return write_modified_line(line, new_ref, new_alt, true, b"flip,swap,GT", annotate, w);
    }

    stats.nerr += 1;
    write_line_annotated(line, b"err", annotate, w)
}

/// Write a VCF line with modified REF (col3), ALT (col4), and optionally
/// swapped GT, followed by FIXREF annotation if requested.
fn write_modified_line(
    line: &[u8],
    new_ref: u8,
    new_alt: u8,
    swap_gt: bool,
    action: &[u8],
    annotate: bool,
    w: &mut dyn Write,
) -> std::io::Result<()> {
    // Write fields 0..2 unchanged.
    let t2 = nth_tab(line, 2).unwrap(); // col3 start
    w.write_all(&line[..t2 + 1])?;
    // Write new REF.
    w.write_all(&[new_ref, b'\t'])?;
    // Write new ALT.
    w.write_all(&[new_alt, b'\t'])?;
    // Write QUAL, FILTER (fields 5..6) unchanged.
    let t4 = nth_tab(line, 4).unwrap();
    let t6 = nth_tab(line, 6).unwrap();
    w.write_all(&line[t4 + 1..t6 + 1])?;
    // INFO field.
    if annotate {
        let info_start = t6 + 1;
        let info_end = nth_tab(line, 7).unwrap_or(line.len());
        let info = &line[info_start..info_end];
        if info == b"." {
            write!(w, "FIXREF={}", String::from_utf8_lossy(action))?;
        } else {
            w.write_all(info)?;
            write!(w, ";FIXREF={}", String::from_utf8_lossy(action))?;
        }
    } else {
        let info_start = t6 + 1;
        let info_end = nth_tab(line, 7).unwrap_or(line.len());
        w.write_all(&line[info_start..info_end])?;
    }
    // FORMAT + samples (field 8+), possibly with swapped GT.
    if let Some(t7) = nth_tab(line, 7) {
        w.write_all(b"\t")?;
        let fmt_samples = &line[t7 + 1..];
        // Strip trailing newline from fmt_samples for processing.
        let fmt_samples = trim_newline(fmt_samples);
        if swap_gt {
            write_swapped_gt_column(fmt_samples, w)?;
        } else {
            w.write_all(fmt_samples)?;
        }
    }
    w.write_all(b"\n")
}

/// Write a VCF line pass-through (no REF/ALT changes).
///
/// Used for records we do not modify (skip, err, check mode, or the `none`
/// action in flip mode where the REF already matches). FIXREF annotation is
/// not emitted for these pass-through records; bcftools +fixref similarly does
/// not tag them.
fn write_line_annotated(
    line: &[u8],
    _action: &[u8],
    _annotate: bool,
    w: &mut dyn Write,
) -> std::io::Result<()> {
    w.write_all(line)?;
    w.write_all(b"\n")
}

// ── Small helpers ─────────────────────────────────────────────────────────────

#[inline]
fn memchr(needle: u8, haystack: &[u8]) -> bool {
    haystack.contains(&needle)
}

#[inline]
fn parse_u64(s: &[u8]) -> Option<u64> {
    if s.is_empty() {
        return None;
    }
    let mut n: u64 = 0;
    for &b in s {
        if b < b'0' || b > b'9' {
            return None;
        }
        n = n.checked_mul(10)?.checked_add((b - b'0') as u64)?;
    }
    Some(n)
}

#[inline]
fn trim_newline(s: &[u8]) -> &[u8] {
    let mut e = s.len();
    if e > 0 && s[e - 1] == b'\n' {
        e -= 1;
    }
    if e > 0 && s[e - 1] == b'\r' {
        e -= 1;
    }
    &s[..e]
}

// ── Stats output ──────────────────────────────────────────────────────────────

fn print_stats(stats: &Stats) {
    let ncmp = stats.nsite - stats.nskip;
    let nactive = ncmp;
    let tot: u64 = stats.count.iter().flatten().sum();

    eprintln!("# ST, substitution types");
    for i in 0u8..4 {
        for j in 0u8..4 {
            if i == j {
                continue;
            }
            let c = stats.count[i as usize][j as usize];
            let pct = if tot > 0 {
                c as f64 * 100.0 / tot as f64
            } else {
                0.0
            };
            eprintln!(
                "ST\t{}>{}	{c}\t{pct:.1}%",
                char::from(int2nt(i)),
                char::from(int2nt(j)),
            );
        }
    }

    let pct = |n: u64, d: u64| {
        if d > 0 {
            n as f64 * 100.0 / d as f64
        } else {
            0.0
        }
    };

    eprintln!("# NS, Number of sites:");
    eprintln!("NS\ttotal        \t{}", stats.nsite);
    eprintln!(
        "NS\tref match    \t{}\t{:.1}%",
        stats.nok,
        pct(stats.nok, ncmp)
    );
    eprintln!(
        "NS\tref mismatch \t{}\t{:.1}%",
        ncmp.saturating_sub(stats.nok),
        pct(ncmp.saturating_sub(stats.nok), ncmp)
    );
    eprintln!(
        "NS\tflipped      \t{}\t{:.1}%",
        stats.nflip,
        pct(stats.nflip, nactive)
    );
    eprintln!(
        "NS\tswapped      \t{}\t{:.1}%",
        stats.nswap,
        pct(stats.nswap, nactive)
    );
    eprintln!(
        "NS\tflip+swap    \t{}\t{:.1}%",
        stats.nflip_swap,
        pct(stats.nflip_swap, nactive)
    );
    eprintln!(
        "NS\tunresolved   \t{}\t{:.1}%",
        stats.nunresolved,
        pct(stats.nunresolved, nactive)
    );
    eprintln!("NS\terrors       \t{}", stats.nerr);
    eprintln!("NS\tskipped      \t{}", stats.nskip);
    eprintln!("NS\tnon-ACGT     \t{}", stats.non_acgt);
    eprintln!("NS\tnon-SNP      \t{}", stats.non_snp);
    eprintln!("NS\tnon-biallelic\t{}", stats.non_biallelic);
}

// ── Public API ────────────────────────────────────────────────────────────────

pub struct FixrefStats {
    pub total: u64,
    pub fixed: u64,
    pub nok: u64,
    pub nswap: u64,
    pub nflip: u64,
    pub nflip_swap: u64,
    pub nunresolved: u64,
    pub nerr: u64,
}

/// Stream a VCF through the fixref logic.
///
/// `fasta_path` must have an adjacent `.fai` index (`<fasta>.fai` or `<base>.fa.fai`).
/// Writes the (possibly modified) VCF to `output` and emits stats to stderr.
pub fn fixref(
    vcf_path: &Path,
    fasta_path: &Path,
    output: &mut dyn Write,
    mode: FixMode,
) -> Result<FixrefStats> {
    let fai_path = {
        let ext = fasta_path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("");
        let with_fai = fasta_path.with_extension(format!("{ext}.fai"));
        if with_fai.exists() {
            with_fai
        } else {
            let mut p = fasta_path.to_path_buf();
            p.set_extension("fai");
            p
        }
    };

    if !fai_path.exists() {
        return Err(RsomicsError::InvalidInput(format!(
            "FASTA index not found: {} — run `samtools faidx {}` first",
            fai_path.display(),
            fasta_path.display()
        )));
    }

    let index: HashMap<Vec<u8>, FaiEntry> = load_fai(&fai_path)?;
    let mut fasta = File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?;
    let vcf_file = File::open(vcf_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", vcf_path.display())))?;

    let mut reader = BufReader::new(vcf_file);
    let mut writer = BufWriter::new(output);
    let mut stats = Stats::default();
    let mut skipped_chroms: HashMap<Vec<u8>, bool> = HashMap::new();
    let mut contig_cache: HashMap<Vec<u8>, Option<Arc<Vec<u8>>>> = HashMap::new();

    let annotate = mode != FixMode::Check;
    let fixref_info_header = b"##INFO=<ID=FIXREF,Number=.,Type=String,Description=\"The change made by bcftools/fixref\">\n";

    let mut buf = Vec::with_capacity(4096);
    loop {
        buf.clear();
        let n = reader
            .read_until(b'\n', &mut buf)
            .map_err(RsomicsError::Io)?;
        if n == 0 {
            break;
        }
        let line = trim_newline(&buf);

        if line.starts_with(b"#") {
            if annotate && line.starts_with(b"#CHROM") {
                writer
                    .write_all(fixref_info_header)
                    .map_err(RsomicsError::Io)?;
            }
            writer.write_all(line).map_err(RsomicsError::Io)?;
            writer.write_all(b"\n").map_err(RsomicsError::Io)?;
            continue;
        }
        if line.is_empty() {
            continue;
        }

        // Extract CHROM and resolve its cached contig.
        let chrom = field(line, 0);
        let contig_seq = contig_cache
            .entry(chrom.to_vec())
            .or_insert_with(|| index.get(chrom).and_then(|e| load_contig(&mut fasta, e)));

        process_line_bytes(
            line,
            contig_seq.as_deref().map(Vec::as_slice),
            mode,
            &mut stats,
            &mut skipped_chroms,
            annotate,
            &mut writer,
        )
        .map_err(RsomicsError::Io)?;
    }

    writer.flush().map_err(RsomicsError::Io)?;
    print_stats(&stats);

    let fixed = stats.nswap + stats.nflip + stats.nflip_swap;
    Ok(FixrefStats {
        total: stats.nsite,
        fixed,
        nok: stats.nok,
        nswap: stats.nswap,
        nflip: stats.nflip,
        nflip_swap: stats.nflip_swap,
        nunresolved: stats.nunresolved,
        nerr: stats.nerr,
    })
}
