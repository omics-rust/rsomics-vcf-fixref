//! rsomics-vcf-fixref — port of bcftools +fixref (MIT, Genome Research Ltd).
//!
//! Checks and fixes VCF REF alleles against a reference FASTA using the
//! four-case fixref algorithm (none / swap / flip / flip+swap). The dbSNP-ID
//! modes (`id`, `stats`) and Illumina TOP conversion (`top`) are out of scope.
//!
//! # Algorithm
//!
//! For each biallelic SNP the reference base `ir` is fetched from the FASTA+FAI.
//! Alleles are encoded A=0 C=1 G=2 T=3; reverse complement is `3 − x`.
//!
//! | FASTA base | action |
//! |------------|--------|
//! | ir == REF  | none   |
//! | ir == ALT  | swap REF↔ALT + swap GT 0↔1 |
//! | ir == rc(REF) | flip both alleles, keep GT |
//! | ir == rc(ALT) | flip + swap REF↔ALT + swap GT |
//! | none       | error  |
//!
//! In `flip` mode A/T and C/G pairs are skipped as unresolvable without a
//! dbSNP anchor. `flip-all` processes them.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::items_after_statements,
    clippy::too_many_lines,
    clippy::too_many_arguments
)]

mod fasta;
mod gt;
mod nucleotide;
mod processor;
mod stats;
mod vcf;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use rsomics_common::{Result, RsomicsError};

use fasta::{FaiEntry, load_contig, load_fai};
use processor::process_line_bytes;
use stats::{Stats, print_stats};
use vcf::{field, trim_newline};

/// The operation to apply for each record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixMode {
    /// Pass through unchanged, emit stats to stderr.
    Check,
    /// Swap/flip based on single-base matching; skip A/T and C/G pairs.
    Flip,
    /// Same as `Flip` but also processes ambiguous A/T and C/G pairs.
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
/// `fasta_path` must have an adjacent `.fai` index. Writes the (possibly
/// modified) VCF to `output` and emits stats to stderr.
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
