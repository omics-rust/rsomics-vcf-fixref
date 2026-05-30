use std::collections::HashMap;
use std::io::Write;

use crate::FixMode;
use crate::gt::write_swapped_gt_column;
use crate::nucleotide::{NON_ACGT, int2nt, nt2int, revint};
use crate::stats::Stats;
use crate::vcf::{field, memchr, nth_tab, parse_u64, trim_newline};

/// Process one VCF data line (raw bytes, newline already stripped by caller).
///
/// Writes the (possibly modified) line to `w` followed by `\n`, updating `stats`.
pub(crate) fn process_line_bytes(
    line: &[u8],
    contig_cache: Option<&[u8]>,
    mode: FixMode,
    stats: &mut Stats,
    skipped_chroms: &mut HashMap<Vec<u8>, bool>,
    annotate: bool,
    w: &mut dyn Write,
) -> std::io::Result<()> {
    stats.nsite += 1;

    let tab3 = nth_tab(line, 2);
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

    if memchr(b',', alt_col) {
        stats.non_biallelic += 1;
        stats.nskip += 1;
        return write_line_annotated(line, b"skip", annotate, w);
    }

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
    // In flip mode bcftools emits FIXREF=skip for unresolvable records; in
    // flip-all mode it emits FIXREF=err (no ambiguous-pair short-circuit in
    // that branch means errors reach a distinct code path).
    let err_tag = if mode == FixMode::Flip {
        b"skip" as &[u8]
    } else {
        b"err"
    };
    write_line_annotated(line, err_tag, annotate, w)
}

fn write_modified_line(
    line: &[u8],
    new_ref: u8,
    new_alt: u8,
    swap_gt: bool,
    action: &[u8],
    annotate: bool,
    w: &mut dyn Write,
) -> std::io::Result<()> {
    let t2 = nth_tab(line, 2).unwrap();
    w.write_all(&line[..t2 + 1])?;
    w.write_all(&[new_ref, b'\t'])?;
    w.write_all(&[new_alt, b'\t'])?;
    let t4 = nth_tab(line, 4).unwrap();
    let t6 = nth_tab(line, 6).unwrap();
    w.write_all(&line[t4 + 1..t6 + 1])?;
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
    if let Some(t7) = nth_tab(line, 7) {
        w.write_all(b"\t")?;
        let fmt_samples = trim_newline(&line[t7 + 1..]);
        if swap_gt {
            write_swapped_gt_column(fmt_samples, w)?;
        } else {
            w.write_all(fmt_samples)?;
        }
    }
    w.write_all(b"\n")
}

/// Write a VCF line pass-through, optionally appending `FIXREF=<action>` to INFO.
pub(crate) fn write_line_annotated(
    line: &[u8],
    action: &[u8],
    annotate: bool,
    w: &mut dyn Write,
) -> std::io::Result<()> {
    if !annotate {
        w.write_all(line)?;
        w.write_all(b"\n")?;
        return Ok(());
    }
    let t6 = match nth_tab(line, 6) {
        Some(i) => i,
        None => {
            w.write_all(line)?;
            w.write_all(b"\n")?;
            return Ok(());
        }
    };
    let info_start = t6 + 1;
    let info_end = nth_tab(line, 7).unwrap_or(line.len());
    let info = &line[info_start..info_end];

    w.write_all(&line[..info_start])?;
    if info == b"." {
        write!(w, "FIXREF={}", String::from_utf8_lossy(action))?;
    } else {
        w.write_all(info)?;
        write!(w, ";FIXREF={}", String::from_utf8_lossy(action))?;
    }
    if info_end < line.len() {
        w.write_all(&line[info_end..])?;
    }
    w.write_all(b"\n")
}

#[inline]
fn base_from_cache(contig: &[u8], pos0: u64) -> Option<u8> {
    let b = *contig.get(pos0 as usize)?;
    let v = nt2int(b);
    if v == NON_ACGT { None } else { Some(v) }
}
