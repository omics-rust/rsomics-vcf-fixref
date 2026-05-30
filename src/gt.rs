use std::io::Write;

/// Write FORMAT+samples (from the 9th column onward) with GT alleles 0↔1 swapped.
///
/// Operates on byte slices to avoid allocations.
pub(crate) fn write_swapped_gt_column(
    fmt_samples: &[u8],
    w: &mut dyn Write,
) -> std::io::Result<()> {
    let tab_pos = fmt_samples.iter().position(|&b| b == b'\t');
    let format_col = match tab_pos {
        Some(p) => &fmt_samples[..p],
        None => fmt_samples,
    };

    let gt_idx = format_col.split(|&b| b == b':').position(|f| f == b"GT");

    let Some(gt_idx) = gt_idx else {
        return w.write_all(fmt_samples);
    };

    let samples_start = match tab_pos {
        Some(p) => p + 1,
        None => return w.write_all(fmt_samples),
    };

    w.write_all(format_col)?;

    let samples_bytes = &fmt_samples[samples_start..];
    for sample in samples_bytes.split(|&b| b == b'\t') {
        w.write_all(b"\t")?;
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
pub(crate) fn write_swapped_gt(gt: &[u8], w: &mut dyn Write) -> std::io::Result<()> {
    let sep = gt
        .iter()
        .position(|&b| b == b'/')
        .map(|p| (p, b'/'))
        .or_else(|| gt.iter().position(|&b| b == b'|').map(|p| (p, b'|')));
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
