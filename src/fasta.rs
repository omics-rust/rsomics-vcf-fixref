use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Arc;

use rsomics_common::{Result, RsomicsError};

pub(crate) struct FaiEntry {
    pub(crate) length: u64,
    /// Byte offset of the first base in the FASTA file.
    pub(crate) offset: u64,
    pub(crate) line_bases: u64,
    pub(crate) line_width: u64,
}

/// Load a .fai index into a name → entry map.
///
/// Keys are `Vec<u8>` for zero-copy lookup against the VCF CHROM field.
pub(crate) fn load_fai(fai_path: &Path) -> Result<HashMap<Vec<u8>, FaiEntry>> {
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
/// Loading once per contig eliminates per-record seek syscalls — the dominant
/// cost for chromosome-sorted VCF inputs.
pub(crate) fn load_contig(fasta: &mut File, entry: &FaiEntry) -> Option<Arc<Vec<u8>>> {
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

