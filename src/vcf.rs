/// Find the byte offset of the `n`-th tab (0-indexed) in `line`.
#[inline]
pub(crate) fn nth_tab(line: &[u8], n: usize) -> Option<usize> {
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

/// Field `n` (0-indexed) of a tab-split VCF line, as a byte slice.
#[inline]
pub(crate) fn field(line: &[u8], n: usize) -> &[u8] {
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

#[inline]
pub(crate) fn trim_newline(s: &[u8]) -> &[u8] {
    let mut e = s.len();
    if e > 0 && s[e - 1] == b'\n' {
        e -= 1;
    }
    if e > 0 && s[e - 1] == b'\r' {
        e -= 1;
    }
    &s[..e]
}

#[inline]
pub(crate) fn memchr(needle: u8, haystack: &[u8]) -> bool {
    haystack.contains(&needle)
}

#[inline]
pub(crate) fn parse_u64(s: &[u8]) -> Option<u64> {
    if s.is_empty() {
        return None;
    }
    let mut n: u64 = 0;
    for &b in s {
        if !b.is_ascii_digit() {
            return None;
        }
        n = n.checked_mul(10)?.checked_add((b - b'0') as u64)?;
    }
    Some(n)
}
