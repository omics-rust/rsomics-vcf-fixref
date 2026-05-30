/// A=0 C=1 G=2 T=3; sentinel for non-ACGT.
pub(crate) const NON_ACGT: u8 = u8::MAX;

#[inline]
pub(crate) fn nt2int(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => NON_ACGT,
    }
}

#[inline]
pub(crate) fn int2nt(x: u8) -> u8 {
    b"ACGT"[x as usize]
}

/// Reverse complement in integer space: A↔T (0↔3), C↔G (1↔2).
#[inline]
pub(crate) fn revint(x: u8) -> u8 {
    3 - x
}
