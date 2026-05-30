use crate::nucleotide::int2nt;

#[derive(Default, Debug)]
pub(crate) struct Stats {
    pub(crate) nsite: u64,
    pub(crate) nok: u64,
    pub(crate) nflip: u64,
    pub(crate) nswap: u64,
    pub(crate) nflip_swap: u64,
    pub(crate) nunresolved: u64,
    pub(crate) nerr: u64,
    pub(crate) nskip: u64,
    pub(crate) non_acgt: u64,
    pub(crate) non_snp: u64,
    pub(crate) non_biallelic: u64,
    /// Substitution type count indexed by `[ref_int][alt_int]`.
    pub(crate) count: [[u64; 4]; 4],
}

pub(crate) fn print_stats(stats: &Stats) {
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
