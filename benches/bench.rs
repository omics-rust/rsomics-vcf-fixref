use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_vcf_fixref(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-fixref");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let vcf = manifest.join("tests/golden/small.vcf");
    let fasta = manifest.join("tests/golden/ref.fa");
    c.bench_function("rsomics-vcf-fixref golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([vcf.to_str().unwrap(), "-f", fasta.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vcf_fixref);
criterion_main!(benches);
