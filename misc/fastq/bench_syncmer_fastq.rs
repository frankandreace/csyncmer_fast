/// Benchmark syncmer throughput on real FASTQ data.
/// Measures raw closed syncmer detection excluding I/O overhead.
///
/// Usage: cargo run --release --example bench_syncmer_fastq -- [-k 31] [-w 1022] <input.fastq>
use packed_seq::{PackedSeqVec, SeqVec};
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut k: usize = 31;
    let mut w: usize = 1022;
    let mut filename: Option<&str> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-k" => { i += 1; k = args[i].parse().unwrap(); }
            "-w" => { i += 1; w = args[i].parse().unwrap(); }
            _ => filename = Some(&args[i]),
        }
        i += 1;
    }
    let filename = filename.expect("Usage: bench_syncmer_fastq [-k 31] [-w 1022] <input.fastq>");
    let syncmer_len = k + w - 1;

    eprintln!("k={k} w={w} syncmer_len={syncmer_len}");

    // mmap the FASTQ file
    let file = std::fs::File::open(filename).expect("cannot open file");
    let mmap = unsafe { memmap2::Mmap::map(&file).expect("mmap failed") };
    let data = &mmap[..];

    // Phase 1: scan FASTQ to collect sequence slices (no copy, just pointers into mmap)
    eprintln!("Scanning FASTQ...");
    let mut seqs: Vec<&[u8]> = Vec::with_capacity(2_000_000);
    let mut total_bp: usize = 0;
    let mut max_len: usize = 0;
    let mut min_len: usize = usize::MAX;

    let mut pos = 0;
    while pos < data.len() {
        if data[pos] != b'@' { break; }
        // line 1: header
        while pos < data.len() && data[pos] != b'\n' { pos += 1; }
        if pos >= data.len() { break; }
        pos += 1;
        // line 2: sequence
        let seq_start = pos;
        while pos < data.len() && data[pos] != b'\n' { pos += 1; }
        let seq_end = pos;
        let slen = seq_end - seq_start;
        if pos >= data.len() { break; }
        pos += 1;
        // line 3: +
        while pos < data.len() && data[pos] != b'\n' { pos += 1; }
        if pos >= data.len() { break; }
        pos += 1;
        // line 4: quality (same length as seq)
        pos += slen;
        if pos < data.len() && data[pos] == b'\n' { pos += 1; }

        seqs.push(&data[seq_start..seq_end]);
        total_bp += slen;
        if slen > max_len { max_len = slen; }
        if slen < min_len { min_len = slen; }
    }

    eprintln!(
        "  {} sequences, {:.2} Gbp, avg {} bp, min {}, max {}",
        seqs.len(),
        total_bp as f64 / 1e9,
        total_bp / seqs.len(),
        min_len,
        max_len
    );

    // Phase 2: run canonical closed syncmer detection (packing included in timing)
    // canonical requires odd syncmer_len; adjust w by +1 if needed
    if (k + w - 1) % 2 == 0 {
        eprintln!("Adjusting w from {w} to {} for odd syncmer_len (canonical requirement)", w + 1);
        w += 1;
    }
    let syncmer_len = k + w - 1;
    eprintln!("Running canonical closed syncmer detection (k={k}, w={w}, syncmer_len={syncmer_len})...");
    eprintln!("  (packing included in timing)");
    let mut positions: Vec<u32> = Vec::new();
    let mut total_syncmers: usize = 0;

    let t0 = Instant::now();
    for seq in &seqs {
        let pseq = PackedSeqVec::from_ascii(seq);
        positions.clear();
        simd_minimizers::canonical_closed_syncmers(k, w).run(pseq.as_slice(), &mut positions);
        total_syncmers += positions.len();
    }
    let elapsed = t0.elapsed().as_secs_f64();

    println!("Sequences:  {}", seqs.len());
    println!("Total bp:   {:.2} Gbp", total_bp as f64 / 1e9);
    println!("Avg length: {} bp", total_bp / seqs.len());
    println!("Syncmers:   {total_syncmers}");
    println!("Time:       {elapsed:.3} s");
    println!("Throughput: {:.2} GB/s", total_bp as f64 / elapsed / 1e9);
    println!(
        "Rate:       {:.1} M syncmers/s",
        total_syncmers as f64 / elapsed / 1e6
    );
}
