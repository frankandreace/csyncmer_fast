/// Benchmark syncmer throughput on real FASTQ data.
/// Measures raw closed syncmer detection excluding I/O and packing overhead.
///
/// Usage: cargo run --release --example bench_syncmer_fastq -- [-k 31] [-w 1022] <input.fastq>
///
/// Two modes:
///   (default)  PackedNSeqVec + run_skip_ambiguous_windows — single call for all reads
///   -perread   Per-read loop with PackedSeqVec + run() — straightforward but slower
use packed_seq::{PackedNSeqVec, PackedSeqVec, SeqVec};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut k: usize = 31;
    let mut w: usize = 1022;
    let mut per_read = false;
    let mut filename: Option<&str> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-k" => { i += 1; k = args[i].parse().unwrap(); }
            "-w" => { i += 1; w = args[i].parse().unwrap(); }
            "-perread" => per_read = true,
            _ => filename = Some(&args[i]),
        }
        i += 1;
    }
    let filename = filename.expect("Usage: bench_syncmer_fastq [-k 31] [-w 1022] [-perread] <input.fastq>");

    // canonical requires odd syncmer_len; adjust w by +1 if needed
    if (k + w - 1) % 2 == 0 {
        eprintln!("Adjusting w from {w} to {} for odd syncmer_len (canonical requirement)", w + 1);
        w += 1;
    }
    let syncmer_len = k + w - 1;
    eprintln!("k={k} w={w} syncmer_len={syncmer_len}");

    // mmap the FASTQ file
    let file = std::fs::File::open(filename).expect("cannot open file");
    let mmap = unsafe { memmap2::Mmap::map(&file).expect("mmap failed") };
    let data = &mmap[..];

    // Phase 1: scan FASTQ to collect sequence slices
    eprintln!("Scanning FASTQ...");
    let mut seqs: Vec<&[u8]> = Vec::with_capacity(2_000_000);
    let mut total_bp: usize = 0;
    let mut max_len: usize = 0;
    let mut min_len: usize = usize::MAX;

    let mut pos = 0;
    while pos < data.len() {
        if data[pos] != b'@' { break; }
        while pos < data.len() && data[pos] != b'\n' { pos += 1; }
        if pos >= data.len() { break; }
        pos += 1;
        let seq_start = pos;
        while pos < data.len() && data[pos] != b'\n' { pos += 1; }
        let seq_end = pos;
        let slen = seq_end - seq_start;
        if pos >= data.len() { break; }
        pos += 1;
        while pos < data.len() && data[pos] != b'\n' { pos += 1; }
        if pos >= data.len() { break; }
        pos += 1;
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

    // PackedNSeqVec path requires k<=96 (par_iter_kmers limit); fall back to per-read
    let (total_syncmers, elapsed) = if per_read || syncmer_len > 96 {
        if syncmer_len > 96 && !per_read {
            eprintln!("syncmer_len={syncmer_len} > 96: falling back to per-read mode (simd-minimizers limit)");
        }
        bench_per_read(&seqs, k, w)
    } else {
        bench_nseq(&seqs, k, w)
    };

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

/// Original per-read approach: pack each read individually, call run() per read.
/// Straightforward but has per-call overhead (hasher construction, TLS access).
fn bench_per_read(seqs: &[&[u8]], k: usize, w: usize) -> (usize, f64) {
    eprintln!("Mode: per-read (packing excluded from timing)");

    // Pre-pack all sequences
    eprintln!("Pre-packing sequences...");
    let packed: Vec<PackedSeqVec> = seqs.iter()
        .map(|s| PackedSeqVec::from_ascii(s))
        .collect();

    let mut positions: Vec<u32> = Vec::new();
    let mut total_syncmers: usize = 0;

    let t0 = Instant::now();
    for pseq in &packed {
        positions.clear();
        simd_minimizers::canonical_closed_syncmers(k, w).run(pseq.as_slice(), &mut positions);
        total_syncmers += positions.len();
    }
    let elapsed = t0.elapsed().as_secs_f64();
    (total_syncmers, elapsed)
}

/// Batch approach: concatenate reads with N-separators into PackedNSeqVec chunks,
/// then call run_skip_ambiguous_windows per chunk. Amortizes per-call overhead.
/// Chunking needed because simd-minimizers uses u32 positions: len*8 must fit in 2^32.
fn bench_nseq(seqs: &[&[u8]], k: usize, w: usize) -> (usize, f64) {
    eprintln!("Mode: PackedNSeqVec batch (packing excluded from timing)");

    // Max bases per chunk: (2^32 / 8) - slack for separators
    const MAX_CHUNK_BP: usize = (1 << 29) - (1 << 20); // ~512 MB with margin

    // Build chunks of N-separated reads, pre-pack each
    eprintln!("Building N-separated chunks...");
    let mut chunks: Vec<PackedNSeqVec> = Vec::new();
    let mut concat = Vec::with_capacity(MAX_CHUNK_BP + 100_000);
    let mut chunk_bp: usize = 0;

    for seq in seqs {
        if chunk_bp > 0 && chunk_bp + seq.len() + 1 > MAX_CHUNK_BP {
            chunks.push(PackedNSeqVec::from_ascii(&concat));
            concat.clear();
            chunk_bp = 0;
        }
        if chunk_bp > 0 {
            concat.push(b'N');
        }
        concat.extend_from_slice(seq);
        chunk_bp += seq.len() + if chunk_bp > 0 { 1 } else { 0 };
    }
    if !concat.is_empty() {
        chunks.push(PackedNSeqVec::from_ascii(&concat));
    }

    eprintln!("  {} chunks", chunks.len());

    let mut positions: Vec<u32> = Vec::new();
    let mut total_syncmers: usize = 0;
    let builder = simd_minimizers::canonical_closed_syncmers(k, w);

    eprintln!("Running canonical closed syncmers...");
    let t0 = Instant::now();
    for chunk in &chunks {
        positions.clear();
        builder.run_skip_ambiguous_windows(chunk.as_slice(), &mut positions);
        total_syncmers += positions.len();
    }
    let elapsed = t0.elapsed().as_secs_f64();

    (total_syncmers, elapsed)
}
