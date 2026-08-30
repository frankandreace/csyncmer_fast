# Tests & Benchmarks

Correctness tests and performance benchmarks for `csyncmer_fast.h`.

## Requirements

- GCC/G++ with C++17 and AVX2 support
- [ntHash](https://github.com/bcgsc/ntHash) library (`-lnthash`)
- zlib (`-lz`)
- A FASTA file for correctness and benchmark tests

A small test file (`test_100kbp.fasta`) is included for quick correctness checks.

## Build

```bash
make          # optimized build
make debug    # with -g, -Wall, -Wextra, assertions
```

## Correctness tests (`test`)

```bash
# Unit tests only (no FASTA file needed)
./test

# Unit tests + cross-implementation validation
./test FASTA_FILE K S
./test test_100kbp.fasta 31 15
```

Without arguments, runs built-in unit tests (base encoding, canonical strand independence, canonical hash values).

With a FASTA file, additionally validates that all implementations agree on syncmer counts:

- **seqhash-based** (6 implementations): naive, circular array, large array, deque, iterator, branchless
- **ntHash128** (3 implementations): generator, deque, naive
- **ntHash32** (11 implementations): naive, 2bit-rescan, 2bit-deque, direct-rescan, fused-deque, fused-rescan, SIMD-rescan, van-herk, TWOSTACK SIMD (count + positions)
- **ntHash64** (5 implementations): naive, iterator, rescan-count, SIMD MW
- **Canonical** (3 implementations): ntHash64 iterator, ntHash64 deque, ntHash32 TWOSTACK SIMD (count + positions)

Different hash sizes produce different syncmer counts (due to tie-breaking). The test verifies agreement within each hash family, and allows ~0.00004% tolerance for the 16-bit TWOSTACK approximation.

## Benchmarks (`benchmark`)

### Quick benchmark

Runs each `csyncmer_fast.h` function once, reports syncmer count and throughput:

```bash
./benchmark --quick FASTA_FILE K S [filter]
```

`filter` selects which implementations to run:
- `all` (default) - 64-bit iterators + all 32-bit SIMD variants
- `32` - 32-bit SIMD only
- `64` - 64-bit iterators only

Example:
```bash
./benchmark --quick ~/data/human.chr19.fasta 31 15
./benchmark --quick ~/data/human.chr19.fasta 31 15 32
```

### Full benchmark

Runs all implementations (including reference/legacy ones) with timing, writes a TSV for plotting:

```bash
./benchmark FASTA_FILE K S OUTPUT.tsv
```

Output lines tagged `[[HASHING ...]]` are pure hashing speed (no syncmer logic). Lines tagged `[[SYNCMERS ...]]` include the full syncmer detection pipeline. The tag also shows the hash type (`syng`, `nth32`, `nth64`, `nth128`) and algorithm variant.

## Reproducing the paper figure

Collects the throughput panels of Figure 1 -- CHM13v2.0 (k=31 s=15) and HiFi
reads (k=1052 s=31) -- and regenerates the plot:

```bash
../paper/run_benchmarks.sh > data.tsv
python3 ../paper/plot_throughput.py data.tsv -o fig_throughput.pdf
```

Dataset paths are set at the top of `run_benchmarks.sh`. Requires the
`simd-minimizers` and `digest` binaries in addition to this directory's
`benchmark` and `misc/fastq/bench_syncmer_fastq`.


## Throughput across (k, s) parameters

Canonical closed-syncmer detection with position output, measured on CHM13v2.0
(3.1 Gbp) for a range of commonly used (k, s) values. Throughput in GB/s;
`w = k - s + 1` is the s-mer window length.

| k, s | w | rescan | twostack | multi-8 | simd-minimizers |
|---|---|---|---|---|---|
| 15, 7 | 9 | 0.278 | 0.462 | 0.527 | 0.521 |
| 15, 11 | 5 | 0.206 | 0.437 | 0.494 | 0.482 |
| 21, 11 | 11 | 0.282 | 0.459 | 0.546 | 0.552 |
| 21, 15 | 7 | 0.241 | 0.447 | 0.519 | 0.497 |
| 31, 15 | 17 | 0.334 | 0.514 | 0.618 | 0.566 |
| 31, 19 | 13 | 0.316 | 0.478 | 0.578 | 0.514 |
| 31, 23 | 9 | 0.280 | 0.477 | 0.529 | 0.519 |
| 41, 21 | 21 | 0.342 | 0.502 | 0.631 | 0.597 |
| 51, 31 | 21 | 0.342 | 0.517 | 0.620 | 0.584 |

Reproduce with:

```bash
./sweep_ks.sh ~/data/chm13v2.0.fa > sweep_chm13.tsv
```

Same measurement on HiFi reads (SRR34765324, 9.9 Gbp):

| k, s | w | rescan | twostack | multi-8 | simd-minimizers |
|---|---|---|---|---|---|
| 15, 11 | 5 | 0.19 | 0.43 | 0.61 | 0.56 |
| 15, 7 | 9 | 0.26 | 0.49 | 0.67 | 0.59 |
| 21, 11 | 11 | 0.27 | 0.46 | 0.70 | 0.61 |
| 21, 15 | 7 | 0.22 | 0.43 | 0.62 | 0.57 |
| 31, 15 | 17 | 0.30 | 0.49 | 0.76 | 0.65 |
| 31, 19 | 13 | 0.30 | 0.49 | 0.69 | 0.60 |
| 31, 23 | 9 | 0.25 | 0.47 | 0.67 | 0.61 |
| 41, 21 | 21 | 0.32 | 0.53 | 0.78 | 0.65 |
| 51, 31 | 21 | 0.32 | 0.53 | 0.82 | 0.63 |
| 285, 31 | 255 | 0.43 | 0.87 | 1.70 | 0.60 |
| 541, 31 | 511 | 0.45 | 0.80 | 1.72 | 0.54 |
| 1052, 31 | 1022 | 0.46 | 0.69 | 1.95 | 0.43 |

```bash
./sweep_ks_hifi.sh ~/data/SRR34765324.20G.fastq > sweep_hifi.tsv
```


Measured on an Intel Core Ultra 5 135H, GCC 15.2.0, rustc 1.92.0
(simd-minimizers 2.3.0, `-C target-cpu=native`), single-threaded.


## Files

| File | Description |
|---|---|
| `test.cpp` | Correctness tests: unit tests + cross-implementation validation |
| `benchmark.cpp` | Performance benchmarks: quick summary + full suite with TSV output |
| `fasta_reader.h` | Shared header-only FASTA parser |
| `Makefile` | Build rules |
| `test_100kbp.fasta` | Small test sequence (100 kbp) |
| `sweep_ks.sh` | (k, s) parameter sweep driver (writes TSV) |
| `results/` | Benchmark output (TSV + plots) |
