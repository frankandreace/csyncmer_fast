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

## Files

| File | Description |
|---|---|
| `test.cpp` | Correctness tests: unit tests + cross-implementation validation |
| `benchmark.cpp` | Performance benchmarks: quick summary + full suite with TSV output |
| `fasta_reader.h` | Shared header-only FASTA parser |
| `Makefile` | Build rules |
| `test_100kbp.fasta` | Small test sequence (100 kbp) |
| `results/` | Benchmark output (TSV + plots) |
