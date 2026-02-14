# Scripts

Utility scripts for benchmarking and testing csyncmer_fast.

## Scripts Overview

| Script | Purpose |
|--------|---------|
| `test.sh` | Full benchmark suite - runs speed test 10x and plots results |
| `test_speed.sh` | Single benchmark run with optional CPU pinning |
| `test_correctness.sh` | Fuzz testing - verifies all algorithms produce identical results |
| `profile.sh` | Deep profiling with perf, FlameGraph, valgrind |
| `plot_result.py` | Generates boxplot charts from benchmark TSV results |

## Usage

All scripts should be run from this directory (`misc/scripts/`).

### Single Speed Benchmark
```bash
./test_speed.sh                          # With CPU pinning (requires sudo)
./test_speed.sh -c                       # Cluster mode (no sudo, no CPU manipulation)
./test_speed.sh -f /path/to/file.fa      # Custom FASTA file
./test_speed.sh -k 31 -s 11              # Custom K and S values
```

Options:
- `-f FILE` - FASTA file to benchmark
- `-k SIZE` - K-mer size (default: 31)
- `-s SIZE` - S-mer size (default: 11)
- `-c` - Cluster mode: skip CPU frequency/SMT manipulation (no sudo required)

### Full Benchmark Suite
```bash
./test.sh
```
Runs `test_speed.sh` 10 times and generates plots in `../tests/results/`.

### Correctness Testing
```bash
./test_correctness.sh
```
Generates 10 random 10kb sequences with random K/S values, verifies all algorithm implementations produce the same syncmer count.

### Profiling
```bash
./profile.sh                             # Use defaults
./profile.sh -f /path/to/file.fa         # Custom FASTA file
```
Runs perf + FlameGraph, cache miss analysis, valgrind massif (memory), and cachegrind (cache simulation).

### Plotting Results
```bash
./plot_result.py benchmark.tsv output_prefix
```
Reads benchmark TSV, generates `output_prefix.png` (all algorithms) and `output_prefix_focus.png` (syncmer implementations only).

## Output

All benchmark results are saved to `misc/tests/results/`.

## Prerequisites

- Built benchmark binary: `make` in `misc/tests/`
- For `test_speed.sh` (without -c): sudo access for CPU frequency/SMT control
- For `profile.sh`: perf, valgrind, FlameGraph (expected at `../../FlameGraph/`)
- For `plot_result.py` (optional): Python with pandas and matplotlib
