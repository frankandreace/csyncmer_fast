# Scripts

Utility scripts for benchmarking and testing csyncmer_fast.

## Scripts Overview

| Script | Purpose |
|--------|---------|
| `test.sh` | Full benchmark suite - runs speed test 10x and plots results |
| `test_speed.sh` | Single benchmark run with CPU pinning for accurate timing |
| `test_speed_cluster.sh` | Simple speed test for cluster environments (no CPU manipulation) |
| `test_correctness.sh` | Fuzz testing - verifies all algorithms produce identical results |
| `profile.sh` | Deep profiling with perf, FlameGraph, valgrind massif/cachegrind |
| `plot_result.py` | Generates boxplot charts from benchmark TSV results |

## Usage

All scripts should be run from this directory (`misc/scripts/`).

### Running the Full Benchmark Suite
```bash
./test.sh
```
Runs `test_speed.sh` 10 times and generates plots in `../../benchmark/results/`.

### Single Speed Benchmark
```bash
./test_speed.sh                          # Use defaults
./test_speed.sh -f /path/to/file.fa      # Custom FASTA file
./test_speed.sh -k 31 -s 11              # Custom K and S values
```
Pins CPU to 2.6GHz, disables hyperthreading, runs benchmark, restores settings.

### Cluster Speed Benchmark
```bash
./test_speed_cluster.sh
```
Same as `test_speed.sh` but without CPU/SMT manipulation (for shared clusters).

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
Runs perf + FlameGraph, cache miss analysis, valgrind massif (memory), and cachegrind (cache simulation). Requires `perf`, `valgrind`, and FlameGraph tools.

### Plotting Results
```bash
./plot_result.py benchmark.tsv output_prefix
```
Reads benchmark TSV, generates `output_prefix.png` (all algorithms) and `output_prefix_focus.png` (syncmer implementations only).

## Prerequisites

- Built benchmark binary: `make` from repo root (produces `misc/bench/benchmark`)
- For `test_speed.sh` and `profile.sh`: sudo access for CPU frequency/SMT control
- For `profile.sh`: perf, valgrind, FlameGraph (expected at `../../FlameGraph/`)
- For `plot_result.py`: Python with pandas and matplotlib
