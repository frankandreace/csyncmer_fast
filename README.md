# csyncmer_fast
Header-only C++ library for fast closed syncmer extraction using ntHash32 and AVX2 SIMD.

### CI/CD
[![Build Docker Image for build and test](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml)

[![C speed benchmark build](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml)

[![python package compilation](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml)

### Project Structure

```
.
├── csyncmer_fast.h          # Main header: fastest ntHash32 implementations
└── misc/
    ├── simd/                # Reusable SIMD building blocks
    │   ├── vec.hpp          # u32x8 AVX2 wrapper with clean API
    │   ├── nthash.hpp       # ntHash32 rolling hash implementation
    │   ├── hash_simd.hpp    # SIMD hash utilities and gather ops
    │   └── sliding_min.hpp  # Sliding window minimum algorithms
    ├── older/               # Archived implementations for reference
    │   ├── syncmer_seqhash.hpp           # Original 64-bit SeqHash syncmers
    │   ├── syncmer_nthash32_variants.hpp # All ntHash32 algorithm variants
    │   └── legacy_infrastructure.hpp     # SeqHash, CircularArray, etc.
    ├── bench/               # Benchmark binary and Makefile
    ├── python/              # Python bindings and packaging
    ├── scripts/             # Benchmark and test scripts
    └── syng/                # Seqhash library (external dependency)
```

### Performance (chr19, K=31, S=15)

| Implementation | Speed | Notes |
|----------------|-------|-------|
| **SIMD_MULTIWINDOW** | ~140-175 MB/s | AVX2 8-lane parallel, fastest |
| **FUSED_RESCAN_BF** | ~110-145 MB/s | Scalar fallback, branch-free |
| SeqHash RESCAN | ~85-90 MB/s | Legacy 64-bit hash |
| SeqHash DEQUE | ~47 MB/s | Legacy deque-based |
| SeqHash NAIVE | ~42 MB/s | O(w) per k-mer baseline |

### Quick Start

Compile the benchmark:
```bash
cd misc/bench
make
```

Run correctness tests:
```bash
./benchmark check /path/to/sequence.fasta 31 15
```

Run performance benchmark:
```bash
./benchmark bench /path/to/sequence.fasta 31 15 output.tsv
```

### API Usage

```cpp
#include "csyncmer_fast.h"

using namespace csyncmer_fast_simd;

const char* seq = "ACGTACGTACGT...";
size_t num_syncmers;

// Fastest: AVX2 SIMD (falls back to scalar if <8 k-mers)
compute_closed_syncmers_nthash32_simd_multiwindow(
    seq, strlen(seq), K, S, &num_syncmers);

// Scalar fallback with branch-free minimum updates
compute_closed_syncmers_nthash32_fused_rescan_branchless(
    seq, strlen(seq), K, S, &num_syncmers);
```

### Key Implementation Details

**ntHash32 Rolling Hash**: 32-bit variant of ntHash optimized for SIMD. Uses direct ASCII lookup tables to skip 2-bit encoding overhead.

**SIMD Multi-Window**: Processes 8 consecutive overlapping windows in parallel using AVX2. Ring buffer with extended mirror area eliminates wrap-around branches.

**RESCAN Algorithm**: O(w) amortized per k-mer. Rescans window only when minimum falls out, with branch-free conditional updates otherwise.

### Closed Syncmers

A k-mer is a closed syncmer iff the minimal LEFTMOST s-mer (s < k) it contains is either at the first or last position (scanning left to right).

### Python Bindings

```bash
cd misc/python
python -m build
```

Or using Docker:
```bash
cd misc/python
docker-compose run python-build
```

### References

- [Curiouscoding: Fast Minimizers](https://curiouscoding.nl/posts/fast-minimizers/)
- [Sliding Window Minimum Algorithm](https://github.com/keegancsmith/Sliding-Window-Minimum)
- [SYNG (Durbin)](https://github.com/richarddurbin/syng/)
- [Strobealign](https://github.com/ksahlin/strobealign)
- [Minimap2](https://github.com/lh3/minimap2)

### Comprehensive Benchmarks

See [csyncmer_fast_benchmark](https://github.com/frankandreace/csyncmer_fast_benchmark) for comparisons against other syncmer libraries.

### TODO

- [x] Add unit tests for rescan
- [x] Implement ntHash32 variants
- [x] Add AVX2 SIMD syncmer computation
- [ ] Pure C implementation (no C++ dependency)
- [ ] ARM NEON support
