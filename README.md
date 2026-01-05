# csyncmer_fast
Header-only pure C library for fast closed syncmer detection using ntHash and AVX2 SIMD.

### CI/CD
[![Build Docker Image for build and test](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml)

[![C speed benchmark build](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml)

[![python package compilation](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml)

### Project Structure

```
.
├── csyncmer_fast.h          # Main header: pure C, header-only
└── misc/
    ├── bench/               # Benchmark suite (benchmark.cpp)
    ├── older/               # Archived implementations
    │   ├── deprecated_simd_mw.h          # Deprecated SIMD multiwindow
    │   ├── syncmer_seqhash.hpp           # Original 64-bit SeqHash syncmers
    │   ├── syncmer_nthash32_variants.hpp # ntHash32 algorithm variants
    │   └── legacy_infrastructure.hpp     # SeqHash, CircularArray, etc.
    ├── simd/                # SIMD utilities (legacy C++)
    ├── python/              # Python bindings
    └── syng/                # External seqhash library
```

### Performance (chr19, K=31, S=15)

| Implementation | Speed | Exact? | Notes |
|----------------|-------|--------|-------|
| **TWOSTACK (32-bit)** | ~350-550 MB/s | ~99.99996% | AVX2 two-stack, fastest |
| **RESCAN (64-bit)** | ~140-210 MB/s | 100% | Scalar, O(1) amortized |
| **Iterator (64-bit)** | ~145-160 MB/s | 100% | Incremental API |

TWOSTACK uses 16-bit hash comparison for speed, causing rare (~4 per 10M) tie-breaking errors.
Use 64-bit RESCAN for exact results.

### Quick Start

Compile the benchmark:
```bash
cd misc/bench
make
```

Run quick performance test:
```bash
./benchmark quick /path/to/sequence.fasta 31 15
```

Run correctness verification:
```bash
./benchmark check /path/to/sequence.fasta 31 15
```

Run full benchmark to file:
```bash
./benchmark bench /path/to/sequence.fasta 31 15 output.tsv
```

### API Usage

```c
#include "csyncmer_fast.h"

const char* seq = "ACGTACGTACGT...";
size_t len = strlen(seq);
size_t K = 31, S = 15;
size_t count;

// TWOSTACK: fastest (~550 MB/s), 32-bit hash, ~99.99996% accurate
csyncmer_compute_twostack_simd_32_count(seq, len, K, S, &count);

// RESCAN: exact results (~200 MB/s), 64-bit hash
csyncmer_compute_fused_rescan_branchless_64(seq, len, K, S, &count);

// Iterator API: incremental processing, exact results
CsyncmerIterator64* iter = csyncmer_iterator_create_64(seq, len, K, S);
size_t pos;
while (csyncmer_iterator_next_64(iter, &pos)) {
    // process syncmer at seq[pos..pos+K)
}
csyncmer_iterator_destroy_64(iter);
```

### Key Implementation Details

**TWOSTACK Algorithm**: Splits sequence into 8 chunks processed in parallel with AVX2. Uses two-stack (prefix-suffix) sliding minimum for O(1) amortized operations. Packs hash (upper 16 bits) + position (lower 16 bits) for SIMD comparison.

**RESCAN Algorithm**: O(1) amortized - only rescans window when minimum falls out (~6% of iterations). Branch-free minimum update using conditional moves. Circular buffer with power-of-2 size for fast modulo.

**ntHash Rolling Hash**: Both 32-bit and 64-bit variants. Uses direct ASCII lookup tables to skip 2-bit encoding overhead.

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

- [simd-minimizers](https://github.com/Daniel-Liu-c0deb0t/simd-minimizers) - TWOSTACK algorithm inspiration
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
- [x] Pure C implementation (no C++ dependency)
- [x] 64-bit ntHash support
- [x] Iterator API for incremental processing
- [ ] ARM NEON support
