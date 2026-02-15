# csyncmer_fast

[![Build Docker Image for build and test](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml)
[![C speed benchmark build](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml)
[![python package compilation](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml)

Header-only pure C library for fast closed syncmer extraction from DNA sequences. Uses forward-strand ntHash rolling hashes and AVX2 SIMD to process sequences at up to ~600 MB/s.

### Quick Start

```c
#include "csyncmer_fast.h"

const char* seq = "ACGTACGTACGT...";
size_t len = strlen(seq);
size_t K = 31, S = 15;

// TWOSTACK: fastest (~600 MB/s), AVX2, ~99.99996% accurate
size_t count = csyncmer_twostack_simd_32_count(seq, len, K, S);

// Canonical TWOSTACK: strand-independent (~550 MB/s), AVX2
size_t canon_count = csyncmer_twostack_simd_32_canonical_count(seq, len, K, S);

// Canonical with positions and strands (~180 MB/s), AVX2
uint32_t positions[10000];
uint8_t strands[10000];  // 0=forward, 1=RC had minimal s-mer
size_t n = csyncmer_twostack_simd_32_canonical_positions(
    seq, len, K, S, positions, strands, 10000);

// Iterator: scalar, portable, exact results (~165 MB/s)
CsyncmerIterator64* iter = csyncmer_iterator_create_64(seq, len, K, S);
size_t pos;
while (csyncmer_iterator_next_64(iter, &pos)) {
    // process syncmer at seq[pos..pos+K)
}
csyncmer_iterator_destroy_64(iter);

// Canonical iterator: strand-independent (~130 MB/s), scalar
CsyncmerIteratorCanonical64* citer = csyncmer_iterator_create_canonical_64(seq, len, K, S);
int strand;
while (csyncmer_iterator_next_canonical_64(citer, &pos, &strand)) {
    // strand: 0=forward, 1=RC had minimal s-mer
}
csyncmer_iterator_destroy_canonical_64(citer);
```

Compile and run the example:
```bash
gcc -std=c11 -o example -march=native example.c
./example
```

### API

**Forward-only (single strand):**

| Function | Output | Speed | Notes |
|----------|--------|-------|-------|
| `csyncmer_twostack_simd_32_count` | Count | ~550-635 MB/s | AVX2, fastest |
| `csyncmer_twostack_simd_32_positions` | Positions | ~250-290 MB/s | AVX2 |
| `csyncmer_iterator_*_64` | Positions | ~145-165 MB/s | Scalar, portable, exact |

**Canonical (strand-independent):**

| Function | Output | Speed | Notes |
|----------|--------|-------|-------|
| `csyncmer_twostack_simd_32_canonical_count` | Count | ~500-570 MB/s | AVX2 |
| `csyncmer_twostack_simd_32_canonical_positions` | Positions + strands | ~170-200 MB/s | AVX2 |
| `csyncmer_iterator_*_canonical_64` | Positions + strands | ~100-130 MB/s | Scalar, portable, exact |

All SIMD implementations use 16-bit hash approximation (~99.99996% accurate, ~4 errors per 10M syncmers).

### Benchmarking

```bash
cd misc/tests
make
./benchmark quick /path/to/sequence.fasta 31 15   # performance test
./benchmark check /path/to/sequence.fasta 31 15  # correctness verification
./benchmark bench /path/to/sequence.fasta 31 15 output.tsv  # full benchmark
```

### Closed Syncmers

A k-mer is a closed syncmer iff the minimal LEFTMOST s-mer (s < k) it contains is either at the first or last position (scanning left to right).

### Key Implementation Details

**TWOSTACK Algorithm** (AVX2): Splits sequence into 8 chunks processed in parallel. Uses two-stack (prefix-suffix) sliding minimum for O(1) amortized operations. Packs hash (upper 16 bits) + position (lower 16 bits) for SIMD comparison. Uses 16-bit hash approximation (~4 errors per 10M syncmers).

**Iterator** (Scalar): Uses RESCAN algorithm internally - O(1) amortized, only rescans window when minimum falls out (~6% of iterations). Branch-free minimum update. Works on any CPU without SIMD.

**ntHash Rolling Hash**: Both 32-bit and 64-bit variants (forward-strand only, not canonical). Uses direct ASCII lookup tables to skip 2-bit encoding overhead.

**Canonical vs Forward-only**: Canonical implementations use `min(forward_hash, rc_hash)` for strand-independent results. Use canonical when you need consistent syncmer positions regardless of which DNA strand is sequenced. Forward-only is faster when strand doesn't matter or you're processing both strands separately.

### Project Structure

```
csyncmer_fast.h              # Main header: pure C, header-only
misc/
  tests/                     # Benchmark/test suite
  code/                      # Reference implementations (32/64/128-bit variants)
    syncmer_seqhash.hpp           # SeqHash (64-bit) implementations
    syncmer_nthash32.hpp          # ntHash32 algorithm variants
    syncmer_nthash64.hpp          # ntHash64 implementations
    syncmer_nthash128.hpp         # ntHash128 implementations
    legacy_infrastructure.hpp     # SeqHash, CircularArray, etc.
    simd/                         # C++ SIMD utilities
  python/                    # Python bindings
  syng/                      # External seqhash library
```

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
- [csyncmer_fast_benchmark](https://github.com/frankandreace/csyncmer_fast_benchmark) - Comparisons against other syncmer libraries

### TODO

- [x] Add unit tests for rescan
- [x] Implement ntHash32 variants
- [x] Add AVX2 SIMD syncmer computation
- [x] Pure C implementation (no C++ dependency)
- [x] 64-bit ntHash support
- [x] Iterator API for incremental processing
- [x] Canonical (strand-independent) support
- [x] AVX2 SIMD canonical implementation
- [ ] ARM NEON support
