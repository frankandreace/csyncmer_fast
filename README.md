# csyncmer_fast

[![Build Docker Image for build and test](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml)
[![C speed benchmark build](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml)
[![python package compilation](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml)

Header-only pure C library for fast closed syncmer extraction from DNA sequences. Uses ntHash rolling hashes with two algorithms: **twostack** (AVX2 SIMD batch processing) and **rescan** (portable scalar iterator).

### Quick Start

```c
#include "csyncmer_fast.h"

const char* seq = "ACGTACGTACGT...";
size_t len = strlen(seq);
size_t K = 31, S = 15;

// twostack: fastest, AVX2, ~99.99996% accurate
size_t count = csyncmer_twostack_simd_32_count(seq, len, K, S);

// Canonical twostack: strand-independent, AVX2
size_t canon_count = csyncmer_twostack_simd_32_canonical_count(seq, len, K, S);

// Canonical with positions and strands, AVX2
uint32_t positions[10000];
uint8_t strands[10000];  // 0=forward, 1=RC had minimal s-mer
size_t n = csyncmer_twostack_simd_32_canonical_positions(
    seq, len, K, S, positions, strands, 10000);

// Iterator: scalar, portable, exact results
CsyncmerIterator64* iter = csyncmer_iterator_create_64(seq, len, K, S);
size_t pos;
while (csyncmer_iterator_next_64(iter, &pos)) {
    // process syncmer at seq[pos..pos+K)
}
csyncmer_iterator_destroy_64(iter);

// Canonical iterator: strand-independent, scalar
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
| `csyncmer_twostack_simd_32_count` | Count | ~1390-1460 MB/s | AVX2, fastest |
| `csyncmer_twostack_simd_32_positions` | Positions | ~780-910 MB/s | AVX2 |
| `csyncmer_iterator_*_64` | Positions | ~380-430 MB/s | Scalar, portable, exact |

**Canonical (strand-independent):**

| Function | Output | Speed | Notes |
|----------|--------|-------|-------|
| `csyncmer_twostack_simd_32_canonical_count` | Count | ~1290 MB/s | AVX2 |
| `csyncmer_twostack_simd_32_canonical_positions` | Positions + strands | ~450-480 MB/s | AVX2 |
| `csyncmer_iterator_*_canonical_64` | Positions + strands | ~285-330 MB/s | Scalar, portable, exact |

All SIMD implementations use 16-bit hash approximation (~99.99996% accurate, ~4 errors per 10M syncmers).
Speeds measured on chr19 (59 MB), best-of-5, Intel Core Ultra 5 135H (4.6 GHz).

Performance is within 10-15% of [simd-minimizers](https://github.com/rust-seq/simd-minimizers) (Rust SIMD), faster for non-canonical and slightly slower for canonical due to `min(fw, rc)` strand tracking overhead vs simd-minimizers' `fw XOR rc`.

### Benchmarking

```bash
cd misc/tests
make
./benchmark --quick /path/to/sequence.fasta 31 15        # quick benchmark (best-of-3)
./benchmark /path/to/sequence.fasta 31 15 output.tsv     # full benchmark suite
```

### Closed Syncmers

A k-mer is a closed syncmer iff the minimal LEFTMOST s-mer (s < k) it contains is either at the first or last position (scanning left to right).

### Key Implementation Details

The library provides two algorithms for syncmer extraction:

**twostack** (AVX2 SIMD, batch): Splits the sequence into 8 chunks processed in parallel. Uses the two-stack (prefix-suffix) sliding minimum algorithm ([Groot Koerkamp & Martayan, SEA 2025](https://curiouscoding.nl/posts/fast-minimizers/); independently implemented in [simd-minimizers](https://github.com/rust-seq/simd-minimizers)) for O(1) amortized operations per element. Packs hash (upper 16 bits) + position (lower 16 bits) into 32-bit values for SIMD comparison. Requires the full sequence upfront.

**rescan** (scalar, streaming iterator): O(1) amortized — maintains a circular buffer and only rescans the window when the current minimum falls out (~6% of iterations). Branch-free minimum update using conditional moves. Works on any CPU without SIMD and processes one syncmer at a time, making it suitable for streaming pipelines.

**Hash precision**: The SIMD twostack path uses 16-bit hash approximation (upper 16 bits of a 32-bit ntHash). The iterator API (`csyncmer_iterator_*_64`) uses 64-bit hashes with no approximation.

**Canonical vs Forward-only**: Canonical implementations use `min(forward_hash, rc_hash)` for strand-independent results. Use canonical when you need consistent syncmer positions regardless of which DNA strand is sequenced.

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

- [Edgar RC (2021) "Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological sequences." PeerJ 9:e10805](https://peerj.com/articles/10805/) - Original syncmer definition
- [Groot Koerkamp & Martayan (SEA 2025) "SimdMinimizers: Computing random minimizers, fast."](https://github.com/rust-seq/simd-minimizers) - twostack sliding minimum algorithm
- [Curiouscoding: Fast Minimizers](https://curiouscoding.nl/posts/fast-minimizers/) - Blog post describing the twostack approach
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
