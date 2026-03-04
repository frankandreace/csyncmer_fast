# csyncmer_fast

[![Build Docker Image for build and test](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml)
[![C speed benchmark build](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml)
[![Compile check](https://github.com/frankandreace/csyncmer_fast/actions/workflows/compile-check.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/compile-check.yml)
[![python package compilation](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml)

![Linux](https://img.shields.io/badge/Linux-supported-brightgreen) ![macOS](https://img.shields.io/badge/macOS-supported-brightgreen) ![Windows](https://img.shields.io/badge/Windows-supported-brightgreen)

Header-only pure C library for fast closed syncmer extraction from DNA sequences using ntHash rolling hashes. Processes human chr19 (59 MB) at up to **1.5 GB/s** with AVX2 SIMD acceleration, and remains fully functional on all platforms via automatic scalar fallback.

### Quick Start

```c
#include "csyncmer_fast.h"

const char* seq = "ACGTACGTACGT...";
size_t len = strlen(seq);
size_t K = 31, S = 15;

// Batch API — forward-only (single strand)
// ~1400 MB/s with AVX2 SIMD, ~400 MB/s scalar fallback, ~99.99996% accurate
size_t count = csyncmer_twostack_simd_32_count(seq, len, K, S);

uint32_t positions[10000];
size_t n = csyncmer_twostack_simd_32_positions(
    seq, len, K, S, positions, 10000);

// Batch API — canonical (strand-independent)
uint8_t strands[10000];  // 0=forward, 1=RC had minimal s-mer
size_t cn = csyncmer_twostack_simd_32_canonical_positions(
    seq, len, K, S, positions, strands, 10000);

// Streaming API — forward-only
// Scalar, portable, exact 64-bit hashes (~400 MB/s)
CsyncmerIterator64* iter = csyncmer_iterator_create_64(seq, len, K, S);
size_t pos;
while (csyncmer_iterator_next_64(iter, &pos)) {
    // process syncmer at seq[pos..pos+K)
}
csyncmer_iterator_destroy_64(iter);

// Streaming API — canonical (strand-independent)
CsyncmerIteratorCanonical64* citer = csyncmer_iterator_create_canonical_64(seq, len, K, S);
int strand;
while (csyncmer_iterator_next_canonical_64(citer, &pos, &strand)) {
    // strand: 0=forward, 1=RC had minimal s-mer
}
csyncmer_iterator_destroy_canonical_64(citer);
```

Compile and run the example:
```bash
cc -std=c11 -O3 -march=native -o example example.c
./example
```

### Performance

Measured on human chr19 (59 MB), best-of-5, Intel Core Ultra 5 135H (4.6 GHz).

**K=31, S=15:**

| Variant | AVX2 SIMD (`-march=native`) | Scalar (portable) |
|---|---|---|
| **Forward-only** | | |
| Count | **1390 MB/s** | 400 MB/s |
| Positions | **790 MB/s** | 400 MB/s |
| **Canonical** | | |
| Count | **1290 MB/s** | 330 MB/s |
| Positions | **500 MB/s** | 330 MB/s |
| **Streaming iterator (64-bit)** | | |
| Forward positions | — | 430 MB/s |
| Canonical positions | — | 330 MB/s |

**K=21, S=11:**

| Variant | AVX2 SIMD (`-march=native`) | Scalar (portable) |
|---|---|---|
| **Forward-only** | | |
| Count | **1460 MB/s** | 400 MB/s |
| Positions | **910 MB/s** | 400 MB/s |
| **Canonical** | | |
| Count | **1290 MB/s** | 370 MB/s |
| Positions | **460 MB/s** | 285 MB/s |
| **Streaming iterator (64-bit)** | | |
| Forward positions | — | 380 MB/s |
| Canonical positions | — | 285 MB/s |

The batch API (`csyncmer_twostack_simd_32_*`) is the recommended entry point. On x86 with AVX2 (`-march=native`), it uses SIMD-accelerated twostack processing (8 parallel chunks, 16-bit hash packing, ~99.99996% accurate). Without AVX2, it automatically falls back to scalar rescan with exact 32-bit hashes. The streaming iterator provides exact 64-bit hashes and processes one syncmer at a time.

Performance is within 10-15% of [simd-minimizers](https://github.com/rust-seq/simd-minimizers) (Rust SIMD), faster for non-canonical and slightly slower for canonical due to `min(fw, rc)` strand tracking overhead vs simd-minimizers' `fw XOR rc`.

### Benchmarking

```bash
cd misc/tests
make
./benchmark --quick /path/to/sequence.fasta 31 15        # quick benchmark (best-of-3)
./benchmark /path/to/sequence.fasta 31 15 output.tsv     # full benchmark suite
```

Performance results are automatically generated in the [CI job summary](../../actions/workflows/c_speed_bench_build.yml) on every push.

### FASTQ Multi-Read Mode

Special code path for FASTQ reads (e.g. HiFi). This multi-read mode processes 8 reads simultaneously in AVX2 lanes, eliminating the ~60% warmup overhead that single-read SIMD has on kilobase-sized reads with kilobase-sized `-w` values.

```bash
cd misc/fastq
make
./bench_syncmer_fastq ~/data/reads.fastq                         # multi-8 monolithic (default)
./bench_syncmer_fastq -twopass ~/data/reads.fastq                # two-pass: hash then twostack (with strand)
./bench_syncmer_fastq -twopass-nostrand ~/data/reads.fastq       # two-pass: hash then twostack (no strand, fastest)
./bench_syncmer_fastq -hashonly ~/data/reads.fastq               # hash pass only (returns s-mer count)
./bench_syncmer_fastq -packonly ~/data/reads.fastq               # pack+interleave only (I/O bound)
./bench_syncmer_fastq -single ~/data/reads.fastq                 # single-read SIMD (baseline)
./bench_syncmer_fastq -k 31 -w 1022 ~/data/reads.fastq          # custom k/w
```

**Modes:**
- **default** (`multi-8-bucketed`): Monolithic pass — hash computation and twostack sliding minimum fused together. Canonical (tracks strand).
- **`-twopass`**: Split into hash-only pass (writes to intermediate buffer) then twostack pass. Canonical. Same syncmer output as default.
- **`-twopass-nostrand`**: Like `-twopass` but skips strand tracking. Fastest mode (~1.9 GB/s). Produces identical syncmer positions.
- **`-hashonly`** / **`-packonly`**: Isolate individual pipeline stages for profiling.
- **`-single`**: One-read-at-a-time SIMD (baseline for measuring multi-read speedup).

All syncmer-producing modes must agree on count (correctness check). The two-pass split enables independent optimization of hash computation and sliding minimum.

A reference Rust benchmark (`bench_syncmer_fastq.rs`) using [simd-minimizers](https://github.com/rust-seq/simd-minimizers) is included for cross-library comparison.

### Closed Syncmers

A k-mer is a closed syncmer iff the minimal LEFTMOST s-mer (s < k) it contains is either at the first or last position (scanning left to right).

### Architecture

The library is organized around two APIs backed by three internal algorithms:

**Batch API** (`csyncmer_twostack_simd_32_*`) — pass the full sequence, get all syncmer positions at once. Internally selects the best available algorithm:
- **twostack** (AVX2 SIMD): splits the sequence into 8 chunks processed in parallel using the two-stack sliding minimum algorithm ([Groot Koerkamp & Martayan, SEA 2025](https://curiouscoding.nl/posts/fast-minimizers/)). Packs hash (upper 16 bits) + position (lower 16 bits) into 32-bit values for SIMD comparison. ~99.99996% accurate (~4 errors per 10M syncmers).
- **rescan** (scalar fallback): O(1) amortized — maintains a circular buffer and only rescans the window when the current minimum falls out (~6% of iterations). Branch-free minimum update using conditional moves. Exact 32-bit hashes.

**Streaming API** (`csyncmer_iterator_*_64`) — create/next/destroy pattern, one syncmer at a time. Uses exact 64-bit ntHash with no approximation. Suitable for streaming pipelines where the full sequence is not available upfront or memory is constrained.

**Canonical vs Forward-only**: Canonical implementations use `min(forward_hash, rc_hash)` for strand-independent results. Use canonical when you need consistent syncmer positions regardless of which DNA strand is sequenced.

### Project Structure

```
csyncmer_fast.h              # Main header: pure C, header-only
example.c                    # Standalone usage example
misc/
  tests/                     # Benchmark/test suite
  code/                      # Reference implementations (32/64/128-bit variants)
    syncmer_seqhash.hpp           # SeqHash (64-bit) implementations
    syncmer_nthash32.hpp          # ntHash32 algorithm variants
    syncmer_nthash64.hpp          # ntHash64 implementations
    syncmer_nthash128.hpp         # ntHash128 implementations
    legacy_infrastructure.hpp     # SeqHash, CircularArray, etc.
    simd/                         # C++ SIMD utilities
  fastq/                     # FASTQ multi-read mode (8 reads in 8 AVX2 lanes)
    csyncmer_fastq.h              # Multi-read 8-lane monolithic SIMD
    csyncmer_fastq_hash.h         # Hash-only pass (two-pass mode)
    csyncmer_fastq_pack.h         # 2-bit DNA packing + interleaving
    csyncmer_fastq_split.h        # Twostack-only pass (two-pass mode, no strand)
    csyncmer_fastq_twostack.h     # Twostack pass with strand tracking
    bench_syncmer_fastq.c         # C benchmark (all modes)
    bench_syncmer_fastq.rs        # Rust benchmark (simd-minimizers, for reference)
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
- [x] Cross-platform compilation (Linux, macOS, Windows)
- [ ] ARM NEON SIMD acceleration

**Note on ARM NEON**: NEON registers are 128-bit (4×u32), half the width of AVX2 (256-bit, 8×u32). Emulating 8 lanes with 2×uint32x4_t doubles every SIMD instruction. While Apple M-series out-of-order execution partially hides this overhead for count-only variants, position collection suffers from the lack of native 256-bit operations (movemask, gather, permute). In our experiments, an 8-lane NEON twostack implementation achieved ~460 MB/s for count (vs ~500 MB/s scalar rescan) and ~270 MB/s for positions (vs ~460 MB/s rescan) — the SIMD overhead exceeds the algorithmic benefit. The scalar rescan fallback already performs well on ARM (~400-500 MB/s), making a dedicated NEON path low priority until wider SIMD extensions (SVE/SVE2) become widespread.
