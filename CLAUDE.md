# csyncmer_fast

Fast closed syncmer detection using ntHash. Pure C header-only library with AVX2 SIMD support.

## Quick Start

```c
#include "csyncmer_fast.h"

// Batch API (fastest)
size_t count;
csyncmer_compute_fused_rescan_branchless_64(seq, len, K, S, &count);

// Iterator API (incremental)
CsyncmerIterator64* iter = csyncmer_iterator_create_64(seq, len, K, S);
size_t pos;
while (csyncmer_iterator_next_64(iter, &pos)) {
    // process syncmer at seq[pos..pos+K)
}
csyncmer_iterator_destroy_64(iter);
```

## Project Structure

```
csyncmer_fast.h          # Main header (pure C, header-only)
misc/
  bench/                 # Benchmark suite (benchmark.cpp)
  simd/                  # SIMD utilities (legacy C++, not used by main header)
  older/                 # Archived C++ implementations
  syng/                  # External syng library (seqhash)
  python/                # Python bindings
```

## Implementations in csyncmer_fast.h

### 64-bit ntHash (recommended)
| Function | Type | Speed | Description |
|----------|------|-------|-------------|
| `csyncmer_compute_fused_rescan_branchless_64` | Batch | ~150-200 MB/s | Scalar, O(1) amortized |
| `csyncmer_compute_simd_multiwindow_64` | Batch | ~130-185 MB/s | AVX2 8-lane parallel |
| `csyncmer_iterator_create/next/destroy_64` | Iterator | ~115-155 MB/s | Incremental, O(1) amortized |

### 32-bit ntHash (legacy)
| Function | Type | Speed | Description |
|----------|------|-------|-------------|
| `csyncmer_compute_fused_rescan_branchless` | Batch | ~110-145 MB/s | Scalar fallback |
| `csyncmer_compute_simd_multiwindow` | Batch | ~140-175 MB/s | AVX2 8-lane parallel |

## Key Algorithms

### RESCAN (Scalar)
- O(1) amortized: only rescans window when minimum falls out (~6% of iterations)
- Branch-free minimum update using conditional moves
- Circular buffer with power-of-2 size for fast modulo

### SIMD Multiwindow
- Processes 8 consecutive k-mers in parallel
- Extended ring buffer (64 + 32 elements) eliminates wrap-around branches
- O(window_size) per 8 k-mers, but parallelism compensates
- Uses lower 32 bits for comparison (64-bit version)

### Iterator
- Shared static lookup tables (~2KB saved per iterator)
- Local variable caching for register allocation
- Same RESCAN algorithm as batch

## Benchmarking

```bash
cd misc/bench
make clean && make
./benchmark check ~/data/human.chr19.fasta 31 15       # correctness
./benchmark quick ~/data/human.chr19.fasta 31 15       # quick perf test
./benchmark quick ~/data/human.chr19.fasta 31 15 64    # 64-bit only
./benchmark bench ~/data/human.chr19.fasta 31 15 out.txt  # full benchmark
```

## Notes

- Different hash sizes (32/64/128-bit) produce different syncmer counts due to different tie-breaking
- 64-bit scalar RESCAN is often faster than 64-bit SIMD due to O(1) vs O(w) complexity
- Iterator is ~75-80% of batch speed (inherent overhead from state save/restore)
