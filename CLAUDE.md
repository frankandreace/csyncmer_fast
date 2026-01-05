# csyncmer_fast

Fast closed syncmer detection using ntHash. Pure C header-only library with AVX2 SIMD support.

## Quick Start

```c
#include "csyncmer_fast.h"

// TWOSTACK: fastest (~550 MB/s), 32-bit hash, ~99.99996% accurate
size_t count;
csyncmer_compute_twostack_simd_32_count(seq, len, K, S, &count);

// RESCAN: exact results (~200 MB/s), 64-bit hash
csyncmer_compute_fused_rescan_branchless_64(seq, len, K, S, &count);

// Iterator API (incremental, exact)
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

### 32-bit ntHash TWOSTACK (fastest)
| Function | Type | Speed | Description |
|----------|------|-------|-------------|
| `csyncmer_compute_twostack_simd_32_count` | Batch | ~350-550 MB/s | Count only, AVX2 |
| `csyncmer_compute_twostack_simd_32` | Batch | ~190-250 MB/s | With positions, AVX2 |

### 64-bit ntHash (recommended for exact results)
| Function | Type | Speed | Description |
|----------|------|-------|-------------|
| `csyncmer_compute_fused_rescan_branchless_64` | Batch | ~140-210 MB/s | Scalar, O(1) amortized |
| `csyncmer_iterator_create/next/destroy_64` | Iterator | ~145-160 MB/s | Incremental, O(1) amortized |

### 32-bit ntHash (internal)
| Function | Type | Speed | Description |
|----------|------|-------|-------------|
| `csyncmer_compute_fused_rescan_branchless` | Batch | ~110-155 MB/s | TWOSTACK fallback |

## Key Algorithms

### RESCAN (Scalar)
- O(1) amortized: only rescans window when minimum falls out (~6% of iterations)
- Branch-free minimum update using conditional moves
- Circular buffer with power-of-2 size for fast modulo

### TWOSTACK (Fastest)
- Two-stack sliding minimum algorithm with O(1) amortized operations
- Splits sequence into 8 chunks processed in parallel with AVX2
- Packs hash (upper 16 bits) + position (lower 16 bits) for SIMD comparison
- Falls back to RESCAN for small inputs or window_size > 64

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

## Optimization Notes

### TWOSTACK Algorithm

The fastest implementation uses the simd-minimizers approach:
- Splits sequence into 8 independent chunks processed in parallel with AVX2
- Uses two-stack (prefix-suffix) sliding minimum algorithm: O(1) amortized
- Packs hash (upper 16 bits) + position (lower 16 bits) for SIMD comparison

**Limitation**: Uses 16-bit hash comparison, which can cause ~0.00004% error rate when two s-mers have identical upper 16 bits but different lower 16 bits. For exact results, use the 64-bit RESCAN implementation.

### RESCAN Algorithm

The 64-bit RESCAN approach provides exact results:
- O(1) amortized: only rescans window when minimum falls out (~6% of iterations)
- Branch-free minimum update using conditional moves
- Circular buffer with power-of-2 size for fast modulo

### Performance Comparison

| Implementation | Speed | Exact? | Use Case |
|----------------|-------|--------|----------|
| TWOSTACK (32-bit) | ~350-550 MB/s | ~99.99996% | High throughput |
| RESCAN (64-bit) | ~140-210 MB/s | 100% | Exact results |
| Iterator (64-bit) | ~145-160 MB/s | 100% | Incremental processing |

## Notes

- Different hash sizes (32/64/128-bit) produce different syncmer counts due to different tie-breaking
- TWOSTACK is fastest but uses 16-bit hash approximation; use RESCAN for exact results
- Iterator is ~75-80% of batch speed (inherent overhead from state save/restore)
- Deprecated implementations (SIMD multiwindow) are archived in `misc/older/deprecated_simd_mw.h`
