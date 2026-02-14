# csyncmer_fast

Fast closed syncmer detection using ntHash. Pure C header-only library with AVX2 SIMD support.

## Quick Start

```c
#include "csyncmer_fast.h"

// TWOSTACK: fastest (~550 MB/s), AVX2, ~99.99996% accurate
size_t count = csyncmer_compute_twostack_simd_32_count(seq, len, K, S);

// Iterator: scalar, portable, exact results (~160 MB/s)
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
  tests/                 # Tests and benchmarks (test.cpp, benchmark.cpp)
  legacy/                # Deprecated implementations (organized by hash type)
  simd/                  # SIMD utilities (legacy C++)
  syng/                  # External syng library (seqhash)
  python/                # Python bindings
```

## Implementations in csyncmer_fast.h

### 32-bit ntHash TWOSTACK (fastest)
| Function | Type | Speed | Description |
|----------|------|-------|-------------|
| `csyncmer_compute_twostack_simd_32_count` | Batch | ~350-550 MB/s | Count only, AVX2 |
| `csyncmer_compute_twostack_simd_32` | Batch | ~190-250 MB/s | With positions, AVX2 |

### 64-bit ntHash (scalar, portable, exact)
| Function | Type | Speed | Description |
|----------|------|-------|-------------|
| `csyncmer_iterator_create/next/destroy_64` | Iterator | ~145-160 MB/s | Scalar, incremental, O(1) amortized |

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

## Testing & Benchmarking

```bash
cd misc/tests
make clean && make

# Correctness tests (uses included 100kbp test file)
./test test_100kbp.fasta 31 15

# Quick performance benchmark
./benchmark --quick ~/data/human.chr19.fasta 31 15
./benchmark --quick ~/data/human.chr19.fasta 31 15 64    # 64-bit only

# Full benchmark suite
./benchmark ~/data/human.chr19.fasta 31 15 out.txt
```

## Optimization Notes

### TWOSTACK Algorithm

The fastest implementation uses the simd-minimizers approach:
- Splits sequence into 8 independent chunks processed in parallel with AVX2
- Uses two-stack (prefix-suffix) sliding minimum algorithm: O(1) amortized
- Packs hash (upper 16 bits) + position (lower 16 bits) for SIMD comparison

**Limitation**: Uses 16-bit hash comparison, which can cause ~0.00004% error rate when two s-mers have identical upper 16 bits but different lower 16 bits. For exact results, use the 64-bit RESCAN implementation.

### RESCAN Algorithm

The 64-bit Iterator uses the RESCAN algorithm internally:
- O(1) amortized: only rescans window when minimum falls out (~6% of iterations)
- Branch-free minimum update using conditional moves
- Circular buffer with power-of-2 size for fast modulo

### Performance Comparison

| Implementation | Speed | Exact? | Use Case |
|----------------|-------|--------|----------|
| TWOSTACK (32-bit) | ~350-550 MB/s | ~99.99996% | High throughput, AVX2 required |
| Iterator (64-bit) | ~145-160 MB/s | 100% | Exact results, scalar, portable |

## Notes

- Different hash sizes (32/64-bit) produce different syncmer counts due to different tie-breaking
- TWOSTACK is fastest but requires AVX2 and uses 16-bit hash approximation
- Iterator is scalar (no SIMD) - works on any CPU including ARM
- Deprecated implementations are archived in `misc/legacy/`
