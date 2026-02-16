# csyncmer_fast - Internal Development Notes

See README.md for user-facing documentation.

## Testing & Benchmarking

```bash
cd misc/tests
make clean && make

# Correctness tests
./test test_100kbp.fasta 31 15

# Quick performance benchmark
./benchmark --quick ~/data/human.chr19.fasta 31 15
./benchmark --quick ~/data/human.chr19.fasta 31 15 64    # 64-bit only

# Full benchmark suite
./benchmark ~/data/human.chr19.fasta 31 15 out.txt
```

## Key Algorithms

### TWOSTACK (Fastest)
- Two-stack sliding minimum with O(1) amortized operations
- Splits sequence into 8 chunks processed in parallel with AVX2
- Packs hash (upper 16 bits) + position (lower 16 bits) for SIMD comparison
- Falls back to RESCAN for small inputs (num_kmers < 64) or window_size > 64.
  Note: the fallback uses full 32-bit hash resolution, unlike the SIMD path's
  16-bit packing — inconsequential for small inputs, theoretically matters for
  large windows on long sequences.
- Canonical variant computes both forward and RC hashes in parallel

### RESCAN (Scalar Iterator)
- O(1) amortized: only rescans window when minimum falls out (~6% of iterations)
- Branch-free minimum update using conditional moves
- Circular buffer with power-of-2 size for fast modulo

### Canonical Implementation
- Pre-computed `rotr7(RC[base])` table eliminates one rotation per iteration
- `prev_remove_base` tracking for correct RC rolling direction
- ILP restructuring: forward and RC computations interleaved

## Implementation Details

### Rolling Hash Formulas

**Forward rolling:**
```c
fw_state = fw_hash ^ f_rot[old_first]
fw_hash = rotl7(fw_state) ^ F[new_last]
```

**RC rolling:**
```c
rc_hash = rotr7(rc_hash) ^ rotr7(RC[old_first]) ^ c_rot[new_last]
```

### SIMD Lookup Tables

```c
CSYNCMER_SIMD_F32[8]       // F[A,C,T,G] repeated 2x for permutevar
CSYNCMER_SIMD_RC32[8]      // RC[A,C,T,G] = F[T,G,A,C] repeated 2x
CSYNCMER_SIMD_RC32_ROTR7[8] // Pre-computed rotr7(RC[base])
```

### Chunk Validity Fix

When `chunk_ends[i] < chunk_starts[i]` (boundary chunks), compute `last_lane_limit` as minimum across all valid lanes to avoid counting invalid k-mers.

## Performance Notes

- 16-bit hash approximation causes ~0.00004% error rate
- Canonical SIMD is 4-5x faster than scalar canonical iterator
- Count-only variants are ~2x faster than position-collecting variants
- Different hash sizes (32/64-bit) produce different syncmer counts due to tie-breaking

## Benchmark Speeds

### Realistic TWOSTACK Speeds

- **~550-600 MB/s** cold (first call, data not in CPU cache)
- **~630-660 MB/s** warm (after repeated calls, data in L2/L3)
- Binary size (.text) has negligible impact: a 7KB minimal binary and the 214KB full benchmark binary produce the same speeds when measured correctly

### Benchmarking Pitfall: Dead Code Elimination

When the return value of `csyncmer_twostack_simd_32_count` is discarded, GCC `-O3` creates a stripped `.isra` clone that eliminates the entire TWOSTACK algorithm (0 AVX2 instructions vs 153 in the real function). This produces bogus speeds of 1200-1500 MB/s. Always use `volatile` or otherwise consume the result when benchmarking.

## TG-count Optimization (Investigated, Not Applicable)

The simd-minimizers library uses TG-count for strand selection:
- Counts (TG - AC) bases in each window
- If TG > AC, the RC strand is "preferred" for consistent minimizer reporting
- This is for **strand preference**, not hash computation

For canonical syncmers, this approach doesn't apply because:
- We need `min(forward_hash, rc_hash)` as the actual hash value
- ntHash is position-dependent, so TG-count can't predict which hash is smaller
- We must compute both hashes to find the minimum s-mer position

**Key insight**: simd-minimizers' "canonical" is a boolean (which strand to report), while our "canonical" is a hash value (minimum of two hashes).

## Code Generation Architecture

csyncmer_fast.h uses two parameterized macros to generate all function variants:

- `CSYNCMER_DEFINE_RESCAN_32(FUNC_NAME, CANONICAL)` — generates scalar RESCAN.
  Count vs positions handled via runtime `if(out_positions)` — compiler eliminates
  dead branch when inlined with NULL.
- `CSYNCMER_DEFINE_TWOSTACK_SIMD_32(FUNC_NAME, CANONICAL, COLLECT_POSITIONS)` —
  generates SIMD TWOSTACK. `COLLECT_POSITIONS` must be compile-time because of
  static `CSYNCMER_THREAD_LOCAL` lane buffers (use `##FUNC_NAME` token pasting
  for unique TLS names per expansion).

All public API functions are thin inline wrappers calling the macro-generated internals.

## CI / Docker

- `misc/python/Dockerfile` builds the Docker image used by `c_speed_bench_build.yml`
- The Docker image MUST include ntHash (built from source) because misc/tests/
  {benchmark,test}.cpp depend on `nthash/nthash.hpp` and `-lnthash` via
  `legacy_infrastructure.hpp`. Removing ntHash from the Dockerfile breaks CI.
- Docker image rebuild triggers on Dockerfile changes; bench CI triggers on
  changes to `csyncmer_fast.h`, `misc/tests/**`, `misc/syng/**`.
  If both trigger on the same push, the bench CI may use the stale image —
  re-run manually after the Docker build completes.

## Project Structure

```
csyncmer_fast.h          # Main header (pure C, header-only)
example.c                # Standalone usage example (compile with gcc -mavx2)
misc/
  tests/                 # Tests and benchmarks (test.cpp, benchmark.cpp)
  code/                  # Reference implementations (32/64/128-bit variants)
    simd/                # C++ SIMD utilities (used by reference implementations)
  syng/                  # External syng library (seqhash)
  python/                # Python bindings + Dockerfile for CI
```
