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
- Non-canonical positions use batch+transpose with `CSYNCMER_UNIQSHUF` LUT
  (SIMD shuffle compaction via `csyncmer_append_filtered`)
- Canonical positions use per-iteration scatter (batch+transpose provides no
  benefit because the RC hash computation dominates, not position collection)
- Count-only variants are ~2x faster than position-collecting variants
- Different hash sizes (32/64-bit) produce different syncmer counts due to tie-breaking

## Benchmark Speeds (chr19, 59 MB, best-of-5)

| Variant | K=21 S=11 | K=31 S=15 |
|---|---|---|
| Non-canonical count | ~1400 MB/s | ~1400 MB/s |
| Non-canonical positions | ~910 MB/s | ~790 MB/s |
| Canonical count | ~1260 MB/s | ~1260 MB/s |
| Canonical positions | ~460 MB/s | ~500 MB/s |

### vs simd-minimizers Rust (same chr19 data)

| Variant | csyncmer_fast | simd-minimizers Rust |
|---|---|---|
| Non-canonical K=21/S=11 | 913 MB/s | 637 MB/s (**+43%**) |
| Non-canonical K=31/S=15 | 791 MB/s | 711 MB/s (**+11%**) |
| Canonical K=21/S=11 | 461 MB/s | 534 MB/s (-14%) |
| Canonical K=31/S=15 | 497 MB/s | 560 MB/s (-11%) |

Non-canonical counts match exactly between the two projects. Canonical counts
differ slightly because of different canonical hash semantics (see below).

### Why Canonical Positions Are Slower Than simd-minimizers

The ~13% gap is **fundamental to the `min(fw, rc)` canonical hash approach**,
not a code quality issue. csyncmer_fast uses `min(forward_hash, rc_hash)` which
requires strand tracking (cmpgt + movemask + bitwise ops propagated through the
TWOSTACK prefix/suffix min — ~14 extra SIMD ops/iteration). simd-minimizers uses
`forward_hash XOR rc_hash` (strand-symmetric), eliminating all strand tracking.

Investigated optimizations (all reverted as not worth the trade-offs):
- **Embed strand in packed value bit 16**: Eliminates explicit strand tracking but
  reduces hash precision to 15 bits, causing count≠positions disagreement (~0.02%)
  and canonical count regression (~6%). Canonical positions improved only ~3-5%.
- **Batch+transpose for canonical positions**: Same approach as non-canonical.
  Marginal improvement because the bottleneck is the RC hash computation in the
  main loop, not the position collection method.

The per-iteration cycle breakdown explains the bound:
- Non-canonical count: ~26 cycles/iter → ~1400 MB/s
- Canonical count: ~31 cycles/iter (+5 for RC hash) → ~1260 MB/s
- Non-canonical positions: ~42 cycles/iter (+16 for batch+transpose) → ~910 MB/s
- Canonical positions: ~77 cycles/iter (+31 for RC hash + strand tracking + scatter) → ~460 MB/s

### Benchmarking Pitfall: Dead Code Elimination

When the return value of `csyncmer_twostack_simd_32_count` is discarded, GCC `-O3` creates a stripped `.isra` clone that eliminates the entire TWOSTACK algorithm (0 AVX2 instructions vs 153 in the real function). This produces bogus speeds of 1200-1500 MB/s. Always use `volatile` or otherwise consume the result when benchmarking.

## Canonical Hash: min(fw,rc) vs XOR (simd-minimizers)

csyncmer_fast and simd-minimizers use fundamentally different canonical hash strategies:

| | csyncmer_fast | simd-minimizers |
|---|---|---|
| **Canonical hash** | `min(fw_hash, rc_hash)` | `fw_hash XOR rc_hash` |
| **Strand tracking** | Required (cmpgt + movemask + bitwise propagation through TWOSTACK) | Not needed (XOR is symmetric) |
| **Sliding min** | Single minimum (leftmost) | Both leftmost AND rightmost simultaneously |
| **Strand selection** | Implicit in min() | TG-count majority vote + blend |
| **Position collection** | Non-canonical: batch+transpose; Canonical: per-iteration scatter | Batch+transpose for both |

These produce different canonical syncmer sets (e.g., 13,537,894 vs 13,505,495
for K=21/S=11 on chr19). Neither is "wrong" — they define canonicality differently.

**Why TG-count can't help csyncmer_fast**: simd-minimizers uses TG-count to
cheaply select a preferred strand (a boolean), then picks left or right minimum
accordingly. This works because XOR hash is strand-symmetric — you don't need to
know the strand to compute the hash. csyncmer_fast needs `min(forward_hash,
rc_hash)` as the actual hash value. Since ntHash is position-dependent, TG-count
can't predict which hash is smaller — both hashes must be computed and compared.

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
