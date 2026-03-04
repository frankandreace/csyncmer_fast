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
- **Legacy 64-bit code uses 32-bit comparison**: `misc/code/syncmer_nthash64.hpp`
  functions (`csyncmer_nthash64_rescan_count` at line ~305, `csyncmer_nthash64_simd_multiwindow`
  at line ~433) compute full 64-bit ntHash but compare only `(uint32_t)hash` (lower 32 bits).
  This produces ~0.5-1.6% different counts from the true 64-bit iterator API in csyncmer_fast.h.
  The test suite prints `[NOTE] 64-bit legacy implementations use 32-bit comparison internally`.
  These are reference/legacy implementations only — the public API in csyncmer_fast.h is correct

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

## FASTQ Benchmark (misc/fastq)

Standalone syncmer throughput benchmark on real FASTQ data (no hash table / I/O overhead).

```bash
cd misc/fastq
make clean && make

# All syncmer-producing modes must agree on count (e.g. 16709730 for SRR34765324.20G)
./bench_syncmer_fastq ~/data/SRR34765324.20G.fastq                   # monolithic (default)
./bench_syncmer_fastq -twopass ~/data/SRR34765324.20G.fastq           # two-pass with strand
./bench_syncmer_fastq -twopass-nostrand ~/data/SRR34765324.20G.fastq  # split two-pass, no strand
./bench_syncmer_fastq -hashonly ~/data/SRR34765324.20G.fastq          # hash pass only (returns s-mer count, not syncmers)
./bench_syncmer_fastq -packonly ~/data/SRR34765324.20G.fastq          # pack+interleave only
./bench_syncmer_fastq -single ~/data/SRR34765324.20G.fastq            # single-sequence API
```

Data: `~/data/SRR34765324.20G.fastq` (cichlid HiFi, 9.9 Gbp, 1.3M reads, avg 7376 bp).
Default params: k=31 w=1022 (K=1052 S=31).

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
  fastq/                 # FASTQ syncmer throughput benchmark (bench_syncmer_fastq)
```

## simd-minimizers Comparison (2026-03-04)

Analyzed `~/tools/simd-minimizers/src/canonical_minimizers_simd.cpp` (~1100 lines) for optimization ideas applicable to the FASTQ twostack. Same algorithmic family: 8-lane AVX2 twostack sliding-window minimum with hash|pos packing.

**Nothing worth copying.** Key findings:

1. **Their sliding min does NOT have our pre-load optimization** — same store→load disambiguation issue we fixed. They store to ring, then load next element, same cache line.

2. **SIMD compaction via Lemire's LUT (`UNIQSHUF[256]`)** — 256-entry shuffle mask table for compressing scattered lane hits into contiguous output. Not useful for syncmers: at 0.3% hit rate, scalar `tzcnt` loop is faster (branch predictor nearly always predicts "no hit").

3. **8x8 transpose for per-lane batching** — transposes 8 registers into per-lane output arrays. Not applicable: we use `movemask` + `tzcnt` to scatter directly.

4. **Fused hash+twostack single pass** — saves hash_buf I/O but loses profiling granularity and mixes machine clear sources between the two passes. Not worth the debugging cost.

5. **`prefix_min = elem` after rescan** (vs our `UINT32_MAX`) — saves 1 `vpminud` per 1022 iterations. Negligible.

6. **Clang being faster for them** is likely C++-specific (templates, STL inlining). Not transferable to our C code — GCC 15 beats Clang 20 by 8% on our workload (5.37s vs 5.80s, Clang emits 1.63x more instructions).

## Machine Clears Optimization (2026-03-04)

**Problem**: VTune uarch-exploration showed **24% of pipeline slots** lost to "Machine Clears: Other Nukes" in the twopass-nostrand pipeline. The twostack function alone wasted ~45% of its slots on clears.

**Root cause**: Store at `ring_ptr`, then load at `ring_ptr+1` (after increment) — same cache line. CPU memory disambiguator sees store followed by load at pointer-derived address, speculatively assumes aliasing, triggers machine clear (~32 cycles per misprediction).

**Fix** (in `csyncmer_fastq_split.h`):

1. **Pre-load suf_min BEFORE the ring store**: `__m256i suf_min_preload = *(ring_ptr + 1)` executes before `*ring_ptr = elem`. On non-wrap iterations (~99.9%), the load reads a value committed long ago — no disambiguation issue. On wrap iterations, reload after rescan.

2. **Anti-4K-aliasing pad**: Added 2048B between hash_buf and ring_buf in work buffer layout to break low-12-bit address collisions (CPU compares only bits [11:0] for store-forwarding).

3. **+1 ring padding element**: Allocated `window_size + 1` elements, initialized to `UINT32_MAX`, so the pre-load at `ring_end - 1` doesn't read out of bounds.

**Results** (GCC 15, pinned to P-core):

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| Other Nukes (overall) | 24.0% | 19.6% | **-18%** |
| Twostack clears (per-function) | ~45% | 36.1% | **-20%** |
| Hash pass clears | ~35% | ~35% | unchanged |
| Best wall time | — | 5.37s / 1.84 GB/s | — |

Disassembly confirms GCC emits load before store in the hot loop (vmovdqa load at 0x7f25, store at 0x7f46).

### Deep dive: remaining clears (2026-03-04)

`machine_clears.disambiguation` exists on Meteor Lake E-cores (cpu_atom) but NOT P-cores (cpu_core). Ran full breakdown on E-core (taskset -c 16):

| Counter | Count | % |
|---------|-------|---|
| any (total) | 905K | 100% |
| page_fault | 311K | 34% |
| memory_ordering | 245K | 27% |
| **disambiguation** | **8K** | **0.9%** |
| mrn_nuke | 5K | 0.5% |
| smc/fp_assist | ~0 | 0% |
| *unaccounted* | ~335K | 37% |

**Pre-load fix confirmed effective**: disambiguation is 8K out of 905K — essentially zero.

Per-function VTune (P-core): twostack 37.6% Other Nukes, hash 33.2%, main 0.0%. Twostack has 0% memory_ordering — ALL its clears are "Other Nukes." But these are NOT disambiguation (E-core confirmed). The 37% unaccounted clears on E-core have no public PMU sub-counter.

Total clears: 1.29M over 33B P-core cycles (0.004%). The per-function 37.6% is a fraction of twostack's own small slot count, not overall impact. Page faults (34%) are from mmap I/O (kernel, not user-space functions).

**Conclusion**: User-space machine clears are solved. Remaining "Other Nukes" in twostack/hash are a P-core microarchitectural effect without a public counter, and total clear count is negligible.

## VTune Cheat Sheet

**Setup**: `source /opt/intel/oneapi/setvars.sh --force 2>/dev/null`

**Collect** (pinned to P-core):
```bash
rm -rf /tmp/vtune_result
taskset -c 0 vtune -collect uarch-exploration -knob sampling-interval=0.1 \
  -result-dir /tmp/vtune_result -- ./bench_syncmer_fastq -twopass-nostrand ~/data/SRR34765324.20G.fastq
```

**Per-function TMA** (use column names with `(%)` suffix, NO `:Self`):
```bash
vtune -report hotspots -r /tmp/vtune_result -format csv -csv-delimiter '|' \
  -group-by function \
  -column "Function,CPU Time,Performance-core (P-core):Retiring(%),Performance-core (P-core):Bad Speculation:Machine Clears:Other Nukes(%)"
```

**Per-source-line TMA** (use `:Self` suffix on ALL TMA columns, filter by source-function):
```bash
vtune -report hotspots -r /tmp/vtune_result -format csv -csv-delimiter ',' \
  -group-by source-line -source-object bench_syncmer_fastq \
  -filter "source-function=csyncmer_twostack_only_multi" \
  -column "Clockticks:Performance-core (P-core),Performance-core (P-core):Retiring:Self,Performance-core (P-core):Front-End Bound:Self,Performance-core (P-core):Bad Speculation:Self,Performance-core (P-core):Bad Speculation:Machine Clears:Other Nukes:Self,Performance-core (P-core):Back-End Bound:Self,Performance-core (P-core):Back-End Bound:Memory Bound:Self,Performance-core (P-core):Back-End Bound:Core Bound:Self,Source Line,Source File"
```

**Key gotchas**:
- Function-level grouping: column names use `(%)` suffix, NO `:Self`
- Source-line grouping: column names MUST have `:Self` suffix
- Filter uses `source-function=` (not `Function=`)
- VTune attributes inlined SIMD intrinsics to their header files (avx2intrin.h etc.), not to caller source lines. Use `source-function` filter + sort by Clockticks to see hot lines.
- Sort output by `Clockticks:Performance-core (P-core)` column for meaningful ordering (CPU Time rounds to 0.0 for individual lines)
- `machine_clears.disambiguation` counter exists only on E-cores (cpu_atom), not P-cores (cpu_core). Run on E-core with `taskset -c 16` for full breakdown.

### Per-source-line TMA analysis (twopass-nostrand, 2026-03-04)

**Twostack** (top lines by clockticks):

| Clockticks | Line | Code | Retiring | Bad Spec | BE Bound |
|-----------|------|------|---------|----------|----------|
| 1394M | 144 | rescan loop | 43% | 37% nukes | 18% (core) |
| 1005M | 111 | main loop back-edge | 25% | 39% nukes | 32% (mem+core) |
| 1000M | 122 | suf_min preload | **59%** | 11% | 26% (mem) |
| 984M | 174 | while(km) scatter | 23% | 65% nukes | 11% |
| 964M | 136 | pos/fs/ls += ones | **54%** | 0% | **41%** (port) |
| 924M | 166 | lane limit check | 34% | 54% nukes | 10% |

- Line 122 (preload): 59% retiring = most efficient hot line, optimization working
- Line 136 (3x vpaddd): 41% back-end/core = port-limited, 3 ALU ops competing
- Line 144 (rescan): 37% nukes from sequential store→load in rescan loop
- Lines 174-178 (scatter): high nuke % but low total clockticks (0.3% hit rate)

**Hash** (essentially 2 lines):

| Clockticks | Line | Code | Retiring | Bad Spec | BE Bound |
|-----------|------|------|---------|----------|----------|
| 2744M | 290 | HASH_ONE_ITER body | 32% | 28% nukes | **40%** (core) |
| 1272M | 289 | loop control | 28% | 43% nukes | 28% (core) |

Hash is **40% core-bound** — port contention from dense SIMD (rotl7/rotr7 + 4x vpermd LUT lookups + XOR chain). This is the fundamental ALU throughput bottleneck.
