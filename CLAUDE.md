# csyncmer_fast

Fast closed syncmer detection using ntHash32.

## Project Structure

```
csyncmer_fast.h          # Main header: _fused_rescan_branchless + _simd_multiwindow only
misc/
  simd/                  # SIMD utilities (vec.hpp, nthash.hpp, hash_simd.hpp, sliding_min.hpp)
  older/                 # Legacy/redundant code (C++)
    legacy_infrastructure.hpp    # SeqHash, CircularArray, NtHash wrappers, iterators
    syncmer_seqhash.hpp          # SeqHash-based syncmer implementations
    syncmer_nthash32_variants.hpp # Redundant ntHash32 variants + hash benchmarks
  bench/                 # Benchmark suite
  scripts/               # Test and profiling scripts
  syng/                  # External syng library
```

## Implementations

### In csyncmer_fast.h (active)
| Function | Description |
|----------|-------------|
| `_fused_rescan_branchless` | Scalar fallback (~110-145 MB/s) |
| `_simd_multiwindow` | SIMD 8-window parallel (~140-175 MB/s) |

### In misc/older/ (archived C++)
| File | Contents |
|------|----------|
| `syncmer_nthash32_variants.hpp` | `_2bit_rescan`, `_2bit_deque`, `_direct_rescan`, `_fused_deque`, `_fused_rescan`, `_vanherk`, `_simd_rescan`, hash benchmarks |
| `syncmer_seqhash.hpp` | SeqHash-based implementations (64-bit hash) |
| `legacy_infrastructure.hpp` | SeqHash, CircularArray, NtHash C wrappers, iterators |

## Key Implementation Details

### ntHash Roll Signature
```cpp
uint32_t roll(uint8_t out_base, uint8_t in_base, uint8_t new_first_base)
```
- `out_base`: base leaving the window (unused in computation, encoded in state)
- `in_base`: base entering the window
- `new_first_base`: first base of the NEW s-mer (for state update)

The state `fw_` encodes `hash ^ f_rot[first_base]`, so rolling requires the new first base, not the outgoing base.

### SIMD_MULTIWINDOW Implementation
- Uses extended ring buffer (64 + 32 elements) to eliminate wrap-around branches
- Direct ASCII lookup tables skip 2-bit encoding step
- Processes 8 consecutive overlapping windows per iteration
- O(w) work per 8 k-mers, but SIMD parallelism makes it faster than O(1) amortized RESCAN

## Benchmarking
```bash
cd misc/bench
make clean && make
./benchmark check ~/data/human.chr19.fasta 31 15  # correctness
./benchmark bench ~/data/human.chr19.fasta 31 15 /tmp/bench.txt  # performance
```

## Performance (chr19, K=31, S=15)
- FUSED_RESCAN_BF: ~110-145 MB/s
- SIMD_MULTIWINDOW: ~140-175 MB/s
