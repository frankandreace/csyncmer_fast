#pragma once

#include "vec.hpp"
#include "nthash.hpp"
#include <array>
#include <cstdint>

namespace csyncmer_simd {

/// SIMD rotl32 - rotate left by R bits (7)
/// Since R=7 is constant, we can optimize this
inline u32x8 simd_rotl7(u32x8 x) {
#if defined(CSYNCMER_SIMD_NEON)
    // NEON: shift left and OR with shift right
    uint32x4_t lo_shl = vshlq_n_u32(x.lo, 7);
    uint32x4_t hi_shl = vshlq_n_u32(x.hi, 7);
    uint32x4_t lo_shr = vshrq_n_u32(x.lo, 25);  // 32 - 7 = 25
    uint32x4_t hi_shr = vshrq_n_u32(x.hi, 25);
    return {vorrq_u32(lo_shl, lo_shr), vorrq_u32(hi_shl, hi_shr)};
#elif defined(CSYNCMER_SIMD_AVX2)
    // AVX2: _mm256_or and shifts
    __m256i shl = _mm256_slli_epi32(x.v, 7);
    __m256i shr = _mm256_srli_epi32(x.v, 25);
    return u32x8{_mm256_or_si256(shl, shr)};
#else
    // Scalar
    u32x8 r;
    for (size_t i = 0; i < 8; ++i) {
        uint32_t v = x.data[i];
        r.data[i] = (v << 7) | (v >> 25);
    }
    return r;
#endif
}

/// SIMD rotl32 by variable amount (per-lane)
inline u32x8 simd_rotl(u32x8 x, u32x8 n) {
#if defined(CSYNCMER_SIMD_NEON)
    // NEON variable shift is tricky, use scalar for now
    auto x_arr = x.to_array();
    auto n_arr = n.to_array();
    std::array<uint32_t, 8> result;
    for (size_t i = 0; i < 8; ++i) {
        uint32_t shift = n_arr[i] & 31;
        result[i] = (x_arr[i] << shift) | (x_arr[i] >> (32 - shift));
    }
    return u32x8::from_array(result);
#elif defined(CSYNCMER_SIMD_AVX2)
    // AVX2: use variable shift (AVX2 has _mm256_sllv_epi32)
    __m256i mask = _mm256_set1_epi32(31);
    __m256i n_masked = _mm256_and_si256(n.v, mask);
    __m256i inv_n = _mm256_sub_epi32(_mm256_set1_epi32(32), n_masked);
    __m256i shl = _mm256_sllv_epi32(x.v, n_masked);
    __m256i shr = _mm256_srlv_epi32(x.v, inv_n);
    return u32x8{_mm256_or_si256(shl, shr)};
#else
    auto x_arr = x.to_array();
    auto n_arr = n.to_array();
    std::array<uint32_t, 8> result;
    for (size_t i = 0; i < 8; ++i) {
        uint32_t shift = n_arr[i] & 31;
        result[i] = (x_arr[i] << shift) | (x_arr[i] >> (32 - shift));
    }
    return u32x8::from_array(result);
#endif
}

/// ASCII to 2-bit encoding: A->0, C->1, T->2, G->3
/// Uses (c >> 1) & 3 for ACGT, treats N and unknown chars as A (0)
/// Matches encode_sequence_2bit behavior for consistency
inline uint8_t ascii_to_2bit(uint8_t c) {
    // Handle N and n as A (0) to match encode_sequence_2bit
    // Also handles lowercase
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'T': case 't': case 'U': case 'u': return 2;
        case 'G': case 'g': return 3;
        default: return 0;  // N, n, and other chars -> A
    }
}

/// Convert array of ASCII bases to 2-bit indices
inline std::array<uint8_t, 8> ascii_to_2bit_array(const std::array<uint8_t, 8>& ascii) {
    std::array<uint8_t, 8> result;
    for (size_t i = 0; i < 8; ++i) {
        result[i] = ascii_to_2bit(ascii[i]);
    }
    return result;
}

/// Lookup table as SIMD vector (4 values -> 8 lanes based on base indices)
/// Optimized with SIMD gather where possible
/// NOTE: indices should be 2-bit encoded (0-3), not ASCII!
inline u32x8 simd_lookup4(const std::array<uint32_t, 4>& table,
                          const std::array<uint8_t, 8>& indices) {
#if defined(CSYNCMER_SIMD_AVX2)
    // AVX2: use gather instruction
    __m256i idx = _mm256_set_epi32(
        indices[7] & 3, indices[6] & 3, indices[5] & 3, indices[4] & 3,
        indices[3] & 3, indices[2] & 3, indices[1] & 3, indices[0] & 3
    );
    __m256i result = _mm256_i32gather_epi32(
        reinterpret_cast<const int*>(table.data()), idx, 4
    );
    return u32x8{result};
#else
    // Scalar/NEON fallback - unroll for better pipelining
    return u32x8::from_array({
        table[indices[0] & 0x3],
        table[indices[1] & 0x3],
        table[indices[2] & 0x3],
        table[indices[3] & 0x3],
        table[indices[4] & 0x3],
        table[indices[5] & 0x3],
        table[indices[6] & 0x3],
        table[indices[7] & 0x3]
    });
#endif
}

/// Gather from a 256-entry table using 8 byte indices (e.g., ASCII lookup)
/// Uses AVX2 gather for efficient parallel lookup
template<size_t N>
inline u32x8 simd_gather(const std::array<uint32_t, N>& table,
                         const std::array<uint8_t, 8>& indices) {
#if defined(CSYNCMER_SIMD_AVX2)
    // AVX2: gather with 32-bit indices (extended from uint8_t)
    __m256i idx = _mm256_set_epi32(
        indices[7], indices[6], indices[5], indices[4],
        indices[3], indices[2], indices[1], indices[0]
    );
    __m256i result = _mm256_i32gather_epi32(
        reinterpret_cast<const int*>(table.data()), idx, 4
    );
    return u32x8{result};
#else
    // Scalar fallback
    return u32x8::from_array({
        table[indices[0]], table[indices[1]], table[indices[2]], table[indices[3]],
        table[indices[4]], table[indices[5]], table[indices[6]], table[indices[7]]
    });
#endif
}

/// Gather from a raw pointer table using 8 byte indices
inline u32x8 simd_gather_ptr(const uint32_t* table,
                              const std::array<uint8_t, 8>& indices) {
#if defined(CSYNCMER_SIMD_AVX2)
    __m256i idx = _mm256_set_epi32(
        indices[7], indices[6], indices[5], indices[4],
        indices[3], indices[2], indices[1], indices[0]
    );
    __m256i result = _mm256_i32gather_epi32(
        reinterpret_cast<const int*>(table), idx, 4
    );
    return u32x8{result};
#else
    return u32x8::from_array({
        table[indices[0]], table[indices[1]], table[indices[2]], table[indices[3]],
        table[indices[4]], table[indices[5]], table[indices[6]], table[indices[7]]
    });
#endif
}

/// 8 parallel ntHash instances for SIMD processing
/// Each lane processes a different chunk of the sequence
class NtHashSimd {
public:
    NtHashSimd() = default;

    explicit NtHashSimd(size_t k) : k_(k) {
        // Pre-compute rotated tables
        f_rot_ = hash::detail::make_f_rot(k);

        // Compute initial hash value (hashing k-1 zeros)
        fw_init_ = hash::detail::make_fw_init(k);

        // Initialize all 8 lanes with the same initial value
        fw_ = u32x8::splat(fw_init_);
    }

    /// Initialize all 8 hashers at their starting positions
    /// Returns the first 8 hashes (one per chunk)
    ///
    /// CORRECT APPROACH: Compute hash directly, then set up state for rolling
    template<typename SeqT>
    u32x8 init_parallel(const SeqT& seq, const std::array<size_t, 8>& starts) {
        // Compute hash directly: H = rotl^(k-1)(F[s[0]]) ^ rotl^(k-2)(F[s[1]]) ^ ... ^ F[s[k-1]]
        fw_ = U32X8_ZERO();
        for (size_t i = 0; i < k_; ++i) {
            std::array<uint8_t, 8> bases;
            for (size_t lane = 0; lane < 8; ++lane) {
                size_t pos = starts[lane] + i;
                bases[lane] = (pos < seq.size()) ? ascii_to_2bit(seq[pos]) : 0;
            }
            u32x8 f_base = simd_lookup4(hash::detail::NTHASH_F, bases);
            fw_ = simd_rotl7(fw_) ^ f_base;
        }
        u32x8 hash = fw_;

        // Set up state for rolling: fw_ = hash ^ f_rot[first_base]
        std::array<uint8_t, 8> first_bases;
        for (size_t lane = 0; lane < 8; ++lane) {
            first_bases[lane] = (starts[lane] < seq.size()) ? ascii_to_2bit(seq[starts[lane]]) : 0;
        }
        fw_ = hash ^ simd_lookup4(f_rot_, first_bases);

        return hash;
    }

    /// Roll: add in_bases, state uses new_first_bases for next roll
    /// After rolling from position i to i+1:
    ///   - out_bases = seq[i] (being removed, not actually used in computation)
    ///   - in_bases = seq[i + k] (being added)
    ///   - new_first_bases = seq[i + 1] (first base of new s-mer, needed for state)
    /// Accepts ASCII bases (A,C,G,T) and converts internally
    u32x8 roll(const std::array<uint8_t, 8>& out_bases,
               const std::array<uint8_t, 8>& in_bases,
               const std::array<uint8_t, 8>& new_first_bases) {
        (void)out_bases;  // Not needed - already encoded in fw_ state
        auto in_idx = ascii_to_2bit_array(in_bases);
        auto new_first_idx = ascii_to_2bit_array(new_first_bases);
        u32x8 f_in = simd_lookup4(hash::detail::NTHASH_F, in_idx);
        u32x8 f_rot_new_first = simd_lookup4(f_rot_, new_first_idx);
        u32x8 fw_out = simd_rotl7(fw_) ^ f_in;
        fw_ = fw_out ^ f_rot_new_first;
        return fw_out;
    }

    /// Convenience: roll with bases from sequence positions
    /// out_positions: positions of bases being removed
    /// in_positions: positions of bases being added
    /// new_first_positions: positions of first bases in new s-mers (= out_positions + 1)
    template<typename SeqT>
    u32x8 roll(const SeqT& seq,
               const std::array<size_t, 8>& out_positions,
               const std::array<size_t, 8>& in_positions,
               const std::array<size_t, 8>& new_first_positions) {
        std::array<uint8_t, 8> out_bases;
        std::array<uint8_t, 8> in_bases;
        std::array<uint8_t, 8> new_first_bases;
        for (size_t lane = 0; lane < 8; ++lane) {
            out_bases[lane] = (out_positions[lane] < seq.size()) ? seq[out_positions[lane]] : 0;
            in_bases[lane] = (in_positions[lane] < seq.size()) ? seq[in_positions[lane]] : 0;
            new_first_bases[lane] = (new_first_positions[lane] < seq.size()) ? seq[new_first_positions[lane]] : 0;
        }
        return roll(out_bases, in_bases, new_first_bases);
    }

    [[nodiscard]] u32x8 hash() const { return fw_; }

private:
    size_t k_ = 0;
    std::array<uint32_t, 4> f_rot_{};
    uint32_t fw_init_ = 0;
    u32x8 fw_;
};

/// State for processing 8 parallel chunks
struct ParallelChunkState {
    std::array<size_t, 8> starts;   // Starting position of each chunk
    std::array<size_t, 8> ends;     // Ending position of each chunk
    std::array<size_t, 8> current;  // Current position in each chunk
    size_t chunk_len;               // Length of each chunk

    /// Initialize for a sequence with num_kmers k-mers
    void init(size_t num_kmers, size_t w) {
        size_t num_windows = (num_kmers > w - 1) ? num_kmers - w + 1 : 0;
        chunk_len = (num_windows + 7) / 8;

        for (size_t lane = 0; lane < 8; ++lane) {
            size_t start_window = lane * chunk_len;
            if (start_window >= num_windows) {
                // Lane is inactive
                starts[lane] = num_kmers;
                ends[lane] = num_kmers;
                current[lane] = num_kmers;
            } else {
                starts[lane] = start_window;
                ends[lane] = std::min((lane + 1) * chunk_len, num_windows) + w - 1;
                current[lane] = start_window;
            }
        }
    }

    /// Get bases at current position from each chunk
    template<typename SeqT>
    std::array<uint8_t, 8> get_bases(const SeqT& seq, size_t offset) const {
        std::array<uint8_t, 8> bases;
        for (size_t lane = 0; lane < 8; ++lane) {
            size_t pos = starts[lane] + offset;
            bases[lane] = (pos < seq.size()) ? seq[pos] : 0;
        }
        return bases;
    }

    /// Check if any lane is still active
    bool any_active(size_t step) const {
        for (size_t lane = 0; lane < 8; ++lane) {
            if (starts[lane] + step < ends[lane]) {
                return true;
            }
        }
        return false;
    }

    /// Get validity mask for step
    uint8_t get_validity_mask(size_t step, size_t num_kmers) const {
        uint8_t mask = 0;
        for (size_t lane = 0; lane < 8; ++lane) {
            size_t pos = starts[lane] + step;
            if (pos < ends[lane] && pos < num_kmers) {
                mask |= (1u << lane);
            }
        }
        return mask;
    }
};

} // namespace csyncmer_simd
