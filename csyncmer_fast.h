#ifndef CSYNCMER_FAST_H
#define CSYNCMER_FAST_H

// Fast closed syncmer detection using ntHash32
// Legacy infrastructure moved to misc/older/legacy_infrastructure.h

#include <stddef.h>
#include <stdint.h>
#include <cstdio>

#ifdef __cplusplus
#include <vector>
#include <array>
#include <algorithm>
#include <cstdint>
#include "misc/simd/vec.hpp"
#include "misc/simd/nthash.hpp"
#include "misc/simd/hash_simd.hpp"
#include "misc/simd/sliding_min.hpp"

namespace csyncmer_fast_simd {

// Redundant ntHash32 variants moved to misc/older/syncmer_nthash32_variants.h
// Keeping only: _fused_rescan_branchless (scalar fallback) and _simd_multiwindow (fastest)

/// Direct ASCII lookup table (256 entries, most unused)
/// Avoids separate encoding step - used by syncmer functions
inline constexpr std::array<uint32_t, 256> make_ascii_hash_table_early() {
    std::array<uint32_t, 256> table{};
    // A=0, C=1, T=2, G=3
    table['A'] = table['a'] = csyncmer_simd::hash::detail::NTHASH_F[0];
    table['C'] = table['c'] = csyncmer_simd::hash::detail::NTHASH_F[1];
    table['T'] = table['t'] = csyncmer_simd::hash::detail::NTHASH_F[2];
    table['G'] = table['g'] = csyncmer_simd::hash::detail::NTHASH_F[3];
    return table;
}

inline constexpr std::array<uint8_t, 256> make_ascii_to_idx_early() {
    std::array<uint8_t, 256> table{};
    table['A'] = table['a'] = 0;
    table['C'] = table['c'] = 1;
    table['T'] = table['t'] = 2;
    table['G'] = table['g'] = 3;
    return table;
}

inline size_t compute_closed_syncmers_nthash32_fused_rescan_branchless(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    size_t* num_syncmers
) {
    if (length < K) {
        *num_syncmers = 0;
        return 0;
    }

    // Direct ASCII tables
    static const auto F_ASCII = make_ascii_hash_table_early();
    static const auto IDX_ASCII = make_ascii_to_idx_early();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    // Power-of-2 circular buffer
    size_t buf_size = 1;
    while (buf_size < window_size) buf_size <<= 1;
    size_t buf_mask = buf_size - 1;
    std::vector<uint32_t> hash_buffer(buf_size);

    size_t syncmer_count = 0;

    // First hash - compute directly (correct approach)
    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_simd::hash::detail::rotl7(hash) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];  // State uses first base of s-mer
    hash_buffer[0] = hash;

    // Fill initial window
    for (size_t i = 1; i < window_size; ++i) {
        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];  // FIXED: use new first base
        hash_buffer[i & buf_mask] = hash;
    }

    // Find minimum in first window
    uint32_t min_hash = UINT32_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < window_size; ++i) {
        if (hash_buffer[i] < min_hash) {
            min_hash = hash_buffer[i];
            min_pos = i;
        }
    }

    // Check first k-mer
    if (min_pos == 0 || min_pos == window_size - 1) {
        syncmer_count++;
    }

    // Main loop with branch-free min update
    for (size_t kmer_idx = 1; kmer_idx < num_smers - window_size + 1; ++kmer_idx) {
        size_t i = kmer_idx + window_size - 1;

        // Compute new hash
        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];  // FIXED: use new first base
        hash_buffer[i & buf_mask] = hash;

        // RESCAN: update minimum
        if (min_pos < kmer_idx) {
            // Minimum fell out - rescan (same as before)
            min_hash = UINT32_MAX;
            for (size_t j = kmer_idx; j <= i; ++j) {
                uint32_t h = hash_buffer[j & buf_mask];
                if (h < min_hash) {
                    min_hash = h;
                    min_pos = j;
                }
            }
        } else {
            // Branch-free update using conditional moves
            // is_smaller is 1 if hash < min_hash, 0 otherwise
            uint32_t is_smaller = (hash < min_hash) ? 1 : 0;
            // Use ternary which compiles to cmov
            min_hash = is_smaller ? hash : min_hash;
            min_pos = is_smaller ? i : min_pos;
        }

        // Check closed syncmer condition
        size_t min_offset = min_pos - kmer_idx;
        if (min_offset == 0 || min_offset == window_size - 1) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_FUSED_RESCAN_BF]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}
#if defined(__AVX2__)
/// SIMD Multi-Window Syncmer Detection
/// Processes 8 consecutive overlapping windows in parallel using AVX2
/// Uses SIMD rescan for O(w) work per 8 k-mers with branch-free loads
inline size_t compute_closed_syncmers_nthash32_simd_multiwindow(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    size_t* num_syncmers
) {
    using namespace csyncmer_simd;

    if (length < K) {
        *num_syncmers = 0;
        return 0;
    }

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = num_smers - window_size + 1;

    if (num_kmers < 8) {
        return compute_closed_syncmers_nthash32_fused_rescan_branchless(
            sequence, length, K, S, num_syncmers);
    }

    // Ring buffer with extended area to eliminate wrap-around branches
    // Main buffer: 64 elements (power of 2)
    // Extended area: 32 elements mirroring positions 0..31
    // This allows contiguous SIMD loads even when buffer wraps
    constexpr size_t BUF_SIZE = 64;
    constexpr size_t BUF_MASK = BUF_SIZE - 1;
    alignas(32) uint32_t buf[BUF_SIZE + 32];

    size_t syncmer_count = 0;
    u32x8 window_size_minus_1 = u32x8::splat(static_cast<uint32_t>(window_size - 1));

    // Direct ASCII lookup for speed
    static const auto F_ASCII = make_ascii_hash_table_early();
    static const auto IDX_ASCII = make_ascii_to_idx_early();
    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);
    const uint8_t* raw_seq = reinterpret_cast<const uint8_t*>(sequence);

    // Initialize first hash directly (not using NtHash class overhead)
    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_simd::hash::detail::rotl7(hash) ^ F_ASCII[raw_seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[raw_seq[0]]];
    buf[0] = hash;
    buf[BUF_SIZE] = hash;  // Mirror

    // Fill initial window + first batch
    size_t init_size = window_size + 7;
    for (size_t i = 1; i < init_size; ++i) {
        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[raw_seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[raw_seq[i]]];
        size_t buf_idx = i & BUF_MASK;
        buf[buf_idx] = hash;
        if (buf_idx < 32) buf[BUF_SIZE + buf_idx] = hash;
    }
    size_t hash_count = init_size;

    // Process 8 k-mers at a time
    size_t k = 0;
    for (; k + 8 <= num_kmers; k += 8) {
        // SIMD minimum finding
        size_t start_idx = k & BUF_MASK;
        u32x8 min_hash = U32X8_MAX();
        u32x8 min_pos = U32X8_ZERO();

        for (size_t j = 0; j < window_size; ++j) {
            u32x8 cur_hash = u32x8::load(&buf[start_idx + j]);
            u32x8 cur_pos = u32x8::splat(static_cast<uint32_t>(j));
            u32x8 is_smaller = min_hash.cmp_gt(cur_hash);
            min_hash = u32x8::blend(is_smaller, cur_hash, min_hash);
            min_pos = u32x8::blend(is_smaller, cur_pos, min_pos);
        }

        // Check closed syncmer condition
        u32x8 is_first = min_pos.cmp_eq(U32X8_ZERO());
        u32x8 is_last = min_pos.cmp_eq(window_size_minus_1);
        u32x8 is_syncmer = is_first | is_last;
        syncmer_count += __builtin_popcount(is_syncmer.movemask());

        // Compute next 8 hashes
        size_t end_hash = k + 8 + 7 + window_size;
        if (end_hash > num_smers) end_hash = num_smers;
        for (; hash_count < end_hash; ++hash_count) {
            hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[raw_seq[hash_count]]];
            size_t buf_idx = hash_count & BUF_MASK;
            buf[buf_idx] = hash;
            if (buf_idx < 32) buf[BUF_SIZE + buf_idx] = hash;
        }
    }

    // Handle remaining k-mers with scalar rescan using buffer
    for (; k < num_kmers; ++k) {
        // Ensure we have hashes for this k-mer's window
        size_t end_hash = k + window_size;
        for (; hash_count < end_hash; ++hash_count) {
            hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[raw_seq[hash_count]]];
            size_t buf_idx = hash_count & BUF_MASK;
            buf[buf_idx] = hash;
            if (buf_idx < 32) buf[BUF_SIZE + buf_idx] = hash;
        }

        uint32_t min_hash = UINT32_MAX;
        uint32_t min_pos = 0;
        for (size_t j = 0; j < window_size; ++j) {
            uint32_t h = buf[(k + j) & BUF_MASK];
            if (h < min_hash) { min_hash = h; min_pos = static_cast<uint32_t>(j); }
        }

        if (min_pos == 0 || min_pos == static_cast<uint32_t>(window_size - 1)) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_SIMD_MULTIWINDOW]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}
#endif  // __AVX2__
} // namespace csyncmer_fast_simd

#endif // __cplusplus

#endif // CSYNCMER_FAST_H
