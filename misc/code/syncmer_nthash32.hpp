#ifndef SYNCMER_NTHASH32_H
#define SYNCMER_NTHASH32_H

// ntHash32-based syncmer implementations (deprecated/legacy)
// Superseded by _fused_rescan_branchless and TWOSTACK in csyncmer_fast.h
//
// This file consolidates:
// - Various ntHash32 syncmer algorithms (rescan, deque, vanherk, etc.)
// - Hash-only benchmarks
// - SIMD multi-window implementation (from deprecated_simd_mw.h)

#include "legacy_infrastructure.hpp"
#include "simd/vec.hpp"
#include "simd/nthash.hpp"
#include "simd/hash_simd.hpp"

namespace csyncmer_nthash32 {

using namespace csyncmer_simd;

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/// Convert ASCII sequence to 2-bit encoding for SIMD processing
/// A=0, C=1, T=2, G=3 (matches ntHash table ordering)
inline std::vector<uint8_t> encode_sequence_2bit(const char* seq, size_t len) {
    std::vector<uint8_t> encoded(len);
    for (size_t i = 0; i < len; ++i) {
        switch (seq[i]) {
            case 'A': case 'a': encoded[i] = 0; break;
            case 'C': case 'c': encoded[i] = 1; break;
            case 'T': case 't': encoded[i] = 2; break;
            case 'G': case 'g': encoded[i] = 3; break;
            default: encoded[i] = 0; break;
        }
    }
    return encoded;
}

/// Direct ASCII lookup table (256 entries, most unused)
inline constexpr std::array<uint32_t, 256> make_ascii_hash_table_local() {
    std::array<uint32_t, 256> table{};
    table['A'] = table['a'] = csyncmer_simd::hash::detail::NTHASH_F[0];
    table['C'] = table['c'] = csyncmer_simd::hash::detail::NTHASH_F[1];
    table['T'] = table['t'] = csyncmer_simd::hash::detail::NTHASH_F[2];
    table['G'] = table['g'] = csyncmer_simd::hash::detail::NTHASH_F[3];
    return table;
}

inline constexpr std::array<uint8_t, 256> make_ascii_to_idx_local() {
    std::array<uint8_t, 256> table{};
    table['A'] = table['a'] = 0;
    table['C'] = table['c'] = 1;
    table['T'] = table['t'] = 2;
    table['G'] = table['g'] = 3;
    return table;
}

// ============================================================================
// SYNCMER IMPLEMENTATIONS
// ============================================================================

/// Naive O(N*W) implementation - compute each hash from scratch, scan each window
inline size_t csyncmer_nthash32_naive(
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

    static const auto F_ASCII = make_ascii_hash_table_local();

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = length - K + 1;

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    // Compute all s-mer hashes from scratch (no rolling hash)
    std::vector<uint32_t> smer_hashes(num_smers);
    for (size_t i = 0; i < num_smers; ++i) {
        uint32_t h = 0;
        for (size_t j = 0; j < S; ++j) {
            h = csyncmer_simd::hash::detail::rotl7(h) ^ F_ASCII[seq[i + j]];
        }
        smer_hashes[i] = h;
    }

    // O(N*W) scan to find syncmers
    size_t syncmer_count = 0;
    for (size_t i = 0; i < num_kmers; ++i) {
        uint32_t min_hash = UINT32_MAX;
        size_t min_pos = 0;
        for (size_t j = 0; j < window_size; ++j) {
            if (smer_hashes[i + j] < min_hash) {
                min_hash = smer_hashes[i + j];
                min_pos = j;
            }
        }
        if (min_pos == 0 || min_pos == window_size - 1) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_NAIVE]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    printf("[NTHASH32_NAIVE]:: HASHED %lu S-MERS\n", num_smers);
    return syncmer_count;
}

/// 2-bit encoding + precompute all hashes + RESCAN
inline size_t csyncmer_nthash32_2bit_rescan(
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

    auto encoded = encode_sequence_2bit(sequence, length);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = length - K + 1;

    std::vector<uint32_t> smer_hashes(num_smers);
    csyncmer_simd::hash::NtHash hasher(S);
    smer_hashes[0] = hasher.init(encoded, 0);
    for (size_t i = 1; i < num_smers; ++i) {
        smer_hashes[i] = hasher.roll(encoded[i - 1], encoded[i + S - 1], encoded[i]);
    }

    size_t syncmer_count = 0;
    size_t rescan_count = 0;

    uint32_t min_hash = UINT32_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < window_size; ++i) {
        if (smer_hashes[i] < min_hash) {
            min_hash = smer_hashes[i];
            min_pos = i;
        }
    }

    if (min_pos == 0 || min_pos == window_size - 1) {
        syncmer_count++;
    }

    for (size_t kmer_idx = 1; kmer_idx < num_kmers; ++kmer_idx) {
        size_t window_start = kmer_idx;
        size_t window_end = kmer_idx + window_size;
        size_t new_pos = window_end - 1;

        if (min_pos < window_start) {
            min_hash = UINT32_MAX;
            for (size_t i = window_start; i < window_end; ++i) {
                if (smer_hashes[i] < min_hash) {
                    min_hash = smer_hashes[i];
                    min_pos = i;
                }
            }
            rescan_count++;
        } else {
            if (smer_hashes[new_pos] < min_hash) {
                min_hash = smer_hashes[new_pos];
                min_pos = new_pos;
            }
        }

        size_t min_offset = min_pos - window_start;
        if (min_offset == 0 || min_offset == window_size - 1) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_2BIT_RESCAN]:: COMPUTED %lu CLOSED SYNCMERS (rescans: %lu)\n",
           syncmer_count, rescan_count);
    return syncmer_count;
}

/// 2-bit encoding + precompute all hashes + DEQUE
inline size_t csyncmer_nthash32_2bit_deque(
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

    auto encoded = encode_sequence_2bit(sequence, length);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    std::vector<uint32_t> smer_hashes(num_smers);
    csyncmer_simd::hash::NtHash hasher(S);
    smer_hashes[0] = hasher.init(encoded, 0);
    for (size_t i = 1; i < num_smers; ++i) {
        smer_hashes[i] = hasher.roll(encoded[i - 1], encoded[i + S - 1], encoded[i]);
    }

    std::vector<size_t> deque(num_smers);
    size_t front = 0, back = 0;
    size_t syncmer_count = 0;

    for (size_t i = 0; i < num_smers; ++i) {
        while (back > front && smer_hashes[deque[back - 1]] > smer_hashes[i]) {
            back--;
        }
        deque[back++] = i;

        if (i >= window_size && deque[front] <= i - window_size) {
            front++;
        }

        if (i >= window_size - 1) {
            size_t min_pos = deque[front];
            size_t kmer_start = i - window_size + 1;
            size_t min_offset = min_pos - kmer_start;

            if (min_offset == 0 || min_offset == window_size - 1) {
                syncmer_count++;
            }
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_2BIT_DEQUE]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}

/// Direct ASCII + precompute all hashes + RESCAN
inline size_t csyncmer_nthash32_direct_rescan(
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

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = length - K + 1;

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    std::vector<uint32_t> smer_hashes(num_smers);

    uint32_t h = 0;
    for (size_t j = 0; j < S; ++j) {
        h = csyncmer_simd::hash::detail::rotl7(h) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = h ^ f_rot[IDX_ASCII[seq[0]]];
    smer_hashes[0] = h;

    for (size_t i = 1; i < num_smers; ++i) {
        uint32_t fw_out = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = fw_out ^ f_rot[IDX_ASCII[seq[i]]];
        smer_hashes[i] = fw_out;
    }

    size_t syncmer_count = 0;
    size_t rescan_count = 0;

    uint32_t min_hash = UINT32_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < window_size; ++i) {
        if (smer_hashes[i] < min_hash) {
            min_hash = smer_hashes[i];
            min_pos = i;
        }
    }

    if (min_pos == 0 || min_pos == window_size - 1) {
        syncmer_count++;
    }

    for (size_t kmer_idx = 1; kmer_idx < num_kmers; ++kmer_idx) {
        size_t window_start = kmer_idx;
        size_t window_end = kmer_idx + window_size;
        size_t new_pos = window_end - 1;

        if (min_pos < window_start) {
            min_hash = UINT32_MAX;
            for (size_t i = window_start; i < window_end; ++i) {
                if (smer_hashes[i] < min_hash) {
                    min_hash = smer_hashes[i];
                    min_pos = i;
                }
            }
            rescan_count++;
        } else {
            if (smer_hashes[new_pos] < min_hash) {
                min_hash = smer_hashes[new_pos];
                min_pos = new_pos;
            }
        }

        size_t min_offset = min_pos - window_start;
        if (min_offset == 0 || min_offset == window_size - 1) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_DIRECT_RESCAN]:: COMPUTED %lu CLOSED SYNCMERS (rescans: %lu)\n",
           syncmer_count, rescan_count);
    return syncmer_count;
}

/// Fused hash + deque (single pass, no large buffer)
inline size_t csyncmer_nthash32_fused_deque(
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

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    size_t buf_size = 1;
    while (buf_size < window_size) buf_size <<= 1;
    size_t buf_mask = buf_size - 1;
    std::vector<uint32_t> hash_buffer(buf_size);

    std::vector<size_t> deque(buf_size);
    size_t front = 0, back = 0;

    size_t syncmer_count = 0;

    uint32_t first_hash = 0;
    for (size_t j = 0; j < S; ++j) {
        first_hash = csyncmer_simd::hash::detail::rotl7(first_hash) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = first_hash ^ f_rot[IDX_ASCII[seq[0]]];

    for (size_t i = 0; i < num_smers; ++i) {
        uint32_t hash;
        if (i == 0) {
            hash = first_hash;
        } else {
            hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        }

        hash_buffer[i & buf_mask] = hash;

        while (back > front && hash_buffer[deque[(back - 1) & buf_mask] & buf_mask] > hash) {
            back--;
        }
        deque[back++ & buf_mask] = i;

        while (front < back && deque[front & buf_mask] + window_size <= i) {
            front++;
        }

        if (i >= window_size - 1) {
            size_t min_pos = deque[front & buf_mask];
            size_t kmer_start = i - window_size + 1;
            size_t min_offset = min_pos - kmer_start;

            if (min_offset == 0 || min_offset == window_size - 1) {
                syncmer_count++;
            }
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_FUSED_DEQUE]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}

/// Fused hash + RESCAN (single pass, circular buffer)
inline size_t csyncmer_nthash32_fused_rescan(
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

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    size_t buf_size = 1;
    while (buf_size < window_size) buf_size <<= 1;
    size_t buf_mask = buf_size - 1;
    std::vector<uint32_t> hash_buffer(buf_size);

    size_t syncmer_count = 0;

    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_simd::hash::detail::rotl7(hash) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];
    hash_buffer[0] = hash;

    for (size_t i = 1; i < window_size; ++i) {
        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;
    }

    uint32_t min_hash = UINT32_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < window_size; ++i) {
        if (hash_buffer[i] < min_hash) {
            min_hash = hash_buffer[i];
            min_pos = i;
        }
    }

    if (min_pos == 0 || min_pos == window_size - 1) {
        syncmer_count++;
    }

    for (size_t kmer_idx = 1; kmer_idx < num_smers - window_size + 1; ++kmer_idx) {
        size_t i = kmer_idx + window_size - 1;

        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;

        if (min_pos < kmer_idx) {
            min_hash = UINT32_MAX;
            for (size_t j = kmer_idx; j <= i; ++j) {
                uint32_t h = hash_buffer[j & buf_mask];
                if (h < min_hash) {
                    min_hash = h;
                    min_pos = j;
                }
            }
        } else {
            if (hash < min_hash) {
                min_hash = hash;
                min_pos = i;
            }
        }

        size_t min_offset = min_pos - kmer_idx;
        if (min_offset == 0 || min_offset == window_size - 1) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_FUSED_RESCAN]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}

/// Fused hash + two-stack sliding minimum (block-based Van Herk/Gil-Werman)
/// O(1) worst-case per element: ~3 comparisons (1 prefix + 1 window + 1 amortized suffix)
/// Eliminates the variable-length rescan loop of fused_rescan.
inline size_t csyncmer_nthash32_fused_twostack(
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

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    if (window_size > 64) {
        return csyncmer_nthash32_fused_rescan(sequence, length, K, S, num_syncmers);
    }

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    // Block size = window_size. Pack (hash << 32) | position for combined comparison.
    uint64_t ring[64];
    uint64_t suffix_min[64];
    uint64_t prefix_min = UINT64_MAX;

    size_t syncmer_count = 0;

    // Initialize first s-mer hash
    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_simd::hash::detail::rotl7(hash) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];

    // Fill initial block (s-mers 0 to window_size-1)
    ring[0] = (static_cast<uint64_t>(hash) << 32) | 0u;
    prefix_min = ring[0];

    for (size_t i = 1; i < window_size; ++i) {
        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        uint64_t packed = (static_cast<uint64_t>(hash) << 32) | static_cast<uint32_t>(i);
        ring[i] = packed;
        if (packed < prefix_min) prefix_min = packed;
    }

    // Check first window (kmer 0): entirely in first block
    {
        size_t min_pos = static_cast<size_t>(prefix_min & 0xFFFFFFFF);
        if (min_pos == 0 || min_pos == window_size - 1) {
            syncmer_count++;
        }
    }

    // Block-unrolled main loop: process W elements per outer iteration.
    // Eliminates per-element branches for block boundary and start_offset wrap.
    for (size_t block_start = window_size; block_start < num_smers; block_start += window_size) {
        // Suffix-min scan (once per block, amortized O(1) per element)
        suffix_min[window_size - 1] = ring[window_size - 1];
        for (size_t j = window_size - 1; j > 0; --j) {
            suffix_min[j - 1] = std::min(ring[j - 1], suffix_min[j]);
        }
        prefix_min = UINT64_MAX;

        size_t block_end = std::min(block_start + window_size, num_smers);
        size_t block_len = block_end - block_start;
        size_t inner_count = std::min(block_len, window_size - 1);

        // Inner loop: window spans previous block (suffix) + current block (prefix)
        // No per-element branch for block boundary or start_offset wrap.
        for (size_t j = 0; j < inner_count; ++j) {
            size_t i = block_start + j;

            hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[seq[i]]];

            uint64_t packed = (static_cast<uint64_t>(hash) << 32) | static_cast<uint32_t>(i);
            ring[j] = packed;
            prefix_min = std::min(prefix_min, packed);

            uint64_t window_min = std::min(suffix_min[j + 1], prefix_min);

            size_t kmer_idx = i - window_size + 1;
            size_t min_pos = static_cast<size_t>(window_min & 0xFFFFFFFF);
            size_t min_offset = min_pos - kmer_idx;
            syncmer_count += (min_offset == 0) | (min_offset == window_size - 1);
        }

        // Last element of full block: window entirely in current block (no suffix needed)
        if (block_len == window_size) {
            size_t i = block_start + window_size - 1;

            hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[seq[i]]];

            uint64_t packed = (static_cast<uint64_t>(hash) << 32) | static_cast<uint32_t>(i);
            ring[window_size - 1] = packed;
            prefix_min = std::min(prefix_min, packed);

            uint64_t window_min = prefix_min;

            size_t kmer_idx = i - window_size + 1;
            size_t min_pos = static_cast<size_t>(window_min & 0xFFFFFFFF);
            size_t min_offset = min_pos - kmer_idx;
            syncmer_count += (min_offset == 0) | (min_offset == window_size - 1);
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_FUSED_TWOSTACK]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}

/// Van Herk/Gil-Werman O(1) block algorithm
inline size_t csyncmer_nthash32_vanherk(
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

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = num_smers - window_size + 1;

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    std::vector<uint32_t> hashes(num_smers);

    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_simd::hash::detail::rotl7(hash) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];
    hashes[0] = hash;

    for (size_t i = 1; i < num_smers; ++i) {
        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hashes[i] = hash;
    }

    size_t num_blocks = (num_smers + window_size - 1) / window_size;

    std::vector<uint64_t> suffix_min(num_smers);
    std::vector<uint64_t> prefix_min(num_smers);

    auto pack = [](uint32_t h, size_t pos) -> uint64_t {
        return (static_cast<uint64_t>(h) << 32) | static_cast<uint32_t>(pos);
    };
    auto unpack_pos = [](uint64_t packed) -> size_t {
        return static_cast<size_t>(packed & 0xFFFFFFFF);
    };

    for (size_t block = 0; block < num_blocks; ++block) {
        size_t block_start = block * window_size;
        size_t block_end = std::min(block_start + window_size, num_smers);

        suffix_min[block_end - 1] = pack(hashes[block_end - 1], block_end - 1);
        for (size_t i = block_end - 1; i > block_start; --i) {
            uint64_t curr = pack(hashes[i - 1], i - 1);
            suffix_min[i - 1] = std::min(curr, suffix_min[i]);
        }

        prefix_min[block_start] = pack(hashes[block_start], block_start);
        for (size_t i = block_start + 1; i < block_end; ++i) {
            uint64_t curr = pack(hashes[i], i);
            prefix_min[i] = std::min(prefix_min[i - 1], curr);
        }
    }

    size_t syncmer_count = 0;

    for (size_t kmer_idx = 0; kmer_idx < num_kmers; ++kmer_idx) {
        size_t start = kmer_idx;
        size_t end = kmer_idx + window_size - 1;

        size_t start_block = start / window_size;
        size_t end_block = end / window_size;

        uint64_t window_min;
        if (start_block == end_block) {
            window_min = suffix_min[start];
        } else {
            window_min = std::min(suffix_min[start], prefix_min[end]);
        }

        size_t min_pos = unpack_pos(window_min);
        size_t min_offset = min_pos - kmer_idx;

        if (min_offset == 0 || min_offset == window_size - 1) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_VANHERK]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}

#if defined(__AVX2__)
/// SIMD argmin: find minimum value and position in circular buffer window
inline std::pair<uint32_t, size_t> simd_argmin_window(
    const uint32_t* buffer,
    size_t buf_mask,
    size_t start,
    size_t count
) {
    alignas(32) uint32_t temp[32];

    for (size_t i = 0; i < count; ++i) {
        temp[i] = buffer[(start + i) & buf_mask];
    }
    for (size_t i = count; i < 32; ++i) {
        temp[i] = UINT32_MAX;
    }

    __m256i v0 = _mm256_load_si256((__m256i*)&temp[0]);
    __m256i v1 = _mm256_load_si256((__m256i*)&temp[8]);
    __m256i v2 = _mm256_load_si256((__m256i*)&temp[16]);
    __m256i v3 = _mm256_load_si256((__m256i*)&temp[24]);

    __m256i min01 = _mm256_min_epu32(v0, v1);
    __m256i min23 = _mm256_min_epu32(v2, v3);
    __m256i min_all = _mm256_min_epu32(min01, min23);

    __m128i lo = _mm256_castsi256_si128(min_all);
    __m128i hi = _mm256_extracti128_si256(min_all, 1);
    __m128i m = _mm_min_epu32(lo, hi);
    m = _mm_min_epu32(m, _mm_shuffle_epi32(m, 0x4E));
    m = _mm_min_epu32(m, _mm_shuffle_epi32(m, 0xB1));
    uint32_t min_val = static_cast<uint32_t>(_mm_cvtsi128_si32(m));

    __m256i target = _mm256_set1_epi32(static_cast<int32_t>(min_val));

    __m256i cmp0 = _mm256_cmpeq_epi32(v0, target);
    int mask0 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp0));
    if (mask0) {
        return {min_val, static_cast<size_t>(__builtin_ctz(mask0))};
    }

    __m256i cmp1 = _mm256_cmpeq_epi32(v1, target);
    int mask1 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp1));
    if (mask1) {
        return {min_val, 8 + static_cast<size_t>(__builtin_ctz(mask1))};
    }

    __m256i cmp2 = _mm256_cmpeq_epi32(v2, target);
    int mask2 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp2));
    if (mask2) {
        return {min_val, 16 + static_cast<size_t>(__builtin_ctz(mask2))};
    }

    for (size_t i = 24; i < count; ++i) {
        if (temp[i] == min_val) return {min_val, i};
    }
    return {min_val, 0};
}

/// Fused hash + SIMD RESCAN
inline size_t csyncmer_nthash32_simd_rescan(
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

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    size_t buf_size = 1;
    while (buf_size < window_size) buf_size <<= 1;
    size_t buf_mask = buf_size - 1;
    std::vector<uint32_t> hash_buffer(buf_size);

    size_t syncmer_count = 0;

    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_simd::hash::detail::rotl7(hash) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];
    hash_buffer[0] = hash;

    for (size_t i = 1; i < window_size; ++i) {
        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;
    }

    auto [min_hash, min_offset] = simd_argmin_window(hash_buffer.data(), buf_mask, 0, window_size);
    size_t min_pos = min_offset;

    if (min_pos == 0 || min_pos == window_size - 1) {
        syncmer_count++;
    }

    for (size_t kmer_idx = 1; kmer_idx < num_smers - window_size + 1; ++kmer_idx) {
        size_t i = kmer_idx + window_size - 1;

        hash = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;

        if (min_pos < kmer_idx) {
            auto [new_min, offset] = simd_argmin_window(
                hash_buffer.data(), buf_mask, kmer_idx, window_size
            );
            min_hash = new_min;
            min_pos = kmer_idx + offset;
        } else {
            if (hash < min_hash) {
                min_hash = hash;
                min_pos = i;
            }
        }

        size_t min_off = min_pos - kmer_idx;
        if (min_off == 0 || min_off == window_size - 1) {
            syncmer_count++;
        }
    }

    *num_syncmers = syncmer_count;
    printf("[NTHASH32_SIMD_RESCAN]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    return syncmer_count;
}

// ============================================================================
// SIMD MULTI-WINDOW IMPLEMENTATION (DEPRECATED)
// Performance: ~192 MB/s (slower than TWOSTACK @ 551 MB/s)
// ============================================================================

// Note: This function requires csyncmer_fast.h to be included first
// for csyncmer_rotl7, csyncmer_init_ascii_hash_table, etc.

static inline size_t csyncmer_nthash32_simd_multiwindow(
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

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = num_smers - window_size + 1;

    // Fall back to scalar for small inputs
    if (num_kmers < 8) {
        return csyncmer_rescan_32_count(
            sequence, length, K, S, num_syncmers);
    }

    // Initialize lookup tables
    uint32_t F_ASCII[256];
    uint8_t IDX_ASCII[256];
    uint32_t f_rot[4];
    csyncmer_init_ascii_hash_table(F_ASCII);
    csyncmer_init_ascii_to_idx(IDX_ASCII);
    csyncmer_make_f_rot_32(S, f_rot);

    // Ring buffer with extended area to eliminate wrap-around branches
    #define BUF_SIZE 64
    #define BUF_MASK (BUF_SIZE - 1)
    uint32_t buf[BUF_SIZE + 32] __attribute__((aligned(32)));

    const uint8_t* raw_seq = (const uint8_t*)sequence;

    size_t syncmer_count = 0;
    __m256i window_size_minus_1 = _mm256_set1_epi32((int)(window_size - 1));
    __m256i zero = _mm256_setzero_si256();
    __m256i max_val = _mm256_set1_epi32((int)UINT32_MAX);

    // Initialize first hash directly
    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_rotl7_32(hash) ^ F_ASCII[raw_seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[raw_seq[0]]];
    buf[0] = hash;
    buf[BUF_SIZE] = hash;  // Mirror

    // Fill initial window + first batch
    size_t init_size = window_size + 7;
    for (size_t i = 1; i < init_size; ++i) {
        hash = csyncmer_rotl7_32(fw) ^ F_ASCII[raw_seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[raw_seq[i]]];
        size_t buf_idx = i & BUF_MASK;
        buf[buf_idx] = hash;
        if (buf_idx < 32) buf[BUF_SIZE + buf_idx] = hash;
    }
    size_t hash_count = init_size;

    // Process 8 k-mers at a time
    size_t k = 0;
    for (; k + 8 <= num_kmers; k += 8) {
        size_t start_idx = k & BUF_MASK;
        __m256i min_hash = max_val;
        __m256i min_pos = zero;

        __m256i sign = _mm256_set1_epi32((int)0x80000000u);

        for (size_t j = 0; j < window_size; ++j) {
            __m256i cur_hash = _mm256_loadu_si256((__m256i*)&buf[start_idx + j]);
            __m256i cur_pos = _mm256_set1_epi32((int)j);
            __m256i min_signed = _mm256_xor_si256(min_hash, sign);
            __m256i cur_signed = _mm256_xor_si256(cur_hash, sign);
            __m256i is_smaller = _mm256_cmpgt_epi32(min_signed, cur_signed);
            min_hash = _mm256_blendv_epi8(min_hash, cur_hash, is_smaller);
            min_pos = _mm256_blendv_epi8(min_pos, cur_pos, is_smaller);
        }

        __m256i is_first = _mm256_cmpeq_epi32(min_pos, zero);
        __m256i is_last = _mm256_cmpeq_epi32(min_pos, window_size_minus_1);
        __m256i is_syncmer = _mm256_or_si256(is_first, is_last);
        syncmer_count += __builtin_popcount(_mm256_movemask_ps(_mm256_castsi256_ps(is_syncmer)));

        size_t end_hash = k + 8 + 7 + window_size;
        if (end_hash > num_smers) end_hash = num_smers;
        for (; hash_count < end_hash; ++hash_count) {
            hash = csyncmer_rotl7_32(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[raw_seq[hash_count]]];
            size_t buf_idx = hash_count & BUF_MASK;
            buf[buf_idx] = hash;
            if (buf_idx < 32) buf[BUF_SIZE + buf_idx] = hash;
        }
    }

    // Handle remaining k-mers with scalar
    for (; k < num_kmers; ++k) {
        size_t end_hash = k + window_size;
        for (; hash_count < end_hash; ++hash_count) {
            hash = csyncmer_rotl7_32(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[raw_seq[hash_count]]];
            size_t buf_idx = hash_count & BUF_MASK;
            buf[buf_idx] = hash;
            if (buf_idx < 32) buf[BUF_SIZE + buf_idx] = hash;
        }

        uint32_t local_min_hash = UINT32_MAX;
        uint32_t local_min_pos = 0;
        for (size_t j = 0; j < window_size; ++j) {
            uint32_t h = buf[(k + j) & BUF_MASK];
            if (h < local_min_hash) {
                local_min_hash = h;
                local_min_pos = (uint32_t)j;
            }
        }

        if (local_min_pos == 0 || local_min_pos == (uint32_t)(window_size - 1)) {
            syncmer_count++;
        }
    }

    #undef BUF_SIZE
    #undef BUF_MASK

    *num_syncmers = syncmer_count;
    return syncmer_count;
}

#endif  // __AVX2__

// ============================================================================
// HASH-ONLY BENCHMARKS (no syncmer detection)
// ============================================================================

/// Benchmark 32-bit ntHash with 2-bit encoding
inline void benchmark_nthash32_2bit(const char* sequence, size_t S) {
    size_t len = strlen(sequence);
    if (len < S) return;

    auto encoded = encode_sequence_2bit(sequence, len);

    csyncmer_simd::hash::NtHash hasher(S);
    size_t count = 0;

    uint32_t h = hasher.init(encoded, 0);
    count++;

    for (size_t i = 1; i <= len - S; ++i) {
        h = hasher.roll(encoded[i - 1], encoded[i + S - 1], encoded[i]);
        count++;
    }

    printf("[BENCHMARK_HASH_2BIT]:: HASHED %lu S-MERS (last hash: %u)\n", count, h);
}

/// Benchmark ntHash32 with direct ASCII lookup (no encoding step)
inline void benchmark_nthash32_direct(const char* sequence, size_t S) {
    size_t len = strlen(sequence);
    if (len < S) return;

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();

    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);
    uint32_t fw_init = csyncmer_simd::hash::detail::make_fw_init(S);

    size_t num_smers = len - S + 1;
    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    uint32_t fw = fw_init;
    for (size_t i = 0; i < S - 1; ++i) {
        uint32_t fw_out = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i]];
        fw = fw_out ^ f_rot[0];
    }

    uint32_t h = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[S - 1]];
    fw = h ^ f_rot[IDX_ASCII[seq[0]]];
    size_t count = 1;

    for (size_t i = 1; i < num_smers; ++i) {
        uint32_t fw_out = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = fw_out ^ f_rot[IDX_ASCII[seq[i - 1]]];
        h = fw_out;
        count++;
    }

    printf("[BENCHMARK_HASH_DIRECT]:: HASHED %lu S-MERS (last=%u)\n", count, h);
}

/// SIMD ntHash32 - 8 hashers processing 8 disjoint chunks
inline void benchmark_nthash32_simd(const char* sequence, size_t S) {
    size_t len = strlen(sequence);
    if (len < S) return;

    static const auto F_ASCII = make_ascii_hash_table_local();
    static const auto IDX_ASCII = make_ascii_to_idx_local();
    auto f_rot = csyncmer_simd::hash::detail::make_f_rot(S);
    uint32_t fw_init = csyncmer_simd::hash::detail::make_fw_init(S);
    std::array<uint32_t, 4> f_rot_arr = {f_rot[0], f_rot[1], f_rot[2], f_rot[3]};

    size_t num_smers = len - S + 1;
    const uint8_t* seq = reinterpret_cast<const uint8_t*>(sequence);

    size_t chunk_len = (num_smers + 7) / 8;

    std::array<size_t, 8> chunk_starts, chunk_ends;
    for (size_t lane = 0; lane < 8; ++lane) {
        chunk_starts[lane] = lane * chunk_len;
        chunk_ends[lane] = std::min((lane + 1) * chunk_len, num_smers);
    }

    std::array<uint32_t, 8> fw_arr;
    for (size_t lane = 0; lane < 8; ++lane) {
        size_t start = chunk_starts[lane];
        if (start >= num_smers) {
            fw_arr[lane] = 0;
            continue;
        }
        uint32_t fw = fw_init;
        for (size_t i = 0; i < S - 1; ++i) {
            uint32_t fw_out = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[start + i]];
            fw = fw_out ^ f_rot[0];
        }
        uint32_t h = csyncmer_simd::hash::detail::rotl7(fw) ^ F_ASCII[seq[start + S - 1]];
        fw = h ^ f_rot[IDX_ASCII[seq[start]]];
        fw_arr[lane] = fw;
    }

    csyncmer_simd::u32x8 fw = csyncmer_simd::u32x8::from_array(fw_arr);
    csyncmer_simd::u32x8 hash_out = csyncmer_simd::u32x8::zero();

    size_t count = std::min(size_t(8), num_smers);
    for (size_t step = 1; step < chunk_len; ++step) {
        for (size_t lane = 0; lane < 8; ++lane) {
            if (chunk_starts[lane] + step < chunk_ends[lane]) count++;
        }
    }

    for (size_t step = 1; step < chunk_len; ++step) {
        std::array<uint8_t, 8> in_bases, out_idx;
        for (size_t lane = 0; lane < 8; ++lane) {
            size_t pos = chunk_starts[lane] + step;
            if (pos < chunk_ends[lane]) {
                in_bases[lane] = seq[pos + S - 1];
                out_idx[lane] = IDX_ASCII[seq[pos - 1]];
            } else {
                in_bases[lane] = 'A';
                out_idx[lane] = 0;
            }
        }

        csyncmer_simd::u32x8 f_in = csyncmer_simd::simd_gather(F_ASCII, in_bases);
        csyncmer_simd::u32x8 f_rot_out = csyncmer_simd::simd_lookup4(f_rot_arr, out_idx);

        hash_out = csyncmer_simd::simd_rotl7(fw) ^ f_in;
        fw = hash_out ^ f_rot_out;
    }

    auto final_arr = hash_out.to_array();
    volatile uint32_t sink = final_arr[0];
    (void)sink;

    printf("[BENCHMARK_HASH_SIMD]:: HASHED %lu S-MERS\n", count);
}

} // namespace csyncmer_nthash32

// Backward compatibility alias
namespace csyncmer_nthash32_variants = csyncmer_nthash32;

#endif // SYNCMER_NTHASH32_H
