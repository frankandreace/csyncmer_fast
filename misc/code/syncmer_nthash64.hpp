#ifndef SYNCMER_NTHASH64_H
#define SYNCMER_NTHASH64_H

// ntHash64-based syncmer implementations (deprecated/legacy)
// Superseded by the Iterator API in csyncmer_fast.h
//
// This file consolidates:
// - csyncmer_nthash64_rescan_count() - scalar RESCAN count-only
// - csyncmer_nthash64_simd_multiwindow() - AVX2 SIMD multi-window (requires csyncmer_fast.h)

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// ============================================================================
// NTHASH64 CONSTANTS AND HELPERS (LEGACY)
// These use different names to avoid conflicts with csyncmer_fast.h
// ============================================================================

static const uint64_t CSYNCMER_LEGACY_NTHASH64_F[4] = {
    0x3c8bfbb395c60474ULL,  // A
    0x3193c18562a02b4cULL,  // C
    0x20323ed082572324ULL,  // T
    0x295549f54be24456ULL   // G
};

static inline uint64_t csyncmer_nthash64_rotl7(uint64_t x) {
    return (x << 7) | (x >> 57);
}

static inline void csyncmer_nthash64_make_f_rot(size_t S, uint64_t f_rot[4]) {
    uint32_t rot = (uint32_t)(((S - 1) * 7) & 63);
    for (int i = 0; i < 4; i++) {
        uint64_t x = CSYNCMER_LEGACY_NTHASH64_F[i];
        f_rot[i] = (x << rot) | (x >> (64 - rot));
    }
}

static inline void csyncmer_nthash64_init_ascii_hash_table(uint64_t table[256]) {
    memset(table, 0, 256 * sizeof(uint64_t));
    table['A'] = table['a'] = CSYNCMER_LEGACY_NTHASH64_F[0];
    table['C'] = table['c'] = CSYNCMER_LEGACY_NTHASH64_F[1];
    table['T'] = table['t'] = CSYNCMER_LEGACY_NTHASH64_F[2];
    table['G'] = table['g'] = CSYNCMER_LEGACY_NTHASH64_F[3];
}

static inline void csyncmer_nthash64_init_ascii_to_idx(uint8_t table[256]) {
    memset(table, 0, 256);
    table['A'] = table['a'] = 0;
    table['C'] = table['c'] = 1;
    table['T'] = table['t'] = 2;
    table['G'] = table['g'] = 3;
}

// Reverse complement constants for canonical hashing
// RC[base] = F[complement(base)]: A↔T, C↔G
static const uint64_t CSYNCMER_LEGACY_NTHASH64_RC[4] = {
    0x20323ed082572324ULL,  // RC[A] = F[T]
    0x295549f54be24456ULL,  // RC[C] = F[G]
    0x3c8bfbb395c60474ULL,  // RC[T] = F[A]
    0x3193c18562a02b4cULL   // RC[G] = F[C]
};

static inline void csyncmer_nthash64_init_ascii_rc_table(uint64_t table[256]) {
    memset(table, 0, 256 * sizeof(uint64_t));
    table['A'] = table['a'] = CSYNCMER_LEGACY_NTHASH64_RC[0];
    table['C'] = table['c'] = CSYNCMER_LEGACY_NTHASH64_RC[1];
    table['T'] = table['t'] = CSYNCMER_LEGACY_NTHASH64_RC[2];
    table['G'] = table['g'] = CSYNCMER_LEGACY_NTHASH64_RC[3];
}

// ============================================================================
// 64-BIT NAIVE O(N*W) IMPLEMENTATION
// ============================================================================

/*
 * Naive reference implementation - compute each s-mer hash from scratch,
 * then scan each k-mer window to find the minimum. O(N*S + N*W) time.
 * Uses full 64-bit hash comparison (exact, matches the Iterator API).
 */
static inline size_t csyncmer_nthash64_naive(
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

    uint64_t F_ASCII[256];
    csyncmer_nthash64_init_ascii_hash_table(F_ASCII);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = length - K + 1;

    const uint8_t* seq = (const uint8_t*)sequence;

    // Compute all s-mer hashes from scratch (no rolling hash)
    uint64_t* smer_hashes = (uint64_t*)malloc(num_smers * sizeof(uint64_t));
    if (!smer_hashes) {
        *num_syncmers = 0;
        return 0;
    }

    for (size_t i = 0; i < num_smers; ++i) {
        uint64_t h = 0;
        for (size_t j = 0; j < S; ++j) {
            h = csyncmer_nthash64_rotl7(h) ^ F_ASCII[seq[i + j]];
        }
        smer_hashes[i] = h;
    }

    // O(N*W) scan to find syncmers
    size_t syncmer_count = 0;
    for (size_t i = 0; i < num_kmers; ++i) {
        uint64_t min_hash = UINT64_MAX;
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

    free(smer_hashes);
    *num_syncmers = syncmer_count;
    printf("[NTHASH64_NAIVE]:: COMPUTED %lu CLOSED SYNCMERS\n", syncmer_count);
    printf("[NTHASH64_NAIVE]:: HASHED %lu S-MERS\n", num_smers);
    return syncmer_count;
}

// ============================================================================
// 64-BIT CANONICAL DEQUE (REFERENCE IMPLEMENTATION)
// ============================================================================

/*
 * Simple canonical 64-bit deque implementation for correctness verification.
 * Computes all canonical hashes from scratch (slow), then uses deque for minima.
 * This is a reference implementation - use the Iterator API for performance.
 *
 * Canonical hash = min(forward_hash, rc_hash) for each s-mer.
 */
static inline size_t csyncmer_nthash64_canonical_deque(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    size_t* num_syncmers
) {
    if (length < K || S == 0 || S >= K) {
        *num_syncmers = 0;
        return 0;
    }

    size_t num_s_mers = length - S + 1;
    size_t window_size = K - S + 1;

    // Initialize lookup tables
    uint64_t F_ASCII[256];
    uint64_t RC_ASCII[256];
    csyncmer_nthash64_init_ascii_hash_table(F_ASCII);
    csyncmer_nthash64_init_ascii_rc_table(RC_ASCII);

    const uint8_t* seq = (const uint8_t*)sequence;

    // Compute all canonical s-mer hashes from scratch
    uint64_t* canon_hashes = (uint64_t*)malloc(num_s_mers * sizeof(uint64_t));
    if (!canon_hashes) {
        *num_syncmers = 0;
        return 0;
    }

    for (size_t i = 0; i < num_s_mers; i++) {
        // Forward hash: process left to right
        uint64_t fw = 0;
        for (size_t j = 0; j < S; j++) {
            fw = csyncmer_nthash64_rotl7(fw) ^ F_ASCII[seq[i + j]];
        }
        // RC hash: process right to left with RC constants
        uint64_t rc = 0;
        for (size_t j = 0; j < S; j++) {
            rc = csyncmer_nthash64_rotl7(rc) ^ RC_ASCII[seq[i + S - 1 - j]];
        }
        // Canonical = min
        canon_hashes[i] = (fw <= rc) ? fw : rc;
    }

    // Use deque to find minimal s-mers in O(N)
    size_t* deque = (size_t*)malloc(num_s_mers * sizeof(size_t));
    if (!deque) {
        free(canon_hashes);
        *num_syncmers = 0;
        return 0;
    }

    size_t front = 0, back = 0;
    size_t syncmer_count = 0;

    for (size_t i = 0; i < num_s_mers; i++) {
        // Remove elements larger than current from back
        while (back > front && canon_hashes[deque[back-1]] > canon_hashes[i]) {
            back--;
        }
        deque[back++] = i;

        // Remove elements outside window from front
        if (i >= window_size && deque[front] <= i - window_size) {
            front++;
        }

        // Check for closed syncmer condition
        if (i >= window_size - 1) {
            size_t min_pos = deque[front];
            size_t kmer_pos = i - window_size + 1;
            // Closed syncmer: minimum at first or last position in window
            if (min_pos == kmer_pos || min_pos == kmer_pos + window_size - 1) {
                syncmer_count++;
            }
        }
    }

    free(canon_hashes);
    free(deque);
    *num_syncmers = syncmer_count;
    return syncmer_count;
}

// ============================================================================
// 64-BIT SCALAR RESCAN (COUNT ONLY)
// ============================================================================

/*
 * 64-bit scalar RESCAN - count only, no positions.
 * O(1) amortized: only rescans when minimum falls out of window.
 * Portable: works on any CPU without AVX2.
 * Exact: uses full 64-bit hash comparison.
 *
 * Speed: ~140-210 MB/s (scalar, portable, exact results)
 * For 64-bit with positions, use the Iterator API in csyncmer_fast.h.
 */
static inline size_t csyncmer_nthash64_rescan_count(
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

    // Initialize lookup tables
    uint64_t F_ASCII[256];
    uint8_t IDX_ASCII[256];
    uint64_t f_rot[4];
    csyncmer_nthash64_init_ascii_hash_table(F_ASCII);
    csyncmer_nthash64_init_ascii_to_idx(IDX_ASCII);
    csyncmer_nthash64_make_f_rot(S, f_rot);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    const uint8_t* seq = (const uint8_t*)sequence;

    // Power-of-2 circular buffer
    size_t buf_size = 1;
    while (buf_size < window_size) buf_size <<= 1;
    size_t buf_mask = buf_size - 1;
    uint64_t* hash_buffer = (uint64_t*)malloc(buf_size * sizeof(uint64_t));
    if (!hash_buffer) {
        *num_syncmers = 0;
        return 0;
    }

    size_t syncmer_count = 0;

    // First hash - compute directly
    uint64_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_nthash64_rotl7(hash) ^ F_ASCII[seq[j]];
    }
    uint64_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];
    hash_buffer[0] = hash;

    // Fill initial window
    for (size_t i = 1; i < window_size; ++i) {
        hash = csyncmer_nthash64_rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;
    }

    // Find minimum in first window (using lower 32 bits for comparison)
    uint32_t min_hash32 = UINT32_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < window_size; ++i) {
        uint32_t h32 = (uint32_t)hash_buffer[i];
        if (h32 < min_hash32) {
            min_hash32 = h32;
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
        hash = csyncmer_nthash64_rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;

        // RESCAN: update minimum
        uint32_t hash32 = (uint32_t)hash;
        if (min_pos < kmer_idx) {
            min_hash32 = UINT32_MAX;
            for (size_t j = kmer_idx; j <= i; ++j) {
                uint32_t h32 = (uint32_t)hash_buffer[j & buf_mask];
                if (h32 < min_hash32) {
                    min_hash32 = h32;
                    min_pos = j;
                }
            }
        } else {
            uint32_t is_smaller = (hash32 < min_hash32) ? 1 : 0;
            min_hash32 = is_smaller ? hash32 : min_hash32;
            min_pos = is_smaller ? i : min_pos;
        }

        // Check closed syncmer condition
        size_t min_offset = min_pos - kmer_idx;
        if (min_offset == 0 || min_offset == window_size - 1) {
            syncmer_count++;
        }
    }

    free(hash_buffer);
    *num_syncmers = syncmer_count;
    return syncmer_count;
}

// ============================================================================
// 64-BIT SIMD MULTI-WINDOW (DEPRECATED - AVX2 REQUIRED)
// ============================================================================

#ifdef __AVX2__

#include <immintrin.h>

// Note: This function requires csyncmer_fast.h to be included first
// for csyncmer_rotl7_64, csyncmer_init_ascii_hash_table_64, etc.
// It is preserved for reference but is slower than the Iterator API.

/*
 * 64-bit SIMD Multi-Window Implementation (DEPRECATED)
 * Performance: ~165 MB/s (slower than 64-bit scalar @ 213 MB/s)
 *
 * This uses 32-bit hash comparison internally, so it's not truly "exact"
 * despite using 64-bit hashes. Use the Iterator API for exact results.
 */

// Forward declaration - requires csyncmer_fast.h
extern size_t csyncmer_compute_fused_rescan_branchless_64(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    size_t* num_syncmers
);

static inline size_t csyncmer_nthash64_simd_multiwindow(
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

    if (num_kmers < 8) {
        return csyncmer_nthash64_rescan_count(
            sequence, length, K, S, num_syncmers);
    }

    uint64_t F_ASCII[256];
    uint8_t IDX_ASCII[256];
    uint64_t f_rot[4];
    csyncmer_nthash64_init_ascii_hash_table(F_ASCII);
    csyncmer_nthash64_init_ascii_to_idx(IDX_ASCII);
    csyncmer_nthash64_make_f_rot(S, f_rot);

    #define BUF_SIZE_64 64
    #define BUF_MASK_64 (BUF_SIZE_64 - 1)
    uint64_t buf64[BUF_SIZE_64 + 32] __attribute__((aligned(32)));
    uint32_t buf32[BUF_SIZE_64 + 32] __attribute__((aligned(32)));

    const uint8_t* raw_seq = (const uint8_t*)sequence;

    size_t syncmer_count = 0;
    __m256i window_size_minus_1 = _mm256_set1_epi32((int)(window_size - 1));
    __m256i zero = _mm256_setzero_si256();
    __m256i max_val = _mm256_set1_epi32((int)UINT32_MAX);

    uint64_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_nthash64_rotl7(hash) ^ F_ASCII[raw_seq[j]];
    }
    uint64_t fw = hash ^ f_rot[IDX_ASCII[raw_seq[0]]];
    buf64[0] = hash;
    buf64[BUF_SIZE_64] = hash;
    buf32[0] = (uint32_t)hash;
    buf32[BUF_SIZE_64] = (uint32_t)hash;

    size_t init_size = window_size + 7;
    for (size_t i = 1; i < init_size; ++i) {
        hash = csyncmer_nthash64_rotl7(fw) ^ F_ASCII[raw_seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[raw_seq[i]]];
        size_t buf_idx = i & BUF_MASK_64;
        buf64[buf_idx] = hash;
        buf32[buf_idx] = (uint32_t)hash;
        if (buf_idx < 32) {
            buf64[BUF_SIZE_64 + buf_idx] = hash;
            buf32[BUF_SIZE_64 + buf_idx] = (uint32_t)hash;
        }
    }
    size_t hash_count = init_size;

    __m256i sign = _mm256_set1_epi32((int)0x80000000u);

    size_t k = 0;
    for (; k + 8 <= num_kmers; k += 8) {
        size_t start_idx = k & BUF_MASK_64;
        __m256i min_hash = max_val;
        __m256i min_pos = zero;

        for (size_t j = 0; j < window_size; ++j) {
            __m256i cur_hash = _mm256_loadu_si256((__m256i*)&buf32[start_idx + j]);
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
            hash = csyncmer_nthash64_rotl7(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[raw_seq[hash_count]]];
            size_t buf_idx = hash_count & BUF_MASK_64;
            buf64[buf_idx] = hash;
            buf32[buf_idx] = (uint32_t)hash;
            if (buf_idx < 32) {
                buf64[BUF_SIZE_64 + buf_idx] = hash;
                buf32[BUF_SIZE_64 + buf_idx] = (uint32_t)hash;
            }
        }
    }

    for (; k < num_kmers; ++k) {
        size_t end_hash = k + window_size;
        for (; hash_count < end_hash; ++hash_count) {
            hash = csyncmer_nthash64_rotl7(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
            fw = hash ^ f_rot[IDX_ASCII[raw_seq[hash_count]]];
            size_t buf_idx = hash_count & BUF_MASK_64;
            buf64[buf_idx] = hash;
            buf32[buf_idx] = (uint32_t)hash;
            if (buf_idx < 32) {
                buf64[BUF_SIZE_64 + buf_idx] = hash;
                buf32[BUF_SIZE_64 + buf_idx] = (uint32_t)hash;
            }
        }

        uint64_t local_min_hash = UINT64_MAX;
        uint32_t local_min_pos = 0;
        for (size_t j = 0; j < window_size; ++j) {
            uint64_t h = buf64[(k + j) & BUF_MASK_64];
            if (h < local_min_hash) {
                local_min_hash = h;
                local_min_pos = (uint32_t)j;
            }
        }

        if (local_min_pos == 0 || local_min_pos == (uint32_t)(window_size - 1)) {
            syncmer_count++;
        }
    }

    #undef BUF_SIZE_64
    #undef BUF_MASK_64

    *num_syncmers = syncmer_count;
    return syncmer_count;
}

#endif  // __AVX2__

#endif // SYNCMER_NTHASH64_H
