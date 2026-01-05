// Deprecated SIMD multiwindow implementations
// These were slower than TWOSTACK and have been moved here for reference.
//
// To use: #include "csyncmer_fast.h" first, then this file
// Requires AVX2 support

#ifndef CSYNCMER_DEPRECATED_SIMD_MW_H
#define CSYNCMER_DEPRECATED_SIMD_MW_H

#ifdef __AVX2__

// ============================================================================
// 32-bit SIMD Multi-Window Implementation (DEPRECATED)
// Performance: ~192 MB/s (slower than TWOSTACK @ 551 MB/s)
// ============================================================================

static inline size_t csyncmer_compute_simd_multiwindow(
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
        return csyncmer_compute_fused_rescan_branchless(
            sequence, length, K, S, num_syncmers);
    }

    // Initialize lookup tables
    uint32_t F_ASCII[256];
    uint8_t IDX_ASCII[256];
    uint32_t f_rot[4];
    csyncmer_init_ascii_hash_table(F_ASCII);
    csyncmer_init_ascii_to_idx(IDX_ASCII);
    csyncmer_make_f_rot(S, f_rot);

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
        hash = csyncmer_rotl7(hash) ^ F_ASCII[raw_seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[raw_seq[0]]];
    buf[0] = hash;
    buf[BUF_SIZE] = hash;  // Mirror

    // Fill initial window + first batch
    size_t init_size = window_size + 7;
    for (size_t i = 1; i < init_size; ++i) {
        hash = csyncmer_rotl7(fw) ^ F_ASCII[raw_seq[i + S - 1]];
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
            hash = csyncmer_rotl7(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
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
            hash = csyncmer_rotl7(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
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

// ============================================================================
// 64-bit SIMD Multi-Window Implementation (DEPRECATED)
// Performance: ~165 MB/s (slower than 64-bit scalar @ 213 MB/s)
// ============================================================================

static inline size_t csyncmer_compute_simd_multiwindow_64(
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
        return csyncmer_compute_fused_rescan_branchless_64(
            sequence, length, K, S, num_syncmers);
    }

    uint64_t F_ASCII[256];
    uint8_t IDX_ASCII[256];
    uint64_t f_rot[4];
    csyncmer_init_ascii_hash_table_64(F_ASCII);
    csyncmer_init_ascii_to_idx(IDX_ASCII);
    csyncmer_make_f_rot_64(S, f_rot);

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
        hash = csyncmer_rotl7_64(hash) ^ F_ASCII[raw_seq[j]];
    }
    uint64_t fw = hash ^ f_rot[IDX_ASCII[raw_seq[0]]];
    buf64[0] = hash;
    buf64[BUF_SIZE_64] = hash;
    buf32[0] = (uint32_t)hash;
    buf32[BUF_SIZE_64] = (uint32_t)hash;

    size_t init_size = window_size + 7;
    for (size_t i = 1; i < init_size; ++i) {
        hash = csyncmer_rotl7_64(fw) ^ F_ASCII[raw_seq[i + S - 1]];
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
            hash = csyncmer_rotl7_64(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
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
            hash = csyncmer_rotl7_64(fw) ^ F_ASCII[raw_seq[hash_count + S - 1]];
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

#endif  // CSYNCMER_DEPRECATED_SIMD_MW_H
