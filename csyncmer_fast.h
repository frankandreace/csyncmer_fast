#ifndef CSYNCMER_FAST_H
#define CSYNCMER_FAST_H

// Fast closed syncmer detection using ntHash32 - Pure C implementation
// Includes AVX2 SIMD acceleration when available

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// ============================================================================
// ntHash32 Constants and Helpers
// ============================================================================

static const uint32_t CSYNCMER_NTHASH32_F[4] = {
    0x95c60474,  // A
    0x62a02b4c,  // C
    0x82572324,  // T
    0x4be24456   // G
};

static inline uint32_t csyncmer_rotl7(uint32_t x) {
    return (x << 7) | (x >> 25);
}

static inline void csyncmer_make_f_rot(size_t S, uint32_t f_rot[4]) {
    // ntHash uses (k-1)*R rotation, matching nthash.hpp::make_f_rot
    uint32_t rot = (uint32_t)(((S - 1) * 7) & 31);
    for (int i = 0; i < 4; i++) {
        uint32_t x = CSYNCMER_NTHASH32_F[i];
        f_rot[i] = (x << rot) | (x >> (32 - rot));
    }
}

// ============================================================================
// ASCII Lookup Tables
// ============================================================================

static inline void csyncmer_init_ascii_hash_table(uint32_t table[256]) {
    memset(table, 0, 256 * sizeof(uint32_t));
    table['A'] = table['a'] = CSYNCMER_NTHASH32_F[0];
    table['C'] = table['c'] = CSYNCMER_NTHASH32_F[1];
    table['T'] = table['t'] = CSYNCMER_NTHASH32_F[2];
    table['G'] = table['g'] = CSYNCMER_NTHASH32_F[3];
}

static inline void csyncmer_init_ascii_to_idx(uint8_t table[256]) {
    memset(table, 0, 256);
    table['A'] = table['a'] = 0;
    table['C'] = table['c'] = 1;
    table['T'] = table['t'] = 2;
    table['G'] = table['g'] = 3;
}

// ============================================================================
// Scalar Implementation: Fused RESCAN with Branch-free Updates
// ============================================================================

static inline size_t csyncmer_compute_fused_rescan_branchless(
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
    uint32_t F_ASCII[256];
    uint8_t IDX_ASCII[256];
    uint32_t f_rot[4];
    csyncmer_init_ascii_hash_table(F_ASCII);
    csyncmer_init_ascii_to_idx(IDX_ASCII);
    csyncmer_make_f_rot(S, f_rot);

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;

    const uint8_t* seq = (const uint8_t*)sequence;

    // Power-of-2 circular buffer
    size_t buf_size = 1;
    while (buf_size < window_size) buf_size <<= 1;
    size_t buf_mask = buf_size - 1;
    uint32_t* hash_buffer = (uint32_t*)malloc(buf_size * sizeof(uint32_t));
    if (!hash_buffer) {
        *num_syncmers = 0;
        return 0;
    }

    size_t syncmer_count = 0;

    // First hash - compute directly
    uint32_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_rotl7(hash) ^ F_ASCII[seq[j]];
    }
    uint32_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];  // State uses first base of s-mer
    hash_buffer[0] = hash;

    // Fill initial window
    for (size_t i = 1; i < window_size; ++i) {
        hash = csyncmer_rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];  // Use new first base
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
        hash = csyncmer_rotl7(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];  // Use new first base
        hash_buffer[i & buf_mask] = hash;

        // RESCAN: update minimum
        if (min_pos < kmer_idx) {
            // Minimum fell out - rescan
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
            uint32_t is_smaller = (hash < min_hash) ? 1 : 0;
            min_hash = is_smaller ? hash : min_hash;
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
// SIMD Implementation: Multi-Window AVX2
// ============================================================================

#ifdef __AVX2__

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
    // Main buffer: 64 elements (power of 2)
    // Extended area: 32 elements mirroring positions 0..31
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
        // SIMD minimum finding
        size_t start_idx = k & BUF_MASK;
        __m256i min_hash = max_val;
        __m256i min_pos = zero;

        // Sign bit for unsigned comparison (AVX2 only has signed cmpgt)
        __m256i sign = _mm256_set1_epi32((int)0x80000000u);

        for (size_t j = 0; j < window_size; ++j) {
            __m256i cur_hash = _mm256_loadu_si256((__m256i*)&buf[start_idx + j]);
            __m256i cur_pos = _mm256_set1_epi32((int)j);
            // Unsigned comparison: flip sign bits before signed compare
            __m256i min_signed = _mm256_xor_si256(min_hash, sign);
            __m256i cur_signed = _mm256_xor_si256(cur_hash, sign);
            __m256i is_smaller = _mm256_cmpgt_epi32(min_signed, cur_signed);
            min_hash = _mm256_blendv_epi8(min_hash, cur_hash, is_smaller);
            min_pos = _mm256_blendv_epi8(min_pos, cur_pos, is_smaller);
        }

        // Check closed syncmer condition
        __m256i is_first = _mm256_cmpeq_epi32(min_pos, zero);
        __m256i is_last = _mm256_cmpeq_epi32(min_pos, window_size_minus_1);
        __m256i is_syncmer = _mm256_or_si256(is_first, is_last);
        syncmer_count += __builtin_popcount(_mm256_movemask_ps(_mm256_castsi256_ps(is_syncmer)));

        // Compute next 8 hashes
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

    // Handle remaining k-mers with scalar rescan using buffer
    for (; k < num_kmers; ++k) {
        // Ensure we have hashes for this k-mer's window
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

#endif  // __AVX2__

// ============================================================================
// ntHash64 Constants and Helpers
// ============================================================================

static const uint64_t CSYNCMER_NTHASH64_F[4] = {
    0x3c8bfbb395c60474ULL,  // A
    0x3193c18562a02b4cULL,  // C
    0x20323ed082572324ULL,  // T
    0x295549f54be24456ULL   // G
};

static inline uint64_t csyncmer_rotl7_64(uint64_t x) {
    return (x << 7) | (x >> 57);
}

static inline void csyncmer_make_f_rot_64(size_t S, uint64_t f_rot[4]) {
    // ntHash uses (k-1)*R rotation
    uint32_t rot = (uint32_t)(((S - 1) * 7) & 63);
    for (int i = 0; i < 4; i++) {
        uint64_t x = CSYNCMER_NTHASH64_F[i];
        f_rot[i] = (x << rot) | (x >> (64 - rot));
    }
}

static inline void csyncmer_init_ascii_hash_table_64(uint64_t table[256]) {
    memset(table, 0, 256 * sizeof(uint64_t));
    table['A'] = table['a'] = CSYNCMER_NTHASH64_F[0];
    table['C'] = table['c'] = CSYNCMER_NTHASH64_F[1];
    table['T'] = table['t'] = CSYNCMER_NTHASH64_F[2];
    table['G'] = table['g'] = CSYNCMER_NTHASH64_F[3];
}

// ============================================================================
// Shared Static Lookup Tables for 64-bit Iterator
// ============================================================================

static uint64_t CSYNCMER_F_ASCII_64[256];
static uint8_t CSYNCMER_IDX_ASCII[256];
static int CSYNCMER_TABLES_INITIALIZED = 0;

static inline void csyncmer_ensure_tables_initialized(void) {
    if (!CSYNCMER_TABLES_INITIALIZED) {
        csyncmer_init_ascii_hash_table_64(CSYNCMER_F_ASCII_64);
        csyncmer_init_ascii_to_idx(CSYNCMER_IDX_ASCII);
        CSYNCMER_TABLES_INITIALIZED = 1;
    }
}

// ============================================================================
// 64-bit Scalar Implementation
// ============================================================================

static inline size_t csyncmer_compute_fused_rescan_branchless_64(
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
    csyncmer_init_ascii_hash_table_64(F_ASCII);
    csyncmer_init_ascii_to_idx(IDX_ASCII);
    csyncmer_make_f_rot_64(S, f_rot);

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
        hash = csyncmer_rotl7_64(hash) ^ F_ASCII[seq[j]];
    }
    uint64_t fw = hash ^ f_rot[IDX_ASCII[seq[0]]];
    hash_buffer[0] = hash;

    // Fill initial window
    for (size_t i = 1; i < window_size; ++i) {
        hash = csyncmer_rotl7_64(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;
    }

    // Find minimum in first window (using lower 32 bits for comparison, matching SIMD)
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
        hash = csyncmer_rotl7_64(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;

        // RESCAN: update minimum (using lower 32 bits for comparison, matching SIMD)
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
// 64-bit Iterator API
// ============================================================================

// Optimized iterator state - shared lookup tables, local variable caching
typedef struct CsyncmerIterator64 {
    const uint8_t* seq;
    uint64_t* hash_buffer;
    uint64_t hash;
    uint64_t fw;
    size_t min_pos;
    size_t kmer_idx;
    size_t num_kmers;
    size_t window_size;
    size_t S;
    size_t buf_mask;
    uint32_t min_hash32;
    uint64_t f_rot[4];
} CsyncmerIterator64;

// Create iterator for closed syncmer detection
// K = total syncmer length, S = smer size (window_size = K - S + 1)
static inline CsyncmerIterator64* csyncmer_iterator_create_64(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S
) {
    if (length < K || S >= K) {
        return NULL;
    }

    // Ensure shared lookup tables are initialized
    csyncmer_ensure_tables_initialized();

    CsyncmerIterator64* iter = (CsyncmerIterator64*)malloc(sizeof(CsyncmerIterator64));
    if (!iter) return NULL;

    iter->seq = (const uint8_t*)sequence;
    iter->S = S;
    iter->window_size = K - S + 1;
    iter->num_kmers = length - K + 1;

    // Initialize S-dependent rotation table
    csyncmer_make_f_rot_64(S, iter->f_rot);

    // Allocate power-of-2 circular buffer
    size_t buf_size = 1;
    while (buf_size < iter->window_size) buf_size <<= 1;
    iter->buf_mask = buf_size - 1;
    iter->hash_buffer = (uint64_t*)malloc(buf_size * sizeof(uint64_t));
    if (!iter->hash_buffer) {
        free(iter);
        return NULL;
    }

    // Use shared static tables for initialization
    const uint64_t* F_ASCII = CSYNCMER_F_ASCII_64;
    const uint8_t* IDX_ASCII = CSYNCMER_IDX_ASCII;

    // Compute first s-mer hash
    uint64_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_rotl7_64(hash) ^ F_ASCII[iter->seq[j]];
    }
    iter->hash = hash;
    iter->fw = hash ^ iter->f_rot[IDX_ASCII[iter->seq[0]]];
    iter->hash_buffer[0] = hash;

    // Fill initial window with remaining s-mer hashes
    for (size_t i = 1; i < iter->window_size; ++i) {
        hash = csyncmer_rotl7_64(iter->fw) ^ F_ASCII[iter->seq[i + S - 1]];
        iter->fw = hash ^ iter->f_rot[IDX_ASCII[iter->seq[i]]];
        iter->hash_buffer[i & iter->buf_mask] = hash;
    }
    iter->hash = hash;

    // Find minimum in first window
    uint32_t min_hash32 = UINT32_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < iter->window_size; ++i) {
        uint32_t h32 = (uint32_t)iter->hash_buffer[i];
        if (h32 < min_hash32) {
            min_hash32 = h32;
            min_pos = i;
        }
    }
    iter->min_hash32 = min_hash32;
    iter->min_pos = min_pos;
    iter->kmer_idx = 0;

    return iter;
}

// Get next syncmer position. Returns 1 if valid, 0 when exhausted.
static inline int csyncmer_iterator_next_64(
    CsyncmerIterator64* iter,
    size_t* pos
) {
    if (!iter) return 0;

    // Cache hot state in local variables for register allocation
    const uint8_t* seq = iter->seq;
    size_t kmer_idx = iter->kmer_idx;
    const size_t num_kmers = iter->num_kmers;
    const size_t window_size = iter->window_size;
    const size_t S = iter->S;
    uint64_t hash = iter->hash;
    uint64_t fw = iter->fw;
    uint32_t min_hash32 = iter->min_hash32;
    size_t min_pos = iter->min_pos;
    const size_t buf_mask = iter->buf_mask;
    uint64_t* hash_buffer = iter->hash_buffer;
    const uint64_t* f_rot = iter->f_rot;

    // Use shared static tables
    const uint64_t* F_ASCII = CSYNCMER_F_ASCII_64;
    const uint8_t* IDX_ASCII = CSYNCMER_IDX_ASCII;

    // Check first k-mer on first call
    if (kmer_idx == 0) {
        size_t min_offset = min_pos;
        if (min_offset == 0 || min_offset == window_size - 1) {
            *pos = 0;
            iter->kmer_idx = 1;
            return 1;
        }
        kmer_idx = 1;
    }

    // Main loop
    while (kmer_idx < num_kmers) {
        size_t i = kmer_idx + window_size - 1;

        // Compute new hash for the new s-mer entering the window
        hash = csyncmer_rotl7_64(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;

        // RESCAN: update minimum
        uint32_t hash32 = (uint32_t)hash;
        if (min_pos < kmer_idx) {
            // Minimum fell out of window - rescan
            min_hash32 = UINT32_MAX;
            for (size_t j = kmer_idx; j <= i; ++j) {
                uint32_t h32 = (uint32_t)hash_buffer[j & buf_mask];
                if (h32 < min_hash32) {
                    min_hash32 = h32;
                    min_pos = j;
                }
            }
        } else {
            // Check if new hash is smaller
            if (hash32 < min_hash32) {
                min_hash32 = hash32;
                min_pos = i;
            }
        }

        // Check closed syncmer condition
        size_t min_offset = min_pos - kmer_idx;
        size_t current_kmer = kmer_idx;
        kmer_idx++;

        if (min_offset == 0 || min_offset == window_size - 1) {
            // Write back modified state
            iter->hash = hash;
            iter->fw = fw;
            iter->min_hash32 = min_hash32;
            iter->min_pos = min_pos;
            iter->kmer_idx = kmer_idx;
            *pos = current_kmer;
            return 1;
        }
    }

    // Write back state on exhaustion
    iter->hash = hash;
    iter->fw = fw;
    iter->min_hash32 = min_hash32;
    iter->min_pos = min_pos;
    iter->kmer_idx = kmer_idx;
    return 0;  // Exhausted
}

// Free iterator
static inline void csyncmer_iterator_destroy_64(CsyncmerIterator64* iter) {
    if (iter) {
        if (iter->hash_buffer) {
            free(iter->hash_buffer);
        }
        free(iter);
    }
}

// ============================================================================
// 64-bit SIMD Implementation (8-lane using lower 32-bits for comparison)
// ============================================================================

#ifdef __AVX2__

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

    // Fall back to scalar for small inputs
    if (num_kmers < 8) {
        return csyncmer_compute_fused_rescan_branchless_64(
            sequence, length, K, S, num_syncmers);
    }

    // Initialize lookup tables
    uint64_t F_ASCII[256];
    uint8_t IDX_ASCII[256];
    uint64_t f_rot[4];
    csyncmer_init_ascii_hash_table_64(F_ASCII);
    csyncmer_init_ascii_to_idx(IDX_ASCII);
    csyncmer_make_f_rot_64(S, f_rot);

    // Ring buffer for 64-bit hashes, but we'll extract lower 32 bits for SIMD comparison
    #define BUF_SIZE_64 64
    #define BUF_MASK_64 (BUF_SIZE_64 - 1)
    uint64_t buf64[BUF_SIZE_64 + 32] __attribute__((aligned(32)));
    // Separate buffer for lower 32 bits (for 8-lane SIMD)
    uint32_t buf32[BUF_SIZE_64 + 32] __attribute__((aligned(32)));

    const uint8_t* raw_seq = (const uint8_t*)sequence;

    size_t syncmer_count = 0;
    __m256i window_size_minus_1 = _mm256_set1_epi32((int)(window_size - 1));
    __m256i zero = _mm256_setzero_si256();
    __m256i max_val = _mm256_set1_epi32((int)UINT32_MAX);

    // Initialize first hash
    uint64_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_rotl7_64(hash) ^ F_ASCII[raw_seq[j]];
    }
    uint64_t fw = hash ^ f_rot[IDX_ASCII[raw_seq[0]]];
    buf64[0] = hash;
    buf64[BUF_SIZE_64] = hash;
    buf32[0] = (uint32_t)hash;
    buf32[BUF_SIZE_64] = (uint32_t)hash;

    // Fill initial window + first batch
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

    // Sign bit for unsigned 32-bit comparison
    __m256i sign = _mm256_set1_epi32((int)0x80000000u);

    // Process 8 k-mers at a time (using lower 32 bits)
    size_t k = 0;
    for (; k + 8 <= num_kmers; k += 8) {
        size_t start_idx = k & BUF_MASK_64;
        __m256i min_hash = max_val;
        __m256i min_pos = zero;

        for (size_t j = 0; j < window_size; ++j) {
            __m256i cur_hash = _mm256_loadu_si256((__m256i*)&buf32[start_idx + j]);
            __m256i cur_pos = _mm256_set1_epi32((int)j);
            // Unsigned comparison via sign flip
            __m256i min_signed = _mm256_xor_si256(min_hash, sign);
            __m256i cur_signed = _mm256_xor_si256(cur_hash, sign);
            __m256i is_smaller = _mm256_cmpgt_epi32(min_signed, cur_signed);
            min_hash = _mm256_blendv_epi8(min_hash, cur_hash, is_smaller);
            min_pos = _mm256_blendv_epi8(min_pos, cur_pos, is_smaller);
        }

        // Check closed syncmer condition
        __m256i is_first = _mm256_cmpeq_epi32(min_pos, zero);
        __m256i is_last = _mm256_cmpeq_epi32(min_pos, window_size_minus_1);
        __m256i is_syncmer = _mm256_or_si256(is_first, is_last);
        syncmer_count += __builtin_popcount(_mm256_movemask_ps(_mm256_castsi256_ps(is_syncmer)));

        // Compute next 8 hashes
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

    // Handle remaining k-mers with scalar (using 32-bit comparison to match SIMD)
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

#ifdef __cplusplus
}
#endif

#endif // CSYNCMER_FAST_H
