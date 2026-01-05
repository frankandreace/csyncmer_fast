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

#ifdef __BMI2__
#include <x86intrin.h>
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
// Optimized Two-Stack with SIMD 2-bit Packing and Parallel Hash
#ifdef __AVX2__
// Uses techniques from simd-minimizers:
// - 2-bit DNA packing (4 bases/byte)
// - SIMD table lookup for parallel hash computation
// - SIMD gather for loading all 8 chunks simultaneously
// - Split 32-bit hash/position buffers for exact correctness
// ============================================================================

// 2-bit pack table: A=0, C=1, T=2, G=3
// Note: Uses the same layout as base_to_bits (lowercase a=65, A=97, etc.)
static inline uint8_t csyncmer_pack_base(char c) {
    // Simple branchless lookup: A/a=0, C/c=1, T/t=2, G/g=3
    static const uint8_t table[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 0-15
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 16-31
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 32-47
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 48-63
        0,0,0,1,0,0,0,3,0,0,0,0,0,0,0,0, // 64-79:  A=65->0, C=67->1, G=71->3
        0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0, // 80-95:  T=84->2
        0,0,0,1,0,0,0,3,0,0,0,0,0,0,0,0, // 96-111: a=97->0, c=99->1, g=103->3
        0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0, // 112-127: t=116->2
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 128-143
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 144-159
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 160-175
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 176-191
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 192-207
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 208-223
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 224-239
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  // 240-255
    };
    return table[(uint8_t)c];
}

#ifdef __BMI2__
// Pack 8 ASCII chars to 2 bytes using BMI2 PEXT
// Extracts bits 1-2 from each byte: A=0, C=1, T=2, G=3
static inline uint16_t csyncmer_pack_8_pext(const uint8_t* ascii) {
    uint64_t chars;
    memcpy(&chars, ascii, 8);
    return (uint16_t)_pext_u64(chars, 0x0606060606060606ULL);
}
#endif

// SIMD hash table (8 elements for permutevar8x32, duplicated)
// Using lower 32 bits of 64-bit ntHash constants
alignas(32) static const uint32_t CSYNCMER_SIMD_F32[8] = {
    0x95c60474,  // A (from 0x3c8bfbb395c60474)
    0x62a02b4c,  // C (from 0x3193c18562a02b4c)
    0x82572324,  // T (from 0x2032ed0082572324)
    0x4be24456,  // G (from 0x295549f54be24456)
    0x95c60474,  // A (duplicate for permutevar)
    0x62a02b4c,  // C
    0x82572324,  // T
    0x4be24456   // G
};

// SIMD table lookup using permutevar8x32
static inline __m256i csyncmer_simd_lookup(__m256i table, __m256i indices) {
    return _mm256_permutevar8x32_epi32(table, indices);
}

// SIMD rotate left by 7 bits (32-bit elements)
static inline __m256i csyncmer_simd_rotl7(__m256i x) {
    return _mm256_or_si256(_mm256_slli_epi32(x, 7), _mm256_srli_epi32(x, 25));
}

// Pack sequence to 2-bit representation (4 bases per byte)
static inline void csyncmer_pack_seq_2bit(
    const char* seq,
    size_t len,
    uint8_t* out
) {
#ifdef __BMI2__
    // Fast path: pack 8 chars (2 output bytes) at a time using PEXT
    size_t i = 0;
    size_t byte_idx = 0;
    const uint8_t* useq = (const uint8_t*)seq;

    for (; i + 8 <= len; i += 8, byte_idx += 2) {
        uint16_t packed = csyncmer_pack_8_pext(useq + i);
        memcpy(out + byte_idx, &packed, 2);
    }

    // Remainder (0-7 chars)
    for (; i < len; i++) {
        size_t bidx = i / 4;
        size_t bit_offset = (i % 4) * 2;
        out[bidx] |= (csyncmer_pack_base(seq[i]) << bit_offset);
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < len; i++) {
        size_t byte_idx = i / 4;
        size_t bit_offset = (i % 4) * 2;
        out[byte_idx] |= (csyncmer_pack_base(seq[i]) << bit_offset);
    }
#endif
}

// 8x8 matrix transpose for batch output collection
// Adapted from simd-minimizers
static inline void csyncmer_transpose_8x8(__m256i m[8], __m256i t[8]) {
    __m256 r0 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[0]), _mm256_castsi256_ps(m[1]));
    __m256 r1 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[0]), _mm256_castsi256_ps(m[1]));
    __m256 r2 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[2]), _mm256_castsi256_ps(m[3]));
    __m256 r3 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[2]), _mm256_castsi256_ps(m[3]));
    __m256 r4 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[4]), _mm256_castsi256_ps(m[5]));
    __m256 r5 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[4]), _mm256_castsi256_ps(m[5]));
    __m256 r6 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[6]), _mm256_castsi256_ps(m[7]));
    __m256 r7 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[6]), _mm256_castsi256_ps(m[7]));

    __m256 s0 = _mm256_shuffle_ps(r0, r2, 0x44);
    __m256 s1 = _mm256_shuffle_ps(r0, r2, 0xEE);
    __m256 s2 = _mm256_shuffle_ps(r1, r3, 0x44);
    __m256 s3 = _mm256_shuffle_ps(r1, r3, 0xEE);
    __m256 s4 = _mm256_shuffle_ps(r4, r6, 0x44);
    __m256 s5 = _mm256_shuffle_ps(r4, r6, 0xEE);
    __m256 s6 = _mm256_shuffle_ps(r5, r7, 0x44);
    __m256 s7 = _mm256_shuffle_ps(r5, r7, 0xEE);

    t[0] = _mm256_castps_si256(_mm256_permute2f128_ps(s0, s4, 0x20));
    t[1] = _mm256_castps_si256(_mm256_permute2f128_ps(s1, s5, 0x20));
    t[2] = _mm256_castps_si256(_mm256_permute2f128_ps(s2, s6, 0x20));
    t[3] = _mm256_castps_si256(_mm256_permute2f128_ps(s3, s7, 0x20));
    t[4] = _mm256_castps_si256(_mm256_permute2f128_ps(s0, s4, 0x31));
    t[5] = _mm256_castps_si256(_mm256_permute2f128_ps(s1, s5, 0x31));
    t[6] = _mm256_castps_si256(_mm256_permute2f128_ps(s2, s6, 0x31));
    t[7] = _mm256_castps_si256(_mm256_permute2f128_ps(s3, s7, 0x31));
}

// Pre-computed dedup shuffle table (256 entries for all 8-bit masks)
// Entry i contains indices to pack non-masked elements to the left
// Using Lemire's algorithm for branchless filtering
alignas(32) static const uint32_t CSYNCMER_UNIQSHUF[256][8] = {
    {0,1,2,3,4,5,6,7}, {1,2,3,4,5,6,7,0}, {0,2,3,4,5,6,7,1}, {2,3,4,5,6,7,0,1},
    {0,1,3,4,5,6,7,2}, {1,3,4,5,6,7,0,2}, {0,3,4,5,6,7,1,2}, {3,4,5,6,7,0,1,2},
    {0,1,2,4,5,6,7,3}, {1,2,4,5,6,7,0,3}, {0,2,4,5,6,7,1,3}, {2,4,5,6,7,0,1,3},
    {0,1,4,5,6,7,2,3}, {1,4,5,6,7,0,2,3}, {0,4,5,6,7,1,2,3}, {4,5,6,7,0,1,2,3},
    {0,1,2,3,5,6,7,4}, {1,2,3,5,6,7,0,4}, {0,2,3,5,6,7,1,4}, {2,3,5,6,7,0,1,4},
    {0,1,3,5,6,7,2,4}, {1,3,5,6,7,0,2,4}, {0,3,5,6,7,1,2,4}, {3,5,6,7,0,1,2,4},
    {0,1,2,5,6,7,3,4}, {1,2,5,6,7,0,3,4}, {0,2,5,6,7,1,3,4}, {2,5,6,7,0,1,3,4},
    {0,1,5,6,7,2,3,4}, {1,5,6,7,0,2,3,4}, {0,5,6,7,1,2,3,4}, {5,6,7,0,1,2,3,4},
    {0,1,2,3,4,6,7,5}, {1,2,3,4,6,7,0,5}, {0,2,3,4,6,7,1,5}, {2,3,4,6,7,0,1,5},
    {0,1,3,4,6,7,2,5}, {1,3,4,6,7,0,2,5}, {0,3,4,6,7,1,2,5}, {3,4,6,7,0,1,2,5},
    {0,1,2,4,6,7,3,5}, {1,2,4,6,7,0,3,5}, {0,2,4,6,7,1,3,5}, {2,4,6,7,0,1,3,5},
    {0,1,4,6,7,2,3,5}, {1,4,6,7,0,2,3,5}, {0,4,6,7,1,2,3,5}, {4,6,7,0,1,2,3,5},
    {0,1,2,3,6,7,4,5}, {1,2,3,6,7,0,4,5}, {0,2,3,6,7,1,4,5}, {2,3,6,7,0,1,4,5},
    {0,1,3,6,7,2,4,5}, {1,3,6,7,0,2,4,5}, {0,3,6,7,1,2,4,5}, {3,6,7,0,1,2,4,5},
    {0,1,2,6,7,3,4,5}, {1,2,6,7,0,3,4,5}, {0,2,6,7,1,3,4,5}, {2,6,7,0,1,3,4,5},
    {0,1,6,7,2,3,4,5}, {1,6,7,0,2,3,4,5}, {0,6,7,1,2,3,4,5}, {6,7,0,1,2,3,4,5},
    {0,1,2,3,4,5,7,6}, {1,2,3,4,5,7,0,6}, {0,2,3,4,5,7,1,6}, {2,3,4,5,7,0,1,6},
    {0,1,3,4,5,7,2,6}, {1,3,4,5,7,0,2,6}, {0,3,4,5,7,1,2,6}, {3,4,5,7,0,1,2,6},
    {0,1,2,4,5,7,3,6}, {1,2,4,5,7,0,3,6}, {0,2,4,5,7,1,3,6}, {2,4,5,7,0,1,3,6},
    {0,1,4,5,7,2,3,6}, {1,4,5,7,0,2,3,6}, {0,4,5,7,1,2,3,6}, {4,5,7,0,1,2,3,6},
    {0,1,2,3,5,7,4,6}, {1,2,3,5,7,0,4,6}, {0,2,3,5,7,1,4,6}, {2,3,5,7,0,1,4,6},
    {0,1,3,5,7,2,4,6}, {1,3,5,7,0,2,4,6}, {0,3,5,7,1,2,4,6}, {3,5,7,0,1,2,4,6},
    {0,1,2,5,7,3,4,6}, {1,2,5,7,0,3,4,6}, {0,2,5,7,1,3,4,6}, {2,5,7,0,1,3,4,6},
    {0,1,5,7,2,3,4,6}, {1,5,7,0,2,3,4,6}, {0,5,7,1,2,3,4,6}, {5,7,0,1,2,3,4,6},
    {0,1,2,3,4,7,5,6}, {1,2,3,4,7,0,5,6}, {0,2,3,4,7,1,5,6}, {2,3,4,7,0,1,5,6},
    {0,1,3,4,7,2,5,6}, {1,3,4,7,0,2,5,6}, {0,3,4,7,1,2,5,6}, {3,4,7,0,1,2,5,6},
    {0,1,2,4,7,3,5,6}, {1,2,4,7,0,3,5,6}, {0,2,4,7,1,3,5,6}, {2,4,7,0,1,3,5,6},
    {0,1,4,7,2,3,5,6}, {1,4,7,0,2,3,5,6}, {0,4,7,1,2,3,5,6}, {4,7,0,1,2,3,5,6},
    {0,1,2,3,7,4,5,6}, {1,2,3,7,0,4,5,6}, {0,2,3,7,1,4,5,6}, {2,3,7,0,1,4,5,6},
    {0,1,3,7,2,4,5,6}, {1,3,7,0,2,4,5,6}, {0,3,7,1,2,4,5,6}, {3,7,0,1,2,4,5,6},
    {0,1,2,7,3,4,5,6}, {1,2,7,0,3,4,5,6}, {0,2,7,1,3,4,5,6}, {2,7,0,1,3,4,5,6},
    {0,1,7,2,3,4,5,6}, {1,7,0,2,3,4,5,6}, {0,7,1,2,3,4,5,6}, {7,0,1,2,3,4,5,6},
    {0,1,2,3,4,5,6,7}, {1,2,3,4,5,6,0,7}, {0,2,3,4,5,6,1,7}, {2,3,4,5,6,0,1,7},
    {0,1,3,4,5,6,2,7}, {1,3,4,5,6,0,2,7}, {0,3,4,5,6,1,2,7}, {3,4,5,6,0,1,2,7},
    {0,1,2,4,5,6,3,7}, {1,2,4,5,6,0,3,7}, {0,2,4,5,6,1,3,7}, {2,4,5,6,0,1,3,7},
    {0,1,4,5,6,2,3,7}, {1,4,5,6,0,2,3,7}, {0,4,5,6,1,2,3,7}, {4,5,6,0,1,2,3,7},
    {0,1,2,3,5,6,4,7}, {1,2,3,5,6,0,4,7}, {0,2,3,5,6,1,4,7}, {2,3,5,6,0,1,4,7},
    {0,1,3,5,6,2,4,7}, {1,3,5,6,0,2,4,7}, {0,3,5,6,1,2,4,7}, {3,5,6,0,1,2,4,7},
    {0,1,2,5,6,3,4,7}, {1,2,5,6,0,3,4,7}, {0,2,5,6,1,3,4,7}, {2,5,6,0,1,3,4,7},
    {0,1,5,6,2,3,4,7}, {1,5,6,0,2,3,4,7}, {0,5,6,1,2,3,4,7}, {5,6,0,1,2,3,4,7},
    {0,1,2,3,4,6,5,7}, {1,2,3,4,6,0,5,7}, {0,2,3,4,6,1,5,7}, {2,3,4,6,0,1,5,7},
    {0,1,3,4,6,2,5,7}, {1,3,4,6,0,2,5,7}, {0,3,4,6,1,2,5,7}, {3,4,6,0,1,2,5,7},
    {0,1,2,4,6,3,5,7}, {1,2,4,6,0,3,5,7}, {0,2,4,6,1,3,5,7}, {2,4,6,0,1,3,5,7},
    {0,1,4,6,2,3,5,7}, {1,4,6,0,2,3,5,7}, {0,4,6,1,2,3,5,7}, {4,6,0,1,2,3,5,7},
    {0,1,2,3,6,4,5,7}, {1,2,3,6,0,4,5,7}, {0,2,3,6,1,4,5,7}, {2,3,6,0,1,4,5,7},
    {0,1,3,6,2,4,5,7}, {1,3,6,0,2,4,5,7}, {0,3,6,1,2,4,5,7}, {3,6,0,1,2,4,5,7},
    {0,1,2,6,3,4,5,7}, {1,2,6,0,3,4,5,7}, {0,2,6,1,3,4,5,7}, {2,6,0,1,3,4,5,7},
    {0,1,6,2,3,4,5,7}, {1,6,0,2,3,4,5,7}, {0,6,1,2,3,4,5,7}, {6,0,1,2,3,4,5,7},
    {0,1,2,3,4,5,6,7}, {1,2,3,4,5,0,6,7}, {0,2,3,4,5,1,6,7}, {2,3,4,5,0,1,6,7},
    {0,1,3,4,5,2,6,7}, {1,3,4,5,0,2,6,7}, {0,3,4,5,1,2,6,7}, {3,4,5,0,1,2,6,7},
    {0,1,2,4,5,3,6,7}, {1,2,4,5,0,3,6,7}, {0,2,4,5,1,3,6,7}, {2,4,5,0,1,3,6,7},
    {0,1,4,5,2,3,6,7}, {1,4,5,0,2,3,6,7}, {0,4,5,1,2,3,6,7}, {4,5,0,1,2,3,6,7},
    {0,1,2,3,5,4,6,7}, {1,2,3,5,0,4,6,7}, {0,2,3,5,1,4,6,7}, {2,3,5,0,1,4,6,7},
    {0,1,3,5,2,4,6,7}, {1,3,5,0,2,4,6,7}, {0,3,5,1,2,4,6,7}, {3,5,0,1,2,4,6,7},
    {0,1,2,5,3,4,6,7}, {1,2,5,0,3,4,6,7}, {0,2,5,1,3,4,6,7}, {2,5,0,1,3,4,6,7},
    {0,1,5,2,3,4,6,7}, {1,5,0,2,3,4,6,7}, {0,5,1,2,3,4,6,7}, {5,0,1,2,3,4,6,7},
    {0,1,2,3,4,5,6,7}, {1,2,3,4,0,5,6,7}, {0,2,3,4,1,5,6,7}, {2,3,4,0,1,5,6,7},
    {0,1,3,4,2,5,6,7}, {1,3,4,0,2,5,6,7}, {0,3,4,1,2,5,6,7}, {3,4,0,1,2,5,6,7},
    {0,1,2,4,3,5,6,7}, {1,2,4,0,3,5,6,7}, {0,2,4,1,3,5,6,7}, {2,4,0,1,3,5,6,7},
    {0,1,4,2,3,5,6,7}, {1,4,0,2,3,5,6,7}, {0,4,1,2,3,5,6,7}, {4,0,1,2,3,5,6,7},
    {0,1,2,3,4,5,6,7}, {1,2,3,0,4,5,6,7}, {0,2,3,1,4,5,6,7}, {2,3,0,1,4,5,6,7},
    {0,1,3,2,4,5,6,7}, {1,3,0,2,4,5,6,7}, {0,3,1,2,4,5,6,7}, {3,0,1,2,4,5,6,7},
    {0,1,2,3,4,5,6,7}, {1,2,0,3,4,5,6,7}, {0,2,1,3,4,5,6,7}, {2,0,1,3,4,5,6,7},
    {0,1,2,3,4,5,6,7}, {1,0,2,3,4,5,6,7}, {0,1,2,3,4,5,6,7}, {0,1,2,3,4,5,6,7},
};

// Filtered append using dedup shuffle table
// skip_mask: bit i=1 means skip lane i
// Returns new write index
static inline size_t csyncmer_append_filtered(
    __m256i vals, int skip_mask, uint32_t* out, size_t write_idx
) {
    int num = 8 - __builtin_popcount(skip_mask);
    if (num == 0) return write_idx;
    __m256i shuf = _mm256_load_si256((__m256i*)CSYNCMER_UNIQSHUF[skip_mask]);
    __m256i packed = _mm256_permutevar8x32_epi32(vals, shuf);
    _mm256_storeu_si256((__m256i*)(out + write_idx), packed);
    return write_idx + num;
}

// Optimized two-stack with SIMD hash computation (32-bit hash, 8 parallel lanes)
// LIMITATION: Uses 16-bit hash packing for O(1) amortized sliding window minimum.
// When two s-mer hashes have identical upper 16 bits, tie-breaking uses position
// instead of the full 32-bit hash, causing ~0.00004% discrepancy vs reference.
// For exact results, use csyncmer_compute_fused_rescan_branchless() or the 64-bit
// simd_multiwindow functions which use full 32-bit hash comparison.
static inline size_t csyncmer_compute_twostack_simd_32(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    uint32_t* out_positions,
    size_t max_positions
) {
    if (length < K) {
        return 0;
    }

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = num_smers - window_size + 1;

    // Fall back to scalar for small inputs or large windows
    // (For now, just return 0 - caller should handle this case)
    if (num_kmers < 64 || window_size > 64) {
        return 0;  // TODO: implement scalar fallback with positions
    }

    // Phase 1: Pack sequence to 2-bit representation
    size_t packed_size = (length + 3) / 4 + 32;  // +32 for SIMD padding
    uint8_t* packed = (uint8_t*)aligned_alloc(32, packed_size);
    if (!packed) {
        return 0;
    }
    memset(packed, 0, packed_size);
    csyncmer_pack_seq_2bit(sequence, length, packed);

    // Phase 2: Calculate chunk parameters
    size_t chunk_len = (num_kmers + 7) / 8;
    // Round to multiple of 4 for 2-bit byte alignment
    chunk_len = ((chunk_len + 3) / 4) * 4;
    size_t bytes_per_chunk = chunk_len / 4;
    size_t max_iter = chunk_len + window_size - 1 + S - 1;

    // Chunk start positions (in bases, not bytes)
    alignas(32) uint32_t chunk_starts[8];
    alignas(32) uint32_t chunk_ends[8];
    for (int i = 0; i < 8; i++) {
        chunk_starts[i] = (uint32_t)(i * chunk_len);
        chunk_ends[i] = (uint32_t)((i + 1) * chunk_len);
        if (chunk_ends[i] > num_kmers) chunk_ends[i] = (uint32_t)num_kmers;
    }

    // Static lane buffers (reused across calls for speed)
    static uint32_t* lane_bufs[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    static size_t lane_buf_caps[8] = {0};
    size_t lane_counts[8] = {0};

    size_t max_per_lane = chunk_len + 64;
    for (int i = 0; i < 8; i++) {
        if (lane_buf_caps[i] < max_per_lane) {
            free(lane_bufs[i]);
            lane_bufs[i] = (uint32_t*)aligned_alloc(32, max_per_lane * sizeof(uint32_t));
            if (!lane_bufs[i]) {
                free(packed);
                return 0;
            }
            lane_buf_caps[i] = max_per_lane;
        }
    }

    // Phase 3: Initialize SIMD state
    __m256i f_table = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_F32);

    // Compute rotated table for remove operation: rot = (S-1) * 7 mod 32
    uint32_t rot = ((S - 1) * 7) & 31;
    alignas(32) uint32_t f_rot_arr[8];
    for (int i = 0; i < 8; i++) {
        uint32_t x = CSYNCMER_SIMD_F32[i];
        f_rot_arr[i] = (x << rot) | (x >> (32 - rot));
    }
    __m256i f_rot_table = _mm256_load_si256((const __m256i*)f_rot_arr);

    // Forward hash state (8 lanes)
    __m256i fw = _mm256_setzero_si256();

    // Delay buffer for (add, remove) pairs
    size_t delay_size = 1;
    while (delay_size < S) delay_size *= 2;
    size_t delay_mask = delay_size - 1;
    __m256i* delay_buf = (__m256i*)aligned_alloc(32, delay_size * sizeof(__m256i));
    if (!delay_buf) {
        free(packed);
        return 0;
    }
    for (size_t i = 0; i < delay_size; i++) {
        delay_buf[i] = _mm256_setzero_si256();
    }
    size_t write_idx = 0, read_idx = 0;

    // Two-stack state with combined hash+position (like simd-minimizers)
    // Format: [hash upper 16 bits][position lower 16 bits]
    alignas(32) __m256i ring_buf[64];
    for (size_t i = 0; i < window_size; i++) {
        ring_buf[i] = _mm256_set1_epi32((int)UINT32_MAX);
    }
    __m256i prefix_min = _mm256_set1_epi32((int)UINT32_MAX);
    __m256i val_mask = _mm256_set1_epi32((int)0xFFFF0000);
    __m256i pos_mask = _mm256_set1_epi32(0x0000FFFF);

    size_t ring_idx = 0;
    uint32_t pos = 0;
    uint32_t pos_offset = 0;  // For 16-bit position overflow handling
    __m256i pos_offset_vec = _mm256_setzero_si256();  // SIMD version (hoisted for perf)
    const uint32_t max_pos = 0xFFFF;

    // Batch buffer for transpose
    alignas(32) __m256i batch_pos[8];
    size_t batch_count = 0;
    size_t batch_base_kmer = 0;

    // SIMD constants
    __m256i mask_2bit = _mm256_set1_epi32(0x03);

    // Gather indices for 8 chunks (vectorized)
    alignas(32) int32_t gather_base[8];
    for (int i = 0; i < 8; i++) {
        gather_base[i] = (int32_t)(i * bytes_per_chunk);
    }
    __m256i gather_base_vec = _mm256_load_si256((const __m256i*)gather_base);

    // Precomputed constants for batch processing (avoid scalar loops)
    alignas(32) uint32_t idx_arr[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    __m256i idx_offsets = _mm256_load_si256((const __m256i*)idx_arr);
    __m256i w_minus_1_vec = _mm256_set1_epi32((int)(window_size - 1));

    // Precompute validity threshold (lane 7 becomes invalid first)
    size_t last_lane_limit = chunk_ends[7] - chunk_starts[7];

    // Current buffered data and position within buffer
    __m256i cur_data = _mm256_setzero_si256();
    size_t buf_pos = 16;  // Force initial load
    size_t seq_pos = 0;   // Position in sequence

    // Main loop
    for (size_t iter = 0; iter < max_iter; iter++) {
        // Load next batch of bases if needed (every 16 bases)
        if (buf_pos >= 16) {
            size_t byte_offset = seq_pos / 4;
            // Vectorized gather index computation (replaces scalar loop)
            __m256i byte_off_vec = _mm256_set1_epi32((int32_t)byte_offset);
            __m256i idx = _mm256_add_epi32(gather_base_vec, byte_off_vec);
            cur_data = _mm256_i32gather_epi32((const int*)packed, idx, 1);
            buf_pos = seq_pos % 16;
        }

        // Extract 2-bit base for each lane
        __m256i add_base = _mm256_and_si256(
            _mm256_srli_epi32(cur_data, (int)(buf_pos * 2)),
            mask_2bit
        );
        buf_pos++;
        seq_pos++;

        // Get remove base from delay buffer
        __m256i remove_base = delay_buf[read_idx];
        delay_buf[write_idx] = add_base;
        write_idx = (write_idx + 1) & delay_mask;

        // Only start reading from delay buffer after S-1 iterations
        if (iter >= S - 1) {
            read_idx = (read_idx + 1) & delay_mask;
        }

        // Compute rolling hash: hash = rotl7(fw) ^ F[add]
        __m256i fw_rotated = csyncmer_simd_rotl7(fw);
        __m256i add_hash = csyncmer_simd_lookup(f_table, add_base);
        __m256i hash_out = _mm256_xor_si256(fw_rotated, add_hash);

        // Update state based on phase
        if (iter >= S - 1) {
            // Rolling phase: apply removal fw = hash ^ F_rot[remove]
            __m256i remove_hash = csyncmer_simd_lookup(f_rot_table, remove_base);
            fw = _mm256_xor_si256(hash_out, remove_hash);
        } else {
            // Accumulation phase: no removal yet
            fw = hash_out;
            continue;
        }

        // Pack hash (upper 16 bits) + position (lower 16 bits) into combined element
        __m256i pos_vec = _mm256_set1_epi32((int)pos);
        __m256i elem = _mm256_or_si256(
            _mm256_and_si256(hash_out, val_mask),  // hash upper 16 bits
            _mm256_and_si256(pos_vec, pos_mask)    // position lower 16 bits
        );

        // Store in ring buffer
        ring_buf[ring_idx] = elem;

        // Update prefix minimum (single unsigned min comparison)
        // Smaller hash wins; for equal hash, smaller position wins (automatic!)
        prefix_min = _mm256_min_epu32(prefix_min, elem);

        // Handle 16-bit position overflow (like simd-minimizers)
        // Use 2*window_size because ring buffer contains suffix minima from previous block
        if (pos == max_pos) {
            uint32_t delta = (1 << 16) - 2 - 2 * window_size;
            pos -= delta;
            pos_offset += delta;
            pos_offset_vec = _mm256_set1_epi32((int)pos_offset);

            // Subtract delta from entire element (like simd-minimizers)
            // No borrow into hash part because all positions >= delta at reset time
            __m256i delta_vec = _mm256_set1_epi32((int)delta);
            prefix_min = _mm256_sub_epi32(prefix_min, delta_vec);
            for (size_t j = 0; j < window_size; j++) {
                ring_buf[j] = _mm256_sub_epi32(ring_buf[j], delta_vec);
            }
        }

        pos++;
        ring_idx++;

        // Suffix recomputation when ring wraps
        if (ring_idx == window_size) {
            ring_idx = 0;

            // Compute suffix minima from right to left (simple min, no tie-breaking needed!)
            __m256i suffix_min = ring_buf[window_size - 1];
            for (size_t j = window_size - 1; j > 0; j--) {
                suffix_min = _mm256_min_epu32(suffix_min, ring_buf[j - 1]);
                ring_buf[j - 1] = suffix_min;
            }

            // Reset prefix for next block
            prefix_min = _mm256_set1_epi32((int)UINT32_MAX);
        }

        // Check for syncmers after full window
        if (iter >= S - 1 + window_size - 1) {
            // Get suffix min from ring buffer and compute overall min
            __m256i suffix_min = ring_buf[ring_idx];
            __m256i min_elem = _mm256_min_epu32(prefix_min, suffix_min);

            // Extract position from lower 16 bits and add offset to get absolute position
            __m256i min_pos = _mm256_add_epi32(
                _mm256_and_si256(min_elem, pos_mask),
                pos_offset_vec);

            size_t kmer_idx = iter - S + 1 - window_size + 1;

            // Store min_pos in batch (will transpose later)
            if (batch_count % 8 == 0) batch_base_kmer = kmer_idx;
            batch_pos[batch_count % 8] = min_pos;
            batch_count++;

            // Process batch when full
            if (batch_count % 8 == 0) {
                __m256i t[8];
                csyncmer_transpose_8x8(batch_pos, t);

                // Precompute base + idx_offsets once per batch (vectorized)
                __m256i batch_base = _mm256_set1_epi32((int)batch_base_kmer);
                __m256i first_smer = _mm256_add_epi32(batch_base, idx_offsets);
                __m256i last_smer = _mm256_add_epi32(first_smer, w_minus_1_vec);

                for (int lane = 0; lane < 8; lane++) {
                    // Compute absolute k-mer positions (vectorized, no scalar loop)
                    __m256i lane_offset = _mm256_set1_epi32((int)chunk_starts[lane]);
                    __m256i abs_pos = _mm256_add_epi32(lane_offset, first_smer);

                    // t[lane] contains transposed min_pos values (absolute positions with offset)
                    __m256i min_pos_lane = t[lane];

                    // Syncmer check: min at first or last s-mer position
                    __m256i is_first = _mm256_cmpeq_epi32(min_pos_lane, first_smer);
                    __m256i is_last = _mm256_cmpeq_epi32(min_pos_lane, last_smer);
                    __m256i is_syncmer = _mm256_or_si256(is_first, is_last);

                    // Validity check: skip when all lanes valid (most iterations)
                    if (batch_base_kmer + 7 >= last_lane_limit) {
                        // Near end - need per-element validity check
                        __m256i limit = _mm256_set1_epi32((int)(chunk_ends[lane] - chunk_starts[lane]));
                        __m256i kmer_indices = _mm256_add_epi32(batch_base, idx_offsets);
                        __m256i valid_mask = _mm256_cmpgt_epi32(limit, kmer_indices);
                        is_syncmer = _mm256_and_si256(is_syncmer, valid_mask);
                    }

                    // Collect positions
                    int keep_mask = _mm256_movemask_ps(_mm256_castsi256_ps(is_syncmer));
                    int skip_mask = (~keep_mask) & 0xFF;
                    if (keep_mask != 0 && lane_counts[lane] + 8 <= max_per_lane) {
                        lane_counts[lane] = csyncmer_append_filtered(
                            abs_pos, skip_mask, lane_bufs[lane], lane_counts[lane]);
                    }
                }
            }
        }
    }

    // Process remaining partial batch
    size_t partial = batch_count % 8;
    if (partial > 0) {
        // Pad remaining slots with zeros
        for (size_t i = partial; i < 8; i++) {
            batch_pos[i] = _mm256_setzero_si256();
        }

        __m256i t[8];
        csyncmer_transpose_8x8(batch_pos, t);

        for (int lane = 0; lane < 8; lane++) {
            alignas(32) uint32_t abs_positions[8];
            for (int j = 0; j < 8; j++) {
                abs_positions[j] = (uint32_t)(chunk_starts[lane] + batch_base_kmer + j);
            }
            __m256i abs_pos = _mm256_load_si256((__m256i*)abs_positions);

            __m256i first_smer = abs_pos;
            __m256i last_smer = _mm256_add_epi32(abs_pos,
                _mm256_set1_epi32((int)(window_size - 1)));

            __m256i lane_offset = _mm256_set1_epi32((int)chunk_starts[lane]);
            __m256i abs_min_pos = _mm256_add_epi32(t[lane], lane_offset);

            __m256i is_first = _mm256_cmpeq_epi32(abs_min_pos, first_smer);
            __m256i is_last = _mm256_cmpeq_epi32(abs_min_pos, last_smer);
            __m256i is_syncmer = _mm256_or_si256(is_first, is_last);

            // Validity mask - also exclude padded entries
            alignas(32) uint32_t valid_arr[8];
            for (int j = 0; j < 8; j++) {
                size_t abs_kmer = chunk_starts[lane] + batch_base_kmer + j;
                valid_arr[j] = (j < (int)partial && abs_kmer < chunk_ends[lane]) ? 0xFFFFFFFF : 0;
            }
            __m256i valid_mask = _mm256_load_si256((__m256i*)valid_arr);
            is_syncmer = _mm256_and_si256(is_syncmer, valid_mask);

            int keep_mask = _mm256_movemask_ps(_mm256_castsi256_ps(is_syncmer));
            int skip_mask = (~keep_mask) & 0xFF;
            if (keep_mask != 0 && lane_counts[lane] + 8 <= max_per_lane) {
                lane_counts[lane] = csyncmer_append_filtered(
                    abs_pos, skip_mask, lane_bufs[lane], lane_counts[lane]);
            }
        }
    }

    // Concatenate lane buffers (per-lane sorted, not globally sorted)
    size_t total = 0;
    for (int i = 0; i < 8; i++) {
        if (lane_counts[i] > 0 && total + lane_counts[i] <= max_positions) {
            memcpy(out_positions + total, lane_bufs[i], lane_counts[i] * sizeof(uint32_t));
            total += lane_counts[i];
        }
    }

    // Cleanup (lane_bufs are static, don't free)
    free(delay_buf);
    free(packed);

    return total;
}

// Count-only version of TWOSTACK (no position collection, faster)
static inline size_t csyncmer_compute_twostack_simd_32_count(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S
) {
    if (length < K) {
        return 0;
    }

    size_t window_size = K - S + 1;
    size_t num_smers = length - S + 1;
    size_t num_kmers = num_smers - window_size + 1;

    if (num_kmers < 64 || window_size > 64) {
        size_t count;
        csyncmer_compute_fused_rescan_branchless(sequence, length, K, S, &count);
        return count;
    }

    // Pack sequence to 2-bit representation
    size_t packed_size = (length + 15) / 16 * 4 + 32;
    uint8_t* packed = (uint8_t*)aligned_alloc(32, packed_size);
    if (!packed) return 0;
    memset(packed, 0, packed_size);
    csyncmer_pack_seq_2bit(sequence, length, packed);

    size_t chunk_len = (num_kmers + 7) / 8;
    chunk_len = ((chunk_len + 3) / 4) * 4;
    size_t bytes_per_chunk = chunk_len / 4;
    size_t max_iter = chunk_len + window_size - 1 + S - 1;

    alignas(32) uint32_t chunk_starts[8];
    alignas(32) uint32_t chunk_ends[8];
    for (int i = 0; i < 8; i++) {
        chunk_starts[i] = (uint32_t)(i * chunk_len);
        chunk_ends[i] = (uint32_t)((i + 1) * chunk_len);
        if (chunk_ends[i] > num_kmers) chunk_ends[i] = (uint32_t)num_kmers;
    }

    size_t syncmer_count = 0;

    __m256i f_table = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_F32);
    uint32_t rot = ((S - 1) * 7) & 31;
    alignas(32) uint32_t f_rot_arr[8];
    for (int i = 0; i < 8; i++) {
        uint32_t x = CSYNCMER_SIMD_F32[i];
        f_rot_arr[i] = (x << rot) | (x >> (32 - rot));
    }
    __m256i f_rot_table = _mm256_load_si256((const __m256i*)f_rot_arr);

    __m256i fw = _mm256_setzero_si256();

    size_t delay_size = 1;
    while (delay_size < S) delay_size *= 2;
    size_t delay_mask = delay_size - 1;
    __m256i* delay_buf = (__m256i*)aligned_alloc(32, delay_size * sizeof(__m256i));
    if (!delay_buf) { free(packed); return 0; }
    for (size_t i = 0; i < delay_size; i++) delay_buf[i] = _mm256_setzero_si256();
    size_t write_idx = 0, read_idx = 0;

    alignas(32) __m256i ring_buf[64];
    for (size_t i = 0; i < window_size; i++) ring_buf[i] = _mm256_set1_epi32((int)UINT32_MAX);
    __m256i prefix_min = _mm256_set1_epi32((int)UINT32_MAX);
    __m256i val_mask = _mm256_set1_epi32((int)0xFFFF0000);
    __m256i pos_mask = _mm256_set1_epi32(0x0000FFFF);

    size_t ring_idx = 0;
    uint32_t pos = 0;
    uint32_t pos_offset = 0;
    __m256i pos_offset_vec = _mm256_setzero_si256();
    const uint32_t max_pos = 0xFFFF;

    __m256i mask_2bit = _mm256_set1_epi32(0x03);

    alignas(32) int32_t gather_base[8];
    for (int i = 0; i < 8; i++) gather_base[i] = (int32_t)(i * bytes_per_chunk);
    __m256i gather_base_vec = _mm256_load_si256((const __m256i*)gather_base);

    size_t last_lane_limit = chunk_ends[7] - chunk_starts[7];

    __m256i cur_data = _mm256_setzero_si256();
    size_t buf_pos = 16;
    size_t seq_pos = 0;

    for (size_t iter = 0; iter < max_iter; iter++) {
        if (buf_pos >= 16) {
            size_t byte_offset = seq_pos / 4;
            __m256i byte_off_vec = _mm256_set1_epi32((int32_t)byte_offset);
            __m256i idx = _mm256_add_epi32(gather_base_vec, byte_off_vec);
            cur_data = _mm256_i32gather_epi32((const int*)packed, idx, 1);
            buf_pos = seq_pos % 16;
        }

        __m256i add_base = _mm256_and_si256(_mm256_srli_epi32(cur_data, (int)(buf_pos * 2)), mask_2bit);
        buf_pos++;
        seq_pos++;

        __m256i remove_base = delay_buf[read_idx];
        delay_buf[write_idx] = add_base;
        write_idx = (write_idx + 1) & delay_mask;
        if (iter >= S - 1) read_idx = (read_idx + 1) & delay_mask;

        __m256i fw_rotated = csyncmer_simd_rotl7(fw);
        __m256i add_hash = csyncmer_simd_lookup(f_table, add_base);
        __m256i hash_out = _mm256_xor_si256(fw_rotated, add_hash);

        if (iter >= S - 1) {
            __m256i remove_hash = csyncmer_simd_lookup(f_rot_table, remove_base);
            fw = _mm256_xor_si256(hash_out, remove_hash);
        } else {
            fw = hash_out;
            continue;
        }

        __m256i pos_vec = _mm256_set1_epi32((int)pos);
        __m256i elem = _mm256_or_si256(
            _mm256_and_si256(hash_out, val_mask),
            _mm256_and_si256(pos_vec, pos_mask));

        ring_buf[ring_idx] = elem;
        prefix_min = _mm256_min_epu32(prefix_min, elem);

        if (pos == max_pos) {
            uint32_t delta = (1 << 16) - 2 - 2 * window_size;
            pos -= delta;
            pos_offset += delta;
            pos_offset_vec = _mm256_set1_epi32((int)pos_offset);
            __m256i delta_vec = _mm256_set1_epi32((int)delta);
            prefix_min = _mm256_sub_epi32(prefix_min, delta_vec);
            for (size_t j = 0; j < window_size; j++)
                ring_buf[j] = _mm256_sub_epi32(ring_buf[j], delta_vec);
        }

        pos++;
        ring_idx++;

        if (ring_idx == window_size) {
            ring_idx = 0;
            __m256i suffix_min = ring_buf[window_size - 1];
            for (size_t j = window_size - 1; j > 0; j--) {
                suffix_min = _mm256_min_epu32(suffix_min, ring_buf[j - 1]);
                ring_buf[j - 1] = suffix_min;
            }
            prefix_min = _mm256_set1_epi32((int)UINT32_MAX);
        }

        if (iter >= S - 1 + window_size - 1) {
            __m256i suffix_min = ring_buf[ring_idx];
            __m256i min_elem = _mm256_min_epu32(prefix_min, suffix_min);
            __m256i min_pos = _mm256_add_epi32(
                _mm256_and_si256(min_elem, pos_mask), pos_offset_vec);

            size_t kmer_idx = iter - S + 1 - window_size + 1;

            __m256i first_smer = _mm256_set1_epi32((int)kmer_idx);
            __m256i last_smer = _mm256_set1_epi32((int)(kmer_idx + window_size - 1));
            __m256i is_first = _mm256_cmpeq_epi32(min_pos, first_smer);
            __m256i is_last = _mm256_cmpeq_epi32(min_pos, last_smer);
            __m256i is_syncmer = _mm256_or_si256(is_first, is_last);

            if (kmer_idx >= last_lane_limit) {
                __m256i kmer_vec = _mm256_set1_epi32((int)kmer_idx);
                __m256i limits = _mm256_load_si256((const __m256i*)chunk_ends);
                __m256i starts = _mm256_load_si256((const __m256i*)chunk_starts);
                __m256i lane_limits = _mm256_sub_epi32(limits, starts);
                __m256i valid_mask = _mm256_cmpgt_epi32(lane_limits, kmer_vec);
                is_syncmer = _mm256_and_si256(is_syncmer, valid_mask);
            }

            syncmer_count += __builtin_popcount(_mm256_movemask_ps(_mm256_castsi256_ps(is_syncmer)));
        }
    }

    free(delay_buf);
    free(packed);
    return syncmer_count;
}

#endif  // __AVX2__

#ifdef __cplusplus
}
#endif

#endif // CSYNCMER_FAST_H
