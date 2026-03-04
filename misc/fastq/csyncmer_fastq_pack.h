#ifndef CSYNCMER_FASTQ_PACK_H
#define CSYNCMER_FASTQ_PACK_H

// Packing, interleaving, and base extraction helpers for multi-read SIMD.
// Transforms 8 separately-packed 2-bit sequences into formats optimized
// for SIMD hash computation.

#include "../../csyncmer_fast.h"

#ifdef __AVX2__

// ============================================================================
// Two-pass multi-read SIMD: separate hash and twostack passes
// ============================================================================
// Splits the monolithic loop into two passes to reduce L1d cache pressure:
//   Pass 1: Hash computation — sequential reads from interleaved packed data,
//           sequential writes to hash_buf. No ring buffer in working set.
//   Pass 2: Twostack sliding minimum — sequential reads from hash_buf,
//           ring buffer fits entirely in L1d (32 KB for w=1022).
//
// Also uses interleaved packing layout: for each group of 16 base positions,
// 8 contiguous int32s (one per lane). Replaces vpgatherdd with vmovdqa.

// Interleave 8 separately-packed reads into contiguous layout.
// Input: packed[lane * ps_uniform + g*4] has lane's group g (16 bases).
// Output: interleaved[g * 32 + lane * 4] = packed[lane * ps_uniform + g*4].
static inline void csyncmer_interleave_packed(
    const uint8_t* packed,
    size_t ps_uniform,
    uint8_t* interleaved,
    size_t n_groups
) {
    alignas(32) int32_t offsets[8];
    for (int i = 0; i < 8; i++)
        offsets[i] = (int32_t)(i * ps_uniform);
    __m256i off_vec = _mm256_load_si256((const __m256i*)offsets);
    for (size_t g = 0; g < n_groups; g++) {
        __m256i byte_off = _mm256_set1_epi32((int32_t)(g * 4));
        __m256i idx = _mm256_add_epi32(off_vec, byte_off);
        __m256i data = _mm256_i32gather_epi32(
            (const int*)packed, idx, 1);
        _mm256_store_si256((__m256i*)(interleaved + g * 32), data);
    }
}

// Extract all bases from interleaved buffer to flat __m256i array.
// Uses constant-shift unrolling: 16 bases per word with shifts 0,2,...,30.
// Output: out[0..count-1] where out[i] has the 2-bit base index (0-3) for
// position (start + i) across all 8 lanes.
static inline void csyncmer_extract_bases(
    const uint8_t* interleaved,
    size_t start,
    size_t count,
    __m256i* out
) {
    __m256i mask = _mm256_set1_epi32(0x03);
    size_t grp = start / 16;
    size_t off = start % 16;
    size_t i = 0;

    // Handle unaligned head
    if (off != 0) {
        __m256i word = _mm256_load_si256(
            (const __m256i*)(interleaved + grp * 32));
        while (off < 16 && i < count) {
            out[i] = _mm256_and_si256(
                _mm256_srli_epi32(word, (int)(off * 2)), mask);
            off++;
            i++;
        }
        grp++;
    }

    // Aligned full groups of 16 (constant shifts, no variable shift)
    while (i + 16 <= count) {
        __m256i w = _mm256_load_si256(
            (const __m256i*)(interleaved + grp * 32));
        out[i +  0] = _mm256_and_si256(w, mask);
        out[i +  1] = _mm256_and_si256(_mm256_srli_epi32(w,  2), mask);
        out[i +  2] = _mm256_and_si256(_mm256_srli_epi32(w,  4), mask);
        out[i +  3] = _mm256_and_si256(_mm256_srli_epi32(w,  6), mask);
        out[i +  4] = _mm256_and_si256(_mm256_srli_epi32(w,  8), mask);
        out[i +  5] = _mm256_and_si256(_mm256_srli_epi32(w, 10), mask);
        out[i +  6] = _mm256_and_si256(_mm256_srli_epi32(w, 12), mask);
        out[i +  7] = _mm256_and_si256(_mm256_srli_epi32(w, 14), mask);
        out[i +  8] = _mm256_and_si256(_mm256_srli_epi32(w, 16), mask);
        out[i +  9] = _mm256_and_si256(_mm256_srli_epi32(w, 18), mask);
        out[i + 10] = _mm256_and_si256(_mm256_srli_epi32(w, 20), mask);
        out[i + 11] = _mm256_and_si256(_mm256_srli_epi32(w, 22), mask);
        out[i + 12] = _mm256_and_si256(_mm256_srli_epi32(w, 24), mask);
        out[i + 13] = _mm256_and_si256(_mm256_srli_epi32(w, 26), mask);
        out[i + 14] = _mm256_and_si256(_mm256_srli_epi32(w, 28), mask);
        out[i + 15] = _mm256_and_si256(_mm256_srli_epi32(w, 30), mask);
        i += 16;
        grp++;
    }

    // Handle tail
    if (i < count) {
        __m256i word = _mm256_load_si256(
            (const __m256i*)(interleaved + grp * 32));
        size_t off2 = 0;
        while (i < count) {
            out[i] = _mm256_and_si256(
                _mm256_srli_epi32(word, (int)(off2 * 2)), mask);
            off2++;
            i++;
        }
    }
}

#endif  // __AVX2__

#endif  // CSYNCMER_FASTQ_PACK_H
