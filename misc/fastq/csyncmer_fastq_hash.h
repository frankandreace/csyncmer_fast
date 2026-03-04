#ifndef CSYNCMER_FASTQ_HASH_H
#define CSYNCMER_FASTQ_HASH_H

// Hash-only benchmark: pack + interleave + hash pass, no twostack.
// For measuring raw hash throughput independently.
//
// Optimization: bases are extracted inline from the interleaved buffer
// using two sliding cursors (add/remove), eliminating the separate
// csyncmer_extract_bases() pass and the all_bases temporary array.
// The inner loop is 4x unrolled for better ILP.

#include "../../csyncmer_fast.h"
#include "csyncmer_fastq_pack.h"

#ifdef __AVX2__

// SIMD pack + interleave: convert ASCII sequences directly to interleaved 2-bit
// format in a single pass. Replaces csyncmer_pack_seq_2bit() +
// csyncmer_interleave_packed(), eliminating the intermediate temp_packed buffer.
//
// Encoding: (c >> 1) & 3  — equivalent to PEXT with mask 0x06.
//   A(0x41)→0, C(0x43)→1, T(0x54)→2, G(0x47)→3, N(0x4E)→3
//
// Each group of 16 base positions across 8 lanes is converted using SSE:
//   load 16 ASCII → vpand+vpsrlw+vpand (2-bit) → vpmaddubsw+vpmaddwd (pack) →
//   vpshufb (collect) → scalar store to interleaved[g*32 + lane*4].
static inline void csyncmer_simd_pack_interleave(
    const char* sequences[8],
    const size_t lengths[8],
    const size_t per_lane_kmers[8],
    uint8_t* interleaved,
    size_t n_groups
) {
    if (n_groups == 0) return;

    const __m128i mask_fe    = _mm_set1_epi8((char)0xFE);
    const __m128i mask_03    = _mm_set1_epi8(0x03);
    const __m128i pack_pairs = _mm_set1_epi16(0x0401);     // maddubs: byte pairs → 4-bit
    const __m128i pack_quads = _mm_set1_epi32(0x00100001);  // madd: 16-bit pairs → 8-bit
    const __m128i extract_shuf = _mm_setr_epi8(
        0, 4, 8, 12,  -1,-1,-1,-1,  -1,-1,-1,-1,  -1,-1,-1,-1);

    // Find min length among valid lanes → number of groups needing no bounds check
    size_t min_valid_len = (size_t)-1;
    int n_valid = 0;
    for (int i = 0; i < 8; i++) {
        if (per_lane_kmers[i] > 0) {
            n_valid++;
            if (lengths[i] < min_valid_len)
                min_valid_len = lengths[i];
        }
    }
    size_t full_groups = (min_valid_len != (size_t)-1) ? min_valid_len / 16 : 0;
    if (full_groups > n_groups) full_groups = n_groups;

    // Only zero tail groups when all 8 lanes are valid (common case);
    // otherwise zero the whole buffer for inactive lane slots.
    if (n_valid == 8) {
        for (size_t g = full_groups; g < n_groups; g++)
            memset(interleaved + g * 32, 0, 32);
    } else {
        memset(interleaved, 0, n_groups * 32);
    }

    // Full groups: all valid lanes have ≥16 bases at this position
    for (size_t g = 0; g < full_groups; g++) {
        uint32_t *out = (uint32_t*)(interleaved + g * 32);
        for (int lane = 0; lane < 8; lane++) {
            if (per_lane_kmers[lane] == 0) continue;
            __m128i chars = _mm_loadu_si128(
                (const __m128i*)(sequences[lane] + g * 16));
            __m128i bases = _mm_and_si128(
                _mm_srli_epi16(_mm_and_si128(chars, mask_fe), 1), mask_03);
            __m128i p16 = _mm_maddubs_epi16(bases, pack_pairs);
            __m128i p32 = _mm_madd_epi16(p16, pack_quads);
            __m128i packed = _mm_shuffle_epi8(p32, extract_shuf);
            out[lane] = (uint32_t)_mm_cvtsi128_si32(packed);
        }
    }

    // Tail groups: per-lane bounds checking
    for (size_t g = full_groups; g < n_groups; g++) {
        uint32_t *out = (uint32_t*)(interleaved + g * 32);
        size_t pos = g * 16;
        for (int lane = 0; lane < 8; lane++) {
            if (per_lane_kmers[lane] == 0 || pos >= lengths[lane]) continue;
            size_t remaining = lengths[lane] - pos;
            __m128i chars;
            if (remaining >= 16) {
                chars = _mm_loadu_si128(
                    (const __m128i*)(sequences[lane] + pos));
            } else {
                alignas(16) char tmp[16];
                memset(tmp, 0, 16);
                memcpy(tmp, sequences[lane] + pos, remaining);
                chars = _mm_load_si128((const __m128i*)tmp);
            }
            __m128i bases = _mm_and_si128(
                _mm_srli_epi16(_mm_and_si128(chars, mask_fe), 1), mask_03);
            __m128i p16 = _mm_maddubs_epi16(bases, pack_pairs);
            __m128i p32 = _mm_madd_epi16(p16, pack_quads);
            __m128i packed = _mm_shuffle_epi8(p32, extract_shuf);
            out[lane] = (uint32_t)_mm_cvtsi128_si32(packed);
        }
    }
}

// Hash-only benchmark: pack + interleave + hash pass only, no twostack.
// Returns total s-mer hashes computed (for throughput measurement).
static inline size_t csyncmer_hash_only_multi(
    const char* sequences[8],
    const size_t lengths[8],
    size_t K,
    size_t S,
    uint8_t* work_buf,
    size_t work_buf_size,
    __m256i **out_hash_buf,        // if non-NULL: store hash_buf ptr, skip free
    size_t *out_per_lane_kmers     // if non-NULL: write 8-element per-lane kmer counts
) {
    if (S == 0 || S >= K) return 0;

    size_t _per_lane_kmers[8];
    size_t *per_lane_kmers = out_per_lane_kmers ? out_per_lane_kmers : _per_lane_kmers;
    size_t max_num_kmers = 0;
    int any_valid = 0;
    for (int i = 0; i < 8; i++) {
        if (sequences[i] && lengths[i] >= K) {
            per_lane_kmers[i] = lengths[i] - K + 1;
            if (per_lane_kmers[i] > max_num_kmers)
                max_num_kmers = per_lane_kmers[i];
            any_valid = 1;
        } else {
            per_lane_kmers[i] = 0;
        }
    }
    if (!any_valid || max_num_kmers < 64) return 0;

    size_t max_len = 0;
    for (int i = 0; i < 8; i++)
        if (lengths[i] > max_len) max_len = lengths[i];
    size_t max_num_smers = max_len - S + 1;
    size_t n_groups = (max_len + 15) / 16;
    size_t interleaved_size = (n_groups * 32 + 31) & ~(size_t)31;

    // Buffer layout: interleaved | 2KB pad | hash_buf
    size_t interleaved_hash_pad = 2048;  // break 4K aliasing between interleaved loads and hash_buf stores
    size_t off_hash = interleaved_size + interleaved_hash_pad;
    size_t needed = off_hash +
        ((max_num_smers * sizeof(__m256i) + 31) & ~(size_t)31);

    int own_buf = 0;
    if (!work_buf || work_buf_size < needed) {
        if (out_hash_buf) return 0;  // caller needs buffer alive
        work_buf = (uint8_t*)aligned_alloc(32, needed);
        if (!work_buf) return 0;
        own_buf = 1;
    }
    uint8_t *interleaved = work_buf;
    __m256i *hash_buf    = (__m256i*)(work_buf + off_hash);

    // SIMD pack + interleave (single pass, no temp_packed buffer)
    csyncmer_simd_pack_interleave(sequences, lengths, per_lane_kmers,
                                  interleaved, n_groups);

    // SIMD tables
    __m256i f_table = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_F32);
    uint32_t rot = ((S - 1) * 7) & 31;
    alignas(32) uint32_t f_rot_arr[8];
    for (int i = 0; i < 8; i++) {
        uint32_t x = CSYNCMER_SIMD_F32[i];
        f_rot_arr[i] = (x << rot) | (x >> (32 - rot));
    }
    __m256i f_rot_table = _mm256_load_si256((const __m256i*)f_rot_arr);
    __m256i rc_table_v = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_RC32);
    __m256i rc_rotr7_table = _mm256_load_si256(
        (const __m256i*)CSYNCMER_SIMD_RC32_ROTR7);
    alignas(32) uint32_t c_rot_arr[8];
    for (int i = 0; i < 8; i++) {
        uint32_t cx = CSYNCMER_SIMD_RC32[i];
        c_rot_arr[i] = (cx << rot) | (cx >> (32 - rot));
    }
    __m256i c_rot_table = _mm256_load_si256((const __m256i*)c_rot_arr);
    __m256i mask = _mm256_set1_epi32(0x03);
    // Extract initial S bases into temp array (init is one-time cost)
    __m256i init_bases[S];  // VLA — S is small (typically 8-31)
    {
        size_t ig = 0;
        int ileft = 16;
        __m256i iw = _mm256_load_si256((const __m256i*)interleaved);
        for (size_t j = 0; j < S; j++) {
            init_bases[j] = _mm256_and_si256(iw, mask);
            iw = _mm256_srli_epi32(iw, 2);
            if (--ileft == 0) {
                ileft = 16;
                ig++;
                iw = _mm256_load_si256(
                    (const __m256i*)(interleaved + ig * 32));
            }
        }
    }

    // Initial forward hash (bases 0..S-1)
    __m256i fw = _mm256_setzero_si256();
    for (size_t j = 0; j < S; j++)
        fw = _mm256_xor_si256(csyncmer_simd_rotl7(fw),
                              csyncmer_simd_lookup(f_table, init_bases[j]));
    __m256i fw_hash_0 = fw;

    // Initial RC hash
    __m256i rc = _mm256_setzero_si256();
    for (size_t j = 0; j < S; j++)
        rc = _mm256_xor_si256(csyncmer_simd_rotl7(rc),
                              csyncmer_simd_lookup(rc_table_v, init_bases[S - 1 - j]));

    {
        __m256i h0 = _mm256_min_epu32(fw_hash_0, rc);
        _mm256_store_si256(&hash_buf[0], h0);
    }

    // Set up rolling state
    fw = _mm256_xor_si256(fw_hash_0,
                          csyncmer_simd_lookup(f_rot_table, init_bases[0]));
    __m256i prev_rm = init_bases[0];

    // ---- Inline extraction cursors (raw variables, no struct) ----
    // Add cursor: base at position (si + S - 1); for si=1 that's position S
    size_t add_grp = S / 16;
    int add_left = 16 - (int)(S % 16);
    __m256i add_word;
    {
        __m256i aw = _mm256_load_si256(
            (const __m256i*)(interleaved + add_grp * 32));
        int add_off = (int)(S % 16);
        add_word = (add_off > 0) ? _mm256_srli_epi32(aw, add_off * 2) : aw;
    }

    // Remove cursor: base at position si; for si=1 that's position 1
    size_t rm_grp = 0;
    int rm_left = 15;
    __m256i rm_word = _mm256_srli_epi32(
        _mm256_load_si256((const __m256i*)interleaved), 2);

    // ---- Rolling hash loop (4x unrolled) ----
    // Cursor state as raw locals so the compiler can keep them in registers.
    // Each cursor: extract base = word & mask, then word >>= 2.
    // Every 16th extraction, reload next interleaved group.

#define EXTRACT_ADD(var) do {                                               \
    var = _mm256_and_si256(add_word, mask);                                  \
    add_word = _mm256_srli_epi32(add_word, 2);                              \
    if (__builtin_expect(--add_left == 0, 0)) {                             \
        add_left = 16; add_grp++;                                           \
        add_word = _mm256_load_si256(                                       \
            (const __m256i*)(interleaved + add_grp * 32));                  \
    }                                                                       \
} while(0)

#define EXTRACT_RM(var) do {                                                \
    var = _mm256_and_si256(rm_word, mask);                                   \
    rm_word = _mm256_srli_epi32(rm_word, 2);                                \
    if (__builtin_expect(--rm_left == 0, 0)) {                              \
        rm_left = 16; rm_grp++;                                             \
        rm_word = _mm256_load_si256(                                        \
            (const __m256i*)(interleaved + rm_grp * 32));                   \
    }                                                                       \
} while(0)

#define HASH_ONE_ITER(SI) do {                                              \
    __m256i add_base, remove_base;                                          \
    EXTRACT_ADD(add_base);                                                  \
    EXTRACT_RM(remove_base);                                                \
                                                                            \
    __m256i fw_hash = _mm256_xor_si256(                                     \
        csyncmer_simd_rotl7(fw),                                            \
        csyncmer_simd_lookup(f_table, add_base));                           \
    fw = _mm256_xor_si256(fw_hash,                                          \
        csyncmer_simd_lookup(f_rot_table, remove_base));                    \
                                                                            \
    __m256i rc_hash = _mm256_xor_si256(                                     \
        _mm256_xor_si256(csyncmer_simd_rotr7(rc),                           \
            csyncmer_simd_lookup(rc_rotr7_table, prev_rm)),                 \
        csyncmer_simd_lookup(c_rot_table, add_base));                       \
    rc = rc_hash;                                                           \
    prev_rm = remove_base;                                                  \
                                                                            \
    __m256i can = _mm256_min_epu32(fw_hash, rc_hash);                       \
    _mm256_store_si256(&hash_buf[(SI)], can);                              \
} while(0)

    for (size_t si = 1; si < max_num_smers; si++) {
        HASH_ONE_ITER(si);
    }
#undef EXTRACT_ADD
#undef EXTRACT_RM
#undef HASH_ONE_ITER

    if (out_hash_buf)
        *out_hash_buf = hash_buf;
    else if (own_buf)
        free(work_buf);
    return max_num_smers;
}

// Pack-only benchmark: measure just pack + interleave cost (no hashing).
// Returns total bases processed (for throughput measurement).
static inline size_t csyncmer_pack_only_multi(
    const char* sequences[8],
    const size_t lengths[8],
    size_t K,
    size_t S,
    uint8_t* work_buf,
    size_t work_buf_size
) {
    if (S == 0 || S >= K) return 0;

    size_t per_lane_kmers[8];
    int any_valid = 0;
    for (int i = 0; i < 8; i++) {
        if (sequences[i] && lengths[i] >= K) {
            per_lane_kmers[i] = lengths[i] - K + 1;
            any_valid = 1;
        } else {
            per_lane_kmers[i] = 0;
        }
    }
    if (!any_valid) return 0;

    size_t max_len = 0;
    for (int i = 0; i < 8; i++)
        if (lengths[i] > max_len) max_len = lengths[i];
    size_t n_groups = (max_len + 15) / 16;
    size_t interleaved_size = (n_groups * 32 + 31) & ~(size_t)31;

    size_t needed = interleaved_size;

    int own_buf = 0;
    if (!work_buf || work_buf_size < needed) {
        work_buf = (uint8_t*)aligned_alloc(32, needed);
        if (!work_buf) return 0;
        own_buf = 1;
    }
    uint8_t *interleaved = work_buf;

    // SIMD pack + interleave (single pass, no temp_packed buffer)
    csyncmer_simd_pack_interleave(sequences, lengths, per_lane_kmers,
                                  interleaved, n_groups);

    size_t total_bp = 0;
    for (int i = 0; i < 8; i++) total_bp += lengths[i];

    if (own_buf) free(work_buf);
    return total_bp;
}

#endif  // __AVX2__

#endif  // CSYNCMER_FASTQ_HASH_H
