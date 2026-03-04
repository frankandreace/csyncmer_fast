#ifndef CSYNCMER_FASTQ_TWOSTACK_H
#define CSYNCMER_FASTQ_TWOSTACK_H

// Twostack sliding-window minimum syncmer detection (multi-read SIMD).
// Contains both monolithic (single-pass) and two-pass (hash + twostack) variants.

#include "../../csyncmer_fast.h"
#include "csyncmer_fastq_pack.h"

#ifdef __AVX2__

// ============================================================================
// Multi-read SIMD: 8 different reads in 8 lanes (optimal for short reads)
// ============================================================================
// Instead of splitting one long sequence into 8 chunks, this puts 8 separate
// reads into the 8 AVX2 lanes. Each lane processes its own read independently.
// This eliminates the 60% warmup overhead for HiFi reads (~15kb).
// Only the canonical+positions variant is generated (the one syng uses).

#define CSYNCMER_DEFINE_TWOSTACK_SIMD_32_MULTI(FUNC_NAME, CANONICAL, COLLECT_POSITIONS) \
static inline void FUNC_NAME(                                                  \
    const char* sequences[8],                                                  \
    const size_t lengths[8],                                                   \
    size_t K,                                                                  \
    size_t S,                                                                  \
    uint32_t* out_positions[8],                                                \
    uint8_t*  out_strands[8],                                                  \
    size_t max_per_read,                                                       \
    size_t out_counts[8],                                                      \
    uint8_t* work_buf,                                                         \
    size_t work_buf_size                                                       \
) {                                                                            \
    for (int i = 0; i < 8; i++) out_counts[i] = 0;                            \
    if (S == 0 || S >= K) return;                                              \
                                                                               \
    size_t window_size = K - S + 1;                                            \
                                                                               \
    /* Compute per-lane num_kmers and max across lanes */                      \
    size_t per_lane_kmers[8];                                                  \
    size_t max_num_kmers = 0;                                                  \
    int any_valid = 0;                                                         \
    for (int i = 0; i < 8; i++) {                                              \
        if (sequences[i] && lengths[i] >= K) {                                 \
            per_lane_kmers[i] = lengths[i] - K + 1;                            \
            if (per_lane_kmers[i] > max_num_kmers)                             \
                max_num_kmers = per_lane_kmers[i];                             \
            any_valid = 1;                                                     \
        } else {                                                               \
            per_lane_kmers[i] = 0;                                             \
        }                                                                      \
    }                                                                          \
    if (!any_valid) return;                                                     \
                                                                               \
    /* Fall back to single-read scalar for very short reads */                 \
    if (max_num_kmers < 64) {                                                  \
        for (int i = 0; i < 8; i++) {                                          \
            if (per_lane_kmers[i] > 0) {                                       \
                if (CANONICAL && COLLECT_POSITIONS) {                           \
                    out_counts[i] = csyncmer_canonical_rescan_32_(             \
                        sequences[i], lengths[i], K, S,                        \
                        out_positions[i], out_strands[i], max_per_read);       \
                }                                                              \
            }                                                                  \
        }                                                                      \
        return;                                                                \
    }                                                                          \
                                                                               \
    /* Phase 1: Pack each read into concatenated buffer */                     \
    /* Use uniform region size based on longest read so gather doesn't        */ \
    /* read past the allocation for shorter lanes.                            */ \
    size_t packed_offsets[8];                                                   \
    size_t max_len = 0;                                                        \
    for (int i = 0; i < 8; i++)                                                \
        if (lengths[i] > max_len) max_len = lengths[i];                        \
    size_t ps_uniform = ((max_len + 3) / 4 + 32 + 31) & ~(size_t)31;          \
    size_t total_packed = 8 * ps_uniform;                                      \
    for (int i = 0; i < 8; i++)                                                \
        packed_offsets[i] = i * ps_uniform;                                    \
                                                                               \
    size_t delay_size = 1;                                                     \
    while (delay_size < S) delay_size *= 2;                                    \
    size_t delay_mask = delay_size - 1;                                        \
                                                                               \
    /* Compute work buffer layout: packed | delay_buf | ring_buf | strand */   \
    size_t off_packed = 0;                                                     \
    size_t off_delay = (total_packed + 32 + 31) & ~(size_t)31;                \
    size_t off_ring = off_delay +                                              \
        ((delay_size * sizeof(__m256i) + 31) & ~(size_t)31);                  \
    size_t off_strand = off_ring +                                             \
        ((window_size * sizeof(__m256i) + 31) & ~(size_t)31);                 \
    size_t needed = off_strand + ((window_size + 31) & ~(size_t)31);          \
                                                                               \
    int own_buf = 0;                                                           \
    if (!work_buf || work_buf_size < needed) {                                 \
        work_buf = (uint8_t*)aligned_alloc(32, needed);                        \
        if (!work_buf) return;                                                 \
        own_buf = 1;                                                           \
    }                                                                          \
    uint8_t *packed = work_buf + off_packed;                                    \
    __m256i *delay_buf = (__m256i*)(work_buf + off_delay);                     \
    __m256i *ring_buf = (__m256i*)(work_buf + off_ring);                       \
    uint8_t *strand_ring = work_buf + off_strand;                              \
                                                                               \
    memset(packed, 0, total_packed + 32);                                       \
    for (int i = 0; i < 8; i++)                                                \
        if (per_lane_kmers[i] > 0)                                             \
            csyncmer_pack_seq_2bit(sequences[i], lengths[i],                   \
                                   packed + packed_offsets[i]);                 \
                                                                               \
    /* Phase 2: Per-lane chunk parameters (each lane = one read) */            \
    alignas(32) uint32_t chunk_starts[8];                                      \
    alignas(32) uint32_t chunk_ends[8];                                        \
    for (int i = 0; i < 8; i++) {                                              \
        chunk_starts[i] = 0;                                                   \
        chunk_ends[i] = (uint32_t)per_lane_kmers[i];                           \
    }                                                                          \
                                                                               \
    size_t max_iter = max_num_kmers + window_size - 1 + S - 1;                \
    size_t max_per_lane = max_num_kmers + 64;                                  \
                                                                               \
    /* Thread-local lane buffers */                                            \
    size_t lane_counts[8] = {0};                                               \
                                                                               \
    /* Multi-read: write directly to out_positions/out_strands, no TLS */      \
    (void)max_per_lane;                                                        \
                                                                               \
    size_t syncmer_count = 0;                                                  \
                                                                               \
    /* Phase 3: Initialize SIMD state (identical to single-read macro) */      \
    __m256i f_table = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_F32);    \
                                                                               \
    uint32_t rot = ((S - 1) * 7) & 31;                                        \
    alignas(32) uint32_t f_rot_arr[8];                                         \
    for (int i = 0; i < 8; i++) {                                              \
        uint32_t x = CSYNCMER_SIMD_F32[i];                                    \
        f_rot_arr[i] = (x << rot) | (x >> (32 - rot));                        \
    }                                                                          \
    __m256i f_rot_table = _mm256_load_si256((const __m256i*)f_rot_arr);        \
                                                                               \
    __m256i rc_table_v = _mm256_setzero_si256();                               \
    __m256i rc_rotr7_table = _mm256_setzero_si256();                           \
    __m256i c_rot_table = _mm256_setzero_si256();                              \
    if (CANONICAL) {                                                           \
        rc_table_v = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_RC32);    \
        rc_rotr7_table = _mm256_load_si256(                                    \
            (const __m256i*)CSYNCMER_SIMD_RC32_ROTR7);                        \
        alignas(32) uint32_t c_rot_arr[8];                                     \
        for (int i = 0; i < 8; i++) {                                          \
            uint32_t cx = CSYNCMER_SIMD_RC32[i];                               \
            c_rot_arr[i] = (cx << rot) | (cx >> (32 - rot));                   \
        }                                                                      \
        c_rot_table = _mm256_load_si256((const __m256i*)c_rot_arr);            \
    }                                                                          \
                                                                               \
    __m256i fw = _mm256_setzero_si256();                                       \
    __m256i rc = _mm256_setzero_si256();                                       \
    __m256i prev_remove_base = _mm256_setzero_si256();                         \
                                                                               \
    for (size_t i = 0; i < delay_size; i++)                                    \
        delay_buf[i] = _mm256_setzero_si256();                                \
    size_t wr_idx = 0, rd_idx = 0;                                            \
                                                                               \
    for (size_t i = 0; i < window_size; i++)                                   \
        ring_buf[i] = _mm256_set1_epi32((int)UINT32_MAX);                     \
    __m256i prefix_min = _mm256_set1_epi32((int)UINT32_MAX);                   \
    __m256i val_mask = _mm256_set1_epi32((int)0xFFFF0000);                     \
    __m256i pos_mask = _mm256_set1_epi32(0x0000FFFF);                          \
                                                                               \
    uint8_t prefix_strand = 0;                                                 \
    if (CANONICAL && COLLECT_POSITIONS) {                                       \
        for (size_t i = 0; i < window_size; i++) strand_ring[i] = 0;          \
    }                                                                          \
                                                                               \
    size_t ring_idx = 0;                                                       \
    uint32_t pos = 0;                                                          \
    uint32_t pos_offset = 0;                                                   \
    __m256i pos_offset_vec = _mm256_setzero_si256();                           \
    const uint32_t max_pos_val = 0xFFFF;                                       \
                                                                               \
    alignas(32) __m256i batch_pos[8];                                          \
    size_t batch_count = 0;                                                    \
    size_t batch_base_kmer = 0;                                                \
    alignas(32) uint32_t idx_arr[8] = {0, 1, 2, 3, 4, 5, 6, 7};              \
    __m256i idx_offsets = _mm256_load_si256((const __m256i*)idx_arr);           \
    __m256i w_minus_1_vec = _mm256_set1_epi32((int)(window_size - 1));         \
                                                                               \
    __m256i mask_2bit = _mm256_set1_epi32(0x03);                               \
    /* Gather base: each lane reads from its own read's packed data */         \
    alignas(32) int32_t gather_base[8];                                        \
    for (int i = 0; i < 8; i++)                                                \
        gather_base[i] = (int32_t)packed_offsets[i];                           \
    __m256i gather_base_vec = _mm256_load_si256((const __m256i*)gather_base);   \
                                                                               \
    /* Validity: min lane length across all active lanes */                    \
    size_t last_lane_limit = max_num_kmers;                                    \
    for (int i = 0; i < 8; i++) {                                              \
        if (per_lane_kmers[i] < last_lane_limit)                               \
            last_lane_limit = per_lane_kmers[i];                               \
    }                                                                          \
                                                                               \
    __m256i cur_data = _mm256_setzero_si256();                                 \
    size_t buf_pos = 16;                                                       \
    size_t seq_pos = 0;                                                        \
                                                                               \
    /* === MAIN LOOP (identical to single-read macro) === */                   \
    for (size_t iter = 0; iter < max_iter; iter++) {                           \
        if (buf_pos >= 16) {                                                   \
            size_t byte_offset = seq_pos / 4;                                  \
            __m256i byte_off_vec = _mm256_set1_epi32((int32_t)byte_offset);    \
            __m256i bv_idx = _mm256_add_epi32(gather_base_vec, byte_off_vec);  \
            cur_data = _mm256_i32gather_epi32((const int*)packed, bv_idx, 1);  \
            buf_pos = seq_pos % 16;                                            \
        }                                                                      \
        __m256i add_base = _mm256_and_si256(                                   \
            _mm256_srli_epi32(cur_data, (int)(buf_pos * 2)), mask_2bit);       \
        buf_pos++;                                                             \
        seq_pos++;                                                             \
                                                                               \
        __m256i remove_base = delay_buf[rd_idx];                               \
        delay_buf[wr_idx] = add_base;                                          \
        wr_idx = (wr_idx + 1) & delay_mask;                                    \
        if (iter >= S - 1) rd_idx = (rd_idx + 1) & delay_mask;                \
                                                                               \
        __m256i fw_rotated = csyncmer_simd_rotl7(fw);                          \
        __m256i add_hash_fw = csyncmer_simd_lookup(f_table, add_base);         \
        __m256i fw_hash = _mm256_xor_si256(fw_rotated, add_hash_fw);           \
                                                                               \
        __m256i hash_out;                                                      \
        uint8_t strand_mask = 0;                                               \
                                                                               \
        if (CANONICAL) {                                                       \
            __m256i rc_hash;                                                   \
            if (iter < S - 1) {                                                \
                fw = fw_hash;                                                  \
                prev_remove_base = remove_base;                                \
                continue;                                                      \
            } else if (iter == S - 1) {                                        \
                rc = _mm256_setzero_si256();                                   \
                for (size_t j = 0; j < S; j++) {                               \
                    size_t bi = (S - 1 - j) & delay_mask;                     \
                    __m256i bj = delay_buf[bi];                                \
                    __m256i rv = csyncmer_simd_lookup(rc_table_v, bj);         \
                    rc = _mm256_xor_si256(csyncmer_simd_rotl7(rc), rv);       \
                }                                                              \
                rc_hash = rc;                                                  \
                __m256i fw_rm = csyncmer_simd_lookup(f_rot_table, remove_base);\
                fw = _mm256_xor_si256(fw_hash, fw_rm);                        \
                prev_remove_base = remove_base;                                \
            } else {                                                           \
                __m256i rc_rotr = csyncmer_simd_rotr7(rc);                     \
                __m256i rc_rm = csyncmer_simd_lookup(                          \
                    rc_rotr7_table, prev_remove_base);                         \
                __m256i rc_ad = csyncmer_simd_lookup(c_rot_table, add_base);   \
                rc_hash = _mm256_xor_si256(                                    \
                    _mm256_xor_si256(rc_rotr, rc_rm), rc_ad);                 \
                __m256i fw_rm = csyncmer_simd_lookup(f_rot_table, remove_base);\
                fw = _mm256_xor_si256(fw_hash, fw_rm);                        \
                rc = rc_hash;                                                  \
                prev_remove_base = remove_base;                                \
            }                                                                  \
            hash_out = _mm256_min_epu32(fw_hash, rc_hash);                     \
            if (COLLECT_POSITIONS) {                                           \
                __m256i cmp = _mm256_cmpgt_epi32(                              \
                    _mm256_xor_si256(fw_hash,                                  \
                        _mm256_set1_epi32((int)0x80000000)),                   \
                    _mm256_xor_si256(rc_hash,                                  \
                        _mm256_set1_epi32((int)0x80000000)));                  \
                strand_mask = (uint8_t)_mm256_movemask_ps(                    \
                    _mm256_castsi256_ps(cmp));                                \
            }                                                                  \
        } else {                                                               \
            if (iter >= S - 1) {                                               \
                __m256i rm = csyncmer_simd_lookup(f_rot_table, remove_base);   \
                fw = _mm256_xor_si256(fw_hash, rm);                           \
            } else {                                                           \
                fw = fw_hash;                                                  \
                continue;                                                      \
            }                                                                  \
            hash_out = fw_hash;                                                \
        }                                                                      \
                                                                               \
        __m256i pos_vec = _mm256_set1_epi32((int)pos);                         \
        __m256i elem = _mm256_or_si256(                                        \
            _mm256_and_si256(hash_out, val_mask),                              \
            _mm256_and_si256(pos_vec, pos_mask));                              \
        ring_buf[ring_idx] = elem;                                             \
                                                                               \
        if (CANONICAL && COLLECT_POSITIONS) {                                  \
            strand_ring[ring_idx] = strand_mask;                               \
            __m256i pcmp = _mm256_cmpgt_epi32(                                \
                _mm256_xor_si256(prefix_min,                                   \
                    _mm256_set1_epi32((int)0x80000000)),                       \
                _mm256_xor_si256(elem,                                         \
                    _mm256_set1_epi32((int)0x80000000)));                      \
            uint8_t umask = (uint8_t)_mm256_movemask_ps(                      \
                _mm256_castsi256_ps(pcmp));                                    \
            prefix_min = _mm256_min_epu32(prefix_min, elem);                   \
            prefix_strand = (prefix_strand & ~umask) |                        \
                            (strand_mask & umask);                             \
        } else {                                                               \
            prefix_min = _mm256_min_epu32(prefix_min, elem);                   \
        }                                                                      \
                                                                               \
        if (pos == max_pos_val) {                                              \
            uint32_t delta = (1 << 16) - 2 - 2 * window_size;                 \
            pos -= delta;                                                      \
            pos_offset += delta;                                               \
            pos_offset_vec = _mm256_set1_epi32((int)pos_offset);              \
            __m256i dv = _mm256_set1_epi32((int)delta);                       \
            prefix_min = _mm256_sub_epi32(prefix_min, dv);                    \
            for (size_t j = 0; j < window_size; j++)                           \
                ring_buf[j] = _mm256_sub_epi32(ring_buf[j], dv);             \
        }                                                                      \
                                                                               \
        pos++;                                                                 \
        ring_idx++;                                                            \
                                                                               \
        if (ring_idx == window_size) {                                         \
            ring_idx = 0;                                                      \
            __m256i smin = ring_buf[window_size - 1];                          \
            if (CANONICAL && COLLECT_POSITIONS) {                              \
                uint8_t sstr = strand_ring[window_size - 1];                   \
                for (size_t j = window_size - 1; j > 0; j--) {                \
                    __m256i prev = ring_buf[j - 1];                            \
                    __m256i sc = _mm256_cmpgt_epi32(                           \
                        _mm256_xor_si256(smin,                                 \
                            _mm256_set1_epi32((int)0x80000000)),               \
                        _mm256_xor_si256(prev,                                 \
                            _mm256_set1_epi32((int)0x80000000)));              \
                    uint8_t sm = (uint8_t)_mm256_movemask_ps(                 \
                        _mm256_castsi256_ps(sc));                              \
                    smin = _mm256_min_epu32(smin, prev);                       \
                    sstr = (sstr & ~sm) | (strand_ring[j - 1] & sm);          \
                    ring_buf[j - 1] = smin;                                    \
                    strand_ring[j - 1] = sstr;                                 \
                }                                                              \
            } else {                                                           \
                for (size_t j = window_size - 1; j > 0; j--) {                \
                    smin = _mm256_min_epu32(smin, ring_buf[j - 1]);            \
                    ring_buf[j - 1] = smin;                                    \
                }                                                              \
            }                                                                  \
            prefix_min = _mm256_set1_epi32((int)UINT32_MAX);                   \
            if (CANONICAL && COLLECT_POSITIONS) prefix_strand = 0;             \
        }                                                                      \
                                                                               \
        /* === SYNCMER CHECK === */                                            \
        if (iter < S - 1 + window_size - 1) continue;                         \
                                                                               \
        __m256i suf_min = ring_buf[ring_idx];                                  \
        __m256i min_elem;                                                      \
        uint8_t overall_strand = 0;                                            \
                                                                               \
        if (CANONICAL && COLLECT_POSITIONS) {                                  \
            uint8_t suf_str = strand_ring[ring_idx];                           \
            __m256i oc = _mm256_cmpgt_epi32(                                  \
                _mm256_xor_si256(prefix_min,                                   \
                    _mm256_set1_epi32((int)0x80000000)),                       \
                _mm256_xor_si256(suf_min,                                      \
                    _mm256_set1_epi32((int)0x80000000)));                      \
            uint8_t om = (uint8_t)_mm256_movemask_ps(                         \
                _mm256_castsi256_ps(oc));                                      \
            min_elem = _mm256_min_epu32(prefix_min, suf_min);                  \
            overall_strand = (prefix_strand & ~om) | (suf_str & om);          \
        } else {                                                               \
            min_elem = _mm256_min_epu32(prefix_min, suf_min);                  \
        }                                                                      \
                                                                               \
        __m256i min_pos_vec = _mm256_add_epi32(                                \
            _mm256_and_si256(min_elem, pos_mask), pos_offset_vec);             \
        size_t kmer_idx = iter - S + 1 - window_size + 1;                     \
                                                                               \
        if (COLLECT_POSITIONS && CANONICAL) {                                  \
            /* Per-lane: kmer_idx is the position within each read */          \
            __m256i fs = _mm256_set1_epi32((int)kmer_idx);                    \
            __m256i ls = _mm256_set1_epi32(                                   \
                (int)(kmer_idx + window_size - 1));                            \
            __m256i isy = _mm256_or_si256(                                    \
                _mm256_cmpeq_epi32(min_pos_vec, fs),                          \
                _mm256_cmpeq_epi32(min_pos_vec, ls));                         \
            /* Mask out lanes where kmer_idx >= their num_kmers */            \
            if (kmer_idx >= last_lane_limit) {                                 \
                __m256i lims = _mm256_load_si256(                              \
                    (const __m256i*)chunk_ends);                               \
                isy = _mm256_and_si256(isy,                                   \
                    _mm256_cmpgt_epi32(lims, fs));                            \
            }                                                                  \
            int km = _mm256_movemask_ps(_mm256_castsi256_ps(isy));            \
            while (km) {                                                       \
                int lane = __builtin_ctz(km);                                  \
                size_t lc = lane_counts[lane];                                 \
                if (lc < max_per_read) {                                       \
                    out_positions[lane][lc] = (uint32_t)kmer_idx;              \
                    out_strands[lane][lc] = (overall_strand >> lane) & 1;      \
                }                                                              \
                lane_counts[lane] = lc + 1;                                    \
                km &= km - 1;                                                  \
            }                                                                  \
        } else if (COLLECT_POSITIONS && !CANONICAL) {                          \
            /* Non-canonical multi-read not needed */                          \
            (void)batch_pos; (void)batch_count; (void)batch_base_kmer;        \
            (void)idx_offsets; (void)w_minus_1_vec;                            \
        } else {                                                               \
            __m256i fs = _mm256_set1_epi32((int)kmer_idx);                    \
            __m256i ls = _mm256_set1_epi32(                                   \
                (int)(kmer_idx + window_size - 1));                            \
            __m256i isy = _mm256_or_si256(                                    \
                _mm256_cmpeq_epi32(min_pos_vec, fs),                          \
                _mm256_cmpeq_epi32(min_pos_vec, ls));                         \
            if (kmer_idx >= last_lane_limit) {                                 \
                __m256i lims = _mm256_load_si256(                              \
                    (const __m256i*)chunk_ends);                               \
                isy = _mm256_and_si256(isy,                                   \
                    _mm256_cmpgt_epi32(lims, fs));                            \
            }                                                                  \
            syncmer_count += __builtin_popcount(                               \
                _mm256_movemask_ps(_mm256_castsi256_ps(isy)));                \
        }                                                                      \
    } /* end main loop */                                                      \
                                                                               \
    /* Output: already written directly to out_positions/out_strands */        \
    if (COLLECT_POSITIONS) {                                                   \
        for (int i = 0; i < 8; i++) {                                          \
            size_t n = lane_counts[i] < max_per_read                           \
                ? lane_counts[i] : max_per_read;                               \
            out_counts[i] = n;                                                 \
        }                                                                      \
    }                                                                          \
                                                                               \
    if (own_buf) free(work_buf);                                               \
}

/* Generate multi-read canonical+positions variant */
CSYNCMER_DEFINE_TWOSTACK_SIMD_32_MULTI(csyncmer_twostack_simd_32_multi_canonical_positions_, 1, 1)

/* Public API wrapper */

static inline void csyncmer_twostack_simd_32_multi_canonical_positions(
    const char* sequences[8],
    const size_t lengths[8],
    size_t K,
    size_t S,
    uint32_t* out_positions[8],
    uint8_t*  out_strands[8],
    size_t max_per_read,
    size_t out_counts[8],
    uint8_t* work_buf,
    size_t work_buf_size
) {
    csyncmer_twostack_simd_32_multi_canonical_positions_(
        sequences, lengths, K, S, out_positions, out_strands,
        max_per_read, out_counts, work_buf, work_buf_size);
}

// Compute work buffer size needed for multi-read SIMD with given parameters.
static inline size_t csyncmer_multi_work_buf_size(size_t max_read_len, size_t K, size_t S) {
    if (S == 0 || S >= K) return 0;
    size_t window_size = K - S + 1;
    size_t ps_uniform = ((max_read_len + 3) / 4 + 32 + 31) & ~(size_t)31;
    size_t total_packed = 8 * ps_uniform;
    size_t delay_size = 1;
    while (delay_size < S) delay_size *= 2;
    size_t off_delay = (total_packed + 32 + 31) & ~(size_t)31;
    size_t off_ring = off_delay + ((delay_size * 32 + 31) & ~(size_t)31);
    size_t off_strand = off_ring + ((window_size * 32 + 31) & ~(size_t)31);
    return off_strand + ((window_size + 31) & ~(size_t)31);
}

#define CSYNCMER_DEFINE_TWOSTACK_SIMD_32_MULTI_TWOPASS(FUNC_NAME, CANONICAL, COLLECT_POSITIONS) \
static inline void FUNC_NAME(                                                  \
    const char* sequences[8],                                                  \
    const size_t lengths[8],                                                   \
    size_t K,                                                                  \
    size_t S,                                                                  \
    uint32_t* out_positions[8],                                                \
    uint8_t*  out_strands[8],                                                  \
    size_t max_per_read,                                                       \
    size_t out_counts[8],                                                      \
    uint8_t* work_buf,                                                         \
    size_t work_buf_size                                                       \
) {                                                                            \
    for (int i = 0; i < 8; i++) out_counts[i] = 0;                            \
    if (S == 0 || S >= K) return;                                              \
                                                                               \
    size_t window_size = K - S + 1;                                            \
                                                                               \
    /* Compute per-lane num_kmers and max across lanes */                      \
    size_t per_lane_kmers[8];                                                  \
    size_t max_num_kmers = 0;                                                  \
    int any_valid = 0;                                                         \
    for (int i = 0; i < 8; i++) {                                              \
        if (sequences[i] && lengths[i] >= K) {                                 \
            per_lane_kmers[i] = lengths[i] - K + 1;                            \
            if (per_lane_kmers[i] > max_num_kmers)                             \
                max_num_kmers = per_lane_kmers[i];                             \
            any_valid = 1;                                                     \
        } else {                                                               \
            per_lane_kmers[i] = 0;                                             \
        }                                                                      \
    }                                                                          \
    if (!any_valid) return;                                                     \
                                                                               \
    /* Fall back to single-read scalar for very short reads */                 \
    if (max_num_kmers < 64) {                                                  \
        for (int i = 0; i < 8; i++) {                                          \
            if (per_lane_kmers[i] > 0) {                                       \
                if (CANONICAL && COLLECT_POSITIONS) {                           \
                    out_counts[i] = csyncmer_canonical_rescan_32_(             \
                        sequences[i], lengths[i], K, S,                        \
                        out_positions[i], out_strands[i], max_per_read);       \
                }                                                              \
            }                                                                  \
        }                                                                      \
        return;                                                                \
    }                                                                          \
                                                                               \
    /* ---- Layout computation ---- */                                         \
    size_t max_len = 0;                                                        \
    for (int i = 0; i < 8; i++)                                                \
        if (lengths[i] > max_len) max_len = lengths[i];                        \
    size_t max_num_smers = max_len - S + 1;                                    \
    size_t ps_uniform = ((max_len + 3) / 4 + 32 + 31) & ~(size_t)31;          \
    size_t total_packed = 8 * ps_uniform;                                      \
    size_t n_groups = (max_len + 15) / 16;                                     \
    size_t interleaved_size = (n_groups * 32 + 31) & ~(size_t)31;              \
                                                                               \
    /* Work buffer layout:                                                     \
     *   temp_packed | interleaved | all_bases | hash_buf | strand_buf |       \
     *   ring_buf    | strand_ring                                             \
     * temp_packed is only needed during packing, then reusable. */            \
    size_t off_packed = 0;                                                     \
    size_t off_interleaved = (total_packed + 32 + 31) & ~(size_t)31;           \
    size_t off_bases = off_interleaved + interleaved_size;                     \
    size_t off_hash = off_bases +                                              \
        ((max_len * sizeof(__m256i) + 31) & ~(size_t)31);                     \
    size_t off_sbuf = off_hash +                                               \
        ((max_num_smers * sizeof(__m256i) + 31) & ~(size_t)31);               \
    size_t off_ring = off_sbuf + ((max_num_smers + 31) & ~(size_t)31);        \
    size_t off_strand = off_ring +                                             \
        ((window_size * sizeof(__m256i) + 31) & ~(size_t)31);                 \
    size_t needed = off_strand + ((window_size + 31) & ~(size_t)31);          \
                                                                               \
    int own_buf = 0;                                                           \
    if (!work_buf || work_buf_size < needed) {                                 \
        work_buf = (uint8_t*)aligned_alloc(32, needed);                        \
        if (!work_buf) return;                                                 \
        own_buf = 1;                                                           \
    }                                                                          \
    uint8_t *temp_packed  = work_buf + off_packed;                             \
    uint8_t *interleaved  = work_buf + off_interleaved;                        \
    __m256i *all_bases    = (__m256i*)(work_buf + off_bases);                   \
    __m256i *hash_buf     = (__m256i*)(work_buf + off_hash);                   \
    uint8_t *strand_buf   = work_buf + off_sbuf;                               \
    __m256i *ring_buf     = (__m256i*)(work_buf + off_ring);                   \
    uint8_t *strand_ring  = work_buf + off_strand;                             \
                                                                               \
    /* ---- Phase 0: Pack & interleave ---- */                                 \
    memset(temp_packed, 0, total_packed + 32);                                 \
    for (int i = 0; i < 8; i++) {                                              \
        if (per_lane_kmers[i] > 0)                                             \
            csyncmer_pack_seq_2bit(sequences[i], lengths[i],                   \
                                   temp_packed + i * ps_uniform);              \
    }                                                                          \
    csyncmer_interleave_packed(temp_packed, ps_uniform,                        \
                               interleaved, n_groups);                         \
                                                                               \
    /* ---- SIMD hash tables ---- */                                           \
    __m256i f_table = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_F32);    \
    uint32_t rot = ((S - 1) * 7) & 31;                                        \
    alignas(32) uint32_t f_rot_arr[8];                                         \
    for (int i = 0; i < 8; i++) {                                              \
        uint32_t x = CSYNCMER_SIMD_F32[i];                                    \
        f_rot_arr[i] = (x << rot) | (x >> (32 - rot));                        \
    }                                                                          \
    __m256i f_rot_table = _mm256_load_si256((const __m256i*)f_rot_arr);        \
                                                                               \
    __m256i rc_table_v = _mm256_setzero_si256();                               \
    __m256i rc_rotr7_table = _mm256_setzero_si256();                           \
    __m256i c_rot_table = _mm256_setzero_si256();                              \
    if (CANONICAL) {                                                           \
        rc_table_v = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_RC32);    \
        rc_rotr7_table = _mm256_load_si256(                                    \
            (const __m256i*)CSYNCMER_SIMD_RC32_ROTR7);                        \
        alignas(32) uint32_t c_rot_arr[8];                                     \
        for (int i = 0; i < 8; i++) {                                          \
            uint32_t cx = CSYNCMER_SIMD_RC32[i];                               \
            c_rot_arr[i] = (cx << rot) | (cx >> (32 - rot));                   \
        }                                                                      \
        c_rot_table = _mm256_load_si256((const __m256i*)c_rot_arr);            \
    }                                                                          \
                                                                               \
                                                                               \
    /* ================================================================ */     \
    /* PASS 1: Hash computation (chunked pre-extraction)                */     \
    /* Extract bases in L1d-sized chunks, tight loop per chunk          */     \
    /* Output: hash_buf[] (sequential write -> L2 streaming)            */     \
    /* ================================================================ */     \
                                                                               \
    /* Extract initial S bases for first s-mer */                              \
    csyncmer_extract_bases(interleaved, 0, S, all_bases);                      \
                                                                               \
    /* Compute initial forward hash for first s-mer (bases 0..S-1) */          \
    __m256i fw = _mm256_setzero_si256();                                       \
    for (size_t j = 0; j < S; j++)                                             \
        fw = _mm256_xor_si256(csyncmer_simd_rotl7(fw),                         \
                              csyncmer_simd_lookup(f_table, all_bases[j]));    \
    __m256i fw_hash_0 = fw;                                                    \
                                                                               \
    /* Compute initial RC hash for first s-mer */                              \
    __m256i rc = _mm256_setzero_si256();                                       \
    if (CANONICAL) {                                                           \
        for (size_t j = 0; j < S; j++)                                         \
            rc = _mm256_xor_si256(csyncmer_simd_rotl7(rc),                     \
                csyncmer_simd_lookup(rc_table_v, all_bases[S - 1 - j]));       \
    }                                                                          \
                                                                               \
    /* Store first s-mer hash */                                               \
    if (CANONICAL) {                                                           \
        hash_buf[0] = _mm256_min_epu32(fw_hash_0, rc);                         \
        if (COLLECT_POSITIONS) {                                               \
            __m256i is_fw = _mm256_cmpeq_epi32(                                \
                hash_buf[0], fw_hash_0);                                       \
            strand_buf[0] = (uint8_t)_mm256_movemask_ps(                       \
                _mm256_castsi256_ps(_mm256_andnot_si256(                       \
                    is_fw, _mm256_set1_epi32(-1))));                           \
        }                                                                      \
    } else {                                                                   \
        hash_buf[0] = fw_hash_0;                                               \
        if (COLLECT_POSITIONS) strand_buf[0] = 0;                              \
    }                                                                          \
                                                                               \
    /* Chunked rolling loop: extract + hash in L1d-sized blocks */             \
    {                                                                          \
        fw = _mm256_xor_si256(fw_hash_0,                                       \
            csyncmer_simd_lookup(f_rot_table, all_bases[0]));                  \
        __m256i prev_rm = all_bases[0];                                        \
        const size_t chunk_sz = 512;                                           \
                                                                               \
        for (size_t cs = 1; cs < max_num_smers; cs += chunk_sz) {              \
            size_t ce = cs + chunk_sz;                                         \
            if (ce > max_num_smers) ce = max_num_smers;                        \
            size_t n_bases = (ce - cs) + S - 1;                                \
            csyncmer_extract_bases(interleaved, cs, n_bases, all_bases);       \
                                                                               \
            for (size_t ci = 0; ci < ce - cs; ci++) {                          \
                size_t si = cs + ci;                                           \
                __m256i add_base = all_bases[ci + S - 1];                      \
                __m256i remove_base = all_bases[ci];                           \
                                                                               \
                __m256i fw_hash = _mm256_xor_si256(                            \
                    csyncmer_simd_rotl7(fw),                                   \
                    csyncmer_simd_lookup(f_table, add_base));                  \
                fw = _mm256_xor_si256(fw_hash,                                 \
                    csyncmer_simd_lookup(f_rot_table, remove_base));           \
                                                                               \
                if (CANONICAL) {                                               \
                    __m256i rc_hash = _mm256_xor_si256(                        \
                        _mm256_xor_si256(csyncmer_simd_rotr7(rc),              \
                            csyncmer_simd_lookup(rc_rotr7_table, prev_rm)),    \
                        csyncmer_simd_lookup(c_rot_table, add_base));          \
                    rc = rc_hash;                                              \
                    prev_rm = remove_base;                                     \
                                                                               \
                    hash_buf[si] = _mm256_min_epu32(fw_hash, rc_hash);         \
                    if (COLLECT_POSITIONS) {                                   \
                        __m256i is_fw = _mm256_cmpeq_epi32(                    \
                            hash_buf[si], fw_hash);                            \
                        strand_buf[si] = (uint8_t)_mm256_movemask_ps(          \
                            _mm256_castsi256_ps(_mm256_andnot_si256(           \
                                is_fw, _mm256_set1_epi32(-1))));               \
                    }                                                          \
                } else {                                                       \
                    hash_buf[si] = fw_hash;                                    \
                    if (COLLECT_POSITIONS) strand_buf[si] = 0;                 \
                }                                                              \
            }                                                                  \
        } /* end chunked hash loop */                                          \
    }                                                                          \
                                                                               \
    /* ================================================================ */     \
    /* PASS 2: Twostack sliding window minimum                          */     \
    /* Working set: ring_buf (32 KB) in L1d, hash_buf reads from L2     */     \
    /* ================================================================ */     \
                                                                               \
    alignas(32) uint32_t chunk_ends[8];                                        \
    for (int i = 0; i < 8; i++)                                                \
        chunk_ends[i] = (uint32_t)per_lane_kmers[i];                           \
                                                                               \
    size_t last_lane_limit = max_num_kmers;                                    \
    for (int i = 0; i < 8; i++) {                                              \
        if (per_lane_kmers[i] < last_lane_limit)                               \
            last_lane_limit = per_lane_kmers[i];                               \
    }                                                                          \
                                                                               \
    /* Initialize ring buffer and twostack state */                            \
    for (size_t i = 0; i < window_size; i++)                                   \
        ring_buf[i] = _mm256_set1_epi32((int)UINT32_MAX);                     \
    __m256i prefix_min = _mm256_set1_epi32((int)UINT32_MAX);                   \
    __m256i val_mask = _mm256_set1_epi32((int)0xFFFF0000);                     \
    __m256i pos_mask = _mm256_set1_epi32(0x0000FFFF);                          \
    uint8_t prefix_strand = 0;                                                 \
    if (CANONICAL && COLLECT_POSITIONS) {                                       \
        for (size_t i = 0; i < window_size; i++) strand_ring[i] = 0;          \
    }                                                                          \
                                                                               \
    size_t ring_idx = 0;                                                       \
    uint32_t pos = 0;                                                          \
    uint32_t pos_offset = 0;                                                   \
    __m256i pos_offset_vec = _mm256_setzero_si256();                           \
    const uint32_t max_pos_val = 0xFFFF;                                       \
                                                                               \
    size_t lane_counts[8] = {0};                                               \
    size_t syncmer_count = 0;                                                  \
    (void)syncmer_count;                                                       \
                                                                               \
    /* Twostack loop: iterate over s-mer hashes */                             \
    for (size_t si = 0; si < max_num_smers; si++) {                            \
        __m256i hash_out = hash_buf[si];                                       \
        uint8_t strand_mask = 0;                                               \
        if (CANONICAL && COLLECT_POSITIONS)                                     \
            strand_mask = strand_buf[si];                                      \
                                                                               \
        /* Pack hash|pos into ring element */                                  \
        __m256i pos_vec = _mm256_set1_epi32((int)pos);                         \
        __m256i elem = _mm256_or_si256(                                        \
            _mm256_and_si256(hash_out, val_mask),                              \
            _mm256_and_si256(pos_vec, pos_mask));                              \
        ring_buf[ring_idx] = elem;                                             \
                                                                               \
        if (CANONICAL && COLLECT_POSITIONS) {                                  \
            strand_ring[ring_idx] = strand_mask;                               \
            /* Update prefix_min with strand tracking (opt #4) */             \
            __m256i old_pm = prefix_min;                                        \
            prefix_min = _mm256_min_epu32(prefix_min, elem);                   \
            __m256i changed = _mm256_andnot_si256(                             \
                _mm256_cmpeq_epi32(prefix_min, old_pm),                        \
                _mm256_set1_epi32(-1));                                        \
            uint8_t umask = (uint8_t)_mm256_movemask_ps(                      \
                _mm256_castsi256_ps(changed));                                \
            prefix_strand = (prefix_strand & ~umask) |                        \
                            (strand_mask & umask);                             \
        } else {                                                               \
            prefix_min = _mm256_min_epu32(prefix_min, elem);                   \
        }                                                                      \
                                                                               \
        /* Position overflow handling */                                       \
        if (pos == max_pos_val) {                                              \
            uint32_t delta = (1 << 16) - 2 - 2 * window_size;                 \
            pos -= delta;                                                      \
            pos_offset += delta;                                               \
            pos_offset_vec = _mm256_set1_epi32((int)pos_offset);              \
            __m256i dv = _mm256_set1_epi32((int)delta);                       \
            prefix_min = _mm256_sub_epi32(prefix_min, dv);                    \
            for (size_t j = 0; j < window_size; j++)                           \
                ring_buf[j] = _mm256_sub_epi32(ring_buf[j], dv);             \
        }                                                                      \
                                                                               \
        pos++;                                                                 \
        ring_idx++;                                                            \
                                                                               \
        /* Rescan when ring wraps */                                           \
        if (ring_idx == window_size) {                                         \
            ring_idx = 0;                                                      \
            __m256i smin = ring_buf[window_size - 1];                          \
            if (CANONICAL && COLLECT_POSITIONS) {                              \
                uint8_t sstr = strand_ring[window_size - 1];                   \
                for (size_t j = window_size - 1; j > 0; j--) {                \
                    __m256i prev = ring_buf[j - 1];                            \
                    __m256i old_smin = smin;                                    \
                    smin = _mm256_min_epu32(smin, prev);                       \
                    __m256i ch = _mm256_andnot_si256(                          \
                        _mm256_cmpeq_epi32(smin, old_smin),                    \
                        _mm256_set1_epi32(-1));                                \
                    uint8_t sm = (uint8_t)_mm256_movemask_ps(                 \
                        _mm256_castsi256_ps(ch));                              \
                    sstr = (sstr & ~sm) | (strand_ring[j - 1] & sm);          \
                    ring_buf[j - 1] = smin;                                    \
                    strand_ring[j - 1] = sstr;                                 \
                }                                                              \
            } else {                                                           \
                for (size_t j = window_size - 1; j > 0; j--) {                \
                    smin = _mm256_min_epu32(smin, ring_buf[j - 1]);            \
                    ring_buf[j - 1] = smin;                                    \
                }                                                              \
            }                                                                  \
            prefix_min = _mm256_set1_epi32((int)UINT32_MAX);                   \
            if (CANONICAL && COLLECT_POSITIONS) prefix_strand = 0;             \
        }                                                                      \
                                                                               \
        /* === SYNCMER CHECK === */                                            \
        if (si < window_size - 1) continue;                                    \
                                                                               \
        __m256i suf_min = ring_buf[ring_idx];                                  \
        __m256i min_elem;                                                      \
        uint8_t overall_strand = 0;                                            \
                                                                               \
        if (CANONICAL && COLLECT_POSITIONS) {                                  \
            uint8_t suf_str = strand_ring[ring_idx];                           \
            __m256i old_pm = prefix_min;                                        \
            min_elem = _mm256_min_epu32(prefix_min, suf_min);                  \
            __m256i ch = _mm256_andnot_si256(                                  \
                _mm256_cmpeq_epi32(min_elem, old_pm),                          \
                _mm256_set1_epi32(-1));                                        \
            uint8_t om = (uint8_t)_mm256_movemask_ps(                         \
                _mm256_castsi256_ps(ch));                                      \
            overall_strand = (prefix_strand & ~om) | (suf_str & om);          \
        } else {                                                               \
            min_elem = _mm256_min_epu32(prefix_min, suf_min);                  \
        }                                                                      \
                                                                               \
        __m256i min_pos_vec = _mm256_add_epi32(                                \
            _mm256_and_si256(min_elem, pos_mask), pos_offset_vec);             \
        size_t kmer_idx = si - window_size + 1;                                \
                                                                               \
        if (COLLECT_POSITIONS && CANONICAL) {                                  \
            __m256i fs = _mm256_set1_epi32((int)kmer_idx);                    \
            __m256i ls = _mm256_set1_epi32(                                   \
                (int)(kmer_idx + window_size - 1));                            \
            __m256i isy = _mm256_or_si256(                                    \
                _mm256_cmpeq_epi32(min_pos_vec, fs),                          \
                _mm256_cmpeq_epi32(min_pos_vec, ls));                         \
            if (kmer_idx >= last_lane_limit) {                                 \
                __m256i lims = _mm256_load_si256(                              \
                    (const __m256i*)chunk_ends);                               \
                isy = _mm256_and_si256(isy,                                   \
                    _mm256_cmpgt_epi32(lims, fs));                            \
            }                                                                  \
            int km = _mm256_movemask_ps(_mm256_castsi256_ps(isy));            \
            while (km) {                                                       \
                int lane = __builtin_ctz(km);                                  \
                size_t lc = lane_counts[lane];                                 \
                if (lc < max_per_read) {                                       \
                    out_positions[lane][lc] = (uint32_t)kmer_idx;              \
                    out_strands[lane][lc] = (overall_strand >> lane) & 1;      \
                }                                                              \
                lane_counts[lane] = lc + 1;                                    \
                km &= km - 1;                                                  \
            }                                                                  \
        } else if (!COLLECT_POSITIONS) {                                       \
            __m256i fs = _mm256_set1_epi32((int)kmer_idx);                    \
            __m256i ls = _mm256_set1_epi32(                                   \
                (int)(kmer_idx + window_size - 1));                            \
            __m256i isy = _mm256_or_si256(                                    \
                _mm256_cmpeq_epi32(min_pos_vec, fs),                          \
                _mm256_cmpeq_epi32(min_pos_vec, ls));                         \
            if (kmer_idx >= last_lane_limit) {                                 \
                __m256i lims = _mm256_load_si256(                              \
                    (const __m256i*)chunk_ends);                               \
                isy = _mm256_and_si256(isy,                                   \
                    _mm256_cmpgt_epi32(lims, fs));                            \
            }                                                                  \
            syncmer_count += __builtin_popcount(                               \
                _mm256_movemask_ps(_mm256_castsi256_ps(isy)));                \
        }                                                                      \
    } /* end pass 2 loop */                                                    \
                                                                               \
    if (COLLECT_POSITIONS) {                                                   \
        for (int i = 0; i < 8; i++) {                                          \
            size_t n = lane_counts[i] < max_per_read                           \
                ? lane_counts[i] : max_per_read;                               \
            out_counts[i] = n;                                                 \
        }                                                                      \
    }                                                                          \
                                                                               \
    if (own_buf) free(work_buf);                                               \
}

/* Generate two-pass multi-read canonical+positions variant */
CSYNCMER_DEFINE_TWOSTACK_SIMD_32_MULTI_TWOPASS(csyncmer_twostack_simd_32_multi_canonical_positions_twopass_, 1, 1)

/* Public API wrapper for two-pass variant */
static inline void csyncmer_twostack_simd_32_multi_canonical_positions_twopass(
    const char* sequences[8],
    const size_t lengths[8],
    size_t K,
    size_t S,
    uint32_t* out_positions[8],
    uint8_t*  out_strands[8],
    size_t max_per_read,
    size_t out_counts[8],
    uint8_t* work_buf,
    size_t work_buf_size
) {
    csyncmer_twostack_simd_32_multi_canonical_positions_twopass_(
        sequences, lengths, K, S, out_positions, out_strands,
        max_per_read, out_counts, work_buf, work_buf_size);
}

// Compute work buffer size for two-pass multi-read SIMD.
static inline size_t csyncmer_multi_work_buf_size_twopass(size_t max_read_len, size_t K, size_t S) {
    if (S == 0 || S >= K) return 0;
    size_t window_size = K - S + 1;
    size_t max_num_smers = max_read_len - S + 1;
    size_t ps_uniform = ((max_read_len + 3) / 4 + 32 + 31) & ~(size_t)31;
    size_t total_packed = 8 * ps_uniform;
    size_t n_groups = (max_read_len + 15) / 16;
    size_t interleaved_size = (n_groups * 32 + 31) & ~(size_t)31;
    size_t off_interleaved = (total_packed + 32 + 31) & ~(size_t)31;
    size_t off_bases = off_interleaved + interleaved_size;
    size_t off_hash = off_bases + ((max_read_len * 32 + 31) & ~(size_t)31);
    size_t off_sbuf = off_hash + ((max_num_smers * 32 + 31) & ~(size_t)31);
    size_t off_ring = off_sbuf + ((max_num_smers + 31) & ~(size_t)31);
    size_t off_strand = off_ring + ((window_size * 32 + 31) & ~(size_t)31);
    return off_strand + ((window_size + 31) & ~(size_t)31);
}

#endif  // __AVX2__

#endif  // CSYNCMER_FASTQ_TWOSTACK_H
