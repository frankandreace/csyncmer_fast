#ifndef CSYNCMER_FASTQ_SPLIT_H
#define CSYNCMER_FASTQ_SPLIT_H

// Split two-pass syncmer detection:
//   Pass 1: csyncmer_hash_only_multi (with out_hash_buf to retain buffer)
//   Pass 2: twostack sliding-window minimum (no strand tracking)
//
// Usage:
//   __m256i *hb;
//   size_t plk[8];
//   size_t mns = csyncmer_hash_only_multi(seqs, lens, K, S, wb, ws, &hb, plk);
//   __m256i *ring = csyncmer_split_ring_buf(wb, mns, S);
//   csyncmer_twostack_only_multi(hb, mns, plk, K, S, ...);

#include "../../csyncmer_fast.h"
#include "csyncmer_fastq_hash.h"

#ifdef __AVX2__

// Compute ring_buf pointer from work buffer layout (interleaved | 2KB pad | hash_buf | 2KB pad | ring_buf).
static inline __m256i* csyncmer_split_ring_buf(
    uint8_t* work_buf, size_t max_num_smers, size_t S
) {
    size_t max_len = max_num_smers + S - 1;
    size_t n_groups = (max_len + 15) / 16;
    size_t interleaved_size = (n_groups * 32 + 31) & ~(size_t)31;
    size_t interleaved_hash_pad = 2048;  // break 4K aliasing between interleaved and hash_buf
    size_t hash_size = (max_num_smers * sizeof(__m256i) + 31) & ~(size_t)31;
    size_t hash_ring_pad = 2048;  // break 4K aliasing between hash_buf and ring_buf
    return (__m256i*)(work_buf + interleaved_size + interleaved_hash_pad + hash_size + hash_ring_pad);
}

// Twostack-only pass: sliding-window minimum on pre-computed hash_buf.
// No strand tracking — positions only.
//
// Optimizations:
//   1. pos_vec YMM with fs/ls derived at use sites (port scheduling)
//   2. Running ring pointer (eliminates index→address math)
//   3. Chunked overflow handling (removes branch from hot loop)
//   4. Software prefetch on hash_buf (masks L2 latency)
//   5. Pre-load suffix min before ring store (disambiguation avoidance)
//   6. Pre-load in rescan backward sweep (same technique)
static inline void csyncmer_twostack_only_multi(
    const __m256i* hash_buf,
    size_t max_num_smers,
    const size_t per_lane_kmers[8],
    size_t K,
    size_t S,
    uint32_t* out_positions[8],
    size_t max_per_read,
    size_t out_counts[8],
    __m256i* ring_buf       // pre-allocated, window_size + 1 elements
) {
    for (int i = 0; i < 8; i++) out_counts[i] = 0;

    size_t window_size = K - S + 1;

    // Find max_num_kmers and last_lane_limit
    size_t max_num_kmers = 0;
    size_t last_lane_limit = (size_t)-1;
    for (int i = 0; i < 8; i++) {
        if (per_lane_kmers[i] > max_num_kmers)
            max_num_kmers = per_lane_kmers[i];
        if (per_lane_kmers[i] > 0 && per_lane_kmers[i] < last_lane_limit)
            last_lane_limit = per_lane_kmers[i];
    }
    if (max_num_kmers == 0) return;
    if (last_lane_limit == (size_t)-1) last_lane_limit = max_num_kmers;

    alignas(32) uint32_t chunk_ends[8];
    for (int i = 0; i < 8; i++)
        chunk_ends[i] = (uint32_t)per_lane_kmers[i];

    // Initialize ring buffer (+1 padding element for pre-load)
    for (size_t i = 0; i <= window_size; i++)
        ring_buf[i] = _mm256_set1_epi32((int)UINT32_MAX);

    __m256i prefix_min = _mm256_set1_epi32((int)UINT32_MAX);
    const __m256i val_mask = _mm256_set1_epi32((int)0xFFFF0000);
    const __m256i pos_mask = _mm256_set1_epi32(0x0000FFFF);
    const __m256i ones = _mm256_set1_epi32(1);

    // [Opt 1] Maintained pos_vec YMM — no scalar→SIMD broadcasts in hot loop.
    // pos_vec: used BEFORE increment (value = si, tracks pos for element packing).
    // After increment, fs/ls are derived at use sites as pos_vec - w_vec / pos_vec - ones,
    // spreading p015 ops between p01-restricted vpcmpeqd for better port scheduling.
    // Comparisons are done in reduced pos space (no pos_offset addition needed).
    __m256i pos_vec = _mm256_setzero_si256();
    const __m256i w_vec = _mm256_set1_epi32((int)window_size);

    uint32_t pos_scalar = 0;     // shadow of pos_vec lane 0 for overflow math
    uint32_t pos_offset = 0;
    __m256i pos_offset_vec = _mm256_setzero_si256();
    const uint32_t max_pos_val = 0xFFFF;

    // [Opt 2] Running ring pointer — replaces ring_idx * 32 + base addressing.
    __m256i *ring_ptr = ring_buf;
    __m256i * const ring_end = ring_buf + window_size;

    const size_t warmup_end = window_size - 1;
    size_t lane_counts[8] = {0};

    // [Opt 3] Chunked outer loop — overflow check removed from inner loop.
    // Each chunk runs at most (max_pos_val - pos_scalar + 1) iterations,
    // guaranteeing pos fits in 16 bits throughout.
    size_t si = 0;
    while (si < max_num_smers) {
        size_t safe_iters = (size_t)(max_pos_val - pos_scalar + 1);
        size_t chunk_end = si + safe_iters;
        if (chunk_end > max_num_smers) chunk_end = max_num_smers;
        size_t chunk_start = si;

        for (; si < chunk_end; si++) {
            // [Opt 4] Software prefetch — hash_buf is ~230 KB (L2 resident).
            _mm_prefetch((const char*)&hash_buf[si + 16], _MM_HINT_T0);

            __m256i hash_out = hash_buf[si];

            // [Opt 5] Pre-load suffix min BEFORE storing to ring.
            // Avoids store→load memory disambiguation machine clears:
            // without this, the CPU sees store(A) then load(A+32) on the same
            // cache line and may speculatively assume aliasing.
            // +1 padding element is allocated so this is safe even at ring_end-1.
            __m256i suf_min_preload = *(ring_ptr + 1);

            // Pack hash|pos into ring element
            __m256i elem = _mm256_or_si256(
                _mm256_and_si256(hash_out, val_mask),
                _mm256_and_si256(pos_vec, pos_mask));

            *ring_ptr = elem;

            // Update prefix_min
            prefix_min = _mm256_min_epu32(prefix_min, elem);

            // Advance ring pointer and pos vector
            ring_ptr++;
            pos_vec = _mm256_add_epi32(pos_vec, ones);

            // [Opt 6] Ring wrap + suffix rescan with pre-load pattern.
            // Pre-loads ring_buf[j-2] before storing ring_buf[j-1] to break
            // store→load disambiguation on shared cache lines.
            if (ring_ptr == ring_end) {
                ring_ptr = ring_buf;
                __m256i smin = ring_buf[window_size - 1];
                __m256i next_elem = ring_buf[window_size - 2];
                for (size_t j = window_size - 1; j > 1; j--) {
                    smin = _mm256_min_epu32(smin, next_elem);
                    next_elem = ring_buf[j - 2];
                    ring_buf[j - 1] = smin;
                }
                smin = _mm256_min_epu32(smin, next_elem);
                ring_buf[0] = smin;
                suf_min_preload = *ring_ptr;
                prefix_min = _mm256_set1_epi32((int)UINT32_MAX);
            }

            // Skip warmup phase
            if (si < warmup_end) continue;

            // Syncmer check: compare in reduced pos space (no pos_offset needed).
            // fs/ls derived from pos_vec at use site — interleaves p015 vpsubd
            // between p01-restricted vpcmpeqd for better port scheduling.
            __m256i min_elem = _mm256_min_epu32(prefix_min, suf_min_preload);
            __m256i min_pos_raw = _mm256_and_si256(min_elem, pos_mask);

            __m256i fs_reduced = _mm256_sub_epi32(pos_vec, w_vec);
            __m256i isy = _mm256_cmpeq_epi32(min_pos_raw, fs_reduced);
            __m256i ls_reduced = _mm256_sub_epi32(pos_vec, ones);
            isy = _mm256_or_si256(isy,
                _mm256_cmpeq_epi32(min_pos_raw, ls_reduced));

            size_t kmer_idx = si - warmup_end;
            if (kmer_idx >= last_lane_limit) {
                __m256i lims = _mm256_load_si256(
                    (const __m256i*)chunk_ends);
                // Recover absolute kmer index for bounds check (rare path)
                __m256i fs_abs = _mm256_add_epi32(fs_reduced, pos_offset_vec);
                isy = _mm256_and_si256(isy,
                    _mm256_cmpgt_epi32(lims, fs_abs));
            }

            int km = _mm256_movemask_ps(_mm256_castsi256_ps(isy));
            while (km) {
                int lane = __builtin_ctz(km);
                size_t lc = lane_counts[lane];
                if (lc < max_per_read)
                    out_positions[lane][lc] = (uint32_t)kmer_idx;
                lane_counts[lane] = lc + 1;
                km &= km - 1;
            }
        }

        // Update scalar pos tracker from iteration count
        pos_scalar += (uint32_t)(si - chunk_start);

        // Overflow adjustment (fires at most once per ~64K iterations)
        if (si < max_num_smers) {
            uint32_t delta = (1 << 16) - 2 - 2 * (uint32_t)window_size;
            __m256i dv = _mm256_set1_epi32((int)delta);
            pos_scalar -= delta;
            pos_vec = _mm256_sub_epi32(pos_vec, dv);
            // pos_offset compensates for pos field reduction; used only
            // in the rare bounds-check path to recover absolute kmer index.
            pos_offset += delta;
            pos_offset_vec = _mm256_set1_epi32((int)pos_offset);
            prefix_min = _mm256_sub_epi32(prefix_min, dv);
            for (size_t j = 0; j < window_size; j++)
                ring_buf[j] = _mm256_sub_epi32(ring_buf[j], dv);
        }
    }

    for (int i = 0; i < 8; i++) {
        size_t n = lane_counts[i] < max_per_read
            ? lane_counts[i] : max_per_read;
        out_counts[i] = n;
    }
}

// Compute work buffer size for split two-pass.
// Layout: interleaved | 2KB pad | hash_buf | 2KB pad | ring_buf
static inline size_t csyncmer_split_work_buf_size(
    size_t max_read_len, size_t K, size_t S
) {
    if (S == 0 || S >= K) return 0;
    size_t window_size = K - S + 1;
    size_t max_num_smers = max_read_len - S + 1;
    size_t n_groups = (max_read_len + 15) / 16;
    size_t interleaved_size = (n_groups * 32 + 31) & ~(size_t)31;
    size_t interleaved_hash_pad = 2048;  // break 4K aliasing between interleaved and hash_buf
    size_t hash_size = (max_num_smers * sizeof(__m256i) + 31) & ~(size_t)31;
    size_t hash_ring_pad = 2048;  // break 4K aliasing between hash_buf and ring_buf
    size_t ring_size = ((window_size + 1) * sizeof(__m256i) + 31) & ~(size_t)31;  // +1 pre-load padding
    return interleaved_size + interleaved_hash_pad + hash_size + hash_ring_pad + ring_size;
}

#endif  // __AVX2__
#endif  // CSYNCMER_FASTQ_SPLIT_H
