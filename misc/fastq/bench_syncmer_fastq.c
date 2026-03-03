/* bench_syncmer_fastq.c - standalone syncmer throughput benchmark
 * Measures raw csyncmer throughput on real FASTQ data,
 * excluding hash table and I/O overhead (ntHash packing is included).
 *
 * Build: cc -O3 -mavx2 -march=native -o bench_syncmer_fastq bench_syncmer_fastq.c -lz
 * Usage: ./bench_syncmer_fastq [-k 31] [-w 1022] [-single] ~/data/reads.fastq
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "csyncmer_fastq.h"

#define N_BUCKETS    256  /* enough for fine-grained bucket_shift values */

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    int k = 31, w = 1022, single_mode = 0, bucket_shift = 11;
    char *filename = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-k") && i+1 < argc) k = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-w") && i+1 < argc) w = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-b") && i+1 < argc) bucket_shift = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-single")) single_mode = 1;
        else filename = argv[i];
    }
    if (!filename) {
        fprintf(stderr, "Usage: %s [-k 31] [-w 1022] [-single] [-b bucket_shift] input.fastq\n", argv[0]);
        fprintf(stderr, "  -single  Use single-sequence API instead of multi (8-lane SIMD)\n");
        fprintf(stderr, "  -b N     Bucket shift for length grouping (default 11 = 2048bp)\n");
        fprintf(stderr, "  default: multi-8 with length bucketing (same as syng)\n");
        return 1;
    }

    int K = w + k - 1;
    int S = k;
    fprintf(stderr, "k=%d w=%d K=%d S=%d  mode=%s  bucket_shift=%d (%d bp)\n", k, w, K, S,
            single_mode ? "single" : "multi-8-bucketed",
            bucket_shift, 1 << bucket_shift);

    /* mmap the FASTQ file */
    int fd = open(filename, O_RDONLY);
    if (fd < 0) { perror("open"); return 1; }
    struct stat st;
    fstat(fd, &st);
    char *map = mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED) { perror("mmap"); return 1; }
    madvise(map, st.st_size, MADV_SEQUENTIAL);

    /* Phase 1: scan FASTQ to collect sequence pointers (no copy) */
    fprintf(stderr, "Scanning FASTQ to collect sequences...\n");
    size_t cap = 2000000, nseqs = 0;
    char **seqs = malloc(cap * sizeof(char*));
    size_t *lens = malloc(cap * sizeof(size_t));
    size_t total_bp = 0, max_len = 0, min_len = (size_t)-1;

    char *p = map, *end = map + st.st_size;
    while (p < end) {
        if (*p != '@') break;
        while (p < end && *p != '\n') p++;
        if (p >= end) break;
        p++;

        char *seq_start = p;
        while (p < end && *p != '\n') p++;
        size_t slen = p - seq_start;
        if (p >= end) break;
        p++;

        while (p < end && *p != '\n') p++;
        if (p >= end) break;
        p++;

        p += slen;
        if (p < end && *p == '\n') p++;

        if (nseqs >= cap) {
            cap *= 2;
            seqs = realloc(seqs, cap * sizeof(char*));
            lens = realloc(lens, cap * sizeof(size_t));
        }
        seqs[nseqs] = seq_start;
        lens[nseqs] = slen;
        total_bp += slen;
        if (slen > max_len) max_len = slen;
        if (slen < min_len) min_len = slen;
        nseqs++;
    }
    fprintf(stderr, "  %zu sequences, %.2f Gbp, avg %zu bp, min %zu, max %zu\n",
            nseqs, total_bp / 1e9, total_bp / nseqs, min_len, max_len);

    /* Allocate output buffers */
    size_t max_per_read = max_len;
    size_t total_syncmers = 0;
    double elapsed;

    if (single_mode) {
        /* Single-sequence API: one read at a time */
        uint32_t *pos_buf = malloc(max_per_read * sizeof(uint32_t));
        uint8_t *strand_buf = malloc(max_per_read * sizeof(uint8_t));

        fprintf(stderr, "Running syncmer detection (single-sequence)...\n");
        double t0 = now();

        for (size_t i = 0; i < nseqs; i++) {
            size_t count = csyncmer_twostack_simd_32_canonical_positions(
                seqs[i], lens[i], (size_t)K, (size_t)S,
                pos_buf, strand_buf, max_per_read);
            total_syncmers += count;
        }

        elapsed = now() - t0;
        free(pos_buf);
        free(strand_buf);
    } else {
        /* Multi API with length bucketing (mirrors syng's strategy) */
        size_t work_size = csyncmer_multi_work_buf_size(max_len, K, S);
        uint8_t *work_buf = aligned_alloc(32, work_size);

        uint32_t *pos_bufs[8];
        uint8_t *strand_bufs[8];
        for (int i = 0; i < 8; i++) {
            pos_bufs[i] = malloc(max_per_read * sizeof(uint32_t));
            strand_bufs[i] = malloc(max_per_read * sizeof(uint8_t));
        }

        /* Bucket state */
        typedef struct { char *seq[8]; size_t len[8]; int n; } Bucket;
        Bucket buckets[N_BUCKETS];
        for (int b = 0; b < N_BUCKETS; b++) buckets[b].n = 0;

        fprintf(stderr, "Running syncmer detection (multi-8 bucketed)...\n");
        double t0 = now();

        for (size_t i = 0; i < nseqs; i++) {
            if ((size_t)K > lens[i]) continue;

            int bi = (int)((lens[i] - K) >> bucket_shift);
            if (bi >= N_BUCKETS) bi = N_BUCKETS - 1;
            Bucket *bk = &buckets[bi];
            bk->seq[bk->n] = seqs[i];
            bk->len[bk->n] = lens[i];
            bk->n++;

            if (bk->n == 8) {
                const char *mseqs[8];
                size_t mlens[8], counts[8] = {0};
                for (int j = 0; j < 8; j++) { mseqs[j] = bk->seq[j]; mlens[j] = bk->len[j]; }

                csyncmer_twostack_simd_32_multi_canonical_positions(
                    mseqs, mlens, (size_t)K, (size_t)S,
                    pos_bufs, strand_bufs, max_per_read, counts,
                    work_buf, work_size);

                for (int j = 0; j < 8; j++) total_syncmers += counts[j];
                bk->n = 0;
            }
        }

        /* Flush partial buckets */
        for (int bi = 0; bi < N_BUCKETS; bi++) {
            Bucket *bk = &buckets[bi];
            if (bk->n == 0) continue;
            const char *mseqs[8] = {NULL};
            size_t mlens[8] = {0}, counts[8] = {0};
            for (int j = 0; j < bk->n; j++) { mseqs[j] = bk->seq[j]; mlens[j] = bk->len[j]; }

            csyncmer_twostack_simd_32_multi_canonical_positions(
                mseqs, mlens, (size_t)K, (size_t)S,
                pos_bufs, strand_bufs, max_per_read, counts,
                work_buf, work_size);

            for (int j = 0; j < bk->n; j++) total_syncmers += counts[j];
        }

        elapsed = now() - t0;
        for (int i = 0; i < 8; i++) { free(pos_bufs[i]); free(strand_bufs[i]); }
        free(work_buf);
    }

    printf("Sequences:  %zu\n", nseqs);
    printf("Total bp:   %.2f Gbp\n", total_bp / 1e9);
    printf("Avg length: %zu bp\n", total_bp / nseqs);
    printf("Syncmers:   %zu\n", total_syncmers);
    printf("Time:       %.3f s\n", elapsed);
    printf("Throughput: %.2f GB/s\n", total_bp / elapsed / 1e9);
    printf("Rate:       %.1f M syncmers/s\n", total_syncmers / elapsed / 1e6);

    free(seqs);
    free(lens);
    munmap(map, st.st_size);
    close(fd);
    return 0;
}
