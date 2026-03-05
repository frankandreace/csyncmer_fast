#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "csyncmer_fast.h"

int main() {
    // Test sequence (random DNA)
    const char* seq =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
        "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
    size_t len = strlen(seq);
    size_t K = 31, S = 15;

    printf("Sequence length: %zu\n", len);
    printf("K=%zu, S=%zu\n\n", K, S);

    size_t max_positions = len - K + 1;
    uint32_t* positions = (uint32_t*)malloc(max_positions * sizeof(uint32_t));

    // 1. AVX2 SIMD twostack (~1400 MB/s, requires -mavx2 or -march=native)
#ifdef __AVX2__
    size_t count1 = csyncmer_twostack_simd_32_count(seq, len, K, S);
    printf("TWOSTACK count:     %zu syncmers (AVX2 SIMD)\n", count1);

    size_t count2 = csyncmer_twostack_simd_32_positions(seq, len, K, S, positions, max_positions);
    printf("TWOSTACK positions: %zu syncmers\n", count2);
    if (count2 > 0)
        printf("  First 3 positions: %u, %u, %u\n", positions[0], positions[1], positions[2]);
#else
    // 2. Scalar rescan (~400 MB/s, portable, exact 32-bit hashes)
    size_t count1 = csyncmer_rescan_32_count(seq, len, K, S, NULL);
    printf("RESCAN count:       %zu syncmers (scalar)\n", count1);

    size_t count2 = csyncmer_rescan_32_positions(seq, len, K, S, positions, max_positions);
    printf("RESCAN positions:   %zu syncmers\n", count2);
    if (count2 > 0)
        printf("  First 3 positions: %u, %u, %u\n", positions[0], positions[1], positions[2]);
#endif

    // 3. Streaming iterator (scalar, portable, exact 64-bit hashes, ~350 MB/s)
    CsyncmerIterator64* iter = csyncmer_iterator_create_64(seq, len, K, S);
    size_t count3 = 0;
    size_t pos;
    while (csyncmer_iterator_next_64(iter, &pos)) {
        count3++;
    }
    csyncmer_iterator_destroy_64(iter);
    printf("Iterator 64-bit:    %zu syncmers (scalar, portable, exact)\n", count3);

    printf("\nNotes:\n");
    printf("- 32-bit and 64-bit counts may differ due to different hash tie-breaking\n");

    free(positions);
    return 0;
}
