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

    // 1. TWOSTACK count only (~1400 MB/s with AVX2, ~400 MB/s scalar fallback)
    size_t count1 = csyncmer_twostack_simd_32_count(seq, len, K, S);
    printf("TWOSTACK count:     %zu syncmers\n", count1);

    // 2. TWOSTACK with positions
    size_t max_positions = len - K + 1;
    uint32_t* positions = (uint32_t*)malloc(max_positions * sizeof(uint32_t));
    size_t count2 = csyncmer_twostack_simd_32_positions(seq, len, K, S, positions, max_positions);
    printf("TWOSTACK positions: %zu syncmers\n", count2);
    if (count2 > 0) {
        printf("  First 3 positions: %u, %u, %u\n",
               positions[0], positions[1], positions[2]);
    }

    // 3. Iterator API (scalar, portable, exact results, ~350 MB/s)
    CsyncmerIterator64* iter = csyncmer_iterator_create_64(seq, len, K, S);
    size_t count3 = 0;
    size_t pos;
    while (csyncmer_iterator_next_64(iter, &pos)) {
        count3++;
    }
    csyncmer_iterator_destroy_64(iter);
    printf("Iterator 64-bit:    %zu syncmers (scalar, portable, exact)\n", count3);

    printf("\nNotes:\n");
    printf("- 32-bit and 64-bit counts differ due to different hash tie-breaking\n");
#ifdef __AVX2__
    printf("- TWOSTACK uses AVX2 SIMD (~1400 MB/s, 16-bit hash approximation)\n");
#else
    printf("- TWOSTACK falls back to scalar RESCAN (~400 MB/s, no AVX2 on this platform)\n");
#endif
    printf("- Iterator: scalar, portable, exact 64-bit hashes (~350 MB/s)\n");

    free(positions);
    return 0;
}
