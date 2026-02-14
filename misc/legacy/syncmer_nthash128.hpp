#ifndef SYNCMER_NTHASH128_H
#define SYNCMER_NTHASH128_H

// Syncmer implementations using ntHash 128-bit
// Superseded by ntHash32/64 implementations in csyncmer_fast.h

#include "legacy_infrastructure.hpp"  // For NtHash wrapper, U128, etc.

// ============================================================================
// NTHASH128-BASED SYNCMER IMPLEMENTATIONS
// ============================================================================

void compute_closed_syncmer_deque_nthash(const char *sequence_input, size_t length, size_t K, size_t S){
     NtHashHandle rolling_hash = nthash_create(sequence_input,strlen(sequence_input), S, 2);
    if (rolling_hash == NULL) {
        return;
    }

    size_t num_s_mers = length - S + 1 ;
    size_t window_size = K - S + 1 ;
    size_t computed_syncmers = 0;

    U128 *s_mer_hashes = (U128 *)malloc(sizeof(U128) * num_s_mers);
    for (size_t i = 0; i < num_s_mers; i++){
        nthash_roll(rolling_hash);
        s_mer_hashes[i] = nthash_get_canonical_hash_128(rolling_hash);
    }

    // Initialize deque
    U128 *deque = (U128 *)malloc(num_s_mers * sizeof(U128));
    U128 front = 0, back = 0;

    // Use deque to find minimal s-mers in O(N)
    for(size_t i = 0; i < num_s_mers; i++) {
        while(back > front  && s_mer_hashes[deque[back-1]] > s_mer_hashes[i]) {
            back--;
        }
        deque[back++] = i;
        if(i >= window_size && deque[front] <= i - window_size) {
                front++;
        }
        // Check for closed syncmer condition
        if(i >= window_size - 1) {
            size_t min_pos = deque[front];
            size_t kmer_pos = i - window_size + 1;
            if(min_pos == kmer_pos || min_pos == kmer_pos + K - S) {
                computed_syncmers++;
            }
        }
    }

    printf("[DEQUE]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ;
    printf("[DEQUE]:: HASHED %lu S-MERS\n", num_s_mers) ;

    free(s_mer_hashes);
    free(deque);
    nthash_destroy(rolling_hash);
}

#endif // SYNCMER_NTHASH128_H
