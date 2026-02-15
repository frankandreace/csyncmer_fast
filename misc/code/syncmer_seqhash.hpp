#ifndef SYNCMER_SEQHASH_H
#define SYNCMER_SEQHASH_H

// Syncmer implementations using SeqHash (64-bit hash)
// Superseded by ntHash32 implementations in csyncmer_fast.h

#include "legacy_infrastructure.hpp"  // For SeqHash, CircularArray, etc.

// ============================================================================
// SEQHASH-BASED SYNCMER IMPLEMENTATIONS
// ============================================================================

void csyncmer_seqhash_rescan_iterator(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer) {
    if(length < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    U64 smer_value ;
    size_t kmer_position ;
    size_t smer_position ;

    U64 computed_syncmer = 0;

    syncmerIterator *syncmer_iterator = syncmerIteratorInitialization(K, S, sequence_input, length) ;

    if (!syncmer_iterator->isDone){
        while(get_next_syncmer(syncmer_iterator, &kmer_position, &smer_position, &smer_value)){
            computed_syncmer++ ;
        }
    }
    *num_syncmer = computed_syncmer ;
    syncmerIteratorDestroy(syncmer_iterator);

    printf("[RESCAN_ITERATOR]:: COMPUTED %llu CLOSED SYNCMERS\n", computed_syncmer) ;
}

void csyncmer_seqhash_rescan_branchless(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer){
    if (length < K){
        printf("SEQUENCE SMALLER THAN K.\n") ;
        return ;
    }

    U64 seed  = 7 ;
    size_t window_size = (size_t)K - (size_t)S + 1 ;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;
    CircularArray *ca = circularArrayCreate( window_size) ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;
    if (ca == NULL){
        fprintf(stderr, "ca is null.\n") ;
        exit (-1) ;
    }
    if (sh == NULL){
        fprintf(stderr, "sh is null.\n") ;
        exit (-1) ;
    }
    if (si == NULL){
        fprintf(stderr, "si is null.\n") ;
        exit (-1) ;
    }

    U64 current_smer = 0 ;
    bool current_orientation;
    size_t current_position = 0 ;
    size_t smer_position = 0;

    // 1 - precompute 1st window
    size_t precompute = 0 ;
    while(precompute < window_size - 1){
        bool x = seqhashNext(si, &current_smer, &current_orientation);
        if(!x){
            exit(1) ;
        }
        circularInsertBranchless(ca,current_smer);
        precompute++;
        computed_smers++;
    }

    // 2 - run until end of given sequence
    while(seqhashNext(si, &current_smer, &current_orientation)){
        computed_smers++ ;
        circularInsertBranchless(ca, current_smer);
        if (is_syncmer(ca, &smer_position)){
            computed_syncmers++ ;
        }
        current_position++ ;
    }
    circularArrayDestroy(ca) ;
    free(sh) ;
    free(si) ;
    printf("[RESCAN_CIRCULAR_ARRAY_BRANCHLESS]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ;
    printf("[RESCAN_CIRCULAR_ARRAY_BRANCHLESS]:: HASHED %lu S-MERS\n", computed_smers) ;
    *num_syncmer = computed_syncmers;
}

void csyncmer_seqhash_rescan_circular(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer){
    if (length < K){
        printf("SEQUENCE SMALLER THAN K.\n") ;
        return ;
    }

    U64 seed  = 7 ;
    size_t window_size = (size_t)K - (size_t)S + 1 ;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;
    CircularArray *ca = circularArrayCreate( window_size) ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;
    if (ca == NULL){
        fprintf(stderr, "ca is null.\n") ;
        exit (-1) ;
    }
    if (sh == NULL){
        fprintf(stderr, "sh is null.\n") ;
        exit (-1) ;
    }
    if (si == NULL){
        fprintf(stderr, "si is null.\n") ;
        exit (-1) ;
    }

    U64 current_smer = 0 ;
    bool current_orientation;
    size_t current_position = 0 ;
    size_t smer_position = 0;

    // 1 - precompute 1st window
    size_t precompute = 0 ;
    while(precompute < window_size - 1){
        bool x = seqhashNext(si, &current_smer, &current_orientation);
        if(!x){
            exit(1) ;
        }
        circularInsert(ca,current_smer);
        precompute++;
        computed_smers++;
    }

    // 2 - run until end of given sequence
    while(seqhashNext(si, &current_smer, &current_orientation)){
        computed_smers++ ;
        circularInsert(ca, current_smer);
        if (is_syncmer(ca, &smer_position)){
            computed_syncmers++ ;
        }
        current_position++ ;
    }
    printf("[RESCAN_CIRCULAR_ARRAY]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ;
    printf("[RESCAN_CIRCULAR_ARRAY]:: HASHED %lu S-MERS\n", computed_smers) ;
    printf("[RESCAN_CIRCULAR_ARRAY]:: RESCANS %lu (%.2f%% of s-mers)\n", ca->rescan_count, 100.0 * ca->rescan_count / computed_smers) ;
    printf("[RESCAN_CIRCULAR_ARRAY]:: CONSECUTIVE RESCANS %lu (%.2f%% of rescans)\n", ca->consecutive_rescan_count, 100.0 * ca->consecutive_rescan_count / ca->rescan_count) ;
    free(ca) ;
    free(sh) ;
    free(si) ;
    *num_syncmer = computed_syncmers;
}

void csyncmer_seqhash_naive(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer) {
    if(length < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    U64 seed  = 7;

    size_t num_s_mers = length - S + 1;
    size_t num_k_mers = length - K + 1;
    U64 *s_mer_hashes = (U64 *)malloc(num_s_mers * sizeof(U64));
    size_t window_size = K - S + 1;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;

    Seqhash *sh = seqhashCreate(S, window_size, seed);
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length);

    U64 smer;
    bool current_orientation;
    size_t s_pos = 0;

    // HASH ALL S-MERS
    while(seqhashNext(si, &smer, &current_orientation)){
        s_mer_hashes[s_pos++] = smer;
        computed_smers++;
    }
    if (computed_smers == 0){ exit(-1) ;}

    size_t min_pos;
    U64 min_smer;

    // USE ARRAY SCAN TO COMPUTE SYNCMERS IN O(N*K)

    for(size_t i = 0; i < num_k_mers; i++){
        min_smer = U64MAX;
        min_pos = -1;
        for (size_t j = 0; j < window_size; j++){
            if (s_mer_hashes[j+i] < min_smer){
                min_pos = j;
                min_smer = s_mer_hashes[j+i];
            }

        }
        if (min_pos == 0 || min_pos == window_size - 1){
            computed_syncmers++;
        }
    }

    // releasing memory
    free(s_mer_hashes) ;
    free(sh);
    free(si);
    *num_syncmer = computed_syncmers;
    printf("[NAIVE]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ;
    printf("[NAIVE]:: HASHED %lu S-MERS\n", computed_smers) ;
}

void csyncmer_seqhash_rescan_array(char *sequence_input, size_t sequence_length, size_t K, size_t S, size_t *num_syncmer){
    #define ARRAYSIZE 17
    if (sequence_length < K){
        printf("SEQUENCE SMALLER THAN K.\n") ;
        return ;
    }

    if(K - S + 1 > (1 << ARRAYSIZE)){
        printf("WINDOW SIZE > BUFFER LENGTH.\n") ;
        return ;
    }

    U64 seed  = 7 ;
    size_t num_s_mers = sequence_length - S + 1 ;
    size_t window_size = K - S + 1 ;
    size_t array_size = (1 << ARRAYSIZE);
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;
    size_t rescan_count = 0 ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, sequence_length) ;

    U64 *hashvector = (U64 *)malloc(sizeof(U64) * array_size) ;
    bool *orientationvector = (bool *)malloc(sizeof(bool) * array_size) ;

    if (sh == NULL){
        fprintf(stderr, "sh is null.\n") ;
        exit (-1) ;
    }
    if (si == NULL){
        fprintf(stderr, "si is null.\n") ;
        exit (-1) ;
    }

    U64 current_smer = 0 ;
    bool current_orientation;

    size_t current_position = 0;

    U64 minimum = U64MAX;
    size_t minimum_position = 0;
    size_t absolute_kmer_position = 0;

    bool first_loop = true ;

    // precompute first w-1 elements
    while(current_position < window_size - 1){

        seqhashNext(si, &current_smer, &current_orientation);
        hashvector[current_position] = current_smer;
        orientationvector[current_position] = current_orientation;

        current_position++;
        computed_smers++;
    }

    size_t len_scan = (num_s_mers < array_size) ? num_s_mers : array_size;

    while (num_s_mers > 0) {

        //scan the first w-1 to update minimum
        for(size_t i = 0; i < window_size - 1; i++){
            if (hashvector[i] < minimum){
                minimum = hashvector[i] ;
                minimum_position = i;
            }
        }

        // finishing computing hashes in the vector
        while(current_position < len_scan){

            seqhashNext(si, &current_smer, &current_orientation);
            hashvector[current_position] = current_smer;
            orientationvector[current_position] = current_orientation;

            current_position++;
            computed_smers++;
        }

        // 2 - find syncmers
        // scan the rest of the array
        for(size_t i = window_size - 1; i < len_scan; i++){

            // rescan last w-1 elements if minimum out of context
            if(minimum_position < i - window_size + 1){
                rescan_count++;
                size_t curr_pos;
                minimum = U64MAX ;
                for(size_t j = 1; j < window_size; j++){
                    curr_pos = i - window_size + j;
                    if ( hashvector[curr_pos] < minimum ){
                        minimum = hashvector[curr_pos] ;
                        minimum_position = curr_pos;
                    }
                }
            }

            //update minimum
            if (hashvector[i] < minimum){
                minimum = hashvector[i] ;
                minimum_position = i;
            }

            //check syncmer condition
            if(i == minimum_position || minimum_position == i - window_size + 1){
                computed_syncmers++;
            }

            //advcance k-mer condition
            absolute_kmer_position++;
        }

        //if I need to loop again, transfer the last w-1 elements at the beginnig to continue the scan
        if (num_s_mers > array_size){
            for(size_t i = 0; i < window_size - 1; i++){
                hashvector[i] = hashvector[array_size - window_size + 1 + i] ;
                orientationvector[i] = orientationvector[array_size - window_size + 1 + i] ;
            }
            num_s_mers = (first_loop) ? num_s_mers - len_scan: num_s_mers - len_scan + window_size - 1 ;
            len_scan = (num_s_mers + window_size - 1 < array_size) ? num_s_mers + window_size - 1 : array_size;
            first_loop = false;
            minimum = U64MAX;
            current_position = window_size - 1;
        }
        else{
            break;
        }
    }

    free(sh) ;
    free(si) ;
    free(hashvector) ;

    printf("[RESCAN_LARGE_ARRAY]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ;
    printf("[RESCAN_LARGE_ARRAY]:: HASHED %lu S-MERS\n", computed_smers) ;
    printf("[RESCAN_LARGE_ARRAY]:: RESCANS %lu (%.2f%% of s-mers)\n", rescan_count, 100.0 * rescan_count / computed_smers) ;
    *num_syncmer = computed_syncmers;
}

void csyncmer_seqhash_deque(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer) {
    if(length < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    U64 seed  = 7 ;
    size_t num_s_mers = length - S + 1 ;
    size_t window_size = K - S + 1 ;

    U64 current_smer = 0 ;
    bool current_orientation;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;

    // Precompute all s-mer hashes
    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;
    U64 *s_mer_hashes = (U64 *)malloc(sizeof(U64) * num_s_mers) ;

    for(size_t s_mer_pos = 0; s_mer_pos < num_s_mers; s_mer_pos++) {
        seqhashNext(si, &current_smer, &current_orientation);
        s_mer_hashes[s_mer_pos] = current_smer;
        computed_smers++;
    }

    // Initialize deque
    size_t *deque = (size_t *)malloc(num_s_mers * sizeof(size_t));
    size_t front = 0, back = 0;

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
    printf("[DEQUE]:: HASHED %lu S-MERS\n", computed_smers) ;
    *num_syncmer = computed_syncmers;

    free(sh);
    free(si);
    free(s_mer_hashes);
    free(deque);
}

// Note: compute_closed_syncmer_deque_nthash() moved to syncmer_nthash128.hpp

#endif // SYNCMER_SEQHASH_H
