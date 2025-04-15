// #include "circular_array.h"
// #include "hashing.h"
#include "csyncmer_fast_iterator.h"

typedef struct
{
    U64 minimum ;
    size_t minimum_position ;
    size_t absolute_kmer_position ; 
} Syncmer;


void compute_closed_syncmers_rescan_iterator(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer) {
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

void compute_closed_syncmers_branchless(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer){
    // if len < K, return
    if (length < K){
        printf("SEQUENCE SMALLER THAN K.\n") ;
        return ;
    }

    U64 seed  = 7 ;
    size_t num_s_mers = length - S + 1 ;
    size_t num_k_mers = length - K + 1 ;
    size_t window_size = (size_t)K - (size_t)S + 1 ;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;
    // initialize seqhashiterator
    // CircularArray *ca = circularArrayCreate(window_size) ;
    CircularArray *ca = circularArrayCreate( window_size) ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;
    // Syncmer *sync = (Syncmer *)malloc(sizeof(Syncmer)) ;
    // if (sync == NULL){
    //     printf("sync is null.\n") ;
    //     exit (-1) ;
    // }
    if (ca == NULL){
        printf("ca is null.\n") ;
        exit (-1) ;
    }
    if (sh == NULL){
        printf("sh is null.\n") ;
        exit (-1) ;
    }
    if (si == NULL){
        printf("si is null.\n") ;
        exit (-1) ;
    }

    U64 current_smer = 0 ;
    size_t current_position = 0 ;
    size_t smer_position = 0;

    // 1 - precompute 1st window
    size_t precompute = 0 ;
    while(precompute < window_size - 1){
        bool x = seqhashNext(si, &current_smer);
        // printf("x is %s, !x is %s\n",  x?"true":"false", !x?"true":"false") ;
        if(!x){
            // printf("seqhashnext is false\n") ;
            exit(1) ;
        }
        // printf("PC-HASHED %llu\n", current_smer) ;
        circularInsertBranchless(ca,current_smer);
        precompute++;
        computed_smers++;
    }

    // 2 - run until end of given sequence
    while(seqhashNext(si, &current_smer)){
        // printf("HASHED %llu\n", current_smer) ;
        computed_smers++ ;
        circularInsertBranchless(ca, current_smer);
        if (is_syncmer(ca, &smer_position)){
            computed_syncmers++ ;
            // printf("SYNCMER AT POS %lu; S-MER AT %lu\n", current_position, current_position + smer_position);
        }
        current_position++ ;
    }
    free(ca) ;
    free(si) ;
    // qfree(sync) ;
    printf("[RESCAN_CIRCULAR_ARRAY_BRANCHLESS]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ; 
    printf("[RESCAN_CIRCULAR_ARRAY_BRANCHLESS]:: HASHED %lu S-MERS\n", computed_smers) ;
    *num_syncmer = computed_syncmers;
}

void compute_closed_syncmers(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer){
    // if len < K, return
    if (length < K){
        printf("SEQUENCE SMALLER THAN K.\n") ;
        return ;
    }

    U64 seed  = 7 ;
    size_t num_s_mers = length - S + 1 ;
    size_t num_k_mers = length - K + 1 ;
    size_t window_size = (size_t)K - (size_t)S + 1 ;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;
    // initialize seqhashiterator
    // CircularArray *ca = circularArrayCreate(window_size) ;
    CircularArray *ca = circularArrayCreate( window_size) ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;
    // Syncmer *sync = (Syncmer *)malloc(sizeof(Syncmer)) ;
    // if (sync == NULL){
    //     printf("sync is null.\n") ;
    //     exit (-1) ;
    // }
    if (ca == NULL){
        printf("ca is null.\n") ;
        exit (-1) ;
    }
    if (sh == NULL){
        printf("sh is null.\n") ;
        exit (-1) ;
    }
    if (si == NULL){
        printf("si is null.\n") ;
        exit (-1) ;
    }

    U64 current_smer = 0 ;
    size_t current_position = 0 ;
    size_t smer_position = 0;

    // 1 - precompute 1st window
    size_t precompute = 0 ;
    while(precompute < window_size - 1){
        bool x = seqhashNext(si, &current_smer);
        // printf("x is %s, !x is %s\n",  x?"true":"false", !x?"true":"false") ;
        if(!x){
            // printf("seqhashnext is false\n") ;
            exit(1) ;
        }
        // printf("PC-HASHED %llu\n", current_smer) ;
        circularInsert(ca,current_smer);
        precompute++;
        computed_smers++;
    }

    // 2 - run until end of given sequence
    while(seqhashNext(si, &current_smer)){
        // printf("HASHED %llu\n", current_smer) ;
        computed_smers++ ;
        circularInsert(ca, current_smer);
        if (is_syncmer(ca, &smer_position)){
            computed_syncmers++ ;
            // printf("SYNCMER AT POS %lu; S-MER AT %lu\n", current_position, current_position + smer_position);
        }
        current_position++ ;
    }
    free(ca) ;
    free(si) ;
    // qfree(sync) ;
    printf("[RESCAN_CIRCULAR_ARRAY]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ; 
    printf("[RESCAN_CIRCULAR_ARRAY]:: HASHED %lu S-MERS\n", computed_smers) ;
    *num_syncmer = computed_syncmers;
}

void compute_closed_syncmers_naive(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer) {
    if(length < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    U64 num_results = 0;

    size_t num_s_mers = length - S + 1;
    size_t num_k_mers = length - K + 1;
    U64 *s_mer_hashes = (U64 *)malloc(num_s_mers * sizeof(U64));
    size_t window_size = K - S + 1;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;

    Seqhash *sh = seqhashCreate(S, window_size, seed);
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length);
    // printf("iterator is %s", si == NULL ? "null" : "ok") ;

    U64 smer;
    size_t s_pos = 0;
    bool hashnext;

    // HASH ALL S-MERS
    // printf("starting hashing\n") ;
    // hashnext = seqhashNext(si, &smer) ;
    // printf ("seqhashnext is %s\n", hashnext ? "true" : "false") ;
    while(seqhashNext(si, &smer)){
        // printf("HASH COMPUTED IS %llu\n", smer) ;
        s_mer_hashes[s_pos++] = smer;
        computed_smers++;
        // hashnext = seqhashNext(si, &smer) ;
        // printf("%llu\t%lu\n", smer, s_pos) ;
    }
    // printf("HASHED %lu s-mers\n", computed_smers) ;
    if (computed_smers == 0){ exit(-1) ;}

    size_t front = 0, back = 0;
    int min_pos;
    U64 min_smer;

    // USE ARRAY SCAN TO COMPUTE SYNCMERS IN O(N*K)

    for(size_t i = 0; i < num_k_mers; i++){
        min_smer = U64MAX;
        min_pos = -1;
        // printf("[") ;
        for (size_t j = 0; j < window_size; j++){
            // printf("%llu ", s_mer_hashes[j+i]);
            if (s_mer_hashes[j+i] < min_smer){
                min_pos = j;
                min_smer = s_mer_hashes[j+i];
            }

        }
        // printf("]\t");
        // printf("MIN IS %llu, pos is %lu\n", min_smer, min_pos) ;
        if (min_pos == 0 || min_pos == window_size - 1){
            // printf("%llu\t%lu\n", min_smer, i+ min_pos);
            // printf("SYNCMER. MIN IS %llu, pos is %lu\n", min_smer, min_pos) ;
            computed_syncmers++;
        }
    }

    // releasing memory
    free(s_mer_hashes) ;
    *num_syncmer = computed_syncmers;
    printf("[NAIVE]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ; 
    printf("[NAIVE]:: HASHED %lu S-MERS\n", computed_smers) ; 
}

void hahsing_speed_benchmark(char *sequence_input, size_t length, size_t K, size_t S){
    //
    U64 seed  = 7;
    size_t num_s_mers = length - S + 1;
    size_t num_k_mers = length - K + 1;
    size_t window_size = (size_t)K - (size_t)S + 1;
    size_t computed_smers = 0 ;

    // initialize seqhashiterator
    Seqhash *sh = seqhashCreate(S, window_size, seed);
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;
    // CircularArray *ca = circularArrayCreate(si->sh->w) ;

    // Syncmer *sync = (Syncmer *)malloc(sizeof(Syncmer)) ;
    U64 current_smer;
    size_t current_position;
    
    // printf("%llu\t%lu\n", current_smer, current_position) ;
    while(seqhashNext(si, &current_smer)){
        // printf("%llu\t%lu\n", current_smer, current_position) ;
        computed_smers++;
        current_position++;
    }
    printf("[HASHING_BENCHMARK]:: HASHED %lu S-MERS\n", computed_smers) ; 

}


void compute_closed_syncmers_rescan(char *sequence_input, size_t sequence_length, size_t K, size_t S, size_t *num_syncmer){
    #define ARRAYSIZE 21
    // if len < K, return
    if (sequence_length < K){
        printf("SEQUENCE SMALLER THAN K.\n") ;
        return ;
    }

    if(K - S + 1 > (1 << ARRAYSIZE)){
        printf("WINDOW SIZE > BUFFER LENGTH.\n") ;
        return ;
    }

    // initialization

    U64 seed  = 7 ;
    size_t num_s_mers = sequence_length - S + 1 ;
    size_t num_k_mers = sequence_length - K + 1 ;
    size_t window_size = K - S + 1 ;
    size_t array_size = (1 << ARRAYSIZE);
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, sequence_length) ;

    U64 *hashvector = (U64 *)malloc(sizeof(U64) * array_size) ;

    // checking pointers are fine
    if (sh == NULL){
        printf("sh is null.\n") ;
        exit (-1) ;
    }
    if (si == NULL){
        printf("si is null.\n") ;
        exit (-1) ;
    }

    // syncmer computation values

    U64 current_smer = 0 ;

    size_t current_position = 0;

    U64 minimum = U64MAX;
    size_t minimum_position;
    size_t absolute_kmer_position = 0;

    bool first_loop = true ;

    // precompute first w-1 elements
    while(current_position < window_size - 1){

        seqhashNext(si, &current_smer);
        hashvector[current_position] = current_smer;
        
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

            seqhashNext(si, &current_smer);
            hashvector[current_position] = current_smer;

            current_position++;
            computed_smers++;
        }

        // 2 - find syncmers
        // scan the rest of the array
        for(size_t i = window_size - 1; i < len_scan; i++){
            
            // rescan last w-1 elements if minimum out of context
            if(minimum_position < i - window_size + 1){
                size_t curr_pos;
                minimum = U64MAX ;
                for(int j = 1; j < window_size; j++){
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

    free(si) ;

    printf("[RESCAN_LARGE_ARRAY]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ; 
    printf("[RESCAN_LARGE_ARRAY]:: HASHED %lu S-MERS\n", computed_smers) ;
    *num_syncmer = computed_syncmers;
}

void compute_closed_syncmers_deque_rayan(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer) {
    if(length < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    U64 seed  = 7 ;
    size_t num_s_mers = length - S + 1 ;
    size_t num_k_mers = length - K + 1 ;
    size_t window_size = K - S + 1 ;

    U64 current_smer = 0 ;
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;

    // Precompute all s-mer hashes
    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;
    U64 *s_mer_hashes = (U64 *)malloc(sizeof(U64) * num_s_mers) ;

    for(size_t s_mer_pos = 0; s_mer_pos < num_s_mers; s_mer_pos++) {
        seqhashNext(si, &current_smer);
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

    free(s_mer_hashes);
    free(deque);
}

// static inline void update_minimum_branchless(U64 candidate, U64 *current_minimum, size_t candidate_position, size_t *current_minimum_position)
// {
//     // U64 smaller = (candidate < *current_minimum);
//     U64 mask = -(candidate < *current_minimum);
//     *current_minimum = (*current_minimum & ~mask) | (candidate & mask);
//     *current_minimum_position  = (*current_minimum_position & ~mask) | (candidate_position & mask);
// }


// void compute_closed_syncmers_rescan_branchless(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer){
//     // #define CIRCULARARRAYSIZE 27
//     // if len < K, return
//     if (length < K){
//         printf("SEQUENCE SMALLER THAN K.\n") ;
//         return ;
//     }

//     U64 seed  = 7 ;
//     size_t num_s_mers = length - S + 1 ;
//     size_t num_k_mers = length - K + 1 ;
//     size_t window_size = K - S + 1 ;
//     size_t circular_array_size = window_size;

//     size_t computed_syncmers = 0 ;
//     size_t computed_smers = 0 ;

//     U64 minimum_value = U64MAX;

//     // initialize 

//     Seqhash *sh = seqhashCreate(S, window_size, seed) ;
//     SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;

//     U64 *hashvector = (U64 *)malloc(sizeof(U64) * circular_array_size) ;

//     if (sh == NULL){
//         printf("sh is null.\n") ;
//         exit (-1) ;
//     }
//     if (si == NULL){
//         printf("si is null.\n") ;
//         exit (-1) ;
//     }

//     U64 current_smer = 0 ;

//     // filling first w-1 elements
//     for (size_t precompute = 0;  precompute < window_size - 1; precompute++){
//         seqhashNext(si, &current_smer) ;
//         hashvector[precompute] = current_smer;
//         computed_smers++ ;
//     }


//     // 2 - find syncmers
//     size_t minimum_position = 0;
//     size_t absolute_kmer_position = 0;

//     //fiund minimum first w-1
//     for(size_t i =0; i < window_size - 1; i++){
//         update_minimum_branchless(hashvector[i], &minimum_value, i, &minimum_position);
//     }

//     // scan the rest
//     size_t leftInWindow;
//     size_t idx = window_size - 2 ; 
//     for(size_t i = window_size - 1; i < num_s_mers; i++){
//         idx = (idx + 1) % window_size;
//         seqhashNext(si, &current_smer) ;

//         hashvector[idx] = current_smer;
//          // rescan if minimum out of context

//         leftInWindow = i + 1 - window_size; // earliest valid pos

//         if(minimum_position < leftInWindow){
//             size_t curr_pos;
//             minimum_value = U64MAX;
//             for(size_t j = 1; j < window_size; j++)
//             {
//                 curr_pos = (i - window_size + j)%window_size;
//                 update_minimum_branchless(hashvector[curr_pos], &minimum_value, (i - window_size + j), &minimum_position);
//             }
                
//         }
//         //verify minimum
//         update_minimum_branchless(hashvector[idx], &minimum_value, i, &minimum_position);

//         // minimum_position = get_pos(minimum_value_position);
//         if (minimum_position == i || minimum_position == i - window_size + 1U) {
//             computed_syncmers++;
//         }

//         absolute_kmer_position++;
//     }
//     free(si) ;

//     *num_syncmer = computed_syncmers;
// }