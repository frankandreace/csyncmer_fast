#include "circular_array.h"
#include "circular_array1.h"
#include "hashing.h"

typedef struct
{
    U64 minimum ;
    size_t minimum_position ;
    size_t absolute_kmer_position ; 
} Syncmer;


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
    Syncmer *sync = (Syncmer *)malloc(sizeof(Syncmer)) ;
    if (sync == NULL){
        printf("sync is null.\n") ;
        exit (-1) ;
    }
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
    free(sync) ;
    printf("COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ; 
    printf("HASHED %lu S-MERS\n", computed_smers) ;
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
    printf("COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ; 
    printf("HASHED %lu S-MERS\n", computed_smers) ; 
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
    printf("HASHED %lu S-MERS\n", computed_smers) ; 

}


void compute_closed_syncmers_rescan(char *sequence_input, size_t length, size_t K, size_t S, size_t *num_syncmer){
    #define CIRCULARARRAYSIZE 27
    // if len < K, return
    if (length < K){
        printf("SEQUENCE SMALLER THAN K.\n") ;
        return ;
    }

    U64 seed  = 7 ;
    size_t num_s_mers = length - S + 1 ;
    size_t num_k_mers = length - K + 1 ;
    size_t window_size = (size_t)K - (size_t)S + 1 ;
    size_t circular_array_size = (1 << CIRCULARARRAYSIZE);
    size_t computed_syncmers = 0 ;
    size_t computed_smers = 0 ;

    // initialize 

    // CircularArray1 *ca = circularArrayCreate1(circular_array_size, window_size) ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length) ;

    printf("CIRCULAR ARRAY SIZE IS %lu \n", circular_array_size) ;

    U64 *hashvector = (U64 *)malloc(sizeof(U64) * circular_array_size) ;

    // if (ca == NULL){
    //     printf("ca is null.\n") ;
    //     exit (-1) ;
    // }
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

    // 1 - fill the 
    size_t precompute = 0 ;
    size_t precompute_elements = num_s_mers < circular_array_size ? num_s_mers : circular_array_size ;

    while(precompute < precompute_elements){
        seqhashNext(si, &current_smer);
        // printf("%llu ", current_smer) ;
        hashvector[precompute] = current_smer;
        precompute++;
        computed_smers++;
    }
    // printf("\n") ;

    // 2 - find syncmers
    U64 minimum = U64MAX;
    size_t minimum_position;
    size_t absolute_kmer_position = 0;

    //precompute first window
    for(size_t i =0; i < window_size - 1; i++){
        if (hashvector[i] < minimum){
            minimum = hashvector[i] ;
            minimum_position = i;
        }
    }
    // printf("AT PRECOMPUTE MIN IS %llu, POS IS %lu\n", minimum, minimum_position) ;
    // scan the rest of the array
    for(size_t i = window_size - 1; i < precompute_elements; i++){
         // rescan if minimum out of context
        // printf("min_pos is %lu, comparison is %lu\n", minimum_position, i - window_size + 1) ;
        if(minimum_position < i - window_size + 1){
            // printf("recompute\n");
            size_t curr_pos;
            minimum = U64MAX;
            for(int j = 1; j < window_size; j++){
                curr_pos = i - window_size + j;
                if ( hashvector[curr_pos] < minimum ){
                    // printf("FOUND %llu AT %lu\n",ca->hashVector[scan_position], scan_position ) ;
                    minimum = hashvector[curr_pos] ;
                    minimum_position = curr_pos;
                  }
            }
            // printf("RESCAN. MIN: %llu; POS: %lu\n", minimum, minimum_position) ;
        }

        //verify minimum
        if (hashvector[i] < minimum){
            minimum = hashvector[i] ;
            minimum_position = i;
        }
        // printf("MIN IS %llu, pos is %lu, i is %lu, start is %lu\n", minimum, minimum_position, i, i - window_size + 1) ;
        if(i == minimum_position || minimum_position == i - window_size + 1){
            // printf("%llu\t%lu\n", minimum, minimum_position);
            computed_syncmers++;
        }

        absolute_kmer_position++;
    }

    // free(ca) ;
    free(si) ;
    // free(sync) ;
    printf("COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers) ; 
    printf("HASHED %lu S-MERS\n", computed_smers) ;
    *num_syncmer = computed_syncmers;
}
