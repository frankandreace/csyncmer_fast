#include <stddef.h>

#include "csyncmer_fast/hashing.h"
#include "circular_array.h"

// wrap functions to compute closed syncmers
typedef struct
{
    SeqhashIterator *si ;
    CircularArray *ca;

    U64 current_smer ;
    U64 current_min_smer ;

    size_t current_kmer_position ;
    size_t smer_position_in_kmer ;
    size_t current_min_smer_position ;

    bool isDone;

} syncmerIterator;

void syncmerIteratorDestroy(syncmerIterator *syncit){
    SeqhashIteratorDestroy(syncit->si);
    circularArrayDestroy(syncit->ca) ;
    free(syncit);
}

syncmerIterator *syncmerIteratorInitialization(int K, int S, char *sequence, size_t sequence_length){
    
    // memory allocation
    syncmerIterator *syncit = (syncmerIterator *)malloc( sizeof(syncmerIterator) ) ;

    // initialization of structures and variables
    U64 seed  = 7 ;
    size_t window_size = K - S + 1 ;

    Seqhash *sh = seqhashCreate(S, window_size, seed) ;
    if (sh == NULL){ printf("Error initializing the seqhash structure.\n") ; exit (-1) ; }

    syncit->si = seqhashIterator(sh, sequence, sequence_length) ;
    if (syncit->si == NULL){ printf("Error initializing the seqhashiterator structure.\n") ; exit (-1) ; }

    syncit->ca = circularArrayCreate( window_size) ;
    if (syncit->ca == NULL){ printf("Error initializing the circularArray structure.\n") ; exit (-1) ; }

    syncit->current_smer = 0 ;
    syncit->current_min_smer = 0 ;

    syncit->smer_position_in_kmer = 0 ;
    syncit->current_min_smer_position = 0 ;
    syncit->current_kmer_position = 0 ;

    syncit->isDone = false;

    // 1 - precompute 1st window
    size_t precompute = 0 ;
    while(precompute < window_size - 1){
        bool x = seqhashNext(syncit->si, &syncit->current_smer);

        if(!x){
            exit(1) ;
        }

        circularInsert(syncit->ca,syncit->current_smer);
        precompute++;
    }

    // 2 - run until end of given sequence
    while(seqhashNext(syncit->si, &syncit->current_smer)){

        circularInsert(syncit->ca, syncit->current_smer);

        if (is_syncmer(syncit->ca, &syncit->smer_position_in_kmer)){
            
            syncit->current_min_smer_position = syncit->current_kmer_position + syncit->smer_position_in_kmer ;
            syncit->current_min_smer = get_smer(syncit->ca);
            return syncit; 

        }

        syncit->current_kmer_position++ ;
    }

    // if no syncmer is found in the sequence (too short to find one), setting IsDone to true.
    syncit->isDone = true;
    return syncit;
}

bool get_next_syncmer(syncmerIterator *syncit, size_t *kmer_pos, size_t *smer_pos, U64 *smer){

    // add elements for return
    if(kmer_pos) *kmer_pos = syncit->current_kmer_position ;
    if(smer_pos) *smer_pos = syncit->current_min_smer_position ;
    if(smer) *smer = syncit->current_min_smer ;

    if(syncit->isDone){return false;}

    // compute next syncmer
    while(seqhashNext(syncit->si, &syncit->current_smer)){

        circularInsert(syncit->ca, syncit->current_smer);

        if (is_syncmer(syncit->ca, &syncit->smer_position_in_kmer)){
            
            syncit->current_min_smer_position = syncit->current_kmer_position + syncit->smer_position_in_kmer ;
            syncit->current_min_smer = get_smer(syncit->ca);
            return true; 

        }

        syncit->current_kmer_position++ ;
    }

    // finished reading sequence
    syncit->isDone = true;
    return true;

}