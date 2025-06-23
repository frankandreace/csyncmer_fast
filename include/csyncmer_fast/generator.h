#include "nthash_wrapper.h"
#include "utils.h"

// 512KB: Assumption of 512Mb M2
// TOT LINE OF CACHE = 1024 * 8 -> 8192
// 8192 - 5125 - 2560
// 507 lines of cache for the rest

#define ARRAY_UINT128_SIZE (1024 * 4 * 5) // 20480 elements ( 4 per lin of cache -> 1024 * 5 = 5125)
#define CACHE_SYNCMER_SIZE (1024 * 5) // 5120 elements (2 per line of cache -> 512*5 = 2560 lines of cache)
#define L2_CACHE_SIZE (512 * 1024) // in bytes
#define CACHE_LINE_SIZE 64 // size of a cache line in bytes (64*8 = 512 bit)
#define ELEMENTS_PER_CACHE_LINE (CACHE_LINE_SIZE / sizeof(U128)) // 4 elements per cache line as U128 is 16 bytes
#define NTHASH_FLAG_U128 2

// to cache syncmers for next()
typedef struct {
    U128 hash_value;
    size_t kmer_position;
    size_t smer_position;
} Syncmer;


typedef struct SyncmerIterator
{
    NtHashHandle rolling_hash;

    size_t k, s, window_size; // k-mer and s-mer length, window size
    size_t smers_remaining; // number of s-mers to process
    size_t kmer_positon;

    U128 *hash_buffer; // s-mer hash buffer
    size_t buffer_size; // size of hash buffer
    U128 current_hash; // current hash value
    size_t current_hash_position; // current_positon in hash buffer
    size_t start_fill;

    Syncmer *cached_syncmer;
    size_t cache_capacity; // size of the cache
    size_t cache_count; // number of syncmers inserted in cache
    size_t cache_index; // pointer to the current cached syncmer

    bool is_done; // finish generation if no more s-mers

} SyncmerIterator;

SyncmerIterator* syncmer_generator_create(const char *sequence_input, const size_t sequence_length, const size_t K, const size_t S){
    
    // ADRESSING CASE SEQUENCE < K
    if (sequence_length < K){
        printf("SEQUENCE < KMER SIZE. EXIT.\n");
    }

    SyncmerIterator* si = (SyncmerIterator *)calloc(1,(sizeof(SyncmerIterator)));
    
    // BUILDING HASHING OBJECT
    si->rolling_hash = nthash_create(sequence_input,strlen(sequence_input), S, NTHASH_FLAG_U128);
    
    // AND VERIFYING IT IS CONSTRUCTED CORRECTLY
    if (si->rolling_hash == NULL){
        fprintf(stderr, "[NT_HASH_BENCH_ERROR]:: Failed to initialize NtHash object. Check sequence/S-mer size. Aborting.\n");
        return;
    }

    si->buffer_size = ARRAY_UINT128_SIZE;
    si->cache_capacity = CACHE_SYNCMER_SIZE;
    si->smers_remaining = sequence_length - S + 1;
    si->s = S;
    si->k = K;
    si->window_size = K - S + 1;
    si->is_done = false;
    si->start_fill = si->window_size - 1;
    si->kmer_positon = 0;

    // Allocating memory for hash vector
    si->hash_buffer = (U128 *)aligned_alloc(64, si->buffer_size * sizeof(U128));

    // Allocating memory for syncmer cache
    si->cached_syncmer = (Syncmer *)aligned_alloc(64, si->cache_capacity * sizeof(Syncmer));

    // Fill 1st window_size -1 elements in the vector
    for(size_t buffer_position = 0; buffer_position < si->window_size - 1; buffer_position++){
        nthash_roll(si->rolling_hash);
        si->hash_buffer[buffer_position] = nthash_get_canonical_hash_128(si->rolling_hash);
    }

    return si;

}

void syncmer_generator_destroy(SyncmerIterator* si){
    if(!si) return;
    nthash_destroy(si->rolling_hash);
    free(si->hash_buffer);
    free(si->cached_syncmer);
    free(si);
}

static void process_chunk_and_cache(SyncmerIterator *si){
    // If in this function, the cache is empty
    si->cache_count = 0;
    si->cache_index = 0;

    // 1 - FILL THE ARRAY OF HASHESH WITH HASHES
    // STOP CONDITION: NO MORE HASHES OR ARRAY FULL
    // NORMAL CASE: K-S hashes are already inserted
    // EDGE CASE: IN THE PREVIOUS CALL, THE SYNCMER BUFFER WAS FILLED AND THERE ARE ALREADY PRE-FILLED MORE THAN K-S HASHES FRO THE PREVIOUS CALL
    // THE FIRST HASH BUFFER POSITION TO BE INSERTED IS STORED IN START_FILL

    size_t buffer_position = si->start_fill;
    size_t end_positon = (si->buffer_size <= si->smers_remaining + si->start_fill) ? si->buffer_size: si->smers_remaining + si->start_fill;
    si->smers_remaining -= end_positon - si->start_fill;


    // Hashing new s-mers till either the end of the vector or the last s-mer in the sequence.
    while(buffer_position < end_positon){
        nthash_roll(si->rolling_hash);
        si->hash_buffer[buffer_position] = nthash_get_canonical_hash_128(si->rolling_hash);
        buffer_position++;
    }

    // 2 - COMPUTE CLOSED SYNCMERS AND ADD THEM IN THE SYNCMER CACHE


    // IF CACHE IS FULL, MOVE THE ELEMENTS TO THE BEGINNING
    if (si->cache_count == si->cache_capacity){
        // set new start to fill
        si->start_fill = si->buffer_size - si->current_hash_position + si->window_size - 1;

        // do memmove 
        memmove(si->hash_buffer, si->hash_buffer + si->current_hash_position - si->window_size + 1, si->start_fill * sizeof(U128));
        // return
        return;
    }


    // IF VECTOR IS PROCESSED, TRANSFER LAST W-1 ELEMENTS TO THE BEGINNING
    if (si->smers_remaining > 0 && end_positon >= si->window_size){
        si->start_fill = si->window_size - 1;
        memmove(si->hash_buffer, si->hash_buffer + si->buffer_size - si->window_size + 1, si->start_fill * sizeof(U128));
    }

    // RETURN
    return;

}

bool syncmer_iterator_next(SyncmerIterator * si, Syncmer * syncmer){
    if(!si || !syncmer || si->is_done) return false;

    // If no syncmer in the buffer, compute again
    if(si->cache_index >= si->cache_count){
        process_chunk_and_cache(si);
    }

    // If no syncmer, process until either it finds more or it consumes the sequence
    while(si->cache_count == 0){
        // If no more s-mers, signal computation ended
        if( si->smers_remaining == 0){
            si->is_done = true;
            return false;
        }
        process_chunk_and_cache(si);
    }

    // return one of the syncmer in the buffer
    *syncmer = si->cached_syncmer[si->cache_index++];

    return true;
}

