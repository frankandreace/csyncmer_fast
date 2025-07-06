#include "hashing.h"
#include "utils.h"
#include <string.h>

/*
 * L2 Cache Budget (512KB total, even if my pc has 10MB):
 * - Hash buffer: 20,480 × 128 bit (16 Bytes) = 320KB (5,120 cache lines)
 * - Syncmer128 cache: 5,120 × 256 bit(32 Bytes) = 160KB (2,560 cache lines)
 * - Total used by the 2 arrays: 480KB (7,680 cache lines)
 * - Free space for the rest (stack/other data) : 32KB (512 cache lines)
 */

// 20,480 U128 elements = 320KB
// 4 elements per 64B cache line = 5,120 cache lines
#define ARRAY_UINT128_SIZE (20480) //20480 40960

// 5,120 Syncmer128 structs = 160KB  
// 2 syncmers per 64B cache line = 2,560 cache lines
#define CACHE_SYNCMER_SIZE (5120)

// Standard cache parameters
#define L2_CACHE_SIZE (512 * 1024)  // 512KB in bytes
#define CACHE_LINE_SIZE 64           // 64 bytes per line
#define ELEMENTS_PER_CACHE_LINE 4    // 64B / 16B = 4 U128s per line
#define NTHASH_FLAG_U128 2

#define SYNG_SEED 7

// to cache syncmers for next()
typedef struct {
    U64 hash_value;
    size_t kmer_position;
    size_t smer_position;
} Syncmer64;


typedef struct SyncmerIteratorS
{
    SeqhashIterator *rolling_hash;

    size_t window_size; // k-mer and s-mer length, window size
    size_t smers_remaining; // number of s-mers to process
    size_t kmer_position;

    U64 *hash_buffer; // s-mer hash buffer
    size_t buffer_size; // size of hash buffer
    size_t start_fill;

    Syncmer64 *cached_syncmer;
    size_t cache_capacity; // size of the cache
    size_t cache_count; // number of syncmers inserted in cache
    size_t cache_index; // pointer to the current cached syncmer
    
    // size_t rescan_count;
    // size_t tot_kmers;
    size_t cache_full;
    size_t num_rolls;
    size_t seq_length;

    bool is_done; // finish generation if no more s-mers

} SyncmerIteratorS;

SyncmerIteratorS* syncmer_generator_createS(char *sequence_input, size_t sequence_length_provided, const size_t K, const size_t S){

    const size_t sequence_length = sequence_length_provided;
    // printf("SEQUENCE LENGTH: %lu\n",sequence_length);

    if (K <= S){
        fprintf(stderr, "K-MER SIZE < S-MER SIZE. EXIT.\n");
        return NULL;
    }
    // ADRESSING CASE SEQUENCE < K
    if (sequence_length < K){
        fprintf(stderr, "SEQUENCE < K-MER SIZE. EXIT.\n");
        return NULL;
    }

    if (S <=2 ){
        fprintf(stderr, "S-MER SIZE HAS TO BE > 2.\n");
        return NULL;
    }

    SyncmerIteratorS* si = (SyncmerIteratorS *)calloc(1,(sizeof(SyncmerIteratorS)));
    si->window_size = K - S + 1;

    if (!si) {
        fprintf(stderr, "Failed to allocate SyncmerIterator structure.\n");
        return NULL;
    }

    // BUILDING HASHING OBJECT
    Seqhash *sh = seqhashCreate(S, si->window_size, SYNG_SEED) ;
    si->rolling_hash = seqhashIterator(sh, sequence_input, sequence_length) ;

    // AND VERIFYING IT IS CONSTRUCTED CORRECTLY
    if (si == NULL){
        fprintf(stderr, "[NT_HASH_BENCH_ERROR]:: Failed to initialize NtHash object. Check sequence/S-mer size. Aborting.\n");
        free(si);
        return NULL;
    }

    si->buffer_size = ARRAY_UINT128_SIZE;
    si->cache_capacity = CACHE_SYNCMER_SIZE;
    si->smers_remaining = sequence_length - S + 1;
    si->is_done = false;
    si->start_fill = si->window_size - 1;
    si->kmer_position = 0;
    si->num_rolls = 0;
    si->seq_length = sequence_length;


    // Allocating memory for hash vector
    si->hash_buffer = (U64 *)aligned_alloc(64, si->buffer_size * sizeof(U64));

    if (!si->hash_buffer) {
        fprintf(stderr, "Could not allocate the hash buffer.\n");
        free(si);
        return NULL;
    }

    // Allocating memory for syncmer cache
    si->cached_syncmer = (Syncmer64 *)aligned_alloc(64, si->cache_capacity * sizeof(Syncmer64));
    if (!si->cached_syncmer) {
        fprintf(stderr, "Could not allocate the cache of syncmer\n");
        SeqhashIteratorDestroy(si->rolling_hash);
        free(si->hash_buffer);
        free(si);
        return NULL;
    }

    // Fill 1st window_size -1 elements in the vector
    for(size_t buffer_position = 0; buffer_position < si->window_size - 1; buffer_position++){
        si->num_rolls++;
        seqhashNext(si->rolling_hash, &si->hash_buffer[buffer_position]);
        // printf("RETURN SEQHASHNEXT: %llu\n", si->hash_buffer[buffer_position]);
    }
    si->smers_remaining -= si->window_size - 1;
    // printf("RETURNING ALL GOOD. (IN THEORY)\n");
    return si;

}

void syncmer_generator_destroyS(SyncmerIteratorS* si){
    if(!si) return;
    SeqhashIteratorDestroy(si->rolling_hash);
    free(si->hash_buffer);
    free(si->cached_syncmer);
    free(si);
}

static void process_chunk_and_cacheS(SyncmerIteratorS *si){
    // If in this function, the cache is empty
    si->cache_count = 0;
    si->cache_index = 0;
    if(si->smers_remaining == 0) return;
    // printf("STARTING. cache count: %lu ; cache_index: %lu.\n",si->cache_count,si->cache_index);
    // 1 - FILL THE ARRAY OF HASHESH WITH HASHES
    // STOP CONDITION: NO MORE HASHES OR ARRAY FULL
    // NORMAL CASE: K-S hashes are already inserted
    // EDGE CASE: IN THE PREVIOUS CALL, THE SYNCMER BUFFER WAS FILLED AND THERE ARE ALREADY PRE-FILLED MORE THAN K-S HASHES FRO THE PREVIOUS CALL
    // THE FIRST HASH BUFFER POSITION TO BE INSERTED IS STORED IN START_FILL
    // U64 synmcer = U64MAX;
    size_t buffer_position = si->start_fill;
    size_t end_positon = (si->buffer_size <= si->smers_remaining + si->start_fill) ? si->buffer_size: si->smers_remaining + si->start_fill;
    // printf("Remaning s-mers to compute: %lu. Num tot s-mers - num_rolls: %lu\n", si->smers_remaining, (si->seq_length - si->num_rolls));
    si->smers_remaining -= end_positon - si->start_fill;
    // printf("NUM ELEMENTS TO FILL: %lu\n",end_positon - buffer_position);
    // printf("Remaning s-mers to compute: %lu. Num tot s-mers - num_rolls: %lu\n",si->smers_remaining, (si->seq_length - si->num_rolls));


    // Hashing new s-mers till either the end of the vector or the last s-mer in the sequence.
    while(buffer_position < end_positon){
        // bool return_seqhashnext;
        seqhashNext(si->rolling_hash, &si->hash_buffer[buffer_position]);
        // printf("RETURN SEQHASHNEXT: %llu\n", si->hash_buffer[buffer_position]);
        // si->hash_buffer[buffer_position] = synmcer;
        buffer_position++;
        si->num_rolls++;
    }
    // printf("FINISHED HASHING.\n");
    // printf("HASH AT POS 0: %lu\n", si->hash_buffer[0]);
    // 2 - COMPUTE CLOSED SYNCMERS AND ADD THEM IN THE SYNCMER CACHE
    // starting at the k-s+1 position to verify syncmer
    size_t min_hash_position = 0;

    // scan first w-1 to get minimum
    for(size_t i = 0; i < si->window_size - 1; i++){
        // printf("HASH_BUFFER[%lu (i)]: %llu <? HASH_BUFFER[%lu (min_hash_pos)]: %llu\n",i,si->hash_buffer[i],min_hash_position, si->hash_buffer[min_hash_position]);
        if (si->hash_buffer[i] < si->hash_buffer[min_hash_position]) {
            min_hash_position = i;
        }
    }


    // scan the rest of the array for syncmers
    for (size_t i = si->window_size - 1; i < end_positon; i++){

        // rescan last w-1 elements if minimum is out of context
        if(min_hash_position < i - si->window_size + 1){
            // si->rescan_count++;
            // the current position is i - window_size + j
            min_hash_position = i - si->window_size + 1;
            for (size_t j = 2; j <= si->window_size; j++){
                if(si->hash_buffer[i - si->window_size + j] < si->hash_buffer[min_hash_position]){
                    min_hash_position = i - si->window_size + j;
                }
            }
        }
        // update minimum
        else if (si->hash_buffer[i] < si->hash_buffer[min_hash_position]){
            min_hash_position = i;
        }

        //check syncmer condition
        if (i == min_hash_position || min_hash_position ==  i - si->window_size + 1){
            // add syncmer to the cache of syncmer
            // printf("MIN HASH POS: %lu; hash: %llu, pos: %lu\n", min_hash_position, si->hash_buffer[min_hash_position]);
            si->cached_syncmer[si->cache_count].hash_value = si->hash_buffer[min_hash_position];
            si->cached_syncmer[si->cache_count].kmer_position = si->kmer_position;
            si->cached_syncmer[si->cache_count].smer_position = si->kmer_position + si->window_size - (i - min_hash_position + 1);

            si->cache_count++;

             // IF CACHE IS FULL, MOVE THE ELEMENTS TO THE BEGINNING
            if (si->cache_count == si->cache_capacity){
                // set new start to fill
                si->start_fill = si->buffer_size - i + si->window_size - 1;
                // do memmove 
                memmove(si->hash_buffer, si->hash_buffer + i - si->window_size + 1, si->start_fill * sizeof(U64));
                si->cache_full++;
                si->kmer_position++;
                // return
                printf("RETURNING WITH CACHE FULL. cache count: %lu ; cache_index: %lu.\n",si->cache_count,si->cache_index);
                return;
            }
        }
        //advance kmer condition
        si->kmer_position++;
    }


    // IF VECTOR IS PROCESSED, TRANSFER LAST W-1 ELEMENTS TO THE BEGINNING
    if (si->smers_remaining > 0 && end_positon >= si->window_size){
        si->start_fill = si->window_size - 1;
        memmove(si->hash_buffer, si->hash_buffer + end_positon - si->window_size + 1, si->start_fill * sizeof(U64));
    }

    // RETURN
    // printf("RETURNING. cache count: %lu ; cache_index: %lu.\n",si->cache_count,si->cache_index);
    return;

}

bool syncmer_iterator_nextS(SyncmerIteratorS * si, Syncmer64 * syncmer){
    if(!si || !syncmer || si->is_done) {
        // if (!si) {printf("SI IS NONE\n");}
        // else if (!syncmer) {printf("SYNCMER IS NONE");}
        // else {printf("SI->IS_DONE IS TRUE");}
        // printf("Either SI, SYNCMER or IS_DONE ARE DONE.\n");
        return false;
    }
    // printf("NEXT, # SMERS REMAINING: %lu\n", si->smers_remaining);
    // If no syncmer in the buffer, compute again
    if(si->cache_index >= si->cache_count){
        // printf("Processing chunk. si->cache_index (%lu) >= si->cache_count (%lu).\n", si->cache_index, si->cache_count);
        process_chunk_and_cacheS(si);
    }
    // printf("NEXT, # SMERS REMAINING AFTER PROCESS_: %lu ; cache_count: %lu\n", si->smers_remaining, si->cache_count);
    // If no syncmer, process until either it finds more or it consumes the sequence
    while(si->cache_count == 0){
        // If no more s-mers, signal computation ended
        if( si->smers_remaining == 0){
            // printf("Had to rescan %lu times out of %lu k-mers.\n", si->rescan_count, si->tot_kmers);
            // printf("Cache was full %lu times.\n", si->cache_full);
            // printf("Num computed s-mers: %lu\n", si->num_rolls);
            si->is_done = true;
            return false;
        }
        // printf("Processing chunk. si->cache_count == 0.\n");
        process_chunk_and_cacheS(si);
    }

    // return one of the syncmer in the buffer
    *syncmer = si->cached_syncmer[si->cache_index++];

    return true;
}

