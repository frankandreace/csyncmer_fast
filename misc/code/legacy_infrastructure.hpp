#ifndef LEGACY_INFRASTRUCTURE_H
#define LEGACY_INFRASTRUCTURE_H

// Legacy infrastructure for SeqHash, CircularArray, NtHash wrappers, etc.
// Used by syncmer_seqhash.h and other legacy implementations

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
#include <nthash/nthash.hpp>
#endif

// ============================================================================
// COMMON TYPES AND UTILITIES
// ============================================================================

#ifndef UNIT_DEFINED
#define UNIT_DEFINED

typedef unsigned long long U64;
typedef __uint128_t U128;
const static U64 U64MAX = 0xffffffffffffffff;
const static __uint128_t U128MAX = (__uint128_t)(-1);

// Helper function to print __uint128_t in hexadecimal (still needed for printf)
void print_uint128_hex_debug(__uint128_t val) {
    uint64_t high = (uint64_t)(val >> 64);
    uint64_t low = (uint64_t)val;
    fprintf(stderr, "%016lX%016lX", high, low); // Print to stderr for debug
}

#endif

// ============================================================================
// SEQHASH IMPLEMENTATION
// ============================================================================

// BASE TO BITS ARRAY
// USING CHAR CONVERTION TO INT
// DEFAULT IS 0
// 65: A ; 97: a
// 67: C ; 99: c
// 71: G ; 103: g
// 84: T ; 116: t
// 85: U ; 117: u
// first 4 elements as it is receiving the sequence already binarized
static const unsigned char base_to_bits_array[256] = {
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0};

static inline char base_to_bits(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0; // Treat Ns and  unknown as 'A'
    }
}

typedef struct
{
  int seed; /* seed */
  int k;    /* kmer */
  int w;    /* window */
  U64 mask; /* 2*k bits */
  int shift1, shift2;
  U64 factor1, factor2;
  U64 patternRC[4]; /* one per base */
} Seqhash;

typedef struct
{
  Seqhash *sh;
  char *s, *sEnd;   /* sequence currently being hashed, end marker */
  U64 h, hRC, hash; /* current k-mer values */
  bool isDone;
} SeqhashIterator;

/*---- HASHING PART ----*/

static inline U64 kHash(Seqhash *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1); }

static inline U64 hashRC(SeqhashIterator *si, bool *isForward)
{
  U64 hashF = kHash(si->sh, si->h);
  U64 hashR = kHash(si->sh, si->hRC);
  if (hashF < hashR)
  {
    *isForward = true;
    return hashF;
  }
  else
  {
    *isForward = false;
    return hashR;
  }
}

static inline U64 advanceHashRC(SeqhashIterator *si, bool *isForward)
{
  Seqhash *sh = si->sh;
  if (si->s < si->sEnd)
  {
    si->h = ((si->h << 2) & sh->mask) | base_to_bits_array[*si->s];
    si->hRC = (si->hRC >> 2) | sh->patternRC[base_to_bits_array[*si->s]];
    ++si->s;
    return hashRC(si, isForward);
  }
  else
    return U64MAX;
}

Seqhash *seqhashCreate(int k, int w, int seed)
{
  assert(sizeof(U64) == 8);
  Seqhash *sh = (Seqhash *)malloc(sizeof(Seqhash));
  sh->k = k;
  if (k < 1 || k >= 32)
  {
    fprintf(stderr, "seqhash k %d must be between 1 and 32\n", k);
    exit(-1);
  }
  sh->w = w;
  if (w < 1)
  {
    fprintf(stderr, "seqhash w %d must be positive\n", w);
    exit(-1);
  }
  sh->seed = seed;
  sh->mask = ((U64)1 << (2 * k)) - 1;
  int i;

  srandom(seed);
  sh->factor1 = (random() << 32) | random() | 0x01;
  sh->shift1 = 64 - 2 * k;
  sh->factor2 = (random() << 32) | random() | 0x01;
  sh->shift2 = 2 * k;
  for (i = 0; i < 4; ++i)
  {
    sh->patternRC[i] = (3 - i);
    sh->patternRC[i] <<= 2 * (k - 1);
  }
  return sh;
}

SeqhashIterator *seqhashIterator(Seqhash *sh, char *s, int len)
{
  assert(s && len >= 0);
  SeqhashIterator *si = (SeqhashIterator *)malloc(sizeof(SeqhashIterator));
  si->sh = sh;
  si->s = s;
  si->sEnd = s + len;
  si->hash = 0;
  si->h = 0;
  si->hRC = 0;
  si->isDone = false;
  bool isForward;
  if (len < sh->k)
  {
    si->isDone = true;
  } // edge case
  else
  {
    int i; /* preinitialise the hashes for the first kmer */
    for (i = 0; i < sh->k; ++i, ++si->s)
    {
      base_to_bits_array[*si->s];
      si->h = (si->h << 2) | base_to_bits_array[*si->s];
      si->hRC = (si->hRC >> 2) | sh->patternRC[base_to_bits_array[*si->s]];
    }
    si->hash = hashRC(si, &isForward);
  }
  return si;
}

static void SeqhashDestroy(Seqhash *sh) { free(sh); }

static void SeqhashIteratorDestroy(SeqhashIterator *si)
{
  // prevent double free
  if (si == NULL)
    return;
  if (si->sh != NULL)
  {
    // Free the seqhash
    SeqhashDestroy(si->sh);
    si->sh = NULL; // prevent double free
  }
  // free the struct
  free(si);
}

bool seqhashNext(SeqhashIterator *si, U64 *kmer, bool *isForward)
{
  if (si->isDone)
  {
    return false; /* we are done */
  }
  if (kmer)
    *kmer = si->hash;

  if (si->s >= si->sEnd)
  {
    si->isDone = true;
  }
  else
    si->hash = advanceHashRC(si, isForward);
  return true;
}

// ============================================================================
// SEQHASH-BASED SYNCMER ITERATOR
// ============================================================================

/*
 * L2 Cache Budget (512KB total, even if my pc has 10MB):
 * - Hash buffer: 20,480 × 128 bit (16 Bytes) = 320KB (5,120 cache lines)
 * - Syncmer128 cache: 5,120 × 256 bit(32 Bytes) = 160KB (2,560 cache lines)
 * - Total used by the 2 arrays: 480KB (7,680 cache lines)
 * - Free space for the rest (stack/other data) : 32KB (512 cache lines)
 */

// 20,480 U64 elements = 160KB
// 8 elements per 64B cache line = 2560 cache lines
#define ARRAY_UINT64_SIZE (40960) //20480 40960

// 5,120 Syncmer64 structs = 160KB
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
    bool smer_at_start;
    bool is_forward;
} Syncmer64;

typedef struct SyncmerIteratorS
{
    SeqhashIterator *rolling_hash;

    size_t window_size; // k-mer and s-mer length, window size
    size_t smers_remaining; // number of s-mers to process
    size_t kmer_position;

    U64 *hash_buffer; // s-mer hash buffer
    bool *orientation_buffer; // s-mer forward or reversed orientation buffer
    size_t buffer_size; // size of hash buffer
    size_t start_fill;

    Syncmer64 *cached_syncmer;
    size_t cache_capacity; // size of the cache
    size_t cache_count; // number of syncmers inserted in cache
    size_t cache_index; // pointer to the current cached syncmer
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

    si->buffer_size = ARRAY_UINT64_SIZE;
    si->cache_capacity = CACHE_SYNCMER_SIZE;
    si->smers_remaining = sequence_length - S + 1;
    si->is_done = false;
    si->start_fill = si->window_size - 1;
    si->kmer_position = 0;
    si->num_rolls = 0;
    si->seq_length = sequence_length;
    bool test;


    // Allocating memory for hash vector
    si->hash_buffer = (U64 *)aligned_alloc(64, si->buffer_size * sizeof(U64));
    si->orientation_buffer = (bool *)aligned_alloc(8, si->buffer_size * sizeof(bool));
    // || !si->orientation_buffer
    if (!si->hash_buffer) {
        fprintf(stderr, "Could not allocate the hash or orientation buffer.\n");
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
        seqhashNext(si->rolling_hash, &si->hash_buffer[buffer_position], &test);
    }
    si->smers_remaining -= si->window_size - 1;
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

    size_t buffer_position = si->start_fill;
    size_t end_positon = (si->buffer_size <= si->smers_remaining + si->start_fill) ? si->buffer_size: si->smers_remaining + si->start_fill;
    // printf("Remaning s-mers to compute: %lu. Num tot s-mers - num_rolls: %lu\n", si->smers_remaining, (si->seq_length - si->num_rolls));
    si->smers_remaining -= end_positon - si->start_fill;
    // printf("NUM ELEMENTS TO FILL: %lu\n",end_positon - buffer_position);
    // printf("Remaning s-mers to compute: %lu. Num tot s-mers - num_rolls: %lu\n",si->smers_remaining, (si->seq_length - si->num_rolls));
    // Hashing new s-mers till either the end of the vector or the last s-mer in the sequence.
    while(buffer_position < end_positon){
        // bool return_seqhashnext;
        seqhashNext(si->rolling_hash, &si->hash_buffer[buffer_position], &si->orientation_buffer[buffer_position]);
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
            // si->cached_syncmer[si->cache_count].smer_at_start = (i == min_hash_position) ? false : true;
            // si->cached_syncmer[si->cache_count].is_forward = si->orientation_buffer[min_hash_position];
            // si->cached_syncmer[si->cache_count].smer_position = si->kmer_position + si->window_size - (i - min_hash_position + 1);

            si->cache_count++;

             // IF CACHE IS FULL, MOVE THE ELEMENTS TO THE BEGINNING
            if (si->cache_count == si->cache_capacity){
                // set new start to fill
                si->start_fill = si->buffer_size - i + si->window_size - 1;
                // do memmove
                memmove(si->hash_buffer, si->hash_buffer + i - si->window_size + 1, si->start_fill * sizeof(U64));
                memmove(si->orientation_buffer, si->orientation_buffer + i - si->window_size + 1, si->start_fill * sizeof(bool));
                si->cache_full++;
                si->kmer_position++;
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
        memmove(si->orientation_buffer, si->orientation_buffer + end_positon - si->window_size + 1, si->start_fill * sizeof(bool));

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

// ============================================================================
// NTHASH WRAPPER
// ============================================================================

// Opaque pointer to hide C++ class details from C
typedef void* NtHashHandle;

#ifdef __cplusplus
extern "C" { // Ensure C linkage for these functions
#endif

/**
 * @brief Creates a new NtHash object for a given sequence and k-mer length.
 * @param sequence The DNA/RNA sequence string.
 * @param k The k-mer length.
 * @param hash_length 1 for 64-bit canonical, 2 for 128-bit canonical.
 * @return An opaque handle to the NtHash object, or NULL on failure.
 */
NtHashHandle nthash_create(const char* sequence, size_t sequence_length, unsigned int k, unsigned int hash_length_flag);

/**
 * @brief Rolls the hash window to the next k-mer.
 * @param handle Handle to the NtHash object.
 * @return 1 if successful (more k-mers), 0 if no more k-mers.
 */
int nthash_roll(NtHashHandle handle);

/**
 * @brief Gets the current forward hash (always 64-bit).
 * @param handle Handle to the NtHash object.
 * @return The 64-bit forward hash.
 */
uint64_t nthash_get_forward_hash(NtHashHandle handle);

/**
 * @brief Gets the current canonical hash.
 * Returns a 128-bit hash. If 64-bit was requested during creation,
 * the higher 64 bits of the __uint128_t will be zero.
 * @param handle Handle to the NtHash object.
 * @return The 128-bit canonical hash.
 */
U128 nthash_get_canonical_hash_128(NtHashHandle handle);

/**
 * @brief Destroys the NtHash object and frees allocated memory.
 * @param handle Opaque handle to the NtHash object.
 */
void nthash_destroy(NtHashHandle handle);

#ifdef __cplusplus
} // extern "C"
#endif

// --- C++ Implementation (Only visible to C++ compilers) ---
// This part contains the actual C++ code that implements the C interface.
// It's placed *after* the extern "C" block and is only compiled if __cplusplus is defined.

#ifdef __cplusplus

// Include necessary C++ headers here.
// These will *not* be included when a C compiler processes this header.
#include <nthash/nthash.hpp> // The actual ntHash C++ library header
#include <vector>             // For std::vector
#include <stdexcept>          // For std::runtime_error (or other exceptions)
#include <iostream>

// Anonymous namespace or internal linkage to prevent symbol clashes
namespace{

// Helper to convert void* handle to C++ NtHash*
inline nthash::NtHash* get_nthash_obj(NtHashHandle handle){
    if (!handle) {
        throw std::runtime_error("Invalid NtHash handle.");
    }
    return static_cast<nthash::NtHash*>(handle);
}
}// anonymous namespace

// Implementations of the C functions, using the extern "C" linkage.
// Marked as 'inline' to suggest the compiler to inline them,
// which is standard for header-only libraries.
extern "C" inline NtHashHandle nthash_create(const char* sequence, size_t sequence_length, unsigned int k, unsigned int hash_length_flag){
    try{
        return new nthash::NtHash(sequence, sequence_length, hash_length_flag, k);
    } catch (const std::exception& e){
        std::cerr << "CONSTRUCTION OF NTHASH FAILED MISERABLY.\n" << std::endl;
        return nullptr;
    }
}

extern "C" inline int nthash_roll(NtHashHandle handle){
    try {
        return get_nthash_obj(handle)->roll() ? 1: 0;
    } catch (...) {return 0; }
}

extern "C" inline uint64_t nthash_get_forward_hash(NtHashHandle handle){
    try {
        return get_nthash_obj(handle)->get_forward_hash();
    } catch (...) {return 0; }
}

extern "C" inline U128 nthash_get_canonical_hash_128(NtHashHandle handle){
    U128 result = 0;
    try {
        const uint64_t* hashes_ptr = get_nthash_obj(handle)->hashes();
        result = ((U128)hashes_ptr[1] << 64) | hashes_ptr[0];
    } catch (...) {}
    return result;
}

extern "C" inline void nthash_destroy(NtHashHandle handle) {
    try {
        // before I deleted get_nthash_obj(handle);
        delete static_cast<nthash::NtHash*>(handle); // Use get_nthash_obj to check handle validity
    } catch (...) { /* ignore errors during destruction */ }
}

#endif // __cplusplus

// ============================================================================
// NTHASH-BASED SYNCMER ITERATOR
// ============================================================================

/*
 * L2 Cache Budget (512KB total, even if my pc has 10MB):
 * - Hash buffer: 20,480 × 128 bit (16 Bytes) = 320KB (5,120 cache lines)
 * - Syncmer128 cache: 5,120 × 256 bit(32 Bytes) = 160KB (2,560 cache lines)
 * - Total used by the 2 arrays: 480KB (7,680 cache lines)
 * - Free space for the rest (stack/other data) : 32KB (512 cache lines)
 */

// 20,480 U128 elements = 320KB
// 4 elements per 64B cache line = 5,120 cache lines
#define ARRAY_UINT128_SIZE (20480) //20480

// to cache syncmers for next()
typedef struct {
    U128 hash_value;
    size_t kmer_position;
    size_t smer_position;
} Syncmer128;


typedef struct SyncmerIterator
{
    NtHashHandle rolling_hash;

    size_t window_size; // k-mer and s-mer length, window size
    size_t smers_remaining; // number of s-mers to process
    size_t kmer_position;

    U128 *hash_buffer; // s-mer hash buffer
    size_t buffer_size; // size of hash buffer
    size_t start_fill;

    Syncmer128 *cached_syncmer;
    size_t cache_capacity; // size of the cache
    size_t cache_count; // number of syncmers inserted in cache
    size_t cache_index; // pointer to the current cached syncmer
    size_t cache_full;
    // size_t num_rolls;
    // size_t seq_length;

    bool is_done; // finish generation if no more s-mers

} SyncmerIterator;

SyncmerIterator* syncmer_generator_create(const char *sequence_input, const size_t K, const size_t S){

    const size_t sequence_length = strlen(sequence_input);
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

    SyncmerIterator* si = (SyncmerIterator *)calloc(1,(sizeof(SyncmerIterator)));
    if (!si) {
        fprintf(stderr, "Failed to allocate SyncmerIterator structure.\n");
        return NULL;
    }

    // BUILDING HASHING OBJECT
    si->rolling_hash = nthash_create(sequence_input, sequence_length, S, NTHASH_FLAG_U128);

    // AND VERIFYING IT IS CONSTRUCTED CORRECTLY
    if (si->rolling_hash == NULL){
        fprintf(stderr, "[NT_HASH_BENCH_ERROR]:: Failed to initialize NtHash object. Check sequence/S-mer size. Aborting.\n");
        free(si);
        return NULL;
    }

    si->buffer_size = ARRAY_UINT128_SIZE;
    si->cache_capacity = CACHE_SYNCMER_SIZE;
    si->smers_remaining = sequence_length - S + 1;
    si->window_size = K - S + 1;
    si->is_done = false;
    si->start_fill = si->window_size - 1;
    si->kmer_position = 0;
    // si->cache_full = 0;
    // si->num_rolls = 0;
    // si->seq_length = sequence_length - S + 1;

    // Allocating memory for hash vector
    si->hash_buffer = (U128 *)aligned_alloc(64, si->buffer_size * sizeof(U128));
    if (!si->hash_buffer) {
        fprintf(stderr, "Could not allocate the hash buffer.\n");
        nthash_destroy(si->rolling_hash);
        free(si);
        return NULL;
    }

    // Allocating memory for syncmer cache
    si->cached_syncmer = (Syncmer128 *)aligned_alloc(64, si->cache_capacity * sizeof(Syncmer128));
    if (!si->cached_syncmer) {
        fprintf(stderr, "Could not allocate the cache of syncmer\n");
        nthash_destroy(si->rolling_hash);
        free(si->hash_buffer);
        free(si);
        return NULL;
    }

    // Fill 1st window_size -1 elements in the vector
    for(size_t buffer_position = 0; buffer_position < si->window_size - 1; buffer_position++){
        nthash_roll(si->rolling_hash);
        // si->num_rolls++;
        si->hash_buffer[buffer_position] = nthash_get_canonical_hash_128(si->rolling_hash);
    }
    si->smers_remaining -= si->window_size - 1;

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

    // If no new s-mers to process, we're done
    if (end_positon <= si->start_fill) {
        si->is_done = true;
        return;
    }

    si->smers_remaining -= end_positon - si->start_fill;

    // Hashing new s-mers till either the end of the vector or the last s-mer in the sequence.
    while(buffer_position < end_positon){
        nthash_roll(si->rolling_hash);
        si->hash_buffer[buffer_position] = nthash_get_canonical_hash_128(si->rolling_hash);
        buffer_position++;
    }

    // 2 - COMPUTE CLOSED SYNCMERS AND ADD THEM IN THE SYNCMER CACHE
    // starting at the k-s+1 position to verify syncmer
    size_t min_hash_position = 0;

    // scan first w-1 to get minimum
    for(size_t i = 0; i < si->window_size - 1; i++){
        if (si->hash_buffer[i] < si->hash_buffer[min_hash_position]) {
            min_hash_position = i;
        }
    }


    // scan the rest of the array for syncmers
    for (size_t i = si->window_size - 1; i < end_positon; i++){

        // rescan last w-1 elements if minimum is out of context
        if(min_hash_position < i - si->window_size + 1){
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

            si->cached_syncmer[si->cache_count].hash_value = si->hash_buffer[min_hash_position];
            si->cached_syncmer[si->cache_count].smer_position = si->kmer_position - (i - min_hash_position);
            si->cached_syncmer[si->cache_count].kmer_position = si->kmer_position;
            si->cache_count++;

             // IF CACHE IS FULL, MOVE THE ELEMENTS TO THE BEGINNING
            if (si->cache_count == si->cache_capacity){
                // set new start to fill
                si->start_fill = si->buffer_size - i + si->window_size - 1;
                // do memmove
                memmove(si->hash_buffer, si->hash_buffer + i - si->window_size + 1, si->start_fill * sizeof(U128));
                si->cache_full++;
                si->kmer_position++;
                // return
                return;
            }
        }
        //advance kmer condition
        si->kmer_position++;
    }


    // IF VECTOR IS PROCESSED, TRANSFER LAST W-1 ELEMENTS TO THE BEGINNING
    if (si->smers_remaining > 0 && end_positon >= si->window_size){
        si->start_fill = si->window_size - 1;
        memmove(si->hash_buffer, si->hash_buffer + end_positon - si->window_size + 1, si->start_fill * sizeof(U128));
    }

    // RETURN
    return;

}

bool syncmer_iterator_next(SyncmerIterator * si, Syncmer128 * syncmer){
    if(!si || !syncmer || si->is_done) {
        // if (!si) {printf("SI IS NONE\n");}
        // else if (!syncmer) {printf("SYNCMER IS NONE");}
        // else {printf("SI->IS_DONE IS TRUE");}
        // printf("Either SI, SYNCMER or IS_DONE ARE DONE.\n");
        return false;
    }
    // If no syncmer in the buffer, compute again
    if(si->cache_index >= si->cache_count){
        // printf("Processing chunk. si->cache_index >= si->cache_count.\n");
        process_chunk_and_cache(si);
    }

    // If no syncmer, process until either it finds more or it consumes the sequence
    while(si->cache_count == 0){
        // If no more s-mers, signal computation ended
        if( si->smers_remaining == 0){
            si->is_done = true;
            return false;
        }
        // printf("Processing chunk. si->cache_count == 0.\n");
        process_chunk_and_cache(si);
    }

    // return one of the syncmer in the buffer
    *syncmer = si->cached_syncmer[si->cache_index++];

    return true;
}

// ============================================================================
// CIRCULAR ARRAY DATA STRUCTURE
// ============================================================================

/*---- circular array structure to handle the hash value of the s-mers in a sliding window ----*/
typedef struct {
    size_t current_position ; // keeps track of the current position in the circula array (for window scanning purposes)
    size_t window_size ; // the number of s-mers in the k-mers (K-S+1)
    size_t buffer_size ; // power of 2 >= window_size for fast modulo
    size_t mask ; // buffer_size - 1, for bitwise AND instead of modulo
    U64 minimum ;
    size_t minimum_position ;
    size_t rescan_count ; // Track number of rescans for performance debugging
    size_t consecutive_rescan_count ; // Track rescans that will trigger immediate next rescan
    U64 hashVector[] ; // contains the hashes of the s-mers for the length of a window w
} CircularArray ;

CircularArray *circularArrayCreate(size_t window_size) {
    // Round up to next power of 2 for fast modulo via bitmask
    size_t buffer_size = 1;
    while (buffer_size < window_size) buffer_size <<= 1;

    CircularArray *ca = (CircularArray *)malloc(sizeof(CircularArray) + buffer_size * sizeof(U64)) ;

    if (ca == NULL){
        fprintf(stderr, "CANNOT INITIALIZE CIRCULAR ARRAY.") ;
        exit (-1) ;
    }

    ca->current_position = 0 ;
    ca->window_size = window_size ;
    ca->buffer_size = buffer_size ;
    ca->mask = buffer_size - 1 ;  // e.g., 32 -> 31 (0x1F)
    ca->minimum = U64MAX;
    ca->minimum_position = window_size + 1;
    ca->rescan_count = 0;
    ca->consecutive_rescan_count = 0;
    return ca ;
}

static void circularArrayDestroy (CircularArray *ca) {if (ca != NULL) free(ca);}

// print array status
void print_status(CircularArray *ca){
    printf("[") ;
    for (size_t i = 0 ; i < ca->window_size ; i++ ) {
        printf("%llu,", ca->hashVector[i]) ;
    }
    printf("]\n") ;
    printf("MIN: %llu ; MIN_P: %lu ; CURR_P : %lu ; WS: %lu\n", ca->minimum, ca->minimum_position, ca->current_position, ca->window_size) ;
}


static inline void update_minimum_branchless(U64 candidate, U64 *current_minimum, size_t candidate_position, size_t *current_minimum_position)
{
    U64 mask = -(candidate < *current_minimum);
    *current_minimum = (*current_minimum & ~mask) | (candidate & mask);
    *current_minimum_position  = (*current_minimum_position & ~mask) | (candidate_position & mask);
}


/*---- perform a re-scan of the entire array when the current minimum is out of context and return the min and position ----*/
void circularScanBranchless(CircularArray *ca){
    // Find minimum using strict < (leftmost wins ties, matching NAIVE behavior)
    ca->minimum = U64MAX;
    ca->minimum_position = 0;

    // Scan positions that will remain after expiring position is gone
    for (size_t i = 1; i < ca->window_size; i++) {
        size_t scan_pos = (ca->current_position + ca->buffer_size - ca->window_size + i) & ca->mask;
        update_minimum_branchless(ca->hashVector[scan_pos], &ca->minimum, scan_pos, &ca->minimum_position);
    }
}

/*---- insert a new element in the circular array----*/
void circularInsertBranchless(CircularArray *ca, U64 value) {
    // Check if minimum is at the position about to expire from the window
    size_t expiring = (ca->current_position + ca->buffer_size - ca->window_size) & ca->mask;
    if (ca->minimum_position == expiring) { circularScanBranchless(ca) ; }

    ca->hashVector[ca->current_position] = value;
    update_minimum_branchless(value, &ca->minimum, ca->current_position, &ca->minimum_position) ;
    ca->current_position = (ca->current_position + 1) & ca->mask;  // Fast wrap via bitmask
}

/*---- perform a re-scan of the entire array when the current minimum is out of context and return the min and position ----*/
void circularScan(CircularArray *ca){
    ca->rescan_count++;  // Track number of rescans

    // Find minimum using strict < (leftmost wins ties, matching NAIVE behavior)
    U64 min_val = U64MAX;
    size_t min_pos = 0;

    // Scan positions that will remain after expiring position is gone
    // These are: (current_position - window_size + 1) through (current_position - 1)
    for (size_t i = 1; i < ca->window_size; i++) {
        size_t scan_pos = (ca->current_position + ca->buffer_size - ca->window_size + i) & ca->mask;
        if (ca->hashVector[scan_pos] < min_val) {  // < for leftmost wins
            min_val = ca->hashVector[scan_pos];
            min_pos = scan_pos;
        }
    }

    ca->minimum = min_val;
    ca->minimum_position = min_pos;

    // Debug: Check if we're setting up for an immediate rescan
    size_t next_expiring = (ca->current_position + ca->buffer_size - ca->window_size + 1) & ca->mask;
    if (ca->minimum_position == next_expiring) {
        ca->consecutive_rescan_count++;
    }
}

/*---- insert a new element in the circular array----*/
void circularInsert(CircularArray *ca, U64 value) {
    // Check if minimum is at the position about to expire from the window
    size_t expiring = (ca->current_position + ca->buffer_size - ca->window_size) & ca->mask;
    if (ca->minimum_position == expiring) { circularScan(ca) ; }

    ca->hashVector[ca->current_position] = value;
    if (value < ca->minimum){
        ca->minimum = value ;
        ca->minimum_position = ca->current_position ;
    }
    ca->current_position = (ca->current_position + 1) & ca->mask;  // Fast wrap via bitmask
}

bool is_syncmer(CircularArray *ca, size_t *position){
    // After insert, current_position has been incremented
    // oldest = (current_position - window_size) & mask
    // newest = (current_position - 1) & mask
    size_t oldest = (ca->current_position + ca->buffer_size - ca->window_size) & ca->mask;
    size_t newest = (ca->current_position + ca->buffer_size - 1) & ca->mask;

    if (ca->minimum_position == oldest) {
        *position = 0;
        return true;
    }
    else if (ca->minimum_position == newest) {
        *position = ca->window_size - 1;
        return true;
    }
    return false;
}

U64 get_smer(CircularArray *ca){
    return ca->hashVector[ca->minimum_position];
}

// ============================================================================
// CIRCULAR ARRAY SYNCMER ITERATOR
// ============================================================================

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
    if (sh == NULL){ fprintf(stderr, "Error initializing the seqhash structure.\n") ; exit (-1) ; }

    syncit->si = seqhashIterator(sh, sequence, sequence_length) ;
    if (syncit->si == NULL){ fprintf(stderr, "Error initializing the seqhashiterator structure.\n") ; exit (-1) ; }

    syncit->ca = circularArrayCreate( window_size) ;
    if (syncit->ca == NULL){ fprintf(stderr, "Error initializing the circularArray structure.\n") ; exit (-1) ; }

    syncit->current_smer = 0 ;
    syncit->current_min_smer = 0 ;

    syncit->smer_position_in_kmer = 0 ;
    syncit->current_min_smer_position = 0 ;
    syncit->current_kmer_position = 0 ;

    syncit->isDone = false;

    // 1 - precompute 1st window
    size_t precompute = 0 ;
    bool current_orientation;
    while(precompute < window_size - 1){
        bool x = seqhashNext(syncit->si, &syncit->current_smer, &current_orientation);

        if(!x){
            exit(1) ;
        }

        circularInsert(syncit->ca,syncit->current_smer);
        precompute++;
    }

    // 2 - run until end of given sequence
    while(seqhashNext(syncit->si, &syncit->current_smer, &current_orientation)){

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
    bool current_orientation;
    // compute next syncmer
    while(seqhashNext(syncit->si, &syncit->current_smer, &current_orientation)){

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

// ============================================================================
// SYNCMER ALGORITHM IMPLEMENTATIONS
// ============================================================================

typedef struct
{
    U64 minimum ;
    size_t minimum_position ;
    size_t absolute_kmer_position ;
} Syncmer;


#endif // LEGACY_INFRASTRUCTURE_H
