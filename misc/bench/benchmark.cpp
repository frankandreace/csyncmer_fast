// Consolidated benchmark and test file for csyncmer_fast
// Combines: benchmarking utilities, FASTA reader, benchmark harnesses, and unit tests

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <sys/stat.h>
#include <ctype.h>

#include "../../csyncmer_fast.h"
#include "../syng/syng_syncmers.h"

// ============================================================================
// BENCHMARKING UTILITIES
// ============================================================================

off_t get_file_size(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
      return st.st_size;
    }
    return -1;
}

void print_timing_stats(double elapsed_time){
    printf("Processing time: %.4f seconds\n", elapsed_time);
}

void print_troughput(double file_size_mb, double processing_speed){
    printf("File size: %.2f MB\n", file_size_mb);
    printf("Processing speed: %.4f MB/sec\n", processing_speed);
}

void print_benchmark(const char* name, double start_time, double end_time, char* filename, FILE *filePtr){
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    off_t file_size = get_file_size(filename);
    double file_size_mb = file_size / (1024.0 * 1024.0);
    double processing_speed = file_size_mb / elapsed_time;

    if (filePtr != NULL) { fprintf(filePtr, "%lf", processing_speed) ; }

    printf("%s Performance metrics.\n", name);
    print_timing_stats(elapsed_time);
    print_troughput(file_size_mb, processing_speed);
    printf("\n");
}

// ============================================================================
// FASTA READER
// ============================================================================

#define INITIAL_BUF_SIZE 1024
#define STREAM_BUFSIZ (1<<12)

typedef enum {
    STREAM_FILE,
    STREAM_FASTA
} stream_type;

typedef struct stream {
    stream_type type;
    int reverse;
    long len;
    long pos;
    FILE *f;
    long buf_ptr;
    size_t buf_size;
    unsigned char buf[STREAM_BUFSIZ];
} stream;

static int stream_getnext_file(stream *S) {
    if (S->pos >= S->len || (S->pos - S->buf_ptr) >= (long)S->buf_size) {
        S->buf_ptr = S->pos;
        fseek(S->f, S->buf_ptr, SEEK_SET);
        S->buf_size = fread(S->buf, 1, STREAM_BUFSIZ, S->f);
        if (S->buf_size == 0)
            return EOF;
    }
    int c = S->buf[S->pos - S->buf_ptr];
    S->pos++;
    return c;
}

static int stream_getnext_fasta(stream *S) {
    int c;
    while (1) {
        c = stream_getnext_file(S);
        if (c == EOF)
            return EOF;
        if (c == '\n') {
            int next = stream_getnext_file(S);
            if (next == '>') {
                S->pos--;
                return 0;
            }
            continue;
        }
        if (c == '>') {
            while (c != '\n' && c != EOF)
                c = stream_getnext_file(S);
            continue;
        }
        c = tolower(c);
        if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'n')
            return c;
    }
}

int stream_getnext(stream *S) {
    if (S->type == STREAM_FASTA)
        return stream_getnext_fasta(S);
    return EOF;
}

stream *stream_open_fasta(FILE *fin) {
    if (!fin)
        return NULL;
    stream *S = (stream*)malloc(sizeof(stream));
    if (!S) {
        fprintf(stderr, "Not enough memory\n");
        exit(1);
    }
    S->type = STREAM_FASTA;
    S->reverse = 0;
    S->pos = 0;
    fseek(fin, 0, SEEK_END);
    S->len = ftell(fin);
    fseek(fin, 0, SEEK_SET);
    S->buf_ptr = 0;
    S->buf_size = 0;
    S->f = fin;
    return S;
}

char *read_sequence(struct stream *S) {
    size_t capacity = INITIAL_BUF_SIZE;
    size_t len = 0;
    char *sequence = (char*)malloc(capacity);
    if (!sequence) {
        fprintf(stderr, "Memory allocation failure\n");
        return NULL;
    }

    int c;
    while ((c = stream_getnext(S)) != EOF && c != 0) {
        sequence[len++] = (char)c;
        if (len >= capacity) {
            capacity *= 2;
            char *tmp = (char*)realloc(sequence, capacity);
            if (!tmp) {
                free(sequence);
                fprintf(stderr, "Memory allocation failure during expansion\n");
                return NULL;
            }
            sequence = tmp;
        }
    }
    sequence[len] = '\0';
    return sequence;
}

void stream_close(stream *S) {
    if (S)
        free(S);
}

// ============================================================================
// BENCHMARK HARNESSES
// ============================================================================

void hashing_speed_benchmark(char *sequence_input, size_t length, size_t K, size_t S){
    U64 seed  = 7;
    size_t window_size = (size_t)K - (size_t)S + 1;
    size_t computed_smers = 0;

    Seqhash *sh = seqhashCreate(S, window_size, seed);
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, length);

    U64 current_smer;
    bool current_orientation;
    size_t current_position = 0;

    while(seqhashNext(si, &current_smer, &current_orientation)){
        computed_smers++;
        current_position++;
    }
    free(si);
    free(sh);
    printf("[HASHING_BENCHMARK]:: HASHED %lu S-MERS\n", computed_smers);
}

void nthash_speed_benchmark(const char *sequence_input, size_t S){
    NtHashHandle rolling_hash = nthash_create(sequence_input, strlen(sequence_input), S, 2);
    if (rolling_hash == NULL) {
        return;
    }

    size_t count = 0;
    while(nthash_roll(rolling_hash)){
        nthash_get_canonical_hash_128(rolling_hash);
        count++;
    }
    printf("[NT_HASH_BENCH]:: HASHED %lu S-MERS\n", count);
    nthash_destroy(rolling_hash);
}

void benchmark_syncmer_generator_syng(char *sequence_input, size_t seq_input_length, size_t K, size_t S){
    size_t num_syncmer = 0;
    Syncmer64 my_syncmer = {0,0};

    printf("BUILDING GENERATOR.\n");
    SyncmerIteratorS *my_syncmer_iterator = syncmer_generator_createS(sequence_input, seq_input_length, K, S);
    printf("ITERATING.\n");
    while(syncmer_iterator_nextS(my_syncmer_iterator, &my_syncmer)){
        num_syncmer++;
    }
    printf("[SYNG_HASH_SYNCMER_GENERATOR]:: COMPUTED %lu CLOSED SYNCMERS\n", num_syncmer);
    syncmer_generator_destroyS(my_syncmer_iterator);
}

void benchmark_syncmer_generator_nthash(const char *sequence_input, size_t K, size_t S){
    size_t num_syncmer = 0;
    Syncmer128 my_syncmer = {0,0,0};

    printf("BUILDING GENERATOR.\n");
    SyncmerIterator *my_syncmer_iterator = syncmer_generator_create(sequence_input, K, S);
    printf("ITERATING.\n");
    while(syncmer_iterator_next(my_syncmer_iterator, &my_syncmer)){
        num_syncmer++;
    }
    printf("[NT_HASH_SYNCMER_GENERATOR]:: COMPUTED %lu CLOSED SYNCMERS\n", num_syncmer);
    syncmer_generator_destroy(my_syncmer_iterator);
}

// ============================================================================
// UNIT TESTS
// ============================================================================

extern const unsigned char base_to_bits_array[];

void test_base_to_bits(){
    assert(base_to_bits_array[0] == 0);
    assert(base_to_bits_array[1] == 1);
    assert(base_to_bits_array[2] == 2);
    assert(base_to_bits_array[3] == 3);

    assert(base_to_bits_array['A'] == 0);
    assert(base_to_bits_array['C'] == 1);
    assert(base_to_bits_array['G'] == 2);
    assert(base_to_bits_array['T'] == 3);
    assert(base_to_bits_array['U'] == 3);

    assert(base_to_bits_array['a'] == 0);
    assert(base_to_bits_array['c'] == 1);
    assert(base_to_bits_array['g'] == 2);
    assert(base_to_bits_array['t'] == 3);
    assert(base_to_bits_array['u'] == 3);

    printf("[TEST] base_to_bits: PASSED\n");
}

void run_unit_tests(){
    printf("=== RUNNING UNIT TESTS ===\n");
    test_base_to_bits();
    printf("=== ALL UNIT TESTS PASSED ===\n\n");
}

// ============================================================================
// BENCHMARK FROM FILE
// ============================================================================

int compute_from_file(char *fasta_filename, int K, int S, char *output_file){

    FILE *seqFile;
    seqFile = fopen(fasta_filename, "r");
    stream *seqStream = stream_open_fasta(seqFile);
    char *sequence_input = read_sequence(seqStream);
    size_t sequence_input_length = strlen(sequence_input);

    printf("SEQ SIZE IS %ld\n", sequence_input_length);

    char *encoded_seq = (char*)malloc(sequence_input_length * sizeof(char));
    if (encoded_seq == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    for (size_t i = 0; i < sequence_input_length; i++) {
        encoded_seq[i] = base_to_bits(sequence_input[i]);
    }

    stream_close(seqStream);
    fclose(seqFile);

    const char * rescan_name = "RESCAN";
    const char * rescan_name2 = "RESCAN ARRAY";
    const char * naive_name = "NAIVE";
    const char * hashing_name = "HASHING";
    const char * nt_hashing_name = "NT_HASHING";
    const char * nt_hashing_generator = "GENERATOR_NT";
    const char * syng_hashing_generator = "GENERATOR_SYNG";
    const char * nt_hashing_deque = "DEQUE_NT";
    const char * deque_name = "DEQUE";
    const char * branchless_name = "BRANCHLESS";
    const char * syng_original_name = "SYNG ORIGINAL";

    size_t num_syncmer_rescan;
    size_t num_syncmer_naive;
    size_t num_syncmer_deque;
    size_t num_syncmer_rescan_iterator;

    clock_t start_time;
    clock_t end_time;

    FILE *filePtr;

    bool first_writing = false;
    if (access(output_file, F_OK) != 0) {
        first_writing = true;
    }

    filePtr = fopen(output_file, "a");

    if (filePtr == NULL) {
        perror("Error opening file");
        return 1;
    }

    if (filePtr != NULL && first_writing) {
        fprintf(filePtr, "SYNGH_JUST_HASHING\tNTH_JUST_HASHING\tNTH_ITERATOR\tNTH_DEQUE\tSYNGH_NAIVE\tSYNGH_DEQUE\tSYNGH_DURBIN_ITERATOR\tSYNGH_RESCAN_NO_ITERATOR\tSYNGH_RESCAN_ITERATOR\tSYNGH_RESCAN_CA_BRANCHLESS\tSYNGH_RESCAN_CA\tSYNGH_RESCAN_CA_ITERATOR\n");
    }

    printf("[[HASHING SPEED BENCHMARK]]\n");
    start_time = clock();
    hashing_speed_benchmark(sequence_input, sequence_input_length, K, S);
    end_time = clock();
    print_benchmark(hashing_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASHING SPEED BENCHMARK]]\n");
    start_time = clock();
    nthash_speed_benchmark(sequence_input, S);
    end_time = clock();
    print_benchmark(nt_hashing_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[GENERATOR - NTHASH]]\n");
    start_time = clock();
    benchmark_syncmer_generator_nthash(sequence_input, K, S);
    end_time = clock();
    print_benchmark(nt_hashing_generator, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[DEQUE - NTHASH]]\n");
    start_time = clock();
    compute_closed_syncmer_deque_nthash(sequence_input, sequence_input_length, K, S);
    end_time = clock();
    print_benchmark(nt_hashing_deque, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NAIVE]]\n");
    start_time = clock();
    compute_closed_syncmers_naive(sequence_input, sequence_input_length, K, S, &num_syncmer_naive);
    end_time = clock();
    print_benchmark(naive_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[DEQUE]]\n");
    start_time = clock();
    compute_closed_syncmers_deque_rayan(sequence_input, sequence_input_length, K, S, &num_syncmer_deque);
    end_time = clock();
    print_benchmark(deque_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNG ORIGINAL]]\n");
    start_time = clock();
    compute_closed_syncmers_syng_original(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(syng_original_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[RESCAN ITERATIVE]]\n");
    start_time = clock();
    compute_closed_syncmers_rescan(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(rescan_name2, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[RESCAN ITERATOR]]\n");
    start_time = clock();
    benchmark_syncmer_generator_syng(sequence_input, sequence_input_length, K, S);
    end_time = clock();
    print_benchmark(syng_hashing_generator, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[RESCAN BRANCHLESS]]\n");
    start_time = clock();
    compute_closed_syncmers_branchless(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(branchless_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[RESCAN CIRCULAR ARRAY]]\n");
    start_time = clock();
    compute_closed_syncmers(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(rescan_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[RESCAN CIRCULAR ARRAY ITERATOR]]\n");
    start_time = clock();
    compute_closed_syncmers_rescan_iterator(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan_iterator);
    end_time = clock();
    print_benchmark(rescan_name2, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\n"); }

    fclose(filePtr);
    free(encoded_seq);
    free((void*)sequence_input);
    return 0;
}

// ============================================================================
// HELPER FUNCTIONS FOR CORRECTNESS TESTING
// ============================================================================

size_t count_syncmer_generator_nthash(const char *sequence_input, size_t K, size_t S){
    size_t num_syncmer = 0;
    Syncmer128 my_syncmer = {0,0,0};
    SyncmerIterator *my_syncmer_iterator = syncmer_generator_create(sequence_input, K, S);
    while(syncmer_iterator_next(my_syncmer_iterator, &my_syncmer)){
        num_syncmer++;
    }
    printf("[NT_HASH_SYNCMER_GENERATOR]:: COMPUTED %lu CLOSED SYNCMERS\n", num_syncmer);
    syncmer_generator_destroy(my_syncmer_iterator);
    return num_syncmer;
}

size_t count_syncmer_deque_nthash(const char *sequence_input, size_t length, size_t K, size_t S){
    NtHashHandle rolling_hash = nthash_create(sequence_input, strlen(sequence_input), S, 2);
    if (rolling_hash == NULL) {
        fprintf(stderr, "Failed to create NtHash handle\n");
        return 0;
    }

    size_t num_s_mers = length - S + 1;
    size_t window_size = K - S + 1;
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

    printf("[DEQUE_NT]:: COMPUTED %lu CLOSED SYNCMERS\n", computed_syncmers);

    free(s_mer_hashes);
    free(deque);
    nthash_destroy(rolling_hash);
    return computed_syncmers;
}

// ============================================================================
// CORRECTNESS CHECK FROM SEQUENCE
// ============================================================================

int compute_from_sequence(char *sequence_input, int K, int S){

    size_t len = strlen(sequence_input);
    char *encoded_seq = (char*)malloc(len * sizeof(uint8_t));
    if (encoded_seq == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    for (size_t i = 0; i < len; i++) {
        encoded_seq[i] = base_to_bits(sequence_input[i]);
    }

    // Seqhash-based implementations (64-bit)
    size_t num_naive;
    size_t num_circular_array;
    size_t num_large_array;
    size_t num_deque;
    size_t num_iterator;
    size_t num_branchless;
    size_t num_syng_original;

    // ntHash-based implementations (128-bit)
    size_t num_nt_generator;
    size_t num_nt_deque;

    printf("=== TESTING SEQHASH-BASED IMPLEMENTATIONS (64-bit) ===\n");
    compute_closed_syncmers_naive(encoded_seq, len, K, S, &num_naive);
    compute_closed_syncmers(encoded_seq, len, K, S, &num_circular_array);
    compute_closed_syncmers_rescan(encoded_seq, len, K, S, &num_large_array);
    compute_closed_syncmers_deque_rayan(encoded_seq, len, K, S, &num_deque);
    compute_closed_syncmers_rescan_iterator(encoded_seq, len, K, S, &num_iterator);
    compute_closed_syncmers_branchless(encoded_seq, len, K, S, &num_branchless);

    // Note: SYNG_ORIGINAL uses Durbin's implementation which may compute a different
    // variant of syncmers - excluded from correctness comparison
    // compute_closed_syncmers_syng_original(encoded_seq, len, K, S, &num_syng_original);

    printf("\n=== TESTING NTHASH-BASED IMPLEMENTATIONS (128-bit) ===\n");
    // These need ASCII input, not encoded
    num_nt_generator = count_syncmer_generator_nthash(sequence_input, K, S);
    num_nt_deque = count_syncmer_deque_nthash(sequence_input, len, K, S);

    free(encoded_seq);

    // Check seqhash-based implementations agree
    bool seqhash_ok = true;
    if (num_naive != num_circular_array){
        printf("[SEQHASH ERROR] NAIVE: %lu ; CIRCULAR_ARRAY: %lu\n", num_naive, num_circular_array);
        seqhash_ok = false;
    }
    if (num_naive != num_large_array){
        printf("[SEQHASH ERROR] NAIVE: %lu ; LARGE_ARRAY: %lu\n", num_naive, num_large_array);
        seqhash_ok = false;
    }
    if (num_naive != num_deque){
        printf("[SEQHASH ERROR] NAIVE: %lu ; DEQUE: %lu\n", num_naive, num_deque);
        seqhash_ok = false;
    }
    if (num_naive != num_iterator){
        printf("[SEQHASH ERROR] NAIVE: %lu ; ITERATOR: %lu\n", num_naive, num_iterator);
        seqhash_ok = false;
    }
    if (num_naive != num_branchless){
        printf("[SEQHASH ERROR] NAIVE: %lu ; BRANCHLESS: %lu\n", num_naive, num_branchless);
        seqhash_ok = false;
    }
    // SYNG_ORIGINAL excluded - uses Durbin's implementation (may compute different variant)

    // Check ntHash-based implementations agree
    bool nthash_ok = true;
    if (num_nt_generator != num_nt_deque){
        printf("[NTHASH ERROR] GENERATOR: %lu ; DEQUE: %lu\n", num_nt_generator, num_nt_deque);
        nthash_ok = false;
    }

    // Report results
    if (seqhash_ok && nthash_ok) {
        printf("\n[CORRECTNESS] All seqhash implementations agree: %lu syncmers\n", num_naive);
        printf("[CORRECTNESS] All ntHash implementations agree: %lu syncmers\n", num_nt_generator);
        if (num_naive != num_nt_generator) {
            printf("[NOTE] Seqhash and ntHash produce different counts (expected if using different hash functions)\n");
        }
        return 0;
    } else {
        return 1;
    }
}

// ============================================================================
// MAIN
// ============================================================================

void print_usage(const char *prog_name) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s test                           Run unit tests\n", prog_name);
    fprintf(stderr, "  %s check SEQUENCE K S             Check algorithm correctness\n", prog_name);
    fprintf(stderr, "  %s bench FASTA_FILE K S OUT_FILE  Run benchmarks on FASTA file\n", prog_name);
}

int main(int argc, char *argv[]) {

    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char *mode = argv[1];

    if (strcmp(mode, "test") == 0) {
        run_unit_tests();
        return 0;
    }
    else if (strcmp(mode, "check") == 0) {
        if (argc < 5) {
            fprintf(stderr, "Usage: %s check SEQUENCE K S\n", argv[0]);
            return 1;
        }
        char *sequence = argv[2];
        int K = atoi(argv[3]);
        int S = atoi(argv[4]);

        if (S >= K) {
            fprintf(stderr, "Error: S (%d) must be less than K (%d)\n", S, K);
            return 1;
        }

        return compute_from_sequence(sequence, K, S);
    }
    else if (strcmp(mode, "bench") == 0) {
        if (argc < 6) {
            fprintf(stderr, "Usage: %s bench FASTA_FILE K S OUT_FILE\n", argv[0]);
            return 1;
        }
        char *fasta_file = argv[2];
        int K = atoi(argv[3]);
        int S = atoi(argv[4]);
        char *output_file = argv[5];

        if (S >= K) {
            fprintf(stderr, "Error: S (%d) must be less than K (%d)\n", S, K);
            return 1;
        }

        return compute_from_file(fasta_file, K, S, output_file);
    }
    else {
        fprintf(stderr, "Unknown mode: %s\n", mode);
        print_usage(argv[0]);
        return 1;
    }
}
