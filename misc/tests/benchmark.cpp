// Performance benchmarks for csyncmer_fast

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <vector>
#include <array>

#include "../../csyncmer_fast.h"
#include "../code/syncmer_seqhash.hpp"
#include "../code/syncmer_nthash32.hpp"
#include "../code/syncmer_nthash128.hpp"
#include "../syng/syng_syncmers.h"
#include "fasta_reader.h"

using namespace csyncmer_nthash32;

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
// FULL BENCHMARK
// ============================================================================

int run_full_benchmark(char *fasta_filename, int K, int S, char *output_file){

    FILE *seqFile;
    seqFile = fopen(fasta_filename, "r");
    if (seqFile == NULL) {
        fprintf(stderr, "Error: Cannot open file '%s'\n", fasta_filename);
        return 1;
    }
    stream *seqStream = stream_open_fasta(seqFile);
    if (seqStream == NULL) {
        fprintf(stderr, "Error: Cannot create stream for '%s'\n", fasta_filename);
        fclose(seqFile);
        return 1;
    }
    char *sequence_input = read_sequence(seqStream);
    if (sequence_input == NULL) {
        fprintf(stderr, "Error: Cannot read sequence from '%s'\n", fasta_filename);
        stream_close(seqStream);
        fclose(seqFile);
        return 1;
    }
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
        fprintf(filePtr, "HASH_SYNGH\t"
                         "HASH_NTH128\t"
                         "SYNCMER_NTH128_ITER\t"
                         "SYNCMER_NTH128_DEQUE\t"
                         "SYNCMER_SYNGH_NAIVE\t"
                         "SYNCMER_SYNGH_DEQUE\t"
                         "SYNCMER_SYNGH_SYNG_ORIG\t"
                         "SYNCMER_SYNGH_RESCAN_ARRAY\t"
                         "SYNCMER_SYNGH_RESCAN_ITER\t"
                         "SYNCMER_SYNGH_RESCAN_BRANCHLESS\t"
                         "SYNCMER_SYNGH_RESCAN_CA\t"
                         "SYNCMER_SYNGH_RESCAN_CA_ITER\t"
                         "HASH_NTH32_2BIT\t"
                         "HASH_NTH32_DIRECT\t"
                         "HASH_NTH32_SIMD\t"
                         "SYNCMER_NTH32_RESCAN_DIRECT\t"
                         "SYNCMER_NTH32_RESCAN_2BIT\t"
                         "SYNCMER_NTH32_DEQUE_2BIT\t"
                         "SYNCMER_NTH32_DEQUE_FUSED\t"
                         "SYNCMER_NTH32_RESCAN_FUSED\t"
                         "SYNCMER_NTH32_RESCAN_SIMD\t"
                         "SYNCMER_NTH32_VANHERK\t"
                         "SYNCMER_NTH64_ITER_POS\t"
                         "SYNCMER_NTH64_CANON_POS\t"
                         "SYNCMER_NTH32_TWOSTACK_COUNT\t"
                         "SYNCMER_NTH32_TWOSTACK_POS\t"
                         "SYNCMER_NTH32_CANON_TWOSTACK_COUNT\t"
                         "SYNCMER_NTH32_CANON_TWOSTACK_POS\n");
    }

    printf("[[HASHING syng]]\n");
    start_time = clock();
    hashing_speed_benchmark(sequence_input, sequence_input_length, K, S);
    end_time = clock();
    print_benchmark(hashing_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[HASHING nth128]]\n");
    start_time = clock();
    nthash_speed_benchmark(sequence_input, S);
    end_time = clock();
    print_benchmark(nt_hashing_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS nth128 iterator]]\n");
    start_time = clock();
    benchmark_syncmer_generator_nthash(sequence_input, K, S);
    end_time = clock();
    print_benchmark(nt_hashing_generator, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS nth128 deque]]\n");
    start_time = clock();
    csyncmer_nthash128_deque(sequence_input, sequence_input_length, K, S);
    end_time = clock();
    print_benchmark(nt_hashing_deque, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng naive]]\n");
    start_time = clock();
    csyncmer_seqhash_naive(sequence_input, sequence_input_length, K, S, &num_syncmer_naive);
    end_time = clock();
    print_benchmark(naive_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng deque]]\n");
    start_time = clock();
    csyncmer_seqhash_deque(sequence_input, sequence_input_length, K, S, &num_syncmer_deque);
    end_time = clock();
    print_benchmark(deque_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng original]]\n");
    start_time = clock();
    compute_closed_syncmers_syng_original(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(syng_original_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng rescan-array]]\n");
    start_time = clock();
    csyncmer_seqhash_rescan_array(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(rescan_name2, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng rescan-iterator]]\n");
    start_time = clock();
    benchmark_syncmer_generator_syng(sequence_input, sequence_input_length, K, S);
    end_time = clock();
    print_benchmark(syng_hashing_generator, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng rescan-branchless]]\n");
    start_time = clock();
    csyncmer_seqhash_rescan_branchless(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(branchless_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng rescan-circular]]\n");
    start_time = clock();
    csyncmer_seqhash_rescan_circular(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan);
    end_time = clock();
    print_benchmark(rescan_name, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS syng rescan-circular-iter]]\n");
    start_time = clock();
    csyncmer_seqhash_rescan_iterator(sequence_input, sequence_input_length, K, S, &num_syncmer_rescan_iterator);
    end_time = clock();
    print_benchmark(rescan_name2, start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    // ntHash32 Benchmarks
    size_t num_nthash32_syncmer;

    printf("[[HASHING nth32 scalar]]\n");
    start_time = clock();
    benchmark_nthash32_2bit(sequence_input, S);
    end_time = clock();
    print_benchmark("NTH32_SCALAR", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[HASHING nth32 direct]]\n");
    start_time = clock();
    benchmark_nthash32_direct(sequence_input, S);
    end_time = clock();
    print_benchmark("NTH32_DIRECT", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[HASHING nth32 simd]]\n");
    start_time = clock();
    benchmark_nthash32_simd(sequence_input, S);
    end_time = clock();
    print_benchmark("NTH32_SIMD", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS nth32 direct-rescan]]\n");
    start_time = clock();
    csyncmer_nthash32_direct_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_DIRECT_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS nth32 2bit-rescan]]\n");
    start_time = clock();
    csyncmer_nthash32_2bit_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_2BIT_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS nth32 2bit-deque]]\n");
    start_time = clock();
    csyncmer_nthash32_2bit_deque(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_2BIT_DEQUE", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS nth32 fused-deque]]\n");
    start_time = clock();
    csyncmer_nthash32_fused_deque(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_FUSED_DEQUE", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[SYNCMERS nth32 fused-rescan]]\n");
    start_time = clock();
    csyncmer_nthash32_fused_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_FUSED_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

#if defined(__AVX2__)
    printf("[[SYNCMERS nth32 simd-rescan]]\n");
    start_time = clock();
    csyncmer_nthash32_simd_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_SIMD_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }
#else
    printf("[[SYNCMERS nth32 simd-rescan]] (skipped - no AVX2)\n");
    if (filePtr != NULL) { fprintf(filePtr, "0\t"); }
#endif

    printf("[[SYNCMERS nth32 van-herk]]\n");
    start_time = clock();
    csyncmer_nthash32_vanherk(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_VANHERK", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    // ntHash64 iterator (scalar, portable, exact)
    printf("[[SYNCMERS nth64 iterator-pos]]\n");
    start_time = clock();
    {
        size_t num_iter = 0;
        size_t pos;
        CsyncmerIterator64* it = csyncmer_iterator_create_64(sequence_input, sequence_input_length, K, S);
        if (it) {
            while (csyncmer_iterator_next_64(it, &pos)) num_iter++;
            csyncmer_iterator_destroy_64(it);
        }
        printf("[NTH64_ITERATOR_POS]:: COMPUTED %zu CLOSED SYNCMERS\n", num_iter);
    }
    end_time = clock();
    print_benchmark("NTH64_ITERATOR_POS", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    // ntHash64 canonical iterator (strand-independent, scalar, exact)
    printf("[[SYNCMERS nth64 canonical-pos]]\n");
    start_time = clock();
    {
        size_t num_iter = 0;
        size_t pos;
        int strand;
        CsyncmerIteratorCanonical64* it = csyncmer_iterator_create_canonical_64(sequence_input, sequence_input_length, K, S);
        if (it) {
            while (csyncmer_iterator_next_canonical_64(it, &pos, &strand)) num_iter++;
            csyncmer_iterator_destroy_canonical_64(it);
        }
        printf("[NTH64_CANONICAL_POS]:: COMPUTED %zu CLOSED SYNCMERS\n", num_iter);
    }
    end_time = clock();
    print_benchmark("NTH64_CANONICAL_POS", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

#if defined(__AVX2__)
    // TWOSTACK SIMD count-only
    printf("[[SYNCMERS nth32 twostack-count]]\n");
    start_time = clock();
    {
        volatile size_t count = csyncmer_twostack_simd_32_count(sequence_input, sequence_input_length, K, S);
        printf("[NTH32_TWOSTACK_COUNT]:: COMPUTED %zu CLOSED SYNCMERS\n", (size_t)count);
    }
    end_time = clock();
    print_benchmark("NTH32_TWOSTACK_COUNT", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    // TWOSTACK SIMD with position collection
    printf("[[SYNCMERS nth32 twostack-pos]]\n");
    {
        size_t max_positions = sequence_input_length;
        uint32_t* positions = (uint32_t*)aligned_alloc(32, max_positions * sizeof(uint32_t));
        if (positions) {
            start_time = clock();
            size_t count = csyncmer_twostack_simd_32_positions(sequence_input, sequence_input_length, K, S, positions, max_positions);
            end_time = clock();
            printf("[NTH32_TWOSTACK_POS]:: COMPUTED %zu CLOSED SYNCMERS\n", count);
            print_benchmark("NTH32_TWOSTACK_POS", start_time, end_time, fasta_filename, filePtr);
            free(positions);
        }
    }
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    // Canonical TWOSTACK SIMD count-only
    printf("[[SYNCMERS nth32 canonical-twostack-count]]\n");
    start_time = clock();
    {
        volatile size_t count = csyncmer_twostack_simd_32_canonical_count(sequence_input, sequence_input_length, K, S);
        printf("[NTH32_CANON_TWOSTACK_COUNT]:: COMPUTED %zu CLOSED SYNCMERS\n", (size_t)count);
    }
    end_time = clock();
    print_benchmark("NTH32_CANON_TWOSTACK_COUNT", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    // Canonical TWOSTACK SIMD with position and strand collection
    printf("[[SYNCMERS nth32 canonical-twostack-pos]]\n");
    {
        size_t max_positions = sequence_input_length;
        uint32_t* positions = (uint32_t*)aligned_alloc(32, max_positions * sizeof(uint32_t));
        uint8_t* strands = (uint8_t*)aligned_alloc(32, max_positions);
        if (positions && strands) {
            start_time = clock();
            size_t count = csyncmer_twostack_simd_32_canonical_positions(
                sequence_input, sequence_input_length, K, S, positions, strands, max_positions);
            end_time = clock();
            printf("[NTH32_CANON_TWOSTACK_POS]:: COMPUTED %zu CLOSED SYNCMERS\n", count);
            print_benchmark("NTH32_CANON_TWOSTACK_POS", start_time, end_time, fasta_filename, filePtr);
        }
        free(positions);
        free(strands);
    }
#endif

    if (filePtr != NULL) { fprintf(filePtr, "\n"); }

    fclose(filePtr);
    free(encoded_seq);
    free((void*)sequence_input);
    return 0;
}

// ============================================================================
// QUICK BENCHMARK
// ============================================================================

int run_quick_benchmark(char *fasta_filename, int K, int S, const char *filter){

    FILE *seqFile = fopen(fasta_filename, "r");
    if (!seqFile) { fprintf(stderr, "Error: Cannot open file\n"); return 1; }
    stream *seqStream = stream_open_fasta(seqFile);
    char *sequence = read_sequence(seqStream);
    size_t length = strlen(sequence);
    stream_close(seqStream);
    fclose(seqFile);

    off_t file_size = get_file_size(fasta_filename);
    double file_size_mb = file_size / (1024.0 * 1024.0);
    size_t num;

    bool run_32 = (strcmp(filter, "32") == 0 || strcmp(filter, "all") == 0);
    bool run_64 = (strcmp(filter, "64") == 0 || strcmp(filter, "all") == 0);

    printf("%-25s %10s %12s\n", "Implementation", "Syncmers", "Speed (MB/s)");
    printf("%-25s %10s %12s\n", "-------------------------", "----------", "------------");

    if (run_64) {
        // Iterator benchmark (scalar, portable, exact)
        {
            clock_t start = clock();
            num = 0;
            size_t pos;
            CsyncmerIterator64* it = csyncmer_iterator_create_64(sequence, length, K, S);
            if (it) {
                while (csyncmer_iterator_next_64(it, &pos)) num++;
                csyncmer_iterator_destroy_64(it);
            }
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH64_ITERATOR_POS", num, file_size_mb / elapsed);
        }

        // Canonical iterator benchmark (strand-independent, scalar, exact)
        {
            clock_t start = clock();
            num = 0;
            size_t pos;
            int strand;
            CsyncmerIteratorCanonical64* it = csyncmer_iterator_create_canonical_64(sequence, length, K, S);
            if (it) {
                while (csyncmer_iterator_next_canonical_64(it, &pos, &strand)) num++;
                csyncmer_iterator_destroy_canonical_64(it);
            }
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH64_CANONICAL_POS", num, file_size_mb / elapsed);
        }
    }

    if (run_32) {
#if defined(__AVX2__)
        // TWOSTACK count-only (no position collection)
        {
            clock_t start = clock();
            size_t count = csyncmer_twostack_simd_32_count(sequence, length, K, S);
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH32_TWOSTACK_COUNT", count, file_size_mb / elapsed);
        }

        // TWOSTACK with position collection
        size_t max_positions = length;
        uint32_t* positions = (uint32_t*)aligned_alloc(32, max_positions * sizeof(uint32_t));
        if (positions) {
            clock_t start = clock();
            size_t count = csyncmer_twostack_simd_32_positions(sequence, length, K, S, positions, max_positions);
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH32_TWOSTACK_POS", count, file_size_mb / elapsed);
            free(positions);
        }

        // Canonical TWOSTACK count-only
        {
            clock_t start = clock();
            size_t count = csyncmer_twostack_simd_32_canonical_count(sequence, length, K, S);
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH32_CANON_TWOSTACK_CNT", count, file_size_mb / elapsed);
        }

        // Canonical TWOSTACK with position and strand collection
        uint8_t* strands = (uint8_t*)aligned_alloc(32, max_positions);
        if (positions && strands) {
            positions = (uint32_t*)aligned_alloc(32, max_positions * sizeof(uint32_t));
            clock_t start = clock();
            size_t count = csyncmer_twostack_simd_32_canonical_positions(
                sequence, length, K, S, positions, strands, max_positions);
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH32_CANON_TWOSTACK_POS", count, file_size_mb / elapsed);
            free(positions);
            free(strands);
        }
#endif
    }

    free((void*)sequence);
    return 0;
}

// ============================================================================
// MAIN
// ============================================================================

void print_usage(const char *prog_name) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s FASTA_FILE K S OUT_FILE        Full benchmark suite\n", prog_name);
    fprintf(stderr, "  %s --quick FASTA_FILE K S [filter] Quick benchmark (filter: 32, 64, all)\n", prog_name);
}

int main(int argc, char *argv[]) {

    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Check for --quick flag
    if (strcmp(argv[1], "--quick") == 0) {
        if (argc < 5) {
            fprintf(stderr, "Usage: %s --quick FASTA_FILE K S [filter]\n", argv[0]);
            fprintf(stderr, "filter: 32, 64, all (default: all)\n");
            return 1;
        }
        char *fasta_file = argv[2];
        int K = atoi(argv[3]);
        int S = atoi(argv[4]);
        const char *filter = (argc > 5) ? argv[5] : "all";

        if (S >= K) {
            fprintf(stderr, "Error: S (%d) must be less than K (%d)\n", S, K);
            return 1;
        }

        return run_quick_benchmark(fasta_file, K, S, filter);
    }

    // Default: full benchmark
    if (argc < 5) {
        fprintf(stderr, "Usage: %s FASTA_FILE K S OUT_FILE\n", argv[0]);
        return 1;
    }

    char *fasta_file = argv[1];
    int K = atoi(argv[2]);
    int S = atoi(argv[3]);
    char *output_file = argv[4];

    if (S >= K) {
        fprintf(stderr, "Error: S (%d) must be less than K (%d)\n", S, K);
        return 1;
    }

    return run_full_benchmark(fasta_file, K, S, output_file);
}
