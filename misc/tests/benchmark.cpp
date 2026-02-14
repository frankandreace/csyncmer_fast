// Performance benchmarks for csyncmer_fast

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <vector>
#include <array>

#include "../../csyncmer_fast.h"
#include "../legacy/syncmer_seqhash.hpp"
#include "../legacy/syncmer_nthash32.hpp"
#include "../legacy/syncmer_nthash128.hpp"
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
                         "SYNCMER_NTH32_VANHERK\n");
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
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    // ntHash32 Benchmarks
    size_t num_nthash32_syncmer;

    printf("[[NTHASH32 SCALAR]]\n");
    start_time = clock();
    benchmark_nthash32_2bit(sequence_input, S);
    end_time = clock();
    print_benchmark("NTH32_SCALAR", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASH32 DIRECT]]\n");
    start_time = clock();
    benchmark_nthash32_direct(sequence_input, S);
    end_time = clock();
    print_benchmark("NTH32_DIRECT", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASH32 SIMD]]\n");
    start_time = clock();
    benchmark_nthash32_simd(sequence_input, S);
    end_time = clock();
    print_benchmark("NTH32_SIMD", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASH32 DIRECT RESCAN]]\n");
    start_time = clock();
    compute_closed_syncmers_nthash32_direct_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_DIRECT_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASH32 2BIT RESCAN]]\n");
    start_time = clock();
    compute_closed_syncmers_nthash32_2bit_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_2BIT_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASH32 2BIT DEQUE]]\n");
    start_time = clock();
    compute_closed_syncmers_nthash32_2bit_deque(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_2BIT_DEQUE", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASH32 FUSED DEQUE]]\n");
    start_time = clock();
    compute_closed_syncmers_nthash32_fused_deque(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_FUSED_DEQUE", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

    printf("[[NTHASH32 FUSED RESCAN]]\n");
    start_time = clock();
    compute_closed_syncmers_nthash32_fused_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_FUSED_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }

#if defined(__AVX2__)
    printf("[[NTHASH32 SIMD RESCAN]]\n");
    start_time = clock();
    compute_closed_syncmers_nthash32_simd_rescan(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_SIMD_RESCAN", start_time, end_time, fasta_filename, filePtr);
    if (filePtr != NULL) { fprintf(filePtr, "\t"); }
#else
    printf("[[NTHASH32 SIMD RESCAN]] (skipped - no AVX2)\n");
    if (filePtr != NULL) { fprintf(filePtr, "0\t"); }
#endif

    printf("[[NTHASH32 VAN HERK]]\n");
    start_time = clock();
    compute_closed_syncmers_nthash32_vanherk(sequence_input, sequence_input_length, K, S, &num_nthash32_syncmer);
    end_time = clock();
    print_benchmark("NTH32_VANHERK", start_time, end_time, fasta_filename, filePtr);

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

    if (run_32) {
#if defined(__AVX2__)
        // TWOSTACK count-only (no position collection)
        {
            clock_t start = clock();
            size_t count = csyncmer_compute_twostack_simd_32_count(sequence, length, K, S);
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH32_TWOSTACK_COUNT", count, file_size_mb / elapsed);
        }

        // TWOSTACK with position collection
        size_t max_positions = length;
        uint32_t* positions = (uint32_t*)aligned_alloc(32, max_positions * sizeof(uint32_t));
        if (positions) {
            clock_t start = clock();
            size_t count = csyncmer_compute_twostack_simd_32(sequence, length, K, S, positions, max_positions);
            double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%-25s %10zu %12.2f\n", "NTH32_TWOSTACK_POS", count, file_size_mb / elapsed);
            free(positions);
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
