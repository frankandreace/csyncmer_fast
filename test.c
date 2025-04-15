#include <string.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>

#include "src/implementations.h"
#include "src/benchmarking.h"
#include "src/fasta_reader.h"

#include "syng/syng_syncmers.h"

/*---- conversion of bases (ascii char) into bits ----*/
static inline uint8_t base_to_bits(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0; // Treat Ns and  unknown as 'A'
    }
  }

  int compute_from_file(char *fasta_filename, int K, int S){

    // OPENING FILE AND READING THE SEQUENCE

    FILE *seqFile;
    seqFile = fopen(fasta_filename,"r");
    stream *seqStream = stream_open_fasta(seqFile) ;
    char *sequence_input = read_sequence(seqStream) ;
    int sequence_input_length = strlen(sequence_input) ;

    printf("SEQ SIZE IS %d\n", sequence_input_length) ;

    // ENCODING ASCII INTO 2-BIT AGCT ENCODING
    uint8_t *encoded_seq = malloc(sequence_input_length * sizeof(uint8_t));
    uint8_t *encoded_seq1 = malloc(sequence_input_length * sizeof(uint8_t));
    uint8_t *encoded_seq2 = malloc(sequence_input_length* sizeof(uint8_t));

    if (encoded_seq == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Loop through the sequence and convert each character
    for (size_t i = 0; i < sequence_input_length; i++) {
        encoded_seq[i] = base_to_bits(sequence_input[i]);
        encoded_seq1[i] = base_to_bits(sequence_input[i]);
        encoded_seq2[i] = base_to_bits(sequence_input[i]);
    }

    for (int i = 0; i < 20; i++){
        printf("%c,", sequence_input[i]) ;
    }  
    printf("\n") ;

    for (int i = 0; i < 20; i++){
        printf("%u,", encoded_seq[i]) ;
    }
    printf("\n") ;

    for (int i = 0; i < 20; i++){
        printf("%u,", encoded_seq1[i]) ;
    }
    printf("\n") ;
    for (int i = 0; i < 20; i++){
        printf("%u,", encoded_seq2[i]) ;
    }
    printf("\n") ;

    stream_close(seqStream);


    // BENCHMARKING TIME
    char * rescan_name = "RESCAN" ;
    char * rescan_name2 = "RESCAN ARRAY" ;
    char * naive_name = "NAIVE" ;
    char * hashing_name = "HASHING" ;
    char * deque_name = "DEQUE" ;
    char * branchless_name = "BRANCHLESS" ;
    char * syng_original_name = "SYNG ORIGINAL" ;

    size_t num_syncmer_rescan;
    size_t num_syncmer_naive;
    size_t num_syncmer_deque;
    size_t num_syncmer_rescan_iterator;

    clock_t start_time;
    clock_t end_time;

    //BENCHMARK FILE
    FILE *filePtr;
    char* filename = "benchmark/benchmark.tsv" ;

    bool first_writing = false;
    if (access(filename, F_OK) != 0) {
        first_writing = true;
    }

    filePtr = fopen(filename, "a");

    if (filePtr == NULL) {
        // If the file could not be opened, print an error message and exit
        perror("Error opening file");
        return 1;
    }

    if (filePtr != NULL && first_writing) { 
        fprintf(filePtr, "HASHING\tNAIVE\tDEQUE\tSYNG_ORIGINAL\tRESCAN\tBRANCHLESS\tRESCAN_CIRCULAR_ARRAY\tRESCAN_CA_ITERATOR\n") ; 
    }

    //benchmark speed for just hashing
    start_time = clock();
    hahsing_speed_benchmark(encoded_seq2, sequence_input_length, K, S) ;
    end_time = clock();
    print_benchmark(hashing_name, start_time, end_time, fasta_filename, filePtr) ;
    if (filePtr != NULL) { fprintf(filePtr, "\t") ; }

    //benchmark speed for naive
    start_time = clock();
    compute_closed_syncmers_naive(encoded_seq1, sequence_input_length, K, S, &num_syncmer_naive) ;
    end_time = clock();
    print_benchmark(naive_name, start_time, end_time, fasta_filename, filePtr) ;
    if (filePtr != NULL) { fprintf(filePtr, "\t") ; }

    //benchmark speed for deque
    start_time = clock();
    compute_closed_syncmers_deque_rayan(encoded_seq, sequence_input_length, K, S, &num_syncmer_deque);
    end_time = clock();
    print_benchmark(deque_name, start_time, end_time, fasta_filename, filePtr) ;
    if (filePtr != NULL) { fprintf(filePtr, "\n") ; }

    //benchmark speed for syng implementation
    start_time = clock();
    compute_closed_syncmers_syng_original(encoded_seq, sequence_input_length, K, S, &num_syncmer_rescan) ;
    end_time = clock();
    print_benchmark(syng_original_name, start_time, end_time, fasta_filename, filePtr) ;
    if (filePtr != NULL) { fprintf(filePtr, "\t") ; }

    //benchmark speed for rescan without circular array
    start_time = clock();
    compute_closed_syncmers_rescan(encoded_seq, sequence_input_length, K, S, &num_syncmer_rescan) ;
    end_time = clock();
    print_benchmark(rescan_name2, start_time, end_time, fasta_filename, filePtr) ;
    if (filePtr != NULL) { fprintf(filePtr, "\t") ; }

    // //benchmark speed for branchless rescan
    // start_time = clock();
    // compute_closed_syncmers_rescan_branchless(encoded_seq, sequence_input_length, K, S, &num_syncmer_rescan) ;
    // end_time = clock();
    // print_benchmark(branchless_name, start_time, end_time, fasta_filename, filePtr) ;
    // if (filePtr != NULL) { fprintf(filePtr, "\t") ; }

    //benchmark speed for syncmer with rescan and ciruclar array
    start_time = clock();
    compute_closed_syncmers(encoded_seq, sequence_input_length, K, S, &num_syncmer_rescan) ;
    end_time = clock();
    print_benchmark(rescan_name, start_time, end_time, fasta_filename, filePtr) ;
    if (filePtr != NULL) { fprintf(filePtr, "\t") ; }

    //benchmark speed for rescan with circular array and iterator
    start_time = clock();
    compute_closed_syncmers_rescan_iterator(encoded_seq, sequence_input_length, K, S, &num_syncmer_rescan_iterator) ;
    end_time = clock();
    print_benchmark(rescan_name2, start_time, end_time, fasta_filename, filePtr) ;
    if (filePtr != NULL) { fprintf(filePtr, "\t") ; }

    fclose(seqFile) ;

    return 0;
    
}

int compute_from_sequence(char *sequence_input, int K, int S){

    size_t kmer_size = (size_t)K ;
    size_t smer_size = (size_t)S ;

    // convert sequence in binary format

    size_t len = strlen(sequence_input);
    uint8_t *encoded_seq = malloc(len * sizeof(uint8_t));
    if (encoded_seq == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // Loop through the sequence and convert each character
    for (size_t i = 0; i < len; i++) {
        encoded_seq[i] = base_to_bits(sequence_input[i]);
    }

    size_t num_syncmer_rescan;
    size_t num_syncmer_naive;
    size_t num_syncmer_rescan_long_array;
    size_t num_syncmer_rescan_deque;
    size_t num_syncmer_rescan_iterator;

    // CHECK CORRECTNESS
    compute_closed_syncmers_naive(encoded_seq, len, K, S, &num_syncmer_naive) ;

    compute_closed_syncmers(encoded_seq, len, K, S, &num_syncmer_rescan) ;

    compute_closed_syncmers_rescan(encoded_seq, len, K, S, &num_syncmer_rescan_long_array) ;

    compute_closed_syncmers_deque_rayan(encoded_seq, len, K, S, &num_syncmer_rescan_deque) ;

    compute_closed_syncmers_rescan_iterator(encoded_seq, len, K, S, &num_syncmer_rescan_iterator) ;


    if (num_syncmer_naive != num_syncmer_rescan){
        printf("NAIVE IS: %lu ; RESCAN IS: %lu\n", num_syncmer_naive, num_syncmer_rescan) ;
        exit(-1) ;
    }
    else if (num_syncmer_naive != num_syncmer_rescan_long_array){
        printf("NAIVE IS: %lu ; RESCAN LONG ARRAY IS: %lu\n", num_syncmer_naive, num_syncmer_rescan_long_array) ;
        exit(-1) ;
    }
    else if (num_syncmer_naive != num_syncmer_rescan_deque){
        printf("NAIVE IS: %lu ; DEQUE IS: %lu\n", num_syncmer_naive, num_syncmer_rescan_deque) ;
        exit(-1) ;
    }
    else if (num_syncmer_naive != num_syncmer_rescan_iterator){
        printf("NAIVE IS: %lu ; RESCAN ITERATOR IS: %lu\n", num_syncmer_naive, num_syncmer_rescan_iterator) ;
        exit(-1) ;
    }
    else{ return 0 ; }
}

/*--------------*/
/* --- MAIN ----*/
/*-------------*/
int main(int argc, char *argv[]) {

    if(argc <5) {
        fprintf(stderr, "Usage: %s SEQUENCE K S F\n", argv[0]);
        return 1;
    }

    char *sequence_input = argv[1];
    int K = atoi(argv[2]);
    int S = atoi(argv[3]);
    int F = atoi(argv[4]);

    // printf("FILE: %s, K: %d, S: %d, F: %d", sequence_input, K, S, F) ;

    if(S >= K) {
        fprintf(stderr, "Error: S (%d) must be less than K (%d)\n", S, K);
        return 1;
    }

    if(F == 0){
        return compute_from_sequence (sequence_input, K, S) ;
    }
    else if (F == 1){
        return compute_from_file(sequence_input, K , S) ;
    }
    else {
        fprintf(stderr, "Error: F must be 0 (sequence) or 1 (file)\n");
        return 1;
    }


    return 0 ;
}