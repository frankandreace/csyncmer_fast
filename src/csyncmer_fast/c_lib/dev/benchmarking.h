#include <time.h>
#include <sys/stat.h>
#include <ctype.h>
#include <stdio.h>

// Get file size in bytes
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
    
    // precomputing useful variables
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    // file size
    off_t file_size = get_file_size(filename);
    double file_size_mb = file_size / (1024.0 * 1024.0);
    double processing_speed = file_size_mb / elapsed_time;

    if (filePtr != NULL) { fprintf(filePtr, "%lf", processing_speed) ; }

    printf("%s Performance metrics.\n", name);

    // printf("Syncmers found: %d\n", syncmer_count);
    // print time stats
    print_timing_stats(elapsed_time);
    print_troughput(file_size_mb, processing_speed);
    // estimated throughput

    printf("\n");
}
