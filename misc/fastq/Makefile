CC = cc
CFLAGS = -O3 -mavx2 -march=native -g
LDFLAGS = -lz
SYNG_DIR = ../syng

SYNG_OBJS = $(SYNG_DIR)/seqhash.o $(SYNG_DIR)/utils_d.o

$(SYNG_DIR)/%.o: $(SYNG_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

all: bench_syncmer_fastq

bench_syncmer_fastq: bench_syncmer_fastq.c csyncmer_fastq.h ../../csyncmer_fast.h $(SYNG_OBJS)
	$(CC) $(CFLAGS) -o $@ $< $(SYNG_OBJS) $(LDFLAGS)

clean:
	rm -f bench_syncmer_fastq $(SYNG_DIR)/seqhash.o $(SYNG_DIR)/utils_d.o

.PHONY: all clean
