CC = cc
CFLAGS = -O3 -mavx2 -march=native -g
LDFLAGS = -lz

bench_syncmer_fastq: bench_syncmer_fastq.c csyncmer_fastq.h ../../csyncmer_fast.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f bench_syncmer_fastq

.PHONY: clean
