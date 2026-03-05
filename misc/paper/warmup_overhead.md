# Warmup overhead in chunk-parallel SIMD twostack

## Problem

The single-sequence SIMD twostack splits one sequence into 8 chunks across 8
AVX2 lanes. Each lane must independently warm up its rolling hash (S-1 iters)
and twostack window (window_size-1 iters) before producing valid syncmer
decisions. Adjacent chunks overlap by K-1 s-mer positions.

## Quantification

Per lane, the warmup is K-1 iterations out of chunk_len + K-1 total iterations:

    warmup fraction = (K-1) / (num_kmers/8 + K-1)

For HiFi reads (avg ~7376 bp) with K=1052, S=31:

    num_kmers = 7376 - 1052 + 1 = 6325
    chunk_len = 6325 / 8 ~ 791
    warmup = K - 1 = 1051
    warmup fraction = 1051 / (791 + 1051) = 57%

With 8 lanes and 7 chunk boundaries, total redundant work is 7*(K-1) = 7357
iterations -- more than the 6325 actual k-mers.

## Solution: multi-read SIMD

Instead of splitting one read into 8 chunks, put 8 different reads into the 8
lanes. Each lane processes a complete read, so warmup occurs once per read (the
minimum possible). No redundant overlap.

The warmup fraction per read is then just (K-1) / (num_kmers + K-1), which for
a 7376 bp read is 1051/7376 = 14%, vs 57% with chunking.
