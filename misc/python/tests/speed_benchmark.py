#!/usr/bin/env python3

from csyncmer_fast import SyncmerIterator, CanonicalSyncmerIterator, count_syncmers, count_syncmers_canonical
import argparse
from time import time


def load_sequence(fasta_file: str):
    sequence: list = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence.append(line.strip())
    return ''.join(sequence)


def benchmark(filename: str, kmer_length: int = 31, smer_length: int = 11):
    sequence = load_sequence(filename)
    sequence_size_in_MB = len(sequence) / (1024 ** 2)

    print(f"Sequence length: {len(sequence)} ({sequence_size_in_MB:.2f} MB)")
    print(f"k={kmer_length}, s={smer_length}")
    print()

    # Forward iterator: collect as list
    start_time = time()
    positions = list(SyncmerIterator(sequence, kmer_length, smer_length))
    elapsed_time = time() - start_time
    speed = sequence_size_in_MB / elapsed_time
    print(f"--- Forward iterator (list) ---")
    print(f"  Syncmers: {len(positions)}, Time: {elapsed_time:.3f}s, Speed: {speed:.1f} MB/s")

    # Forward iterator: batch numpy
    start_time = time()
    arr = SyncmerIterator(sequence, kmer_length, smer_length).get_all_positions()
    elapsed_time = time() - start_time
    speed = sequence_size_in_MB / elapsed_time
    print(f"--- Forward batch (numpy) ---")
    print(f"  Syncmers: {len(arr)}, Time: {elapsed_time:.3f}s, Speed: {speed:.1f} MB/s")

    # SIMD count (forward)
    start_time = time()
    count = count_syncmers(sequence, kmer_length, smer_length)
    elapsed_time = time() - start_time
    speed = sequence_size_in_MB / elapsed_time
    print(f"--- SIMD count (forward) ---")
    print(f"  Syncmers: {count}, Time: {elapsed_time:.3f}s, Speed: {speed:.1f} MB/s")

    # Canonical iterator: collect as list
    start_time = time()
    results = list(CanonicalSyncmerIterator(sequence, kmer_length, smer_length))
    elapsed_time = time() - start_time
    speed = sequence_size_in_MB / elapsed_time
    print(f"--- Canonical iterator (list) ---")
    print(f"  Syncmers: {len(results)}, Time: {elapsed_time:.3f}s, Speed: {speed:.1f} MB/s")

    # Canonical iterator: batch numpy
    start_time = time()
    pos_arr, strand_arr = CanonicalSyncmerIterator(sequence, kmer_length, smer_length).get_all_positions()
    elapsed_time = time() - start_time
    speed = sequence_size_in_MB / elapsed_time
    print(f"--- Canonical batch (numpy) ---")
    print(f"  Syncmers: {len(pos_arr)}, Time: {elapsed_time:.3f}s, Speed: {speed:.1f} MB/s")

    # SIMD count (canonical)
    start_time = time()
    count_c = count_syncmers_canonical(sequence, kmer_length, smer_length)
    elapsed_time = time() - start_time
    speed = sequence_size_in_MB / elapsed_time
    print(f"--- SIMD count (canonical) ---")
    print(f"  Syncmers: {count_c}, Time: {elapsed_time:.3f}s, Speed: {speed:.1f} MB/s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fasta", type=str, help="input fasta file", required=True)
    parser.add_argument("-k", "--kmer_size", type=int, help="k-mer length", default=31)
    parser.add_argument("-s", "--smer_size", type=int, help="s-mer length", default=11)
    args = parser.parse_args()

    benchmark(args.input_fasta, args.kmer_size, args.smer_size)
