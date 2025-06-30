#!/usr/bin/env python3

from csyncmer_fast import SyncmerIterator
import argparse
from time import time

def load_sequence(fasta_file: str):
    sequence: list = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence.append(line.strip())

        return ''.join(sequence)

def benchmark(filename: str, kmer_length:int = 31, smer_length: int = 11):
    sequence = load_sequence(filename)
    sequence_size_in_MB = len(sequence) / (1024 ** 2)

    ## RETURNING THEM IN A LIST USING PYTHON ITERATOR
    print(f"Running on seq {sequence[:10]}...")
    start_time = time()
    syncmers: list = list(SyncmerIterator(sequence, kmer_length,smer_length))
    elapsed_time = time() - start_time

    print("----- RETURN ALL AS LIST OF TUPLE USING PYTHON ITERATOR-----")
    print(f"{syncmers[0]}")
    print(f"ELAPSED TIME:{elapsed_time}")
    print(f"SEQUENCE SIZE:{sequence_size_in_MB}")
    speed_MB_sec = sequence_size_in_MB / elapsed_time
    num_syncmer_computed: int = len(syncmers)
    print(f"SPEED MB/S {speed_MB_sec}")
    print(f'{speed_MB_sec:.3f}\t{num_syncmer_computed}\t{len(sequence)}')
    
    #ITERATING OVER ELEMENTS SYNCMER BY SYNCMER
    print(f"Running on seq {sequence[:10]}...")
    start_time = time()
    tot_syncmer = 0
    for syncmer in SyncmerIterator(sequence, kmer_length,smer_length):
        tot_syncmer+=1
    elapsed_time = time() - start_time

    print("----- LOOP OVER ITERATOR 1 ELEMENT AT A TIME -----")
    print(f"ELAPSED TIME:{elapsed_time}")
    print(f"SEQUENCE SIZE:{sequence_size_in_MB}")
    speed_MB_sec = sequence_size_in_MB / elapsed_time
    print(f"SPEED MB/S {speed_MB_sec}")
    print(f'{speed_MB_sec:.3f}\t{tot_syncmer}\t{len(sequence)}')

    #DIRECTLY RETURN THEM WITHOUT USING AN ITERATOR
    print(f"Running on seq {sequence[:10]}...")
    start_time = time()
    tot_syncmer = 0
    syncmer_iterator = SyncmerIterator(sequence,kmer_length,smer_length)
    hash_values, kmer_positions, smer_positions = syncmer_iterator.get_all_syncmers()
    elapsed_time = time() - start_time

    print("----- DIRECT RETURN ALL AS 3 NUMPY ARRAYS -----")
    print(f"{hash_values[0]},{kmer_positions[0]},{smer_positions[0]}")
    print(f"ELAPSED TIME:{elapsed_time}")
    print(f"SEQUENCE SIZE:{sequence_size_in_MB}")
    speed_MB_sec = sequence_size_in_MB / elapsed_time
    print(f"SPEED MB/S {speed_MB_sec}")
    print(f'{speed_MB_sec:.3f}\t{len(hash_values)}\t{len(sequence)}')


# my_syncmers = list(SyncmerIterator(sequence, kmer_length, smer_length)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fasta", type=str, help="input fasta file", required=True)
    # parser.add_argument("-o", "--output", type=str, help="output benchmark tsv", required=True)
    parser.add_argument("-k", "--kmer_size", type=int, help="k-mer length", default=31)
    parser.add_argument("-s", "--smer_size", type=int, help="s-mer length", default=11)
    args = parser.parse_args()

    benchmark(args.input_fasta, args.kmer_size, args.smer_size) #, args.output
