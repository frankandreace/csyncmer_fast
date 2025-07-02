# csyncmer_fast
Header only library for fast syncmer extraction from a binary sequence

### CI/CD
[![Build Docker Image for build and test](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/build-docker.yml)

[![C speed benchmark build](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/c_speed_bench_build.yml)

[![python package compilation](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml/badge.svg)](https://github.com/frankandreace/csyncmer_fast/actions/workflows/python-build.yml)  

### How to compile the speed benchmark binary
Please do: 
```
make
```

### How to run speed bench
Please do
```
./test.sh
```
In the benchmarks folder there will be the tsv table with the speed for the 10 rounds of tests and 2 plots of the benchmarked speed of the different implementations.

### How to download 



### Benchmark against current other libraries that compute closed syncmers
Please see [the github directory of the comprehensive tests](https://github.com/frankandreace/csyncmer_fast_benchmark).

### HOW TO USE
If you want to give a sequence (GCAAGTGACAATTCCTGAGAATAAT) in command line with k=5 and s=2, please use:

```
./test GCAAGTGACAATTCCTGAGAATAAT 5 2 0
```
If want to pass a fasta file (data/chr19_bit.fa), k=31 and s=11, please use:

```
./test data/chr19_bit.fa 31 11 1
```

### Basic information
This library aims at providing a header only library, in C, to compute closed syncmers in a given sequence.
On my machine, using an intel i9 CPU set to 2.6 GHz and hyperthreading disabled, on one core, I get:

- ~460 MB/S for HASHING
- ~90 MB/S for RE-SCAN USING A LARGE ARRAY WHERE TO STORE HASHES
- ~85 MB/S for RE-SCAN WITH A CIRCULAR ARRAY STRUCTURE
- ~80 MB/S for RE-SCAN ITERATOR WITH CIRCULAR ARRAY STRUCT
- ~47 MB/S for DEQUE IMPLEMENTATION
- ~42 MB/S for NAIVE IMPLEMENTATION

No significative change by using an array of length "window_size" and using a 256 MB array for the rescan computation.

#### Input
At the moment I am converting ASCII of the Fasta char into the binary encoding.

The library has a main interface that is composed of a function that accepts a binary sequence (A is 00, C is 01, G is 10 and T is 11) , a k-mer lenght value K and a s-mer length value S.

#### Output
The functions returns, as an iterator, the computed closed syncmers in the sequence, one after the other.

### Closed Syncmers
A k-mer is a closed filter iff the minimal LEFTMOST canonical s-mer (s < k) that it contains is either the first or the last present in the k-mer (by scanning it left to right).

### Hash function
I am using Richard Durbin's hash function used in SYNG. Please see below for the link.

### Some of existing code

Strobealign code: uses 64-bit smers 
https://github.com/ksahlin/strobealign/blob/71866c31b578e5166c83aaf1fde79d238246490d/src/randstrobes.cpp#L57

Minimap2 code: uses 64-bit smers
https://github.com/lh3/minimap2/blob/c2f07ff2ac8bdc5c6768e63191e614ea9012bd5d/sketch.c#L145-L192

Curiouscoding blog:
https://curiouscoding.nl/posts/fast-minimizers/

Sliding window minimum algorithm explanation:
https://github.com/keegancsmith/Sliding-Window-Minimum?tab=readme-ov-file#sliding-window-minimum-algorithm

SYNG: 
https://github.com/richarddurbin/syng/


### TO DO LIST:

- [ ] Add unit tests for rescan
- [x] Write nextflow pipeline for benchmark
- [x] Add simd_minimizer benchmark 
- [x] Restructure folder as benchmarking folder
- [ ] Integrate Prof Sadakane AVX hashing
- [ ] Integrate Prof Sadakane AVX syncmer computation
- [x] Rescan function implementation
- [x] Unit test it
- [x] Wrap it for time and throughput estimation
- [x] Test on file for speed
- [x] Add Makefile
- [x] Rescan iterator implementation
- [x] Run speed test on them
- [x] Python script to plot
