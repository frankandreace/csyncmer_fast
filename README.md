# csyncmer_fast
Header only library for fast syncmer extraction from a binary sequence

### How to compile
Please do: 
```
gcc -march=native -O3 -o test test.c
```

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
On my machine, with 2.6 GH and hyperthreading disabled I get:

460 MB/S for HASHING

42 MB/S for NAIVE

82 MB/S for RE-SCAN WITH CIRCULAR ARRAY STRUCT

85 MB/S for RE-SCAN WITHOUT STRUCT AND WITH LARGE ARRAY

No significative change by using an array of length "window_size" and using a 256 MB array for the rescan computation.

#### Input
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

- [x] Rescan function implementation
- [x] Unit test it
- [x] Wrap it for time and throughput estimation
- [x] Test on file for speed
- [ ] Add Makefile
- [ ] Rescan iterator implementation
- [ ] Integrate Prof Sadakane AVX hashing
- [ ] Integrate Prof Sadakane AVX syncmer computation
- [ ] Run speed test on them
- [ ] Python script to plot
- [ ] Update README