# csyncmer_fast
Header only library for fast syncmer extraction from a binary sequence

### How to compile
Please do: 
```
gcc -march=native -O3 -o test test.c
```

### Basic information
This library aims at providing a header only library, in C, to compute closed syncmers in a given sequence.

#### Input
The library has a main interface that is composed of a function that accepts a binary sequence (A is 00, C is 01, G is 10 and T is 11) , a k-mer lenght value K and a s-mer length value S.

#### Output
The functions returns, as an iterator, the computed closed syncmers in the sequence, one after the other.

### Closed Syncmers
A k-mer is a closed filter iff the minimal LEFTMOST canonical s-mer (s < k) that it contains is either the first or the last present in the k-mer (by scanning it left to right).

### Re-scan implementation


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
- [ ] Updarte readme