#ifndef CSYNCMER_FASTQ_H
#define CSYNCMER_FASTQ_H

// Multi-read SIMD extension for csyncmer_fast.
// Processes 8 separate reads in 8 AVX2 lanes (optimal for short reads / FASTQ).
// Includes the upstream single-read API via csyncmer_fast.h.

#include "../../csyncmer_fast.h"
#include "csyncmer_fastq_pack.h"
#include "csyncmer_fastq_twostack.h"
#include "csyncmer_fastq_hash.h"

#endif  // CSYNCMER_FASTQ_H
