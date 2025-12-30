/*  File: seqhash.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: header file for seqhash package - minimizers and moshers
 * Exported functions: see below
 * HISTORY:
 * Last edited: Jul 24 14:38 2023 (rd109)
 * Created: Mon Mar  5 08:43:45 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "utils_d.h"

typedef struct {
  int seed ;			/* seed */
  int k ;			/* kmer */
  int w ;			/* window */
  U64 mask ;			/* 2*k bits */
  int shift1, shift2 ;
  U64 factor1, factor2 ;
  U64 patternRC[4] ;		/* one per base */
} SeqhashD ;

typedef struct {
  SeqhashD *sh ;
  char *s, *sEnd ;     		/* sequence currently being hashed, end marker */
  U64 h, hRC ;			/* current k-mer values */
  U64 *hash ;			/* buffer of length w holding hashes for current window */
  bool *isForward ;		/* buffer of length w holding isForward for current window */
  int base ;			/* start of buf in sequence */
  int iStart, iMin ;		/* position in buf of start of current window, next min */
  U64 min ;                     /* needed for syncmers */
  bool isDone ;
} SeqhashIteratorD ;

SeqhashD *seqhashCreateD (int k, int w, int seed) ;
static void seqhashDestroyD (SeqhashD *sh) { free (sh) ; }

// for all iterators sequence must continue to exist through the life of the iterator
// all the *next functions return any/all of kmer, pos, isF - get hash from seqhash(sh,kmer)

// simple iterator to return all kmers
SeqhashIteratorD *seqhashIteratorD (SeqhashD *sh, char *s, int len) ;
bool seqhashNextD (SeqhashIteratorD *si, U64 *kmer, int *pos, bool *isF) ;

static void seqhashIteratorDestroyD (SeqhashIteratorD *si)
{ free (si->hash) ; free (si->isForward) ; free (si) ; }

// (closed) syncmer extracts w-mers that end with a minimal kmer
// these provide a cover, and have good distribution properties
SeqhashIteratorD *syncmerIteratorD (SeqhashD *sh, char *s, int len) ;
bool syncmerNextD (SeqhashIteratorD *si, U64 *kmer, size_t *pos, bool *isF) ;

// utilities
static inline U64 kHashD (SeqhashD *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1) ; }


/******* end of file ********/
