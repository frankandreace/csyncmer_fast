#ifndef HASHING_H
#define HASHING_H

#include "utils.h"

// BASE TO BITS ARRAY
// USING CHAR CONVERTION TO INT
// DEFAULT IS 0
// 65: A ; 97: a
// 67: C ; 99: c
// 71: G ; 103: g
// 84: T ; 116: t
// 85: U ; 117: u
// first 4 elements as it is receiving the sequence already binarized
static const unsigned char base_to_bits_array[256] = {
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0};

static inline char base_to_bits(char base) {
    switch(base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0; // Treat Ns and  unknown as 'A'
    }
}

/*---- STRUCTURES ----*/

typedef struct
{
  int seed; /* seed */
  int k;    /* kmer */
  int w;    /* window */
  U64 mask; /* 2*k bits */
  int shift1, shift2;
  U64 factor1, factor2;
  U64 patternRC[4]; /* one per base */
} Seqhash;

typedef struct
{
  Seqhash *sh;
  char *s, *sEnd;   /* sequence currently being hashed, end marker */
  U64 h, hRC, hash; /* current k-mer values */
  bool isDone;
} SeqhashIterator;

/*---- HAHSING PART ----*/

static inline U64 kHash(Seqhash *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1); }

static inline U64 hashRC(SeqhashIterator *si, bool *isForward)
{
  U64 hashF = kHash(si->sh, si->h);
  U64 hashR = kHash(si->sh, si->hRC);
  if (hashF < hashR)
  {
    *isForward = true;
    return hashF;
  }
  else
  {
    *isForward = false;
    return hashR;
  }
}

static inline U64 advanceHashRC(SeqhashIterator *si, bool *isForward)
{
  Seqhash *sh = si->sh;
  if (si->s < si->sEnd)
  {
    si->h = ((si->h << 2) & sh->mask) | base_to_bits_array[*si->s];
    si->hRC = (si->hRC >> 2) | sh->patternRC[base_to_bits_array[*si->s]];
    ++si->s;
    return hashRC(si, isForward);
  }
  else
    return U64MAX;
}

Seqhash *seqhashCreate(int k, int w, int seed)
{
  assert(sizeof(U64) == 8);
  Seqhash *sh = (Seqhash *)malloc(sizeof(Seqhash));
  sh->k = k;
  if (k < 1 || k >= 32)
  {
    printf("seqhash k %d must be between 1 and 32\n", k);
    exit(-1);
  }
  sh->w = w;
  if (w < 1)
  {
    printf("seqhash w %d must be positive\n", w);
    exit(-1);
  }
  sh->seed = seed;
  sh->mask = ((U64)1 << (2 * k)) - 1;
  int i;

  srandom(seed);
  sh->factor1 = (random() << 32) | random() | 0x01;
  sh->shift1 = 64 - 2 * k;
  sh->factor2 = (random() << 32) | random() | 0x01;
  sh->shift2 = 2 * k;
  for (i = 0; i < 4; ++i)
  {
    sh->patternRC[i] = (3 - i);
    sh->patternRC[i] <<= 2 * (k - 1);
  }
  return sh;
}

SeqhashIterator *seqhashIterator(Seqhash *sh, char *s, int len)
{
  assert(s && len >= 0);
  SeqhashIterator *si = (SeqhashIterator *)malloc(sizeof(SeqhashIterator));
  si->sh = sh;
  si->s = s;
  si->sEnd = s + len;
  si->hash = 0;
  si->h = 0;
  si->hRC = 0;
  si->isDone = false;
  bool isForward;
  if (len < sh->k)
  {
    si->isDone = true;
  } // edge case
  else
  {
    int i; /* preinitialise the hashes for the first kmer */
    for (i = 0; i < sh->k; ++i, ++si->s)
    {
      base_to_bits_array[*si->s];
      si->h = (si->h << 2) | base_to_bits_array[*si->s];
      si->hRC = (si->hRC >> 2) | sh->patternRC[base_to_bits_array[*si->s]];
    }
    si->hash = hashRC(si, &isForward);
  }
  return si;
}

static void SeqhashDestroy(Seqhash *sh) { free(sh); }

static void SeqhashIteratorDestroy(SeqhashIterator *si)
{
  // prevent double free
  if (si == NULL)
    return;
  if (si->sh != NULL)
  {
    // Free the seqhash
    SeqhashDestroy(si->sh);
    si->sh = NULL; // prevent double free
  }
  // free the struct
  free(si);
}

bool seqhashNext(SeqhashIterator *si, U64 *kmer, bool *isForward)
{
  if (si->isDone)
  {
    return false; /* we are done */
  }
  if (kmer)
    *kmer = si->hash;

  if (si->s >= si->sEnd)
  {
    si->isDone = true;
  }
  else
    si->hash = advanceHashRC(si, isForward);
  return true;
}

#endif