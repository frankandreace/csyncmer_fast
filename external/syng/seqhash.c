/*  File: seqhash.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: seqhash package - uses random bit patterns for bases and rotate/XOR
 	compile with -DTEST to test, and -DDEBUG to debug, -DNDEBUG turns off asserts
	see test main() at end for standard usage pattern
 * Exported functions: see seqhash.h
 * HISTORY:
 * Last edited: Jan 14 09:28 2025 (rd109)
 * Created: Sat Feb 24 19:20:18 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "seqhash.h"

SeqhashD *seqhashCreateD (int k, int w, int seed)
{
  assert (sizeof (U64) == 8) ;
  SeqhashD *sh = new0 (1, SeqhashD) ;
  sh->k = k ; if (k < 1 || k >= 32) die ("seqhash k %d must be between 1 and 32\n", k) ;
  sh->w = w ; if (w < 1) die ("seqhash w %d must be positive\n", w) ;
  sh->seed = seed ;
  sh->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  
  srandom (seed) ;
  sh->factor1 = (random() << 32) | random() | 0x01 ;
  sh->shift1 = 64 - 2*k ;
  sh->factor2 = (random() << 32) | random() | 0x01 ;
  sh->shift2 = 2*k ;
  for (i = 0 ; i < 4 ; ++i) { sh->patternRC[i] = (3-i) ; sh->patternRC[i] <<= 2*(k-1) ; }
  return sh ;
}

#include <stdio.h>


/************** basic hash functions *************/

static inline U64 hashRCD (SeqhashIteratorD *si, bool *isForward)
{ U64 hashF = kHashD (si->sh, si->h) ;
  U64 hashR = kHashD (si->sh, si->hRC) ;
#ifdef DEBUG
  printf ("hashRC: h %lx hRC %lx hashF %lx hashR %lx\n", si->h, si->hRC, hashF, hashR) ;
#endif
  if (hashF < hashR) { *isForward = true ; return hashF ; }
  else { *isForward = false ; return hashR ; }
}

static inline U64 advanceHashRCD (SeqhashIteratorD *si, bool *isForward)
{ SeqhashD *sh = si->sh ;
  if (si->s < si->sEnd)
    { si->h = ((si->h << 2) & sh->mask) | *(si->s) ;
      si->hRC = (si->hRC >> 2) | sh->patternRC[(int)*(si->s)] ;
      ++si->s ;
      return hashRCD (si, isForward) ;
    }
  else
    return U64MAXD;
}

/************ iterators to run across a sequence, returning (a subset of) hashes ***********/

/*************** this basic one returns all the hashes *********************/
/*************** and its creator is used as a base by the others ***********/

SeqhashIteratorD *seqhashIteratorD (SeqhashD *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashIteratorD *si = new0 (1, SeqhashIteratorD) ;
  si->sh = sh ;
  si->s = s ; si->sEnd = s + len ;
  si->hash = new0 (sh->w, U64) ;
  si->isForward = new0 (sh->w, bool) ;
  if (len < sh->k)
    si->isDone = true ; // edge case
  else
    { int i ;			/* preinitialise the hashes for the first kmer */
      for (i = 0 ; i < sh->k ; ++i, ++si->s)
	{ si->h = (si->h << 2) | *si->s ;
	  si->hRC = (si->hRC >> 2) | sh->patternRC[(int)*(si->s)] ;
	}
      *si->hash = hashRCD (si, si->isForward) ;
    }
  return si ;
}

bool seqhashNextD (SeqhashIteratorD *si, U64 *kmer, int *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

  if (kmer) *kmer = si->h ;
  if (pos) *pos = si->iStart ;
  if (isF) *isF = *si->isForward ;

  if (si->s >= si->sEnd)
    si->isDone = true ;
  else
    { *si->hash = advanceHashRCD (si, si->isForward) ;
      ++si->iStart ;
    }
  
  return true ;
}



/************ same for closed syncmer ***********/

SeqhashIteratorD *syncmerIteratorD (SeqhashD *sh, char *s, int len)
{
  SeqhashIteratorD *si = seqhashIteratorD (sh, s, len) ;
  if (len < sh->w + sh->k) si->isDone = true ; // because we are looking for w-mers not k-mers here
  if (si->isDone) return si ;
    
  /* store first w hashes in hash and set ->min */
  si->min = si->hash[0] ;
  int i ;
  for (i = 1 ; i < sh->w ; ++i)
    { si->hash[i] = advanceHashRCD (si, &si->isForward[i]) ;
      if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }

  // si->iStart = 0 ; // from initialisation
  if (si->hash[0] == si->min || si->hash[sh->w-1] == si->min) return si ; // we are done
  while (true)
    { U64 x = advanceHashRCD (si, &si->isForward[si->iStart]) ;
      if (si->s >= si->sEnd) { si->isDone = true ; return si ; }
      si->hash[si->iStart++] = x ;
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ; return si ; }
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	return si ;
    }
  die ("syncmer initialisation failure") ;
}

bool syncmerNextD (SeqhashIteratorD *si, U64 *kmer, size_t *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, min %" PRIx64 "\n", si->base, si->iStart, si->min) ;
  int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hash[j]) ;
  printf ("\n") ;
#endif

  if (kmer) *kmer = si->hash[si->iStart] ;
  if (pos) *pos = si->base + si->iStart ;
  if (isF) *isF = si->isForward[si->iStart] ;

  if (si->hash[si->iStart] == si->min) // need to find new min - could use a heap, but not so bad to search here
    { int i ;
      si->min = si->hash[si->iStart] = U64MAXD;
      for (i = 0 ; i < si->sh->w ; ++i) if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }
  
  while (true) // move forwards to the next minimum
    { U64 x = advanceHashRCD (si, &si->isForward[si->iStart]) ;
      if (si->s >= si->sEnd) { si->isDone = true ; return true ; }
      si->hash[si->iStart++] = x ;
      if (si->iStart == si->sh->w) { si->base += si->sh->w ; si->iStart = 0 ; }
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ;
	  //	  printf (" syncmerNext   at_end %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
	  return true ;
	}
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	{
	  // printf (" syncmerNext at_start %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
	  return true ;
	}
    }
}


void seqhashForCompilerHappiness (SeqhashD *sh, SeqhashIteratorD *sit)
{ seqhashDestroyD (sh) ; seqhashIteratorDestroyD (sit) ; }

/**************** end of file ****************/
