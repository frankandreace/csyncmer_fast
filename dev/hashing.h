#include "csyncmer_fast/utils.h"

/*---- STRUCTURES ----*/

typedef struct {
    int seed ;			/* seed */
    int k ;			/* kmer */
    int w ;			/* window */
    U64 mask ;			/* 2*k bits */
    int shift1, shift2 ;
    U64 factor1, factor2 ;
    U64 patternRC[4] ;		/* one per base */
  } Seqhash ;

typedef struct {
Seqhash *sh ;
char *s, *sEnd ;     		/* sequence currently being hashed, end marker */
U64 h, hRC, hash ;			/* current k-mer values */
bool isDone;
} SeqhashIterator ;


/*---- HAHSING PART ----*/

static inline U64 kHash (Seqhash *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1) ; }

static inline U64 hashRC (SeqhashIterator *si)
{ U64 hashF = kHash (si->sh, si->h) ;
  U64 hashR = kHash (si->sh, si->hRC) ;
  if (hashF < hashR) { return hashF ; }
  else { return hashR ; }
}

static inline U64 advanceHashRC (SeqhashIterator *si)
{ Seqhash *sh = si->sh ;
  if (si->s < si->sEnd)
    { si->h = ((si->h << 2) & sh->mask) | *si->s ;
      si->hRC = (si->hRC >> 2) | sh->patternRC[(int)*si->s] ;
      // printf("ADVRC. H: %llu, HRC: %llu\n", si->h, si->hRC) ;
      ++si->s ;
      return hashRC (si) ;
    }
  else
    return U64MAX ;
}

Seqhash *seqhashCreate (int k, int w, int seed)
{
  assert (sizeof (U64) == 8) ;
  // printf("malloc seqhash\n") ;
  Seqhash *sh = (Seqhash *)malloc(sizeof(Seqhash)) ;
  // printf("end seqhash is %s \n", sh == NULL ? "null" : "working") ;

  sh->k = k ; if (k < 1 || k >= 32) {printf ("seqhash k %d must be between 1 and 32\n", k) ; exit (-1) ;}
  sh->w = w ; if (w < 1) {printf ("seqhash w %d must be positive\n", w) ; exit (-1) ;} 
  sh->seed = seed ;
  sh->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  
  srandom (seed) ;
  sh->factor1 = (random() << 32) | random() | 0x01 ;
  sh->shift1 = 64 - 2*k ;
  sh->factor2 = (random() << 32) | random() | 0x01 ;
  sh->shift2 = 2*k ;
  for (i = 0 ; i < 4 ; ++i) { 
    sh->patternRC[i] = (3-i) ; 
    sh->patternRC[i] <<= 2*(k-1) ; 
  }
  // printf("end seqhash creation \n") ;
  return sh ;
}

SeqhashIterator *seqhashIterator (Seqhash *sh, char *s, int len)
{
    assert (s && len >= 0) ;
    SeqhashIterator *si =  (SeqhashIterator *)malloc(sizeof(SeqhashIterator)) ;
    si->sh = sh ;
    si->s = s ; si->sEnd = s + len ;
    si->hash = 0;
    si->h = 0;
    si->hRC = 0;
    si->isDone = false;
    if (len < sh->k) { si->isDone = true ;} // edge case
    else
      { 
        int i ;			/* preinitialise the hashes for the first kmer */
        // printf("H is %llu, HRC is %llu\n", si->h, si->hRC) ;
        for (i = 0 ; i < sh->k ; ++i, ++si->s)
        { 
          si->h = (si->h << 2) | *si->s ;
          si->hRC = (si->hRC >> 2) | sh->patternRC[(int)*si->s] ;
          // printf("H is %llu, HRC is %llu, i is %d\n", si->h, si->hRC,i) ;
        }
        // U64 x = hashRC (si) ;
        si->hash = hashRC (si) ;
      }
    return si ;
}

bool seqhashNext (SeqhashIterator *si, U64 *kmer)
{
    if (si->isDone){
      return false ; /* we are done */
    }
    if (kmer) *kmer = si->hash ;

    if (si->s >= si->sEnd){
      si->isDone = true ;
    }
    else si->hash = advanceHashRC (si) ;
    return true ;
}