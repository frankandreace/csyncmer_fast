#include "seqhash.h"

/* seqhash takes bases as 0-3 codes (A,C,G,T), which syng's seqIO produces while
   reading (dna2index4Conv); non-ACGT -> 0 as there.  Encode outside timed regions. */
static inline void syng_encode_bases(const char *s, size_t len, char *out) {
    static unsigned char to_index[256];
    if (!to_index['C']) {
        to_index['C'] = to_index['c'] = 1;
        to_index['G'] = to_index['g'] = 2;
        to_index['T'] = to_index['t'] = 3;
    }
    for (size_t i = 0; i < len; i++) out[i] = (char)to_index[(unsigned char)s[i]];
}

void compute_closed_syncmers_syng_original(char *sequence_input, int len, int K, int S, size_t *num_syncmers) {
  /* Durbin's SYNG implementation of closed syncmers enumeration from a sequence.
     sequence_input must be 0-3 encoded (see syng_encode_bases), not ASCII. */
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        *num_syncmers = 0;
        return;
    }
    // Durbin's function works with the parameters as follows:
    // its K is normal S
    // its W is the window size, i.e. K-S+1
    U64 seed  = 7;
    size_t window_size = (size_t)K - (size_t)S + 1;

    SeqhashD *sh = seqhashCreateD(S, window_size, seed);
    SeqhashIteratorD *si = syncmerIteratorD(sh, sequence_input, len);
    size_t count = 0;

    U64 kmer;
    int pos;
    bool isF;
    while (syncmerNextD(si, &kmer, &pos, &isF))
        count++;

    printf("COUNT IS %zu\n", count);
    seqhashIteratorDestroyD(si);
    seqhashDestroyD(sh);

    *num_syncmers = count;
}
