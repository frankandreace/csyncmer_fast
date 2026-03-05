#include "seqhash.h"

void compute_closed_syncmers_syng_original(char *sequence_input, int len, int K, int S, size_t *num_syncmers) {
  /* Durbin's SYNG implementation of closed syncmers enumeration from a sequence */
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
