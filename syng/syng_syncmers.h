#include "seqhash.h"

void compute_closed_syncmers_syng_original(char *sequence_input, int len, int K, int S, size_t *num_syncmers) {
    //, MinimizerResult *results, int *num_results
  /* Durbin's SYNG implementation of closed syncmers enrumeration from a sequence*/
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }
    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    size_t syncmer_count = 0;

    // Durbin's function works with the parameters as follows:
    // its K is normal S
    // its W is the window size, i.e. K-S+1
    size_t window_size = (U64)K - (U64)S + 1;

    SeqhashD *sh = seqhashCreateD(S, window_size, seed);

    // initializing the syncmer iterator
    SeqhashIteratorD *si = syncmerIteratorD(sh, sequence_input, len); 
    size_t count = 0 ;

    if (si->iMin != U64MAXD) {
        U64 kmer ;
        size_t s_pos = U64MAXD;
        size_t last_k_pos;
        bool isF ;

        while (syncmerNextD(si, &kmer, &s_pos, &isF))
        {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            last_k_pos = s_pos ;
            count++;
        }
        // printf("I SEE k_pos == %lu, last_k_pos == %lu\n", k_pos, last_k_pos) ;
        if(s_pos != U64MAXD && last_k_pos != s_pos) {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            syncmer_count++;
        }
    }
    printf("COUNT IS %lu\n", count) ;
    // releasing memory
    seqhashIteratorDestroyD(si);
    seqhashDestroyD(sh);

    *num_syncmers = syncmer_count;
}