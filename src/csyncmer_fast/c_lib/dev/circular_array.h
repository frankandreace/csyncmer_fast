#include "csyncmer_fast/utils.h"

/*---- circular array structure to handle the hash value of the s-mers in a sliding window ----*/
typedef struct {
    size_t current_position ; // keeps track of the current position in the circula array (for window scanning purposes)
    size_t window_size ; // the number of s-mers in the k-mers (K-S+1)
    U64 minimum ;
    size_t minimum_position ;
    U64 hashVector[] ; // contains the hashes of the s-mers for the length of a window w
} CircularArray ;

CircularArray *circularArrayCreate(size_t window_size) {
    CircularArray *ca = (CircularArray *)malloc(sizeof(CircularArray) + window_size * sizeof(U64)) ;

    if (ca == NULL){
        fprintf(stderr, "CANNOT INITIALIZE CIRCULAR ARRAY.") ;
        exit (-1) ;
    }

    ca->current_position = 0 ;
    ca->window_size = window_size ;
    ca->minimum = U64MAX;
    ca->minimum_position = window_size + 1;
    return ca ;
}

static void circularArrayDestroy (CircularArray *ca) {if (ca != NULL) free(ca);}

// print array status
void print_status(CircularArray *ca){
    printf("[") ;
    for (size_t i = 0 ; i < ca->window_size ; i++ ) {
        printf("%llu,", ca->hashVector[i]) ;
    }
    printf("]\n") ;
    printf("MIN: %llu ; MIN_P: %lu ; CURR_P : %lu ; WS: %lu\n", ca->minimum, ca->minimum_position, ca->current_position, ca->window_size) ;
}


static inline void update_minimum_branchless(U64 candidate, U64 *current_minimum, size_t candidate_position, size_t *current_minimum_position)
{
    // U64 smaller = (candidate < *current_minimum);
    U64 mask = -(candidate < *current_minimum);
    *current_minimum = (*current_minimum & ~mask) | (candidate & mask);
    *current_minimum_position  = (*current_minimum_position & ~mask) | (candidate_position & mask);
}


/*---- perform a re-scan of the entire array when the current minimum is out of context and returnt the min and position ----*/
void circularScanBranchless(CircularArray *ca){

    size_t scan_position ;
    ca->minimum = U64MAX ;
 
    for (size_t i = 1 ; i < ca->window_size ; i++ ) {
      scan_position = (i + ca->current_position) % ca->window_size ;
      update_minimum_branchless(ca->hashVector[scan_position], &ca->minimum, scan_position, &ca->minimum_position) ;
    }
}

/*---- insert a new element in the circular array----*/
void circularInsertBranchless(CircularArray *ca, U64 value) {

    if (ca->minimum_position == ca->current_position) { circularScanBranchless(ca) ; }

    ca->hashVector[ca->current_position++] = value;

    update_minimum_branchless(value, &ca->minimum, ca->current_position - 1, &ca->minimum_position) ;
    if (ca->current_position == ca->window_size) {ca->current_position = 0;} // go back to zero if at the end of the array
}

/*---- perform a re-scan of the entire array when the current minimum is out of context and returnt the min and position ----*/
void circularScan(CircularArray *ca){

    size_t scan_position ;
    // U64 current_minimum = U64MAX ;
    size_t current_minimum_position = (1 + ca->current_position) % ca->window_size; 
 
    for (size_t i = 1 ; i < ca->window_size ; i++ ) {
      scan_position = (i + ca->current_position) % ca->window_size ;
      if ( ca->hashVector[scan_position] < ca->hashVector[current_minimum_position] ){
        // current_minimum = ca->hashVector[scan_position] ;
        current_minimum_position = scan_position;
      }
    }

    ca->minimum = ca->hashVector[current_minimum_position] ;
    ca->minimum_position = current_minimum_position ;
}

/*---- insert a new element in the circular array----*/
void circularInsert(CircularArray *ca, U64 value) {
    // if minimum out of scope, recompute
    if (ca->minimum_position == ca->current_position) { circularScan(ca) ; }

    ca->hashVector[ca->current_position++] = value;
    if (value < ca->minimum){
        ca->minimum = value ;
        ca->minimum_position = ca->current_position - 1 ;
    }
    if (ca->current_position == ca->window_size) {ca->current_position = 0;} // go back to zero if at the end of the array
}

bool is_syncmer(CircularArray *ca, size_t *position){
    //verify if it is at the begin or end of the window
    if (ca->minimum_position == ca->current_position) {
        *position = 0;
        return  true;
    }
    else if (ca->minimum_position == (ca->current_position + ca->window_size - 1) % ca->window_size)
    {
        *position = ca->window_size - 1;
        return  true;
    }
    return false;
}

U64 get_smer(CircularArray *ca){
    return ca->hashVector[ca->minimum_position];
}
