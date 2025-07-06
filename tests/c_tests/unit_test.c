#include <assert.h>
#include <stdio.h>

#include "csyncmer_fast/iterator_syng.h"

extern const unsigned char  base_to_bits_array[];

char* SMALL_SEQUENCE = "ACGTAACGTATACGTA";
U64 KMER_LENGTH = 7;
U64 SMER_LENGTH = 4;

void test_base_to_bits(){
    assert(base_to_bits_array[0] == 0);
    assert(base_to_bits_array[1] == 1);
    assert(base_to_bits_array[2] == 2);
    assert(base_to_bits_array[3] == 3);

    assert(base_to_bits_array['A'] == 0);
    assert(base_to_bits_array['C'] == 1);
    assert(base_to_bits_array['G'] == 2);
    assert(base_to_bits_array['T'] == 3);
    assert(base_to_bits_array['U'] == 3);

    assert(base_to_bits_array['a'] == 0);
    assert(base_to_bits_array['c'] == 1);
    assert(base_to_bits_array['g'] == 2);
    assert(base_to_bits_array['t'] == 3);
    assert(base_to_bits_array['u'] == 3);

}

void test_closed_syncmer_next(){
    
    Syncmer64 my_syncmer = {0,0,0};
    SyncmerIteratorS *my_syncmer_iterator = syncmer_generator_createS(SMALL_SEQUENCE, strlen(SMALL_SEQUENCE), KMER_LENGTH, SMER_LENGTH);
    bool return_flag;
    
    return_flag = syncmer_iterator_nextS(my_syncmer_iterator, &my_syncmer);
    assert(return_flag == true);
    // printf("NEW_SYNCMER: %llu\t%lu\t%lu\n",my_syncmer.hash_value, my_syncmer.kmer_position, my_syncmer.smer_position);
    assert(my_syncmer.hash_value == 5);
    assert(my_syncmer.kmer_position == 1);
    assert(my_syncmer.smer_position == 4);

    return_flag = syncmer_iterator_nextS(my_syncmer_iterator, &my_syncmer);
    assert(return_flag == true);
    assert(my_syncmer.hash_value == 5);
    assert(my_syncmer.kmer_position == 4);
    assert(my_syncmer.smer_position == 4);

    return_flag = syncmer_iterator_nextS(my_syncmer_iterator, &my_syncmer);
    assert(return_flag == true);
    assert(my_syncmer.hash_value == 52);
    assert(my_syncmer.kmer_position == 6);
    assert(my_syncmer.smer_position == 6);

    return_flag = syncmer_iterator_nextS(my_syncmer_iterator, &my_syncmer);
    assert(return_flag == true);
    assert(my_syncmer.hash_value == 52);
    assert(my_syncmer.kmer_position == 7);
    assert(my_syncmer.smer_position == 10);

    // CHECKING THAT THE RETURN FLAG HAS CHANGED TO FALSE AND THE SYNCMER HAS NOT BEEN CHANGED
    return_flag = syncmer_iterator_nextS(my_syncmer_iterator, &my_syncmer);
    assert(return_flag == false);
    assert(my_syncmer.hash_value == 52);
    assert(my_syncmer.kmer_position == 7);
    assert(my_syncmer.smer_position == 10);

}


int main(int argc, char **argv){

    // 1 - TESTING NUCLEOBASE TO BIT ENCODING 
    test_base_to_bits();

    test_closed_syncmer_next();

    printf("All tests run succesfully.");

    return 0;
}
