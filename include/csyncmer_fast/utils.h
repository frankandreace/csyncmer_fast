#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h> 

// useful variables

#ifndef UNIT_DEFINED
#define UNIT_DEFINED

typedef unsigned long long U64 ;
typedef __uint128_t U128;
const static U64 U64MAX = 0xffffffffffffffff ;


const static __uint128_t U128MAX = (__uint128_t)(-1) ;

// Helper function to print __uint128_t in hexadecimal (still needed for printf)
void print_uint128_hex_debug(__uint128_t val) {
    uint64_t high = (uint64_t)(val >> 64);
    uint64_t low = (uint64_t)val;
    fprintf(stderr, "%016lX%016lX", high, low); // Print to stderr for debug
}

#endif

