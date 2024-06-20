#ifndef MATGEN_H
#define MATGEN_H

#define MAX_CHAIN_LEN 64

#define POW2(BITPOS) (1UL << (BITPOS))

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

void
gencrs( uint64_t* rowPtr,
        double* valArr, 
        uint64_t* colInd, 
        double* jump_rate_Array, 
        int* bit, 
        uint64_t dim, 
        int n);

#endif
