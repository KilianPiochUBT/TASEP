#ifndef SOLV_H
#define SOLV_H

#define POW2(BITPOS) (1UL << (BITPOS))
#define MAX_CHAIN_LEN 64

#include <stdint.h>

#ifdef __has_include
#  if __has_include(<omp.h>)
#    include <omp.h>
#    define HAVE_OMP 1
#  elif __has_include(<experimental/omp.h>)
#    include <experimental/omp.h>
#    define HAVE_OMP 1
#  else
#	 undef HAVE_OMP
#  endif
#endif

#if defined(__linux__)

    #define USE_OMP 1

#elif defined(__APPLE__)

    #define USE_OMP 0

#else
#   error "Unsupported Operating System"
#endif


void 
rk4_step(   int step, 
            uint64_t dim, 
            double *k, 
            double *x, 
            double *x1,
            uint64_t * rowPtr, 
            uint64_t* colInd, 
            double *valArr, 
            double h);
 

uint_fast64_t
gen_meq_row(int n, 
            uint_fast64_t k, 
            double * jump_rate_Array, 
            double *valArr, 
            uint_fast64_t *colInd);

void
rk4_step_reduced(   int step, 
                    uint64_t valCount, 
                    uint_fast64_t currIndex, 
                    double *k, 
                    double *x, 
                    double *x1, 
                    uint64_t* colInd, 
                    double *valArr, 
                    double h);

#endif
