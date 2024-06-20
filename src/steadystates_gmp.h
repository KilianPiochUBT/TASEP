
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "lib/io.h"

#if defined(__linux__)

    #define OS_TYPE 'L'
    #include <getopt.h>
    #include <time.h>
    #include <gmp.h>

#elif defined(__APPLE__)

    #define OS_TYPE 'M'
    #include <unistd.h>
    #include <getopt.h>
    #include <sys/stat.h>
    #include <sys/time.h>

#else
#   error "Unsupported Operating System"
#endif


struct option   long_opt[] =
    {
      { "nodes"             , required_argument , NULL , 0},
      { "alpha"             , required_argument , NULL , 1},
      { "beta"              , required_argument , NULL , 2},
      { "output"            , required_argument , NULL , 3},
      { "precision"         , required_argument , NULL , 4},
      { "print-precision"   , required_argument , NULL , 5},
      { "pp"                , required_argument , NULL , 5},
      { "help"              , no_argument       , NULL , 6},
      { "verbose"           , no_argument       , NULL , 7},
      {NULL                 , 0                 , NULL , 0}
   };



uint_fast8_t 
allocfill_B(mpf_t*** __B, uint32_t n);

uint_fast8_t
allocfill_Z(mpf_t ** __Z, mpf_t ** __B, uint32_t n, double alpha, double beta);

uint_fast8_t
allocfill_Rho(mpf_t ** __Rho, mpf_t ** __B, mpf_t * __Z, uint32_t n, double beta);

uint_fast64_t
fact(uint_fast64_t n);

uint_fast8_t
calculate_Z_i(mpf_t * Z_i, double alpha, double beta, int i);
