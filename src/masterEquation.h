#if defined(__linux__)

    #define OS_TYPE 'L'
    #include <getopt.h>
    #include <sys/time.h>
    #include <sys/stat.h>

#elif defined(__APPLE__)

    #define OS_TYPE 'M'
    #include <unistd.h>
    #include <getopt.h>
    #include <sys/stat.h>
    #include <sys/time.h>

#else
#   error "Unsupported Operating System"
#endif

#include <string.h> 
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <assert.h> 
#include <math.h>



#include "lib/matvec.h"
#include "lib/matgen.h"
#include "lib/solv.h"
#include "lib/graphlib.h"
#include "lib/io.h"

#define POW2(BITPOS) (1UL << (BITPOS))



double * 
enterProbs(int num);

static void 
usage(const char *progname);

struct option   long_opt[] =
    {
      {"time-var"     , no_argument       , NULL , 0},
      {"function-type", required_argument , NULL , 1},
      {"offset"       , required_argument , NULL , 2},
      {"amplitude"    , required_argument , NULL , 3},
      {"period"       , required_argument , NULL , 4},
      {"constant"     , required_argument , NULL , 5},
      {"save-select"  , required_argument , NULL , 6},
      {"thread-num"   , required_argument , NULL , 7},
      {"alpha"        , required_argument , NULL , 8},
      {"beta"         , required_argument , NULL , 9},
      {"gamma"      , required_argument , NULL , 10},
      {NULL           , 0                 , NULL , 0}
   };

struct timeval start, stop;
struct timeval perf_start, perf_stop;
