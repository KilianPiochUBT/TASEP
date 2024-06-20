#ifndef MF_FULL_H
#define MF_FULL_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>

#include "lib/io.h"
#include "lib/solv.h"
#include "lib/matvec.h"

#define POW2(BITPOS) (1UL << (BITPOS))

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


struct option   
long_opt[] =
   {
    { "alpha"           ,   required_argument , NULL , 0},
    { "beta"            ,   required_argument , NULL , 1},
    { "gamma"           ,   required_argument , NULL , 2},
    { "num-sites"       ,   required_argument , NULL , 3},
    { "input"           ,   required_argument , NULL , 4},
    { "output"          ,   required_argument , NULL , 5},
    { "trajectories"    ,   required_argument , NULL , 6},
    { "method"          ,   required_argument , NULL , 7},
    { "parallel-threads",   required_argument , NULL , 8},
    { "start"           ,   required_argument , NULL , 9},
    { "stop"            ,   required_argument , NULL , 10},
    { "steps"           ,   required_argument , NULL , 11},
    { "verbose"         ,   no_argument       , NULL , 12},
    { "help"            ,   required_argument , NULL , 13},
    { "prune"           ,   no_argument       , NULL , 14},
    { "depth"           ,   required_argument , NULL , 15},
    {"time-var"         ,   no_argument       , NULL , 16},
    {"function-type"    ,   required_argument , NULL , 17},
    {"offset"           ,   required_argument , NULL , 18},
    {"amplitude"        ,   required_argument , NULL , 19},
    {"period"           ,   required_argument , NULL , 20},
    {"constant"         ,   required_argument , NULL , 21},
    { NULL              ,    0                 , NULL , 0}
   };


typedef struct
{
    //index of desitination-vertex
    uint_fast64_t to;
    uint_fast64_t from;
    double weight;

}edge_t;

/*
connected holds the number of edges connected to vertex 
*/
typedef struct
{
    uint_fast64_t index;
    uint_fast32_t connected;
    edge_t * E;
}vertex_t;



/*
We assume, the verteces to be ordered in order from lowest index to highest.
The same is true for the edges, were we assume them to be ordered from lowest
to highest for the originating indeces ("edge_t.from")
*/
typedef struct
{
    vertex_t * V;
    uint_fast32_t num_vert;
    uint_fast32_t num_edge;
}network_t;

typedef struct
{
    network_t subnet;
    
    uint_fast64_t num;
    
    double value;
    double v_delta;

    uint_fast8_t *bit_array;

}network_state_t;

//type defines how closures are structured, by size (type = 0)
// or by number of neighbours (type = 1)
typedef struct
{
    uint_fast8_t type;
    //gives the number of nodes included in the closure
    uint_fast64_t size;
    uint_fast32_t dos;
    uint_fast32_t overlap;

}network_closure_t;


// used in case of line networks
typedef struct
{
    //Primary indentifiers
    uint_fast32_t width;
    uint_fast32_t num;
    uint_fast32_t first;

    //value of Expected Value (\tau_{first}...\tau_{first+width})
    double value;
    double orig_value;
    double v_delta;

    double k1, k2, k3, k4;

    //explicit bit representation of state 
    uint_fast8_t    *bit_array;

    //just for testing
    uint_fast8_t pruned; 
    uint_fast8_t approx_deps;

    double max_value;

} state_t;

typedef struct
{
    uint_fast32_t close_at;

    uint_fast32_t u_cl_num;    
    uint_fast32_t u_closures[64];
    uint_fast32_t u_cl_pos[64];
    
    uint_fast32_t l_closures[64];
    uint_fast32_t l_cl_pos[64];
    uint_fast32_t l_cl_num;

} closure_t;

typedef struct 
{
    int     time_var_found;
    int     funtyp_found;
    int     offset_found;
    int     period_found;
    int     amplit_found;
    int     constf_found;
} tv_set_inputs_t;

typedef struct 
{
    char*   funTyp_filename;
    char*   offset_filename;
    char*   period_filename;
    char*   amplit_filename;
    char*   constf_filename;
} tv_filenames_t;

typedef struct
{
    int*    funTyp_array; 
    double* offset_array;
    double* period_array;
    double* amplit_array;
    double* constf_array;
} tv_arrays_t;

static void 
usage(const char *progname);

int
gen_bit_array(state_t * state);


int
malloc_state_array(state_t **** st_array, size_t chain_length);
int
malloc_state_array_verbose(FILE * F, state_t **** st_array, size_t chain_length, closure_t closure);

int
init_state_array(state_t *** st_array, size_t chain_length, closure_t closure);
int
init_state_array_verbose(FILE * F, state_t *** st_array, size_t chain_length, closure_t closure);

int
preserve_called(state_t *** state_array, uint_fast32_t i, uint_fast32_t j, uint_fast32_t k, size_t chain_length);
int
preserve_called_depth(state_t *** state_array, uint_fast32_t i, uint_fast32_t j, uint_fast32_t k, size_t chain_length, uint_fast32_t depth);
int
approx_state_deps(state_t *** state_array, uint_fast32_t i, uint_fast32_t j, uint_fast32_t k);
int
prune_4332(state_t *** state_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j);
int
close_4_33_2(state_t *** state_array, uint_fast32_t chain_length);
double
closure_change(state_t *** state_array, closure_t * closure, state_t * state, uint_fast32_t j);
double
get_cl_change(state_t *** state_array, double * jr_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j, closure_t * closure);

static inline
int
free_state_array_full(state_t *** st_array, size_t chain_length, uint_fast32_t close_at);

void
print_state(state_t * state);


uint_fast8_t
get_bit(uint_fast32_t num, state_t * state);

///Get pointer to state state_arr[i][k]
state_t *
get_state(state_t *** state_arr, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j);

double
get_state_value(state_t *** state_arr, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j);

double 
get_change(state_t *** state_array, double * jr_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j, closure_t * closure);
double 
get_approx_change(state_t *** state_array, double * jr_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j);
double 
get_4332_change(state_t *** state_array, double * jr_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j);

int
solve_rk4_step(state_t *** state_array, double * jr_array, double h, size_t chain_length, closure_t * closure);

int
solve_rk4_step_verbose(FILE * F, state_t *** state_array, double * jr_array, double h, size_t chain_length, closure_t * closure);

int
solve_euler_step(state_t *** state_array, double * jr_array, double h, size_t chain_length, closure_t * closure);

int
solve_euler_step_verbose(FILE * F, state_t *** state_array, double * jr_array, double h, size_t chain_length, closure_t * closure);

int
solve_rk4(FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays);

int
solve_rk4_verbose(FILE * logfile, FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays);

int
solve_euler(FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays);

int
solve_euler_verbose(FILE * logfile, FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays);



int
save_state_fvalue(FILE * F, state_t *** state_array, size_t chain_length, double t, uint_fast32_t close_at);

int
save_state_singles(FILE * F, state_t *** state_array, size_t chain_length, double t);

int 
save_state_full(FILE * F, state_t *** state_array, size_t chain_length, double t);

int
save_state_col_nums(FILE * F, size_t chain_length);

uint_fast64_t
prune_array(state_t *** state_array, double * jr_array, size_t chain_length, closure_t * closure);

uint_fast32_t
prune_state(        state_t *** state_array, 
                    double * jr_array, 
                    uint_fast32_t chain_length, 
                    uint_fast32_t i, 
                    uint_fast32_t k, 
                    uint_fast32_t j,
                    closure_t * closure);

uint_fast32_t
unprune_closure_states(state_t *** state_array, double * jr_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j, closure_t * closure);

uint_fast32_t
closure_prune(state_t *** state_array, closure_t * closure, state_t * state, uint_fast32_t j);

uint_fast32_t
closure_prune_check(state_t *** state_array, closure_t * closure, state_t * state, uint_fast32_t j);

double
check_consistency(state_t *** state_array, uint_fast32_t state_size, uint_fast64_t chain_length, closure_t closure);

uint_fast8_t
closure_rhs(state_t *** state_array, closure_t closure, state_t *state, uint_fast32_t j, uint_fast64_t i_c, uint_fast64_t k_c, uint_fast64_t j_c);

uint_fast8_t
is_on_rhs(uint_fast64_t * lhs, uint_fast64_t * candidate, state_t *** state_array, closure_t closure, uint_fast64_t chain_length);

void
save_jacobi_matrix(state_t *** state_array, uint_fast64_t chain_length, closure_t closure, double * jr_array);

double
rhs_jacobi_eval(uint_fast64_t * lhs_indices, uint_fast64_t * rhs_indices, state_t *** state_array, closure_t closure, uint_fast64_t chain_length, double * jr_array);

#endif
