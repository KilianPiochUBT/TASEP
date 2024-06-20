#include "mf_full.h"


static void usage(const char *progname)
{
  printf(   "Usage: %s\n"
            "---Switches for hop rate control [Steady State Analysis]---\n"
            "   -a --alpha              Entry hop rate\n"
            "   -b --beta               Exit hop rate\n"
            "   -g --gamma              Inside flow rate\n"
            "---Switches for chain config---\n"
            "   -n --num-sites          Number of sites in the chain\n"
            "---Switches for I/O---\n"
            "   -i --input              Give file of hop rates. Must contain n+1 entries of type float seperated by whitespace\n"
            "   -o --output             Give filename of output file. Results will be saved to that location instead of stdout\n"
            "   --trajectories          Enable recording of trajectories of ALL state values. Results will be printed to <file> specified by -o\n"
            "---Switches for solver-options---\n"
            "   -m --method             Method used by solver [rk4|e]\n"
            "   -p --parallel-threads   Specify number of threads the calculation is allowed to occupy\n"
            "   -t --start              Start time\n"
            "   -T --stop               Stop time\n"
            "   -s --steps              Number of timesteps\n"
            "   -m --solver             Solving method to use [e(uler)|rk4]"
            "   -v --verbose            Verbose output\n"
            "---Switches for time variant simulations---\n"
            "--amplitude <char* filename> reads amplitudes from given file and sets functions, when using <lin> these serve as the gradient\n"
            "--period <char* filename> reads period multipliers from given file an sets jump rate functions. (default is 1)\n"
            "--offset <char* filename> reads offsets from given file and sets jump rate functions. (default is 0)\n"
            "--constant <char* filename> reads constants from given file and sets jump rate functions. (default is 0)\n"
            "--thread-num <char* number_of_threads_to_use> sets number of threads used during parallel portion of the code. Only applies if -f is also selected. (default is max Threads)\n"
            "--save-select <char* shorts for what to save> shorts are: 'a': all 'e': expected values 'p': probabilities 'c': covariance Matrices 'm': multiplicative expected values\n"
            "\n"
            "---Help---\n"
            "   -h --help -?            Display this dialog\n", progname);
}

//calculate and state description S_{k,n}
//with binary representation of k
//if k is too large to fit in n bits of data, 2 is returned. 
//if an allocation error is encountered -1 is returned. 
int 
gen_bit_array(state_t * state)
{   
    uint_fast32_t bitset;
    for(unsigned i = 0; i < state->width; i++){
        bitset = state->num;
        if((bitset >> i)&1)
            state->bit_array[i] = 1;
        else 
            state->bit_array[i] = 0;

    }

    return 0;
} 

int 
net_gen_bit_array(network_state_t * nstate)
{
    uint_fast64_t bitset;

    for(unsigned i = 0; i < nstate->subnet.num_vert; i++)
    {
        bitset = nstate->num;
        if((bitset >> i)&1)
            nstate->bit_array[i] = 1;
        else         
            nstate->bit_array[i] = 0;
    }    
    return 0;
}


//this is probably overkill
//expects pointer of type state_t**** 
int
malloc_state_array(state_t **** state_array, size_t chain_length)
{
        //Array elements will be accessed in the following way://

    /*                                             j = 
                                            <----first---->
    state_array ->  [0]     ->  [0]     ->  [0][1]...[n-1-i]
                |    k = num|
           /\   |           ->  [k-1]   ->  [0][1]...[n-1-i]
            w   | !!!k<POW2(i+1)!!!
            i   ->  [1]
            d   |
            t   |
            h   |
                ->  [2]
               .
               .
           =   .
                ->  [n-3]
                |
                |
            i   |
           \/   ->  [n-2]
                |
                |
                ->  [n=chain_length-1] 

    */

    //sanity check
    //because k = {0..POW2(i+1)-1}, i = {0..chain_length-1}

    
        *state_array = (state_t ***) malloc(chain_length * sizeof(state_t**));
        if(state_array == NULL)
            return -1;

        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            (*state_array)[i] = (state_t **)malloc((size_t)POW2(i+1) * sizeof(state_t*));
            if((*state_array)[i] == NULL)
                return -1;

                for(uint_fast32_t k = 0; k < POW2(i+1); k++)
                {   
                    (*state_array)[i][k] = (state_t *)malloc((size_t)(chain_length-i)*sizeof(state_t));
                    if((*state_array)[i][k] == NULL)
                        return -1;

                    for(uint_fast32_t j = 0; j < chain_length-i; j++){
                        (*state_array)[i][k][j].bit_array = (uint_fast8_t *) malloc((size_t)(i+1) * sizeof(uint_fast8_t));
                        if((*state_array)[i][k][j].bit_array == NULL)
                            return -1;
                    }

            }
        }
    
    return 0;

}

int
malloc_state_array_verbose(FILE * F, state_t **** state_array, size_t chain_length, closure_t closure)
{
    
        *state_array = (state_t ***) malloc(chain_length * sizeof(state_t**));
        if(state_array == NULL)
            return -1;
    
        for(uint_fast32_t i = 0; i < closure.close_at; i++)
        {
            fprintf(F, "Allocating state array (%"PRIuFAST32"/%ld)\r", i+1, chain_length);
            (*state_array)[i] = (state_t **)malloc((size_t)POW2(i+1) * sizeof(state_t*));
            if((*state_array)[i] == NULL)
                return -1;
         
                for(uint_fast32_t k = 0; k < POW2(i+1); k++)
                {   
                    (*state_array)[i][k] = (state_t *)malloc((size_t)(chain_length-i)*sizeof(state_t));
                    if((*state_array)[i][k] == NULL)
                        return -1;
    
                    for(uint_fast32_t j = 0; j < chain_length-i; j++){
                        (*state_array)[i][k][j].bit_array = (uint_fast8_t *) malloc((size_t)(i+1) * sizeof(uint_fast8_t));
                        if((*state_array)[i][k][j].bit_array == NULL)
                            return -1;
                    }
    
            }
        }
    
    
    fprintf(F, "Finished allocating state arrays\n");
    return 0;
}


//set starting value to "naive" starting value.
//the remaining values are set according to array position
int
init_state(state_t * state, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j)
{
    state->width = i+1;
    state->num = k;
    state->first = j;
    state->value = pow(0.5, i+1);
    state->orig_value = pow(0.5, i+1);
    state->v_delta = 0;
    state->pruned = 0;
    state->approx_deps = 0;
    state->max_value = 1.0;
    gen_bit_array(state);

    state->pruned = 1;
    return 0;
}

int
init_state_custom(state_t * state, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j, double * array)
{
    state->width = i+1;
    state->num = k;
    state->first = j;
   //state->orig_value = pow(0.5, i+1);
    gen_bit_array(state);
    if(state->bit_array[0]){
        state->value = array[j];
        state->orig_value = array[j];
    } else {
        state->value = 1-array[j];
        state->orig_value = 1-array[j];
    }
    for(uint_fast32_t l = j+1; l < j+i+1; l++){
        if(state->bit_array[l] == 1){
            state->value*=array[l];
            state->orig_value*=array[l];
        } else {
            state->value*=(1-array[l]);
            state->orig_value*=(1-array[l]);
        }
    }
    //printf("%lf ", state->value);
    state->v_delta = 0;
    state->pruned = 0;
    state->approx_deps = 0;
    state->max_value = 1.0;

    state->pruned = 1;
    return 0;

}

int 
prepare_time_var_arrays(tv_arrays_t * tv_arrays, tv_filenames_t * tv_filenames, tv_set_inputs_t * tv_inputs, int n)
{
    int err = 0;
    FILE * F;
    if(tv_inputs->funtyp_found){
        if(fopen(tv_filenames->funTyp_filename, "r") == NULL){
            printf("ERROR: Cannot open file containing function Types. Aborting\n");
            return -1;
        }
        tv_arrays->funTyp_array = (int *)malloc(((size_t)n+1) * sizeof(int));
        if(tv_arrays->funTyp_array == NULL){
            printf("ERROR: Not enough space to allocate function type array. Aborting\n");
            return -1;
        }
        F = fopen(tv_filenames->funTyp_filename, "r");
        char funTypSelector[256];
        for(int i = 0; i < n+1; i++){
            err = fscanf(F, "%255s", funTypSelector);
            if(err == EOF)
                return err;
            switch(funTypSelector[0]){
                case 's':
                    tv_arrays->funTyp_array[i] = 0;
                    break;
                case 'c':
                    tv_arrays->funTyp_array[i] = 1;
                    break;
                case 'l':
                    tv_arrays->funTyp_array[i] = 2;
                    break;
                default:
                    printf("ERROR: Non identifiable function type fonud in"
                           "%s. Read %s. Aborting!\n", tv_filenames->funTyp_filename, funTypSelector);
            }
        }
        printf("Read function types successfully\n");
    }
    else{
        tv_arrays->funTyp_array = (int *)malloc(((size_t)n+1) * sizeof(int));
        if(tv_arrays->funTyp_array == NULL){
            printf("ERROR: Not enough space to allocate function type array. Aborting\n");
            return -1;
        }
        for(int i = 0; i < n+1; i++){
            tv_arrays->funTyp_array[i] = 0;
        }
        
        printf("Wrote function type array successfully\n");
    }

    if(tv_inputs->offset_found){
        //read offset values
        if(fopen(tv_filenames->offset_filename, "r") == NULL){
            printf("ERROR: Cannot open file containing offset values. Aborting\n");
            return -1;
        }
        tv_arrays->offset_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->offset_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        F = fopen(tv_filenames->offset_filename, "r");
        for(int i = 0; i < n+1; i++){
            err = fscanf(F, "%lf, ", &tv_arrays->offset_array[i]);
            if(err == EOF)
                return err;
        }

                printf("Read offset values successfully\n");
    }
    else{
        tv_arrays->offset_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->offset_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        for(int i = 0; i < n+1; i++)
            tv_arrays->offset_array[i] = 0.;

            printf("Wrote offset values sucessfully\n");
        }

    if(tv_inputs->amplit_found){
        //read amplitude values
        if(fopen(tv_filenames->amplit_filename, "r") == NULL){
            printf("ERROR: Cannot open file containing amplitude values. Aborting\n");
            return -1;
        }
        tv_arrays->amplit_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->amplit_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        F = fopen(tv_filenames->amplit_filename, "r");
        for(int i = 0; i < n+1; i++){
            err = fscanf(F, "%lf, ", &tv_arrays->amplit_array[i]);
            if(err == EOF)
                return err;
        }

        printf("Read amplitude values correctly\n");
    }
    else{
        tv_arrays->amplit_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->amplit_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        for(int i = 0; i < n+1; i++)
            tv_arrays->amplit_array[i] = 1.;
        printf("Wrote amplitude values correctly\n");
    }

    if(tv_inputs->period_found){
        //read amplitude values
        if(fopen(tv_filenames->period_filename, "r") == NULL){
            printf("ERROR: Cannot open file containing amplitude values. Aborting\n");
            return -1;
        }
        tv_arrays->period_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->period_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        F = fopen(tv_filenames->period_filename, "r");
        for(int i = 0; i < n+1; i++){
            err = fscanf(F, "%lf, ", &tv_arrays->period_array[i]);
            if(err == EOF)
                return -1;
        }

        printf("Read period values correctly\n");
    }
    else{
        tv_arrays->period_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->period_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        for(int i = 0; i < n+1; i++)
            tv_arrays->period_array[i] = 1.;

        printf("Wrote period values correctly\n");
    }

    if(tv_inputs->constf_found){
        //read amplitude values
        if(fopen(tv_filenames->constf_filename, "r") == NULL){
            printf("ERROR: Cannot open file containing amplitude values. Aborting\n");
            return -1;
        }
        tv_arrays->constf_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->constf_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        F = fopen(tv_filenames->constf_filename, "r");
        for(int i = 0; i < n+1; i++){
            err = fscanf(F, "%lf, ", &tv_arrays->constf_array[i]);
            if(err == EOF)
                return err;
        }

        printf("Read constant values correctly\n");
    }
    else{
        tv_arrays->constf_array = (double *)malloc(((size_t)n+1) * sizeof(double));
        if(tv_arrays->constf_array == NULL){
            printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
            return -1;
        }
        for(int i = 0; i < n+1; i++)
            tv_arrays->constf_array[i] = 1.;

        printf("Wrote constant values correctly\n");
    }
    return 0;
}

static inline 
int
free_state_array_full(state_t *** state_array, size_t chain_length, uint_fast32_t close_at)
{
    /*if(close_at > 0 && close_at < chain_length){
        for(uint_fast32_t i = 0; i < close_at; i++)
        {
            for(uint_fast32_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length - i; j++){
                    if(state_array[i][k][j].pruned == 0)
                        free(state_array[i][k][j].bit_array);
                }
                free(state_array[i][k]);
            }
            free(state_array[i]);
        }
        free(state_array);
    }
    else{*/
        for(uint_fast32_t i = 0; i < close_at; i++)
        {
            for(uint_fast32_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length - i; j++){
                    if(state_array[i][k][j].pruned == 0)
                        free(state_array[i][k][j].bit_array);
                }
                free(state_array[i][k]);
            }
            free(state_array[i]);
        }
        free(state_array);
    //}
    return 0;
}

void
print_state(state_t * state)
{

    printf("Width: %lu, State number: %"PRIuFAST32", First Element: %"PRIuFAST32"\n", state->width, state->num, state->first);
    printf("Value: %lf\n", state->value);
    for(uint_fast32_t i = 0; i < state->width; i++)
    {
        printf("%"PRIuFAST8" ", state->bit_array[i]);
    }
    printf("\n");
}

uint_fast8_t
get_bit(uint_fast32_t num, state_t * state)
{
    
    return state->bit_array[num];

}

//get porinter at position [i][k][j] from state_arr
state_t *
get_state(state_t *** state_array, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j){
    
    return &(state_array[i][k][j]);
}

int
free_pruned_entries(state_t *** state_array, uint_fast32_t chain_length, closure_t closure){

    for(uint_fast32_t i = 0; i < closure.close_at; i++)
    {
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                if(state_array[i][k][j].pruned == 1){
                    free(state_array[i][k][j].bit_array);
                }
            }
        }
    }
    return 0;
}

uint_fast64_t
prune_array(state_t *** state_array, double * jr_array, size_t chain_length, closure_t * closure)
{
    uint_fast64_t num_pruned = 0, num_pruned_prev = 0;
    printf("Pruning unused Point-Function array entries\n");
    state_array[0][0][chain_length-1].pruned = 0;

    for(uint_fast32_t i = 0; i < chain_length; i++)
        state_array[0][1][i].pruned = 0;

    num_pruned = chain_length+1;

    while(num_pruned > num_pruned_prev || num_pruned_prev == 0){
        if(num_pruned != 0)
            num_pruned_prev = num_pruned;
        for(uint_fast32_t i = 0; i < closure->close_at; i++)
        {
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    if(state_array[i][k][j].pruned == 0)
                        num_pruned += prune_state(state_array, jr_array, chain_length, i, k, j, closure);
                }
            }
        }
    }
    return num_pruned;
}

int
approx_state_deps(state_t *** state_array, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j)
{
    // very simple first implementation: only singles will be called. 
    state_t * state = get_state(state_array, i, k, j);
    size_t state_chain_l = state_array[i][k][j].width;


    for(size_t b = 0; b < state_chain_l; b++)
    {
        if(!get_bit(b, state))
        {
            state_array[0][0][state->first+b].pruned = 0;
        }
    }
    return 0;
}


double 
pseudo_prune_array(state_t *** state_array, uint_fast32_t chain_length, uint_fast32_t close_at){
    state_t * state;
    if(close_at > 0){
        for(uint_fast32_t i = 0; i < close_at; i++)
        {
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state = get_state(state_array, i, k, j);
                    state->pruned = 0;
                }
            }
        }
    }
    else{
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state = get_state(state_array, i, k, j);
                    state->pruned = 0;
                }
            }
        }
    }
    return 0;
}

double 
get_change( state_t *** state_array, 
            double * jr_array, 
            uint_fast32_t chain_length, 
            uint_fast32_t i, 
            uint_fast32_t k, 
            uint_fast32_t j,
            closure_t * closure){
    
           
    state_t * state;
    double change = 0, change_temp = 0;

    uint_fast32_t cl_type = 0;

    state = get_state(state_array, i, k, j);
    uint_fast32_t le = state->width;
    
    
    if(state_array[i][k][j].width+1 == closure->close_at){
        change = get_cl_change(state_array, jr_array, chain_length, i, k, j, closure);
        return change;
    }

        
        
    if(state->pruned == 0){
        if((k >> i)&1){
            //state is result of jump in
            if(j+i == chain_length-1){
                change += jr_array[chain_length] * (state_array[i][k-POW2(i)][j].value);
            }else {
                change += jr_array[j+i+1] * (state_array[i+1][k+POW2(i)][j].value);
            }
        } else { 
            //state is result of no(!) jump in
            if(j+i == chain_length-1){
                //state stays the same, no jump from left position in chain
                change -= jr_array[chain_length] * (state_array[i][k][j].value);
            } else {
                change -= jr_array[j+i+1] * (state_array[i+1][k+POW2(i+1)][j].value);
            }
        }

        if(!(get_bit(0, state))){
            //state is result of jump out
            if(j > 0){
                change += jr_array[j] * (state_array[i+1][2*(k+1)][j-1].value);
            }else{
                change += jr_array[0] * (state_array[i][k+1][j].value);
            }
        } else { 
            //state is result of no(!) jump out
            if(j>0){ 
                change -= jr_array[j] * (state_array[i+1][k*2][j-1].value); 
            } else {
                change -= jr_array[0] * (state_array[i][k][j].value);
            }
        }

        //jumps inside the state
        for(uint_fast32_t l = 0; l < i; l++)
        {
            if( !(get_bit(l+1, state)) && (get_bit(l, state)) )
            {
                change += jr_array[j+l+1] * (state_array[i][k + POW2(l)][j].value);
            }

            if( get_bit(l+1, state) && !(get_bit(l, state)) )
            {
                change -= jr_array[j+l+1] * (state_array[i][k][j].value);
            }
        }
    }
    return change;
}

uint_fast32_t
prune_state(state_t *** state_array, 
            double * jr_array, 
            uint_fast32_t chain_length, 
            uint_fast32_t i, 
            uint_fast32_t k, 
            uint_fast32_t j,
            closure_t * closure)
{   
    uint_fast32_t num_pruned = 0;

    state_t * state;

    uint_fast32_t cl_type = 0;

    state = get_state(state_array, i, k, j);
    uint_fast32_t le = state->width;
    
    
    if(state_array[i][k][j].width+1 == closure->close_at){
        num_pruned =  unprune_closure_states(state_array, jr_array, chain_length, i, k, j, closure);
        return num_pruned;
    } else if(state_array[i][k][j].width+1 > closure->close_at){
        return 0;
    }

        
        
    if(state->pruned == 0){
        if((k >> i)&1){
            //state is result of jump in
            if(j+i == chain_length-1){
                if(state_array[i][k-POW2(i)][j].pruned == 1){
                    state_array[i][k-POW2(i)][j].pruned = 0;
                    num_pruned++;
                }
            }else {
                if(state_array[i+1][k+POW2(i)][j].pruned == 1){
                    state_array[i+1][k+POW2(i)][j].pruned = 0;
                    num_pruned++;
                }
            }
        } else { 
            //state is result of no(!) jump in
            if(j+i != chain_length-1){
                if(state_array[i+1][k+POW2(i+1)][j].pruned == 1){
                    state_array[i+1][k+POW2(i+1)][j].pruned = 0;
                    num_pruned++;
                }
            }
        }

        if(!get_bit(0,state)){
            //state is result of jump out
            if(j > 0){
                if(state_array[i+1][2*(k+1)][j-1].pruned == 1){
                    state_array[i+1][2*(k+1)][j-1].pruned = 0;
                    num_pruned++;
                }
            }else{
                if(state_array[i][k+1][j].pruned == 1){
                    state_array[i][k+1][j].pruned = 0;
                    num_pruned++;
                }
            }
        } else { 
            //state is result of no(!) jump out
            if(j > 0){ 
                if(state_array[i+1][k*2][j-1].pruned == 1){ 
                    state_array[i+1][k*2][j-1].pruned = 0;
                    num_pruned++; 
                }
            }
        }

        //jumps inside the state
        for(uint_fast32_t l = 0; l < i; l++)
        {
            if(  !(get_bit(l+1, state)) && (get_bit(l, state)) )
            {
                if(state_array[i][k + POW2(l)][j].pruned == 1){
                    state_array[i][k + POW2(l)][j].pruned = 0;
                    num_pruned++;
                }
            }

        }
    }

    return num_pruned;
}

uint_fast32_t
unprune_closure_states(state_t *** state_array, double * jr_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j, closure_t * closure)
{
    //Preliminary stuff
    state_t * state = get_state(state_array, i, k, j);
    uint_fast32_t num_pruned = 0;
    if(state->width > closure->close_at)
        return num_pruned;

    //Actually calculate closed change
    if((k >> i)&1){
        //state is result of jump in
        if(j+i == chain_length-1){
            if(state_array[i][k-POW2(i)][j].pruned == 1){
                state_array[i][k-POW2(i)][j].pruned = 0;
                num_pruned++;
            }
        }else {
            if(closure_prune_check(state_array, closure, &state_array[i+1][k+POW2(i)][j], j) > 1){
                num_pruned+=closure_prune(state_array, closure, &state_array[i+1][k+POW2(i)][j], j);
            }
        }
    } else { 
        //state is result of no(!) jump in
        if(j+i != chain_length-1){
            if(closure_prune_check(state_array, closure, &state_array[i+1][k+POW2(i+1)][j], j) > 1){
                num_pruned += closure_prune(state_array, closure, &state_array[i+1][k+POW2(i+1)][j], j);
            }
        }
    }

    if(!get_bit(0,state)){
        //state is result of jump out
        if(j > 0){
            if(closure_prune_check(state_array, closure, &state_array[i+1][2*(k+1)][j-1], j-1) > 1){
                num_pruned+=closure_prune(state_array, closure, &state_array[i+1][2*(k+1)][j-1], j-1);
            }
        }else{
            if(state_array[i][k+1][j].pruned == 1){
                state_array[i][k+1][j].pruned = 0;
                num_pruned++;
            }
        }
    } else { 
        //state is result of no(!) jump out
        if(j > 0){ 
            if(closure_prune_check(state_array, closure, &state_array[i+1][k*2][j-1], j-1) > 1){ 
                num_pruned+=closure_prune(state_array, closure, &state_array[i+1][k*2][j-1], j-1); 
            }
        }
    }

    //jumps inside the state
    for(uint_fast32_t l = 0; l < i; l++)
    {
        if(  !(get_bit(l+1, state)) && (get_bit(l, state)) )
        {
            if(state_array[i][k + POW2(l)][j].pruned == 1){
                state_array[i][k + POW2(l)][j].pruned = 0;
                num_pruned++;
            }
        }

    }

    return num_pruned;
}

uint_fast32_t
closure_prune(state_t *** state_array, closure_t * closure, state_t * state, uint_fast32_t j)
{   
    //Change this: Too slow (probably)
    uint_fast32_t * index_u, * index_l;
    index_u = (uint_fast32_t *)malloc((size_t)(closure->u_cl_num) * sizeof(uint_fast32_t));
    index_l = (uint_fast32_t *)malloc((size_t)(closure->l_cl_num) * sizeof(uint_fast32_t));
    //------------------------------//

    uint_fast32_t pruned = 0;

    state->pruned = 0;

    if(index_l == NULL || index_u == NULL)
        err_exit("ERROR: Could not allocate index arrays for arbitrary closure->", 1);

    for(uint_fast32_t l = 0; l < closure->u_cl_num; l++){
        index_u[l] = 0;
        for(uint_fast32_t m = 0; m < closure->u_closures[l]; m++){
            index_u[l] += POW2(m) * state->bit_array[closure->u_cl_pos[l]+m];
        }
    }

    for(uint_fast32_t l = 0; l < closure->l_cl_num; l++){
        index_l[l] = 0;
        for(uint_fast32_t m = 0; m < closure->l_closures[l]; m++){
            index_l[l] += POW2(m) * state->bit_array[closure->l_cl_pos[l]+m];
        }
    }

    double change = 1.;

    for(uint_fast32_t i = 0; i < closure->u_cl_num; i++){
        if(state_array[closure->u_closures[i]-1][index_u[i]][closure->u_cl_pos[i]+j].pruned == 1){
            state_array[closure->u_closures[i]-1][index_u[i]][closure->u_cl_pos[i]+j].pruned = 0; 
            pruned++;
        }
    }

    for(uint_fast32_t i = 0; i < closure->l_cl_num; i++){
        if(state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].pruned == 1){
            state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].pruned = 0;
            pruned++;
        }
    }
    free(index_l);
    free(index_u);
    return pruned;
}

uint_fast32_t
closure_prune_check(state_t *** state_array, closure_t * closure, state_t * state, uint_fast32_t j)
{   
    //Change this: Too slow (probably)
    uint_fast32_t * index_u, * index_l;
    index_u = (uint_fast32_t *)malloc((size_t)(closure->u_cl_num) * sizeof(uint_fast32_t));
    index_l = (uint_fast32_t *)malloc((size_t)(closure->l_cl_num) * sizeof(uint_fast32_t));
    //------------------------------//

    state->pruned = 0;

    uint_fast32_t pruned = 1;

    if(index_l == NULL || index_u == NULL)
        err_exit("ERROR: Could not allocate index arrays for arbitrary closure->", 1);

    for(uint_fast32_t l = 0; l < closure->u_cl_num; l++){
        index_u[l] = 0;
        for(uint_fast32_t m = 0; m < closure->u_closures[l]; m++){
            index_u[l] += POW2(m) * state->bit_array[closure->u_cl_pos[l]+m];
        }
    }

    for(uint_fast32_t l = 0; l < closure->l_cl_num; l++){
        index_l[l] = 0;
        for(uint_fast32_t m = 0; m < closure->l_closures[l]; m++){
            index_l[l] += POW2(m) * state->bit_array[closure->l_cl_pos[l]+m];
        }
    }

    double change = 1.;

    for(uint_fast32_t i = 0; i < closure->u_cl_num; i++){
        if(state_array[closure->u_closures[i]-1][index_u[i]][closure->u_cl_pos[i]+j].pruned == 1){
            pruned++;
        }
    }

    for(uint_fast32_t i = 0; i < closure->l_cl_num; i++){
        if(state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].pruned == 1){
            pruned++;
        }
    }
    free(index_l);
    free(index_u);
    return pruned;
}

uint_fast8_t
check_closure(closure_t * closure)
{
    //TODO: VERY IMPORTANT
    //check for well definedness of closure
    return 0;
}

double
get_cl_change(state_t *** state_array, double * jr_array, uint_fast32_t chain_length, uint_fast32_t i, uint_fast32_t k, uint_fast32_t j, closure_t * closure)
{
    //Preliminary stuff
    state_t * state = get_state(state_array, i, k, j);

    if(state->width > closure->close_at)
        return 0;

    double change = 0;

    //Actually calculate closed change

        //state is result of jump in
        if((k >> i)&1){
            //state is result of jump in
            if(j+i == chain_length-1){
                change += jr_array[chain_length] * (state_array[i][k-POW2(i)][j].value);
            }else {
                change += jr_array[j+i+1] * closure_change(state_array, closure, &state_array[i+1][k+POW2(i)][j], j);
            }
        } else { 
            //state is result of no(!) jump in
            if(j+i == chain_length-1){
                //state stays the same, no jump from left position in chain
                change -= jr_array[chain_length] * (state_array[i][k][j].value);
            } else {
                change -= jr_array[j+i+1] * closure_change(state_array, closure, &state_array[i+1][k+POW2(i+1)][j], j);
            }
        }

        if(!(get_bit(0, state))){
            //state is result of jump out
            if(j > 0){
                change += jr_array[j] * closure_change(state_array, closure, &state_array[i+1][2*(k+1)][j-1], j-1);
            }else{
                change += jr_array[0] * (state_array[i][k+1][j].value);
            }
        } else { 
            //state is result of no(!) jump out
            if(j>0){ 
                change -= jr_array[j] * closure_change(state_array, closure, &state_array[i+1][k*2][j-1], j-1); 
            } else {
                change -= jr_array[0] * (state_array[i][k][j].value);
            }
        }

        //jumps inside the state
        for(uint_fast32_t l = 0; l < i; l++)
        {
            if( !(get_bit(l+1, state)) && (get_bit(l, state)) )
            {
                change += jr_array[j+l+1] * (state_array[i][k + POW2(l)][j].value);
            }

            if( get_bit(l+1, state) && !(get_bit(l, state)) )
            {
                change -= jr_array[j+l+1] * (state_array[i][k][j].value);
            }
        }
    
    return change;
}

double
closure_change(state_t *** state_array, closure_t * closure, state_t * state, uint_fast32_t j)
{
    //Change this: Too slow (probably)
    uint_fast32_t * index_u, * index_l;
    index_u = (uint_fast32_t *)malloc((size_t)(closure->u_cl_num) * sizeof(uint_fast32_t));
    index_l = (uint_fast32_t *)malloc((size_t)(closure->l_cl_num) * sizeof(uint_fast32_t));
    //------------------------------//
    
    if(index_l == NULL || index_u == NULL)
        err_exit("ERROR: Could not allocate index arrays for arbitrary closure->", 1);

    for(uint_fast32_t l = 0; l < closure->u_cl_num; l++){
        index_u[l] = 0;
        for(uint_fast32_t m = 0; m < closure->u_closures[l]; m++){
            index_u[l] += POW2(m) * state->bit_array[closure->u_cl_pos[l]+m];
        }
    }

    for(uint_fast32_t l = 0; l < closure->l_cl_num; l++){
        index_l[l] = 0;
        for(uint_fast32_t m = 0; m < closure->l_closures[l]; m++){
            index_l[l] += POW2(m) * state->bit_array[closure->l_cl_pos[l]+m];
        }
    }

    double change = 1.;

    for(uint_fast32_t i = 0; i < closure->u_cl_num; i++){
        if(state_array[closure->u_closures[i]-1][index_u[i]][closure->u_cl_pos[i]+j].value < 10e-20)
            change = 0;
        else 
            change = change * state_array[closure->u_closures[i]-1][index_u[i]][closure->u_cl_pos[i]+j].value; 
    }

    for(uint_fast32_t i = 0; i < closure->l_cl_num; i++){
        if(state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].value < 10e-20)
            change = 0;
        else 
            change = change / state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].value;
    }
    free(index_l);
    free(index_u);

    return change;
}

 
int 
solve_rk4_step(state_t *** state_array, double * jr_array, double h, size_t chain_length, closure_t * closure){

   
    if(closure->close_at > 0){
        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);
                    if(state->pruned == 0){
                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                        state->orig_value = state->value;
                        state->k1 = value_delta;
                        state->v_delta = (h * 0.5 * state->k1);
                        state->k1 = h * state->k1;
                    }           
                }
            }
        }

        //update values
        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    
                    state_t *state = get_state(state_array, i,k,j);
                    if(state->pruned == 0)
                        state->value = state->orig_value+state->v_delta;
                }
            }
        }

        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);
                    if(state->pruned == 0){

                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);

                        state->k2 = value_delta;
                        state->v_delta = (h * 0.5 * state->k2);
                        state->k2 = h * state->k2;
                    }           
                }
            }
        }

        //update values
        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {   
                    state_t *state = get_state(state_array, i,k,j);
                    if(state->pruned == 0)
                        state->value = state->orig_value+state->v_delta;
                }
            }
        }

        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);

                    if(state->pruned == 0){
                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                        
                        state->k3 = value_delta;
                        state->v_delta = (h * state->k3);
                        state->k3 = h * state->k3;
                    }           
                }
            }
        }

        //update values
        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t *state = get_state(state_array, i,k,j);
                    if(state->pruned == 0)
                        state->value = state->orig_value+state->v_delta;
                }
            }
        }

        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);

                    if(state->pruned == 0){

                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);

                        state->k4 = value_delta;
                        state->k4 = h * state->k4;

                        state->v_delta = (state->k1+2*state->k2+2*state->k3+state->k4)/6.;
                    }           
                }
            }
        }

        for(uint_fast32_t i = 0; i < closure->close_at-1; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);
                        if(state->pruned == 0){
                        state->value = state->orig_value+state->v_delta;

                        if(fabs(state->value) > fabs(state->max_value))
                            state->max_value = state->value;
                        if(state->value < 0)
                            state->value = 0;
                    }
                }
            }
        }
    }
    else{
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);

                    if(state->pruned == 0){

                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                        state->orig_value = state->value;


                        state->k1 = value_delta;
                        state->value = (h * 0.5 * state->k1) + state->orig_value;
                        state->k1 = h * state->k1;
                    }           
                }
            }
        }
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp prallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {   
                    state_t *state = get_state(state_array, i,k,j);
                    if(state->pruned == 0)
                        state->value = state->orig_value+state->v_delta;
                }
            }
        }
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);

                    if(state->pruned == 0){

                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);

                        state->k2 = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                        state->value = (h * 0.5 * state->k2) + state->orig_value;
                        state->k2 = h * state->k2;
                    }           
                }
            }
        }
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {   
                    state_t *state = get_state(state_array, i,k,j);
                    if(state->pruned == 0)
                        state->value = state->orig_value+state->v_delta;
                }
            }
        }
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);

                    if(state->pruned == 0){

                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);

                        state->k3 = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                        state->value = (h * state->k3) + state->orig_value;
                        state->k3 = h * state->k3;
                    }           
                }
            }
        }
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {   
                    state_t *state = get_state(state_array, i,k,j);
                    if(state->pruned == 0)
                        state->value = state->orig_value+state->v_delta;
                }
            }
        }
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);

                    if(state->pruned == 0){

                        double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);

                        state->k4 = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                        state->k4 = h * state->k4;

                        state->value = state->orig_value;

                        state->v_delta = (state->k1+2*state->k2+2*state->k3+state->k4)/6.;
                    }           
                }
            }
        }

        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            #pragma omp parallel for
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    state_t * state = get_state(state_array, i, k, j);
                    if(state->pruned == 0){
                        state->value += state->v_delta;

                        if(fabs(state->value) > fabs(state->max_value))
                            state->max_value = state->value;
                        if(state->value < 0)
                            state->value = 0;
                    }
                }
            }
        }
    }
    
    return 0;

}


int
solve_euler_step(state_t *** state_array, double * jr_array, double h, size_t chain_length, closure_t * closure){

    for(uint_fast32_t i = 0; i < chain_length; i++)
    {
        #pragma omp parallel for
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                state_t * state = get_state(state_array, i, k, j);
                if(state->pruned == 0){
                    double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                    state->v_delta = h * value_delta;
                }
            }
        }
    }

    for(uint_fast32_t i = 0; i < chain_length; i++)
    {
        #pragma omp parallel for
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                state_t * state = get_state(state_array, i, k, j);
                if(state->pruned == 0)
                    state->value += state->v_delta;
                if(state->value < 0)
                    state->value = 0;
            }
        }
    }

    
    return 0;

}

int
solve_euler_step_verbose(FILE * F, state_t *** state_array, double * jr_array, double h, size_t chain_length, closure_t * closure){

    for(uint_fast32_t i = 0; i < chain_length; i++)
    {
        #pragma omp parallel for
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                double value_delta = get_change(state_array, jr_array, chain_length, i, k, j, closure);
                state_t * state = get_state(state_array, i, k, j);
                if(state->pruned == 0)
                    state->v_delta = h * value_delta;
            }
        }
    }

    for(uint_fast32_t i = 0; i < chain_length; i++)
    {
        #pragma omp parallel for
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                state_t * state = get_state(state_array, i, k, j);
                if(state->pruned == 0)
                     state->value += state->v_delta;
                if(state->value < 0)
                    state->value = 0;
            }
        }
    }

    
    return 0;

}


int
solve_rk4(FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays){
    double h = (T - t) / (double)s;
    int err = 0;
    if(tra_steps < 0){
        for(uint_fast64_t i = 0; i < s; i++){
            //Stop gap solution -> refactor all of these to function pointers at some point
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);  
            solve_rk4_step(state_array, jr_array, h, chain_length, closure);
        }
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    } else {
        
        save_state_col_nums(F, chain_length);
        
        for(uint_fast64_t i = 0; i < s; i++){
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);
            solve_rk4_step(state_array, jr_array, h, chain_length, closure);

            if(i%tra_steps == 0)
                save_state_fvalue(F, state_array, chain_length, (double)i*h, closure->close_at);
        }
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    }

    return 0;
}

int
solve_rk4_verbose(FILE * logfile, FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays){
    double h = (T - t) / (double)s;
    int err = 0;
    if(tra_steps < 0){
        for(uint_fast64_t i = 0; i < s; i++){
            fprintf(logfile, "Calulating step (%"PRIuFAST64"/%"PRIuFAST64")\r", i, s);
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);
            solve_rk4_step(state_array, jr_array, h, chain_length, closure);
        }
        fprintf(logfile, "\n");
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    } else {
        
        save_state_col_nums(F, chain_length);
        
        for(uint_fast64_t i = 0; i < s; i++){
            fprintf(logfile, "Calulating step (%"PRIuFAST64"/%"PRIuFAST64")\r", i, s);
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);
            solve_rk4_step(state_array, jr_array, h, chain_length, closure);

            if(i%tra_steps == 0)
                save_state_fvalue(F, state_array, chain_length, (double)i*h, closure->close_at);
        }
        fprintf(logfile, "\n");
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    }
    return 0;
}

int
solve_euler(FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays){
    double h = (T - t) / (double)s;
    int err = 0;
    if(tra_steps < 0){
        for(uint_fast64_t i = 0; i < s; i++){
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);
            solve_euler_step(state_array, jr_array, h, chain_length, closure);
        }
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    } else {
        
        save_state_col_nums(F, chain_length);

        for(uint_fast64_t i = 0; i < s; i++){
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);
            solve_euler_step(state_array, jr_array, h, chain_length, closure);

            if(i%tra_steps == 0)
                save_state_fvalue(F, state_array, chain_length, (double)i*h, closure->close_at);
        }
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    }

    return 0;
}

int
solve_euler_verbose(FILE * logfile, FILE * F, int tra_steps, state_t *** state_array, double * jr_array, size_t chain_length, double t, double T, uint_fast64_t s, closure_t * closure, int time_var, tv_arrays_t * tv_arrays){
    double h = (T - t) / (double)s;
    int err = 0;
    if(tra_steps < 0){
        for(uint_fast64_t i = 0; i < s; i++){
            fprintf(logfile, "Calulating step (%"PRIuFAST64"/%"PRIuFAST64")\r", i, s);
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);
            solve_euler_step_verbose(logfile, state_array, jr_array, h, chain_length, closure);
        }
        fprintf(logfile, "\n");
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    } else {

        save_state_col_nums(F, chain_length);

        for(uint_fast64_t i = 0; i < s; i++){
            fprintf(logfile, "Calulating step (%"PRIuFAST64"/%"PRIuFAST64")\r", i, s);
            if(time_var)
                timeVariantProbs(jr_array, chain_length + 1, tv_arrays->period_array, tv_arrays->offset_array, tv_arrays->amplit_array, i*h, tv_arrays->funTyp_array, tv_arrays->constf_array);
            solve_euler_step_verbose(logfile, state_array, jr_array, h, chain_length, closure);

            if(i%tra_steps == 0)
                save_state_fvalue(F, state_array, chain_length, (double)i*h, closure->close_at);
        }
        fprintf(logfile, "\n");
        if(err)
            err_exit("ERROR: Issue while saving results of simulation", err);
    }
    return 0;
}

int
init_state_array(state_t *** state_array, size_t chain_length, closure_t closure){
    //#pragma omp parallel for not thread safe

    for(uint_fast32_t i = 0; i < closure.close_at; i++)
    {
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                init_state(&state_array[i][k][j], i, k, j);
            }
        }
    }
    return 0;
}

int
init_state_array_custom(state_t *** state_array, size_t chain_length, uint_fast32_t close_at, double * array){
    //#pragma omp parallel for not thread safe

    if(close_at > 0){

        for(uint_fast32_t i = 0; i < close_at; i++)
        {
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    init_state_custom(&state_array[i][k][j], i, k, j,array);
                }
            }
        }
    }
    else{
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    init_state_custom(&state_array[i][k][j], i, k, j,array);
                }
            }
        }
    }
    return 0;
}

int
init_state_array_verbose(FILE * F, state_t *** state_array, size_t chain_length, closure_t closure){
    fprintf(F, "Setting values for state array for chain-length of %ld\n", chain_length);
    
    for(uint_fast32_t i = 0; i < closure.close_at; i++)
    {
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                init_state(&state_array[i][k][j], i, k, j);
            }
        }
    }

    return 0;
}

int 
save_state_pruned(FILE * F, state_t *** state_array, size_t chain_length, double t, closure_t closure){

    fprintf(F, "%lf:\t", t);

    for(uint_fast32_t i = 0; i < closure.close_at; i++)
    {
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                if(state_array[i][k][j].pruned == 0){
                    for(unsigned int u = 0; u < j; u++){
                        fprintf(F, "x");
                    }
                    for(unsigned int u = j; u < j+i+1; u++){
                        fprintf(F, "%"PRIuFAST8"", state_array[i][k][j].bit_array[u-j]);
                    }
                    for(unsigned int u = j+i+1; u < chain_length; u++){
                        fprintf(F, "x");
                    }
                    //fprintf(F, " %lf ", state_array[i][k][j].value);
                    fprintf(F, " %d ", state_array[i][k][j].pruned);
                }
            }
            fprintf(F, "\n");
        }
    }
    fprintf(F, "\n\n");
    return 0;
}

int 
save_state_full(FILE * F, state_t *** state_array, size_t chain_length, double t){

    fprintf(F, "%lf:\t", t);

    for(uint_fast32_t i = 0; i < chain_length; i++)
    {
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {
                for(unsigned int u = 0; u < j; u++){
                    fprintf(F, "x");
                }
                for(unsigned int u = j; u < j+i+1; u++){
                    fprintf(F, "%"PRIuFAST8"", state_array[i][k][j].bit_array[u-j]);
                }
                for(unsigned int u = j+i+1; u < chain_length; u++){
                    fprintf(F, "x");
                }
                //fprintf(F, " %lf ", state_array[i][k][j].value);
                fprintf(F, " %d ", state_array[i][k][j].pruned);
            }
            fprintf(F, "\n");
        }
    }
    fprintf(F, "\n");
    
    return 0;
}


int 
save_state_fvalue(FILE * F, state_t *** state_array, size_t chain_length, double t, uint_fast32_t close_at){
    
    fprintf(F, "%8.8lf:\t", t);
    if(close_at > 0){
        for(uint_fast32_t i = 0; i < close_at; i++)
        {
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    fprintf(F, "%8.8lf\t", state_array[i][k][chain_length-j-1].value);
                }
                fprintf(F, "\t");
            }
        }
        fprintf(F, "\n");
    }
    else{
        for(uint_fast32_t i = 0; i < chain_length; i++)
        {
            for(uint_fast64_t k = 0; k < POW2(i+1); k++)
            {
                for(uint_fast32_t j = 0; j < chain_length -i; j++)
                {
                    fprintf(F, "%8.8lf\t", state_array[i][k][chain_length-j-1].value);
                }
                fprintf(F, "\t");
            }
        }
        fprintf(F, "\n");   
    }
    return 0;
}

int
save_state_col_nums(FILE * F, size_t chain_length){
    
    uint_fast64_t ctr = 1;
    fprintf(F, "#%"PRIuFAST64"\t\t\t", ctr);
    ctr++;
    for(uint_fast32_t i = 0; i < chain_length; i++)
    {
        for(uint_fast64_t k = 0; k < POW2(i+1); k++)
        {
            for(uint_fast32_t j = 0; j < chain_length -i; j++)
            {   
                if(ctr < 1000)
                    fprintf(F, "%"PRIuFAST64"\t\t\t", ctr);
                else
                    fprintf(F, "%"PRIuFAST64"\t\t", ctr);
                ctr++;
            }
            fprintf(F, "\t");
        }
    }
   
    fprintf(F, "\n");
    return 0;
}

int
save_state_singles(FILE * F, state_t *** state_array, size_t chain_length, double t)
{
    for(uint_fast64_t j = 0; j < chain_length; j++)
    {
        fprintf(F, "%ld, %lf\n",j, state_array[0][1][chain_length-j-1].value);
    }
    fprintf(F, "\n");
    return 0;
}

uint_fast64_t
count_pruned(state_t *** state_array, size_t chain_length, closure_t closure){
    uint_fast64_t pr = 0, npr = 0;

    for(size_t i = 0; i < closure.close_at; i++){
        for(uint_fast64_t k = 0; k < POW2(i+1); k++){
            for(size_t j = 0; j < chain_length - i; j++){
                if(state_array[i][k][j].pruned == 0 && i < closure.close_at-1)
                    npr++;
                else
                    pr++;
            }
        }
    }

    printf("Pruned: %"PRIuFAST64", Not Pruned: %"PRIuFAST64", Ratio: %lf\n", pr, npr, (double)pr/(double)npr);
    return npr;
}

uint_fast64_t
count_pruned_unclosed(state_t *** state_array, size_t chain_length, closure_t closure){
    uint_fast64_t pr = 0, npr = 0;

    for(size_t i = 0; i < closure.close_at; i++){
        for(uint_fast64_t k = 0; k < POW2(i+1); k++){
            for(size_t j = 0; j < chain_length - i; j++){
                //also count point functions which will be approximated by closures
                if(state_array[i][k][j].pruned == 0 && i < closure.close_at)
                    npr++;
                else
                    pr++;
            }
        }
    }

    printf("Pruned: %"PRIuFAST64", Not Pruned: %"PRIuFAST64", Ratio: %lf\n", pr, npr, (double)pr/(double)npr);
    return npr;
}

double
check_consistency(state_t *** state_array, uint_fast32_t state_size, uint_fast64_t chain_length, closure_t closure)
{
    //check if inputs are sane
    if(state_size > chain_length)    
        err_exit("ERROR: Consistency check cannot be performed for states longer than chain length", -1);
    
    if(state_size >= closure.close_at)
        err_exit("ERROR: Consistency check cannot be performed for states longer than closure state length", -1);
    
    uint_fast32_t *unusable_positions;
    unusable_positions = (uint_fast32_t*)malloc((size_t)(chain_length-state_size)*sizeof(uint_fast32_t));

    for(uint_fast32_t i = 0; i < chain_length-state_size; i++){
        unusable_positions[i] = 0;
    }

    uint_fast8_t unusable = 0;

    //preliminary check if consistency can be calculated or if some of the values are pruned.
    printf("Unusable positions for consistency checks of point functions of order %ld are: ", state_size+1);
    for(size_t j = 0; j < chain_length-state_size; j++){
        for(uint_fast64_t k = 0; k < POW2(state_size+1); k++){           
            if(state_array[state_size][k][j].pruned == 1){
                if(unusable_positions[j] != 1){
                    unusable_positions[j] = 1;
                    printf("%ld, ", j);
                    unusable = 1;
                }
            }else{
                unusable_positions[j] = 0;
            }
        }
    }
    if(!unusable)
        printf("none");

    printf("\n");
    double state_sum = 0;

    for(size_t j = 0; j < chain_length-state_size; j++){
        if(!unusable_positions[j]){    
            for(uint_fast64_t k = 0; k < POW2(state_size+1); k++){
                state_sum+=state_array[state_size][k][j].value;
            }
            printf("Index low: %ld, summands: %ld, sum val: %lf, error: %1.10g\n", j, POW2(state_size+1), state_sum, 1.-state_sum);
        }
        state_sum = 0;
    }
    free(unusable_positions);
    return 0;
}

uint_fast8_t
closure_rhs(state_t *** state_array, closure_t closure, state_t *state, uint_fast32_t j, uint_fast64_t i_c, uint_fast64_t k_c, uint_fast64_t j_c)
{
    uint_fast32_t * index_u, * index_l;
    index_u = (uint_fast32_t *)malloc((size_t)(closure.u_cl_num) * sizeof(uint_fast32_t));
    index_l = (uint_fast32_t *)malloc((size_t)(closure.l_cl_num) * sizeof(uint_fast32_t));
    //------------------------------//
    
    if(index_l == NULL || index_u == NULL)
        err_exit("ERROR: Could not allocate index arrays for arbitrary closure->", 1);

    for(uint_fast32_t l = 0; l < closure.u_cl_num; l++){
        index_u[l] = 0;
        for(uint_fast32_t m = 0; m < closure.u_closures[l]; m++){
            index_u[l] += POW2(m) * state->bit_array[closure.u_cl_pos[l]+m];
        }
    }

    for(uint_fast32_t l = 0; l < closure.l_cl_num; l++){
        index_l[l] = 0;
        for(uint_fast32_t m = 0; m < closure.l_closures[l]; m++){
            index_l[l] += POW2(m) * state->bit_array[closure.l_cl_pos[l]+m];
        }
    }

    for(uint_fast32_t i = 0; i < closure.u_cl_num; i++){
        if(i_c == closure.u_closures[i]-1 && k_c == index_u[i] && j_c == closure.u_cl_pos[i]+j){
            free(index_l);
            free(index_u);
            return 1;
        }
    }

    for(uint_fast32_t i = 0; i < closure.l_cl_num; i++){
        if(i_c == closure.l_closures[i]-1 && k_c == index_l[i] && j_c == closure.l_cl_pos[i]+j){
            free(index_l);
            free(index_u);
            return 1;
        }
    }
    free(index_l);
    free(index_u);

    return 0;
}

uint_fast8_t
is_on_rhs(uint_fast64_t * lhs, uint_fast64_t * candidate, state_t *** state_array, closure_t closure, uint_fast64_t chain_length)
{

    uint_fast64_t i,i_c,k,k_c,j,j_c;
    i = lhs[0];
    k = lhs[1];
    j = lhs[2];

    i_c = candidate[0];
    k_c = candidate[1];
    j_c = candidate[2];

    //rhs is closure
    if(lhs[0]+2 == closure.close_at)
    {
        if((k >> i)&1){
            //state is result of jump in
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k-POW2(i) && j_c == j)
                    return 1;
            }else {
                if(closure_rhs(state_array, closure, &state_array[i+1][k+POW2(i)][j], j, i_c, k_c, j_c))
                    return 1;
            }
        } else { 
            //state is result of no(!) jump in
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k && j_c == j)
                    return 1;
            } else {
                if(closure_rhs(state_array, closure, &state_array[i+1][k+POW2(i+1)][j], j, i_c, k_c, j_c))
                    return 1;
            }
        }

        if(!(get_bit(0, &state_array[i][k][j]))){
            //state is result of jump out
            if(j > 0){
                if(closure_rhs(state_array, closure, &state_array[i+1][2*(k+1)][j-1], j-1, i_c, k_c, j_c))
                    return 1;
            }else{
                if(i_c == i && k_c == k+1 && j_c == j)
                    return 1;
            }
        } else { 
            //state is result of no(!) jump out
            if(j>0){ 
                if(closure_rhs(state_array, closure, &state_array[i+1][k*2][j-1], j-1, i_c, k_c, j_c))
                    return 1;
            } else {
                if(i_c == i && k_c == k && j_c == j)
                    return 1;
            }
        }

        //jumps inside the state
        for(uint_fast32_t l = 0; l < i; l++)
        {
            if( !(get_bit(l+1, &state_array[i][k][j])) && (get_bit(l, &state_array[i][k][j])) )
            {
                if(i_c == i && k_c == k+POW2(l) && j_c == j)
                    return 1;
            }

            if( get_bit(l+1, &state_array[i][k][j]) && !(get_bit(l, &state_array[i][k][j])) )
            {
                if(i_c == i && j_c == j && k_c == k)
                    return 1;
            }
        }
    }
    else
    {
        //run through steps of get_change to check if rhs candidate would be called by state encoded in lhs
        //if that is the case, we know that candidate appears on rhs and as such the derivative
        //will appear in the jacobian
        if((k >> i)&1){
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k-POW2(i) && j_c == j)
                    return 1;
            }else{
                if(i_c == i+1 && k_c == k+POW2(i) && j_c == j)
                    return 1;
            }
        } else {
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k && j_c == j)
                    return 1;
            } else {
                if(i_c == i+1 && k_c == k+POW2(i+1) && j_c == j)
                    return 1;
            }
        }

        if(!(get_bit(0, &state_array[i][k][j]))){
            if(j > 0){
                if(i_c == i+1 && k_c == 2*(k+1) && j_c == j-1)
                    return 1;
            } else {
                if(i_c == i && k_c == k+1 && j_c == j)
                    return 1;
            }
        } else {
            if(j>0){
                if(i_c == i+1 && k_c == 2*k && j_c == j-1)
                    return 1;
            } else {
                if(i_c == i && k_c == k && j_c == j)
                    return 1;
            }
        }

        for(uint_fast32_t l = 0; l < i; l++)
        {
            if( !(get_bit(l+1, &state_array[i][k][j])) && get_bit(l, &state_array[i][k][j]) )
            {
                if(i_c == i && k_c == k+POW2(l) && j_c == j)
                    return 1;
            }

            if( get_bit(l+1, &state_array[i][k][j]) && !(get_bit(l, &state_array[i][k][j])) )
            {
                if(i_c == i && k_c == k && j_c == j)
                    return 1;
            }
        }
    }
    return 0;
}

double 
closure_jacobian(state_t *** state_array, closure_t * closure, state_t * state, uint_fast64_t j, uint_fast64_t i_c, uint_fast64_t k_c, uint_fast64_t j_c)
{
    uint_fast32_t * index_u, * index_l;
    index_u = (uint_fast32_t *)malloc((size_t)(closure->u_cl_num) * sizeof(uint_fast32_t));
    index_l = (uint_fast32_t *)malloc((size_t)(closure->l_cl_num) * sizeof(uint_fast32_t));
    //------------------------------//
    
    if(index_l == NULL || index_u == NULL)
        err_exit("ERROR: Could not allocate index arrays for arbitrary closure->", 1);

    for(uint_fast32_t l = 0; l < closure->u_cl_num; l++){
        index_u[l] = 0;
        for(uint_fast32_t m = 0; m < closure->u_closures[l]; m++){
            index_u[l] += POW2(m) * state->bit_array[closure->u_cl_pos[l]+m];
        }
    }

    for(uint_fast32_t l = 0; l < closure->l_cl_num; l++){
        index_l[l] = 0;
        for(uint_fast32_t m = 0; m < closure->l_closures[l]; m++){
            index_l[l] += POW2(m) * state->bit_array[closure->l_cl_pos[l]+m];
        }
    }

    double change = 1.;

    for(uint_fast32_t i = 0; i < closure->u_cl_num; i++){
        if(state_array[closure->u_closures[i]-1][index_u[i]][closure->u_cl_pos[i]+j].value < 10e-20)
            change = 0;
        else if(!(i_c == closure->u_closures[i]-1 && k_c == index_u[i] && j_c == closure->u_cl_pos[i]+j))
            change = change * state_array[closure->u_closures[i]-1][index_u[i]][closure->u_cl_pos[i]+j].value; 
    }

    for(uint_fast32_t i = 0; i < closure->l_cl_num; i++){
        if(state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].value < 10e-20)
            change = 0;
        else if(i_c == closure->l_closures[i]-1 && k_c == index_l[i] && j_c == closure->l_cl_pos[i]+j){
            change = (-1.)*change / (state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].value * 
                                state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].value);
        } else
            change = change / state_array[closure->l_closures[i]-1][index_l[i]][closure->l_cl_pos[i]+j].value;
    }
    free(index_l);
    free(index_u);

    return change;    
}

double
rhs_jacobi_eval(uint_fast64_t * lhs_indices, uint_fast64_t * rhs_indices, state_t *** state_array, closure_t closure, uint_fast64_t chain_length, double * jr_array)
{
    double diff_val = 0;

    uint_fast64_t i,i_c,k,k_c,j,j_c;
    i = lhs_indices[0];
    k = lhs_indices[1];
    j = lhs_indices[2];

    i_c = rhs_indices[0];
    k_c = rhs_indices[1];
    j_c = rhs_indices[2];

    //rhs is closure
    if(i+2 == closure.close_at)
    {
        if((k >> i)&1){
            //state is result of jump in
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k-POW2(i) && j_c == j)
                    diff_val += jr_array[chain_length];
            }else {
                if(closure_rhs(state_array, closure, &state_array[i+1][k+POW2(i)][j], j, i_c, k_c, j_c))
                    diff_val += closure_jacobian(state_array, &closure, &state_array[i+1][k+POW2(i)][j], j, i_c, k_c, j_c);
            }
        } else { 
            //state is result of no(!) jump in
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k && j_c == j)
                    diff_val -= jr_array[chain_length];
            } else {
                if(closure_rhs(state_array, closure, &state_array[i+1][k+POW2(i+1)][j], j, i_c, k_c, j_c))
                    diff_val -= closure_jacobian(state_array, &closure, &state_array[i+1][k+POW2(i+1)][j], j, i_c, k_c, j_c);
            }
        }

        if(!(get_bit(0, &state_array[i][k][j]))){
            //state is result of jump out
            if(j > 0){
                if(closure_rhs(state_array, closure, &state_array[i+1][2*(k+1)][j-1], j-1, i_c, k_c, j_c))
                    diff_val += closure_jacobian(state_array,&closure, &state_array[i+1][2*(k+1)][j-1], j-1, i_c, k_c, j_c);
            }else{
                if(i_c == i && k_c == k+1 && j_c == j)
                    diff_val += jr_array[0];
            }
        } else { 
            //state is result of no(!) jump out
            if(j>0){ 
                if(closure_rhs(state_array, closure, &state_array[i+1][k*2][j-1], j-1, i_c, k_c, j_c))
                    diff_val -= closure_jacobian(state_array, &closure, &state_array[i+1][k*2][j-1], j-1, i_c, k_c, j_c);
            } else {
                if(i_c == i && k_c == k && j_c == j)
                    diff_val -= jr_array[0];
            }
        }

        //jumps inside the state
        for(uint_fast32_t l = 0; l < i; l++)
        {
            if( !(get_bit(l+1, &state_array[i][k][j])) && (get_bit(l, &state_array[i][k][j])) )
            {
                if(i_c == i && k_c == k+POW2(l) && j_c == j)
                    diff_val += jr_array[j+l+1];
            }

            if( get_bit(l+1, &state_array[i][k][j]) && !(get_bit(l, &state_array[i][k][j])) )
            {
                if(i_c == i && j_c == j && k_c == k)
                    diff_val -= jr_array[j+l+1];
            }
        }
    }
    else
    {
        //run through steps of get_change to check if rhs candidate would be called by state encoded in lhs
        //if that is the case, we know that candidate appears on rhs and as such the derivative
        //will appear in the jacobian
        if((k >> i)&1){
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k-POW2(i) && j_c == j)
                    diff_val += jr_array[chain_length];
            }else{
                if(i_c == i+1 && k_c == k+POW2(i) && j_c == j)
                    diff_val += jr_array[j+i+1];
            }
        } else {
            if(j+i == chain_length-1){
                if(i_c == i && k_c == k && j_c == j)
                    diff_val -= jr_array[chain_length];
            } else {
                if(i_c == i+1 && k_c == k+POW2(i+1) && j_c == j)
                    diff_val -= jr_array[j+i+1];
            }
        }

        if(!(get_bit(0, &state_array[i][k][j]))){
            if(j > 0){
                if(i_c == i+1 && k_c == 2*(k+1) && j_c == j-1)
                    diff_val += jr_array[j];
            } else {
                if(i_c == i && k_c == k+1 && j_c == j)
                    diff_val += jr_array[0];
            }
        } else {
            if(j>0){
                if(i_c == i+1 && k_c == 2*k && j_c == j-1)
                    diff_val -= jr_array[j];
            } else {
                if(i_c == i && k_c == k && j_c == j)
                    diff_val -= jr_array[0];
            }
        }

        for(uint_fast32_t l = 0; l < i; l++)
        {
            if( !(get_bit(l+1, &state_array[i][k][j])) && get_bit(l, &state_array[i][k][j]) )
            {
                if(i_c == i && k_c == k+POW2(l) && j_c == j)
                    diff_val += jr_array[j+l+1];
            }

            if( get_bit(l+1, &state_array[i][k][j]) && !(get_bit(l, &state_array[i][k][j])) )
            {
                if(i_c == i && k_c == k && j_c == j)
                    diff_val -= jr_array[j+l+1];
            }
        }
    }

    return diff_val;
}

void
save_jacobi_matrix(state_t *** state_array, uint_fast64_t chain_length, closure_t closure, double * jr_array)
{
    //first generate enumeration array
    uint_fast64_t num_vars = count_pruned(state_array, chain_length, closure);


    uint_fast64_t ** enum_state_array;
    enum_state_array = (uint_fast64_t**)malloc((size_t)num_vars*sizeof(uint_fast64_t*));
    if(enum_state_array == NULL){
        fprintf(stderr, "ERROR: Could not allocate array for jacobian matrix\n");
        return;
    }
    for(uint_fast64_t i = 0; i < num_vars; i++){
        enum_state_array[i]=(uint_fast64_t*)malloc((size_t)3*sizeof(uint_fast64_t));
        if(enum_state_array[i] == NULL){
            fprintf(stderr, "ERROR: Could not allocate array for jacobian matrix\n");
            return;
        }
    }


    uint_fast64_t num_vars_counted = 0;
    for(size_t i = 0; i < closure.close_at; i++){
        for(uint_fast64_t k = 0; k < POW2(i+1); k++){
            for(size_t j = 0; j < chain_length-i; j++){
                if(state_array[i][k][j].pruned == 0 && i < closure.close_at-1){
                    enum_state_array[num_vars_counted][0]=i;
                    enum_state_array[num_vars_counted][1]=k;
                    enum_state_array[num_vars_counted][2]=j;
                    num_vars_counted++;
                }
            }
        }
    }

    if(num_vars_counted != num_vars){
        fprintf(stderr, "ERROR: Missmatch in number of unpruned entries during generation of jacobian matrix: %ld vs %ld\n", num_vars, num_vars_counted);
        for(uint_fast64_t i = 0; i < num_vars; i++){
            free(enum_state_array[i]);
        }
        free(enum_state_array);
        return;
    }

    //calculate Jacobian

    //Allocate sparse arrays first
    num_vars = count_pruned(state_array, chain_length, closure);
    double * val;
    uint_fast64_t * col;    
    uint_fast64_t * row;

    val = (double*)malloc((size_t)num_vars*sizeof(double));
    col = (uint_fast64_t*)malloc((size_t)num_vars*sizeof(uint_fast64_t));
    row = (uint_fast64_t*)malloc((size_t)num_vars*sizeof(uint_fast64_t));

    if(val == NULL || col == NULL || row == NULL){
        fprintf(stderr, "ERROR: Cannot allocate arrays for jacobian matrix\n");
        for(uint_fast64_t i = 0; i < num_vars; i++){
            free(enum_state_array[i]);
        }
        free(enum_state_array);
        (val!=NULL) ? free(val):printf("sparse array for values couldn't be allocated");
        (row!=NULL) ? free(row):printf("sparse array for column couldn't be allocated");
        (col!=NULL) ? free(col):printf("sparse array for row couldn't be allocated");
    }


    uint_fast64_t sparse_counter = 0, arr_size = num_vars;
    for(uint_fast64_t i = 0; i < num_vars; i++){
        for(uint_fast64_t k = 0; k < num_vars; k++){
            if(is_on_rhs(enum_state_array[i], enum_state_array[k], state_array, closure, chain_length)){
                val[sparse_counter] = rhs_jacobi_eval(enum_state_array[i], enum_state_array[k], state_array, closure, chain_length, jr_array);
                row[sparse_counter] = i;
                col[sparse_counter] = k;
                sparse_counter++;
                if(sparse_counter >= arr_size){
                    val = reallocarray(val, 2*arr_size, sizeof(double));
                    col = reallocarray(col, 2*arr_size, sizeof(uint_fast64_t));
                    row = reallocarray(row, 2*arr_size, sizeof(uint_fast64_t));
                    arr_size *=2;
                }
            }
        }
    }

    FILE * F;
    F = fopen("testfile.dat", "w");

    for(uint_fast64_t i = 0; i < sparse_counter; i++){
        fprintf(F, "%ld %ld %lf\n", row[i]+1, col[i]+1, val[i]);
    }

    fclose(F);

    uint_fast32_t counter = 0;
    F = fopen("maple_testfile", "w");
    for(uint_fast32_t i = 0; i < num_vars; i++){
        for(uint_fast32_t k = 0; k < num_vars; k++){
            if (row[counter] == i && col[counter] == k)
            {
                fprintf(F, "%1.17g ", val[counter]);
                counter++;
            }
            else
                fprintf(F, "0 ");
        }
        fprintf(F, "\n");
    }

    for(uint_fast64_t i = 0; i < num_vars_counted; i++){
        free(enum_state_array[i]);
    }
    free(enum_state_array);
    free(row);
    free(val);
    free(col);
    return;
}

int is_equiv_rep(state_t state, int equiv_type)
{
    uint_fast8_t equiv_type_arr[2];
    equiv_type_arr[0] = equiv_type&1;
    equiv_type_arr[1] = (equiv_type>>1)&1;
    if(state.width>1){
        if(state.bit_array[0] == equiv_type_arr[0] && state.bit_array[state.width-1] == equiv_type_arr[1])
            return 1;
    }else{
        if(state.bit_array[0] == equiv_type_arr[0])
            return 1;
    }
    return 0;
}

double equivalence_class_val(closure_t closure, uint_fast8_t equiv_type, state_t state, state_t *** state_array)
{
    if(equiv_type > 3)
        err_exit("Equivalence type not allowed! Acceprable Classes are: 00, 01, 10, 11", -1);
    
       if(closure.close_at == state.width+1){
        //closure equiv class
    }
    else if(closure.close_at > state.width+1){
        if(is_equiv_rep(state, equiv_type))
            return state.value;

        uint_fast8_t eq_low = equiv_type&1;
        uint_fast8_t eq_high = (equiv_type>>1)&1;

        uint_fast32_t new_b = state.num;

        if(state.bit_array[0] != eq_low)
            if(eq_low == 0)
                new_b -= 1;
            if(eq_low == 1)
                new_b += 1;
        if(state.width > 1){
            if(state.bit_array[state.width-1] != eq_high)
                if(eq_high == 0)
                    new_b -= POW2(state.width-1);
                if(eq_high == 1);
                    new_b += POW2(state.width-1);
        }
        //return state_array[state.width-1][new_b][state.first].value;
        //DEBUG
        state_t temp = state_array[state.width-1][new_b][state.first];
        return temp.value;
    } else if(closure.close_at < state.width)
        return 0;
}

int main(int argc, char *argv[])
{   

    /* CMD OPTIONS*/
    int opt;
    int option_index = 0;

    int alpha_found = 0; 
    char * alpha_temp; 
    double alpha;

    int beta_found = 0; 
    char * beta_temp; 
    double beta;

    char * gamma_temp; 
    double gamma = 1;
    
    int chain_length_found = 0;
    char * n_temp;
    size_t chain_length = -1;

    int i_found = 0; 
    char * i_filename;
    int o_found = 0; 
    char * o_filename;

    int tra_found = 0;
    int tra_steps = 0;
    char * tra_temp;

    int solv_found = 0;
    char * solver; 

    int thread_num_found = 0; 
    uint_fast8_t thread_num = 1;

    int t_found = 0, T_found = 0; 
    char * t_temp, *T_temp;
    double t,T;

    int s_found = 0; 
    char * s_temp;
    uint_fast64_t steps; 
    
    uint_fast8_t prune = 0;

    char * depth_temp;
    uint_fast8_t depth_found = 0; 
    uint_fast32_t depth = 0;

    int verbose = 0;
	
	char * model_order_temp;
	int model_order;

    /*NUMERICAL STUFF*/
    const double LOW_D_BOUND = 1e-30;
    const double HIGH_D_BOUND = __DBL_MAX__;

    state_t *** state_array;

    double * jr_array;
    
    
    tv_set_inputs_t tv_inputs;

    tv_inputs.amplit_found = 0;
    tv_inputs.time_var_found = 0;    
    tv_inputs.funtyp_found = 0;
    tv_inputs.offset_found = 0;
    tv_inputs.period_found = 0;
    tv_inputs.constf_found = 0;

    tv_filenames_t tv_filenames;
    tv_arrays_t tv_arrays;
    
    /*MISC STUFF*/
    FILE * F;
    int err; //error state
 
    while ((opt = getopt_long(argc, argv, "a:b:c:g:n:i:o:t:T:s:m:p:d:vh?", long_opt, &option_index)) != -1){
        switch(opt){
        case 'a':
        case  0 :
            alpha_temp = strdup(optarg);
            alpha = atof(alpha_temp);
            free(alpha_temp);
            alpha_found = 1;
            break;
        case 'b':
        case  1 :
            beta_temp = strdup(optarg);
            beta = atof(beta_temp);
            free(beta_temp);
            beta_found = 1;
            break;
        case 'g':
        case  2 :
            gamma_temp = strdup(optarg);
            gamma = atof(gamma_temp);
            free(gamma_temp);
            break;
        case 'n':
        case  3 :
            n_temp = strdup(optarg);
            chain_length = atoi(n_temp);
            free(n_temp);
            chain_length_found = 1;
            break;
        case 'i':
        case  4 :
            i_filename = strdup(optarg);
            i_found = 1;
            break;
        case 'o':
        case  5 :
            o_filename = strdup(optarg);
            o_found = 1;
            break;
        case  6 :
            tra_found = 1;
            tra_temp = strdup(optarg);
            tra_steps = atoi(tra_temp);
            break;
        case 'm':
        case  7 :
            solv_found = 1;
            solver = strdup(optarg);
            break;
        case 'p':
        case  8 :
            thread_num_found = 1;
            char * temp;
            temp = strdup(optarg);
            thread_num = atoi(temp);
            free(temp);
            break;
        case 't':
        case  9 :
            t_temp = strdup(optarg);
            t = atof(t_temp);
            free(t_temp);
            t_found = 1;
            break;
        case 'T':
        case  10 :
            T_temp = strdup(optarg);
            T = atof(T_temp);
            free(T_temp);
            T_found = 1;
            break;
        case 's':
        case  11 :
            s_temp = strdup(optarg);
            steps = atoi(s_temp);
            free(s_temp);
            s_found = 1;
            break;
        case 'v':
        case  12 :
            verbose = 1;
            break;
        case 'h':
        case  13:
            usage(argv[0]);
            exit(0);
        case 14:
            prune = 1;
            break;
        case  15:
        case 'd':
            depth_temp = strdup(optarg);
            depth = atoi(depth_temp);
            free(depth_temp);
            depth_found = 1;
            break;
        case 16:
            tv_inputs.time_var_found = 1;
            break;
        case 17:
            tv_inputs.funtyp_found = 1;
            tv_filenames.funTyp_filename = strdup(optarg);
            break;
        case 18:
            tv_inputs.offset_found = 1;
            tv_filenames.offset_filename = strdup(optarg);
            break;
        case 19:
            tv_inputs.amplit_found = 1;
            tv_filenames.amplit_filename = strdup(optarg);
            break;
        case 20:
            tv_inputs.period_found = 1;
            tv_filenames.period_filename = strdup(optarg);
            break;
        case 21:
            tv_inputs.constf_found = 1;
            tv_filenames.constf_filename = strdup(optarg);
            break;
		case 'c':
			model_order_temp = strdup(optarg);
			model_order=atoi(model_order_temp);
			free(model_order_temp);
			break;
        default: 
            err_exit("ERROR: Unrecocgnized switch", -1);
            break;
        }

    }
    
    /* Check for correctness of inputs in roughly the same order as "long_opt[]" */

    if(alpha_found && beta_found){
        printf( "Found entry and exit rate\n"
                "\tentry rate:\t%lf\n"
                "\tmain rate:\t%lf\n"
                "\texit rate:\t%lf\n", alpha, gamma, beta);              
    } else if(i_found) {
        printf("Input file set to %s\n", i_filename);
        if(fopen(i_filename, "r") == NULL)
            err_exit("ERROR: Could not open input file", 1);
    } else if(!tv_inputs.time_var_found){
        err_exit("ERROR: Missing necessary inputs to configure hop rates", 1);
    }

    if(!chain_length_found)
        err_exit("ERROR: Missing necessary input for chain length (-n --num-sites)", 1);
    
    jr_array = (double *) malloc((size_t) (chain_length+1) * sizeof(double) );
    if(jr_array == NULL)
        err_exit("ERROR: Could not allocate memory for jump rate array", 1);

    if(i_found){
        F = fopen(i_filename, "r");
        for(unsigned i = 0; i < chain_length+1; i++){
            err = fscanf(F, "%lf", &jr_array[chain_length-i]);
            if(err == EOF)
                err_exit("ERROR: Reached early eof while reading from input file", 1);
        }
        fclose(F);
    } else {
        for(size_t i = 1; i < chain_length; i++){
            jr_array[i] = gamma;
        }
        jr_array[chain_length] = alpha;
        jr_array[0] = beta;
    }
    if(o_found){
        if(fopen(o_filename, "w") == NULL)
            err_exit("ERROR: Output file cannot be written to", 1);
        F = fopen(o_filename, "w");
    } else if(tra_found) {
        err_exit("ERROR: Cannot save temporal data without an output file", 1);
    } else {
        F = stdout;
    }

    if(tra_found){
        if(tra_steps < 0){
            printf("CAUTION: Negative number of steps in between snapshots. No snapshots will be saved\n");
            tra_steps = -1;
        } else {
            printf("Snapshots will be saved to file: %s every %d steps\n", o_filename, tra_steps);
        }
    } else {
        printf("No snapshots will be saved\n");
        tra_steps = -1;
    }

    if(solv_found){
        switch(solver[0]){
            case 'e':
                break;
            case 'r':
                break;
            default:
                err_exit("ERROR: Unknown solver type", 1);
                break;
        }
    } else {
        printf("Caution: No Solver specified. Default is Runge-Kutta 4\n");
        solver = "r";
    }

    if(tv_inputs.time_var_found){
        err = prepare_time_var_arrays(&tv_arrays, &tv_filenames, &tv_inputs, chain_length);
        if(err)
            err_exit("ERROR: Time var input files not read correctly.", err);    
    }

    #ifdef HAVE_OMP
        if(thread_num_found){
            if(thread_num > omp_get_max_threads()){
                printf("Caution: value specified by [-p|--parallel-threads] exceeds avaliable number of threads on system. Value will be decreased to match\n");
                thread_num = omp_get_max_threads();
            }

            printf("Threads used during calculation: %d/%d, OMP return: %d\n", thread_num, omp_get_max_threads(), omp_get_num_threads());       
        }
        omp_set_num_threads(thread_num);
    #endif

    if(!t_found || !T_found || !s_found)
        err_exit("ERROR: Missing necessary input for either times or step number", 1);

    if(t > HIGH_D_BOUND || t < 0)
        err_exit("ERROR: t does not have a reasonable value", 1);
    
    if(T > HIGH_D_BOUND || T < LOW_D_BOUND)
        err_exit("ERROR: T does not have a reasonable value", 1);
    
    if(T < t)
        err_exit("ERROR: T must be larger than t", 1);
    
    if(steps < 0 )
        err_exit("ERROR: Step number too low", 1);

    if(verbose){
        printf("Verbose output\n");
    }

    struct timeval start, stop;

    FILE * logfile = stdout;

    state_t  pruned_state_standin;

    if(prune){

        pruned_state_standin.bit_array = (uint_fast8_t *) malloc((size_t)POW2(chain_length)*sizeof(uint_fast8_t));
        if(pruned_state_standin.bit_array == NULL){
            err_exit("ERROR: Could not allocate memory for pruned state standin", 1);
        }
        pruned_state_standin.width = chain_length,
        pruned_state_standin.num = 0,
        pruned_state_standin.first = 0,
        pruned_state_standin.value = 0,
        pruned_state_standin.v_delta = 0,
        pruned_state_standin.pruned = 0;

        gen_bit_array(&pruned_state_standin);


        pruned_state_standin.pruned = 1; 
    }
    closure_t closure;

    closure.u_cl_num = 2;
    if(model_order > 1){
		closure.l_cl_num = 1;
	}else{
		closure.l_cl_num = 0;
	}
	fprintf(logfile, "Model Order is: %d\n", model_order);
	closure.close_at =      model_order+1;
    closure.u_closures[0] = model_order;
    closure.u_closures[1] = model_order;
    closure.l_closures[0] = model_order-1;

    closure.u_cl_pos[0] = 0;
    closure.u_cl_pos[1] = 1;
    
    closure.l_cl_pos[0] = 1;


    check_closure(&closure);

    depth_found = 1;
    depth = closure.close_at;


    if(verbose){
        
        //double probs[10] = {1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
        //double probs[10] = {0.4,0.6,0.4,0.6,0.4,0.6,0.4,0.6,0.4,0.6};
    
        fprintf(logfile, "Allocating state array\n");
        gettimeofday(&start, NULL);
        err = malloc_state_array_verbose(logfile, &state_array, chain_length, closure);
        if(err)
            err_exit("ERROR: Could not allocate state array", err);
        
        //err = init_state_array_custom(state_array, chain_length, closure.close_at + 1, probs);
        err = init_state_array_verbose(logfile, state_array, chain_length, closure);
        if(err)
            err_exit("ERROR: Could not initialize state array", err);
        
        
        if(depth_found){
            prune_array(state_array, jr_array, chain_length, &closure);
            count_pruned(state_array, chain_length, closure);
            free_pruned_entries(state_array, chain_length, closure);
        }
        
        gettimeofday(&stop, NULL);
        double secs = (double)(stop.tv_usec - start.tv_usec)/1000000 + (double)(stop.tv_sec - start.tv_sec);
        fprintf(logfile, "Finished after %lf seconds\n", secs);

        fprintf(logfile, "Starting simulation\n");

        //DEBUG TESTS
        for(uint_fast32_t i = 0; i < closure.close_at; i++){
            for(uint_fast64_t k = 0; k < POW2(i+1);k++){
                for(uint_fast32_t j= 0; j < chain_length-i; j++){
                    if(state_array[i][k][j].pruned == 0)           
                        equivalence_class_val( closure, (uint_fast8_t)3, state_array[i][k][j], state_array);
                }
            }
        }
        //END TESTS
        
        switch(solver[0]){
            default:
                err_exit("Unknown solver", 1);
                break;
            case 'e':
                gettimeofday(&start, NULL);
                err = solve_euler_verbose(logfile, F, tra_steps, state_array, jr_array, chain_length, t, T, steps, &closure, tv_inputs.time_var_found, &tv_arrays);
                if(err)
                    err_exit("ERROR: Issue during simulation", err);
                
                gettimeofday(&stop, NULL);
                break;
            case 'r':
                gettimeofday(&start, NULL);
                err = solve_rk4_verbose(logfile, F, tra_steps, state_array, jr_array, chain_length, t, T, steps, &closure, tv_inputs.time_var_found, &tv_arrays);
                if(err)
                    err_exit("ERROR: Issue during simulation", err);
                
                gettimeofday(&stop, NULL);
                break;
            }

        fprintf(logfile, "Finished simulation\n");
        secs = (double)(stop.tv_usec - start.tv_usec)/1000000 + (double)(stop.tv_sec - start.tv_sec);
        fprintf(logfile, "Time elapsed: %lf seconds\n", secs);
        //save_state_pruned(logfile, state_array, chain_length, T);

        save_state_singles(logfile, state_array, chain_length, T);
        save_state_singles(F, state_array, chain_length, T);

        check_consistency(state_array, closure.close_at-2, chain_length, closure);

    } else {
        err = malloc_state_array(&state_array, chain_length);
        if(err)
            err_exit("ERROR: Could not allocate state array", err);
        err = init_state_array(state_array, chain_length, closure);
        if(err)
            err_exit("ERROR: Could not initialize state array", err);
        
        if(depth_found){
            prune_array(state_array, jr_array, chain_length, &closure);
            count_pruned(state_array, chain_length, closure);
            free_pruned_entries(state_array, chain_length, closure);
        } else {
            pseudo_prune_array(state_array, chain_length, closure.close_at);
        }

        switch(solver[0]){
            case 'e':
                err = solve_euler(F, tra_steps, state_array, jr_array, chain_length, t, T, steps, &closure, tv_inputs.time_var_found, &tv_arrays);
                if(err)
                    err_exit("ERROR: Issue during simulation", err);
                break;
            
            case 'r':
                err = solve_rk4(F, tra_steps, state_array, jr_array, chain_length, t, T, steps, &closure, tv_inputs.time_var_found, &tv_arrays);
                if(err)
                    err_exit("ERROR: Issue during simulation", err);
                break;
            default:
                err_exit("Unknown solver", 1);
        }

        save_state_singles(F, state_array, chain_length, T);
    }
    
    save_jacobi_matrix(state_array, chain_length, closure, jr_array);    



    if(prune)
        free(pruned_state_standin.bit_array);


    if(i_found)
        free(i_filename);
    
    if(o_found){
        free(o_filename);
        fclose(F);
    }

    if(tra_found)   
        free(tra_temp);
    
    if(solv_found)
        free(solver);

    free(jr_array);
    
    free_state_array_full(state_array, chain_length, closure.close_at);

    return 0;

}
