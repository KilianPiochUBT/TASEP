#include "masterEquation.h"

#define uninitialized_var(x) x = x


static inline uint64_t log2_u64(uint64_t n)
{
  int result = 0;
  while (n >>= 1) ++result;
  return result;
}

static void usage(const char *progname)
{
  printf("Usage: %s\n-i jump rate file\n-n chain length\n-o output file name\n-t start time\n-T endtime\n"
  "-s steps\n-f fast option (caution: this requires more memory)\n"
  "-a space between snapshots\n"
  "--time-var sets jump rates to time variant. CAUTION: This cannot be used in combination with -f\n"
  "--function-type <char* filename with sin,cos,lin in order> sets function type of jump rate (when using <lin> --period option is not used)\n"
  "--amplitude <char* filename> reads amplitudes from given file and sets functions, when using <lin> these serve as the gradient\n"
  "--period <char* filename> reads period multipliers from given file an sets jump rate functions. (default is 1)\n"
  "--offset <char* filename> reads offsets from given file and sets jump rate functions. (default is 0)\n"
  "--constant <char* filename> reads constants from given file and sets jump rate functions. (default is 0)\n"
  "--thread-num <char* number_of_threads_to_use> sets number of threads used during parallel portion of the code. Only applies if -f is also selected. (default is max Threads)\n"
  "--save-select <char* shorts for what to save> shorts are: 'a': all 'e': expected values 'p': probabilities 'c': covariance Matrices 'm': multiplicative expected values\n"
  "--alpha <float> entry flow rate\n"
  "--beta <float> exit flow rate\n"
  "--gamma <float> flow rate inside lattice\n"
  "-h\n", progname);
}

double
key_fit(double x, uint_fast32_t id, uint_fast64_t key, uint_fast8_t offset, uint_fast32_t size)
{
    uint_fast8_t add_bool = 1;
    id = id >> offset;
    for(uint_fast32_t i = 0; i < size; i++){
        if(((key >> i)&1) != ((id >> i)&1))
            add_bool = 0;
    }
    if(add_bool)
        return x;
    else 
        return 0.;
}

double
high_order_corr(double * x, uint_fast64_t dim, uint_fast64_t key, uint_fast8_t key_size, uint_fast8_t offset, uint_fast32_t size)
{
    if(key_size+offset > size)
        exit(-1);
    
    if(offset == 0)
        exit(-1);
    
    uint_fast64_t key_left_ex, key_right_ex, key_lr_ex;
    key_left_ex = key+POW2(key_size);
    key_right_ex = (key << 1)+1;
    key_lr_ex = (key<<1)+1+POW2(key_size+1);

    double lr_ex = 0, no_ex = 0, left_ex = 0, right_ex = 0;

    for(uint_fast32_t i = 0; i < dim; i++){
        lr_ex += key_fit(x[i], i, key_lr_ex, offset-1, key_size+2);
        no_ex += key_fit(x[i], i, key, offset, key_size);
        left_ex += key_fit(x[i], i, key_left_ex, offset, key_size+1);
        right_ex += key_fit(x[i], i, key_right_ex, offset-1, key_size+1);
    }

    double ret = 0;

    ret = lr_ex*no_ex-left_ex*right_ex;
    ret = ret / sqrt(left_ex*(no_ex-left_ex)*right_ex*(no_ex-right_ex));

}

double
error_guess(double * x, uint_fast64_t dim, uint_fast64_t key, uint_fast8_t key_size, uint_fast8_t offset, uint_fast32_t size)
{
    if(key_size + offset > size || offset == 0)
        exit(-1);
    uint_fast64_t key_left_ex, key_right_ex, key_lr_ex;
    key_left_ex = key+POW2(key_size);
    key_right_ex = (key << 1)+1;
    key_lr_ex = (key<<1)+1+POW2(key_size+1);

    double lr_ex = 0, no_ex = 0, left_ex = 0, right_ex = 0;

    for(uint_fast32_t i = 0; i < dim; i++){
        lr_ex += key_fit(x[i], i, key_lr_ex, offset-1, key_size+2);
        no_ex += key_fit(x[i], i, key, offset, key_size);
        left_ex += key_fit(x[i], i, key_left_ex, offset, key_size+1);
        right_ex += key_fit(x[i], i, key_right_ex, offset-1, key_size+1);
    }

    double ret = 0;

    ret = (lr_ex/no_ex - (left_ex/no_ex)*(right_ex/no_ex))*no_ex;

    return ret;
}

double
higher_order_pf_approx(double *x, uint_fast64_t dim, uint_fast64_t key, uint_fast8_t key_size, uint_fast8_t offset, uint_fast8_t size)
{
    if(key_size + offset > size)
        exit(-1);

    uint_fast64_t key_left_ex, key_right_ex;
    key_left_ex = key+POW2(key_size);
    key_right_ex = (key << 1)+1;

    double lr_ex = 0, no_ex = 0, left_ex = 0, right_ex = 0;

    for(uint_fast32_t i = 0; i < dim; i++){
        no_ex += key_fit(x[i], i, key, offset, key_size);
        left_ex += key_fit(x[i], i, key_left_ex, offset, key_size+1);
        right_ex += key_fit(x[i], i, key_right_ex, offset-1, key_size+1);
    }

    double ret = 0;

    ret = (left_ex*right_ex)/no_ex;

    return ret;

}

double
higher_order_point_func(double *x, uint_fast64_t dim, uint_fast64_t key, uint_fast8_t key_size, uint_fast8_t offset, uint_fast8_t size)
{
    if(key_size + offset > size)
        exit(-1);

    double pf_val = 0;

    for(uint_fast64_t i = 0; i < dim; i++)
    {
        pf_val += key_fit(x[i], i, key, offset, key_size);
    }

    return pf_val;
}

int main(int argc, char *argv[])
{

    int opt;
    int option_index = 0;
 
    int if_found = 0;
    char* i_filename;
    
    int of_found = 0;
    char* o_filename;
    
    int t_found = 0;
    double t = 0;
    char * t_temp;

    int T_found = 0;
    double T = 0;
    char *T_temp;

    int n_found = 0;
    int n;
    char* n_temp;

    int s_found = 0;
    unsigned int steps;
    char * s_temp;

    int f_found = 0;

    int a_found = 0;
    char* a_temp;
    int animation_resolution = 0;
    double * matrix_data;
    double * multExp_matrix;
    double * exp_values;
    double * prob_values;
    uint64_t animation_data_pointer = 0, animated_step_check;
    unsigned int prob_value_data_pointer = 0;


    double checksum, h, stop_condition = 1e-7;
    double * k4, * k3, * k2, * k1, * x, * x1;
    uint64_t steps_at_stop = 0;
    int early_stop = 0;

    uint64_t xk = 0, xi = 0, temp;
    double cov = 0;

    int     time_var_found = 0;
    int     funtyp_found = 0;
    int     offset_found = 0;
    int     period_found = 0;
    int     amplit_found = 0;
    int     constf_found = 0;
    char*   funTyp_filename;
    char*   offset_filename;
    char*   period_filename;
    char*   amplit_filename;
    char*   constf_filename;
    int*    funTyp_array; 
    double* offset_array;
    double* period_array;
    double* amplit_array;
    double* constf_array;

    double alpha, beta, gamma;
    int alpha_found = 0, beta_found = 0, gamma_found = 0;
    char * alpha_temp, * beta_temp, * gamma_temp;    
    
    int save_select_found = 0;
    char * save_selection;
    int num_threads = 1;
    
    #ifdef HAVE_OMP 
        num_threads = omp_get_max_threads();
    #endif
    
    char* num_thread_temp;

    FILE* F;
    int err;
    uint64_t dim;
    double startVecValue = 0.;
    double * jump_rate_Array;
    uint64_t * colInd;
    double * valArr;
    int*bit;

    double secs;

    double * expValueAnimation;
    uint64_t expValuePtr = 0;
    
    while ((opt = getopt_long(argc, argv, "i:v:o:t:T:n:s:a:fh?", long_opt, &option_index)) != -1){
        switch(opt){
        case 'i':
            if_found = 1;
            i_filename = strdup(optarg);
            break;
        case 'o':
            of_found = 1;
            o_filename = strdup(optarg);
            break;
        case 't':
            t_found = 1;
            t_temp = strdup(optarg);
            t = atof(t_temp);
            free(t_temp);
            break;
        case 'T':
            T_found = 1;
            T_temp = strdup(optarg);
            T = atof(T_temp);
            free(T_temp);
            break;
        case 'n':
            n_found = 1;
            n_temp = strdup(optarg);
            n = atoi(n_temp);
            free(n_temp);
            break;
        case 's':
            s_found = 1;
            s_temp = strdup(optarg);
            steps = atoi(s_temp);
            free(s_temp);
            break;
        case 'f':
            f_found = 1;
            break;
        case 'a':
            a_found = 1;
            a_temp = strdup(optarg);
            animation_resolution = atoi(a_temp);
            free(a_temp);
            break;
        case 0:
            time_var_found = 1;
            break;
        case 1:
            funtyp_found = 1;
            funTyp_filename = strdup(optarg);
            break;
        case 2:
            offset_found = 1;
            offset_filename = strdup(optarg);
            break;
        case 3:
            amplit_found = 1;
            amplit_filename = strdup(optarg);
            break;
        case 4:
            period_found = 1;
            period_filename = strdup(optarg);
            break;
        case 5:
            constf_found = 1;
            constf_filename = strdup(optarg);
            break;
        case 6:
            save_select_found = 1; 
            save_selection = strdup(optarg);
            break;
        case 7:
            num_thread_temp = strdup(optarg);
            num_threads = atoi(num_thread_temp);
            free(num_thread_temp);
            break;
        case 8:
            alpha_temp = strdup(optarg);
            alpha = atof(alpha_temp);
            alpha_found = 1;
            free(alpha_temp);
            break;
        case 9:
            beta_temp = strdup(optarg);
            beta = atof(beta_temp);
            beta_found = 1;
            free(beta_temp);
            break;
        case 10:
            gamma_temp = strdup(optarg);
            gamma = atof(gamma_temp);
            gamma_found = 1;
            free(alpha_temp);
            break;
        default:
            usage(argv[0]);
            exit(-1); 
        }
    }

    if(f_found && time_var_found){
        err_exit("ERROR: -f cannot be selected in combination with --time-var", -1);
    }

    if(!n_found){
        printf("Enter length: ");
        err = fscanf(stdin, "%d", &n);
        if(err == EOF)
            err_exit("ERROR: Issue while reading n from stdin", err);
    }
    dim = POW2(n);
    x = (double *) malloc(((size_t)dim) * sizeof(double));
    if(x == NULL)
    {
        err_exit("ERROR: Not enough memory to allocate x array", -1);
    }

    printf("Allocated result-vector\n");

    if(alpha_found || beta_found || gamma_found){
        if(!alpha_found){
            fprintf(stderr, "ERROR: No value for alpha found. Please ammend option --alpha to command line options. Exiting\n");
            exit(-1);
        }
        else{
            fprintf(stdout,"Found value for\n   alpha: %lf\n",alpha);
        }
        if(!gamma_found){
            fprintf(stderr, "CAUTION: No value for gamma found. If you want to specify gamma, please ammend option --gamma to command line options.\n");
            gamma = 1;
            gamma_found = 1;
            fprintf(stdout,"   gamma: %lf\n",gamma);
        }
        else{
            fprintf(stdout,"   gamma: %lf\n",gamma);
        }
        if(!beta_found){
            fprintf(stderr, "ERROR: No value for beta found. Please ammend option --beta to command line options. Exiting\n");
            exit(-1);
        }
        else{
            fprintf(stdout,"   beta: %lf\n",beta);
        }
    }
    else{
        if(!if_found && !time_var_found){
            printf("Enter %d probabilities: ", n+1);
            jump_rate_Array = enterProbs(n+1);
        }
    }
    if(!T_found){
        printf("No end time found. Please enter: ");
        err = scanf("%lf", &T);
        if(err == EOF)
            err_exit("ERROR: Issue while reading end time from stdin", err);
    }

    if(!t_found){        
        printf("No start time found. Please enter: ");
        err = scanf("%lf", &t);
        if(err == EOF)
            err_exit("ERROR: Issue while reading start time from stdin", err);
    }

    double d = POW2(n);
    startVecValue = 1./d;
    printf("Generated inital vector with identical values(%.30g)\n", startVecValue);
    

    //////CHANGED
    for(unsigned int i = 0; i < dim; i++)
        x[i] = startVecValue;

    //x[3]=1;

    if(!of_found){
        o_filename = (char*)malloc(19*sizeof(char));
        of_found = 1;
        printf("WARNING: No output file found. Generating default output file ('generic_output.vec').\n");
        snprintf(o_filename, 19, "generic_output.vec");
        if(fopen(o_filename, "w") == NULL){
            printf("ERROR: Cannot open generic output file. Check permissions.\n");
            exit(-1);
        }
    }

    if(!s_found){
        printf("No step number found. Please Enter: ");
        err = scanf("%u", &steps);
        if(err == EOF)
            err_exit("ERROR: Issue while reading step number from stdin", err);
    }   


    if(if_found){
        if(fopen(i_filename, "r") != NULL){
            F = fopen(i_filename, "r");
            int ctr = 0;
            jump_rate_Array = (double*)malloc(((size_t)n + 1) * sizeof(double));
            if(jump_rate_Array == NULL){
                printf("ERROR: Not enough memory to allocate probability array. ABORTING\n");
                exit(-1);
            }
            for(ctr = 0; ctr < n+1; ctr++){
                err = fscanf(F, "%lf", &jump_rate_Array[ctr]);
                if(err == EOF)  
                    err_exit("ERROR: Issue while reading from jump rate file name", err);
            }
        }
        fclose(F);
        if(alpha_found || beta_found || gamma_found)
            fprintf(stdout,"CAUTION: Switch '-i' overrides values given for alpha, beta and gamma.\n");
    }
    else if(alpha_found && beta_found){
        jump_rate_Array = (double*)malloc(((size_t)n + 1) * sizeof(double));
        if(jump_rate_Array == NULL){
            printf("ERROR: Not enough memory to allocate probability array. ABORTING\n");
            exit(-1);
        }
        jump_rate_Array[0] = alpha;
        jump_rate_Array[n] = beta;
        for(int_fast16_t i = 1; i < n; i++)
            jump_rate_Array[i] = gamma;
    }


    bit = (int *) malloc((n) * sizeof(int));
    if(bit == NULL)
    {
        fprintf(stderr, "ERROR: Not enough memory to allocate bit array.\n");
        exit(-1);
    }


    x1 = (double *) malloc(((size_t)dim) * sizeof(double));
    if(x1 == NULL)
    {
        fprintf(stderr, "ERROR: Not enough memory to allocate x1 array.\n");
        exit(-1);
    }

    k1 = (double *) malloc(((size_t)dim) * sizeof(double));
    if(k1 == NULL)
    {
        fprintf(stderr, "ERROR: Not enough memory to allocate k1 array.\n");
        exit(-1);
    }

    k2 = (double *) malloc(((size_t)dim) * sizeof(double));
    if(k2 == NULL)
    {
        fprintf(stderr, "ERROR: Not enough memory to allocate k2 array.\n");
        exit(-1);
    }

    k3 = (double *) malloc(((size_t)dim) * sizeof(double));
    if(k3 == NULL)
    {
        fprintf(stderr, "ERROR: Not enough memory to allocate k3 array.\n");
        exit(-1);
    }

    k4 = (double *) malloc(((size_t)dim) * sizeof(double));
    if(k4 == NULL)
    {
        fprintf(stderr, "ERROR: Not enough memory to allocate k4 array.\n");
        exit(-1);
    }

    printf("Allocated calculation vectors\n");

    //Creating a new directory to save the pictures and matrices during calculation
    if(a_found){
        printf("Allocatng memory for animation data\n");
        matrix_data = (double *)malloc(((size_t)(n*n*steps))* sizeof(double));
        if(matrix_data == NULL){
            printf("ERROR: Not Enough memory to allocate animation data array\n");
            exit(-1);
        }
        multExp_matrix = (double *)malloc(((size_t)(n*n*steps))* sizeof(double));
        if(multExp_matrix == NULL){
            printf("ERROR: Not Enough memory to allocate animation data array\n");
            exit(-1);
        }
        exp_values = (double *)malloc(((size_t)n) * sizeof(double));
        if(exp_values == NULL){
            printf("ERROR: Not enough memory to allocate expected values Array\n");
            return(EXIT_FAILURE);
        }

        expValueAnimation = (double *)malloc(((size_t)n * steps) * sizeof(double));
        if(expValueAnimation == NULL){
            printf("Not enough memory to allocate animation data array for expected Values\n");
            exit(-1);
        }

        prob_values = (double *)malloc(((size_t) dim * (steps / animation_resolution)) * sizeof(double));
        if(prob_values == NULL){
            printf("ERROR: Not enough memory to allocate probability values Array\n");
            exit(-1);
        }
    }
    
    #ifdef HAVE_OMP
        if(num_threads <= omp_get_max_threads())
            omp_set_num_threads(num_threads);
        else{
            omp_set_num_threads(1);
        }
    #endif
    
    printf("Number of threads to be used: %d\n", num_threads);

    int save_expected_value_data = 0, save_probability_data = 0, save_covariance_data = 0, save_multiplicative_expVal_data = 0;
    if(save_select_found){
        int savetypes = strlen(save_selection);
        for(int i = 0; i < savetypes; i++){
            switch(save_selection[i]){
                case 'a':
                    save_expected_value_data = 1;
                    save_probability_data = 1;
                    save_covariance_data = 1;
                    break;
                case 'e':
                    save_expected_value_data = 1;
                    break;
                case 'p':
                    save_probability_data = 1;
                    break;
                case 'c':
                    save_covariance_data = 1;
                    break;
                case 'm':
                    save_multiplicative_expVal_data = 1;
                    break;
                default:
                    printf("Could not parse save selction\n");
            }
        }
    }
    else{
        save_expected_value_data = 1;
        save_probability_data = 1;
        save_covariance_data = 1;
        save_multiplicative_expVal_data = 1;
    }


    h = (T - t)/((long double) steps);
    
 
    gettimeofday(&start, NULL);
    double prozent;
    if(f_found){
        printf("Starting Calculation!\n");
        uint64_t size = 2 * POW2(n) + POW2(n-2)*(n-1); //exact number of Values
        uint64_t * rowPtr;
        rowPtr = (uint64_t *) malloc((dim + 1) * sizeof (long unsigned));
        if(rowPtr == NULL)
        {
            fprintf(stderr, "ERROR: Not enough memory to allocate rowPtr.\nTry rerunning program without 'f' flag enabled");
            exit(-1);
        }
        
        valArr = (double *) malloc((size) * sizeof(double));
        if(valArr == NULL)
        {
            fprintf(stderr, "ERROR: Not enough memory to allocate valArr.\nTry rerunning program without 'f' flag enabled");
            exit(-1);
        }

        colInd = (uint64_t *) malloc((size) * sizeof(long unsigned));
        if(colInd == NULL)
        {
            fprintf(stderr, "ERROR: Not enough memory to allocate colInd.\nTry rerunning program without 'f' flag enabled");
            exit(-1);
        }

        printf("Generating Matrix:\n");
        
        gettimeofday(&perf_start, NULL);
        
		gencrs(rowPtr, valArr, colInd, jump_rate_Array, bit, dim, n);//generate matrix in memory
        free(bit);

        gettimeofday(&perf_stop, NULL);
        secs = (double)(perf_stop.tv_usec - perf_start.tv_usec)/1000000 + (double)(perf_stop.tv_sec - perf_start.tv_sec);

        printf("Generating Matrix: 100.000000%%\nTime elapsed: %f seconds\nCalculating result vector..\n\n", secs);

        for(uint64_t i = 0; i < steps; i++){   
            checksum = 0;    
            rk4_step(1, dim, k1, x, x1, rowPtr, colInd, valArr, h);
			
			rk4_step(2, dim, k2, x, x1, rowPtr, colInd, valArr, h);
			
			rk4_step(3, dim, k3, x, x1, rowPtr, colInd, valArr, h);
			
			rk4_step(4, dim, k4, x, x1, rowPtr, colInd, valArr, h);
            

            for(unsigned int k = 0; k < dim; k++){
                x[k] += (k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.;
                checksum += fabs((k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.);
            }
    
            t += h;

            if(a_found){
                    animated_step_check = i % animation_resolution;
                    if(animated_step_check == 0){
                        //Calculate exprected values of x
                        
						exp_val_eval(dim, exp_values, x);
                        
						for(int ctr = 0; ctr < n; ctr++){
                            expValueAnimation[(n*expValuePtr) + ctr] = exp_values[ctr];
                        }

						expValuePtr++;

						cov_matrix_eval(n, dim, exp_values, x, matrix_data, multExp_matrix, &animation_data_pointer);

                        for(uint64_t ctr = 0; ctr < dim; ctr++){
                            prob_values[prob_value_data_pointer * dim + ctr] = x[ctr];
                        }
                        prob_value_data_pointer++;
                    }
                }
            
            printf("Sim Time: %lf ------ Delta %.5g ------ Stop condition %d%%\r", t, checksum, early_stop);
            
            if(checksum < stop_condition*h){
                early_stop++;
                steps_at_stop = i; 
            } 
            else{
                early_stop = 0;
            }
            if(early_stop > 100){
                break;
            }


        }
        free(rowPtr);
    }
    else{

        colInd = (uint64_t *) malloc((n) * sizeof(long unsigned));
        if(colInd == NULL)
        {
            fprintf(stderr, "ERROR: Not enough memory to allocate col.\n");
            exit(-1);
        }

        valArr = (double *) malloc((n) * sizeof(double));
        if(valArr == NULL)
        {
            fprintf(stderr, "ERROR: Not enough memory to allocate val.\n");
            exit(-1);
        }

        if(time_var_found){
            //read function type
            if(funtyp_found){
                if(fopen(funTyp_filename, "r") == NULL){
                    printf("ERROR: Cannot open file containing function Types. Aborting\n");
                    exit(-1);
                }
                funTyp_array = (int *)malloc(((size_t)n+1) * sizeof(int));
                if(funTyp_array == NULL){
                    printf("ERROR: Not enough space to allocate function type array. Aborting\n");
                    exit(-1);
                }
                F = fopen(funTyp_filename, "r");
                char funTypSelector[256];
                for(int i = 0; i < n+1; i++){
                    err = fscanf(F, "%255s", funTypSelector);
                    if(err == EOF)
                        err_exit("ERROR: Issue while reading function type file: Reached early EOF", err);
                    switch(funTypSelector[0]){
                        case 's':
                            funTyp_array[i] = 0;
                            break;
                        case 'c':
                            funTyp_array[i] = 1;
                            break;
                        case 'l':
                            funTyp_array[i] = 2;
                            break;
                        default:
                            printf("ERROR: Non identifiable function type fonud in"
                                   "%s. Read %s. Aborting!\n", funTyp_filename, funTypSelector);
                    }
                }
                printf("Read function types successfully\n");
            }
            else{
                funTyp_array = (int *)malloc(((size_t)n+1) * sizeof(int));
                if(funTyp_array == NULL){
                    printf("ERROR: Not enough space to allocate function type array. Aborting\n");
                    exit(-1);
                }
                for(int i = 0; i < n+1; i++){
                    funTyp_array[i] = 0;
                }
                
                printf("Wrote function type array successfully\n");
            }

            if(offset_found){
                //read offset values
                if(fopen(offset_filename, "r") == NULL){
                    printf("ERROR: Cannot open file containing offset values. Aborting\n");
                    exit(-1);
                }
                offset_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(offset_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                F = fopen(offset_filename, "r");
                for(int i = 0; i < n+1; i++){
                    err = fscanf(F, "%lf, ", &offset_array[i]);
                    if(err == EOF)
                        err_exit("ERROR: Issue while reading offset filename", err);
                }

                printf("Read offset values successfully\n");
            }
            else{
                offset_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(offset_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                for(int i = 0; i < n+1; i++)
                    offset_array[i] = 0.;

                printf("Wrote offset values sucessfully\n");
            }

            if(amplit_found){
                //read amplitude values
                if(fopen(amplit_filename, "r") == NULL){
                    printf("ERROR: Cannot open file containing amplitude values. Aborting\n");
                    exit(-1);
                }
                amplit_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(amplit_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                F = fopen(amplit_filename, "r");
                for(int i = 0; i < n+1; i++){
                    err = fscanf(F, "%lf, ", &amplit_array[i]);
                    if(err == EOF)
                        err_exit("ERROR: Issue while reading amplitude-file", err);
                }

                printf("Read amplitude values correctly\n");
            }
            else{
                amplit_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(amplit_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                for(int i = 0; i < n+1; i++)
                    amplit_array[i] = 1.;

                printf("Wrote amplitude values correctly\n");
            }

            if(period_found){
                //read amplitude values
                if(fopen(period_filename, "r") == NULL){
                    printf("ERROR: Cannot open file containing amplitude values. Aborting\n");
                    exit(-1);
                }
                period_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(period_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                F = fopen(period_filename, "r");
                for(int i = 0; i < n+1; i++){
                    err = fscanf(F, "%lf, ", &period_array[i]);
                    if(err == EOF)
                        err_exit("ERROR: Issue while reading period-file", err);
                }

                printf("Read period values correctly\n");
            }
            else{
                period_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(period_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                for(int i = 0; i < n+1; i++)
                    period_array[i] = 1.;

                printf("Wrote period values correctly\n");
            }

            if(constf_found){
                //read amplitude values
                if(fopen(constf_filename, "r") == NULL){
                    printf("ERROR: Cannot open file containing amplitude values. Aborting\n");
                    exit(-1);
                }
                constf_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(constf_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                F = fopen(constf_filename, "r");
                for(int i = 0; i < n+1; i++){
                    err = fscanf(F, "%lf, ", &constf_array[i]);
                    if(err == EOF)
                        err_exit("ERROR: Issue while reading const-file", err);
                }

                printf("Read constant values correctly\n");
            }
            else{
                constf_array = (double *)malloc(((size_t)n+1) * sizeof(double));
                if(constf_array == NULL){
                    printf("ERROR: Not enough space to allocate offset value array. Aborting\n");
                    exit(-1);
                }
                for(int i = 0; i < n+1; i++)
                    constf_array[i] = 1.;

                printf("Wrote constant values correctly\n");
            }

            //Calculation

            printf("Starting Calculation\n");

            uint_fast64_t valCount = 0;

            for(uint64_t i = 0; i < steps; i++){
                double prozent = 100* (double) i / (double) steps;
                printf("\r%lf%%, step %" PRIu64" of %u", prozent, i, steps);
                timeVariantProbs(jump_rate_Array, n+1, period_array, offset_array, amplit_array, t, funTyp_array, constf_array);
				for(uint_fast64_t k = 0; k < dim; k++){
				    valCount = gen_meq_row(n, k, jump_rate_Array, valArr, colInd);
                    rk4_step_reduced(1, valCount, k, &k1[k], x, x1, colInd, valArr, h);
                    rk4_step_reduced(2, valCount, k, &k2[k], x, x1, colInd, valArr, h);
                    rk4_step_reduced(3, valCount, k, &k3[k], x, x1, colInd, valArr, h);
                    rk4_step_reduced(4, valCount, k, &k4[k], x, x1, colInd, valArr, h);
                }
                for(unsigned int k = 0; k < dim; k++){
                    x[k] += (k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.;
                    checksum = fabs((k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.);
                }
                t += h;

                if(a_found){
                    animated_step_check = i % animation_resolution;
                    if(animated_step_check == 0){
                        //Calculate exprected values of x
                        
						exp_val_eval(dim, exp_values, x);	
					    
						for(int ctr = 0; ctr < n; ctr++){
                            expValueAnimation[(n*expValuePtr) + ctr] = exp_values[ctr];
                        }

                        expValuePtr++;
						
						cov_matrix_eval(n, dim, exp_values, x, matrix_data, multExp_matrix, &animation_data_pointer);
                       
                        for(uint64_t ctr = 0; ctr < dim; ctr++){
                            prob_values[prob_value_data_pointer * dim + ctr] = x[ctr];
                        }
                        prob_value_data_pointer++;
                    }
                }
            }
            printf("\n");

        }

        else{
            //TODO TEST THIS!!!
            uint_fast64_t valCount = 0;
			for(uint64_t i = 0; i < steps; i++){
                checksum = 0;
            	
			 	for(uint_fast64_t k = 0; k < dim; k++){
				    valCount = gen_meq_row(n, k, jump_rate_Array, valArr, colInd);
                    rk4_step_reduced(1, valCount, k, &k1[k], x, x1, colInd, valArr, h);
                    rk4_step_reduced(2, valCount, k, &k2[k], x, x1, colInd, valArr, h);
                    rk4_step_reduced(3, valCount, k, &k3[k], x, x1, colInd, valArr, h);
                    rk4_step_reduced(4, valCount, k, &k4[k], x, x1, colInd, valArr, h);
                }

                double checksum = 0;
                for(unsigned int k = 0; k < dim; k++){
                    x[k] += (k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.;
                    checksum += fabs((k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.);
                }

                t += h;

                
                if(a_found){
                    animated_step_check = i % animation_resolution;
                    if(animated_step_check == 0){
                        //Calculate exprected values of x
                       	exp_val_eval(dim, exp_values, x);

						for(int ctr = 0; ctr < n; ctr++){
                            expValueAnimation[(n*expValuePtr) + ctr] = exp_values[ctr];
                        }

                        expValuePtr++;

						cov_matrix_eval(n, dim, exp_values, x, matrix_data, multExp_matrix, &animation_data_pointer);
                        
                        for(uint64_t ctr = 0; ctr < dim; ctr++){
                            prob_values[prob_value_data_pointer * dim + ctr] = x[ctr];
                        }
                        prob_value_data_pointer++;
                    }
                }

                printf("Simulation time elapsed: %lf ------ Delta %.5g ------ Stop condition reached %d%%\r", t, checksum, early_stop);
                if(checksum < stop_condition * h){
                    early_stop++;
                    steps_at_stop = i;
                }
                else{
                    early_stop = 0;
                }
                if(early_stop > 100){
                    break;
                }

            }
        }        
    }
    free(k1);
    free(k2);
    free(k3);
    free(k4);    
    free(x1);
    free(valArr);
    free(colInd);

    gettimeofday(&stop, NULL);
    secs = (double)(stop.tv_usec - start.tv_usec)/1000000 + (double)(stop.tv_sec - start.tv_sec);
    if(fopen("runtimes.txt", "a") != NULL){
      F = fopen("runtimes.txt", "a");
      fprintf(F, "%i %f\n", n, secs);
      fclose(F);
    }
    printf("\n####Finished Calculation####\n");
    printf("Real Time elapsed: %f\n", secs);
    printf("Simulation time elapsed: %lf\n", t);

    if(early_stop){
        printf("Convergence after %" PRIu64 " steps\nConfidience |\u0394x| = %.20g\n", steps_at_stop, stop_condition*h);
    }
     
    
    checksum = 0;
    for(uint64_t i = 0; i < dim; i++)
        checksum += x[i];
    
    printf("Checksum over result vector: %lf\n", checksum);
    
    
    double covMatrix[n][n];
    int size;
    
    if(of_found){
        //probabilities
        if(save_probability_data){
            if(fopen(o_filename, "w") != NULL){
                F = fopen(o_filename, "w");
                printf("Writing result vector to file '%s'\n", o_filename);
                for(unsigned int i = 0; i < dim; i++){
                    prozent = 100 * (double)(i+1)/(double)(dim);
                    printf("-->%f%%\r", prozent);
                    //fwrite(&i, sizeof(unsigned int), 1, F);
                    unsigned int k = i;
                    for (int j = 0;j<n; j++){
                        if(((k>>(n-j-1))&1)==1)
                            fprintf(F,"%c", '1');
                        else
                            fprintf(F,"%c", '0');
                    } //fprintf(F,"%c",'0' + (k & 1));
                    fprintf(F, "\t\t%.20g\n", x[i]);
                }
            printf("\n");
            fclose(F);
            }
            size = strlen(o_filename);
            size += 14;
            char anim_prob_txt[size];
            snprintf(anim_prob_txt, size, "%s_prb.txt", o_filename);
            
            if(fopen(anim_prob_txt, "w") != NULL){
                F = fopen(anim_prob_txt, "w");


            
                for(uint64_t i = 0; i < prob_value_data_pointer; i++){
                    double time = h * animation_resolution * i; 
                    fprintf(F, "%f, ", time);
                    for(uint64_t k = 0; k < dim - 1; k++){
                        fprintf(F, "%f, ",prob_values[i * dim + k]);
                    }
                    fprintf(F, "%f\n", prob_values[i * dim + dim - 1]);
                }
                fclose(F);
            }
            else
            {
                printf("ERROR: Cant write to %s\n", anim_prob_txt);
            }
            
        }

        //Exp Values
        if(save_expected_value_data){
            int size = strlen(o_filename);
            size += 20; 
            char expected_values_filename[size];

            snprintf(expected_values_filename, size, "%s_exp_values_fin.vec", o_filename);   
            double * y = condense(x, dim);

            if(fopen(expected_values_filename, "w") != NULL){
                F = fopen(expected_values_filename, "w");
                for(unsigned int i = 0; i < log2_u64(dim) - 1; i++){
                    fprintf(F, "%lf, ", y[n-1-i]);
                }

                fprintf(F, "%lf", y[0]);
            }
            fclose(F);

            char gnuPlot_data_expvals_filename[size];
            snprintf(gnuPlot_data_expvals_filename, size, "%s_exp_values_fin.dat", o_filename);
            if(fopen(gnuPlot_data_expvals_filename, "w") != NULL){
                F = fopen(gnuPlot_data_expvals_filename, "w");
                for(int i = 0; i < n ; i++){
                    fprintf(F, "%d Pos%d %lf\n",i,i , y[n-1-i]);
                }
            }
            fclose(F);

            if(a_found){
                size = strlen(o_filename);
                size += 20;

                char anim_exp_txt[size];
                snprintf(anim_exp_txt,  size, "%s_exp_timed.txt", o_filename);
                if(fopen(anim_exp_txt, "w") != NULL){
                    F = fopen(anim_exp_txt, "w");

                    for(unsigned int i = 0; i < expValuePtr; i++){
                        double time = h * animation_resolution * i; 
                        fprintf(F, "%f, ", time);
                        for(int k = 0; k < n - 1; k++){
                            fprintf(F, "%lf, ", expValueAnimation[(n*i) + n - 1 - k]);
                        }
                        fprintf(F, "%lf\n", expValueAnimation[(n*i)]);
                    }
                    fclose(F);
                }
            }
            free(y);
        }

        //Covarinace 
        if(save_covariance_data){
            size -= 5;
            char cavrianceMatrix_ppm_filename[size];
            snprintf(cavrianceMatrix_ppm_filename, size, "%s_cov.mat", o_filename);
            double * y = condense(x, dim);

            for(int i = 0; i < n; i++){
                for(int k = 0; k < n; k++){
                    for(uint64_t s = 0; s < dim; s++){
                        temp = s;
                        xi = (temp >> i)&1;
                        temp = s;
                        xk = (temp >> k)&1;
                        cov += (xk - y[k])*(xi - y[i]) * x[s];
                    }
                    covMatrix[n-1-i][n-1-k] = cov;
                    cov = 0;
                }
            }

            if(fopen(cavrianceMatrix_ppm_filename, "w") != NULL){
                F = fopen(cavrianceMatrix_ppm_filename, "w");
                for(int i = 0; i < n; i++){
                    for(int k = 0; k < n; k++){
                        if(k < n-1){
                            fprintf(F, "%lf, ", covMatrix[i][k]);
                        }
                        else{
                            fprintf(F, "%lf\n", covMatrix[i][k]);
                        }
                    }
                }
                fclose(F);
            }
            free(y);
            //create Covariance picture

            size = strlen(o_filename);
            size += 9;
            char of[size];
            snprintf(of, size, "%s_cov.ppm", o_filename);

            screen_t *S;
            S = create_screen(n, n);
            clear_screen(S);
            pixel_t pixel;
            int r,g,b;
            
            double intensity = 0.25;
            for(int i = 0; i < n; i++){
                for(int k = 0; k < n; k++){
                    cov = covMatrix[i][k];
                    if(cov < 0 && fabs(cov) > 0.0001){
                        r = 0;
                        g = 0;
                        b = (-1)*pow(fabs(cov), intensity) * 255;
                    }
                    else if(cov > 0 && fabs(cov) > 0.0001){
                        r = pow(cov, intensity) * 255;
                        g = 0;
                        b = 0;
                    }
                    else{
                        r = 0;
                        g = fabs(pow(fabs(cov), intensity) * 255);
                        b = 0;
                    }
                    pixel.r = r;
                    pixel.b = b;
                    pixel.g = g;
                    set_pixel(S, i, k, pixel);
                }
            }
            
            write_ppm(S, of);

            if(a_found){
                int size = strlen(o_filename);
                size += 17; 
                char animation_directory_name[size];
                
                printf("Creating directory for covariance animation data\n");
                snprintf(animation_directory_name, size, "%s_covariance_data", o_filename);
                int dir_check = mkdir(animation_directory_name, 0777);
                if(dir_check == -1){
                    printf("Cauton: Directory <%s> already exists, Data may be overriten or lost! Continue? [Y/N]  ", animation_directory_name);
                    char sel = 'N';
                    err = scanf("%c", &sel);
                    if(err == EOF)
                        err_exit("ERROR: Issue while reading from stdin", err);
                    if(sel == 'N' || sel == 'n'){
                        exit(-1);
                    }
                }
                size = strlen(animation_directory_name);
                size += strlen(o_filename);
                size += 20;
                char anim_of_ppm[size];
                snprintf(anim_of_ppm,   size, "%s/%s_fin.ppm", animation_directory_name, o_filename);
                write_ppm(S, anim_of_ppm);
                printf("Wrote final covariance matrix ppm file to <%s>\n", anim_of_ppm);

                char anim_of_txt[size];
                snprintf(anim_of_txt,   size, "%s/%s_all.txt", animation_directory_name, o_filename);
                
                if(fopen(anim_of_txt, "w") != NULL){
                    F = fopen(anim_of_txt, "w");
                    int matrices = animation_data_pointer / n;
                    for(int j = 0; j < matrices; j++){
                        for(int i = 0; i < n; i++){
                            for(int k = 0; k < n; k++){
                                cov = matrix_data[((j+1) * n * n) - 1 - k - (i * n)];
                                if(cov < 0 && fabs(cov) > 0.0001){
                                    r = 0;
                                    g = 0;
                                    b = (-1)*pow(fabs(cov), intensity) * 255;
                                }
                                else if(cov > 0 && fabs(cov) > 0.0001){
                                    r = pow(cov, intensity) * 255;
                                    g = 0;
                                    b = 0;
                                }
                                else{
                                    r = 0;
                                    g = fabs(pow(fabs(cov), intensity) * 255);
                                    b = 0;
                                }
                                pixel_t pixel;
										  pixel.r = r;
										  pixel.g = g;
										  pixel.b = b;
                                set_pixel(S, i, k, pixel);
                            }
                            if((i+1) % n != 0){

                            }
                            else{
                                double time = h * animation_resolution * j; 
                                snprintf(anim_of_ppm, size, "%s/%s_%lf.ppm", animation_directory_name, o_filename, time);
                                write_ppm(S, anim_of_ppm);

                                fprintf(F, "%d Matrix:\n", j);
                                for(int _i = 0; _i < n; _i++){
                                    for(int _k = 0; _k < n; _k++){
                                        if(_k < n-1){
                                            fprintf(F, "%lf, ", matrix_data[(j * n * n) + n-1-_k + ((n-1-_i) * n)]);
                                        }
                                        else{
                                            fprintf(F, "%lf\n", matrix_data[(j * n * n) + n-1-_k + ((n-1-_i) * n)]);
                                        }
                                    }
                                }
                                fprintf(F, "\n\n");              
                            }
                        }
                    }
                    fclose(F);
                }        
            }
            destroy_screen(S);
        }

        //Multiplicative expected values
        if(save_multiplicative_expVal_data){
            size += 1;
            char mevalMatrix_ppm_filename[size];
            snprintf(mevalMatrix_ppm_filename, size, "%s_mulExpVal.mat", o_filename);
            double * y = condense(x, dim);
            double multExpvalMatrix [n][n];
            double meval = 0;
            for(int i = 0; i < n; i++){
                for(int k = 0; k < n; k++){
                    for(uint64_t s = 0; s < dim; s++){
                        temp = s;
                        xi = (temp >> i)&1;
                        temp = s;
                        xk = (temp >> k)&1;
                        meval += (xk - y[k])*(xi - y[i]) * x[s];
                    }
                    multExpvalMatrix[n-1-i][n-1-k] = y[k] * y[i] + meval;
                    meval = 0;
                }
            }
            if(fopen(mevalMatrix_ppm_filename, "w") != NULL){
                F = fopen(mevalMatrix_ppm_filename, "w");
                for(int i = 0; i < n; i++){
                    for(int k = 0; k < n; k++){
                        if(k < n-1){
                            fprintf(F, "%lf, ", multExpvalMatrix[i][k]);
                        }
                        else{
                            fprintf(F, "%lf\n", multExpvalMatrix[i][k]);
                        }
                    }
                }
                fclose(F);
            }
            free(y);
            

            char of[size];
            snprintf(of, size, "%s_mulExpVal.ppm", o_filename);

            screen_t *S;
            S = create_screen(n, n);
            clear_screen(S);
            pixel_t pixel;
            int r,g,b;
            
            double intensity = 0.25;
            for(int i = 0; i < n; i++){
                for(int k = 0; k < n; k++){
                    meval = multExpvalMatrix[i][k];
                    r = (-1)*pow(fabs(meval), intensity) * 255;
                    g = (-1)*pow(fabs(meval), intensity) * 255;
                    b = (-1)*pow(fabs(meval), intensity) * 255;
                    
                    pixel.r = r;
                    pixel.b = b;
                    pixel.g = g;
                    set_pixel(S, i, k, pixel);
                }
            }
            
            write_ppm(S, of);

            if(a_found){
                int size = strlen(o_filename);
                size += 17; 
                char directory_name[size];
                printf("Creating directory for mulitplicative expected values animation data\n");
                snprintf(directory_name, size, "%s_multExpval_data", o_filename);
                int dir_check = mkdir(directory_name, 0777);
                if(dir_check == -1){
                    printf("Cauton: Directory <%s> already exists, Data may be overriten or lost! Continue? [Y/N]  ", directory_name);
                    char sel = 'N';
                    err = scanf("%c", &sel);
                    if(err == EOF)
                        err_exit("ERROR: Issue while reading from stdin", err);
                    if(sel == 'N' || sel == 'n'){
                        exit(-1);
                    }
                }
                size = strlen(directory_name);
                size += strlen(o_filename);
                size += 20;
                
                char anim_of_ppm[size];
                snprintf(anim_of_ppm,   size, "%s/%s_fin.ppm", directory_name, o_filename);
                write_ppm(S, anim_of_ppm);
                printf("Wrote final multiplicative expected values matrix ppm file to <%s>\n", anim_of_ppm);

                char anim_of_txt[size];
                snprintf(anim_of_txt,   size, "%s/%s_all.txt", directory_name, o_filename);
                
                if(fopen(anim_of_txt, "w") != NULL){
                    F = fopen(anim_of_txt, "w");
                    int matrices = animation_data_pointer / n;
                    for(int j = 0; j < matrices; j++){
                        for(int i = 0; i < n; i++){
                            for(int k = 0; k < n; k++){
                                meval = multExp_matrix[((j+1) * n * n) - 1 - k - (i * n)];
                                r = (-1)*pow(fabs(meval), intensity) * 255;
                                g = (-1)*pow(fabs(meval), intensity) * 255;
                                b = (-1)*pow(fabs(meval), intensity) * 255;
                                pixel.r = r;
                                pixel.b = b;
                                pixel.g = g;
                                //pixel = {r,g,b};
                                set_pixel(S, i, k, pixel);
                            }
                            if((i+1) % n != 0){

                            }
                            else{
                                double time = h * animation_resolution * j; 
                                snprintf(anim_of_ppm, size, "%s/%s_%lf.ppm", directory_name, o_filename, time);
                                write_ppm(S, anim_of_ppm);

                                fprintf(F, "%d Matrix:\n", j);
                                for(int _i = 0; _i < n; _i++){
                                    for(int _k = 0; _k < n; _k++){
                                        if(_k < n-1){
                                            fprintf(F, "%lf, ", multExp_matrix[(j * n * n) + n-1-_k + ((n-1-_i) * n)]);
                                        }
                                        else{
                                            fprintf(F, "%lf\n", multExp_matrix[(j * n * n) + n-1-_k + ((n-1-_i) * n)]);
                                        }
                                    }
                                }
                                fprintf(F, "\n\n");              
                            }
                        }
                    }
                    fclose(F);
                }        
            }
            destroy_screen(S);
        }
    }

    int mask_size, mask;
    mask_size = 2;
    mask = 0;
    printf("Calculating Higher Order Correlations\n");
    size = strlen(o_filename);
    size += 11;
    char high_order_filename[size];

    snprintf(high_order_filename, size, "%s_highcorr", o_filename);
    F = NULL;
    if(fopen(high_order_filename, "w") != NULL)
    {   
        F = fopen(high_order_filename, "w");
        for(uint_fast32_t i = 1; i < n-mask_size; i++)
            fprintf(F,"\t\t\t%ld", i);
        fprintf(F,"\n");
        for(uint_fast32_t k = 0; k < POW2(mask_size); k++){

            for(uint_fast32_t i = 0; i < mask_size; i++)
                mask>>(mask_size-1-i) & 1UL ? fprintf(F,"1") : fprintf(F,"0"); 
            fprintf(F,"\t");
            for(uint_fast32_t i = 1; i < n-mask_size; i++){
                double corr = high_order_corr(x, dim, mask, mask_size, n-mask_size-i, n);
                fprintf(F,"%1.4g\t",i , corr);
            }
            fprintf(F,"\n");
            mask++;
        }
    }

    fclose(F);

    //ERROR GUESSING
    mask = POW2(mask_size)-1;

    printf("Calculating error guess\n");
    size = strlen(o_filename);
    size+=13;
    char error_guess_filename[size];

    snprintf(error_guess_filename, size, "%s_error_guess", o_filename);

    F=NULL;
    if(fopen(error_guess_filename, "w")!= NULL)
    {
        F = fopen(error_guess_filename, "w");
        for(uint_fast32_t i = 1; i < n-mask_size; i++){
            double e_guess = error_guess(x, dim, mask, mask_size, n-mask_size-i, n);
            fprintf(F,"Error guess for t_%ld...t_%ld:\t %lf\n", i-1, i+mask_size, e_guess);
        }
    }

    fclose(F);

    printf("Calculating closures\n");
    size = strlen(o_filename);
    size+=10;
    char closure_approx_filename[size];

    snprintf(closure_approx_filename, size, "%s_closures", o_filename);

    F=NULL;
    if(fopen(closure_approx_filename, "w")!= NULL)
    {
        F = fopen(closure_approx_filename, "w");
        for(uint_fast32_t i = 1; i < n-mask_size; i++){
            double e_guess = higher_order_pf_approx(x, dim, mask, mask_size, n-mask_size-i, n);
            fprintf(F,"Value for closure of t_%ld...t_%ld:\t %lf\n", i-1, i+mask_size, e_guess);
        }
    }

    fclose(F);

    //Higher order Point Function Values
    mask_size+=2;
    mask = POW2(mask_size)-1;
    size = strlen(o_filename);
    size += 22;
    char ho_pf_fn[size];
    snprintf(ho_pf_fn, size, "%s_high_order_point_func", o_filename);

    F=NULL;
    if(fopen(ho_pf_fn, "w")!=NULL)
    {
        F = fopen(ho_pf_fn, "w");
        for(uint_fast64_t i = 0; i < n-mask_size+1; i++)
        {
            double pf_val = higher_order_point_func(x, dim, mask, mask_size, n-mask_size-i, n);
            fprintf(F, "t_%ld...t_%ld = %lf\n", i, mask_size+i-1, pf_val);
        }
    }

    fclose(F);

    if(a_found){
        size = strlen(o_filename);
        size += 16;

        char anim_high_order_filename[size];
        snprintf(anim_high_order_filename,  size, "%s_highcorr_timed.txt", o_filename);
        if(fopen(anim_high_order_filename, "w") != NULL){
            F = fopen(anim_high_order_filename, "w");

            for(uint_fast64_t i = 0; i < prob_value_data_pointer; i++){
                prozent = 100. * (double)(i+1)/(double)prob_value_data_pointer;
                printf("-->%f%%\r", prozent);

                double time = h * animation_resolution * i; 
                fprintf(F, "%f ", time);
                for(uint_fast32_t k = 1; k < n-mask_size; k++){
                    double corr = high_order_corr(&prob_values[i*dim], dim, mask, mask_size, n-mask_size-k, n);
                    fprintf(F,"%1.7g ", corr);
                }
                fprintf(F, "\n");
            }
            fclose(F);
        }
    }

    free(x);
    free(jump_rate_Array);
    if(if_found)
        free(i_filename);
    if(of_found)
        free(o_filename);
    if(a_found){
        free(matrix_data);
        free(multExp_matrix);
        free(exp_values);
        free(expValueAnimation);
        free(prob_values);
    }
    if(funtyp_found){
        free(funTyp_filename);
    }
    if(offset_found){
        free(offset_filename);
    }
    if(period_found){
        free(period_filename);
    }
    if(amplit_found){
        free(amplit_filename);
    }
    if(constf_found){
        free(constf_filename);
    }
    if(time_var_found){
        free(funTyp_array);
        free(offset_array);
        free(period_array);
        free(amplit_array);
        free(constf_array);
    }
    if(save_select_found)
        free(save_selection);
}
