#include "steadystates_gmp.h"

static void 
usage(const char *progname)
{
    printf( "Usage of %s: \n"
            "-n, --nodes                Number of nodes in chain\n"
            "-a, --alpha                Entry hop rate \u03b1\n"
            "-b, --beta                 Exit jop rate \u03b2\n"
            "-i,                        First Node Index for correlation calculation\n"
            "-j                         Second Node Index for correlatoin calculation\n"
            "-o, --output               Name of output file to save Expected occupancy values to\n"
            "-p, --precision            Set minimal bit-width for numbers used in calculation\n"
            "--pp --print-precision     Set precision for printed results\n"
            "-v, --verbose              Print Values of B and Z aswell\n"
            "-h, --help                 Display this dialog\n", progname);
    return;
}

void eval_B_gmp(uint_fast64_t n, uint_fast64_t p, mpf_t * array_entry){
    
    //Avoid unnecesary divisions, because n >= p
    uint_fast32_t np_diff = n-p;

    
    if(np_diff == 0){
        //if 0, we know B = 1
        mpf_set_ui(*array_entry, 1UL);
        return;
    }
    if(np_diff == 1){
        mpf_set_ui(*array_entry, (unsigned int)(n-1));
        return;
    }

    mpf_t dividend, divisor, quotient, remainder;

    mpf_init(dividend);
    mpf_set_ui(dividend, 1UL);
    mpf_init(divisor);
    mpf_set_ui(divisor, 1UL);
    mpf_init(quotient);
    mpf_set_ui(quotient, 0UL);
    mpf_init(remainder);
    mpf_set_ui(remainder, 0UL);

    for(unsigned long k = 0; k < np_diff-1; k++){
        mpf_mul_ui(dividend, dividend, n + np_diff - 1 - k);
        mpf_mul_ui(divisor, divisor, np_diff - k);
    }
    
    mpf_mul_ui(dividend, dividend, p);

    mpf_div(quotient, dividend, divisor);
     
    mpf_clear(divisor);
    mpf_clear(dividend);    
    mpf_clear(remainder);

    mpf_set(*array_entry, quotient);

    mpf_clear(quotient);

    return;
}
    

uint_fast8_t 
allocfill_B(mpf_t *** B, uint32_t n){
    *B = (mpf_t **)malloc((size_t) n * sizeof(mpf_t *));
    
    if(*B == NULL)
        return 2;

    for(uint_fast32_t i = 0; i < n; i++){
        (*B)[i] = (mpf_t*)malloc((size_t) (n+1) * sizeof(mpf_t));
        for(unsigned int k = 0; k < n+1; k++)
            mpf_init_set_ui((*B)[i][k], 0UL);

        if((*B)[i] == NULL)
            return 2;
    }

    for(int_fast32_t i = 0; i < n; i++){
        for(int_fast32_t k = 0; k <= n; k++){
            if(k > 0 && k <= i+1)
                eval_B_gmp(i+1, k, (*B)[i][k]);
            
        }
        //Todo Verbose flag
        //fprintf(stdout, "\rB [%ld / %d]", i+1, n);
    }

    return 0;
}

uint_fast8_t
allocfill_Z(mpf_t ** Z, mpf_t ** B, uint32_t n, double alpha, double beta){
    
    *Z = (mpf_t *) malloc((size_t) (n+1) * sizeof(mpf_t));

    for(unsigned int k = 0; k < n+1; k++)
        mpf_init_set_d((*Z)[k], 0.);

    if(*Z == NULL)
        return 2;

    mpf_t power_alpha;
    mpf_init_set_d(power_alpha, 0.);

    mpf_t power_beta;
    mpf_init_set_d(power_beta, 0.);
    

    for(int_fast32_t i = 0; i < n; i++){
        for(int_fast32_t p = 0; p <= i+1; p++){
            for(int_fast32_t q = 0; q <= p; q++){
                mpf_set_d(power_beta, (1./beta));
                mpf_set_d(power_alpha, (1./alpha));

                mpf_pow_ui(power_alpha, power_alpha, q);
                mpf_pow_ui(power_beta, power_beta, p-q);

                mpf_mul(power_alpha, power_alpha, power_beta);
                mpf_mul(power_alpha, power_alpha, B[i][p]);

                mpf_add((*Z)[i+1], (*Z)[i+1], power_alpha);
            }
        }
        //flag
        //fprintf(stdout, "\rZ [%ld / %d]", i+1, n);
    }
    
    mpf_set_ui((*Z)[0], 1UL);
    
    mpf_clear(power_alpha);
    mpf_clear(power_beta);

    return 0;
}

uint_fast8_t
allocfill_Rho(mpf_t ** Rho, mpf_t ** B, mpf_t * Z, uint32_t n, double beta){
    
    *Rho = (mpf_t *) malloc((size_t) n * sizeof(mpf_t));
    for(unsigned int k = 0; k < n; k++)
        mpf_init((*Rho)[k]);

    if(*Rho == NULL)
        return 2;

    mpf_t power_beta;
    mpf_init(power_beta);

    mpf_t acc;
    mpf_init(acc);
    mpf_set_ui(acc, 0UL);

    for(int_fast32_t i = 0; i < n-1 ; i++){
        mpf_set_ui((*Rho)[i], 0UL);
        for(int_fast32_t k = 0; k < n - i - 1; k++){
            mpf_mul(acc, B[k][1], Z[n-1-k]);
            mpf_add((*Rho)[i], (*Rho)[i], acc);
            
            mpf_set_ui(acc, 1UL);
            mpf_set_d(power_beta, 1./beta);
            mpf_pow_ui(power_beta, power_beta, k+2);
            mpf_mul(acc, acc, power_beta);
            mpf_mul(acc, acc, B[n-i-2][k+1]);
            mpf_mul(acc, acc, Z[i]);
            mpf_add((*Rho)[i], (*Rho)[i], acc);
        }
        mpf_div((*Rho)[i], (*Rho)[i], Z[n]);
        //flag
        //fprintf(stdout, "\rB [%ld / %d]", i+1, n);
    }
    
    mpf_set_d(power_beta, beta);
    mpf_mul(acc, power_beta, Z[n]);
    mpf_div((*Rho)[n-1], Z[n-1], acc);
    //fprintf(stdout, "\rB [%d / %d]", n, n);
    mpf_clear(power_beta);
    mpf_clear(acc);

    return 0;
}

double two_point_function(uint_fast32_t i, uint_fast32_t j, double alpha, double beta, int N)
{
    assert(i<j);
    assert(j<N);
    
    mpf_t fact_1, fact_2, fact_3, fact_4;
    mpf_t p, Z_N, Z_N1p, Z_j1, beta_mpf;
    mpf_t tau_1, tau_2;

    mpf_t first_sum_acc, second_sum_acc;


    //inefficient but whatever
    mpf_t ** B, * Z, *Rho;
    mpf_init(tau_2);
    allocfill_B(&B, j-1);
    allocfill_Z(&Z, B, j-1, alpha, beta);
    allocfill_Rho(&Rho, B, Z, j-1, beta);
    mpf_set(tau_2, Rho[i-1]);
    for(uint_fast32_t k = 0; k < j-1; k++){
        for(int l = 0; l < j; l++)
            mpf_clear(B[k][l]);    
        free(B[k]);
        mpf_clear(Z[k]);
        mpf_clear(Rho[k]);
    }
    free(B);
    free(Z);
    free(Rho);

    //precalculate everything possible
    calculate_Z_i(&Z_N, alpha, beta, N);
    calculate_Z_i(&Z_j1, alpha, beta, j-1);
    ////first sum
    //temp vars
    mpz_t fact_part1, fact_part2, fact_part3;
    mpz_init_set_ui(fact_part1 , 0Ul);        
    mpz_init_set_ui(fact_part2, 0Ul);
    mpz_init_set_ui(fact_part3 , 0Ul);

    mpf_t f_temp1, f_temp2;

    mpf_init_set_ui(first_sum_acc, 0UL);

    for(unsigned long p = 0; p < N-j; p++){
        //first factorial
        mpz_fac_ui(fact_part1, 2UL*p);
        mpz_fac_ui(fact_part2, p);
        mpz_fac_ui(fact_part3, p+1);

        mpz_mul(fact_part2, fact_part2, fact_part3);

        mpf_init(f_temp1);
        mpf_init(f_temp2);

        mpf_set_z(f_temp1, fact_part1);
        mpf_set_z(f_temp2, fact_part2);

        mpf_init(fact_1);
        mpf_div(fact_1, f_temp1, f_temp2);

        //tau_n-p+j-1
        mpf_init(tau_1);
        allocfill_B(&B, N-p-1);
        allocfill_Z(&Z, B, N-p-1,alpha, beta);
        allocfill_Rho(&Rho, B, Z, N-p-1, beta);
        mpf_set(tau_1, Rho[i-1]);
        for(uint_fast32_t k = 0; k < N-p-1; k++){
            for(int l = 0; l < N-p; l++)
                mpf_clear(B[k][l]);
            free(B[k]);
            mpf_clear(Z[k]);
            mpf_clear(Rho[k]);
        }
        free(B);
        free(Z);
        free(Rho);
        //Z Quotient

        calculate_Z_i(&Z_N1p, alpha, beta, N-1-p);
        
        mpf_div(Z_N1p, Z_N1p, Z_N);
        
        mpf_set_d(f_temp1, 1.0);
        mpf_mul(f_temp1, Z_N1p, f_temp1);
        mpf_mul(f_temp1, tau_1, f_temp1);
        mpf_mul(f_temp1, fact_1, f_temp1);
        mpf_add(first_sum_acc, first_sum_acc, f_temp1);
    }


    mpf_init_set_ui(second_sum_acc, 0UL);
    mpf_init_set_ui(fact_2,0UL);    
    mpf_init_set_ui(fact_3,0UL);
    //second sum
    for(unsigned long p = 2; p <= N-j+1; p++)
    {
        //setup
        mpz_init(fact_part1);
        mpz_init(fact_part2);
        mpz_init(fact_part3);
        
        //init factorials
        mpz_fac_ui(fact_part1, 2*(N-j)-p);
        mpz_fac_ui(fact_part2, N-j);
        mpz_fac_ui(fact_part3, N-j+1-p);

        //multiply factorial expressions together
        mpz_mul_ui(fact_part1, fact_part1, p-1);
        mpz_mul(fact_part2, fact_part2, fact_part3);

        //transcribe to float
        mpf_set_z(fact_2, fact_part1);
        mpf_set_z(fact_3, fact_part2);

        //division
        mpf_div(fact_2, fact_2, fact_3);

        //beta exponent
        mpf_init_set_d(beta_mpf,1./beta);
        mpf_pow_ui(beta_mpf, beta_mpf, p);

        //final product
        mpf_mul(fact_2, fact_2, beta_mpf);
        
        //accumulation
        mpf_add(second_sum_acc, second_sum_acc, fact_2);
    }

    mpf_mul(second_sum_acc, second_sum_acc, Z_j1);
    mpf_mul(second_sum_acc, second_sum_acc, tau_2);
    mpf_div(second_sum_acc, second_sum_acc, Z_N);

    mpf_add(second_sum_acc, second_sum_acc, first_sum_acc);

    double return_value; 
    return_value = mpf_get_d(second_sum_acc);
    return return_value;
}


uint_fast8_t calculate_Z_i(mpf_t * Z_i, double alpha, double beta, int i)
{
    mpf_init_set_d(*Z_i, 0.);

    mpf_t ** B, *Z;
    allocfill_B(&B, i);
    allocfill_Z(&Z, B, i, alpha, beta);

    mpf_set(Z_i, Z[i-1]);

    return 0;
}

int main(int argc, char *argv[])
{
    int opt;
    int option_index;

    int_fast8_t of_found = 0;
    char * of_filename;
    
    int_fast8_t n_found = 0;
    uint32_t n;
    char * n_temp;

    int_fast8_t alpha_found = 0;
    double alpha;
    char * alpha_temp;

    int_fast8_t beta_found = 0;
    double beta;
    char * beta_temp;

    int_fast8_t precision_found = 0; 
    uint_fast64_t precision = 64; 
    char * precision_temp;

    int_fast8_t pp_found= 0;
    uint_fast64_t pp = 5;
    char * pp_temp;

    uint_fast8_t err_code = 0;

    uint_fast8_t verbose = 0;

    uint_fast32_t tp_i, tp_j;
    char * tp_temp;
    uint_fast8_t tp_i_found, tp_j_found;

    FILE * F;

    while((opt = getopt_long(argc, argv, "n:a:b:o:p:i:j:vh?", long_opt, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'n':
        case  0 :
            n_found = 1;
            n_temp = strdup(optarg);
            n = atoi(n_temp);
            free(n_temp);
        break;

        case 'a':
        case  1 :
            alpha_found = 1;
            alpha_temp = strdup(optarg);
            alpha = atof(alpha_temp);
            free(alpha_temp);
        break;

        case 'b':
        case  2 :
            beta_found = 1;
            beta_temp = strdup(optarg);
            beta = atof(beta_temp);
            free(beta_temp);
        break;

        case 'o':
        case  3 :
            of_found = 1;
            of_filename = strdup(optarg);
        break;

        case 'p':
        case  4 :
            precision_found = 1;
            precision_temp = strdup(optarg);
            precision = atoi(precision_temp);
            free(precision_temp);
            break;

        case  5 :
            pp_found = 1;
            pp_temp = strdup(optarg);
            pp = atoi(pp_temp);
            free(pp_temp);
            break;
        case  7 :
        case 'v':
            verbose = 1;
            break;
        case 'i':
            tp_temp = strdup(optarg);
            tp_i = atoi(tp_temp);
            tp_i_found = 1;
            break;
        case 'j':
            tp_temp = strdup(optarg);
            tp_j = atoi(tp_temp);
            tp_j_found = 1;
            break;
        case  6 :
        case 'h':
        default :
            usage(argv[0]);
            exit(1);
        }
    }

    /* PRELIMINARY CHECKS FOR VALUE CORRECTNESS*/

    if(!n_found)
        err_exit("Number of nodes not specified. Rerun with -n | --nodes", 1);

    if(!alpha_found)
        err_exit("Alpha value not specified. Rerun with -a | --alpha", 1);
    
    if(!beta_found)
        err_exit("Beta value not specified. Rerun with -b | --beta", 1);
    
    if(!of_found){
        //Give Warning
        fprintf(stderr, "CAUTION: No output file specified. Results will be printed to 'stdout'\n");
    }

    if(!pp_found){
        //fprintf(stdout, "INFO: Results will be given with precision of %ld\n", pp);
    }

    /*END OF PRELIMINARY CHECKS*/
    if(precision_found){
        if(precision < (uint_fast64_t) mpf_get_default_prec()){
            fprintf(stderr, "CAUTION: Precision set to less than default value\n");
            fprintf(stderr, "CAUTION: Precision set to %"PRIuFAST64"\n", precision);
            fprintf(stderr, "CAUTION: Standard precision is %"PRIu64"\n", (uint64_t)mpf_get_default_prec());
        }
        else{
            fprintf(stdout, "Precision set to %"PRIuFAST64" bits\n", precision);
        }
    }

    mpf_set_default_prec(precision); 

    mpf_t ** B;
    mpf_t * Rho;
    mpf_t * Z;

    //B array
    err_code = allocfill_B(&B,n);
    if(err_code)
        err_exit("B array could not be created/filled.", err_code);
    
    fprintf(stdout, "\nAllocated and filled B array\n");

    //Z array
    err_code = allocfill_Z(&Z, B, n, alpha, beta);
    if(err_code)
        err_exit("Z array could not be created/filles.", err_code);

    fprintf(stdout, "\nAllocated and filled Z array\n");

    //rho_array
    err_code = allocfill_Rho(&Rho, B, Z, n, beta);
    if(err_code)
        err_exit("Rho array could not be created/filled.", err_code);

    fprintf(stdout, "\nAllocated and filled Rho array\n");

    //save to specified file

    if(of_found){
        if(fopen(of_filename, "w") == NULL)
            err_exit("Could not open <output> file", 2);
        
        F = fopen(of_filename, "w");
        for(unsigned int k = 0; k < n; k++){
            mpf_out_str(F, 10, 5, Rho[k]);
            fprintf(F, ",");
        }
        fclose(F);


        int size = strlen(of_filename);
        size += 5;
        char data_file[size];
        snprintf(data_file, size, "%s.dat" ,of_filename);
        if(fopen(data_file, "a") == NULL)
            err_exit("Could not open <output>.dat file", 2);
        
        F = fopen(data_file, "w");
        for(unsigned int k = 0; k < n; k++){
            fprintf(F, "%d Pos%d ",k,k);
            mpf_out_str(F, 10, pp, Rho[k]);
            fprintf(F, "\n");
        }
        fclose(F);
        
    }
    else{
        
        for(unsigned int k = 0; k < n; k++){
            fprintf(stdout, "  \u27e8\u03c4_%d\u27e9: ",k+1);
            mpf_out_str(stdout, 10, pp, Rho[k]);
            fprintf(stdout,"\n");
        } 
        fprintf(stdout, "\n");
        
        if(verbose)
        {
            for(unsigned int k = 0; k < n; k++){
                for(unsigned int j = 0; j <= n; j++){
                    fprintf(stdout, "  B_%d,%d: ",k+1,j+1);
                    mpf_out_str(stdout, 10, pp, B[k][j]);
                }
                fprintf(stdout,"\n");
            } 
            fprintf(stdout, "\n");
        
            for(unsigned int k = 0; k < n; k++){
                fprintf(stdout, "  Z_%d: ",k+1);
                mpf_out_str(stdout, 10, pp, Z[k]);
                fprintf(stdout,"\n");
            } 
            fprintf(stdout, "\n");
        }
    }

    if(tp_i_found && tp_j_found){
        double test = two_point_function(tp_i, tp_j, alpha, beta, n);
        double rho_i = mpf_get_d(Rho[tp_i-1]);
        double rho_j = mpf_get_d(Rho[tp_j-1]);
        double sigma_i = sqrt(rho_i-rho_i*rho_i);
        double sigma_j = sqrt(rho_j-rho_j*rho_j);
        double corr = (test-rho_i*rho_j)/(sigma_i*sigma_j);
        printf("Correlation of nodes %d,%d: %lf\n",tp_i, tp_j, corr);
    }

    for(uint_fast32_t k = 0; k < n; k++){
        mpf_clear(Z[k]);
        mpf_clear(Rho[k]);
        for(uint_fast32_t l = 0; l < n+1; l++)
            mpf_clear(B[k][l]);
    }
    mpf_clear(Z[n]);

    for(uint32_t i = 0; i < n; i++)
        free(B[i]);
    free(B);
    free(Z);
    free(Rho);
    if(of_found)
        free(of_filename);
    
    return 0;
}
