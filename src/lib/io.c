#include "io.h"

char DATE_TIME[128];

void
timedatestring(void){
    time_t t = time(NULL);
    struct tm * __tm = localtime(&t);
    assert(strftime(DATE_TIME, sizeof(DATE_TIME), "%c", __tm));
    return;
}

void
err_exit(char * msg, int_fast32_t exit_code){
    timedatestring();
    
    fprintf(stderr, "[%s]:%s\n", DATE_TIME, msg);
    exit(exit_code);
}

double * enterProbs(int num){
    double *ret = (double*) malloc((size_t)num * sizeof(double));
    if(ret == NULL){
        printf("ERROR: Can not allocate memory for probability array. Aborting\n");
        exit(-1);
    }
    for(int i = 0; i < num; i++){
        int err = scanf("%lf", &ret[i]);
        if(err == EOF)
            err_exit("ERROR: Could not read from stdin.", err);
    }

    return ret;
}
