#ifndef MATVEC_H
#define MATVEC_H

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double * 
condense(   double *x, 
            uint64_t sizeofx);

void 
timeVariantProbs(   double * jump_rate_Array, 
                    int chainLength, double * period, 
                    double * offset, 
                    double * amplitude ,
                    double time, 
                    int* funType, 
                    double* constants);

void
exp_val_eval(   uint64_t dim, 
                double * exp_values, 
                double * x);

void
cov_matrix_eval(int n, 
                uint64_t dim, 
                double * exp_values, 
                double * x, 
                double * matrix_data, 
                double * multExp_matrix, 
                uint64_t *animation_data_pointer);

#endif
