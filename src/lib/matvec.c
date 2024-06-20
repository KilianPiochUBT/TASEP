#include "matvec.h"

static inline uint64_t log2_u64(uint64_t n)
{
  int result = 0;
  while (n >>= 1) ++result;
  return result;
}

double * condense(double *x, uint64_t sizeofx)
{
    int c_size = log2_u64(sizeofx);
    double *condensed_x;
    condensed_x = (double *) malloc(((size_t)c_size) * sizeof(double));
    if(condensed_x == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate condensed value array.\n");
        return NULL;
    }

    for(int i = 0; i < c_size; i++)
    {
        condensed_x[i] = 0.;
    }
    
    
    for(uint64_t k = 0; k < sizeofx; k++)
    {
        for(int i = 0; i < c_size; i++){
            uint64_t temp = k;
            if((temp >> i)&1){
                condensed_x[i] += x[k]; 
                //printf("Add %d||%d\n", i, temp);
            }
            else{
                //printf("     No %d||%d\n", i, temp);
            }
        }
        //printf("\n\n");
    }
    return condensed_x;
}

void timeVariantProbs(double * jump_rate_Array, int chainLength, double * period, double * offset, double * amplitude ,double time, int* funType, double* constants){
    for(int i = 0; i < chainLength; i++){
        switch(funType[i]){
            case 0:
                jump_rate_Array[i] = amplitude[i] * sin(time * period[i] + offset[i]) + constants[i];
                if(jump_rate_Array[i] < 0){
                    jump_rate_Array[i] = 0.;
                }
                break;
            case 1:
                jump_rate_Array[i] = amplitude[i] * cos(time * period[i] + offset[i]) + constants[i];
                if(jump_rate_Array[i] < 0){
                    jump_rate_Array[i] = 0.;
                }
                break;
            case 2:
                jump_rate_Array[i] = amplitude[i] * time + constants[i];
                if(jump_rate_Array[i] < 0.){
                    jump_rate_Array[i] = 0.;
            }
        }
    }
	
	return;
}

void
exp_val_eval(uint64_t dim, double * exp_values, double * x){
	//Calculate exprected values of x
	int c_size = log2_u64(dim);

	for(int i = 0; i < c_size; i++)
	{
    	exp_values[i] = 0.;
	}

	for(uint64_t k = 0; k < dim; k++)
	{
    	for(int i = 0; i < c_size; i++){
        	uint64_t temp = k;
        	if((temp >> i)&1)
            	exp_values[i] += x[k]; 
    	}
	}

	return;

}

void
cov_matrix_eval(int n, uint64_t dim, double * exp_values, double * x, double * matrix_data, double * multExp_matrix, uint64_t *animation_data_pointer){
	int temp; 
	int xi, xk;
	double cov;
	for(int row = 0; row < n; row++){
   		for(int col = 0; col < n; col++){
       		for(uint64_t s = 0; s < dim; s++){
       	    	temp = s;
     	    	xi = (temp >> row)&1;
       	    	temp = s;
      	     	xk = (temp >> col)&1;

        	   	cov += (xk - exp_values[col])*(xi - exp_values[row]) * x[s];
       		}

       		matrix_data[((*animation_data_pointer) * n) + col] = cov;
       		multExp_matrix[((*animation_data_pointer) * n) + col] = exp_values[row] * exp_values[col] + cov;

			cov = 0;
   		}
   		(*animation_data_pointer)++;
   	}
	
	return;

}
