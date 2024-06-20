#include "solv.h"

void 
rk4_step(int step, uint64_t dim, double *k, double *x, double *x1,uint64_t * rowPtr, uint64_t* colInd, double *valArr, double h){
	
    if(USE_OMP==1){ /*use parallel computation on linux machines switch to HAVE_OMP in future*/
        #pragma omp parallel for 
        for(uint64_t i = 0; i < dim; i++){
           	double _T = 0;
           	k[i] = 0;
           	int addMax = rowPtr[i+1] - rowPtr[i];
           	for(int j = 0; j < addMax; j++)
           	    _T += valArr[rowPtr[i] + j] * x[colInd[rowPtr[i] + j]];
           	switch(step){
                case 1:
                case 2:
                    x1[i] = (0.5 * h * _T) + x[i];
                    break;
                case 3:
                    x1[i] = (h * _T) + x[i];
                    break;
                default:
                    break;
               }

	    	k[i] = h * _T;
   	    }
    } 
    else{
        for(uint64_t i = 0; i < dim; i++){
           	double _T = 0;
           	k[i] = 0;
           	int addMax = rowPtr[i+1] - rowPtr[i];
           	for(int j = 0; j < addMax; j++)
           	    _T += valArr[rowPtr[i] + j] * x[colInd[rowPtr[i] + j]];
           	switch(step){
                case 1:
                case 2:
                    x1[i] = (0.5 * h * _T) + x[i];
                    break;
                case 3:
                    x1[i] = (h * _T) + x[i];
                    break;
                default:
                    break;
               }

	    	k[i] = h * _T;
   	    }
    }

	return;
}


void
rk4_step_reduced(int step, uint64_t valCount, uint_fast64_t currIndex, double *k, double *x, double *x1, uint64_t* colInd, double *valArr, double h){
	
   	double _T = 0;
    
   	for(int j = 0; j < valCount; j++)
   	    _T += valArr[j] * x[colInd[j]];
   	switch(step){
        case 1:
        case 2:
            *x1 = (0.5 * h * _T) + x[currIndex];
            break;
        case 3:
            *x1 = (h * _T) + x[currIndex];
            break;
        default:
            break;
       }

    *k = h * _T;


	return;
}

uint64_t
gen_meq_row(int n, uint_fast64_t k, double * jump_rate_Array, double *valArr, uint_fast64_t *colInd){
	
	uint64_t offsetDiagCounter = 0;
    if(k >= POW2(n-1)){
        offsetDiagCounter = k - POW2(n-1);
    }
  	//#pragma omp parallel for // Not yet implemented thread safely 
	uint64_t valCounter = 0;
	uint64_t bitSet = k;
	int bit[MAX_CHAIN_LEN];
	for(int j = 0; j < n; j++){
        if((bitSet >> j)&1)
            bit[j] = 1;
        else
            bit[j] = 0;
    }

    if(bit[n-1]){
        colInd[0] = offsetDiagCounter;
        valArr[0] = jump_rate_Array[0];
        valCounter++;
        offsetDiagCounter++;
    }

    colInd[valCounter] = k;
    valArr[valCounter] = 0.;
    if(bit[0]){
        valArr[valCounter] -= jump_rate_Array[n];
    }
    if(!bit[n-1]){
        valArr[valCounter] -= jump_rate_Array[0];
    }
    for(int j = 0; j < n - 1; j++){
        if(!bit[j] && bit[j+1]){
            valArr[valCounter] -= jump_rate_Array[n-1-j];
        }
    }
    valCounter++;
    if(!bit[0]){
        colInd[valCounter] = k + 1;
        valArr[valCounter] = jump_rate_Array[n];
        valCounter++;
    }
    for(int j = 0; j < n - 1; j++){
        if(bit[j] && !bit[j+1]){
            colInd[valCounter] = k + POW2(j);
            valArr[valCounter] = jump_rate_Array[n-1-j];
            valCounter++;
        }
    }

	
	return valCounter;
}
