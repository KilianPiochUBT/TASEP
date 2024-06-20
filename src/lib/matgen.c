#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

#include "matgen.h"

void
gencrs(uint64_t* rowPtr,double* valArr, uint64_t* colInd, double* jump_rate_Array, int* bit, uint64_t dim, int n){

	rowPtr[0] = 0;
    uint64_t offsetDiagCounter = 0;
    double prozent;
	uint64_t valCounter;
	uint64_t bitSet;
	//#pragma omp parallel for
	for(uint64_t i = 0; i < dim; i++){
		prozent = 100. * (double)(i+1)/(double)dim;
        printf("-->%f%%\r", prozent);
        valCounter = 0;
        bitSet = i;

		for(int k = 0; k < n; k++){
            if((bitSet >> k)&1)
                bit[k] = 1;
            else
                bit[k] = 0;
        }

        if(bit[n-1]){
            colInd[rowPtr[i]] = offsetDiagCounter;
            valArr[rowPtr[i]] = jump_rate_Array[0];
            valCounter++;
            offsetDiagCounter++;
        }

		colInd[rowPtr[i]+valCounter] = i;
        valArr[rowPtr[i]+valCounter] = 0.;

		if(bit[0]){
            valArr[rowPtr[i]+valCounter] -= jump_rate_Array[n];
        }

		if(!bit[n-1]){
            valArr[rowPtr[i]+valCounter] -= jump_rate_Array[0];
        }

		for(int k = 0; k < n - 1; k++){
            if(!bit[k] && bit[k+1]){
                valArr[rowPtr[i]+valCounter] -= jump_rate_Array[n-1-k];
            }
        }

		valCounter++;

		if(!bit[0]){
            colInd[rowPtr[i]+valCounter] = i + 1;
            valArr[rowPtr[i]+valCounter] = jump_rate_Array[n];
            valCounter++;
        }

		for(int k = 0; k < n - 1; k++){
            if(bit[k] && !bit[k+1]){
                colInd[rowPtr[i]+valCounter] = i + POW2(k);
                valArr[rowPtr[i]+valCounter] = jump_rate_Array[n-1-k];
                valCounter++;
            }
        }

		rowPtr[i + 1] = rowPtr[i] + valCounter;

	}
	
	return;
}

