#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

void 
timedatestring();

void
err_exit(char * msg, int_fast32_t exit_code);

double * 
enterProbs(int num);


#endif
