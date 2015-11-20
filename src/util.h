#ifndef _UTIL_H_
#define _UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ABS(a) (((a)>0)?(a):(-(a)))

void affiche(int m, int n, double *a, int lda, FILE*flux);
double * alloc_matrix(int m, int n);
void init_matrix(int m, int n, double *a);
void init_rand_matrix(int m, int n, double *a);
void cache_buster();

#endif
