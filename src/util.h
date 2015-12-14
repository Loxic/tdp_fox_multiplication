#ifndef _UTIL_H_
#define _UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ABS(a) (((a)>0)?(a):(-(a)))

void affiche(const int m,const int n,const double *a,const int lda, FILE*flux);
double * alloc_matrix(const int m,const int n);
void init_matrix(const int m,const int n, double *a);
void init_rand_matrix(int m, int n, double *a);
void cache_buster();
double *parse_matrix_file(const char * filename, int * size);
void write_matrix_file(const char *filename, double *matrix, int size);

#endif
