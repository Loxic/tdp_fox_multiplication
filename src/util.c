#include "util.h"

void affiche(int m, int n, double *a, int lda, FILE *flux){
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      fprintf(flux, "[%f]", a[i+j*lda]);
    }
    fprintf(flux, "\n");
  }
  fprintf(flux, "\n");
}

double * alloc_matrix(int m, int n){
  double *r = malloc(sizeof(*r)*m*n);
  return r;
}

void init_matrix(int m, int n, double *a){
  for(int i=0; i<m*n; i++){
    a[i] = 0.0f;
  }
}

void init_rand_matrix(int m, int n, double *a){
  srand(time(NULL));
  for(int i=0; i<m*n; i++){
    a[i]=rand()%100;
  }
}

void cache_buster(){
  int size = 20*1024*1024;
  char *c = malloc(size);
  for(int i=0; i<size; i++){
    c[i] = i;
  }
}
