#include "util.h"



void affiche(const int m,const int n,const double *a,const int lda, FILE *flux){
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      fprintf(flux, "[%f]", a[i+j*lda]);
    }
    fprintf(flux, "\n");
  }
  fprintf(flux, "\n");
}

double * alloc_matrix(const int m,const int n){
  double *r = malloc(sizeof(*r)*m*n);
  return r;
}

void init_matrix(const int m,const int n, double *a){
  for(int i=0; i<m*n; i++){
    a[i] = 0.0f;
  }
}

void init_rand_matrix(const int m,const int n, double *a){
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

double *parse_matrix_file(const char * filename, int * size){
  double * matrix;
  FILE * fd;
  fd = fopen(filename, "r");
  fscanf(fd,"%d\n",size);

  matrix = malloc((*size)*(*size)*sizeof(double));
  for(int i = 0; i < *size; i++) {
    for(int j = 0; j < *size; j++) {
      if(fscanf(fd,"%lf",&matrix[i+j*(*size)])) {
      }
    } 
  }
  return matrix;
}

void write_matrix_file(const char *filename, double *matrix, int size){
  FILE *out = fopen(filename, "w");
  fprintf(out, "%d\n", size);
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      fprintf(out, "%lf ", matrix[i*size+j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);
}
