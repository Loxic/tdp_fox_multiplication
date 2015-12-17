#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"
#include "mul_fox.h"

#define DEFAULT_SIZE 10

int main(int argc, char**argv){
  MPI_Init(NULL, NULL);
  srand(time(NULL));

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int n = DEFAULT_SIZE;
  if(argc > 1)
    n = atoi(argv[1]);

  if(!is_nbr_proc_valid(n)){
    fprintf(stdout, "Number of proc isn't a multiple of matrix size. (Default matrix size is 10[x10], this can be changed by passing wanted size as program parameter)\n");
    MPI_Finalize();
    return 0;
  }

  double *a=NULL,*b=NULL,*c=NULL,*c_bis=NULL;
  if(myrank==0){
    a = malloc(sizeof(*a)*n*n);
    b = malloc(sizeof(*b)*n*n);
    c_bis = malloc(sizeof(*c_bis)*n*n);
    for(int i=0; i<n*n; i++){
      a[i] = rand()%100;
      b[i] = rand()%100;
      c_bis[i] = 0.0f;
    }
  }

  MPI_Comm proc_grid;
  mul_fox_grid_comm_create(&proc_grid);

  c = mul_fox(0, a, b, n, proc_grid);
  
  if(myrank==0){
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f, a, n, b, n, 1.0f, c_bis, n);

    affiche(n,n,c,n,stdout);
    affiche(n,n,c_bis,n,stdout);

    for(int i=0; i<n*n; i++){
      if(c_bis[i]!=c[i]){
	fprintf(stdout, "Error on matrix comparison at index %d !\n", i);
      }
    }
    fprintf(stdout, "Matrix comparison finished.\n");
  }

  MPI_Finalize();
}
