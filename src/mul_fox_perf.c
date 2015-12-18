#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"
#include "mul_fox.h"

int main(int argc, char**argv){
  if(argc<3){
    fprintf(stdout, "Usage is %s [matrix size] [nb smoothing iteration].\n", argv[0]);
      return 0;
  }

  MPI_Init(NULL, NULL);
  
  int myrank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  int n = atoi(argv[1]);
  int itr = atoi(argv[2]);

  if(!is_nbr_proc_valid(n)){
    fprintf(stdout, "Number of proc isn't a multiple of matrix size.\n");
    MPI_Finalize();
    return 0;
  }

  double *a=NULL,*b=NULL,*c=NULL;
  if(myrank==0){
    a = malloc(sizeof(*a)*n*n);
    b = malloc(sizeof(*b)*n*n);
    for(int i=0; i<n*n; i++){
      a[i] = rand()%100;
      b[i] = rand()%100;
    }
  }

  MPI_Comm proc_grid;
  mul_fox_grid_comm_create(&proc_grid);

  double smooth_time=0.0f;
  double max_time;

  for(int i=0; i<itr; i++){
    double s_time = MPI_Wtime();
    c = mul_fox(0, a, b, n, proc_grid);
    smooth_time += MPI_Wtime() - s_time;
    free(c);
    cache_buster();
  }
  smooth_time /= itr;
  MPI_Reduce(&smooth_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(myrank==0) fprintf(stdout, "%d %d %lf\n", size, n, max_time);

  MPI_Finalize();
  return 0;
}
