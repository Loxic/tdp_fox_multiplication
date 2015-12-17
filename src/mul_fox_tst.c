#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "mul_fox.h"

#define DEFAULT_OUTPUT "out.dat"
#define ROOT 0

int main(int argc, char**argv){
  MPI_Init(NULL, NULL);

  int myrank;
  int root_proc = ROOT;
  double *a = NULL, *b = NULL;
  double *c = NULL;
  int n;
  char *output_file = DEFAULT_OUTPUT;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(myrank==root_proc){
    if(argc>2){
      a = parse_matrix_file(argv[1], &n);
      b = parse_matrix_file(argv[2], &n);
    }
    else{
      fprintf(stdout, "No input filenames. Usage is %s [Matrix A] [Matrix B] [Output].\n", argv[0]);
      MPI_Finalize();
      return 0;
    }
    if(argc>3){
      output_file = argv[3];
    }
  }
  MPI_Bcast(&n, 1, MPI_INT, root_proc, MPI_COMM_WORLD);

  if(!is_nbr_proc_valid(n)){
    fprintf(stdout, "Number of proc isn't a multiple of matrix size.\n");
    MPI_Finalize();
    return 0;
  }

  MPI_Comm proc_grid;
  mul_fox_grid_comm_create(&proc_grid);

  c = mul_fox(root_proc, a, b, n, proc_grid);

  if(myrank==root_proc) affiche(n,n,c,n,stdout);
  if(myrank==root_proc) write_matrix_file(output_file, c, n);

  MPI_Finalize();
  return 0;
}
