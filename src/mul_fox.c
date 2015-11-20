#include "mul_fox.h"

void mul_fox_grid_comm_create(MPI_Comm *grid){
  int dim[2];
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int n = sqrt(size);
  dim[0]=n;
  dim[1]=n;
  int period[2]={1,0}
  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, 1, grid);
}

double * mul_fox(const double * a_local, const double * b_local, int n, const MPI_Comm proc_grid){
  int size, myrank;
  int mycoord[2];
  MPI_Cart_coords(proc_grid, myrank, 2, mycoord);
  int dimB_remain[2]={1,0};
  int dimA_remain[2]={0,1};
  MPI_Comm b_comm, a_comm;
  MPI_Status status;
  MPI_Cart_sub(proc_grid, dimB_remain, &b_comm);
  MPI_Cart_sub(proc_grid, dimA_remain, &a_comm);
  MPI_Comm_rank(proc_grid, &myrank);
  MPI_Comm_size(proc_grid, &size);
  int n_proc = sqrt(size);
  int m = n/n_proc;

  double *a_buffer = malloc(sizeof(*a_buffer)*m*m);
  for(int i=0; i<m*m; i++)
    a_buffer[i] = a_local[i];

  double *b_buffer = malloc(sizeof(*b_buffer)*m*m);
  for(int i=0; i<m*m; i++)
    b_buffer[i] = b_local[i];

  double *c_local = malloc(sizeof(*c_local)*m*m);
  for(int i=0; i<m*m; i++)
    c_local[i] = 0.0f;

  for(int k=0; k<m; k++){
    int a_coord[1]={(mycoord[0]+k)%n_proc};
    int a_rank;
    int b_from, b_to;
    int a_myrank;
    double *tmp_ptr;
    MPI_Cart_rank(a_comm, a_coord, &a_rank);
    MPI_Comm_rank(a_comm, &a_myrank);
    cblas_dgemm(order, transa, transb, m, m, m, 1.0f, a_buffer, m, b_buffer, m, 1.0f, c_local, m);
    if(a_myrank==a_rank)
      tmp_ptr=a_local;
    else
      tmp_ptr=a_buffer;
    MPI_Bcast(tmp_ptr, m*m, MPI_DOUBLE, a_rank, a_comm);
    MPI_Cart_shift(proc_grid, 0, -1, &b_from, &b_to);
    MPI_Sendrecv_replace(b_buffer, m*m, MPI_DOUBLE, b_to, 42, b_from, 42, proc_grid, &status);
  }

  free(a_buffer);
  free(b_buffer);
  return c_local;
}
