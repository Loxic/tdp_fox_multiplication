#include "mul_fox.h"

void mul_fox_grid_comm_create(MPI_Comm *grid){
  int dim[2];
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int n = sqrt(size);
  dim[0]=n;
  dim[1]=n;
  int period[2]={1,0};
  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, 1, grid);
}

void ab_scatter(int root_proc, const double *a, const double *b, double *a_local, double *b_local, int n, const MPI_Comm proc_grid){
  int myrank;
  MPI_Comm_rank(proc_grid, &myrank);

  int size;
  MPI_Comm_size(proc_grid, &size);
  int m = n/sqrt(size);
  size = sqrt(size)*sqrt(size);
  int *sendcounts, *displs;
  sendcounts = malloc(sizeof(*sendcounts)*size);
  displs = malloc(sizeof(*displs)*size);
  for(int i=0; i<size; i++)
    sendcounts[i]=1;
  for(int i=0; i<size; i++){
    int col_displ = i/sqrt(size);
    displs[i]=m*i + col_displ*m*n;
  }

  MPI_Datatype scatter_type;
  MPI_Type_vector(m, m, n, MPI_DOUBLE, &scatter_type);
  MPI_Type_commit(&scatter_type);
  fprintf(stderr, "lol %d\n", myrank);
  MPI_Scatterv((void *)a, sendcounts, displs, scatter_type, (void *)a_local, m*m, MPI_DOUBLE, root_proc, proc_grid);
  fprintf(stderr, "lol2 %d\n", myrank);
  MPI_Scatterv((void *)b, sendcounts, displs, scatter_type, b_local, m*m, MPI_DOUBLE, root_proc, proc_grid);

  MPI_Type_free(&scatter_type);

}

double * mul_fox(int root_proc, const double * a, const double * b, int n, const MPI_Comm proc_grid){
  int size, myrank;
  MPI_Comm_rank(proc_grid, &myrank);
  MPI_Comm_size(proc_grid, &size);
  int mycoord[2];
  MPI_Cart_coords(proc_grid, myrank, 2, mycoord);
  int dimB_remain[2]={1,0};
  int dimA_remain[2]={0,1};
  MPI_Comm b_comm, a_comm;
  MPI_Status status;
  MPI_Cart_sub(proc_grid, dimB_remain, &b_comm);
  MPI_Cart_sub(proc_grid, dimA_remain, &a_comm);
  int m = n/sqrt(size);

  double *a_local = malloc(sizeof(*a_local)*m*m);
  double *b_local = malloc(sizeof(*a_local)*m*m);

  ab_scatter(root_proc, a, b, a_local, b_local, n, proc_grid);

  double *a_buffer = malloc(sizeof(*a_buffer)*m*m);
  for(int i=0; i<m*m; i++)
    a_buffer[i] = a_local[i];

  double *c_local = malloc(sizeof(*c_local)*m*m);
  for(int i=0; i<m*m; i++)
    c_local[i] = 0.0f;

  for(int k=0; k<m; k++){
    int a_coord[1]={(mycoord[0]+k)%(int)sqrt(size)};
    int a_rank;
    int b_from, b_to;
    int a_myrank;
    double *tmp_ptr;

    MPI_Cart_rank(a_comm, a_coord, &a_rank);
    MPI_Cart_rank(a_comm, mycoord, &a_myrank);

    //MPI_Comm_rank(a_comm, &a_myrank);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0f, a_buffer, m, b_local, m, 1.0f, c_local, m);

    if(a_myrank==a_rank)
      tmp_ptr=a_local;
    else
      tmp_ptr=a_buffer;

    MPI_Bcast(tmp_ptr, m*m, MPI_DOUBLE, a_rank, a_comm);
    MPI_Cart_shift(proc_grid, 0, -1, &b_from, &b_to);
    MPI_Sendrecv_replace(b_local, m*m, MPI_DOUBLE, b_to, 42, b_from, 42, proc_grid, &status);

  }

  free(a_buffer);
  free(a_local);
  free(b_local);
  return c_local;
}
