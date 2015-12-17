#include "mul_fox.h"

#define RANK_DEBUG 16

int is_nbr_proc_valid(int n){
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size == ((int)sqrt(size))*((int)sqrt(size))
    && n%(int)sqrt(size) == 0;
}

void mul_fox_grid_comm_create(MPI_Comm *grid){
  int dim[2];
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int n = sqrt(size);
  dim[0]=n;
  dim[1]=n;
  int period[2]={1,1};
  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, 0, grid);
}

void ab_scatter(int root_proc, const double *a, const double *b, double *a_local, double *b_local, int n, const MPI_Comm proc_grid){
  int size;
  MPI_Comm_size(proc_grid, &size);
  int m = n/sqrt(size);
  int *sendcounts, *displs;
  sendcounts = malloc(sizeof(*sendcounts)*size);
  displs = malloc(sizeof(*displs)*size);
  for(int i=0; i<size; i++)
    sendcounts[i]=1;
  for(int i=0; i<size; i++){
    displs[i]=m*(i%(int)sqrt(size)) 
      + (i/((int)sqrt(size)))*m*n;
  }
  MPI_Datatype scatter_type2, scatter_type;
  MPI_Type_vector(m, m, n, MPI_DOUBLE, &scatter_type2);
  MPI_Type_create_resized(scatter_type2, 0, sizeof(double), &scatter_type);
  MPI_Type_commit(&scatter_type);

  MPI_Scatterv((void *)a, sendcounts, displs, scatter_type, a_local, m*m, MPI_DOUBLE, root_proc, proc_grid);
  MPI_Scatterv((void *)b, sendcounts, displs, scatter_type, b_local, m*m, MPI_DOUBLE, root_proc, proc_grid);

  free(displs);
  free(sendcounts);
  MPI_Type_free(&scatter_type);
  MPI_Type_free(&scatter_type2);
}

void c_gather(int root_proc, double *c, double *c_local, int n, const MPI_Comm proc_grid){
  int size;
  MPI_Comm_size(proc_grid, &size);
  int m = n/sqrt(size);
  int *recvcounts, *displs;
  recvcounts = malloc(sizeof(*recvcounts)*size);
  displs = malloc(sizeof(*displs)*size);
  for(int i=0; i<size; i++)
    recvcounts[i]=1;
  for(int i=0; i<size; i++){
    displs[i]=m*(i%(int)sqrt(size)) 
      + (i/((int)sqrt(size)))*m*n;
  }
  MPI_Datatype gather_type2, gather_type;
  MPI_Type_vector(m, m, n, MPI_DOUBLE, &gather_type2);
  MPI_Type_create_resized(gather_type2, 0, sizeof(double), &gather_type);
  MPI_Type_commit(&gather_type);

  MPI_Gatherv(c_local, m*m, MPI_DOUBLE, c, recvcounts, displs, gather_type, root_proc, proc_grid);

  free(displs);
  free(recvcounts);
  MPI_Type_free(&gather_type);
  MPI_Type_free(&gather_type2);
}

double * mul_fox(int root_proc, const double * a, const double * b, int n, const MPI_Comm proc_grid){
  int size, myrank;
  MPI_Comm_rank(proc_grid, &myrank);
  MPI_Comm_size(proc_grid, &size);
  int mycoord[2];
  MPI_Cart_coords(proc_grid, myrank, 2, mycoord);
  int dimB_remain[2]={0,1};
  int dimA_remain[2]={1,0};
  MPI_Comm b_comm, a_comm;
  MPI_Status status;
  MPI_Cart_sub(proc_grid, dimB_remain, &b_comm);
  MPI_Cart_sub(proc_grid, dimA_remain, &a_comm);
  int proc_n = sqrt(size);
  int m = n/proc_n;

  double *a_local = malloc(sizeof(*a_local)*m*m);
  double *b_local = malloc(sizeof(*b_local)*m*m);

  double *c = NULL;
  if(myrank==root_proc)
    c = malloc(sizeof(*c)*n*n);

  ab_scatter(root_proc, a, b, a_local, b_local, n, proc_grid);

  double *a_buffer = malloc(sizeof(*a_buffer)*m*m);
  for(int i=0; i<m*m; i++)
    a_buffer[i] = a_local[i];

  double *c_local = malloc(sizeof(*c_local)*m*m);
  for(int i=0; i<m*m; i++)
    c_local[i] = 0.0f;

  for(int k=0; k<proc_n; k++){
    int a_coord[1];
    int a_rank;
    int b_from, b_to;
    int a_myrank;
    double *tmp_ptr;

    a_coord[0]=(mycoord[1]+k)%proc_n;
    MPI_Cart_rank(a_comm, a_coord, &a_rank);
    //MPI_Cart_rank(a_comm, mycoord, &a_myrank);
    MPI_Comm_rank(a_comm, &a_myrank);

    if(a_myrank==a_rank)
      tmp_ptr=a_local;
    else
      tmp_ptr=a_buffer;
    MPI_Bcast(tmp_ptr, m*m, MPI_DOUBLE, a_rank, a_comm);
    //if(myrank==RANK_DEBUG) fprintf(stderr, "Information: k=%d | a_coord[0]=%d | a_rank=%d | a_myrank=%d | mycoord[1]=%d\n", k, a_coord[0], a_rank, a_myrank, mycoord[1]);
    //if(myrank==RANK_DEBUG) affiche(m,m,tmp_ptr,m,stdout);

    if(m > 10)
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0f, tmp_ptr, m, b_local, m, 1.0f, c_local, m);
    else
      cblas_dgemm_scalaire(CblasNoTrans, CblasNoTrans, m, m, m, 1.0f, tmp_ptr, m, b_local, m, c_local, m);

    MPI_Cart_shift(proc_grid, 1, -1, &b_from, &b_to);
    MPI_Sendrecv_replace(b_local, m*m, MPI_DOUBLE, b_to, 42, b_from, 42, proc_grid, &status);

  }

  free(a_buffer);
  free(a_local);
  free(b_local);

  c_gather(root_proc, c, c_local, n, proc_grid);
  free(c_local);
  return c;
}
