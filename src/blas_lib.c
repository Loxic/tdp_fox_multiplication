#include "blas_lib.h"

double cblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY){
  double r = 0.0f;
  for(int i=0; i<N/4; i++){
    r += X[i*incX]+Y[i*incY];
    r += X[(i+1)*incX]+Y[(i+1)*incY];
    r += X[(i+2)*incX]+Y[(i+2)*incY];
    r += X[(i+3)*incX]+Y[(i+3)*incY];
  }
  for(int i = 0; i < N%4; i++) {
    r += X[(N-1-i)*incX]+Y[(N-1-i)*incX];
  }
  return r;
}

void cblas_daxpy(const int N,
	   const double alpha, const double *X, const int incX,
	   double *Y, const int incY) {
  for(int i=0; i < N/4; i+=4) {
    Y[i*incX] += alpha*X[i*incY];
    Y[(i+1)*incX] += alpha*X[(i+1)*incY];
    Y[(i+2)*incX] += alpha*X[(i+2)*incY];
    Y[(i+3)*incX] += alpha*X[(i+3)*incY];
  }
  for(int i=0 ; i < N%4; i++) {
    Y[(N-1-i)*incX] += alpha*X[(N-1-i)*incY];
  }
}

void cblas_dgemv(const enum CBLAS_ORDER order,
		 const enum CBLAS_TRANSPOSE Trans, const int M, const int N,
		 const double alpha, const double *A, const int lda,
		 const double *X, const int incX,
		 const double beta, double *Y, const int incY) {	
  int transStepA = 1;
  int ldaBis = lda;
  if(Trans == CblasTrans) {
    transStepA = lda;
    ldaBis = 1;
  }

  for(int i=0; i < N; i++){
    Y[i*incY] *= beta;
  }
  for(int i=0; i < N; i++) {
    for(int j=0; j < M; j++) {
      Y[i*incY] += alpha*A[i*transStepA+j*ldaBis]*X[j*incX];
    }
  }
}


void cblas_dger(const enum CBLAS_ORDER order,
		const int M, const int N,
		const double alpha,
		const double *X, const int incX,
		const double *Y, const int incY,
		double *A, const int lda) {	
  for(int i=0; i < M; i++) {
    for(int j=0; j < N; j++) {
      A[i+j*lda] += alpha*X[i*incX]*Y[j*incY];
    }
  }
}

void cblas_dgemm_scalaire(/*const enum CBLAS_ORDER Order, */
			  const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB,
			  const int M,
			  const int N,
			  const int K,
			  const double alpha,
			  const double *A,
			  const int lda,
			  const double *B,
			  const int ldb,
			  double *C,
			  const int ldc){
  int transStepA = 1;
  int transStepB = 1;
  int ldaBis = lda;
  int ldbBis = ldb;
  if(TransA == CblasTrans ||
     TransA == CblasConjTrans){
    transStepA = lda;
    ldaBis = 1;
  }
  if(TransB == CblasTrans ||
     TransB == CblasConjTrans){
    transStepB = ldb;
    ldbBis = 1;
  }
  if(((TransA == CblasNoTrans && TransB != CblasNoTrans)
      || (TransA != CblasNoTrans && TransB == CblasNoTrans))
     && (M != (N != K)))
    return ;

  for(int i=0; i<M; i++){
    for(int j=0; j<K; j++){
      for(int k=0; k<N; k++){
	C[i+j*ldc] += alpha * B[k*transStepB+j*ldbBis] * A[i*transStepA+k*ldaBis];
      }
    }
  }  
}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){
  char *env_var = getenv("MYLIB_NUM_THREADS");
  int n_threads = DEFAULT_THREADS;
  if(!env_var){
    //fprintf(stderr, "Unable to access environment variable MYLIB_NUM_THREADS.\n");
    //exit(-1);
    //fprintf(stderr, "[WARNING]: couldn't charge environment variable MYLIB_NUM_THREADS, continuing with %d threads.\n", DEFAULT_THREADS);
  }else{
    n_threads = atoi(env_var);
  }

  stack_t *task_stack = stack__create();
  int M_nbr = M/BLOCK_SIZE;
  if(M_nbr==0)
    M_nbr++;
  int N_nbr = N/BLOCK_SIZE;
  if(N_nbr==0)
    N_nbr++;
  int K_nbr = K/BLOCK_SIZE;
  if(K_nbr==0)
    K_nbr++;
  int stepA=1;
  int stepB=1;
  int ldaBis=lda;
  int ldbBis=ldb;

  pthread_t *threads = malloc(sizeof(*threads)*n_threads);
  struct tasks_data_t tasks_data;
  tasks_data.tasks_nbr = M_nbr*K_nbr;
  pthread_mutex_init(&(tasks_data.tasks_nbr_lock), NULL);
  tasks_data.stack = task_stack;
  pthread_mutex_init(&(tasks_data.stack_lock), NULL);

  if(TransA == CblasTrans ||
     TransA == CblasConjTrans){
    stepA = ldaBis;
    ldaBis = 1;
  }
  if(TransB == CblasTrans ||
     TransB == CblasConjTrans){
    stepB = ldbBis;
    ldbBis = 1;
  }
  if(((TransA == CblasNoTrans && TransB != CblasNoTrans)
      || (TransA != CblasNoTrans && TransB == CblasNoTrans))
     && (M != (N != K)))
    return ;
    
  if(beta != 1.0f)
    for(int i=0; i<M*K; i++)
      C[i]*=beta;

  for(int i=0; i<n_threads; i++){
    pthread_create(threads+i, NULL, dgemm_thread_launch, (void *)(&tasks_data));
  }

  for(int i=0; i<M_nbr; i++){
    for(int j=0; j<K_nbr; j++){
      struct dgemm_task_t *task = malloc(sizeof(*task));

      int M_size=BLOCK_SIZE;
      int K_size=BLOCK_SIZE;
      if(j==K_nbr-1){
	if(K/BLOCK_SIZE==0)
	  K_size = K;
	else
	  K_size = BLOCK_SIZE+K%BLOCK_SIZE;
      }
      if(i==M_nbr-1){
	if(M/BLOCK_SIZE==0)
	  M_size = M;
	else
	  M_size = BLOCK_SIZE+M%BLOCK_SIZE;
      }

      task->TransA = TransA;
      task->TransB = TransB;
      task->M      = M_size;
      task->N      = N;
      task->K      = K_size;
      task->alpha  = alpha;
      task->A      = (A+i*BLOCK_SIZE*stepA);
      task->A_inc  = ldaBis*BLOCK_SIZE;
      task->lda    = lda;
      task->B      = (B+j*ldbBis*BLOCK_SIZE);
      task->B_inc  = BLOCK_SIZE*stepB;
      task->ldb    = ldb;
      task->C      = (C+i*BLOCK_SIZE+j*ldc*BLOCK_SIZE);
      task->ldc    = ldc;

      pthread_mutex_lock(&(tasks_data.stack_lock));
      stack__push(task, tasks_data.stack);
      pthread_mutex_unlock(&(tasks_data.stack_lock));
    }
  }
  for(int i=0; i<n_threads; i++){
    pthread_join(threads[i], NULL);
  }
  pthread_mutex_destroy(&(tasks_data.tasks_nbr_lock));
  pthread_mutex_destroy(&(tasks_data.stack_lock));
  stack__destroy(task_stack);
  free(threads);
}

void * dgemm_thread_launch(void *arg){
  struct tasks_data_t *tasks_data = (struct tasks_data_t*) arg;

  while(1){
    pthread_mutex_lock(&(tasks_data->tasks_nbr_lock));
    if(tasks_data->tasks_nbr == 0){
      pthread_mutex_unlock(&(tasks_data->tasks_nbr_lock));
      break;
    }
    pthread_mutex_unlock(&(tasks_data->tasks_nbr_lock));

    pthread_mutex_lock(&(tasks_data->stack_lock));//        stack lock
    if(!stack__is_empty(tasks_data->stack)){
      struct dgemm_task_t *data = (struct dgemm_task_t *)stack__pop(tasks_data->stack);
      pthread_mutex_unlock(&(tasks_data->stack_lock));//    stack unlock
      pthread_mutex_lock(&(tasks_data->tasks_nbr_lock));//  nbr lock
      tasks_data->tasks_nbr--;
      pthread_mutex_unlock(&(tasks_data->tasks_nbr_lock));//nbr unlock

      int N_nbr = data->N/BLOCK_SIZE;
      if(N_nbr==0)
	N_nbr++;
      int N_size = BLOCK_SIZE;
      for(int k=0; k<N_nbr; k++){
	if(k==N_nbr-1){
	  if(data->N/BLOCK_SIZE==0)
	    N_size = data->N;
	  else
	    N_size = BLOCK_SIZE+(data->N%BLOCK_SIZE);
	}
	
	cblas_dgemm_scalaire(data->TransA,
			     data->TransB,
			     data->M,
			     N_size,
			     data->K,
			     data->alpha,
			     (data->A+(data->A_inc*k)), data->lda,
			     (data->B+(data->B_inc*k)), data->ldb,
			     data->C, data->ldc);
      }
      free(data);
      continue;
    }
    pthread_mutex_unlock(&(tasks_data->stack_lock));//   stack unlock
    pthread_yield();
  }
  return NULL;
}

void cblas_dscal(const int N, const double alpha, double *X, const int incX) {
  for(int i = 0; i < N; i++) {
    X[i*incX] *= alpha;
  }
}


void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb) {
  //Side == L, TransA == Nope

  // B = alpha*B
  if(alpha != 1) {
    for(int i = 0 ; i < M; i++) {
      for(int j = 0 ; j < N; j++) {
	B[i+ldb*j] *= alpha;
      }
    }
  }

  
  // Matrix U (Diag == CblasNonUnit)
  if(Uplo == CblasUpper) {
    if(Diag != CblasNonUnit) {
      fprintf(stderr,"CblasUnit is not implemented for upper matrix");
      return;
    }
    for(int j = 0 ; j < N; j++) {
      for(int k = M-1; k >= 0; k--) {
	if(B[k+ldb*j]!=0) {
	  B[k+ldb*j] = B[k+ldb*j]/A[k+lda*k];
	  for(int i=0; i<k; i++)
	    B[i+j*ldb] -= B[k+ldb*j]*A[i+k*lda];
	}
      }
    }
  }
  //Matrix L (Diag == CblasUnit)
  else if(Uplo == CblasLower) {
    if(Diag != CblasUnit) {
      fprintf(stderr,"CblasNonUnit is not implemented for lower matrix");
      return;
    }
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < M; k++) {
	if(B[k+ldb*j]!=0) {
	  //B[k+ldb*j] = B[k+ldb*j]/A[k+lda*k];
	  for(int i = k+1; i < M; i++) {
	    B[i+ldb*j] -= B[k+ldb*j]*A[i+lda*k];
	  }
	}
      }
    }
  }
}
