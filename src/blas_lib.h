#ifndef BLAS_LIB_H_
#define BLAS_LIB_H_

#include <stdio.h>
#include <stdlib.h>

#include <pthread.h>

#include "stack.h"

#define BLOCK_SIZE 16
#define DEFAULT_THREADS (24*12)

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

//Structure that represents the problems data to be sent to the threads
struct tasks_data_t{
  int tasks_nbr;
  pthread_mutex_t tasks_nbr_lock;
  stack_t *stack;
  pthread_mutex_t stack_lock;
};

//Structure representing a bloc of data to compute
struct dgemm_task_t{
  enum CBLAS_TRANSPOSE TransA;
  enum CBLAS_TRANSPOSE TransB;
  int M;
  int N;
  int K;
  double alpha;
  double *A;
  int A_inc;
  int lda;
  double *B;
  int B_inc;
  int ldb;
  double *C;
  int ldc;
};

//Thread function used to compute dgemm
void *dgemm_thread_launch(void *arg);

//Sequentiel dgemm
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
			  const int ldc);

//Multi-core dgemm computation
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

void cblas_daxpy(const int N,
	   const double alpha, const double *X, const int incX,
	   double *Y, const int incY);

double cblas_ddot(const int N, const double *X, const int incX, const double *Y, const int incY);

void cblas_dgemv(const enum CBLAS_ORDER order,
		 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
		 const double alpha, const double *A, const int lda,
		 const double *X, const int incX,
		 const double beta, double *Y, const int incY);


void cblas_dger(const enum CBLAS_ORDER order,
		const int M, const int N,
		const double alpha,
		const double *X, const int incX,
		const double *Y, const int incY,
		double *A, const int lda);

void cblas_dscal(const int N, const double alpha, double *X, const int incX);

void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);



#endif
