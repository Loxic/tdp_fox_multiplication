#ifndef _MUL_FOX_H_
#define _MUL_FOX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "blas_lib.h"

void mul_fox_grid_comm_create(MPI_Comm *grid);
double * mul_fox(int root_proc, const double * a, const double * b, int n, const MPI_Comm proc_grid);

#endif
