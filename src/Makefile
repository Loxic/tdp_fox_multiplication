CC=mpicc
CFLAGS=-std=c99 -Wall -g
LDFLAGS=-lm


mul_fox_tst:
	$(CC) $(CFLAGS) -o $@ mul_fox_tst.c mul_fox.c stack.c util.c blas_lib.c $(LDFLAGS)
