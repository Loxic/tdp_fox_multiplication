CC=mpicc
CFLAGS=-std=c99 -Wall -g
LDFLAGS=-lm
REQ=mul_fox.o stack.o util.o blas_lib.o
TARGETS=mul_fox_tst mul_fox_correctness mul_fox_perf

.PHONY: clean

all: $(TARGETS)

mul_fox_correctness: mul_fox_correctness.o $(REQ)
mul_fox_tst: mul_fox_tst.o $(REQ)
mul_fox_perf: mul_fox_perf.o $(REQ)

$(TARGETS):
	$(CC) -o $@ $^ $(LDFLAGS)

%.o:%.c %.h
	$(CC) $(CFLAGS) -c $< 

clean:
	rm *.o *~ $(TARGETS)
