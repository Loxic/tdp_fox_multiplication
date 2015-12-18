To compile binaries :
   $> cd src/
   $> make

Binaries and usage:
   mul_fox_tst
      Tests the fox multiplication algorithm. Takes two .dat files and produces a out.dat file or other if specified. The number of nodes assigned to this program, as with all others, must be compatible with the size of the matrixes. (ex: 10x10 matrixes -> 1,4,25,100 nodes accepted)
      $> mpirun mul_fox_tst [A.dat filename] [B.dat filename] ([C.dat filename])
   
   mul_fox_correctness
      Validates the correctness of the implementation by comparing the fox product of randomly filled pair of matrixes, with a classic dgemm matrix product on the same matrixes. The size of the matrixes compared can be parametered.
      $> mpirun mul_fox_correctness ([n size of matrix])

   mul_fox_perf
      Times the computation of the fox algorithm product on two randomly filled matrixes.
      $> mpirun mul_fox_perf [matrix size] [smoothing iterations]
