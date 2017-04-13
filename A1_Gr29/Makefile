CC         =  mpicc
CCFLAGS    =  -O3
LIBS       =  -lmpi -lm

matmult: matmult.c
	$(CC) $(CCFLAGS) -o matmult matmult.c $(LIBS)

clean:
	$(RM) matmult
