CC         =  mpicc
CCFLAGS    =  -O3 -g
LIBS       =  -lmpi -lm

matmult: matmult.c
	$(CC) $(CCFLAGS) -o matmult matmult.c $(LIBS)

test: test.c
	gcc $(CCFLAGS) -o test test.c

integral2d:             prog.c
	$(CC) $(CCFLAGS) -o prog prog.c $(LIBS)


clean:
	$(RM) prog
	$(RM) matmult
	$(RM) prog
