CC         =  mpicc
CCFLAGS    =  -O3 -Wall
LIBS       =  -lmpi -lm



integral2d:             prog.c
	$(CC) $(CCFLAGS) -o prog prog.c $(LIBS)


clean:
	$(RM) prog
