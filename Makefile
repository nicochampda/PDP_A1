CC         =  mpicc
CCFLAGS    =  -O3
LIBS       =  -lmpi -lm



integral2d:             prog.c
	$(CC) $(CCFLAGS) -o prog prog.c $(LIBS)


clean:
	$(RM) prog
