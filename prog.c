#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

//#define SIZE 8
//#define PROC 4

int main(int argc, char *argv[]) {
	
	if (argc != 2){
		printf("Wrong number of arguments given");
		return -1;
	}
	const int SIZE = atoi(argv[1]);
	MPI_Init(&argc, &argv);
	int p, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	char i;

	char a[SIZE*SIZE];
	char aa[SIZE*SIZE];
	const int PROC = p;
	const int SPROC = sqrt(PROC);
	const int BLOCKS = SIZE/SPROC;
	if (SIZE%SPROC!=0){
		printf("SIZE has to be divisible by SPROC");
		MPI_Finalize();
		exit(-1);
	}

	if (rank == 0) {
		for (int j=0; j<SIZE*SIZE; j++) {
			a[j] = (char)j;
			int y= SIZE*SIZE-j;
			aa[j] = (char)y;
		}
	}
	char b[BLOCKS*BLOCKS];
	char bb[BLOCKS*BLOCKS];
	for (int j=0; j<BLOCKS*BLOCKS; j++){
		b[j]=0;
		bb[j]=0;
	}
	MPI_Datatype blocktype;
	MPI_Datatype blocktype2;
	MPI_Type_vector(BLOCKS,BLOCKS,SIZE,MPI_CHAR,&blocktype2);
	MPI_Type_create_resized(blocktype2,0,sizeof(char),&blocktype);
	MPI_Type_commit(&blocktype);

	int disps[PROC];
	int counts[PROC];
	for (int j=0; j<SPROC; j++) {
		for (int k=0; k<SPROC; k++) {
			disps[j*SPROC+k] = j*SIZE*BLOCKS+k*BLOCKS;
			counts [j*SPROC+k] = 1;
		}
	}
	MPI_Scatterv(a, counts, disps, blocktype, b, BLOCKS*BLOCKS,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Scatterv(aa, counts, disps, blocktype, bb, BLOCKS*BLOCKS,MPI_CHAR,0,MPI_COMM_WORLD);
	for (int proc=0;proc<PROC;proc++){
		if (proc == rank) {
			printf("Rank = %d\n", rank);
			if (rank == 0) {
				printf("Global matrix: \n");
				for (int k=0;k<SIZE;k++){
					for (int j=0; j<SIZE; j++){
						printf("%3d",(int)a[k*SIZE+j]);
					}
					printf("\n");
				}
				printf("\n");
				for (int k=0;k<SIZE;k++){
					for(int j=0; j<SIZE; j++){	
						printf("%3d",(int)aa[k*SIZE+j]);
					}
					printf("\n");
				}
			}
			printf("Local Matrix:\n");
			for (int j=0; j<BLOCKS; j++){
				for (int k=0; k<BLOCKS; k++){
					printf("%3d",(int)b[j*BLOCKS+k]);
				}
				printf("\n");
			}
			printf("\n");
			for (int j=0; j<BLOCKS; j++){
				for (int k=0; k<BLOCKS; k++){
					printf("%3d",(int)bb[j*BLOCKS+k]);
				}
			
			printf("\n");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}

