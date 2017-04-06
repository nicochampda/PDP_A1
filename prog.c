#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define SIZE 8
#define PROC 4

int main(int argc, char **argv) {

	MPI_Init(&argc, &argv);
	int p, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	char i;

	char a[SIZE*SIZE];
	const int BLOCKS = SIZE/sqrt(PROC);

	if (rank == 0) {
		for (int j=0; j<SIZE*SIZE; j++) {
			a[j] = (char)j;
		}
	}
	if (p != PROC) {
		printf("Error");
		MPI_Finalize();
		exit(-1);
	}
	char b[BLOCKS*BLOCKS];
	for (int j=0; j<BLOCKS*BLOCKS; j++){
		b[j]=0;
	}
	MPI_Datatype blocktype;
	MPI_Datatype blocktype2;
	MPI_Type_vector(BLOCKS,BLOCKS,SIZE,MPI_CHAR,&blocktype2);
	MPI_Type_create_resized(blocktype2,0,sizeof(char),&blocktype);
	MPI_Type_commit(&blocktype);

	int disps[PROC];
	int counts[PROC];
	for (int j=0; j<2; j++) {
		for (int k=0; k<2; k++) {
			disps[j*2+k] = j*SIZE*BLOCKS+k*BLOCKS;
			counts [j*2+k] = 1;
		}
	}
	MPI_Scatterv(a, counts, disps, blocktype, b, BLOCKS*BLOCKS,MPI_CHAR,0,MPI_COMM_WORLD);
	for (int proc=0;proc<p;proc++){
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
			}
			printf("Local Matrix:\n");
			for (int j=0; j<BLOCKS; j++){
				for (int k=0; k<BLOCKS; k++){
					printf("%3d",(int)b[j*BLOCKS+k]);
				}
				printf("\n");
			}
			printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}

