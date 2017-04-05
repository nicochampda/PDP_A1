#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
  
    int rank, size, mat_size, block_size, i, j, p_x, dest;
    double A[mat_size][mat_size];
    double B[mat_size][mat_size];
    double C[mat_size][mat_size];
  
    MPI_Status status;

    MPI_Init(&argc, &argv);       
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    mat_size = atoi(argv[1]);
	/* Generate random Matrices */
    if (rank == 0){
	    for (i = 0; i < mat_size; i++){
			for (j = 0; i < mat_size; j++){
				A[i][j] = i;   /* to be Replace by random */
				B[i][j] = j;
			}
		}
		
		/* Distributes Blocks from A and B */
		p_x = sqrt(size);
		p_y = sqrt(size);
		
		block_size = mat_size / p_x;
		
		for (i = 0; i < p_x; i++){
			for (j = 0; j < p_y; j++){
				//MPI_Isend() A,B blocks
			}
		}
	}
	
	if (rank > 0){
		/* Collect blocks */
		//MPI_Recv () blocks 
	
		/* Apply fox algorithm on blocks */
	
	
		/* Send result */
		//MPI_send()
	}
	
	if (rank == 0){   //put that part in the first if ?
			//MPI_recv()
			//sum up results
	}
  MPI_Finalize();
  
  return 0;
}
