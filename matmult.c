#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Function to perform the step 2 of Fox's algorithm
void Block_matmul(int *subA, int *SubB, int *SubC, int block_size){
     int i,j,k;
     for(i=0;i<block_size;i++){
        for(k=0;k<block_size;k++){
           for(j=0;j<block_size;j++){ // efficient matrix multiplication
              subC[i][j]+=subA[i][k]*subB[k][j];
           } 
        }
     }     
}





int main(int argc, char *argv[]) {
  
    int rank, size, mat_size, block_size, i, j,k, p_x,p_y, dest, offset;
    double A[mat_size][mat_size];
    double B[mat_size][mat_size];
    double C[mat_size][mat_size];
    double *subA;
    double *subB;
    double *subC;
    int dims[2] = {0,0}, periods[2] = {1,1}, reorder = 1; 
    
    // create communicators
    MPI_Comm row_comm;
    MPI_Comm col_comm;
    MPI_Comm block_comm;

  
    MPI_Datatype block;
    MPI_Status status;


    MPI_Init(&argc, &argv);       
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    
     /* Create a virtual 2D-grid topology */
    MPI_Dims_create(size, 2, dims);  
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &block_comm); 
    MPI_Comm_rank(block_comm, &rank); 

    //MPI_Comm_split(MPI_COMM_WORLD, rank / 4, 0, &row_comm);
    //MPI_Comm_split(MPI_COMM_WORLD, rank % 4, 0, &col_comm);
    //get row and column rank
    //
    
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
	    offset = 0;
	
            for (i = 0; i < p_x; i++){
		for (j = 0; j < p_y; j++){
                    //for matrix A
		    MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
		    MPI_Send(&block_size, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
       		    MPI_Send(&a[offset][0], block_size*mat_size, MPI_DOUBLE, i+j, 1, MPI_COMM_WORLD);

                    //for matrix B

      		}
       		offset = offset + block_size;
       	    }
		
	}
	
	if (rank > 0){
		/* Collect blocks */
		//MPI_Recv () blocks 
	
		/* Apply fox algorithm on blocks */
              
              //1. broadcast Aim
              //2. call block_matmul
              //3. shift the B blocks


               fox();
	
	
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
