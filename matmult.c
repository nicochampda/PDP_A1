#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>


// Function to perform the step 2 of Fox's algorithm
void Block_matmul(double **subA, double **subB, double **subC, int block_size){
     int i,j,k;
     for(i = 0; i < block_size; i++){
          for(k = 0; k < block_size; k++){
              for(j = 0; j < block_size; j++){ // efficient matrix multiplication
                  subC[i][j] += subA[i][k] * subB[k][j];
              } 
          }
     }     
}

void PrMat(double **matrix, int row, int col){   
    int i, j;
    printf("\n");
    for (i = 0; i < row; i++){
        for (j = 0; j < col; j++){
            printf("%f", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}





int main(int argc, char *argv[]) {
  
    int rank, size, mat_size, block_size, i, j,k, p_x,p_y, dest, offset;
    double **subA;
    double **subB;
    double **subC;
    double B_blocks[block_size][block_size];
    double *arows;

    int dims[2] = {0,0}, periods[2] = {1,1}, reorder = 1; 
    
    // create communicators
    /*MPI_Comm row_comm;
    MPI_Comm col_comm;
    MPI_Comm block_comm;    
    MPI_Datatype block;*/


    mat_size = atoi(argv[1]);
    double A[mat_size][mat_size];
    double B[mat_size][mat_size];
    double C[mat_size][mat_size];

    MPI_Status status;
    MPI_Status status2;
    MPI_Status status3;

    MPI_Request request1;
    MPI_Request request2;
    MPI_Request request3;

    MPI_Init(&argc, &argv);       
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    
     /* Create a virtual 2D-grid topology */
   /* MPI_Dims_create(size, 2, dims);  
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &block_comm); 
    MPI_Comm_rank(block_comm, &rank); 
*/
    //MPI_Comm_split(MPI_COMM_WORLD, rank / 4, 0, &row_comm);
    //MPI_Comm_split(MPI_COMM_WORLD, rank % 4, 0, &col_comm);
    //get row and column rank
    //
    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;
    MPI_Type_vector(block_size,block_size,mat_size,MPI_DOUBLE,&blocktype2);
    MPI_Type_create_resized(blocktype2,0,sizeof(double),&blocktype);
    MPI_Type_commit(&blocktype);

    p_x = sqrt(size);
    p_y = sqrt(size);

    block_size = mat_size / p_x;
    
	/* Generate random Matrices */
    if (rank == 0){

        for (i = 0; i < mat_size; i++){
	    for (j = 0; j < mat_size; j++){
       	        A[i][j] = i+j;   /* to be Replace by random */
	        B[i][j] = 2.0;
	    }
        }
        for (i = 0; i < mat_size; i++){
	   printf("\n");
		 for (j = 0; j < mat_size; j++){
		    printf("%.1f\t", A[i][j]);
		}
	}
	   printf("\n");

	    /* Distributes Blocks from A and B */
	    
	    offset = 0;
	
            for (i = 0; i < p_x; i++){
		for (j = 0; j < p_y; j++){
                    //for matrix A
                    if ((i+j) != 0){
		    	MPI_Isend(&mat_size, 1, MPI_INT, i*p_x + j, 1, MPI_COMM_WORLD, &request1);
		    	MPI_Isend(&block_size, 1, MPI_INT, i*p_x + j, 2, MPI_COMM_WORLD, &request2);
			for (k = 0; k < block_size;k++){
       		    		MPI_Isend(&A[offset + k][0], mat_size, MPI_DOUBLE, i*p_x + j, 100+k, MPI_COMM_WORLD, &request3);
		}
		}
                    //for matrix B

      		}
       		offset = offset + block_size;
       	    }
		
	}
	int disps[p_x*p_y];
        int counts[p_x*p_y];

        for (i = 0; i < p_x; i++){
                for (j = 0; j < p_y; j++){
                    disps[i*p_x+j]=i*mat_size*block_size+j*block_size;
                    counts[i*p_x+j]=1;       
                 }
        }


        MPI_Scatterv(B,counts,disps,blocktype,B_blocks,block_size*block_size,MPI_DOUBLE,0,MPI_COMM_WORLD);


	if (rank > 0){
		/* Collect A rows */
               MPI_Recv(&mat_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);  
	       MPI_Recv(&block_size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status2);  
	       double arows[mat_size][block_size];
		for (k = 0; k < block_size;k++){
			MPI_Recv(&arows[k], mat_size, MPI_DOUBLE, 0, 100+k, MPI_COMM_WORLD, &status3);  
		}

		sleep(rank);
		printf("mat %i block %i rank %i\n",mat_size, block_size, rank);
		for(i = 0; i < block_size; i++){
			printf("\n");
			for(j = 0; j < mat_size; j++){
				printf("%.1f\t", arows[i][j]);
			}
		} 
	
		printf("\n");
		/* Collect B blocks */
		

		/* Apply fox algorithm on blocks */
              
              //1. broadcast Aim
              //2. call block_matmul
              //3. shift the B blocks


               //fox();
	
	
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
