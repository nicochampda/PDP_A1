#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>


// Function to perform the step 2 of Fox's algorithm
void Block_matmul(double **subA, double **subB, double **subC, int block_size, int block_nbr){
    int i,j,k;
    for(i = 0; i < block_size; i++){
        for(k = 0; k < block_size; k++){
            for(j = 0; j < block_size; j++){ // efficient matrix multiplication
                subC[i][j] += subA[i][k + block_nbr * block_size] * subB[k][j];
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

    int rank, size, col_rank, mat_size, block_size, i, j,k, p_x,p_y, dest, offset;
    double **A_rows;
    double **B_blocks;
    double **C_blocks;

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
    MPI_Status status4;
    

    MPI_Request request1;
    MPI_Request request2;
    MPI_Request request3;
    MPI_Request request4;
    

    MPI_Init(&argc, &argv);       
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    
    p_x = sqrt(size);
    p_y = sqrt(size);
    /* Create a virtual 2D-grid topology */
    /* MPI_Dims_create(size, 2, dims);  
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &block_comm); 
    MPI_Comm_rank(block_comm, &rank); 
    */

    MPI_Comm col_comm;
    //MPI_Comm_split(MPI_COMM_WORLD, rank / 4, 0, &row_comm);
    int color = rank %  p_x;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &col_comm);
    MPI_Comm_rank(col_comm, &col_rank);

    printf("col rank %i\n", col_rank);
    

    block_size = mat_size / p_x;

   
    
	/* Generate random Matrices */
    if (rank == 0){
        printf("mat %i block %i rank %i\n",mat_size, block_size, rank);
        for (i = 0; i < mat_size; i++){
	        for (j = 0; j < mat_size; j++){
                A[i][j] = i+j;   /* to be Replace by random */
                if(i == j)
                    B[i][j] = 1;
                else
                    B[i][j] = 0;
	        }
        }

        for (i = 0; i < mat_size; i++){
	        printf("\n");
		    for (j = 0; j < mat_size; j++){
		        printf("%.1f\t", A[i][j]);
		    }
	    }
	    printf("\n");
        for (i = 0; i < mat_size; i++){
	        printf("\n");
		    for (j = 0; j < mat_size; j++){
		        printf("%.1f\t", B[i][j]);
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

                //for matrix B

                             for (k = 0; k < block_size;k++){
       		    	    MPI_Isend(&B[i*block_size + k][j*block_size], block_size, MPI_DOUBLE, i*p_x + j, 200+k, MPI_COMM_WORLD, &request4);
		            }
                             
		        }


                

      	    }
       		offset = offset + block_size;
       	}
		
    }

	if (rank > 0){
	    /* Collect A rows */
        MPI_Recv(&mat_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);  
	    MPI_Recv(&block_size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status2);  
	    A_rows = (double **)malloc(block_size * sizeof(double *));
        for(i = 0; i < block_size; i++){
            A_rows[i] = (double *)malloc(mat_size * sizeof(double));
        }

		for (k = 0; k < block_size;k++){
			MPI_Recv(A_rows[k], mat_size, MPI_DOUBLE, 0, 100+k, MPI_COMM_WORLD, &status3);  
		}

          /* Collect B blocks */

           B_blocks= (double **)malloc(block_size * sizeof(double *));
        for(i = 0; i < block_size; i++){
            B_blocks[i] = (double *)malloc(block_size * sizeof(double));
        }

                for (k = 0; k < block_size;k++){
			MPI_Recv(B_blocks[k], block_size, MPI_DOUBLE, 0, 200+k, MPI_COMM_WORLD, &status4);  
		}
                
         

		/*Printing part */
		sleep(rank*2);
		printf("mat %i block %i rank %i\n",mat_size, block_size, rank);
		for(i = 0; i < block_size; i++){
		    printf("\n");
			for(j = 0; j < mat_size; j++){
				printf("%.1f\t", A_rows[i][j]);
			}   
		} 
            printf("\n");
	    printf("\n");
		for(i = 0; i < block_size; i++){
		    printf("\n");
			for(j = 0; j < block_size; j++){
				printf("%.1f\t", B_blocks[i][j]);
			}
		} 
	
		printf("\n");
		
		p_y = mat_size / block_size;
	    C_blocks = (double **)malloc(block_size * sizeof(double *));
        for(i = 0; i < block_size; i++){
            C_blocks[i] = (double *)malloc(block_size * sizeof(double));
        }
		/* Apply fox algorithm on blocks */
		for(i = 0; i < 1; i++){
			Block_matmul(A_rows, B_blocks, C_blocks, block_size, col_rank);
			//MPI_Send(&B_blocks, 1, blocktype, (col_rank + 1) % p_y, 1234, col_comm);
			//MPI_Recv(&B_blocks, 1, blocktype, (col_rank - 1) % p_y, 1234, col_comm, &status);
						
		}
		printf("\n");
	    printf("\n");
		for(i = 0; i < block_size; i++){
		    printf("\n");
			for(j = 0; j < block_size; j++){
				printf("%.1f\t", C_blocks[i][j]);
			}
		} 
              
        printf("\n");
              //1. broadcast Aim
              //2. call block_matmul
              //3. shift the B blocks

        /*Free all allocations */
        for(i = 0; i < block_size; i++){
            free(A_rows[i]);
        }
        free(A_rows);

               //fox();
	
	
		/* Send result */
		//MPI_send()
    }
	
	if (rank == 0){   //put that part in the first if ?
			//MPI_recv()
			//sum up results
    }
    MPI_Comm_free(&col_comm);
    MPI_Finalize();
  
    return 0;
}
