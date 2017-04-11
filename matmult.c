#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

/* Function to perform the step 2 of Fox's algorithm */
void Block_matmul(double *subA, double *subB, double *subC, int block_size, int mat_size, int block_nbr){
    int i,j,k;
    for(i = 0; i < block_size; i++){
        for(k = 0; k < block_size; k++){
            for(j = 0; j < block_size; j++){ // efficient matrix multiplication
                subC[i*block_size + j] += subA[i*mat_size + k + block_nbr * block_size] * subB[k * block_size + j];
            }
        }
    }     
}

/* Function to Print matrix */
/*void PrMat(int mat_size, double* matrix){   
    int i, j;
    printf("\n");
    for (i = 0; i < mat_size; i++){
        for (j = 0; j < mat_size; j++){
            printf("%.1f\t", matrix[i*mat_size + j]);
        }
        printf("\n");
    }
    printf("\n");
}*/





int main(int argc, char *argv[]) {

    int rank, nprocs,row_rank, col_rank, mat_size, block_size, i, j,k, k1, p_x,p_y, dest, offset;
    double *A_rows;
    double *B_blocks;
    double *C_blocks;

    double begin, end; // for time measurements

    mat_size = atoi(argv[1]);
    double A[mat_size * mat_size];
    double B[mat_size * mat_size];
    double C[mat_size * mat_size];

    MPI_Status status;
    MPI_Status status2;
    MPI_Status status3;
    MPI_Status status4;
    MPI_Status status5;
    

    MPI_Request request1;
    MPI_Request request2;
    MPI_Request request3;
    MPI_Request request4;
    MPI_Request request5;

    MPI_Init(&argc, &argv);       
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    
    p_x = sqrt(nprocs);
    p_y = sqrt(nprocs);
    block_size = mat_size / p_x;
   
    MPI_Comm row_comm;
    MPI_Comm col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank / p_x, rank, &row_comm);
    int color = rank %  p_x;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &col_comm);
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_rank(row_comm, &row_rank);
    
    MPI_Datatype rawstype, blocktype, blockselect, blockselect2;

    

    MPI_Type_contiguous(block_size * mat_size,  MPI_DOUBLE, &rawstype);
    MPI_Type_commit(&rawstype);
    MPI_Type_contiguous(block_size * block_size, MPI_DOUBLE, &blocktype);
    MPI_Type_commit(&blocktype);
    MPI_Type_vector(block_size, block_size, mat_size, MPI_DOUBLE, &blockselect2);
    MPI_Type_commit(&blockselect2);
    MPI_Type_create_resized(blockselect2, 0, sizeof(double), &blockselect);
    MPI_Type_commit(&blockselect);

/**********************************************************************
 * Master processor work 
 * *******************************************************************/
   
    
	/* Generate random Matrices */
    if (rank == 0){
        printf("mat %i block %i rank %i\n",mat_size, block_size, rank);
        for (i = 0; i < mat_size; i++){
	        for (j = 0; j < mat_size; j++){
                A[i * mat_size + j] = i+j;   /* to be Replace by random */
                B[i * mat_size + j] = i*j;
	        }
        }

        for (i = 0; i < mat_size; i++){
	        printf("\n");
		    for (j = 0; j < mat_size; j++){
		        printf("%.1f\t", A[i * mat_size + j]);
		    }
	    }

	    printf("\n");
        for (i = 0; i < mat_size; i++){
	        printf("\n");
		    for (j = 0; j < mat_size; j++){
		        printf("%.1f\t", B[i * mat_size + j]);
		    }
	    }
	    printf("\n");


	    /* Distributes Blocks from A and B */
	    offset = 0;
	
        for (i = 0; i < p_x; i++){
		    for (j = 0; j < p_y; j++){
                
                if ((i+j) != 0){
                    //for matrix A
		            MPI_Isend(&mat_size, 1, MPI_INT, i*p_x + j, 1, MPI_COMM_WORLD, &request1);
		    	    MPI_Isend(&block_size, 1, MPI_INT, i*p_x + j, 2, MPI_COMM_WORLD, &request2);

       		    	MPI_Isend(&A[offset], 1, rawstype, i*p_x + j, 666, MPI_COMM_WORLD, &request3);

                    //for matrix B

                    MPI_Isend(&B[i*block_size*mat_size + j*block_size], 1, blockselect, i*p_x + j, 999, MPI_COMM_WORLD, &request4);
		                                         
		        }

      	    }
       		offset = offset + block_size * mat_size;
       	}

		/* Blocks for P0 don't need to be send */
	    A_rows = (double *)malloc(block_size * mat_size * sizeof(double));
        for (i = 0; i < mat_size * block_size; i++)
            A_rows[i] = A[i];

        B_blocks= (double *)malloc(block_size * block_size * sizeof(double));
        for(i = 0; i < block_size; i++){
            for(j = 0; j < block_size; j++){
                B_blocks[i * block_size + j] = B[i * mat_size + j];
            }
        }
    }

/**********************************************************************
 * Collecting Blocks of initial matrix
 * *******************************************************************/

	if (rank > 0){ 

	    /* Collect A rows */
        MPI_Recv(&mat_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);  
	    MPI_Recv(&block_size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status2);  
	    
        A_rows = (double *)malloc(block_size * mat_size * sizeof(double));
		MPI_Recv(A_rows, 1, rawstype, 0, 666, MPI_COMM_WORLD, &status3);  

          /* Collect B blocks */

        B_blocks= (double *)malloc(block_size * block_size * sizeof(double));
		MPI_Recv(B_blocks, 1, blocktype, 0, 999, MPI_COMM_WORLD, &status4);    
    }

    
   /* Printing part */
    sleep(rank);
    printf("mat %i block %i rank %i\n",mat_size, block_size, rank);
    for(i = 0; i < block_size; i++){
        printf("\n");
        for(j = 0; j < mat_size; j++){
            printf("%.1f\t", A_rows[i * mat_size + j]);
        }   
    } 
    printf("\n");
    printf("\n");
    for(i = 0; i < block_size; i++){
        printf("\n");
        for(j = 0; j < block_size; j++){
            printf("%.1f\t", B_blocks[i * block_size + j]);
        }
    } 
	
    printf("\n"); 
    
/**********************************************************************
 * FOX Algorithm
 * *******************************************************************/
    begin = MPI_Wtime();

    
    C_blocks = (double *)malloc(block_size * block_size * sizeof(double *));

    /* Apply fox algorithm on blocks */
    for(i = 0; i < p_x; i++){
        MPI_Isend(B_blocks, 1, blocktype, (col_rank - 1 + p_y) % p_y, 1000 + i, col_comm, &request4);
        Block_matmul(A_rows, B_blocks, C_blocks, block_size, mat_size, (col_rank + i) % p_y);
        MPI_Recv(B_blocks, 1, blocktype, (col_rank + 1 + p_y) % p_y, 1000 + i, col_comm, &status);
    }

    end = MPI_Wtime();

    
    sleep(5 + rank);
    printf("\n");
    printf("C_blocks in proc %i\n",rank);
    for(i = 0; i < block_size; i++){
        printf("\n");
        for(j = 0; j < block_size; j++){
            printf("%.1f\t", C_blocks[i*block_size + j]);
        }
    } 
          
    printf("\n"); 

    /* Send C_blocks in processors(i,j) back to root processor */
    
    if(rank>0){

        for(i=0; i<block_size;i++){
           for (j=0; j<block_size; j++){
              MPI_Isend(&C_blocks[i*block_size*mat_size + j*block_size],1,blockselect,0,2000+ rank,MPI_COMM_WORLD,&request5);
           }
        }

    }
    
    /* Reception by root processor of all C_blocks and storing it C */	
   if (rank == 0){ 
       for(i=0; i<block_size;i++){
           for (j=0; j<block_size; j++){
               // C[i*mat_size + j] = C_blocks[i*block_size + j];   
           }
        }


       for (k=1; k< nprocs ; k++){
           MPI_Recv(&C_blocks,1,blocktype,k,2000+k,MPI_COMM_WORLD,&status5);
           for (i=0; i<block_size;i++){
               for (j=0; j<block_size; j++){
                 //  C[((k / p_x)*block_size + i)*mat_size + ((k%p_x)*block_size + j)] = C_blocks[i*block_size + j];
                }
           }
        }  

   printf("C : \n");
  //  PrMat(mat_size, C);
   printf("\n Fox's algorithm time: %g s\n\n", end - begin);
    
   }

    /*Free all allocations */
    /*for(i = 0; i < block_size; i++){
        free(C_blocks[i]);
    }*/
    free(A_rows);
    free(B_blocks);
    free(C_blocks);
    
    MPI_Comm_free(&row_comm);   
    MPI_Comm_free(&col_comm);
    MPI_Finalize();
  
    return 0;
}
