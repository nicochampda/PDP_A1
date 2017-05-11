#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#define randmin 0
#define randmax 50



/* Function to perform the step 2 of Fox's algorithm */
void Block_matmul(double *subA, double *subB, double *subC, int block_size){
    int i, j, k;
    for(i = 0; i < block_size; i++){
        for(j = 0; j < block_size; j++){
            for(k = 0; k < block_size; k++){ 
                subC[i*block_size + j] += subA[i*block_size + k] * subB[k*block_size + j];
            }
        }
    }     
}

/* Function to Print matrix */
void PrMat(int mat_size, double matrix[mat_size * mat_size]){   
    int i, j;
    printf("\n");
    for (i = 0; i < mat_size; i++){
        for (j = 0; j < mat_size; j++){
            printf("%.1f\t", matrix[i*mat_size + j]);
        }
        printf("\n");
    }
    printf("\n");
}

/* Swap pointers */
void swap(double **arr1, double **arr2){
    double *temp;
    temp = *arr1;
    *arr1 = *arr2;
    *arr2 = temp;
}



int main(int argc, char *argv[]) {

    int rank, nprocs,row_rank, col_rank, mat_size = 0, block_size = 0, i, j, p_x, dest, Print = 0;
    double *cur_A_blocks;
    double *init_A_blocks;

    double *cur_B_blocks;
    double *next_B_blocks;
    
    double *C_blocks;

    double begin, mid, end;

    if (argc < 2 || argc > 3){
        printf("usage : ./mpirun -np nprocs ./matmult mat_size [Print 0/1]\n");
        exit(0);
    }

    mat_size = atoi(argv[1]);
    if (argc == 3) Print = atoi(argv[2]);

    double *A = NULL;
    double *B = NULL;
    double *C = NULL;

    MPI_Request req_A_send;
    MPI_Request req_A_recv;
    MPI_Request req_A_bcast;
    MPI_Request req_B_send;
    MPI_Request req_B_recv;
    MPI_Request req_B_init;
    MPI_Request req_C_send;


    MPI_Init(&argc, &argv);       
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    p_x = sqrt(nprocs);

    if (nprocs != p_x * p_x){
        printf("nprocs must be a perfect square number\n");
        MPI_Finalize();
        exit(0);
    }

    block_size = mat_size / p_x;

    if (mat_size%p_x != 0){
        printf("maxtrix size must be divisible by sqrt(nprocs)\n");
        MPI_Finalize();
        exit(0);
    }


    MPI_Request req_arr[nprocs];
    /* Define Row and Column communicators */
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank / p_x, rank, &row_comm);
    MPI_Comm_rank(row_comm, &row_rank);

    MPI_Comm col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank % p_x, rank, &col_comm);
    MPI_Comm_rank(col_comm, &col_rank);
    

    /* Define new type to send datas */
    MPI_Datatype blocktype, blockselect;

    MPI_Type_contiguous(block_size * block_size, MPI_DOUBLE, &blocktype);
    MPI_Type_commit(&blocktype);

    MPI_Type_vector(block_size, block_size, mat_size, MPI_DOUBLE, &blockselect);
    MPI_Type_commit(&blockselect);

/**********************************************************************
 * Master processor work 
 * *******************************************************************/
   
    srand(time(NULL)); 
	/* Generate random Matrices */
    if (rank == 0){
        A = (double *)malloc(sizeof(double) * mat_size * mat_size);
        B = (double *)malloc(sizeof(double) * mat_size * mat_size);
        C = (double *)malloc(sizeof(double) * mat_size * mat_size);

        for (i = 0; i < mat_size; i++){
	        for (j = 0; j < mat_size; j++){
                A[i * mat_size + j] = (rand()/(double)RAND_MAX )*(randmax-randmin) + randmin;  
                B[i * mat_size + j] = (rand()/(double)RAND_MAX )*(randmax-randmin) + randmin; 
	        }
        }

        begin = MPI_Wtime();
	    
        /* Distributes Blocks from A and B */
        for (i = 0; i < p_x; i++){
            int i_b_m = i * block_size * mat_size;
            dest = i * p_x;
		    for (j = 0; j < p_x; j++){
                MPI_Isend(&A[i_b_m + j*block_size], 1, blockselect, dest + j, 3, MPI_COMM_WORLD, &req_A_send);
                MPI_Isend(&B[i_b_m + j*block_size], 1, blockselect, dest + j, 4, MPI_COMM_WORLD, &req_B_init);

      	    }
       	}
    }

/**********************************************************************
 * Collecting Blocks of initial matrix
 * *******************************************************************/
    /* Collect A Blocks */
    init_A_blocks = (double *)malloc(block_size * block_size * sizeof(double));
    MPI_Irecv(init_A_blocks, 1, blocktype, 0, 3, MPI_COMM_WORLD, &req_A_recv);
    
    
    /* Collect B blocks */
    cur_B_blocks = (double *)malloc(block_size * block_size * sizeof(double));
    MPI_Irecv(cur_B_blocks, 1, blocktype, 0, 4, MPI_COMM_WORLD, &req_B_recv);

/**********************************************************************
 * FOX Algorithm
 * *******************************************************************/

    cur_A_blocks = (double *)malloc(block_size * block_size * sizeof(double));
    next_B_blocks = (double *)malloc(block_size * block_size * sizeof(double));
    /* Initialisation of the C blocks */
    C_blocks = (double *)malloc(block_size * block_size * sizeof(double));

    for (i = 0; i < block_size * block_size; i++)
        C_blocks[i] = 0;

    MPI_Wait(&req_A_recv, MPI_STATUS_IGNORE);
    MPI_Wait(&req_B_recv, MPI_STATUS_IGNORE);

    double mid1 = MPI_Wtime();
    /* Apply fox algorithm on blocks */
    for(i = 0; i < p_x; i++){
        /* Broadcast the usefull block of A */
        if (row_rank == (col_rank + i)%p_x){
            memcpy(cur_A_blocks, init_A_blocks, block_size * block_size*sizeof(double));
        }

        MPI_Bcast(cur_A_blocks, 1, blocktype, (col_rank + i)%p_x, row_comm);

        /* Shift B blocks upward */
        MPI_Irecv(next_B_blocks, 1, blocktype, (col_rank + 1) % p_x, 100 + i, col_comm, &req_B_recv);
        MPI_Isend(cur_B_blocks, 1, blocktype, (col_rank - 1 + p_x) % p_x, 100 + i, col_comm, &req_B_send);
        
        Block_matmul(cur_A_blocks, cur_B_blocks, C_blocks, block_size);
        
        MPI_Wait(&req_B_recv, MPI_STATUS_IGNORE);
        MPI_Wait(&req_B_send, MPI_STATUS_IGNORE);
        
        swap(&cur_B_blocks, &next_B_blocks);
    }


    /* Send C_blocks in processors(i,j) back to root processor */
    
    if(rank == 0) 
        printf("\n Fox's algorithm time: %f s\n\n", MPI_Wtime() - mid1);
    MPI_Isend(C_blocks, 1, blocktype, 0, 50, MPI_COMM_WORLD, &req_C_send);
    
    /* Waiting of all C blocks */	
    if (rank == 0){ 
        for (i = 0; i < p_x; i++){
            int i_b_m = i * block_size * mat_size;
            dest = i * p_x;
		    for (j = 0; j < p_x; j++){
                MPI_Irecv(&C[i_b_m + j*block_size], 1, blockselect, dest + j, 50, MPI_COMM_WORLD, &req_arr[dest + j]);
            }
        }
                
        MPI_Waitall(nprocs, req_arr, MPI_STATUS_IGNORE);
        end = MPI_Wtime();
            
        if (Print == 1){
            printf("A : \n");
            PrMat(mat_size, A);
            printf("B : \n");
            PrMat(mat_size, B);
            printf("A * B = C : \n");
            PrMat(mat_size, C);
        }

        free(A);
        free(B);
        free(C);
        printf("\n Total run time: %g s\n\n", end - begin);
    } 

    MPI_Wait(&req_C_send, MPI_STATUSES_IGNORE);

    free(cur_A_blocks);
    free(init_A_blocks);
    free(cur_B_blocks);
    free(next_B_blocks);
    free(C_blocks);

    MPI_Type_free(&blockselect);
    MPI_Type_free(&blocktype);

    MPI_Comm_free(&row_comm);   
    MPI_Comm_free(&col_comm);
    MPI_Finalize();
  
    return 0;
}
