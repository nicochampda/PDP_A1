#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define randmin 0
#define randmax 50




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
void PrMat(int mat_size, double matrix[mat_size * mat_size]){   
    int i, j;
    printf("\n");
    for (i = 0; i < mat_size; i++){
        for (j = 0; j < mat_size; j++){
            printf("%.3f\t", matrix[i*mat_size + j]);
        }
        printf("\n");
    }
    printf("\n");
}





int main(int argc, char *argv[]) {

    int rank, nprocs,row_rank, col_rank, mat_size, block_size, i, j, p_x, dest, Print;
    double *A_rows;
    double *B_blocks;
    double *C_blocks;

    double begin = 0, end; // for time measurements

    if (argc < 2 || argc > 3){
        printf("usage : ./mpirun -np nprocs ./matmult mat_size [Print 0/1]\n");
        exit(0);
    }

    mat_size = atoi(argv[1]);
    if (argc == 3) Print = atoi(argv[2]);

    double A[mat_size * mat_size];
    double B[mat_size * mat_size];
    double C[mat_size * mat_size];

    MPI_Request request3;
    MPI_Request request4;
    MPI_Request request5;


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
    MPI_Datatype rowstype, blocktype, blockselect;

    MPI_Type_contiguous(block_size * mat_size,  MPI_DOUBLE, &rowstype);
    MPI_Type_commit(&rowstype);

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
        for (i = 0; i < mat_size; i++){
	        for (j = 0; j < mat_size; j++){
                A[i * mat_size + j] =( rand()/(double)RAND_MAX )*(randmax-randmin) + randmin;  
                B[i * mat_size + j] =( rand()/(double)RAND_MAX )*(randmax-randmin) + randmin; 
	        }
        }


        begin = MPI_Wtime();
	    
        
        /* Distributes Blocks from A and B */
        for (i = 0; i < p_x; i++){
            MPI_Isend(&A[i * block_size * mat_size], 1, rowstype, i, 3, col_comm, &request3);
		    for (j = 0; j < p_x; j++){
                dest = i * p_x + j;
                MPI_Isend(&B[i*block_size*mat_size + j*block_size], 1, blockselect, dest, 4, MPI_COMM_WORLD, &request4);
      	    }
       	}
    }

/**********************************************************************
 * Collecting Blocks of initial matrix
 * *******************************************************************/
    MPI_Bcast(&mat_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    A_rows = (double *)malloc(block_size * mat_size * sizeof(double));
    /* Collect A rows in the first column of processors */
    if (row_rank == 0){
        MPI_Recv(A_rows, 1, rowstype, 0, 3, col_comm, MPI_STATUS_IGNORE);  
    }

    /* Broadcast A_rows on rows of processors */
    MPI_Bcast(A_rows, 1, rowstype, 0, row_comm);

    /* Collect B blocks */
    B_blocks= (double *)malloc(block_size * block_size * sizeof(double));
    MPI_Recv(B_blocks, 1, blocktype, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    
    
/**********************************************************************
 * FOX Algorithm
 * *******************************************************************/

    
    C_blocks = (double *)malloc(block_size * block_size * sizeof(double *));

    /* Apply fox algorithm on blocks */
    for(i = 0; i < p_x; i++){
        MPI_Isend(B_blocks, 1, blocktype, (col_rank - 1 + p_x) % p_x, 100 + i, col_comm, &request4);
        Block_matmul(A_rows, B_blocks, C_blocks, block_size, mat_size, (col_rank + i) % p_x);
        MPI_Recv(B_blocks, 1, blocktype, (col_rank + 1) % p_x, 100 + i, col_comm, MPI_STATUS_IGNORE);
    }

    
    /* Send C_blocks in processors(i,j) back to root processor */
    
    MPI_Isend(C_blocks, 1, blocktype, 0, 50, MPI_COMM_WORLD, &request5);
    
    /* Reception by root processor of all C_blocks and storing it C */	
    if (rank == 0){ 
        for (i = 0; i < p_x; i++){
            for (j = 0; j < p_x; j++){
                dest = i * p_x + j;
                MPI_Irecv(&C[i * block_size * mat_size + j * block_size], 1, blockselect, dest, 50, MPI_COMM_WORLD, &req_arr[dest]);
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

        printf("\n Fox's algorithm time: %g s\n\n", end - begin);
    } 


    free(A_rows);
    free(B_blocks);
    free(C_blocks);
    
    MPI_Comm_free(&row_comm);   
    MPI_Comm_free(&col_comm);
    MPI_Finalize();
  
    return 0;
}
