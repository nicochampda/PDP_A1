#include <stdio.h>
#include <stdlib.h>

/* Function to Print matrix */
void PrMat(int mat_size, double matrix[mat_size][mat_size]){   
    int i, j;
    printf("\n");
    for (i = 0; i < mat_size; i++){
        for (j = 0; j < mat_size; j++){
            printf("%.1f\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    int mat_size = 8;
    double A[mat_size][mat_size];
    double B[mat_size][mat_size];
    double C[mat_size][mat_size];
    
    int i, j, k;
    /* Fill matrices */

    for (i = 0; i < mat_size; i++){
        for (j = 0; j < mat_size; j++){
            A[i][j] = i+j;
            B[i][j] = i*j;
            C[i][j] = 0;
        }
    }

    /* Compute product */
    for(i = 0; i < mat_size; i++){
        for(k = 0; k < mat_size; k++){
            for(j = 0; j < mat_size; j++){ 
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    /* print A, B, C */
    printf("A : \n");
    PrMat(mat_size, A);
    printf("B : \n");
    PrMat(mat_size, B);
    printf("C : \n");
    PrMat(mat_size, C);

return 0;
}
