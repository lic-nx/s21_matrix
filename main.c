#include "s21_matrix.h"

void main() {
    matrix_t A, B, C;
    s21_create_matrix(2,2, &A);
    s21_create_matrix(2,2, &B);
    s21_create_matrix(2,2, &C);
    for(int i = 0; i < A.rows; i++){
        for(int j = 0; j < A.columns;j++){
            A.matrix[i][j] = 1+i+j;
            B.matrix[i][j] = 1+i+j;
            printf("%f ", A.matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    s21_mult_matrix(&A, &B, &C);
    for(int i = 0; i < C.rows; i++){
        for(int j = 0; j < C.columns;j++){
            printf("%f ", C.matrix[i][j]);
        }
        printf("\n");
    }
}
