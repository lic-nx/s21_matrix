#include "s21_matrix.h"

void main() {
    matrix_t A, B, C;
    s21_create_matrix(6,6, &A);
    s21_create_matrix(6,6, &B);
    s21_create_matrix(6,6, &C);
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6;j++){
            A.matrix[i][j] = 1+i+j;
            B.matrix[i][j] = 1+i+j;
            printf("%f ", A.matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    s21_sum_matrix(&A, &B, &C);
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6;j++){
            printf("%f ", C.matrix[i][j]);
        }
        printf("\n");
    }
}
