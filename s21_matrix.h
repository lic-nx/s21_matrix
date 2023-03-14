#include <math.h>
#include <stdlib.h>
#include "stdio.h"

#define SUCCSES 1
#define FAILURE 0

#define OK 0
#define Err 1 // Ошибка, некорректная матрица
#define ErrOper 2 // Ошибка вычисления 
// (несовпадающие размеры матриц; матрица, для которой нельзя провести вычисления и т.д.)

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int is_not_NULL(matrix_t *matrix);
int s21_eq_matrix(matrix_t *A, matrix_t *B);

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);


int s21_transpose(matrix_t *A, matrix_t *result);

//int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
double revers_determinant(matrix_t *A, double result);