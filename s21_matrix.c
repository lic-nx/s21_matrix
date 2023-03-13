#include "s21_matrix.h"

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
  }
  A->columns = 0;
  A->rows = 0;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int returned = OK, cont = 0;
  s21_remove_matrix(result);
  if (rows > 0 && columns > 0) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = malloc(rows * sizeof(double));  // malloc calloc
    returned = Err;
    if (result->matrix != NULL) {
      for (int i = 0; i < rows && cont == 0; i++) {  // посмотреть кто зануляет
        result->matrix[i] = malloc(columns * sizeof(double));
        if (result->matrix[i] == NULL) {
          cont = 1;
        }
      }
      returned = OK;
    }
  }
  return returned;
}

int is_not_NULL(matrix_t *matrix) {
  int returned = 1;
  if (matrix != NULL && matrix->matrix != NULL && matrix->rows > 0 &&
      matrix->columns > 0) {
    returned = 1;
  } else {
    returned = 0;
  }
  return returned;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int returned = SUCCSES;
  if (A->rows == B->rows && A->columns == B->columns && is_not_NULL(A) &&
      is_not_NULL(B)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-6 && returned == 0) {
          returned = 0;
        } else {
          returned = 1;
        }
      }
    }
  } else {
    returned = FAILURE;
  }
  return returned;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int returned = OK;
  if (is_not_NULL(A) && is_not_NULL(B)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      returned = 2;
    }
  } else {returned = 1;}
  return returned;
}
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int returned = OK;
  if (is_not_NULL(A) && is_not_NULL(B)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
          //printf("%f ", result->matrix[i][j]);
        }
      }
    } else {
      returned = 2;
    }
  } else {returned = 1;}
  return returned;
}


int s21_mult_number(matrix_t *A, double number, matrix_t *result){
  int returned = OK;
   if (is_not_NULL(A)) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
          //printf("%f ", result->matrix[i][j]);
        }
      }
  } else {returned = 1;}
  return returned;
}
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result){
  int returned = OK;
  return returned;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result){
 int returned = OK;
 s21_create_matrix(A->rows, B->columns, result);
 if(A->rows == B->columns){
   int res = 0; 
   for(int j = 0; j< A->rows; j++){
    for(int i = 0; i< A->columns; i++){
      res += A->matrix[j][i] * B->matrix[i][j];
    }
   }
 }
 else{
   returned = ErrOper;
 }
 return returned;
}