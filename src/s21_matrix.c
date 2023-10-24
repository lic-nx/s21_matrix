#include "s21_matrix.h"

int egual(matrix_t A, matrix_t B) {
  int returned = 0;
  if (is_Emty(&A) && is_Emty(&B))
    if (A.columns == B.columns && A.rows == B.rows) {
      returned = 1;
    }
  return returned;
}

int is_Emty(matrix_t *matrix) {
  int returned = 0;
  if (matrix != NULL && (matrix->matrix != 0 && matrix->matrix != NULL) &&
      matrix->rows >= 1 && matrix->columns >= 1) {
    returned = 1;
  } else {
    matrix->columns = 0;
    matrix->matrix = NULL;
    matrix->rows = 0;
  }
  return returned;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
  }
  free(A->matrix);
  A->columns = 0;
  A->rows = 0;
  A->matrix = NULL;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int returned = OK;
  if (rows > 0 && columns > 0 && rows < 9999 && columns < 9999) {
    if (returned == OK) {
      result->columns = columns;
      result->rows = rows;
      result->matrix = calloc(rows, sizeof(double));  // malloc calloc

      if (result->matrix != NULL) {
        for (int i = 0; i < rows && returned != ErrOper; i++) {
          result->matrix[i] = calloc(columns, sizeof(double));
          if (result->matrix[i] == NULL) {
            result->rows = i;
            returned = ErrOper;
            s21_remove_matrix(result);
          }
        }
      } else {
        returned = ErrOper;
      }
    }
  } else {
    returned = 1;
  }
  return returned;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int returned = SUCCESS;

  if (egual(*A, *B)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabsl(A->matrix[i][j] - B->matrix[i][j]) < MINFOREQ &&
            returned == SUCCESS) {
        } else {
          returned = FAILURE;
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
  if (egual(*A, *B)) {
    returned = s21_create_matrix(A->rows, A->columns, result);

    for (int i = 0; i < A->rows && returned == OK; i++) {
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  } else {
    returned = 2;
  }

  return returned;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int returned = OK;
  if (egual(*A, *B)) {
    returned = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && returned == OK; i++) {
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  } else {
    returned = 2;
  }
  return returned;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int returned = OK;
  if (is_Emty(A)) {
    // if (result != NULL &&
    //     (A->rows != result->rows || A->columns != result->columns)) {
    //   s21_remove_matrix(result);
    // }
    returned = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] += A->matrix[i][j] * number;
      }
    }
  } else {
    returned = 1;
  }

  return returned;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int returned = OK;
  if (A->columns == B->rows && is_Emty(A) && is_Emty(B)) {
    returned = s21_create_matrix(A->rows, B->columns, result);
    if (returned == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          result->matrix[i][j] = 0;
          for (int k = 0; k < B->rows; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else {
      returned = ErrOper;
    }
  } else {
    returned = 1;
  }
  return returned;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int returned = OK;
  returned = s21_create_matrix(A->columns, A->rows, result);
  if (returned == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return returned;
}

// минор
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int returned = OK;
  if (A->columns != A->rows) {
    returned = 1;
  }

  if (returned == OK) {
    if (A->columns > 1 && A->rows > 1) {
      returned = s21_create_matrix(A->rows, A->columns, result);
      double mull = 0;
      int Pow = 1;
      matrix_t minor = {.matrix = NULL, .columns = 0, .rows = 0};
      s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          Pow = pow((-1), i + j);
          s21_minor(A, i, j, &minor);
          s21_determinant(&minor, &mull);
          result->matrix[i][j] = mull;
          result->matrix[i][j] *= Pow;
        }
      }

      s21_remove_matrix(&minor);
    } else {
      returned = 1;
    }
  }
  return returned;
}

int s21_minor(matrix_t *A, int i, int j, matrix_t *result) {
  int err = OK;
  // err = s21_create_matrix(A->columns - 1, A->rows - 1, result);
  if (err == OK) {
    for (int rows = 0, m = 0; rows < A->rows; rows++, m++) {
      for (int columns = 0, n = 0; columns < A->rows; columns++, n++) {
        if (columns == j) columns++;
        if (rows == i) rows++;
        if (rows != i && columns != j && columns < A->columns &&
            rows < A->rows) {
          result->matrix[m][n] = A->matrix[rows][columns];
        }
      }
    }
  }
  return err;
}

double revers_determinant(matrix_t *A) {
  double flag = 0.0;
  if (A->rows == 1) {
    flag = A->matrix[0][0];
  } else if (A->rows == 2 && A->columns == 2) {
    flag =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else if (A->rows > 2) {
    matrix_t tmp = {0};
    s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
    for (int i = 0; i < A->columns; i++) {
      s21_minor(A, 0, i, &tmp);
      if (i % 2) {
        flag -= A->matrix[0][i] * revers_determinant(&tmp);
      } else {
        flag += A->matrix[0][i] * revers_determinant(&tmp);
      }
    }
    s21_remove_matrix(&tmp);
  }
  return flag;
}

int s21_determinant(matrix_t *A, double *result) {
  int returned = OK;
  *result = 0;
  if (A->columns != A->rows || !result) {
    returned = 1;
  } else if (A->columns == 1 && A->rows == 1) {
    *result = A->matrix[0][0];
  } else if (is_Emty(A) && A->columns == A->rows) {
    *result = revers_determinant(A);

  } else {
    returned = 1;
  }
  return returned;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int returned = OK;
  double determ;
  returned = s21_determinant(A, &determ);
  if (returned == OK && determ != 0.0) {
    // returned = s21_create_matrix(A->columns, A->rows, result);
    if (returned == OK && A->rows > 1) {
      matrix_t Atrans = {.matrix = NULL, .columns = 0, .rows = 0},
               minor = {.matrix = NULL, .columns = 0, .rows = 0};
      returned = s21_calc_complements(A, &minor);
      if (returned == OK) {
        // s21_create_matrix(minor.columns, minor.rows, &Atrans);
        returned = s21_transpose(&minor, &Atrans);
        if (returned == OK) s21_mult_number(&Atrans, 1 / determ, result);
        s21_remove_matrix(&Atrans);
      }

      s21_remove_matrix(&minor);
    } else {
      returned = 2;
    }
  } else {
    returned = 2;
  }
  return returned;
}
