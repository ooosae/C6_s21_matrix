#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int rc = OK;
  if (rows >= 1 && columns >= 1) {
    if (result) {
      result->rows = rows;
      result->columns = columns;
      result->matrix = (double **)calloc(rows, sizeof(double *));
      if (result->matrix) {
        for (int i = 0; i < rows; ++i) {
          result->matrix[i] = (double *)calloc(columns, sizeof(double));
          if (!result->matrix[i]) {
            s21_remove_matrix(result);
            rc = EXIT_CALC_ERROR;
            break;
          }
          for (int j = 0; j < columns; ++j) result->matrix[i][j] = 0.0;
        }
      } else {
        s21_remove_matrix(result);
        rc = EXIT_CALC_ERROR;
      }
    } else {
      rc = EXIT_INCORRECT_MATRIX_ERROR;
    }
  } else {
    rc = EXIT_CALC_ERROR;
  }

  return rc;
}

void s21_remove_matrix(matrix_t *A) {
  if (!A || !A->matrix) return;
  for (int i = 0; i < A->rows; ++i) free(A->matrix[i]);
  free(A->matrix);
  A->rows = 0;
  A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int rc = SUCCESS;
  if (A && B) {
    if (A->rows == B->rows && A->columns == B->columns && A->rows > 0 &&
        A->columns > 0) {
      if (A->matrix && B->matrix) {
        for (int i = 0; i < A->rows; ++i)
          for (int j = 0; j < A->columns; ++j)
            if (fabs(A->matrix[i][j] - B->matrix[i][j]) > eps) {
              rc = FAILURE;
              break;
            }
      } else {
        rc = FAILURE;
      }
    } else {
      rc = FAILURE;
    }
  } else {
    rc = FAILURE;
  }
  return rc;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int rc = OK;
  if (A && B && result && A->rows > 0 && A->columns > 0 && B->rows > 0 &&
      B->columns > 0 && A->matrix && B->matrix) {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < A->columns; ++j)
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    } else {
      rc = EXIT_CALC_ERROR;
    }
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }

  return rc;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int rc = OK;
  if (A && B && result && A->rows > 0 && A->columns > 0 && B->rows > 0 &&
      B->columns > 0 && A->matrix && B->matrix) {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < A->columns; ++j)
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    } else {
      rc = EXIT_CALC_ERROR;
    }
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }

  return rc;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int rc = OK;
  if (A && A->rows > 0 && A->columns > 0 && A->matrix) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; ++i)
      for (int j = 0; j < A->columns; ++j)
        result->matrix[i][j] = A->matrix[i][j] * number;
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }
  return rc;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int rc = OK;
  if (A && B && result && A->rows > 0 && A->columns > 0 && B->rows > 0 &&
      B->columns > 0 && A->matrix && B->matrix) {
    if (A->columns == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < B->columns; ++j)
          for (int k = 0; k < A->columns; ++k)
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
    } else {
      rc = EXIT_CALC_ERROR;
    }
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }
  return rc;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int rc = OK;
  if (A && result && A->rows > 0 && A->columns > 0 && A->matrix) {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; ++i)
      for (int j = 0; j < A->columns; ++j)
        result->matrix[j][i] = A->matrix[i][j];
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }
  return rc;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int rc = OK;
  if (A && result && A->rows > 0 && A->columns > 0 && A->matrix) {
    if (A->rows == A->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; ++i)
        for (int j = 0; j < A->columns; ++j) {
          double det = 0.0;
          matrix_t minor = {0};
          get_minor_ij(i + 1, j + 1, A, &minor);
          s21_determinant(&minor, &det);
          s21_remove_matrix(&minor);
          result->matrix[i][j] = pow((-1), i + j) * det;
        }
    } else {
      rc = EXIT_CALC_ERROR;
    }
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }
  return rc;
}

int s21_determinant(matrix_t *A, double *result) {
  int rc = OK;
  if (A && result && A->rows > 0 && A->columns > 0 && A->matrix) {
    if (A->rows == A->columns) {
      *result = 0;
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows == 2) {
        *result = A->matrix[0][0] * A->matrix[1][1] -
                  A->matrix[0][1] * A->matrix[1][0];
      } else {
        matrix_t minor = {0};
        for (int i = 0; i < A->rows; i++) {
          double res_minor = 0.0;
          get_minor_ij(1, i + 1, A, &minor);
          s21_determinant(&minor, &res_minor);
          s21_remove_matrix(&minor);
          *result += pow((-1), i) * A->matrix[0][i] * res_minor;
        }
      }
    } else {
      rc = EXIT_CALC_ERROR;
    }
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }
  return rc;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  double det = 0;
  int rc = s21_determinant(A, &det);
  if (A && result && A->rows > 0 && A->columns > 0 && A->matrix) {
    if (A->rows == A->columns) {
      matrix_t tmp_calc = {0};
      if (!s21_calc_complements(A, &tmp_calc)) {
        matrix_t tmp_tr = {0};
        if (!s21_transpose(&tmp_calc, &tmp_tr)) {
          s21_mult_number(&tmp_tr, 1 / det, result);
        }
        s21_remove_matrix(&tmp_tr);
      }
      s21_remove_matrix(&tmp_calc);
    } else {
      rc = EXIT_CALC_ERROR;
    }
  } else {
    rc = EXIT_INCORRECT_MATRIX_ERROR;
  }
  return rc;
}

int get_minor_ij(int row, int column, matrix_t *A, matrix_t *result) {
  int rc = EXIT_INCORRECT_MATRIX_ERROR;
  if (A && A->matrix && result && A->rows > 0 && A->columns > 0) {
    if (row > 0 && column > 0 && A->rows == A->columns) {
      rc = s21_create_matrix(A->rows - 1, A->columns - 1, result);
      if (rc == OK)
        for (int i = 0; i < A->rows; ++i) {
          int tmp_row = i;
          if (i > row - 1) --tmp_row;
          for (int j = 0; j < A->columns; ++j) {
            int tmp_col = j;
            if (j > column - 1) --tmp_col;
            if (i != row - 1 && j != column - 1)
              result->matrix[tmp_row][tmp_col] = A->matrix[i][j];
          }
        }
    } else {
      rc = EXIT_CALC_ERROR;
    }
  }
  return rc;
}
