/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xtrsm.cpp
 *
 * Code generation for function 'xtrsm'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xtrsm.h"

/* Function Definitions */
void xtrsm(int m, int n, const emxArray_real_T *A, int lda, emxArray_real_T *B,
           int ldb)
{
  int j;
  int jBcol;
  int k;
  int kAcol;
  double x;
  double y;
  int i;
  if ((n == 0) || ((B->size[0] == 0) || (B->size[1] == 0))) {
  } else {
    for (j = 1; j <= n; j++) {
      jBcol = ldb * (j - 1) - 1;
      for (k = m; k > 0; k--) {
        kAcol = lda * (k - 1) - 1;
        if (B->data[k + jBcol] != 0.0) {
          x = B->data[k + jBcol];
          y = A->data[k + kAcol];
          B->data[k + jBcol] = x / y;
          for (i = 1; i < k; i++) {
            B->data[i + jBcol] -= B->data[k + jBcol] * A->data[i + kAcol];
          }
        }
      }
    }
  }
}

/* End of code generation (xtrsm.cpp) */
