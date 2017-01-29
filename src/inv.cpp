/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * inv.cpp
 *
 * Code generation for function 'inv'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "inv.h"
#include "CVfast_emxutil.h"

/* Function Declarations */
static void invNxN(const emxArray_real_T *x, emxArray_real_T *y);

/* Function Definitions */
static void invNxN(const emxArray_real_T *x, emxArray_real_T *y)
{
  int n;
  int i6;
  int yk;
  emxArray_real_T *A;
  int b_n;
  emxArray_int32_T *ipiv;
  int k;
  emxArray_int32_T *p;
  int j;
  int mmj;
  int c;
  int ix;
  double smax;
  int i7;
  int jy;
  double s;
  int ijA;
  n = x->size[0];
  i6 = y->size[0] * y->size[1];
  y->size[0] = x->size[0];
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, i6, (int)sizeof(double));
  yk = x->size[0] * x->size[1];
  for (i6 = 0; i6 < yk; i6++) {
    y->data[i6] = 0.0;
  }

  emxInit_real_T(&A, 2);
  i6 = A->size[0] * A->size[1];
  A->size[0] = x->size[0];
  A->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)A, i6, (int)sizeof(double));
  yk = x->size[0] * x->size[1];
  for (i6 = 0; i6 < yk; i6++) {
    A->data[i6] = x->data[i6];
  }

  yk = x->size[0];
  if (yk < 1) {
    b_n = 0;
  } else {
    b_n = yk;
  }

  emxInit_int32_T1(&ipiv, 2);
  i6 = ipiv->size[0] * ipiv->size[1];
  ipiv->size[0] = 1;
  ipiv->size[1] = b_n;
  emxEnsureCapacity((emxArray__common *)ipiv, i6, (int)sizeof(int));
  if (b_n > 0) {
    ipiv->data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      ipiv->data[k - 1] = yk;
    }
  }

  if (x->size[0] < 1) {
    b_n = 0;
  } else {
    if (x->size[0] - 1 <= x->size[0]) {
      i6 = x->size[0] - 1;
    } else {
      i6 = x->size[0];
    }

    for (j = 0; j + 1 <= i6; j++) {
      mmj = n - j;
      c = j * (n + 1);
      if (mmj < 1) {
        yk = -1;
      } else {
        yk = 0;
        if (mmj > 1) {
          ix = c;
          smax = std::abs(A->data[c]);
          for (k = 1; k + 1 <= mmj; k++) {
            ix++;
            s = std::abs(A->data[ix]);
            if (s > smax) {
              yk = k;
              smax = s;
            }
          }
        }
      }

      if (A->data[c + yk] != 0.0) {
        if (yk != 0) {
          ipiv->data[j] = (j + yk) + 1;
          ix = j;
          yk += j;
          for (k = 1; k <= n; k++) {
            smax = A->data[ix];
            A->data[ix] = A->data[yk];
            A->data[yk] = smax;
            ix += n;
            yk += n;
          }
        }

        i7 = c + mmj;
        for (jy = c + 1; jy + 1 <= i7; jy++) {
          A->data[jy] /= A->data[c];
        }
      }

      yk = (n - j) - 1;
      b_n = c + n;
      jy = c + n;
      for (k = 1; k <= yk; k++) {
        smax = A->data[jy];
        if (A->data[jy] != 0.0) {
          ix = c + 1;
          i7 = mmj + b_n;
          for (ijA = 1 + b_n; ijA + 1 <= i7; ijA++) {
            A->data[ijA] += A->data[ix] * -smax;
            ix++;
          }
        }

        jy += n;
        b_n += n;
      }
    }

    b_n = x->size[0];
  }

  emxInit_int32_T1(&p, 2);
  i6 = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = b_n;
  emxEnsureCapacity((emxArray__common *)p, i6, (int)sizeof(int));
  if (b_n > 0) {
    p->data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      p->data[k - 1] = yk;
    }
  }

  for (k = 0; k < ipiv->size[1]; k++) {
    if (ipiv->data[k] > 1 + k) {
      yk = p->data[ipiv->data[k] - 1];
      p->data[ipiv->data[k] - 1] = p->data[k];
      p->data[k] = yk;
    }
  }

  emxFree_int32_T(&ipiv);
  for (k = 0; k + 1 <= n; k++) {
    c = p->data[k] - 1;
    y->data[k + y->size[0] * (p->data[k] - 1)] = 1.0;
    for (j = k; j + 1 <= n; j++) {
      if (y->data[j + y->size[0] * c] != 0.0) {
        for (jy = j + 1; jy + 1 <= n; jy++) {
          y->data[jy + y->size[0] * c] -= y->data[j + y->size[0] * c] * A->
            data[jy + A->size[0] * j];
        }
      }
    }
  }

  emxFree_int32_T(&p);
  if ((x->size[0] != 0) && (!((y->size[0] == 0) || (y->size[1] == 0)))) {
    for (j = 1; j <= n; j++) {
      yk = n * (j - 1);
      for (k = n - 1; k + 1 > 0; k--) {
        b_n = n * k;
        if (y->data[k + yk] != 0.0) {
          smax = y->data[k + yk];
          s = A->data[k + b_n];
          y->data[k + yk] = smax / s;
          for (jy = 0; jy + 1 <= k; jy++) {
            y->data[jy + yk] -= y->data[k + yk] * A->data[jy + b_n];
          }
        }
      }
    }
  }

  emxFree_real_T(&A);
}

void inv(const emxArray_real_T *x, emxArray_real_T *y)
{
  int i5;
  int loop_ub;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    i5 = y->size[0] * y->size[1];
    y->size[0] = x->size[0];
    y->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)y, i5, (int)sizeof(double));
    loop_ub = x->size[0] * x->size[1];
    for (i5 = 0; i5 < loop_ub; i5++) {
      y->data[i5] = x->data[i5];
    }
  } else {
    invNxN(x, y);
  }
}

/* End of code generation (inv.cpp) */
