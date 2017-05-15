/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * mldivide.cpp
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "mldivide.h"
#include "MAVEfast_emxutil.h"
#include "xtrsm.h"
#include "xzgetrf.h"
#include "xgeqp3.h"

/* Function Definitions */
void mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *Y)
{
  emxArray_real_T *b_A;
  unsigned int unnamed_idx_0;
  emxArray_int32_T *jpvt;
  int mn;
  int n;
  int minmn;
  emxArray_real_T *tau;
  int rankR;
  int maxmn;
  double tol;
  int j;
  int k;
  emxArray_real_T *b_B;
  int i;
  if (B->size[1] == 0) {
    unnamed_idx_0 = (unsigned int)A->size[1];
    mn = Y->size[0] * Y->size[1];
    Y->size[0] = (int)unnamed_idx_0;
    Y->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)Y, mn, (int)sizeof(double));
  } else {
    emxInit_real_T(&b_A, 2);
    emxInit_int32_T(&jpvt, 2);
    if (A->size[0] == A->size[1]) {
      n = A->size[1];
      mn = b_A->size[0] * b_A->size[1];
      b_A->size[0] = A->size[0];
      b_A->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)b_A, mn, (int)sizeof(double));
      minmn = A->size[0] * A->size[1];
      for (mn = 0; mn < minmn; mn++) {
        b_A->data[mn] = A->data[mn];
      }

      xzgetrf(A->size[1], A->size[1], b_A, A->size[1], jpvt, &minmn);
      rankR = B->size[1];
      mn = Y->size[0] * Y->size[1];
      Y->size[0] = B->size[0];
      Y->size[1] = B->size[1];
      emxEnsureCapacity((emxArray__common *)Y, mn, (int)sizeof(double));
      minmn = B->size[0] * B->size[1];
      for (mn = 0; mn < minmn; mn++) {
        Y->data[mn] = B->data[mn];
      }

      for (minmn = 0; minmn + 1 < n; minmn++) {
        if (jpvt->data[minmn] != minmn + 1) {
          maxmn = jpvt->data[minmn] - 1;
          for (mn = 0; mn + 1 <= rankR; mn++) {
            tol = Y->data[minmn + Y->size[0] * mn];
            Y->data[minmn + Y->size[0] * mn] = Y->data[maxmn + Y->size[0] * mn];
            Y->data[maxmn + Y->size[0] * mn] = tol;
          }
        }
      }

      for (j = 1; j <= rankR; j++) {
        minmn = n * (j - 1);
        for (k = 0; k + 1 <= n; k++) {
          maxmn = n * k;
          if (Y->data[k + minmn] != 0.0) {
            for (i = k + 1; i + 1 <= n; i++) {
              Y->data[i + minmn] -= Y->data[k + minmn] * b_A->data[i + maxmn];
            }
          }
        }
      }

      xtrsm(A->size[1], B->size[1], b_A, A->size[1], Y, A->size[1]);
    } else {
      mn = b_A->size[0] * b_A->size[1];
      b_A->size[0] = A->size[0];
      b_A->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)b_A, mn, (int)sizeof(double));
      minmn = A->size[0] * A->size[1];
      for (mn = 0; mn < minmn; mn++) {
        b_A->data[mn] = A->data[mn];
      }

      emxInit_real_T1(&tau, 1);
      xgeqp3(b_A, tau, jpvt);
      rankR = 0;
      if (b_A->size[0] < b_A->size[1]) {
        minmn = b_A->size[0];
        maxmn = b_A->size[1];
      } else {
        minmn = b_A->size[1];
        maxmn = b_A->size[0];
      }

      if (minmn > 0) {
        tol = (double)maxmn * std::abs(b_A->data[0]) * 2.2204460492503131E-16;
        while ((rankR < minmn) && (std::abs(b_A->data[rankR + b_A->size[0] *
                 rankR]) >= tol)) {
          rankR++;
        }
      }

      minmn = b_A->size[1];
      maxmn = B->size[1];
      mn = Y->size[0] * Y->size[1];
      Y->size[0] = minmn;
      Y->size[1] = maxmn;
      emxEnsureCapacity((emxArray__common *)Y, mn, (int)sizeof(double));
      minmn *= maxmn;
      for (mn = 0; mn < minmn; mn++) {
        Y->data[mn] = 0.0;
      }

      emxInit_real_T(&b_B, 2);
      mn = b_B->size[0] * b_B->size[1];
      b_B->size[0] = B->size[0];
      b_B->size[1] = B->size[1];
      emxEnsureCapacity((emxArray__common *)b_B, mn, (int)sizeof(double));
      minmn = B->size[0] * B->size[1];
      for (mn = 0; mn < minmn; mn++) {
        b_B->data[mn] = B->data[mn];
      }

      maxmn = b_A->size[0];
      minmn = b_A->size[0];
      mn = b_A->size[1];
      if (minmn <= mn) {
        mn = minmn;
      }

      for (j = 0; j + 1 <= mn; j++) {
        if (tau->data[j] != 0.0) {
          for (k = 0; k + 1 <= B->size[1]; k++) {
            tol = b_B->data[j + b_B->size[0] * k];
            for (i = j + 1; i + 1 <= maxmn; i++) {
              tol += b_A->data[i + b_A->size[0] * j] * b_B->data[i + b_B->size[0]
                * k];
            }

            tol *= tau->data[j];
            if (tol != 0.0) {
              b_B->data[j + b_B->size[0] * k] -= tol;
              for (i = j + 1; i + 1 <= maxmn; i++) {
                b_B->data[i + b_B->size[0] * k] -= b_A->data[i + b_A->size[0] *
                  j] * tol;
              }
            }
          }
        }
      }

      emxFree_real_T(&tau);
      for (k = 0; k + 1 <= B->size[1]; k++) {
        for (i = 0; i + 1 <= rankR; i++) {
          Y->data[(jpvt->data[i] + Y->size[0] * k) - 1] = b_B->data[i +
            b_B->size[0] * k];
        }

        for (j = rankR - 1; j + 1 > 0; j--) {
          Y->data[(jpvt->data[j] + Y->size[0] * k) - 1] /= b_A->data[j +
            b_A->size[0] * j];
          for (i = 0; i + 1 <= j; i++) {
            Y->data[(jpvt->data[i] + Y->size[0] * k) - 1] -= Y->data[(jpvt->
              data[j] + Y->size[0] * k) - 1] * b_A->data[i + b_A->size[0] * j];
          }
        }
      }

      emxFree_real_T(&b_B);
    }

    emxFree_int32_T(&jpvt);
    emxFree_real_T(&b_A);
  }
}

/* End of code generation (mldivide.cpp) */
