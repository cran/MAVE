/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * eig.cpp
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "eig.h"
#include "CVfast_emxutil.h"
#include "schur.h"
#include "relop.h"
#include "xzggev.h"
#include "CVfast_rtwutil.h"

/* Function Definitions */
void eig(const emxArray_real_T *A, emxArray_creal_T *V, emxArray_creal_T *D)
{
  int i8;
  boolean_T p;
  int kend;
  emxArray_creal_T *At;
  boolean_T exitg2;
  int i;
  int exitg1;
  emxArray_creal_T *alpha1;
  emxArray_creal_T *beta1;
  double absxk;
  int n;
  int c;
  int coltop;
  double colnorm;
  double scale;
  double t;
  double alpha1_re;
  double alpha1_im;
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
    i8 = V->size[0] * V->size[1];
    V->size[0] = A->size[0];
    V->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(creal_T));
    i = A->size[0] * A->size[1];
    for (i8 = 0; i8 < i; i8++) {
      V->data[i8].re = A->data[i8];
      V->data[i8].im = 0.0;
    }

    i8 = D->size[0] * D->size[1];
    D->size[0] = A->size[0];
    D->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)D, i8, (int)sizeof(creal_T));
    i = A->size[0] * A->size[1];
    for (i8 = 0; i8 < i; i8++) {
      D->data[i8].re = A->data[i8];
      D->data[i8].im = 0.0;
    }
  } else if ((A->size[0] == 1) && (A->size[1] == 1)) {
    i8 = V->size[0] * V->size[1];
    V->size[0] = 1;
    V->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(creal_T));
    for (i8 = 0; i8 < 1; i8++) {
      V->data[0].re = 1.0;
      V->data[0].im = 0.0;
    }

    i8 = D->size[0] * D->size[1];
    D->size[0] = A->size[0];
    D->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)D, i8, (int)sizeof(creal_T));
    i = A->size[0] * A->size[1];
    for (i8 = 0; i8 < i; i8++) {
      D->data[i8].re = A->data[i8];
      D->data[i8].im = 0.0;
    }
  } else {
    p = (A->size[0] == A->size[1]);
    if (p) {
      kend = 0;
      exitg2 = false;
      while ((!exitg2) && (kend <= A->size[1] - 1)) {
        i = 0;
        do {
          exitg1 = 0;
          if (i <= kend) {
            if (!(A->data[i + A->size[0] * kend] == A->data[kend + A->size[0] *
                  i])) {
              p = false;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            kend++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }

    emxInit_creal_T1(&At, 2);
    if (p) {
      schur(A, V, At);
      i8 = D->size[0] * D->size[1];
      D->size[0] = At->size[0];
      D->size[1] = At->size[1];
      emxEnsureCapacity((emxArray__common *)D, i8, (int)sizeof(creal_T));
      i = At->size[0] * At->size[1];
      for (i8 = 0; i8 < i; i8++) {
        D->data[i8] = At->data[i8];
      }

      absxk = At->data[0].re;
      D->data[0].re = absxk;
      D->data[0].im = 0.0;
      for (kend = 1; kend + 1 <= At->size[0]; kend++) {
        absxk = D->data[kend + D->size[0] * kend].re;
        D->data[kend + D->size[0] * kend].re = absxk;
        D->data[kend + D->size[0] * kend].im = 0.0;
        D->data[kend + D->size[0] * (kend - 1)].re = 0.0;
        D->data[kend + D->size[0] * (kend - 1)].im = 0.0;
        for (i = 1; i <= kend; i++) {
          D->data[(i + D->size[0] * kend) - 1].re = 0.0;
          D->data[(i + D->size[0] * kend) - 1].im = 0.0;
        }
      }
    } else {
      i8 = At->size[0] * At->size[1];
      At->size[0] = A->size[0];
      At->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)At, i8, (int)sizeof(creal_T));
      i = A->size[0] * A->size[1];
      for (i8 = 0; i8 < i; i8++) {
        At->data[i8].re = A->data[i8];
        At->data[i8].im = 0.0;
      }

      emxInit_creal_T(&alpha1, 1);
      emxInit_creal_T(&beta1, 1);
      xzggev(At, &i, alpha1, beta1, V);
      n = A->size[0];
      c = (A->size[0] - 1) * A->size[0];
      for (coltop = 0; coltop + 1 <= c + 1; coltop += n) {
        colnorm = 0.0;
        if (n == 1) {
          colnorm = rt_hypotd_snf(V->data[coltop].re, V->data[coltop].im);
        } else {
          scale = 2.2250738585072014E-308;
          kend = coltop + n;
          for (i = coltop; i + 1 <= kend; i++) {
            absxk = std::abs(V->data[i].re);
            if (absxk > scale) {
              t = scale / absxk;
              colnorm = 1.0 + colnorm * t * t;
              scale = absxk;
            } else {
              t = absxk / scale;
              colnorm += t * t;
            }

            absxk = std::abs(V->data[i].im);
            if (absxk > scale) {
              t = scale / absxk;
              colnorm = 1.0 + colnorm * t * t;
              scale = absxk;
            } else {
              t = absxk / scale;
              colnorm += t * t;
            }
          }

          colnorm = scale * std::sqrt(colnorm);
        }

        i8 = coltop + n;
        for (kend = coltop; kend + 1 <= i8; kend++) {
          absxk = V->data[kend].re;
          scale = V->data[kend].im;
          if (scale == 0.0) {
            V->data[kend].re = absxk / colnorm;
            V->data[kend].im = 0.0;
          } else if (absxk == 0.0) {
            V->data[kend].re = 0.0;
            V->data[kend].im = scale / colnorm;
          } else {
            V->data[kend].re = absxk / colnorm;
            V->data[kend].im = scale / colnorm;
          }
        }
      }

      i8 = D->size[0] * D->size[1];
      D->size[0] = alpha1->size[0];
      D->size[1] = alpha1->size[0];
      emxEnsureCapacity((emxArray__common *)D, i8, (int)sizeof(creal_T));
      i = alpha1->size[0] * alpha1->size[0];
      for (i8 = 0; i8 < i; i8++) {
        D->data[i8].re = 0.0;
        D->data[i8].im = 0.0;
      }

      for (i = 0; i < alpha1->size[0]; i++) {
        alpha1_re = alpha1->data[i].re;
        alpha1_im = alpha1->data[i].im;
        absxk = beta1->data[i].re;
        t = beta1->data[i].im;
        if (t == 0.0) {
          if (alpha1_im == 0.0) {
            D->data[i + D->size[0] * i].re = alpha1_re / absxk;
            D->data[i + D->size[0] * i].im = 0.0;
          } else if (alpha1_re == 0.0) {
            D->data[i + D->size[0] * i].re = 0.0;
            D->data[i + D->size[0] * i].im = alpha1_im / absxk;
          } else {
            D->data[i + D->size[0] * i].re = alpha1_re / absxk;
            D->data[i + D->size[0] * i].im = alpha1_im / absxk;
          }
        } else if (absxk == 0.0) {
          if (alpha1_re == 0.0) {
            D->data[i + D->size[0] * i].re = alpha1_im / t;
            D->data[i + D->size[0] * i].im = 0.0;
          } else if (alpha1_im == 0.0) {
            D->data[i + D->size[0] * i].re = 0.0;
            D->data[i + D->size[0] * i].im = -(alpha1_re / t);
          } else {
            D->data[i + D->size[0] * i].re = alpha1_im / t;
            D->data[i + D->size[0] * i].im = -(alpha1_re / t);
          }
        } else {
          colnorm = std::abs(absxk);
          scale = std::abs(t);
          if (colnorm > scale) {
            scale = t / absxk;
            absxk += scale * t;
            D->data[i + D->size[0] * i].re = (alpha1_re + scale * alpha1_im) /
              absxk;
            D->data[i + D->size[0] * i].im = (alpha1_im - scale * alpha1_re) /
              absxk;
          } else if (scale == colnorm) {
            if (absxk > 0.0) {
              absxk = 0.5;
            } else {
              absxk = -0.5;
            }

            if (t > 0.0) {
              scale = 0.5;
            } else {
              scale = -0.5;
            }

            D->data[i + D->size[0] * i].re = (alpha1_re * absxk + alpha1_im *
              scale) / colnorm;
            D->data[i + D->size[0] * i].im = (alpha1_im * absxk - alpha1_re *
              scale) / colnorm;
          } else {
            scale = absxk / t;
            absxk = t + scale * absxk;
            D->data[i + D->size[0] * i].re = (scale * alpha1_re + alpha1_im) /
              absxk;
            D->data[i + D->size[0] * i].im = (scale * alpha1_im - alpha1_re) /
              absxk;
          }
        }
      }

      emxFree_creal_T(&beta1);
      emxFree_creal_T(&alpha1);
    }

    emxFree_creal_T(&At);
  }
}

/* End of code generation (eig.cpp) */
