/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * schur.cpp
 *
 * Code generation for function 'schur'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "schur.h"
#include "xscal.h"
#include "xzlarf.h"
#include "xzlarfg.h"
#include "xdlanv2.h"
#include "MAVEfast_emxutil.h"
#include "xdhseqr.h"
#include "xgehrd.h"
#include "MAVEfast_rtwutil.h"

/* Function Definitions */
void schur(const emxArray_real_T *A, emxArray_creal_T *V, emxArray_creal_T *T)
{
  emxArray_real_T *b_A;
  int n;
  int iajm1;
  int ia;
  emxArray_real_T *Vr;
  emxArray_real_T *tau;
  int nh;
  int j;
  emxArray_real_T *h;
  int i;
  int iaii;
  emxArray_real_T *work;
  int itau;
  int varargin_1[2];
  double r;
  double b;
  double c;
  double d;
  double s;
  double rt1i;
  double t1_re;
  double t1_im;
  double mu1_im;
  double mu1_re;
  emxInit_real_T(&b_A, 2);
  n = A->size[0];
  iajm1 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)b_A, iajm1, (int)sizeof(double));
  ia = A->size[0] * A->size[1];
  for (iajm1 = 0; iajm1 < ia; iajm1++) {
    b_A->data[iajm1] = A->data[iajm1];
  }

  emxInit_real_T(&Vr, 2);
  emxInit_real_T1(&tau, 1);
  xgehrd(b_A, tau);
  iajm1 = Vr->size[0] * Vr->size[1];
  Vr->size[0] = b_A->size[0];
  Vr->size[1] = b_A->size[1];
  emxEnsureCapacity((emxArray__common *)Vr, iajm1, (int)sizeof(double));
  ia = b_A->size[0] * b_A->size[1];
  for (iajm1 = 0; iajm1 < ia; iajm1++) {
    Vr->data[iajm1] = b_A->data[iajm1];
  }

  if (A->size[0] != 0) {
    nh = A->size[0] - 1;
    for (j = A->size[0]; j > 1; j--) {
      ia = (j - 1) * n;
      for (i = 1; i < j; i++) {
        Vr->data[(ia + i) - 1] = 0.0;
      }

      iajm1 = ia - n;
      for (i = j; i + 1 <= n; i++) {
        Vr->data[ia + i] = Vr->data[iajm1 + i];
      }

      for (i = n; i + 1 <= n; i++) {
        Vr->data[ia + i] = 0.0;
      }
    }

    for (i = 1; i <= n; i++) {
      Vr->data[i - 1] = 0.0;
    }

    Vr->data[0] = 1.0;
    for (j = A->size[0]; j + 1 <= n; j++) {
      ia = j * n;
      for (i = 1; i <= n; i++) {
        Vr->data[(ia + i) - 1] = 0.0;
      }

      Vr->data[ia + j] = 1.0;
    }

    if (!(A->size[0] - 1 < 1)) {
      for (j = A->size[0]; j - 1 < nh; j++) {
        ia = n + (j - 1) * n;
        for (i = 0; i < nh; i++) {
          Vr->data[(ia + i) + 1] = 0.0;
        }

        Vr->data[ia + j] = 1.0;
      }

      emxInit_real_T1(&work, 1);
      itau = A->size[0] - 2;
      varargin_1[0] = Vr->size[1];
      iajm1 = work->size[0];
      work->size[0] = varargin_1[0];
      emxEnsureCapacity((emxArray__common *)work, iajm1, (int)sizeof(double));
      ia = varargin_1[0];
      for (iajm1 = 0; iajm1 < ia; iajm1++) {
        work->data[iajm1] = 0.0;
      }

      for (i = A->size[0] - 1; i >= 1; i--) {
        iaii = ((A->size[0] + i) + (i - 1) * n) + 1;
        if (i < n - 1) {
          Vr->data[iaii - 1] = 1.0;
          xzlarf(n - i, nh - i, iaii, tau->data[itau], Vr, iaii + n, n, work);
          xscal(nh - i, -tau->data[itau], Vr, iaii + 1);
        }

        Vr->data[iaii - 1] = 1.0 - tau->data[itau];
        for (j = 1; j < i; j++) {
          Vr->data[(iaii - j) - 1] = 0.0;
        }

        itau--;
      }

      emxFree_real_T(&work);
    }
  }

  emxFree_real_T(&tau);
  emxInit_real_T(&h, 2);
  eml_dlahqr(b_A, Vr);
  iajm1 = h->size[0] * h->size[1];
  h->size[0] = b_A->size[0];
  h->size[1] = b_A->size[1];
  emxEnsureCapacity((emxArray__common *)h, iajm1, (int)sizeof(double));
  ia = b_A->size[0] * b_A->size[1];
  for (iajm1 = 0; iajm1 < ia; iajm1++) {
    h->data[iajm1] = b_A->data[iajm1];
  }

  if ((b_A->size[0] == 0) || (b_A->size[1] == 0) || (3 >= b_A->size[0])) {
  } else {
    ia = 4;
    if (b_A->size[0] - 4 < b_A->size[1] - 1) {
      iaii = b_A->size[0] - 3;
    } else {
      iaii = b_A->size[1];
    }

    for (j = 1; j <= iaii; j++) {
      for (i = ia; i <= b_A->size[0]; i++) {
        h->data[(i + h->size[0] * (j - 1)) - 1] = 0.0;
      }

      ia++;
    }
  }

  emxFree_real_T(&b_A);
  iajm1 = T->size[0] * T->size[1];
  T->size[0] = h->size[0];
  T->size[1] = h->size[1];
  emxEnsureCapacity((emxArray__common *)T, iajm1, (int)sizeof(creal_T));
  ia = h->size[0] * h->size[1];
  for (iajm1 = 0; iajm1 < ia; iajm1++) {
    T->data[iajm1].re = h->data[iajm1];
    T->data[iajm1].im = 0.0;
  }

  iajm1 = V->size[0] * V->size[1];
  V->size[0] = Vr->size[0];
  V->size[1] = Vr->size[1];
  emxEnsureCapacity((emxArray__common *)V, iajm1, (int)sizeof(creal_T));
  ia = Vr->size[0] * Vr->size[1];
  for (iajm1 = 0; iajm1 < ia; iajm1++) {
    V->data[iajm1].re = Vr->data[iajm1];
    V->data[iajm1].im = 0.0;
  }

  for (iajm1 = 0; iajm1 < 2; iajm1++) {
    varargin_1[iajm1] = h->size[iajm1];
  }

  iaii = varargin_1[0];
  if (varargin_1[1] < varargin_1[0]) {
    iaii = varargin_1[1];
  }

  for (iajm1 = 0; iajm1 < 2; iajm1++) {
    varargin_1[iajm1] = Vr->size[iajm1];
  }

  emxFree_real_T(&Vr);
  ia = varargin_1[0];
  if (varargin_1[1] < varargin_1[0]) {
    ia = varargin_1[1];
  }

  if (iaii <= ia) {
    ia = iaii;
  }

  if (ia != 0) {
    for (iaii = ia - 1; iaii + 1 >= 2; iaii--) {
      if (h->data[iaii + h->size[0] * (iaii - 1)] != 0.0) {
        r = h->data[(iaii + h->size[0] * (iaii - 1)) - 1];
        b = h->data[(iaii + h->size[0] * iaii) - 1];
        c = h->data[iaii + h->size[0] * (iaii - 1)];
        d = h->data[iaii + h->size[0] * iaii];
        xdlanv2(&r, &b, &c, &d, &s, &rt1i, &t1_re, &t1_im, &mu1_im, &mu1_re);
        mu1_re = s - h->data[iaii + h->size[0] * iaii];
        r = rt_hypotd_snf(rt_hypotd_snf(mu1_re, rt1i), h->data[iaii + h->size[0]
                          * (iaii - 1)]);
        if (rt1i == 0.0) {
          mu1_re /= r;
          mu1_im = 0.0;
        } else if (mu1_re == 0.0) {
          mu1_re = 0.0;
          mu1_im = rt1i / r;
        } else {
          mu1_re /= r;
          mu1_im = rt1i / r;
        }

        s = h->data[iaii + h->size[0] * (iaii - 1)] / r;
        for (j = iaii - 1; j + 1 <= ia; j++) {
          t1_re = T->data[(iaii + T->size[0] * j) - 1].re;
          t1_im = T->data[(iaii + T->size[0] * j) - 1].im;
          c = T->data[(iaii + T->size[0] * j) - 1].re;
          d = T->data[(iaii + T->size[0] * j) - 1].im;
          r = T->data[(iaii + T->size[0] * j) - 1].im;
          b = T->data[(iaii + T->size[0] * j) - 1].re;
          T->data[(iaii + T->size[0] * j) - 1].re = (mu1_re * c + mu1_im * d) +
            s * T->data[iaii + T->size[0] * j].re;
          T->data[(iaii + T->size[0] * j) - 1].im = (mu1_re * r - mu1_im * b) +
            s * T->data[iaii + T->size[0] * j].im;
          r = mu1_re * T->data[iaii + T->size[0] * j].re - mu1_im * T->data[iaii
            + T->size[0] * j].im;
          b = mu1_re * T->data[iaii + T->size[0] * j].im + mu1_im * T->data[iaii
            + T->size[0] * j].re;
          T->data[iaii + T->size[0] * j].re = r - s * t1_re;
          T->data[iaii + T->size[0] * j].im = b - s * t1_im;
        }

        for (i = 0; i + 1 <= iaii + 1; i++) {
          t1_re = T->data[i + T->size[0] * (iaii - 1)].re;
          t1_im = T->data[i + T->size[0] * (iaii - 1)].im;
          r = mu1_re * T->data[i + T->size[0] * (iaii - 1)].re - mu1_im *
            T->data[i + T->size[0] * (iaii - 1)].im;
          b = mu1_re * T->data[i + T->size[0] * (iaii - 1)].im + mu1_im *
            T->data[i + T->size[0] * (iaii - 1)].re;
          c = T->data[i + T->size[0] * iaii].re;
          d = T->data[i + T->size[0] * iaii].im;
          T->data[i + T->size[0] * (iaii - 1)].re = r + s * c;
          T->data[i + T->size[0] * (iaii - 1)].im = b + s * d;
          c = T->data[i + T->size[0] * iaii].re;
          d = T->data[i + T->size[0] * iaii].im;
          r = T->data[i + T->size[0] * iaii].im;
          b = T->data[i + T->size[0] * iaii].re;
          T->data[i + T->size[0] * iaii].re = (mu1_re * c + mu1_im * d) - s *
            t1_re;
          T->data[i + T->size[0] * iaii].im = (mu1_re * r - mu1_im * b) - s *
            t1_im;
        }

        for (i = 0; i + 1 <= ia; i++) {
          t1_re = V->data[i + V->size[0] * (iaii - 1)].re;
          t1_im = V->data[i + V->size[0] * (iaii - 1)].im;
          r = mu1_re * V->data[i + V->size[0] * (iaii - 1)].re - mu1_im *
            V->data[i + V->size[0] * (iaii - 1)].im;
          b = mu1_re * V->data[i + V->size[0] * (iaii - 1)].im + mu1_im *
            V->data[i + V->size[0] * (iaii - 1)].re;
          c = V->data[i + V->size[0] * iaii].re;
          d = V->data[i + V->size[0] * iaii].im;
          V->data[i + V->size[0] * (iaii - 1)].re = r + s * c;
          V->data[i + V->size[0] * (iaii - 1)].im = b + s * d;
          c = V->data[i + V->size[0] * iaii].re;
          d = V->data[i + V->size[0] * iaii].im;
          r = V->data[i + V->size[0] * iaii].im;
          b = V->data[i + V->size[0] * iaii].re;
          V->data[i + V->size[0] * iaii].re = (mu1_re * c + mu1_im * d) - s *
            t1_re;
          V->data[i + V->size[0] * iaii].im = (mu1_re * r - mu1_im * b) - s *
            t1_im;
        }

        T->data[iaii + T->size[0] * (iaii - 1)].re = 0.0;
        T->data[iaii + T->size[0] * (iaii - 1)].im = 0.0;
      }
    }
  }

  emxFree_real_T(&h);
}

/* End of code generation (schur.cpp) */
