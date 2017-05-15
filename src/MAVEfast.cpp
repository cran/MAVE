/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * MAVEfast.cpp
 *
 * Code generation for function 'MAVEfast'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "MAVEfast_emxutil.h"
#include "power.h"
#include "mldivide.h"
#include "eye.h"
#include "repmat.h"
#include "exp.h"
#include "sort1.h"
#include "diag.h"
#include "eig.h"
#include "kron.h"
#include "sum.h"
#include "inv.h"
#include "strcmp.h"
#include "mean.h"
#include "std.h"
#include "norm.h"
#include "rdivide.h"
#include "sqrt.h"
#include "abs.h"
#include "quantile.h"
#include "upper.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);
static void unifD(const emxArray_real_T *x, double m, emxArray_real_T *I);

/* Function Definitions */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = rtNaN;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

static void unifD(const emxArray_real_T *x, double m, emxArray_real_T *I)
{
  int p;
  emxArray_real_T *y;
  double k;
  int b_k;
  int i;
  int loop_ub;
  emxArray_real_T *varargin_2;
  int xoffset;
  boolean_T empty_non_axis_sizes;
  int nrows;
  int result;
  emxArray_real_T *Xremain;
  emxArray_int32_T *Ii;
  emxArray_int32_T *b_Ii;
  emxArray_int32_T *iidx;
  emxArray_real_T *b_x;
  emxArray_boolean_T *b;
  emxArray_real_T *b_Xremain;
  emxArray_real_T *b_y;
  emxArray_real_T *c_Xremain;
  emxArray_real_T *d_Xremain;
  emxArray_real_T *e_Xremain;
  emxArray_real_T *b_I;
  emxArray_real_T *f_Xremain;
  double s;
  int nrowx;
  int ncolx;

  /*  m is the number of points where teh gradients are calculated */
  m = std::floor(m);
  p = x->size[1];
  emxInit_real_T(&y, 2);
  if (std::ceil((double)x->size[0] / m) == 1.0) {
    if (x->size[0] < 1) {
      b_k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)y, b_k, (int)sizeof(double));
    } else {
      b_k = x->size[0];
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)((double)b_k - 1.0) + 1;
      emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
      loop_ub = (int)((double)b_k - 1.0);
      for (b_k = 0; b_k <= loop_ub; b_k++) {
        y->data[y->size[0] * b_k] = 1.0 + (double)b_k;
      }
    }

    b_k = I->size[0] * I->size[1];
    I->size[0] = 1;
    I->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)I, b_k, (int)sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (b_k = 0; b_k < loop_ub; b_k++) {
      I->data[b_k] = y->data[b_k];
    }
  } else {
    k = std::floor((double)x->size[0] / m);
    if (x->size[0] < 1) {
      b_k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)y, b_k, (int)sizeof(double));
    } else {
      b_k = x->size[0];
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)((double)b_k - 1.0) + 1;
      emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
      loop_ub = (int)((double)b_k - 1.0);
      for (b_k = 0; b_k <= loop_ub; b_k++) {
        y->data[y->size[0] * b_k] = 1.0 + (double)b_k;
      }
    }

    emxInit_real_T1(&varargin_2, 1);
    b_k = varargin_2->size[0];
    varargin_2->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)varargin_2, b_k, (int)sizeof(double));
    loop_ub = y->size[1];
    for (b_k = 0; b_k < loop_ub; b_k++) {
      varargin_2->data[b_k] = y->data[y->size[0] * b_k];
    }

    if (!((x->size[0] == 0) || (x->size[1] == 0))) {
      xoffset = x->size[0];
    } else if (!(varargin_2->size[0] == 0)) {
      xoffset = varargin_2->size[0];
    } else {
      xoffset = x->size[0];
      if (!(xoffset >= 0)) {
        xoffset = 0;
      }
    }

    empty_non_axis_sizes = (xoffset == 0);
    if (empty_non_axis_sizes || (!((x->size[0] == 0) || (x->size[1] == 0)))) {
      nrows = x->size[1];
    } else {
      nrows = 0;
    }

    if (empty_non_axis_sizes || (!(varargin_2->size[0] == 0))) {
      result = 1;
    } else {
      result = 0;
    }

    emxInit_real_T(&Xremain, 2);
    b_k = Xremain->size[0] * Xremain->size[1];
    Xremain->size[0] = xoffset;
    Xremain->size[1] = nrows + result;
    emxEnsureCapacity((emxArray__common *)Xremain, b_k, (int)sizeof(double));
    for (b_k = 0; b_k < nrows; b_k++) {
      for (i = 0; i < xoffset; i++) {
        Xremain->data[i + Xremain->size[0] * b_k] = x->data[i + xoffset * b_k];
      }
    }

    for (b_k = 0; b_k < result; b_k++) {
      for (i = 0; i < xoffset; i++) {
        Xremain->data[i + Xremain->size[0] * (b_k + nrows)] = varargin_2->data[i
          + xoffset * b_k];
      }
    }

    b_k = I->size[0] * I->size[1];
    I->size[0] = 0;
    I->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)I, b_k, (int)sizeof(double));
    emxInit_int32_T1(&Ii, 1);
    emxInit_int32_T1(&b_Ii, 1);
    emxInit_int32_T1(&iidx, 1);
    emxInit_real_T(&b_x, 2);
    emxInit_boolean_T(&b, 2);
    emxInit_real_T(&b_Xremain, 2);
    emxInit_real_T(&b_y, 2);
    emxInit_real_T(&c_Xremain, 2);
    emxInit_real_T(&d_Xremain, 2);
    emxInit_real_T(&e_Xremain, 2);
    emxInit_real_T(&b_I, 2);
    emxInit_real_T(&f_Xremain, 2);
    while (Xremain->size[0] > 1) {
      if (1 > p) {
        loop_ub = 0;
        result = 0;
      } else {
        loop_ub = p;
        result = p;
      }

      b_k = d_Xremain->size[0] * d_Xremain->size[1];
      d_Xremain->size[0] = 1;
      d_Xremain->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)d_Xremain, b_k, (int)sizeof(double));
      for (b_k = 0; b_k < loop_ub; b_k++) {
        d_Xremain->data[d_Xremain->size[0] * b_k] = Xremain->data[Xremain->size
          [0] * b_k];
      }

      repmat(d_Xremain, (double)Xremain->size[0], b_x);
      loop_ub = Xremain->size[0];
      b_k = c_Xremain->size[0] * c_Xremain->size[1];
      c_Xremain->size[0] = loop_ub;
      c_Xremain->size[1] = result;
      emxEnsureCapacity((emxArray__common *)c_Xremain, b_k, (int)sizeof(double));
      for (b_k = 0; b_k < result; b_k++) {
        for (i = 0; i < loop_ub; i++) {
          c_Xremain->data[i + c_Xremain->size[0] * b_k] = Xremain->data[i +
            Xremain->size[0] * b_k] - b_x->data[i + b_x->size[0] * b_k];
        }
      }

      power(c_Xremain, b_x);
      b_sum(b_x, varargin_2);
      e_sort(varargin_2, iidx);
      b_k = b_Ii->size[0];
      b_Ii->size[0] = iidx->size[0];
      emxEnsureCapacity((emxArray__common *)b_Ii, b_k, (int)sizeof(int));
      loop_ub = iidx->size[0];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        b_Ii->data[b_k] = iidx->data[b_k];
      }

      if (k <= Xremain->size[0]) {
        s = k;
      } else {
        s = Xremain->size[0];
      }

      if (1.0 > s) {
        loop_ub = 0;
      } else {
        loop_ub = (int)s;
      }

      b_k = Ii->size[0];
      Ii->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)Ii, b_k, (int)sizeof(int));
      for (b_k = 0; b_k < loop_ub; b_k++) {
        Ii->data[b_k] = b_Ii->data[b_k];
      }

      if (1 > p) {
        result = 0;
      } else {
        result = p;
      }

      b_k = b_x->size[0] * b_x->size[1];
      b_x->size[0] = Ii->size[0];
      b_x->size[1] = result;
      emxEnsureCapacity((emxArray__common *)b_x, b_k, (int)sizeof(double));
      for (b_k = 0; b_k < result; b_k++) {
        nrows = Ii->size[0];
        for (i = 0; i < nrows; i++) {
          b_x->data[i + b_x->size[0] * b_k] = Xremain->data[(Ii->data[i] +
            Xremain->size[0] * b_k) - 1];
        }
      }

      b_k = e_Xremain->size[0] * e_Xremain->size[1];
      e_Xremain->size[0] = loop_ub;
      e_Xremain->size[1] = result;
      emxEnsureCapacity((emxArray__common *)e_Xremain, b_k, (int)sizeof(double));
      for (b_k = 0; b_k < result; b_k++) {
        for (i = 0; i < loop_ub; i++) {
          e_Xremain->data[i + e_Xremain->size[0] * b_k] = Xremain->data
            [(b_Ii->data[i] + Xremain->size[0] * b_k) - 1];
        }
      }

      b_k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = e_Xremain->size[1];
      emxEnsureCapacity((emxArray__common *)y, b_k, (int)sizeof(double));
      if ((loop_ub == 0) || (result == 0)) {
        b_k = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, b_k, (int)sizeof(double));
        result = y->size[1];
        for (b_k = 0; b_k < result; b_k++) {
          y->data[y->size[0] * b_k] = 0.0;
        }
      } else {
        for (i = 0; i + 1 <= result; i++) {
          xoffset = i * loop_ub;
          s = b_x->data[xoffset];
          for (b_k = 2; b_k <= loop_ub; b_k++) {
            s += b_x->data[(xoffset + b_k) - 1];
          }

          y->data[i] = s;
        }
      }

      if (1 > p) {
        result = 0;
      } else {
        result = p;
      }

      b_k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, b_k, (int)sizeof(double));
      nrows = y->size[0] * y->size[1];
      for (b_k = 0; b_k < nrows; b_k++) {
        b_y->data[b_k] = y->data[b_k] / (double)loop_ub;
      }

      repmat(b_y, (double)loop_ub, b_x);
      b_k = b_Xremain->size[0] * b_Xremain->size[1];
      b_Xremain->size[0] = Ii->size[0];
      b_Xremain->size[1] = result;
      emxEnsureCapacity((emxArray__common *)b_Xremain, b_k, (int)sizeof(double));
      for (b_k = 0; b_k < result; b_k++) {
        loop_ub = Ii->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_Xremain->data[i + b_Xremain->size[0] * b_k] = Xremain->data
            [(Ii->data[i] + Xremain->size[0] * b_k) - 1] - b_x->data[i +
            b_x->size[0] * b_k];
        }
      }

      power(b_Xremain, b_x);
      b_sum(b_x, varargin_2);
      e_sort(varargin_2, iidx);
      b_k = varargin_2->size[0];
      varargin_2->size[0] = iidx->size[0];
      emxEnsureCapacity((emxArray__common *)varargin_2, b_k, (int)sizeof(double));
      loop_ub = iidx->size[0];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        varargin_2->data[b_k] = iidx->data[b_k];
      }

      s = Xremain->data[(b_Ii->data[(int)varargin_2->data[0] - 1] +
                         Xremain->size[0] * p) - 1];
      if (!((I->size[0] == 0) || (I->size[1] == 0))) {
        nrows = I->size[1];
      } else {
        nrows = 0;
      }

      b_k = b_I->size[0] * b_I->size[1];
      b_I->size[0] = 1;
      b_I->size[1] = nrows + 1;
      emxEnsureCapacity((emxArray__common *)b_I, b_k, (int)sizeof(double));
      for (b_k = 0; b_k < nrows; b_k++) {
        b_I->data[b_I->size[0] * b_k] = I->data[b_k];
      }

      b_I->data[b_I->size[0] * nrows] = s;
      b_k = I->size[0] * I->size[1];
      I->size[0] = b_I->size[0];
      I->size[1] = b_I->size[1];
      emxEnsureCapacity((emxArray__common *)I, b_k, (int)sizeof(double));
      loop_ub = b_I->size[1];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        result = b_I->size[0];
        for (i = 0; i < result; i++) {
          I->data[i + I->size[0] * b_k] = b_I->data[i + b_I->size[0] * b_k];
        }
      }

      b_k = iidx->size[0];
      iidx->size[0] = Ii->size[0];
      emxEnsureCapacity((emxArray__common *)iidx, b_k, (int)sizeof(int));
      loop_ub = Ii->size[0];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        iidx->data[b_k] = Ii->data[b_k];
      }

      nrowx = Xremain->size[0];
      ncolx = Xremain->size[1];
      if (iidx->size[0] == 1) {
        nrows = Xremain->size[0] - 1;
        for (result = 0; result + 1 <= ncolx; result++) {
          for (i = iidx->data[0]; i < nrowx; i++) {
            Xremain->data[(i + Xremain->size[0] * result) - 1] = Xremain->data[i
              + Xremain->size[0] * result];
          }
        }
      } else {
        b_k = b->size[0] * b->size[1];
        b->size[0] = 1;
        b->size[1] = Xremain->size[0];
        emxEnsureCapacity((emxArray__common *)b, b_k, (int)sizeof(boolean_T));
        loop_ub = Xremain->size[0];
        for (b_k = 0; b_k < loop_ub; b_k++) {
          b->data[b_k] = false;
        }

        for (b_k = 1; b_k <= iidx->size[0]; b_k++) {
          b->data[iidx->data[b_k - 1] - 1] = true;
        }

        xoffset = 0;
        for (b_k = 1; b_k <= b->size[1]; b_k++) {
          xoffset += b->data[b_k - 1];
        }

        nrows = Xremain->size[0] - xoffset;
        i = 0;
        for (b_k = 1; b_k <= nrowx; b_k++) {
          if ((b_k > b->size[1]) || (!b->data[b_k - 1])) {
            for (result = 0; result + 1 <= ncolx; result++) {
              Xremain->data[i + Xremain->size[0] * result] = Xremain->data[(b_k
                + Xremain->size[0] * result) - 1];
            }

            i++;
          }
        }
      }

      if (1 > nrows) {
        loop_ub = 0;
      } else {
        loop_ub = nrows;
      }

      xoffset = Xremain->size[1];
      b_k = f_Xremain->size[0] * f_Xremain->size[1];
      f_Xremain->size[0] = loop_ub;
      f_Xremain->size[1] = xoffset;
      emxEnsureCapacity((emxArray__common *)f_Xremain, b_k, (int)sizeof(double));
      for (b_k = 0; b_k < xoffset; b_k++) {
        for (i = 0; i < loop_ub; i++) {
          f_Xremain->data[i + f_Xremain->size[0] * b_k] = Xremain->data[i +
            Xremain->size[0] * b_k];
        }
      }

      b_k = Xremain->size[0] * Xremain->size[1];
      Xremain->size[0] = f_Xremain->size[0];
      Xremain->size[1] = f_Xremain->size[1];
      emxEnsureCapacity((emxArray__common *)Xremain, b_k, (int)sizeof(double));
      loop_ub = f_Xremain->size[1];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        result = f_Xremain->size[0];
        for (i = 0; i < result; i++) {
          Xremain->data[i + Xremain->size[0] * b_k] = f_Xremain->data[i +
            f_Xremain->size[0] * b_k];
        }
      }
    }

    emxFree_real_T(&f_Xremain);
    emxFree_real_T(&b_I);
    emxFree_real_T(&e_Xremain);
    emxFree_real_T(&d_Xremain);
    emxFree_real_T(&c_Xremain);
    emxFree_real_T(&b_y);
    emxFree_real_T(&b_Xremain);
    emxFree_boolean_T(&b);
    emxFree_real_T(&b_x);
    emxFree_int32_T(&iidx);
    emxFree_real_T(&varargin_2);
    emxFree_int32_T(&b_Ii);
    emxFree_int32_T(&Ii);
    emxFree_real_T(&Xremain);
  }

  emxFree_real_T(&y);
}

void MAVEfast(emxArray_real_T *x, const emxArray_real_T *y, emxArray_char_T
              *method, emxArray_real_T *which_dim, emxArray_real_T *BB1D,
              emxArray_real_T *ky)
{
  emxArray_real_T *qx;
  emxArray_real_T *a;
  int p;
  int n;
  int i8;
  int loop_ub;
  int br;
  int i9;
  emxArray_real_T *yi;
  int k;
  int nx;
  int vstride;
  int m;
  emxArray_real_T *b_yi;
  double pp;
  int k1;
  emxArray_real_T *ss;
  int ic;
  emxArray_real_T *V;
  emxArray_creal_T *Vc;
  emxArray_creal_T *Dc;
  int ar;
  int ib;
  int ia;
  emxArray_real_T *b;
  emxArray_real_T *yj;
  emxArray_real_T *b_x;
  emxArray_char_T *b_method;
  emxArray_real_T *b_p;
  emxArray_real_T *ky1;
  emxArray_real_T *ky2;
  emxArray_real_T *DD;
  emxArray_real_T *b_qx;
  emxArray_real_T *C;
  emxArray_real_T *U;
  emxArray_real_T *s;
  emxArray_real_T *dc;
  emxArray_int32_T *iidx;
  emxArray_boolean_T *c_x;
  emxArray_int32_T *ii;
  emxArray_real_T *b_s;
  emxArray_real_T *b_Dc;
  emxArray_boolean_T *c_yi;
  emxArray_real_T *c_qx;
  emxArray_real_T *b_ky;
  emxArray_real_T *b_V;
  emxArray_real_T *c_V;
  int nm1d2;
  emxArray_real_T *BB;
  emxArray_real_T *r0;
  emxArray_real_T *onexi;
  emxArray_real_T *B;
  int iter;
  emxArray_real_T *D;
  emxArray_real_T *Ifast;
  emxArray_real_T *b_B;
  emxArray_real_T *xfast;
  emxArray_real_T *c_B;
  emxArray_real_T *xij;
  boolean_T empty_non_axis_sizes;
  emxArray_real_T *dxij;
  emxArray_real_T *abi;
  emxArray_real_T *dd;
  emxArray_real_T *kxijy;
  emxArray_real_T *ddx;
  emxArray_real_T *tmp;
  emxArray_real_T *B0;
  emxArray_real_T *d_B;
  emxArray_creal_T *R;
  emxArray_real_T *b_C;
  boolean_T guard1 = false;
  emxArray_real_T *c_C;
  emxArray_real_T *d_C;
  emxArray_real_T *e_C;
  emxArray_real_T *d_x;
  emxArray_int32_T *r1;
  emxArray_real_T *b_y;
  emxArray_real_T *c_y;
  emxArray_real_T *d_y;
  cell_wrap_0 reshapes[2];
  emxArray_real_T *e_y;
  emxArray_real_T *f_y;
  emxArray_real_T *g_y;
  emxArray_real_T *b_dd;
  emxArray_real_T *b_abi;
  emxArray_real_T *f_C;
  emxArray_real_T *d_V;
  emxArray_real_T *b_U;
  emxArray_real_T *b_xfast;
  emxArray_real_T *g_C;
  emxArray_real_T *e_V;
  emxArray_real_T *e_B;
  emxArray_int32_T *b_Ifast;
  emxArray_int32_T *c_Ifast;
  emxArray_real_T *b_dc;
  emxArray_real_T *c_abi;
  emxArray_int32_T *d_Ifast;
  emxArray_int32_T *e_Ifast;
  boolean_T exitg1;
  double ip;
  double b_m;
  int K;
  int b_iter;
  int j;
  double b_nm1d2;
  boolean_T exitg2;
  boolean_T b_guard1 = false;
  int Ifast_idx_0;
  double h2;
  unsigned int a_idx_0;
  double d;
  double ndbl;
  double apnd;
  double cdiff;
  double absa;
  double absb;
  emxInit_real_T(&qx, 2);
  emxInit_real_T(&a, 2);

  /*  */
  p = x->size[1];
  n = x->size[0];
  mean(x, qx);
  repmat(qx, (double)x->size[0], a);
  i8 = x->size[0] * x->size[1];
  emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
  loop_ub = x->size[1];
  for (i8 = 0; i8 < loop_ub; i8++) {
    br = x->size[0];
    for (i9 = 0; i9 < br; i9++) {
      x->data[i9 + x->size[0] * i8] -= a->data[i9 + a->size[0] * i8];
    }
  }

  i8 = a->size[0] * a->size[1];
  a->size[0] = x->size[1];
  a->size[1] = x->size[0];
  emxEnsureCapacity((emxArray__common *)a, i8, (int)sizeof(double));
  loop_ub = x->size[0];
  for (i8 = 0; i8 < loop_ub; i8++) {
    br = x->size[1];
    for (i9 = 0; i9 < br; i9++) {
      a->data[i9 + a->size[0] * i8] = x->data[i8 + x->size[0] * i9];
    }
  }

  emxInit_real_T(&yi, 2);
  if ((a->size[1] == 1) || (x->size[0] == 1)) {
    i8 = yi->size[0] * yi->size[1];
    yi->size[0] = a->size[0];
    yi->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
    loop_ub = a->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = x->size[1];
      for (i9 = 0; i9 < br; i9++) {
        yi->data[i8 + yi->size[0] * i9] = 0.0;
        vstride = a->size[1];
        for (k1 = 0; k1 < vstride; k1++) {
          yi->data[i8 + yi->size[0] * i9] += a->data[i8 + a->size[0] * k1] *
            x->data[k1 + x->size[0] * i9];
        }
      }
    }
  } else {
    k = a->size[1];
    nx = a->size[0];
    vstride = x->size[1];
    i8 = yi->size[0] * yi->size[1];
    yi->size[0] = nx;
    yi->size[1] = vstride;
    emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
    m = a->size[0];
    i8 = yi->size[0] * yi->size[1];
    emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
    loop_ub = yi->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = yi->size[0];
      for (i9 = 0; i9 < br; i9++) {
        yi->data[i9 + yi->size[0] * i8] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (x->size[1] == 0)) {
    } else {
      vstride = a->size[0] * (x->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i8 = nx + m;
        for (ic = nx; ic + 1 <= i8; ic++) {
          yi->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i8 = br + k;
        for (ib = br; ib + 1 <= i8; ib++) {
          if (x->data[ib] != 0.0) {
            ia = ar;
            i9 = nx + m;
            for (ic = nx; ic + 1 <= i9; ic++) {
              ia++;
              yi->data[ic] += x->data[ib] * a->data[ia - 1];
            }
          }

          ar += m;
        }

        br += k;
        nx += m;
      }
    }
  }

  emxInit_real_T(&b_yi, 2);
  pp = rt_powd_snf((double)n, 3.0);
  eye((double)p, a);
  i8 = b_yi->size[0] * b_yi->size[1];
  b_yi->size[0] = yi->size[0];
  b_yi->size[1] = yi->size[1];
  emxEnsureCapacity((emxArray__common *)b_yi, i8, (int)sizeof(double));
  loop_ub = yi->size[0] * yi->size[1];
  for (i8 = 0; i8 < loop_ub; i8++) {
    b_yi->data[i8] = yi->data[i8] / (double)n + a->data[i8] / pp;
  }

  emxInit_real_T(&ss, 2);
  emxInit_real_T(&V, 2);
  emxInit_creal_T1(&Vc, 2);
  emxInit_creal_T1(&Dc, 2);
  inv(b_yi, ss);
  eig(ss, Vc, Dc);
  i8 = V->size[0] * V->size[1];
  V->size[0] = Vc->size[0];
  V->size[1] = Vc->size[1];
  emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
  loop_ub = Vc->size[0] * Vc->size[1];
  emxFree_real_T(&b_yi);
  for (i8 = 0; i8 < loop_ub; i8++) {
    V->data[i8] = Vc->data[i8].re;
  }

  emxInit_real_T(&b, 2);
  i8 = b->size[0] * b->size[1];
  b->size[0] = Dc->size[0];
  b->size[1] = Dc->size[1];
  emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
  loop_ub = Dc->size[0] * Dc->size[1];
  for (i8 = 0; i8 < loop_ub; i8++) {
    b->data[i8] = Dc->data[i8].re;
  }

  c_sqrt(b);
  emxInit_real_T(&yj, 2);
  if ((V->size[1] == 1) || (b->size[0] == 1)) {
    i8 = yj->size[0] * yj->size[1];
    yj->size[0] = V->size[0];
    yj->size[1] = b->size[1];
    emxEnsureCapacity((emxArray__common *)yj, i8, (int)sizeof(double));
    loop_ub = V->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = b->size[1];
      for (i9 = 0; i9 < br; i9++) {
        yj->data[i8 + yj->size[0] * i9] = 0.0;
        vstride = V->size[1];
        for (k1 = 0; k1 < vstride; k1++) {
          yj->data[i8 + yj->size[0] * i9] += V->data[i8 + V->size[0] * k1] *
            b->data[k1 + b->size[0] * i9];
        }
      }
    }
  } else {
    k = V->size[1];
    nx = V->size[0];
    vstride = b->size[1];
    i8 = yj->size[0] * yj->size[1];
    yj->size[0] = nx;
    yj->size[1] = vstride;
    emxEnsureCapacity((emxArray__common *)yj, i8, (int)sizeof(double));
    m = V->size[0];
    i8 = yj->size[0] * yj->size[1];
    emxEnsureCapacity((emxArray__common *)yj, i8, (int)sizeof(double));
    loop_ub = yj->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = yj->size[0];
      for (i9 = 0; i9 < br; i9++) {
        yj->data[i9 + yj->size[0] * i8] = 0.0;
      }
    }

    if ((V->size[0] == 0) || (b->size[1] == 0)) {
    } else {
      vstride = V->size[0] * (b->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i8 = nx + m;
        for (ic = nx; ic + 1 <= i8; ic++) {
          yj->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i8 = br + k;
        for (ib = br; ib + 1 <= i8; ib++) {
          if (b->data[ib] != 0.0) {
            ia = ar;
            i9 = nx + m;
            for (ic = nx; ic + 1 <= i9; ic++) {
              ia++;
              yj->data[ic] += b->data[ib] * V->data[ia - 1];
            }
          }

          ar += m;
        }

        br += k;
        nx += m;
      }
    }
  }

  i8 = b->size[0] * b->size[1];
  b->size[0] = V->size[1];
  b->size[1] = V->size[0];
  emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
  loop_ub = V->size[0];
  for (i8 = 0; i8 < loop_ub; i8++) {
    br = V->size[1];
    for (i9 = 0; i9 < br; i9++) {
      b->data[i9 + b->size[0] * i8] = V->data[i8 + V->size[0] * i9];
    }
  }

  if ((yj->size[1] == 1) || (b->size[0] == 1)) {
    i8 = ss->size[0] * ss->size[1];
    ss->size[0] = yj->size[0];
    ss->size[1] = b->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i8, (int)sizeof(double));
    loop_ub = yj->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = b->size[1];
      for (i9 = 0; i9 < br; i9++) {
        ss->data[i8 + ss->size[0] * i9] = 0.0;
        vstride = yj->size[1];
        for (k1 = 0; k1 < vstride; k1++) {
          ss->data[i8 + ss->size[0] * i9] += yj->data[i8 + yj->size[0] * k1] *
            b->data[k1 + b->size[0] * i9];
        }
      }
    }
  } else {
    k = yj->size[1];
    nx = yj->size[0];
    vstride = b->size[1];
    i8 = ss->size[0] * ss->size[1];
    ss->size[0] = nx;
    ss->size[1] = vstride;
    emxEnsureCapacity((emxArray__common *)ss, i8, (int)sizeof(double));
    m = yj->size[0];
    i8 = ss->size[0] * ss->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i8, (int)sizeof(double));
    loop_ub = ss->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = ss->size[0];
      for (i9 = 0; i9 < br; i9++) {
        ss->data[i9 + ss->size[0] * i8] = 0.0;
      }
    }

    if ((yj->size[0] == 0) || (b->size[1] == 0)) {
    } else {
      vstride = yj->size[0] * (b->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i8 = nx + m;
        for (ic = nx; ic + 1 <= i8; ic++) {
          ss->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i8 = br + k;
        for (ib = br; ib + 1 <= i8; ib++) {
          if (b->data[ib] != 0.0) {
            ia = ar;
            i9 = nx + m;
            for (ic = nx; ic + 1 <= i9; ic++) {
              ia++;
              ss->data[ic] += b->data[ib] * yj->data[ia - 1];
            }
          }

          ar += m;
        }

        br += k;
        nx += m;
      }
    }
  }

  i8 = a->size[0] * a->size[1];
  a->size[0] = x->size[0];
  a->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)a, i8, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i8 = 0; i8 < loop_ub; i8++) {
    a->data[i8] = x->data[i8];
  }

  emxInit_real_T(&b_x, 2);
  if ((x->size[1] == 1) || (ss->size[0] == 1)) {
    i8 = b_x->size[0] * b_x->size[1];
    b_x->size[0] = x->size[0];
    b_x->size[1] = ss->size[1];
    emxEnsureCapacity((emxArray__common *)b_x, i8, (int)sizeof(double));
    loop_ub = x->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = ss->size[1];
      for (i9 = 0; i9 < br; i9++) {
        b_x->data[i8 + b_x->size[0] * i9] = 0.0;
        vstride = x->size[1];
        for (k1 = 0; k1 < vstride; k1++) {
          b_x->data[i8 + b_x->size[0] * i9] += x->data[i8 + x->size[0] * k1] *
            ss->data[k1 + ss->size[0] * i9];
        }
      }
    }

    i8 = x->size[0] * x->size[1];
    x->size[0] = b_x->size[0];
    x->size[1] = b_x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
    loop_ub = b_x->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = b_x->size[0];
      for (i9 = 0; i9 < br; i9++) {
        x->data[i9 + x->size[0] * i8] = b_x->data[i9 + b_x->size[0] * i8];
      }
    }
  } else {
    k = x->size[1];
    nx = x->size[0];
    vstride = ss->size[1];
    i8 = x->size[0] * x->size[1];
    x->size[0] = nx;
    x->size[1] = vstride;
    emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
    m = a->size[0];
    i8 = x->size[0] * x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
    loop_ub = x->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = x->size[0];
      for (i9 = 0; i9 < br; i9++) {
        x->data[i9 + x->size[0] * i8] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (ss->size[1] == 0)) {
    } else {
      vstride = a->size[0] * (ss->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i8 = nx + m;
        for (ic = nx; ic + 1 <= i8; ic++) {
          x->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i8 = br + k;
        for (ib = br; ib + 1 <= i8; ib++) {
          if (ss->data[ib] != 0.0) {
            ia = ar;
            i9 = nx + m;
            for (ic = nx; ic + 1 <= i9; ic++) {
              ia++;
              x->data[ic] += ss->data[ib] * a->data[ia - 1];
            }
          }

          ar += m;
        }

        br += k;
        nx += m;
      }
    }
  }

  emxFree_real_T(&b_x);
  emxInit_char_T(&b_method, 2);
  i8 = b_method->size[0] * b_method->size[1];
  b_method->size[0] = method->size[0];
  b_method->size[1] = method->size[1];
  emxEnsureCapacity((emxArray__common *)b_method, i8, (int)sizeof(char));
  loop_ub = method->size[0] * method->size[1];
  for (i8 = 0; i8 < loop_ub; i8++) {
    b_method->data[i8] = method->data[i8];
  }

  upper(b_method, method);
  sort(which_dim);
  emxFree_char_T(&b_method);
  if (which_dim->data[0] != p) {
    emxInit_real_T(&b_p, 2);
    i8 = b_p->size[0] * b_p->size[1];
    b_p->size[0] = 1;
    b_p->size[1] = 1 + which_dim->size[1];
    emxEnsureCapacity((emxArray__common *)b_p, i8, (int)sizeof(double));
    b_p->data[0] = p;
    loop_ub = which_dim->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_p->data[b_p->size[0] * (i8 + 1)] = which_dim->data[which_dim->size[0] *
        i8];
    }

    i8 = which_dim->size[0] * which_dim->size[1];
    which_dim->size[0] = 1;
    which_dim->size[1] = b_p->size[1];
    emxEnsureCapacity((emxArray__common *)which_dim, i8, (int)sizeof(double));
    loop_ub = b_p->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      which_dim->data[which_dim->size[0] * i8] = b_p->data[b_p->size[0] * i8];
    }

    emxFree_real_T(&b_p);
  }

  emxInit_real_T(&ky1, 2);
  emxInit_real_T(&ky2, 2);
  emxInit_real_T(&DD, 2);
  emxInit_real_T(&b_qx, 2);
  emxInit_real_T(&C, 2);
  emxInit_real_T(&U, 2);
  emxInit_real_T1(&s, 1);
  emxInit_real_T1(&dc, 1);
  emxInit_int32_T1(&iidx, 1);
  emxInit_boolean_T1(&c_x, 1);
  emxInit_int32_T1(&ii, 1);
  emxInit_real_T(&b_s, 2);
  emxInit_real_T(&b_Dc, 2);
  emxInit_boolean_T(&c_yi, 2);
  emxInit_real_T(&c_qx, 2);
  emxInit_real_T(&b_ky, 2);
  emxInit_real_T(&b_V, 2);
  emxInit_real_T(&c_V, 2);
  if (b_strcmp(method) || c_strcmp(method)) {
    i8 = ky->size[0] * ky->size[1];
    ky->size[0] = y->size[0];
    ky->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)ky, i8, (int)sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      ky->data[i8] = y->data[i8];
    }

    i8 = ky1->size[0] * ky1->size[1];
    ky1->size[0] = y->size[0];
    ky1->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)ky1, i8, (int)sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      ky1->data[i8] = y->data[i8];
    }

    i8 = ky2->size[0] * ky2->size[1];
    ky2->size[0] = 1;
    ky2->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)ky2, i8, (int)sizeof(double));
    ky2->data[0] = 1.0;
    i8 = DD->size[0] * DD->size[1];
    DD->size[0] = 1;
    DD->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)DD, i8, (int)sizeof(double));
    DD->data[0] = 1.0;
  } else {
    if (99 <= n) {
      nm1d2 = 99;
    } else {
      nm1d2 = n;
    }

    pp = std::floor(rt_powd_snf((double)n, 0.6));
    if ((nm1d2 >= pp) || rtIsNaN(pp)) {
      pp = nm1d2;
    }

    if (rtIsNaN(pp)) {
      i8 = qx->size[0] * qx->size[1];
      qx->size[0] = 1;
      qx->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof(double));
      qx->data[0] = rtNaN;
    } else if (pp < 1.0) {
      i8 = qx->size[0] * qx->size[1];
      qx->size[0] = 1;
      qx->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof(double));
    } else if (rtIsInf(pp) && (1.0 == pp)) {
      i8 = qx->size[0] * qx->size[1];
      qx->size[0] = 1;
      qx->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof(double));
      qx->data[0] = rtNaN;
    } else {
      i8 = qx->size[0] * qx->size[1];
      qx->size[0] = 1;
      qx->size[1] = (int)(pp - 1.0) + 1;
      emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof(double));
      loop_ub = (int)(pp - 1.0);
      for (i8 = 0; i8 <= loop_ub; i8++) {
        qx->data[qx->size[0] * i8] = 1.0 + (double)i8;
      }
    }

    i8 = c_qx->size[0] * c_qx->size[1];
    c_qx->size[0] = 1;
    c_qx->size[1] = qx->size[1];
    emxEnsureCapacity((emxArray__common *)c_qx, i8, (int)sizeof(double));
    loop_ub = qx->size[0] * qx->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      c_qx->data[i8] = qx->data[i8] / (pp + 1.0);
    }

    quantile(y, c_qx, b_qx);
    if ((b_qx->size[0] == 0) || (b_qx->size[1] == 0)) {
      vstride = 0;
    } else {
      nx = b_qx->size[0];
      vstride = b_qx->size[1];
      if (nx >= vstride) {
        vstride = nx;
      }
    }

    nx = b_qx->size[0] * b_qx->size[1];
    i8 = qx->size[0] * qx->size[1];
    qx->size[0] = 1;
    qx->size[1] = vstride;
    emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof(double));
    for (k = 0; k + 1 <= nx; k++) {
      qx->data[k] = b_qx->data[k];
    }

    /* change for C */
    /* qx = unique(qx);%change for C */
    b_repmat(y, (double)qx->size[1], yi);
    repmat(qx, (double)n, yj);
    i8 = c_yi->size[0] * c_yi->size[1];
    c_yi->size[0] = yi->size[0];
    c_yi->size[1] = yi->size[1];
    emxEnsureCapacity((emxArray__common *)c_yi, i8, (int)sizeof(boolean_T));
    loop_ub = yi->size[0] * yi->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      c_yi->data[i8] = (yi->data[i8] - yj->data[i8] < 0.0);
    }

    b_abs(c_yi, ky);
    mean(ky, qx);
    repmat(qx, (double)n, a);
    i8 = ky->size[0] * ky->size[1];
    emxEnsureCapacity((emxArray__common *)ky, i8, (int)sizeof(double));
    nm1d2 = ky->size[0];
    vstride = ky->size[1];
    loop_ub = nm1d2 * vstride;
    for (i8 = 0; i8 < loop_ub; i8++) {
      ky->data[i8] -= a->data[i8];
    }

    b_mean(ky, s);
    c_repmat(s, (double)ky->size[1], ky1);
    i8 = ky1->size[0] * ky1->size[1];
    ky1->size[0] = ky->size[0];
    ky1->size[1] = ky->size[1];
    emxEnsureCapacity((emxArray__common *)ky1, i8, (int)sizeof(double));
    loop_ub = ky->size[0] * ky->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      ky1->data[i8] = ky->data[i8] - ky1->data[i8];
    }

    mean(ky1, qx);
    repmat(qx, (double)n, a);
    i8 = ky1->size[0] * ky1->size[1];
    emxEnsureCapacity((emxArray__common *)ky1, i8, (int)sizeof(double));
    nm1d2 = ky1->size[0];
    vstride = ky1->size[1];
    loop_ub = nm1d2 * vstride;
    for (i8 = 0; i8 < loop_ub; i8++) {
      ky1->data[i8] -= a->data[i8];
    }

    if (!((ky->size[0] == 0) || (ky->size[1] == 0))) {
      nx = ky->size[0];
    } else if (!((ky1->size[0] == 0) || (ky1->size[1] == 0))) {
      nx = ky1->size[0];
    } else {
      nx = ky->size[0];
      if (!(nx >= 0)) {
        nx = 0;
      }

      if (ky1->size[0] > nx) {
        nx = ky1->size[0];
      }
    }

    empty_non_axis_sizes = (nx == 0);
    if (empty_non_axis_sizes || (!((ky->size[0] == 0) || (ky->size[1] == 0)))) {
      nm1d2 = ky->size[1];
    } else {
      nm1d2 = 0;
    }

    if (empty_non_axis_sizes || (!((ky1->size[0] == 0) || (ky1->size[1] == 0))))
    {
      vstride = ky1->size[1];
    } else {
      vstride = 0;
    }

    i8 = b_ky->size[0] * b_ky->size[1];
    b_ky->size[0] = nx;
    b_ky->size[1] = nm1d2 + vstride;
    emxEnsureCapacity((emxArray__common *)b_ky, i8, (int)sizeof(double));
    for (i8 = 0; i8 < nm1d2; i8++) {
      for (i9 = 0; i9 < nx; i9++) {
        b_ky->data[i9 + b_ky->size[0] * i8] = ky->data[i9 + nx * i8];
      }
    }

    for (i8 = 0; i8 < vstride; i8++) {
      for (i9 = 0; i9 < nx; i9++) {
        b_ky->data[i9 + b_ky->size[0] * (i8 + nm1d2)] = ky1->data[i9 + nx * i8];
      }
    }

    i8 = ky->size[0] * ky->size[1];
    ky->size[0] = b_ky->size[0];
    ky->size[1] = b_ky->size[1];
    emxEnsureCapacity((emxArray__common *)ky, i8, (int)sizeof(double));
    loop_ub = b_ky->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = b_ky->size[0];
      for (i9 = 0; i9 < br; i9++) {
        ky->data[i9 + ky->size[0] * i8] = b_ky->data[i9 + b_ky->size[0] * i8];
      }
    }

    /* clear yi ky1; */
    i8 = a->size[0] * a->size[1];
    a->size[0] = ky->size[1];
    a->size[1] = ky->size[0];
    emxEnsureCapacity((emxArray__common *)a, i8, (int)sizeof(double));
    loop_ub = ky->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = ky->size[1];
      for (i9 = 0; i9 < br; i9++) {
        a->data[i9 + a->size[0] * i8] = ky->data[i8 + ky->size[0] * i9];
      }
    }

    if ((a->size[1] == 1) || (ky->size[0] == 1)) {
      i8 = C->size[0] * C->size[1];
      C->size[0] = a->size[0];
      C->size[1] = ky->size[1];
      emxEnsureCapacity((emxArray__common *)C, i8, (int)sizeof(double));
      loop_ub = a->size[0];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = ky->size[1];
        for (i9 = 0; i9 < br; i9++) {
          C->data[i8 + C->size[0] * i9] = 0.0;
          vstride = a->size[1];
          for (k1 = 0; k1 < vstride; k1++) {
            C->data[i8 + C->size[0] * i9] += a->data[i8 + a->size[0] * k1] *
              ky->data[k1 + ky->size[0] * i9];
          }
        }
      }
    } else {
      k = a->size[1];
      nx = a->size[0];
      vstride = ky->size[1];
      i8 = C->size[0] * C->size[1];
      C->size[0] = nx;
      C->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)C, i8, (int)sizeof(double));
      m = a->size[0];
      i8 = C->size[0] * C->size[1];
      emxEnsureCapacity((emxArray__common *)C, i8, (int)sizeof(double));
      loop_ub = C->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = C->size[0];
        for (i9 = 0; i9 < br; i9++) {
          C->data[i9 + C->size[0] * i8] = 0.0;
        }
      }

      if ((a->size[0] == 0) || (ky->size[1] == 0)) {
      } else {
        vstride = a->size[0] * (ky->size[1] - 1);
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          i8 = nx + m;
          for (ic = nx; ic + 1 <= i8; ic++) {
            C->data[ic] = 0.0;
          }

          nx += m;
        }

        br = 0;
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          ar = 0;
          i8 = br + k;
          for (ib = br; ib + 1 <= i8; ib++) {
            if (ky->data[ib] != 0.0) {
              ia = ar;
              i9 = nx + m;
              for (ic = nx; ic + 1 <= i9; ic++) {
                ia++;
                C->data[ic] += ky->data[ib] * a->data[ia - 1];
              }
            }

            ar += m;
          }

          br += k;
          nx += m;
        }
      }
    }

    eig(C, Vc, Dc);

    /* change for C */
    i8 = V->size[0] * V->size[1];
    V->size[0] = Vc->size[0];
    V->size[1] = Vc->size[1];
    emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
    loop_ub = Vc->size[0] * Vc->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      V->data[i8] = Vc->data[i8].re;
    }

    /* change for C */
    /* change for C */
    /* clear C; */
    i8 = b_Dc->size[0] * b_Dc->size[1];
    b_Dc->size[0] = Dc->size[0];
    b_Dc->size[1] = Dc->size[1];
    emxEnsureCapacity((emxArray__common *)b_Dc, i8, (int)sizeof(double));
    loop_ub = Dc->size[0] * Dc->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_Dc->data[i8] = Dc->data[i8].re;
    }

    diag(b_Dc, s);
    c_abs(s, dc);
    c_sort(dc, iidx);
    nm1d2 = 2;
    if (dc->size[0] != 1) {
      nm1d2 = 1;
    }

    i8 = s->size[0];
    s->size[0] = dc->size[0];
    emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
    loop_ub = dc->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      s->data[i8] = dc->data[i8];
    }

    if (nm1d2 <= 1) {
      i8 = dc->size[0];
    } else {
      i8 = 1;
    }

    if ((!(dc->size[0] == 0)) && (i8 > 1)) {
      vstride = 1;
      k = 1;
      while (k <= nm1d2 - 1) {
        vstride *= dc->size[0];
        k = 2;
      }

      for (j = 0; j + 1 <= vstride; j++) {
        for (k = 1; k < i8; k++) {
          s->data[j + k * vstride] += s->data[j + (k - 1) * vstride];
        }
      }
    }

    pp = sum(dc);
    i8 = c_x->size[0];
    c_x->size[0] = s->size[0];
    emxEnsureCapacity((emxArray__common *)c_x, i8, (int)sizeof(boolean_T));
    loop_ub = s->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      c_x->data[i8] = (s->data[i8] / pp >= 0.99);
    }

    nx = c_x->size[0];
    vstride = 0;
    i8 = ii->size[0];
    ii->size[0] = c_x->size[0];
    emxEnsureCapacity((emxArray__common *)ii, i8, (int)sizeof(int));
    nm1d2 = 1;
    exitg2 = false;
    while ((!exitg2) && (nm1d2 <= nx)) {
      b_guard1 = false;
      if (c_x->data[nm1d2 - 1]) {
        vstride++;
        ii->data[vstride - 1] = nm1d2;
        if (vstride >= nx) {
          exitg2 = true;
        } else {
          b_guard1 = true;
        }
      } else {
        b_guard1 = true;
      }

      if (b_guard1) {
        nm1d2++;
      }
    }

    if (c_x->size[0] == 1) {
      if (vstride == 0) {
        i8 = ii->size[0];
        ii->size[0] = 0;
        emxEnsureCapacity((emxArray__common *)ii, i8, (int)sizeof(int));
      }
    } else {
      i8 = ii->size[0];
      if (1 > vstride) {
        ii->size[0] = 0;
      } else {
        ii->size[0] = vstride;
      }

      emxEnsureCapacity((emxArray__common *)ii, i8, (int)sizeof(int));
    }

    i8 = s->size[0];
    s->size[0] = ii->size[0];
    emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
    loop_ub = ii->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      s->data[i8] = ii->data[i8];
    }

    nx = s->size[0];
    vstride = (int)s->data[0];
    if (s->size[0] > 1) {
      for (nm1d2 = 1; nm1d2 + 1 <= nx; nm1d2++) {
        if ((int)s->data[nm1d2] < vstride) {
          vstride = (int)s->data[nm1d2];
        }
      }
    }

    nm1d2 = V->size[0];
    i8 = b_V->size[0] * b_V->size[1];
    b_V->size[0] = nm1d2;
    b_V->size[1] = iidx->size[0];
    emxEnsureCapacity((emxArray__common *)b_V, i8, (int)sizeof(double));
    loop_ub = iidx->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      for (i9 = 0; i9 < nm1d2; i9++) {
        b_V->data[i9 + b_V->size[0] * i8] = V->data[i9 + V->size[0] *
          (iidx->data[i8] - 1)];
      }
    }

    i8 = V->size[0] * V->size[1];
    V->size[0] = b_V->size[0];
    V->size[1] = b_V->size[1];
    emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
    loop_ub = b_V->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = b_V->size[0];
      for (i9 = 0; i9 < br; i9++) {
        V->data[i9 + V->size[0] * i8] = b_V->data[i9 + b_V->size[0] * i8];
      }
    }

    if (1 > vstride) {
      loop_ub = 0;
      br = 0;
    } else {
      loop_ub = vstride;
      br = vstride;
    }

    nm1d2 = V->size[0];
    i8 = c_V->size[0] * c_V->size[1];
    c_V->size[0] = nm1d2;
    c_V->size[1] = br;
    emxEnsureCapacity((emxArray__common *)c_V, i8, (int)sizeof(double));
    for (i8 = 0; i8 < br; i8++) {
      for (i9 = 0; i9 < nm1d2; i9++) {
        c_V->data[i9 + c_V->size[0] * i8] = V->data[i9 + V->size[0] * i8];
      }
    }

    i8 = V->size[0] * V->size[1];
    V->size[0] = c_V->size[0];
    V->size[1] = c_V->size[1];
    emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
    br = c_V->size[1];
    for (i8 = 0; i8 < br; i8++) {
      vstride = c_V->size[0];
      for (i9 = 0; i9 < vstride; i9++) {
        V->data[i9 + V->size[0] * i8] = c_V->data[i9 + c_V->size[0] * i8];
      }
    }

    if ((ky->size[1] == 1) || (V->size[0] == 1)) {
      i8 = U->size[0] * U->size[1];
      U->size[0] = ky->size[0];
      U->size[1] = V->size[1];
      emxEnsureCapacity((emxArray__common *)U, i8, (int)sizeof(double));
      br = ky->size[0];
      for (i8 = 0; i8 < br; i8++) {
        vstride = V->size[1];
        for (i9 = 0; i9 < vstride; i9++) {
          U->data[i8 + U->size[0] * i9] = 0.0;
          nm1d2 = ky->size[1];
          for (k1 = 0; k1 < nm1d2; k1++) {
            U->data[i8 + U->size[0] * i9] += ky->data[i8 + ky->size[0] * k1] *
              V->data[k1 + V->size[0] * i9];
          }
        }
      }
    } else {
      k = ky->size[1];
      nx = ky->size[0];
      vstride = V->size[1];
      i8 = U->size[0] * U->size[1];
      U->size[0] = nx;
      U->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)U, i8, (int)sizeof(double));
      m = ky->size[0];
      i8 = U->size[0] * U->size[1];
      emxEnsureCapacity((emxArray__common *)U, i8, (int)sizeof(double));
      br = U->size[1];
      for (i8 = 0; i8 < br; i8++) {
        vstride = U->size[0];
        for (i9 = 0; i9 < vstride; i9++) {
          U->data[i9 + U->size[0] * i8] = 0.0;
        }
      }

      if ((ky->size[0] == 0) || (V->size[1] == 0)) {
      } else {
        vstride = ky->size[0] * (V->size[1] - 1);
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          i8 = nx + m;
          for (ic = nx; ic + 1 <= i8; ic++) {
            U->data[ic] = 0.0;
          }

          nx += m;
        }

        br = 0;
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          ar = 0;
          i8 = br + k;
          for (ib = br; ib + 1 <= i8; ib++) {
            if (V->data[ib] != 0.0) {
              ia = ar;
              i9 = nx + m;
              for (ic = nx; ic + 1 <= i9; ic++) {
                ia++;
                U->data[ic] += V->data[ib] * ky->data[ia - 1];
              }
            }

            ar += m;
          }

          br += k;
          nx += m;
        }
      }
    }

    i8 = s->size[0];
    s->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
    for (i8 = 0; i8 < loop_ub; i8++) {
      s->data[i8] = dc->data[i8];
    }

    d_sqrt(s);

    /* ky1 = bsxfun(@(x,c)x./c, U, s'); */
    i8 = b_s->size[0] * b_s->size[1];
    b_s->size[0] = 1;
    b_s->size[1] = s->size[0];
    emxEnsureCapacity((emxArray__common *)b_s, i8, (int)sizeof(double));
    loop_ub = s->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_s->data[b_s->size[0] * i8] = s->data[i8];
    }

    repmat(b_s, (double)U->size[0], a);
    rdivide(U, a, ky1);
    b_diag(s, a);
    i8 = b->size[0] * b->size[1];
    b->size[0] = V->size[1];
    b->size[1] = V->size[0];
    emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
    loop_ub = V->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = V->size[1];
      for (i9 = 0; i9 < br; i9++) {
        b->data[i9 + b->size[0] * i8] = V->data[i8 + V->size[0] * i9];
      }
    }

    if ((a->size[1] == 1) || (b->size[0] == 1)) {
      i8 = ky2->size[0] * ky2->size[1];
      ky2->size[0] = a->size[0];
      ky2->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)ky2, i8, (int)sizeof(double));
      loop_ub = a->size[0];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = b->size[1];
        for (i9 = 0; i9 < br; i9++) {
          ky2->data[i8 + ky2->size[0] * i9] = 0.0;
          vstride = a->size[1];
          for (k1 = 0; k1 < vstride; k1++) {
            ky2->data[i8 + ky2->size[0] * i9] += a->data[i8 + a->size[0] * k1] *
              b->data[k1 + b->size[0] * i9];
          }
        }
      }
    } else {
      k = a->size[1];
      nx = a->size[0];
      vstride = b->size[1];
      i8 = ky2->size[0] * ky2->size[1];
      ky2->size[0] = nx;
      ky2->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)ky2, i8, (int)sizeof(double));
      m = a->size[0];
      i8 = ky2->size[0] * ky2->size[1];
      emxEnsureCapacity((emxArray__common *)ky2, i8, (int)sizeof(double));
      loop_ub = ky2->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = ky2->size[0];
        for (i9 = 0; i9 < br; i9++) {
          ky2->data[i9 + ky2->size[0] * i8] = 0.0;
        }
      }

      if ((a->size[0] == 0) || (b->size[1] == 0)) {
      } else {
        vstride = a->size[0] * (b->size[1] - 1);
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          i8 = nx + m;
          for (ic = nx; ic + 1 <= i8; ic++) {
            ky2->data[ic] = 0.0;
          }

          nx += m;
        }

        br = 0;
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          ar = 0;
          i8 = br + k;
          for (ib = br; ib + 1 <= i8; ib++) {
            if (b->data[ib] != 0.0) {
              ia = ar;
              i9 = nx + m;
              for (ic = nx; ic + 1 <= i9; ic++) {
                ia++;
                ky2->data[ic] += b->data[ib] * a->data[ia - 1];
              }
            }

            ar += m;
          }

          br += k;
          nx += m;
        }
      }
    }

    if ((ky1->size[1] == 1) || (ky2->size[0] == 1)) {
      i8 = ky->size[0] * ky->size[1];
      ky->size[0] = ky1->size[0];
      ky->size[1] = ky2->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i8, (int)sizeof(double));
      loop_ub = ky1->size[0];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = ky2->size[1];
        for (i9 = 0; i9 < br; i9++) {
          ky->data[i8 + ky->size[0] * i9] = 0.0;
          vstride = ky1->size[1];
          for (k1 = 0; k1 < vstride; k1++) {
            ky->data[i8 + ky->size[0] * i9] += ky1->data[i8 + ky1->size[0] * k1]
              * ky2->data[k1 + ky2->size[0] * i9];
          }
        }
      }
    } else {
      k = ky1->size[1];
      nx = ky1->size[0];
      vstride = ky2->size[1];
      i8 = ky->size[0] * ky->size[1];
      ky->size[0] = nx;
      ky->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)ky, i8, (int)sizeof(double));
      m = ky1->size[0];
      i8 = ky->size[0] * ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i8, (int)sizeof(double));
      loop_ub = ky->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = ky->size[0];
        for (i9 = 0; i9 < br; i9++) {
          ky->data[i9 + ky->size[0] * i8] = 0.0;
        }
      }

      if ((ky1->size[0] == 0) || (ky2->size[1] == 0)) {
      } else {
        vstride = ky1->size[0] * (ky2->size[1] - 1);
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          i8 = nx + m;
          for (ic = nx; ic + 1 <= i8; ic++) {
            ky->data[ic] = 0.0;
          }

          nx += m;
        }

        br = 0;
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          ar = 0;
          i8 = br + k;
          for (ib = br; ib + 1 <= i8; ib++) {
            if (ky2->data[ib] != 0.0) {
              ia = ar;
              i9 = nx + m;
              for (ic = nx; ic + 1 <= i9; ic++) {
                ia++;
                ky->data[ic] += ky2->data[ib] * ky1->data[ia - 1];
              }
            }

            ar += m;
          }

          br += k;
          nx += m;
        }
      }
    }

    /* clear U V; */
    i8 = b->size[0] * b->size[1];
    b->size[0] = ky2->size[1];
    b->size[1] = ky2->size[0];
    emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
    loop_ub = ky2->size[0];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = ky2->size[1];
      for (i9 = 0; i9 < br; i9++) {
        b->data[i9 + b->size[0] * i8] = ky2->data[i8 + ky2->size[0] * i9];
      }
    }

    if ((ky2->size[1] == 1) || (b->size[0] == 1)) {
      i8 = DD->size[0] * DD->size[1];
      DD->size[0] = ky2->size[0];
      DD->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)DD, i8, (int)sizeof(double));
      loop_ub = ky2->size[0];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = b->size[1];
        for (i9 = 0; i9 < br; i9++) {
          DD->data[i8 + DD->size[0] * i9] = 0.0;
          vstride = ky2->size[1];
          for (k1 = 0; k1 < vstride; k1++) {
            DD->data[i8 + DD->size[0] * i9] += ky2->data[i8 + ky2->size[0] * k1]
              * b->data[k1 + b->size[0] * i9];
          }
        }
      }
    } else {
      k = ky2->size[1];
      nx = ky2->size[0];
      vstride = b->size[1];
      i8 = DD->size[0] * DD->size[1];
      DD->size[0] = nx;
      DD->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)DD, i8, (int)sizeof(double));
      m = ky2->size[0];
      i8 = DD->size[0] * DD->size[1];
      emxEnsureCapacity((emxArray__common *)DD, i8, (int)sizeof(double));
      loop_ub = DD->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = DD->size[0];
        for (i9 = 0; i9 < br; i9++) {
          DD->data[i9 + DD->size[0] * i8] = 0.0;
        }
      }

      if ((ky2->size[0] == 0) || (b->size[1] == 0)) {
      } else {
        vstride = ky2->size[0] * (b->size[1] - 1);
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          i8 = nx + m;
          for (ic = nx; ic + 1 <= i8; ic++) {
            DD->data[ic] = 0.0;
          }

          nx += m;
        }

        br = 0;
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          ar = 0;
          i8 = br + k;
          for (ib = br; ib + 1 <= i8; ib++) {
            if (b->data[ib] != 0.0) {
              ia = ar;
              i9 = nx + m;
              for (ic = nx; ic + 1 <= i9; ic++) {
                ia++;
                DD->data[ic] += b->data[ib] * ky2->data[ia - 1];
              }
            }

            ar += m;
          }

          br += k;
          nx += m;
        }
      }
    }
  }

  emxFree_real_T(&c_V);
  emxFree_real_T(&b_V);
  emxFree_real_T(&b_ky);
  emxFree_real_T(&c_qx);
  emxFree_boolean_T(&c_yi);
  emxFree_real_T(&b_Dc);
  emxFree_real_T(&b_s);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&c_x);
  emxFree_real_T(&b_qx);
  if (p == 0) {
    i8 = BB1D->size[0];
    BB1D->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)BB1D, i8, (int)sizeof(double));
    BB1D->data[0] = 1.0;
  } else {
    emxInit_real_T2(&BB, 3);
    i8 = BB->size[0] * BB->size[1] * BB->size[2];
    BB->size[0] = p;
    BB->size[1] = p;
    BB->size[2] = p;
    emxEnsureCapacity((emxArray__common *)BB, i8, (int)sizeof(double));
    loop_ub = p * p * p;
    for (i8 = 0; i8 < loop_ub; i8++) {
      BB->data[i8] = 0.0;
    }

    emxInit_real_T(&r0, 2);
    eye((double)p, r0);
    loop_ub = r0->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      br = r0->size[0];
      for (i9 = 0; i9 < br; i9++) {
        BB->data[(i9 + BB->size[0] * i8) + BB->size[0] * BB->size[1] * (p - 1)] =
          r0->data[i9 + r0->size[0] * i8];
      }
    }

    emxFree_real_T(&r0);
    emxInit_real_T(&onexi, 2);
    i8 = onexi->size[0] * onexi->size[1];
    onexi->size[0] = n;
    onexi->size[1] = p + 1;
    emxEnsureCapacity((emxArray__common *)onexi, i8, (int)sizeof(double));
    loop_ub = n * (p + 1);
    for (i8 = 0; i8 < loop_ub; i8++) {
      onexi->data[i8] = 1.0;
    }

    emxInit_real_T(&B, 2);
    eye((double)p, B);
    i8 = yi->size[0] * yi->size[1];
    yi->size[0] = p;
    yi->size[1] = p;
    emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
    loop_ub = p * p;
    for (i8 = 0; i8 < loop_ub; i8++) {
      yi->data[i8] = 0.0;
    }

    if (p > 1) {
      iter = 0;
      emxInit_real_T(&D, 2);
      emxInit_real_T(&Ifast, 2);
      emxInit_real_T(&xfast, 2);
      emxInit_real_T(&xij, 2);
      emxInit_real_T1(&dxij, 1);
      emxInit_real_T(&abi, 2);
      emxInit_real_T(&dd, 2);
      emxInit_real_T(&kxijy, 2);
      emxInit_real_T(&ddx, 2);
      emxInit_real_T(&tmp, 2);
      emxInit_real_T(&B0, 2);
      emxInit_real_T1(&d_B, 1);
      emxInit_creal_T(&R, 1);
      emxInit_real_T(&b_C, 2);
      emxInit_real_T(&c_C, 2);
      emxInit_real_T(&d_C, 2);
      emxInit_real_T1(&e_C, 1);
      emxInit_int32_T(&r1, 2);
      emxInit_real_T(&b_y, 2);
      emxInit_real_T(&c_y, 2);
      emxInit_real_T(&d_y, 2);
      emxInitMatrix_cell_wrap_0(reshapes);
      emxInit_real_T(&e_y, 2);
      emxInit_real_T(&f_y, 2);
      emxInit_real_T(&g_y, 2);
      emxInit_real_T(&b_dd, 2);
      emxInit_real_T(&b_abi, 2);
      emxInit_real_T(&f_C, 2);
      emxInit_real_T(&d_V, 2);
      emxInit_real_T(&b_U, 2);
      emxInit_real_T(&b_xfast, 2);
      emxInit_real_T(&g_C, 2);
      emxInit_real_T1(&e_V, 1);
      emxInit_real_T(&e_B, 2);
      emxInit_int32_T(&b_Ifast, 2);
      emxInit_int32_T(&c_Ifast, 2);
      emxInit_real_T1(&b_dc, 1);
      emxInit_real_T(&c_abi, 2);
      emxInit_int32_T(&d_Ifast, 2);
      emxInit_int32_T(&e_Ifast, 2);
      exitg1 = false;
      while ((!exitg1) && (iter <= which_dim->size[1] - 1)) {
        ip = which_dim->data[iter];
        nm1d2 = (int)std::floor((double)n / 2.0);
        if (which_dim->data[iter] <= nm1d2) {
          b_m = which_dim->data[iter];
        } else {
          b_m = nm1d2;
        }

        K = ((which_dim->data[iter] <= 10.0) * 5 + ((which_dim->data[iter] <=
               20.0) * (which_dim->data[iter] > 10.0) << 1)) + (which_dim->
          data[iter] > 20.0);
        if (d_strcmp(method)) {
          K = 1;
        }

        if (1.0 > which_dim->data[iter]) {
          loop_ub = 0;
        } else {
          loop_ub = (int)which_dim->data[iter];
        }

        vstride = B->size[0];
        i8 = e_B->size[0] * e_B->size[1];
        e_B->size[0] = vstride;
        e_B->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)e_B, i8, (int)sizeof(double));
        for (i8 = 0; i8 < loop_ub; i8++) {
          for (i9 = 0; i9 < vstride; i9++) {
            e_B->data[i9 + e_B->size[0] * i8] = B->data[i9 + B->size[0] * i8];
          }
        }

        i8 = B->size[0] * B->size[1];
        B->size[0] = e_B->size[0];
        B->size[1] = e_B->size[1];
        emxEnsureCapacity((emxArray__common *)B, i8, (int)sizeof(double));
        loop_ub = e_B->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          br = e_B->size[0];
          for (i9 = 0; i9 < br; i9++) {
            B->data[i9 + B->size[0] * i8] = e_B->data[i9 + e_B->size[0] * i8];
          }
        }

        for (b_iter = 0; b_iter < K; b_iter++) {
          if (d_strcmp(method)) {
            i8 = V->size[0] * V->size[1];
            V->size[0] = y->size[0];
            V->size[1] = y->size[1];
            emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
            loop_ub = y->size[0] * y->size[1];
            for (i8 = 0; i8 < loop_ub; i8++) {
              V->data[i8] = y->data[i8];
            }

            b_m = 1.0;
          } else if ((x->size[1] == 1) || (B->size[0] == 1)) {
            i8 = V->size[0] * V->size[1];
            V->size[0] = x->size[0];
            V->size[1] = B->size[1];
            emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
            loop_ub = x->size[0];
            for (i8 = 0; i8 < loop_ub; i8++) {
              br = B->size[1];
              for (i9 = 0; i9 < br; i9++) {
                V->data[i8 + V->size[0] * i9] = 0.0;
                vstride = x->size[1];
                for (k1 = 0; k1 < vstride; k1++) {
                  V->data[i8 + V->size[0] * i9] += x->data[i8 + x->size[0] * k1]
                    * B->data[k1 + B->size[0] * i9];
                }
              }
            }
          } else {
            k = x->size[1];
            nx = x->size[0];
            vstride = B->size[1];
            i8 = V->size[0] * V->size[1];
            V->size[0] = nx;
            V->size[1] = vstride;
            emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
            m = x->size[0];
            i8 = V->size[0] * V->size[1];
            emxEnsureCapacity((emxArray__common *)V, i8, (int)sizeof(double));
            loop_ub = V->size[1];
            for (i8 = 0; i8 < loop_ub; i8++) {
              br = V->size[0];
              for (i9 = 0; i9 < br; i9++) {
                V->data[i9 + V->size[0] * i8] = 0.0;
              }
            }

            if ((x->size[0] == 0) || (B->size[1] == 0)) {
            } else {
              vstride = x->size[0] * (B->size[1] - 1);
              nx = 0;
              while ((m > 0) && (nx <= vstride)) {
                i8 = nx + m;
                for (ic = nx; ic + 1 <= i8; ic++) {
                  V->data[ic] = 0.0;
                }

                nx += m;
              }

              br = 0;
              nx = 0;
              while ((m > 0) && (nx <= vstride)) {
                ar = 0;
                i8 = br + k;
                for (ib = br; ib + 1 <= i8; ib++) {
                  if (B->data[ib] != 0.0) {
                    ia = ar;
                    i9 = nx + m;
                    for (ic = nx; ic + 1 <= i9; ic++) {
                      ia++;
                      V->data[ic] += B->data[ib] * x->data[ia - 1];
                    }
                  }

                  ar += m;
                }

                br += k;
                nx += m;
              }
            }
          }

          if (100 <= n) {
            nm1d2 = 100;
          } else {
            nm1d2 = n;
          }

          pp = std::floor(std::sqrt((double)n));
          if ((nm1d2 >= pp) || rtIsNaN(pp)) {
            b_nm1d2 = nm1d2;
          } else {
            b_nm1d2 = pp;
          }

          unifD(V, b_nm1d2, Ifast);
          i8 = b_Ifast->size[0] * b_Ifast->size[1];
          b_Ifast->size[0] = Ifast->size[0];
          b_Ifast->size[1] = Ifast->size[1];
          emxEnsureCapacity((emxArray__common *)b_Ifast, i8, (int)sizeof(int));
          loop_ub = Ifast->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            br = Ifast->size[0];
            for (i9 = 0; i9 < br; i9++) {
              b_Ifast->data[i9 + b_Ifast->size[0] * i8] = (int)Ifast->data[i9 +
                Ifast->size[0] * i8];
            }
          }

          Ifast_idx_0 = Ifast->size[0] * Ifast->size[1];
          loop_ub = V->size[1];
          i8 = U->size[0] * U->size[1];
          U->size[0] = Ifast_idx_0;
          U->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)U, i8, (int)sizeof(double));
          for (i8 = 0; i8 < loop_ub; i8++) {
            for (i9 = 0; i9 < Ifast_idx_0; i9++) {
              U->data[i9 + U->size[0] * i8] = V->data[(b_Ifast->data[i9] +
                V->size[0] * i8) - 1];
            }
          }

          i8 = c_Ifast->size[0] * c_Ifast->size[1];
          c_Ifast->size[0] = Ifast->size[0];
          c_Ifast->size[1] = Ifast->size[1];
          emxEnsureCapacity((emxArray__common *)c_Ifast, i8, (int)sizeof(int));
          loop_ub = Ifast->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            br = Ifast->size[0];
            for (i9 = 0; i9 < br; i9++) {
              c_Ifast->data[i9 + c_Ifast->size[0] * i8] = (int)Ifast->data[i9 +
                Ifast->size[0] * i8];
            }
          }

          Ifast_idx_0 = Ifast->size[0] * Ifast->size[1];
          loop_ub = x->size[1];
          i8 = xfast->size[0] * xfast->size[1];
          xfast->size[0] = Ifast_idx_0;
          xfast->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)xfast, i8, (int)sizeof(double));
          for (i8 = 0; i8 < loop_ub; i8++) {
            for (i9 = 0; i9 < Ifast_idx_0; i9++) {
              xfast->data[i9 + xfast->size[0] * i8] = x->data[(c_Ifast->data[i9]
                + x->size[0] * i8) - 1];
            }
          }

          Ifast_idx_0 = Ifast->size[0] * Ifast->size[1];
          b_std(V, qx);
          pp = 1.2 * c_mean(qx);
          pp /= rt_powd_snf((double)n, 1.0 / (b_m + 4.0));
          h2 = 2.0 * pp * pp;
          if (e_strcmp(method) || d_strcmp(method) || c_strcmp(method) || (ip >
               5.0)) {
            i8 = yi->size[0] * yi->size[1];
            yi->size[0] = p;
            yi->size[1] = p;
            emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
            loop_ub = p * p;
            for (i8 = 0; i8 < loop_ub; i8++) {
              yi->data[i8] = 0.0;
            }

            for (j = 0; j < Ifast_idx_0; j++) {
              i8 = xij->size[0] * xij->size[1];
              xij->size[0] = n;
              xij->size[1] = p;
              emxEnsureCapacity((emxArray__common *)xij, i8, (int)sizeof(double));
              for (nx = 0; nx < p; nx++) {
                loop_ub = x->size[0] - 1;
                i8 = e_Ifast->size[0] * e_Ifast->size[1];
                e_Ifast->size[0] = Ifast->size[0];
                e_Ifast->size[1] = Ifast->size[1];
                emxEnsureCapacity((emxArray__common *)e_Ifast, i8, (int)sizeof
                                  (int));
                br = Ifast->size[1];
                for (i8 = 0; i8 < br; i8++) {
                  vstride = Ifast->size[0];
                  for (i9 = 0; i9 < vstride; i9++) {
                    e_Ifast->data[i9 + e_Ifast->size[0] * i8] = (int)Ifast->
                      data[i9 + Ifast->size[0] * i8];
                  }
                }

                pp = x->data[(e_Ifast->data[j] + x->size[0] * nx) - 1];
                for (i8 = 0; i8 <= loop_ub; i8++) {
                  xij->data[i8 + xij->size[0] * nx] = x->data[i8 + x->size[0] *
                    nx] - pp;
                }
              }

              i8 = dxij->size[0];
              dxij->size[0] = n;
              emxEnsureCapacity((emxArray__common *)dxij, i8, (int)sizeof(double));
              for (i8 = 0; i8 < n; i8++) {
                dxij->data[i8] = 0.0;
              }

              for (nx = 0; nx < (int)b_m; nx++) {
                loop_ub = V->size[0];
                i8 = d_Ifast->size[0] * d_Ifast->size[1];
                d_Ifast->size[0] = Ifast->size[0];
                d_Ifast->size[1] = Ifast->size[1];
                emxEnsureCapacity((emxArray__common *)d_Ifast, i8, (int)sizeof
                                  (int));
                br = Ifast->size[1];
                for (i8 = 0; i8 < br; i8++) {
                  vstride = Ifast->size[0];
                  for (i9 = 0; i9 < vstride; i9++) {
                    d_Ifast->data[i9 + d_Ifast->size[0] * i8] = (int)Ifast->
                      data[i9 + Ifast->size[0] * i8];
                  }
                }

                pp = V->data[(d_Ifast->data[j] + V->size[0] * nx) - 1];
                i8 = e_V->size[0];
                e_V->size[0] = loop_ub;
                emxEnsureCapacity((emxArray__common *)e_V, i8, (int)sizeof
                                  (double));
                for (i8 = 0; i8 < loop_ub; i8++) {
                  e_V->data[i8] = V->data[i8 + V->size[0] * nx] - pp;
                }

                b_power(e_V, s);
                i8 = dxij->size[0];
                emxEnsureCapacity((emxArray__common *)dxij, i8, (int)sizeof
                                  (double));
                loop_ub = dxij->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  dxij->data[i8] += s->data[i8];
                }
              }

              i8 = s->size[0];
              s->size[0] = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
              loop_ub = dxij->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                s->data[i8] = dxij->data[i8];
              }

              f_sort(s);
              if ((h2 >= s->data[(int)(2.0 * b_m) - 1]) || rtIsNaN(s->data[(int)
                   (2.0 * b_m) - 1])) {
                pp = h2;
              } else {
                pp = s->data[(int)(2.0 * b_m) - 1];
              }

              i8 = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)dxij, i8, (int)sizeof(double));
              loop_ub = dxij->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                dxij->data[i8] = -dxij->data[i8] / pp;
              }

              b_exp(dxij);
              loop_ub = xij->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = xij->size[0];
                for (i9 = 0; i9 < br; i9++) {
                  onexi->data[i9 + onexi->size[0] * i8] = xij->data[i9 +
                    xij->size[0] * i8];
                }
              }

              c_repmat(dxij, (double)p + 1.0, a);
              i8 = C->size[0] * C->size[1];
              C->size[0] = onexi->size[1];
              C->size[1] = onexi->size[0];
              emxEnsureCapacity((emxArray__common *)C, i8, (int)sizeof(double));
              loop_ub = onexi->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = onexi->size[1];
                for (i9 = 0; i9 < br; i9++) {
                  C->data[i9 + C->size[0] * i8] = onexi->data[i8 + onexi->size[0]
                    * i9] * a->data[i8 + a->size[0] * i9];
                }
              }

              /* abi = inv(xk*onexi+eye(p+1)/n)*(xk*ky1); */
              if ((C->size[1] == 1) || (onexi->size[0] == 1)) {
                i8 = b_C->size[0] * b_C->size[1];
                b_C->size[0] = C->size[0];
                b_C->size[1] = onexi->size[1];
                emxEnsureCapacity((emxArray__common *)b_C, i8, (int)sizeof
                                  (double));
                loop_ub = C->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = onexi->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    b_C->data[i8 + b_C->size[0] * i9] = 0.0;
                    vstride = C->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      b_C->data[i8 + b_C->size[0] * i9] += C->data[i8 + C->size
                        [0] * k1] * onexi->data[k1 + onexi->size[0] * i9];
                    }
                  }
                }
              } else {
                k = C->size[1];
                nx = C->size[0];
                vstride = onexi->size[1];
                i8 = b_C->size[0] * b_C->size[1];
                b_C->size[0] = nx;
                b_C->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)b_C, i8, (int)sizeof
                                  (double));
                m = C->size[0];
                i8 = b_C->size[0] * b_C->size[1];
                emxEnsureCapacity((emxArray__common *)b_C, i8, (int)sizeof
                                  (double));
                loop_ub = b_C->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = b_C->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    b_C->data[i9 + b_C->size[0] * i8] = 0.0;
                  }
                }

                if ((C->size[0] == 0) || (onexi->size[1] == 0)) {
                } else {
                  vstride = C->size[0] * (onexi->size[1] - 1);
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      b_C->data[ic] = 0.0;
                    }

                    nx += m;
                  }

                  br = 0;
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (onexi->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          b_C->data[ic] += onexi->data[ib] * C->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                    nx += m;
                  }
                }
              }

              eye((double)p + 1.0, a);
              if ((C->size[1] == 1) || (ky1->size[0] == 1)) {
                i8 = f_y->size[0] * f_y->size[1];
                f_y->size[0] = C->size[0];
                f_y->size[1] = ky1->size[1];
                emxEnsureCapacity((emxArray__common *)f_y, i8, (int)sizeof
                                  (double));
                loop_ub = C->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = ky1->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    f_y->data[i8 + f_y->size[0] * i9] = 0.0;
                    vstride = C->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      f_y->data[i8 + f_y->size[0] * i9] += C->data[i8 + C->size
                        [0] * k1] * ky1->data[k1 + ky1->size[0] * i9];
                    }
                  }
                }
              } else {
                k = C->size[1];
                nx = C->size[0];
                vstride = ky1->size[1];
                i8 = f_y->size[0] * f_y->size[1];
                f_y->size[0] = nx;
                f_y->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)f_y, i8, (int)sizeof
                                  (double));
                m = C->size[0];
                i8 = f_y->size[0] * f_y->size[1];
                emxEnsureCapacity((emxArray__common *)f_y, i8, (int)sizeof
                                  (double));
                loop_ub = f_y->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = f_y->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    f_y->data[i9 + f_y->size[0] * i8] = 0.0;
                  }
                }

                if ((C->size[0] == 0) || (ky1->size[1] == 0)) {
                } else {
                  vstride = C->size[0] * (ky1->size[1] - 1);
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      f_y->data[ic] = 0.0;
                    }

                    nx += m;
                  }

                  br = 0;
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (ky1->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          f_y->data[ic] += ky1->data[ib] * C->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                    nx += m;
                  }
                }
              }

              i8 = g_C->size[0] * g_C->size[1];
              g_C->size[0] = b_C->size[0];
              g_C->size[1] = b_C->size[1];
              emxEnsureCapacity((emxArray__common *)g_C, i8, (int)sizeof(double));
              loop_ub = b_C->size[0] * b_C->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                g_C->data[i8] = b_C->data[i8] + a->data[i8] / (double)n;
              }

              mldivide(g_C, f_y, abi);
              i8 = abi->size[1];
              if ((i8 == 1) || (DD->size[0] == 1)) {
                loop_ub = abi->size[1];
                i8 = c_abi->size[0] * c_abi->size[1];
                c_abi->size[0] = p;
                c_abi->size[1] = loop_ub;
                emxEnsureCapacity((emxArray__common *)c_abi, i8, (int)sizeof
                                  (double));
                for (i8 = 0; i8 < loop_ub; i8++) {
                  for (i9 = 0; i9 < p; i9++) {
                    c_abi->data[i9 + c_abi->size[0] * i8] = abi->data[i9 +
                      abi->size[0] * i8];
                  }
                }

                i8 = g_y->size[0] * g_y->size[1];
                g_y->size[0] = c_abi->size[0];
                g_y->size[1] = DD->size[1];
                emxEnsureCapacity((emxArray__common *)g_y, i8, (int)sizeof
                                  (double));
                loop_ub = c_abi->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = DD->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    g_y->data[i8 + g_y->size[0] * i9] = 0.0;
                    vstride = c_abi->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      g_y->data[i8 + g_y->size[0] * i9] += c_abi->data[i8 +
                        c_abi->size[0] * k1] * DD->data[k1 + DD->size[0] * i9];
                    }
                  }
                }
              } else {
                i8 = abi->size[1];
                vstride = DD->size[1];
                i9 = g_y->size[0] * g_y->size[1];
                g_y->size[0] = p;
                g_y->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)g_y, i9, (int)sizeof
                                  (double));
                i9 = g_y->size[0] * g_y->size[1];
                emxEnsureCapacity((emxArray__common *)g_y, i9, (int)sizeof
                                  (double));
                loop_ub = g_y->size[1];
                for (i9 = 0; i9 < loop_ub; i9++) {
                  br = g_y->size[0];
                  for (k1 = 0; k1 < br; k1++) {
                    g_y->data[k1 + g_y->size[0] * i9] = 0.0;
                  }
                }

                if (DD->size[1] != 0) {
                  vstride = p * (DD->size[1] - 1);
                  for (nx = 0; nx <= vstride; nx += p) {
                    i9 = nx + p;
                    for (ic = nx; ic + 1 <= i9; ic++) {
                      g_y->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (nx = 0; nx <= vstride; nx += p) {
                    ar = 0;
                    i9 = br + i8;
                    for (ib = br; ib + 1 <= i9; ib++) {
                      if (DD->data[ib] != 0.0) {
                        ia = ar;
                        k1 = nx + p;
                        for (ic = nx; ic + 1 <= k1; ic++) {
                          ia++;
                          g_y->data[ic] += DD->data[ib] * abi->data[(ia - 1) % p
                            + abi->size[0] * ((ia - 1) / p)];
                        }
                      }

                      ar += p;
                    }

                    br += i8;
                  }
                }
              }

              loop_ub = abi->size[1];
              i8 = b->size[0] * b->size[1];
              b->size[0] = loop_ub;
              b->size[1] = p;
              emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
              for (i8 = 0; i8 < p; i8++) {
                for (i9 = 0; i9 < loop_ub; i9++) {
                  b->data[i9 + b->size[0] * i8] = abi->data[i8 + abi->size[0] *
                    i9];
                }
              }

              if ((g_y->size[1] == 1) || (b->size[0] == 1)) {
                i8 = c_C->size[0] * c_C->size[1];
                c_C->size[0] = g_y->size[0];
                c_C->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)c_C, i8, (int)sizeof
                                  (double));
                loop_ub = g_y->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = b->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    c_C->data[i8 + c_C->size[0] * i9] = 0.0;
                    vstride = g_y->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      c_C->data[i8 + c_C->size[0] * i9] += g_y->data[i8 +
                        g_y->size[0] * k1] * b->data[k1 + b->size[0] * i9];
                    }
                  }
                }
              } else {
                k = g_y->size[1];
                nx = g_y->size[0];
                vstride = b->size[1];
                i8 = c_C->size[0] * c_C->size[1];
                c_C->size[0] = nx;
                c_C->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)c_C, i8, (int)sizeof
                                  (double));
                m = g_y->size[0];
                i8 = c_C->size[0] * c_C->size[1];
                emxEnsureCapacity((emxArray__common *)c_C, i8, (int)sizeof
                                  (double));
                loop_ub = c_C->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = c_C->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    c_C->data[i9 + c_C->size[0] * i8] = 0.0;
                  }
                }

                vstride = g_y->size[0] * (b->size[1] - 1);
                for (nx = 0; nx <= vstride; nx += m) {
                  i8 = nx + m;
                  for (ic = nx; ic + 1 <= i8; ic++) {
                    c_C->data[ic] = 0.0;
                  }
                }

                br = 0;
                for (nx = 0; nx <= vstride; nx += m) {
                  ar = 0;
                  i8 = br + k;
                  for (ib = br; ib + 1 <= i8; ib++) {
                    if (b->data[ib] != 0.0) {
                      ia = ar;
                      i9 = nx + m;
                      for (ic = nx; ic + 1 <= i9; ic++) {
                        ia++;
                        c_C->data[ic] += b->data[ib] * g_y->data[ia - 1];
                      }
                    }

                    ar += m;
                  }

                  br += k;
                }
              }

              i8 = yi->size[0] * yi->size[1];
              emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
              nm1d2 = yi->size[0];
              vstride = yi->size[1];
              loop_ub = nm1d2 * vstride;
              for (i8 = 0; i8 < loop_ub; i8++) {
                yi->data[i8] += c_C->data[i8];
              }
            }

            eig(yi, Vc, Dc);
            i8 = B->size[0] * B->size[1];
            B->size[0] = Vc->size[0];
            B->size[1] = Vc->size[1];
            emxEnsureCapacity((emxArray__common *)B, i8, (int)sizeof(double));
            loop_ub = Vc->size[0] * Vc->size[1];
            for (i8 = 0; i8 < loop_ub; i8++) {
              B->data[i8] = Vc->data[i8].re;
            }

            /* Change for C */
            c_diag(Dc, R);
            i8 = s->size[0];
            s->size[0] = R->size[0];
            emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
            loop_ub = R->size[0];
            for (i8 = 0; i8 < loop_ub; i8++) {
              s->data[i8] = R->data[i8].re;
            }

            c_sort(s, iidx);
            i8 = s->size[0];
            s->size[0] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
            loop_ub = iidx->size[0];
            for (i8 = 0; i8 < loop_ub; i8++) {
              s->data[i8] = iidx->data[i8];
            }

            loop_ub = B->size[0];
            i8 = yi->size[0] * yi->size[1];
            yi->size[0] = loop_ub;
            yi->size[1] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
            br = iidx->size[0];
            for (i8 = 0; i8 < br; i8++) {
              for (i9 = 0; i9 < loop_ub; i9++) {
                yi->data[i9 + yi->size[0] * i8] = B->data[i9 + B->size[0] *
                  (iidx->data[i8] - 1)];
              }
            }

            if (1.0 > ip) {
              loop_ub = 0;
            } else {
              loop_ub = (int)ip;
            }

            br = B->size[0];
            i8 = B->size[0] * B->size[1];
            B->size[0] = br;
            B->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)B, i8, (int)sizeof(double));
            for (i8 = 0; i8 < loop_ub; i8++) {
              for (i9 = 0; i9 < br; i9++) {
                B->data[i9 + B->size[0] * i8] = yi->data[i9 + yi->size[0] * i8];
              }
            }
          } else {
            i8 = dd->size[0] * dd->size[1];
            dd->size[0] = (int)(b_m * (double)p);
            dd->size[1] = (int)(b_m * (double)p);
            emxEnsureCapacity((emxArray__common *)dd, i8, (int)sizeof(double));
            loop_ub = (int)(b_m * (double)p) * (int)(b_m * (double)p);
            for (i8 = 0; i8 < loop_ub; i8++) {
              dd->data[i8] = 0.0;
            }

            i8 = dc->size[0];
            dc->size[0] = (int)(b_m * (double)p);
            emxEnsureCapacity((emxArray__common *)dc, i8, (int)sizeof(double));
            loop_ub = (int)(b_m * (double)p);
            for (i8 = 0; i8 < loop_ub; i8++) {
              dc->data[i8] = 0.0;
            }

            i8 = D->size[0] * D->size[1];
            D->size[0] = (int)b_m;
            D->size[1] = (int)b_m;
            emxEnsureCapacity((emxArray__common *)D, i8, (int)sizeof(double));
            loop_ub = (int)b_m * (int)b_m;
            for (i8 = 0; i8 < loop_ub; i8++) {
              D->data[i8] = 0.0;
            }

            /* Change for C */
            for (j = 0; j < Ifast_idx_0; j++) {
              loop_ub = x->size[1];
              i8 = b_xfast->size[0] * b_xfast->size[1];
              b_xfast->size[0] = 1;
              b_xfast->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_xfast, i8, (int)sizeof
                                (double));
              for (i8 = 0; i8 < loop_ub; i8++) {
                b_xfast->data[b_xfast->size[0] * i8] = xfast->data[j +
                  xfast->size[0] * i8];
              }

              repmat(b_xfast, (double)n, xij);
              i8 = xij->size[0] * xij->size[1];
              xij->size[0] = x->size[0];
              xij->size[1] = x->size[1];
              emxEnsureCapacity((emxArray__common *)xij, i8, (int)sizeof(double));
              loop_ub = x->size[0] * x->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                xij->data[i8] = x->data[i8] - xij->data[i8];
              }

              loop_ub = V->size[1];
              i8 = b_U->size[0] * b_U->size[1];
              b_U->size[0] = 1;
              b_U->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_U, i8, (int)sizeof(double));
              for (i8 = 0; i8 < loop_ub; i8++) {
                b_U->data[b_U->size[0] * i8] = U->data[j + U->size[0] * i8];
              }

              repmat(b_U, (double)n, a);
              i8 = d_V->size[0] * d_V->size[1];
              d_V->size[0] = V->size[0];
              d_V->size[1] = V->size[1];
              emxEnsureCapacity((emxArray__common *)d_V, i8, (int)sizeof(double));
              loop_ub = V->size[0] * V->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                d_V->data[i8] = V->data[i8] - a->data[i8];
              }

              power(d_V, a);
              b_sum(a, dxij);
              i8 = s->size[0];
              s->size[0] = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
              loop_ub = dxij->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                s->data[i8] = dxij->data[i8];
              }

              f_sort(s);
              if ((h2 >= s->data[(int)(2.0 * b_m) - 1]) || rtIsNaN(s->data[(int)
                   (2.0 * b_m) - 1])) {
                pp = h2;
              } else {
                pp = s->data[(int)(2.0 * b_m) - 1];
              }

              i8 = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)dxij, i8, (int)sizeof(double));
              loop_ub = dxij->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                dxij->data[i8] = -dxij->data[i8] / pp;
              }

              b_exp(dxij);
              c_repmat(dxij, (double)p + 1.0, yj);
              if ((xij->size[1] == 1) || (B->size[0] == 1)) {
                i8 = d_y->size[0] * d_y->size[1];
                d_y->size[0] = xij->size[0];
                d_y->size[1] = B->size[1];
                emxEnsureCapacity((emxArray__common *)d_y, i8, (int)sizeof
                                  (double));
                loop_ub = xij->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = B->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    d_y->data[i8 + d_y->size[0] * i9] = 0.0;
                    vstride = xij->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      d_y->data[i8 + d_y->size[0] * i9] += xij->data[i8 +
                        xij->size[0] * k1] * B->data[k1 + B->size[0] * i9];
                    }
                  }
                }
              } else {
                k = xij->size[1];
                nx = xij->size[0];
                vstride = B->size[1];
                i8 = d_y->size[0] * d_y->size[1];
                d_y->size[0] = nx;
                d_y->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)d_y, i8, (int)sizeof
                                  (double));
                m = xij->size[0];
                i8 = d_y->size[0] * d_y->size[1];
                emxEnsureCapacity((emxArray__common *)d_y, i8, (int)sizeof
                                  (double));
                loop_ub = d_y->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = d_y->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    d_y->data[i9 + d_y->size[0] * i8] = 0.0;
                  }
                }

                if ((xij->size[0] == 0) || (B->size[1] == 0)) {
                } else {
                  vstride = xij->size[0] * (B->size[1] - 1);
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      d_y->data[ic] = 0.0;
                    }

                    nx += m;
                  }

                  br = 0;
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (B->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          d_y->data[ic] += B->data[ib] * xij->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                    nx += m;
                  }
                }
              }

              if (!((d_y->size[0] == 0) || (d_y->size[1] == 0))) {
                nx = d_y->size[0];
              } else if (!(n == 0)) {
                nx = n;
              } else {
                nx = d_y->size[0];
                if (!(nx >= 0)) {
                  nx = 0;
                }
              }

              empty_non_axis_sizes = (nx == 0);
              if (empty_non_axis_sizes || (!((d_y->size[0] == 0) || (d_y->size[1]
                     == 0)))) {
                nm1d2 = d_y->size[1];
              } else {
                nm1d2 = 0;
              }

              if (empty_non_axis_sizes || (!(n == 0))) {
                vstride = 1;
              } else {
                vstride = 0;
              }

              i8 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
              reshapes[1].f1->size[0] = nx;
              reshapes[1].f1->size[1] = vstride;
              emxEnsureCapacity((emxArray__common *)reshapes[1].f1, i8, (int)
                                sizeof(double));
              loop_ub = nx * vstride;
              for (i8 = 0; i8 < loop_ub; i8++) {
                reshapes[1].f1->data[i8] = 1.0;
              }

              i8 = onexi->size[0] * onexi->size[1];
              onexi->size[0] = nx;
              onexi->size[1] = nm1d2 + reshapes[1].f1->size[1];
              emxEnsureCapacity((emxArray__common *)onexi, i8, (int)sizeof
                                (double));
              for (i8 = 0; i8 < nm1d2; i8++) {
                for (i9 = 0; i9 < nx; i9++) {
                  onexi->data[i9 + onexi->size[0] * i8] = d_y->data[i9 + nx * i8];
                }
              }

              loop_ub = reshapes[1].f1->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = reshapes[1].f1->size[0];
                for (i9 = 0; i9 < br; i9++) {
                  onexi->data[i9 + onexi->size[0] * (i8 + nm1d2)] = reshapes[1].
                    f1->data[i9 + reshapes[1].f1->size[0] * i8];
                }
              }

              i8 = C->size[0] * C->size[1];
              C->size[0] = onexi->size[1];
              C->size[1] = onexi->size[0];
              emxEnsureCapacity((emxArray__common *)C, i8, (int)sizeof(double));
              loop_ub = onexi->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = onexi->size[1];
                for (i9 = 0; i9 < br; i9++) {
                  C->data[i9 + C->size[0] * i8] = onexi->data[i8 + onexi->size[0]
                    * i9] * yj->data[i8 + yj->size[0] * i9];
                }
              }

              /* abi = inv(xk*onexi+eye(size(B,2)+1)/n)*(xk*ky1)*ky2; */
              if ((C->size[1] == 1) || (onexi->size[0] == 1)) {
                i8 = d_C->size[0] * d_C->size[1];
                d_C->size[0] = C->size[0];
                d_C->size[1] = onexi->size[1];
                emxEnsureCapacity((emxArray__common *)d_C, i8, (int)sizeof
                                  (double));
                loop_ub = C->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = onexi->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    d_C->data[i8 + d_C->size[0] * i9] = 0.0;
                    vstride = C->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      d_C->data[i8 + d_C->size[0] * i9] += C->data[i8 + C->size
                        [0] * k1] * onexi->data[k1 + onexi->size[0] * i9];
                    }
                  }
                }
              } else {
                k = C->size[1];
                nx = C->size[0];
                vstride = onexi->size[1];
                i8 = d_C->size[0] * d_C->size[1];
                d_C->size[0] = nx;
                d_C->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)d_C, i8, (int)sizeof
                                  (double));
                m = C->size[0];
                i8 = d_C->size[0] * d_C->size[1];
                emxEnsureCapacity((emxArray__common *)d_C, i8, (int)sizeof
                                  (double));
                loop_ub = d_C->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = d_C->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    d_C->data[i9 + d_C->size[0] * i8] = 0.0;
                  }
                }

                if ((C->size[0] == 0) || (onexi->size[1] == 0)) {
                } else {
                  vstride = C->size[0] * (onexi->size[1] - 1);
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      d_C->data[ic] = 0.0;
                    }

                    nx += m;
                  }

                  br = 0;
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (onexi->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          d_C->data[ic] += onexi->data[ib] * C->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                    nx += m;
                  }
                }
              }

              eye((double)B->size[1] + 1.0, a);
              if ((C->size[1] == 1) || (ky1->size[0] == 1)) {
                i8 = e_y->size[0] * e_y->size[1];
                e_y->size[0] = C->size[0];
                e_y->size[1] = ky1->size[1];
                emxEnsureCapacity((emxArray__common *)e_y, i8, (int)sizeof
                                  (double));
                loop_ub = C->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = ky1->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    e_y->data[i8 + e_y->size[0] * i9] = 0.0;
                    vstride = C->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      e_y->data[i8 + e_y->size[0] * i9] += C->data[i8 + C->size
                        [0] * k1] * ky1->data[k1 + ky1->size[0] * i9];
                    }
                  }
                }
              } else {
                k = C->size[1];
                nx = C->size[0];
                vstride = ky1->size[1];
                i8 = e_y->size[0] * e_y->size[1];
                e_y->size[0] = nx;
                e_y->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)e_y, i8, (int)sizeof
                                  (double));
                m = C->size[0];
                i8 = e_y->size[0] * e_y->size[1];
                emxEnsureCapacity((emxArray__common *)e_y, i8, (int)sizeof
                                  (double));
                loop_ub = e_y->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = e_y->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    e_y->data[i9 + e_y->size[0] * i8] = 0.0;
                  }
                }

                if ((C->size[0] == 0) || (ky1->size[1] == 0)) {
                } else {
                  vstride = C->size[0] * (ky1->size[1] - 1);
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      e_y->data[ic] = 0.0;
                    }

                    nx += m;
                  }

                  br = 0;
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (ky1->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          e_y->data[ic] += ky1->data[ib] * C->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                    nx += m;
                  }
                }
              }

              i8 = f_C->size[0] * f_C->size[1];
              f_C->size[0] = d_C->size[0];
              f_C->size[1] = d_C->size[1];
              emxEnsureCapacity((emxArray__common *)f_C, i8, (int)sizeof(double));
              loop_ub = d_C->size[0] * d_C->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                f_C->data[i8] = d_C->data[i8] + a->data[i8] / (double)n;
              }

              mldivide(f_C, e_y, a);
              if ((a->size[1] == 1) || (ky2->size[0] == 1)) {
                i8 = abi->size[0] * abi->size[1];
                abi->size[0] = a->size[0];
                abi->size[1] = ky2->size[1];
                emxEnsureCapacity((emxArray__common *)abi, i8, (int)sizeof
                                  (double));
                loop_ub = a->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = ky2->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    abi->data[i8 + abi->size[0] * i9] = 0.0;
                    vstride = a->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      abi->data[i8 + abi->size[0] * i9] += a->data[i8 + a->size
                        [0] * k1] * ky2->data[k1 + ky2->size[0] * i9];
                    }
                  }
                }
              } else {
                k = a->size[1];
                nx = a->size[0];
                vstride = ky2->size[1];
                i8 = abi->size[0] * abi->size[1];
                abi->size[0] = nx;
                abi->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)abi, i8, (int)sizeof
                                  (double));
                m = a->size[0];
                i8 = abi->size[0] * abi->size[1];
                emxEnsureCapacity((emxArray__common *)abi, i8, (int)sizeof
                                  (double));
                loop_ub = abi->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = abi->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    abi->data[i9 + abi->size[0] * i8] = 0.0;
                  }
                }

                if ((a->size[0] == 0) || (ky2->size[1] == 0)) {
                } else {
                  vstride = a->size[0] * (ky2->size[1] - 1);
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      abi->data[ic] = 0.0;
                    }

                    nx += m;
                  }

                  br = 0;
                  nx = 0;
                  while ((m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (ky2->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          abi->data[ic] += ky2->data[ib] * a->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                    nx += m;
                  }
                }
              }

              i8 = yi->size[0] * yi->size[1];
              yi->size[0] = xij->size[1];
              yi->size[1] = xij->size[0];
              emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
              loop_ub = xij->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = xij->size[1];
                for (i9 = 0; i9 < br; i9++) {
                  yi->data[i9 + yi->size[0] * i8] = xij->data[i8 + xij->size[0] *
                    i9] * yj->data[i8 + yj->size[0] * i9];
                }
              }

              loop_ub = abi->size[1];
              i8 = b_abi->size[0] * b_abi->size[1];
              b_abi->size[0] = 1;
              b_abi->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_abi, i8, (int)sizeof
                                (double));
              for (i8 = 0; i8 < loop_ub; i8++) {
                b_abi->data[b_abi->size[0] * i8] = abi->data[((int)(b_m + 1.0) +
                  abi->size[0] * i8) - 1];
              }

              repmat(b_abi, (double)n, b);
              i8 = b->size[0] * b->size[1];
              b->size[0] = ky->size[0];
              b->size[1] = ky->size[1];
              emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
              loop_ub = ky->size[0] * ky->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                b->data[i8] = ky->data[i8] - b->data[i8];
              }

              if ((yi->size[1] == 1) || (b->size[0] == 1)) {
                i8 = kxijy->size[0] * kxijy->size[1];
                kxijy->size[0] = yi->size[0];
                kxijy->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)kxijy, i8, (int)sizeof
                                  (double));
                loop_ub = yi->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = b->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    kxijy->data[i8 + kxijy->size[0] * i9] = 0.0;
                    vstride = yi->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      kxijy->data[i8 + kxijy->size[0] * i9] += yi->data[i8 +
                        yi->size[0] * k1] * b->data[k1 + b->size[0] * i9];
                    }
                  }
                }
              } else {
                k = yi->size[1];
                nx = yi->size[0];
                vstride = b->size[1];
                i8 = kxijy->size[0] * kxijy->size[1];
                kxijy->size[0] = nx;
                kxijy->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)kxijy, i8, (int)sizeof
                                  (double));
                m = yi->size[0];
                i8 = kxijy->size[0] * kxijy->size[1];
                emxEnsureCapacity((emxArray__common *)kxijy, i8, (int)sizeof
                                  (double));
                loop_ub = kxijy->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = kxijy->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    kxijy->data[i9 + kxijy->size[0] * i8] = 0.0;
                  }
                }

                if (b->size[1] != 0) {
                  vstride = yi->size[0] * (b->size[1] - 1);
                  for (nx = 0; nx <= vstride; nx += m) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      kxijy->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (nx = 0; nx <= vstride; nx += m) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (b->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          kxijy->data[ic] += b->data[ib] * yi->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                  }
                }
              }

              if ((yi->size[1] == 1) || (xij->size[0] == 1)) {
                i8 = ddx->size[0] * ddx->size[1];
                ddx->size[0] = yi->size[0];
                ddx->size[1] = xij->size[1];
                emxEnsureCapacity((emxArray__common *)ddx, i8, (int)sizeof
                                  (double));
                loop_ub = yi->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = xij->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    ddx->data[i8 + ddx->size[0] * i9] = 0.0;
                    vstride = yi->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      ddx->data[i8 + ddx->size[0] * i9] += yi->data[i8 +
                        yi->size[0] * k1] * xij->data[k1 + xij->size[0] * i9];
                    }
                  }
                }
              } else {
                k = yi->size[1];
                nx = yi->size[0];
                vstride = xij->size[1];
                i8 = ddx->size[0] * ddx->size[1];
                ddx->size[0] = nx;
                ddx->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)ddx, i8, (int)sizeof
                                  (double));
                m = yi->size[0];
                i8 = ddx->size[0] * ddx->size[1];
                emxEnsureCapacity((emxArray__common *)ddx, i8, (int)sizeof
                                  (double));
                loop_ub = ddx->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = ddx->size[0];
                  for (i9 = 0; i9 < br; i9++) {
                    ddx->data[i9 + ddx->size[0] * i8] = 0.0;
                  }
                }

                if (xij->size[1] != 0) {
                  vstride = yi->size[0] * (xij->size[1] - 1);
                  for (nx = 0; nx <= vstride; nx += m) {
                    i8 = nx + m;
                    for (ic = nx; ic + 1 <= i8; ic++) {
                      ddx->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (nx = 0; nx <= vstride; nx += m) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (xij->data[ib] != 0.0) {
                        ia = ar;
                        i9 = nx + m;
                        for (ic = nx; ic + 1 <= i9; ic++) {
                          ia++;
                          ddx->data[ic] += xij->data[ib] * yi->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                  }
                }
              }

              for (k1 = 0; k1 < (int)b_m; k1++) {
                pp = ((1.0 + (double)k1) - 1.0) * (double)p + 1.0;
                d = (1.0 + (double)k1) * (double)p;
                if (d < pp) {
                  i8 = qx->size[0] * qx->size[1];
                  qx->size[0] = 1;
                  qx->size[1] = 0;
                  emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof
                                    (double));
                } else if (pp == pp) {
                  i8 = qx->size[0] * qx->size[1];
                  qx->size[0] = 1;
                  qx->size[1] = (int)(d - pp) + 1;
                  emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof
                                    (double));
                  loop_ub = (int)(d - pp);
                  for (i8 = 0; i8 <= loop_ub; i8++) {
                    qx->data[qx->size[0] * i8] = pp + (double)i8;
                  }
                } else {
                  ndbl = std::floor((d - pp) + 0.5);
                  apnd = pp + ndbl;
                  cdiff = apnd - d;
                  absa = std::abs(pp);
                  absb = std::abs(d);
                  if (absa >= absb) {
                    absb = absa;
                  }

                  if (std::abs(cdiff) < 4.4408920985006262E-16 * absb) {
                    ndbl++;
                    apnd = d;
                  } else if (cdiff > 0.0) {
                    apnd = pp + (ndbl - 1.0);
                  } else {
                    ndbl++;
                  }

                  if (ndbl >= 0.0) {
                    nx = (int)ndbl;
                  } else {
                    nx = 0;
                  }

                  i8 = qx->size[0] * qx->size[1];
                  qx->size[0] = 1;
                  qx->size[1] = nx;
                  emxEnsureCapacity((emxArray__common *)qx, i8, (int)sizeof
                                    (double));
                  if (nx > 0) {
                    qx->data[0] = pp;
                    if (nx > 1) {
                      qx->data[nx - 1] = apnd;
                      nm1d2 = (nx - 1) / 2;
                      for (k = 1; k < nm1d2; k++) {
                        qx->data[k] = pp + (double)k;
                        qx->data[(nx - k) - 1] = apnd - (double)k;
                      }

                      if (nm1d2 << 1 == nx - 1) {
                        qx->data[nm1d2] = (pp + apnd) / 2.0;
                      } else {
                        qx->data[nm1d2] = pp + (double)nm1d2;
                        qx->data[nm1d2 + 1] = apnd - (double)nm1d2;
                      }
                    }
                  }
                }

                loop_ub = abi->size[1];
                i8 = s->size[0];
                s->size[0] = loop_ub;
                emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
                for (i8 = 0; i8 < loop_ub; i8++) {
                  s->data[i8] = abi->data[k1 + abi->size[0] * i8];
                }

                if ((kxijy->size[1] == 1) || (s->size[0] == 1)) {
                  i8 = e_C->size[0];
                  e_C->size[0] = kxijy->size[0];
                  emxEnsureCapacity((emxArray__common *)e_C, i8, (int)sizeof
                                    (double));
                  loop_ub = kxijy->size[0];
                  for (i8 = 0; i8 < loop_ub; i8++) {
                    e_C->data[i8] = 0.0;
                    br = kxijy->size[1];
                    for (i9 = 0; i9 < br; i9++) {
                      e_C->data[i8] += kxijy->data[i8 + kxijy->size[0] * i9] *
                        s->data[i9];
                    }
                  }
                } else {
                  k = kxijy->size[1];
                  a_idx_0 = (unsigned int)kxijy->size[0];
                  i8 = e_C->size[0];
                  e_C->size[0] = (int)a_idx_0;
                  emxEnsureCapacity((emxArray__common *)e_C, i8, (int)sizeof
                                    (double));
                  m = kxijy->size[0];
                  vstride = e_C->size[0];
                  i8 = e_C->size[0];
                  e_C->size[0] = vstride;
                  emxEnsureCapacity((emxArray__common *)e_C, i8, (int)sizeof
                                    (double));
                  for (i8 = 0; i8 < vstride; i8++) {
                    e_C->data[i8] = 0.0;
                  }

                  nx = 0;
                  while (nx <= 0) {
                    for (ic = 1; ic <= m; ic++) {
                      e_C->data[ic - 1] = 0.0;
                    }

                    nx = m;
                  }

                  br = 0;
                  nx = 0;
                  while (nx <= 0) {
                    ar = 0;
                    i8 = br + k;
                    for (ib = br; ib + 1 <= i8; ib++) {
                      if (s->data[ib] != 0.0) {
                        ia = ar;
                        for (ic = 0; ic + 1 <= m; ic++) {
                          ia++;
                          e_C->data[ic] += s->data[ib] * kxijy->data[ia - 1];
                        }
                      }

                      ar += m;
                    }

                    br += k;
                    nx = m;
                  }
                }

                i8 = r1->size[0] * r1->size[1];
                r1->size[0] = 1;
                r1->size[1] = qx->size[1];
                emxEnsureCapacity((emxArray__common *)r1, i8, (int)sizeof(int));
                loop_ub = qx->size[0] * qx->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  r1->data[i8] = (int)qx->data[i8];
                }

                i8 = b_dc->size[0];
                b_dc->size[0] = qx->size[1];
                emxEnsureCapacity((emxArray__common *)b_dc, i8, (int)sizeof
                                  (double));
                loop_ub = qx->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  b_dc->data[i8] = dc->data[(int)qx->data[qx->size[0] * i8] - 1]
                    + e_C->data[i8];
                }

                loop_ub = r1->size[1];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  dc->data[r1->data[r1->size[0] * i8] - 1] = b_dc->data[(*(int (*)
                    [2])r1->size)[0] * i8];
                }
              }

              if (1.0 > b_m) {
                loop_ub = 0;
                br = 0;
              } else {
                loop_ub = (int)b_m;
                br = (int)b_m;
              }

              vstride = abi->size[1];
              i8 = a->size[0] * a->size[1];
              a->size[0] = loop_ub;
              a->size[1] = vstride;
              emxEnsureCapacity((emxArray__common *)a, i8, (int)sizeof(double));
              for (i8 = 0; i8 < vstride; i8++) {
                for (i9 = 0; i9 < loop_ub; i9++) {
                  a->data[i9 + a->size[0] * i8] = abi->data[i9 + abi->size[0] *
                    i8];
                }
              }

              vstride = abi->size[1];
              i8 = b->size[0] * b->size[1];
              b->size[0] = vstride;
              b->size[1] = br;
              emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
              for (i8 = 0; i8 < br; i8++) {
                for (i9 = 0; i9 < vstride; i9++) {
                  b->data[i9 + b->size[0] * i8] = abi->data[i8 + abi->size[0] *
                    i9];
                }
              }

              i8 = abi->size[1];
              if ((i8 == 1) || (b->size[0] == 1)) {
                i8 = tmp->size[0] * tmp->size[1];
                tmp->size[0] = a->size[0];
                tmp->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)tmp, i8, (int)sizeof
                                  (double));
                loop_ub = a->size[0];
                for (i8 = 0; i8 < loop_ub; i8++) {
                  br = b->size[1];
                  for (i9 = 0; i9 < br; i9++) {
                    tmp->data[i8 + tmp->size[0] * i9] = 0.0;
                    vstride = a->size[1];
                    for (k1 = 0; k1 < vstride; k1++) {
                      tmp->data[i8 + tmp->size[0] * i9] += a->data[i8 + a->size
                        [0] * k1] * b->data[k1 + b->size[0] * i9];
                    }
                  }
                }
              } else {
                i8 = abi->size[1];
                vstride = b->size[1];
                i9 = tmp->size[0] * tmp->size[1];
                tmp->size[0] = loop_ub;
                tmp->size[1] = vstride;
                emxEnsureCapacity((emxArray__common *)tmp, i9, (int)sizeof
                                  (double));
                i9 = tmp->size[0] * tmp->size[1];
                emxEnsureCapacity((emxArray__common *)tmp, i9, (int)sizeof
                                  (double));
                br = tmp->size[1];
                for (i9 = 0; i9 < br; i9++) {
                  vstride = tmp->size[0];
                  for (k1 = 0; k1 < vstride; k1++) {
                    tmp->data[k1 + tmp->size[0] * i9] = 0.0;
                  }
                }

                if ((loop_ub == 0) || (b->size[1] == 0)) {
                } else {
                  vstride = loop_ub * (b->size[1] - 1);
                  nx = 0;
                  while ((loop_ub > 0) && (nx <= vstride)) {
                    i9 = nx + loop_ub;
                    for (ic = nx; ic + 1 <= i9; ic++) {
                      tmp->data[ic] = 0.0;
                    }

                    nx += loop_ub;
                  }

                  br = 0;
                  nx = 0;
                  while ((loop_ub > 0) && (nx <= vstride)) {
                    ar = 0;
                    i9 = br + i8;
                    for (ib = br; ib + 1 <= i9; ib++) {
                      if (b->data[ib] != 0.0) {
                        ia = ar;
                        k1 = nx + loop_ub;
                        for (ic = nx; ic + 1 <= k1; ic++) {
                          ia++;
                          tmp->data[ic] += b->data[ib] * a->data[ia - 1];
                        }
                      }

                      ar += loop_ub;
                    }

                    br += i8;
                    nx += loop_ub;
                  }
                }
              }

              kron(tmp, ddx, a);
              i8 = dd->size[0] * dd->size[1];
              emxEnsureCapacity((emxArray__common *)dd, i8, (int)sizeof(double));
              nm1d2 = dd->size[0];
              vstride = dd->size[1];
              loop_ub = nm1d2 * vstride;
              for (i8 = 0; i8 < loop_ub; i8++) {
                dd->data[i8] += a->data[i8];
              }

              i8 = D->size[0] * D->size[1];
              emxEnsureCapacity((emxArray__common *)D, i8, (int)sizeof(double));
              nm1d2 = D->size[0];
              vstride = D->size[1];
              loop_ub = nm1d2 * vstride;
              for (i8 = 0; i8 < loop_ub; i8++) {
                D->data[i8] += tmp->data[i8];
              }
            }

            eye((double)dc->size[0], a);
            i8 = b_dd->size[0] * b_dd->size[1];
            b_dd->size[0] = dd->size[0];
            b_dd->size[1] = dd->size[1];
            emxEnsureCapacity((emxArray__common *)b_dd, i8, (int)sizeof(double));
            loop_ub = dd->size[0] * dd->size[1];
            for (i8 = 0; i8 < loop_ub; i8++) {
              b_dd->data[i8] = dd->data[i8] + a->data[i8] / (double)n;
            }

            inv(b_dd, a);
            if ((a->size[1] == 1) || (dc->size[0] == 1)) {
              i8 = d_B->size[0];
              d_B->size[0] = a->size[0];
              emxEnsureCapacity((emxArray__common *)d_B, i8, (int)sizeof(double));
              loop_ub = a->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                d_B->data[i8] = 0.0;
                br = a->size[1];
                for (i9 = 0; i9 < br; i9++) {
                  d_B->data[i8] += a->data[i8 + a->size[0] * i9] * dc->data[i9];
                }
              }
            } else {
              k = a->size[1];
              a_idx_0 = (unsigned int)a->size[0];
              i8 = d_B->size[0];
              d_B->size[0] = (int)a_idx_0;
              emxEnsureCapacity((emxArray__common *)d_B, i8, (int)sizeof(double));
              m = a->size[0];
              vstride = d_B->size[0];
              i8 = d_B->size[0];
              d_B->size[0] = vstride;
              emxEnsureCapacity((emxArray__common *)d_B, i8, (int)sizeof(double));
              for (i8 = 0; i8 < vstride; i8++) {
                d_B->data[i8] = 0.0;
              }

              if (a->size[0] != 0) {
                nx = 0;
                while ((m > 0) && (nx <= 0)) {
                  for (ic = 1; ic <= m; ic++) {
                    d_B->data[ic - 1] = 0.0;
                  }

                  nx = m;
                }

                br = 0;
                nx = 0;
                while ((m > 0) && (nx <= 0)) {
                  ar = 0;
                  i8 = br + k;
                  for (ib = br; ib + 1 <= i8; ib++) {
                    if (dc->data[ib] != 0.0) {
                      ia = ar;
                      for (ic = 0; ic + 1 <= m; ic++) {
                        ia++;
                        d_B->data[ic] += dc->data[ib] * a->data[ia - 1];
                      }
                    }

                    ar += m;
                  }

                  br += k;
                  nx = m;
                }
              }
            }

            i8 = B0->size[0] * B0->size[1];
            B0->size[0] = p;
            B0->size[1] = (int)b_m;
            emxEnsureCapacity((emxArray__common *)B0, i8, (int)sizeof(double));
            for (k = 0; k + 1 <= d_B->size[0]; k++) {
              B0->data[k] = d_B->data[k];
            }

            if ((B0->size[1] == 1) || (D->size[0] == 1)) {
              i8 = b_y->size[0] * b_y->size[1];
              b_y->size[0] = B0->size[0];
              b_y->size[1] = D->size[1];
              emxEnsureCapacity((emxArray__common *)b_y, i8, (int)sizeof(double));
              loop_ub = B0->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = D->size[1];
                for (i9 = 0; i9 < br; i9++) {
                  b_y->data[i8 + b_y->size[0] * i9] = 0.0;
                  vstride = B0->size[1];
                  for (k1 = 0; k1 < vstride; k1++) {
                    b_y->data[i8 + b_y->size[0] * i9] += B0->data[i8 + B0->size
                      [0] * k1] * D->data[k1 + D->size[0] * i9];
                  }
                }
              }
            } else {
              k = B0->size[1];
              nx = B0->size[0];
              vstride = D->size[1];
              i8 = b_y->size[0] * b_y->size[1];
              b_y->size[0] = nx;
              b_y->size[1] = vstride;
              emxEnsureCapacity((emxArray__common *)b_y, i8, (int)sizeof(double));
              m = B0->size[0];
              i8 = b_y->size[0] * b_y->size[1];
              emxEnsureCapacity((emxArray__common *)b_y, i8, (int)sizeof(double));
              loop_ub = b_y->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = b_y->size[0];
                for (i9 = 0; i9 < br; i9++) {
                  b_y->data[i9 + b_y->size[0] * i8] = 0.0;
                }
              }

              if (D->size[1] != 0) {
                vstride = B0->size[0] * (D->size[1] - 1);
                for (nx = 0; nx <= vstride; nx += m) {
                  i8 = nx + m;
                  for (ic = nx; ic + 1 <= i8; ic++) {
                    b_y->data[ic] = 0.0;
                  }
                }

                br = 0;
                for (nx = 0; nx <= vstride; nx += m) {
                  ar = 0;
                  i8 = br + k;
                  for (ib = br; ib + 1 <= i8; ib++) {
                    if (D->data[ib] != 0.0) {
                      ia = ar;
                      i9 = nx + m;
                      for (ic = nx; ic + 1 <= i9; ic++) {
                        ia++;
                        b_y->data[ic] += D->data[ib] * B0->data[ia - 1];
                      }
                    }

                    ar += m;
                  }

                  br += k;
                }
              }
            }

            i8 = b->size[0] * b->size[1];
            b->size[0] = B0->size[1];
            b->size[1] = B0->size[0];
            emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
            loop_ub = B0->size[0];
            for (i8 = 0; i8 < loop_ub; i8++) {
              br = B0->size[1];
              for (i9 = 0; i9 < br; i9++) {
                b->data[i9 + b->size[0] * i8] = B0->data[i8 + B0->size[0] * i9];
              }
            }

            if ((b_y->size[1] == 1) || (b->size[0] == 1)) {
              i8 = c_y->size[0] * c_y->size[1];
              c_y->size[0] = b_y->size[0];
              c_y->size[1] = b->size[1];
              emxEnsureCapacity((emxArray__common *)c_y, i8, (int)sizeof(double));
              loop_ub = b_y->size[0];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = b->size[1];
                for (i9 = 0; i9 < br; i9++) {
                  c_y->data[i8 + c_y->size[0] * i9] = 0.0;
                  vstride = b_y->size[1];
                  for (k1 = 0; k1 < vstride; k1++) {
                    c_y->data[i8 + c_y->size[0] * i9] += b_y->data[i8 +
                      b_y->size[0] * k1] * b->data[k1 + b->size[0] * i9];
                  }
                }
              }
            } else {
              k = b_y->size[1];
              nx = b_y->size[0];
              vstride = b->size[1];
              i8 = c_y->size[0] * c_y->size[1];
              c_y->size[0] = nx;
              c_y->size[1] = vstride;
              emxEnsureCapacity((emxArray__common *)c_y, i8, (int)sizeof(double));
              m = b_y->size[0];
              i8 = c_y->size[0] * c_y->size[1];
              emxEnsureCapacity((emxArray__common *)c_y, i8, (int)sizeof(double));
              loop_ub = c_y->size[1];
              for (i8 = 0; i8 < loop_ub; i8++) {
                br = c_y->size[0];
                for (i9 = 0; i9 < br; i9++) {
                  c_y->data[i9 + c_y->size[0] * i8] = 0.0;
                }
              }

              vstride = b_y->size[0] * (b->size[1] - 1);
              for (nx = 0; nx <= vstride; nx += m) {
                i8 = nx + m;
                for (ic = nx; ic + 1 <= i8; ic++) {
                  c_y->data[ic] = 0.0;
                }
              }

              br = 0;
              for (nx = 0; nx <= vstride; nx += m) {
                ar = 0;
                i8 = br + k;
                for (ib = br; ib + 1 <= i8; ib++) {
                  if (b->data[ib] != 0.0) {
                    ia = ar;
                    i9 = nx + m;
                    for (ic = nx; ic + 1 <= i9; ic++) {
                      ia++;
                      c_y->data[ic] += b->data[ib] * b_y->data[ia - 1];
                    }
                  }

                  ar += m;
                }

                br += k;
              }
            }

            eig(c_y, Vc, Dc);
            i8 = B->size[0] * B->size[1];
            B->size[0] = Vc->size[0];
            B->size[1] = Vc->size[1];
            emxEnsureCapacity((emxArray__common *)B, i8, (int)sizeof(double));
            loop_ub = Vc->size[0] * Vc->size[1];
            for (i8 = 0; i8 < loop_ub; i8++) {
              B->data[i8] = Vc->data[i8].re;
            }

            /* Change for C */
            c_diag(Dc, R);
            i8 = s->size[0];
            s->size[0] = R->size[0];
            emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
            loop_ub = R->size[0];
            for (i8 = 0; i8 < loop_ub; i8++) {
              s->data[i8] = R->data[i8].re;
            }

            c_sort(s, iidx);
            i8 = s->size[0];
            s->size[0] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)s, i8, (int)sizeof(double));
            loop_ub = iidx->size[0];
            for (i8 = 0; i8 < loop_ub; i8++) {
              s->data[i8] = iidx->data[i8];
            }

            loop_ub = B->size[0];
            i8 = yi->size[0] * yi->size[1];
            yi->size[0] = loop_ub;
            yi->size[1] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)yi, i8, (int)sizeof(double));
            br = iidx->size[0];
            for (i8 = 0; i8 < br; i8++) {
              for (i9 = 0; i9 < loop_ub; i9++) {
                yi->data[i9 + yi->size[0] * i8] = B->data[i9 + B->size[0] *
                  (iidx->data[i8] - 1)];
              }
            }

            if (1.0 > ip) {
              loop_ub = 0;
            } else {
              loop_ub = (int)ip;
            }

            br = B->size[0];
            i8 = B->size[0] * B->size[1];
            B->size[0] = br;
            B->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)B, i8, (int)sizeof(double));
            for (i8 = 0; i8 < loop_ub; i8++) {
              for (i9 = 0; i9 < br; i9++) {
                B->data[i9 + B->size[0] * i8] = yi->data[i9 + yi->size[0] * i8];
              }
            }
          }
        }

        nm1d2 = (int)which_dim->data[iter];
        loop_ub = B->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          br = B->size[0];
          for (i9 = 0; i9 < br; i9++) {
            BB->data[(i9 + BB->size[0] * i8) + BB->size[0] * BB->size[1] *
              (nm1d2 - 1)] = B->data[i9 + B->size[0] * i8];
          }
        }

        if (d_strcmp(method)) {
          i8 = (int)((1.0 + (-1.0 - ((double)p - 1.0))) / -1.0);
          for (nm1d2 = 0; nm1d2 < i8; nm1d2++) {
            vstride = (p - nm1d2) - 1;
            if (1 > vstride) {
              loop_ub = -1;
            } else {
              loop_ub = vstride - 1;
            }

            br = yi->size[0] - 1;
            for (i9 = 0; i9 <= loop_ub; i9++) {
              for (k1 = 0; k1 <= br; k1++) {
                BB->data[(k1 + BB->size[0] * i9) + BB->size[0] * BB->size[1] *
                  (vstride - 1)] = yi->data[k1 + yi->size[0] * i9];
              }
            }
          }

          exitg1 = true;
        } else {
          /* B = B(:,1:(ip-1)); */
          iter++;
        }
      }

      emxFree_int32_T(&e_Ifast);
      emxFree_int32_T(&d_Ifast);
      emxFree_real_T(&c_abi);
      emxFree_real_T(&b_dc);
      emxFree_int32_T(&c_Ifast);
      emxFree_int32_T(&b_Ifast);
      emxFree_real_T(&e_B);
      emxFree_real_T(&e_V);
      emxFree_real_T(&g_C);
      emxFree_real_T(&b_xfast);
      emxFree_real_T(&b_U);
      emxFree_real_T(&d_V);
      emxFree_real_T(&f_C);
      emxFree_real_T(&b_abi);
      emxFree_real_T(&b_dd);
      emxFree_real_T(&g_y);
      emxFree_real_T(&f_y);
      emxFree_real_T(&e_y);
      emxFreeMatrix_cell_wrap_0(reshapes);
      emxFree_real_T(&d_y);
      emxFree_real_T(&c_y);
      emxFree_real_T(&b_y);
      emxFree_int32_T(&r1);
      emxFree_real_T(&e_C);
      emxFree_real_T(&d_C);
      emxFree_real_T(&c_C);
      emxFree_real_T(&b_C);
      emxFree_creal_T(&R);
      emxFree_real_T(&d_B);
      emxFree_real_T(&B0);
      emxFree_real_T(&tmp);
      emxFree_real_T(&ddx);
      emxFree_real_T(&kxijy);
      emxFree_real_T(&dd);
      emxFree_real_T(&abi);
      emxFree_real_T(&dxij);
      emxFree_real_T(&xij);
      emxFree_real_T(&xfast);
      emxFree_real_T(&Ifast);
      emxFree_real_T(&D);
    }

    emxFree_real_T(&onexi);
    nm1d2 = 0;
    emxInit_real_T1(&b_B, 1);
    emxInit_real_T1(&c_B, 1);
    while (nm1d2 <= p - 1) {
      loop_ub = BB->size[0];
      br = BB->size[1];
      i8 = b->size[0] * b->size[1];
      b->size[0] = loop_ub;
      b->size[1] = br;
      emxEnsureCapacity((emxArray__common *)b, i8, (int)sizeof(double));
      for (i8 = 0; i8 < br; i8++) {
        for (i9 = 0; i9 < loop_ub; i9++) {
          b->data[i9 + b->size[0] * i8] = BB->data[(i9 + BB->size[0] * i8) +
            BB->size[0] * BB->size[1] * nm1d2];
        }
      }

      guard1 = false;
      if (ss->size[1] == 1) {
        guard1 = true;
      } else {
        i8 = BB->size[0];
        if (i8 == 1) {
          guard1 = true;
        } else {
          k = ss->size[1];
          i8 = BB->size[1];
          nx = ss->size[0];
          i9 = B->size[0] * B->size[1];
          B->size[0] = nx;
          B->size[1] = i8;
          emxEnsureCapacity((emxArray__common *)B, i9, (int)sizeof(double));
          m = ss->size[0];
          i8 = B->size[0] * B->size[1];
          emxEnsureCapacity((emxArray__common *)B, i8, (int)sizeof(double));
          loop_ub = B->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            br = B->size[0];
            for (i9 = 0; i9 < br; i9++) {
              B->data[i9 + B->size[0] * i8] = 0.0;
            }
          }

          if (ss->size[0] == 0) {
          } else {
            i8 = BB->size[1];
            if (i8 == 0) {
            } else {
              i8 = BB->size[1] - 1;
              vstride = ss->size[0] * i8;
              nx = 0;
              while ((m > 0) && (nx <= vstride)) {
                i8 = nx + m;
                for (ic = nx; ic + 1 <= i8; ic++) {
                  B->data[ic] = 0.0;
                }

                nx += m;
              }

              br = 0;
              nx = 0;
              while ((m > 0) && (nx <= vstride)) {
                ar = 0;
                i8 = br + k;
                for (ib = br; ib + 1 <= i8; ib++) {
                  if (b->data[ib] != 0.0) {
                    ia = ar;
                    i9 = nx + m;
                    for (ic = nx; ic + 1 <= i9; ic++) {
                      ia++;
                      B->data[ic] += b->data[ib] * ss->data[ia - 1];
                    }
                  }

                  ar += m;
                }

                br += k;
                nx += m;
              }
            }
          }
        }
      }

      if (guard1) {
        i8 = B->size[0] * B->size[1];
        B->size[0] = ss->size[0];
        B->size[1] = b->size[1];
        emxEnsureCapacity((emxArray__common *)B, i8, (int)sizeof(double));
        loop_ub = ss->size[0];
        for (i8 = 0; i8 < loop_ub; i8++) {
          br = b->size[1];
          for (i9 = 0; i9 < br; i9++) {
            B->data[i8 + B->size[0] * i9] = 0.0;
            vstride = ss->size[1];
            for (k1 = 0; k1 < vstride; k1++) {
              B->data[i8 + B->size[0] * i9] += ss->data[i8 + ss->size[0] * k1] *
                b->data[k1 + b->size[0] * i9];
            }
          }
        }
      }

      i8 = B->size[1];
      for (nx = 0; nx < i8; nx++) {
        loop_ub = B->size[0];
        i9 = b_B->size[0];
        b_B->size[0] = loop_ub;
        emxEnsureCapacity((emxArray__common *)b_B, i9, (int)sizeof(double));
        for (i9 = 0; i9 < loop_ub; i9++) {
          b_B->data[i9] = B->data[i9 + B->size[0] * nx];
        }

        pp = norm(b_B) + 1.0E-20;
        vstride = B->size[0];
        i9 = c_B->size[0];
        c_B->size[0] = vstride;
        emxEnsureCapacity((emxArray__common *)c_B, i9, (int)sizeof(double));
        for (i9 = 0; i9 < vstride; i9++) {
          c_B->data[i9] = B->data[i9 + B->size[0] * nx] / pp;
        }

        loop_ub = c_B->size[0];
        for (i9 = 0; i9 < loop_ub; i9++) {
          B->data[i9 + B->size[0] * nx] = c_B->data[i9];
        }
      }

      loop_ub = B->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = B->size[0];
        for (i9 = 0; i9 < br; i9++) {
          BB->data[(i9 + BB->size[0] * i8) + BB->size[0] * BB->size[1] * nm1d2] =
            B->data[i9 + B->size[0] * i8];
        }
      }

      nm1d2++;
    }

    emxFree_real_T(&c_B);
    emxFree_real_T(&b_B);
    emxFree_real_T(&B);
    i8 = a->size[0] * a->size[1];
    a->size[0] = x->size[0];
    a->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)a, i8, (int)sizeof(double));
    loop_ub = x->size[0] * x->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      a->data[i8] = x->data[i8];
    }

    inv(ss, b);
    emxInit_real_T(&d_x, 2);
    if ((x->size[1] == 1) || (b->size[0] == 1)) {
      i8 = d_x->size[0] * d_x->size[1];
      d_x->size[0] = x->size[0];
      d_x->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)d_x, i8, (int)sizeof(double));
      loop_ub = x->size[0];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = b->size[1];
        for (i9 = 0; i9 < br; i9++) {
          d_x->data[i8 + d_x->size[0] * i9] = 0.0;
          vstride = x->size[1];
          for (k1 = 0; k1 < vstride; k1++) {
            d_x->data[i8 + d_x->size[0] * i9] += x->data[i8 + x->size[0] * k1] *
              b->data[k1 + b->size[0] * i9];
          }
        }
      }

      i8 = x->size[0] * x->size[1];
      x->size[0] = d_x->size[0];
      x->size[1] = d_x->size[1];
      emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
      loop_ub = d_x->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = d_x->size[0];
        for (i9 = 0; i9 < br; i9++) {
          x->data[i9 + x->size[0] * i8] = d_x->data[i9 + d_x->size[0] * i8];
        }
      }
    } else {
      k = x->size[1];
      nx = x->size[0];
      vstride = b->size[1];
      i8 = x->size[0] * x->size[1];
      x->size[0] = nx;
      x->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
      m = a->size[0];
      i8 = x->size[0] * x->size[1];
      emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
      loop_ub = x->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        br = x->size[0];
        for (i9 = 0; i9 < br; i9++) {
          x->data[i9 + x->size[0] * i8] = 0.0;
        }
      }

      if ((a->size[0] == 0) || (b->size[1] == 0)) {
      } else {
        vstride = a->size[0] * (b->size[1] - 1);
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          i8 = nx + m;
          for (ic = nx; ic + 1 <= i8; ic++) {
            x->data[ic] = 0.0;
          }

          nx += m;
        }

        br = 0;
        nx = 0;
        while ((m > 0) && (nx <= vstride)) {
          ar = 0;
          i8 = br + k;
          for (ib = br; ib + 1 <= i8; ib++) {
            if (b->data[ib] != 0.0) {
              ia = ar;
              i9 = nx + m;
              for (ic = nx; ic + 1 <= i9; ic++) {
                ia++;
                x->data[ic] += b->data[ib] * a->data[ia - 1];
              }
            }

            ar += m;
          }

          br += k;
          nx += m;
        }
      }
    }

    emxFree_real_T(&d_x);

    /*  for ip = 1:p */
    /*       cv(ip) = CVm(x*BB(:,1:ip, ip), ky); */
    /*  end */
    nx = BB->size[0] * BB->size[1] * BB->size[2];
    i8 = BB1D->size[0];
    BB1D->size[0] = (int)rt_powd_snf((double)p, 3.0);
    emxEnsureCapacity((emxArray__common *)BB1D, i8, (int)sizeof(double));
    for (k = 0; k + 1 <= nx; k++) {
      BB1D->data[k] = BB->data[k];
    }

    emxFree_real_T(&BB);
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&qx);
  emxFree_creal_T(&Dc);
  emxFree_creal_T(&Vc);
  emxFree_real_T(&dc);
  emxFree_real_T(&s);
  emxFree_real_T(&U);
  emxFree_real_T(&C);
  emxFree_real_T(&yj);
  emxFree_real_T(&yi);
  emxFree_real_T(&DD);
  emxFree_real_T(&ky2);
  emxFree_real_T(&ky1);
  emxFree_real_T(&V);
  emxFree_real_T(&ss);
}

/* End of code generation (MAVEfast.cpp) */
