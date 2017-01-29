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
#include "CVfast.h"
#include "MAVEfast.h"
#include "CVfast_emxutil.h"
#include "power.h"
#include "inv.h"
#include "eye.h"
#include "repmat.h"
#include "exp.h"
#include "sort1.h"
#include "diag.h"
#include "eig.h"
#include "kron.h"
#include "sum.h"
#include "strcmp.h"
#include "mean.h"
#include "std.h"
#include "rdivide.h"
#include "sqrt.h"
#include "abs.h"
#include "quantile.h"
#include "CVfast_rtwutil.h"
#include <stdio.h>
/* Function Declarations */
static void unifD(const emxArray_real_T *x, double m, emxArray_real_T *I);

/* Function Definitions */
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
  int I_count;
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
    I->size[0] = (int)m;
    I->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)I, b_k, (int)sizeof(double));
    loop_ub = (int)m;
    for (b_k = 0; b_k < loop_ub; b_k++) {
      I->data[b_k] = 0.0;
    }

    /* Jan 13 */
    I_count = -1;

    /* Jan 13 */
    emxInit_int32_T(&Ii, 1);
    emxInit_int32_T(&b_Ii, 1);
    emxInit_int32_T(&iidx, 1);
    emxInit_real_T(&b_x, 2);
    emxInit_boolean_T(&b, 2);
    emxInit_real_T(&b_Xremain, 2);
    emxInit_real_T(&b_y, 2);
    emxInit_real_T(&c_Xremain, 2);
    emxInit_real_T(&d_Xremain, 2);
    emxInit_real_T(&e_Xremain, 2);
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
      sum(b_x, varargin_2);
      sort(varargin_2, iidx);
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
      sum(b_x, varargin_2);
      sort(varargin_2, iidx);
      b_k = varargin_2->size[0];
      varargin_2->size[0] = iidx->size[0];
      emxEnsureCapacity((emxArray__common *)varargin_2, b_k, (int)sizeof(double));
      loop_ub = iidx->size[0];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        varargin_2->data[b_k] = iidx->data[b_k];
      }

      /* I = [I, Xremain(Ii(J(1)),p+1)]; */
      I_count++;

      /* Jan 13 */
      I->data[I_count] = Xremain->data[(b_Ii->data[(int)varargin_2->data[0] - 1]
        + Xremain->size[0] * p) - 1];

      /* Jan 13 */
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

void MAVEfast(emxArray_real_T *x, const emxArray_real_T *y, const
              emxArray_char_T *method, emxArray_real_T *BB1D, emxArray_real_T
              *ky)
{
  emxArray_real_T *qx;
  emxArray_real_T *a;
  int p;
  int n;
  int i15;
  int loop_ub;
  int br;
  int i16;
  emxArray_real_T *ss;
  int k;
  int vstride;
  int nx;
  int m;
  emxArray_real_T *b_ss;
  double pp;
  int i17;
  emxArray_real_T *V;
  int ic;
  emxArray_creal_T *Vc;
  emxArray_creal_T *Dc;
  int ar;
  int ib;
  int ia;
  emxArray_real_T *b;
  emxArray_real_T *yj;
  emxArray_real_T *b_x;
  boolean_T empty_non_axis_sizes;
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
  emxArray_boolean_T *c_ss;
  emxArray_real_T *c_qx;
  emxArray_real_T *b_ky;
  emxArray_real_T *b_V;
  emxArray_real_T *c_V;
  boolean_T guard2 = false;
  int b_m;
  emxArray_real_T *BB;
  emxArray_real_T *r1;
  emxArray_real_T *onexi;
  emxArray_real_T *B;
  int ip;
  emxArray_real_T *D;
  emxArray_real_T *Ifast;
  emxArray_real_T *Vfast;
  emxArray_real_T *xfast;
  emxArray_real_T *xij;
  emxArray_real_T *dxij;
  emxArray_real_T *abi;
  emxArray_real_T *dd;
  emxArray_real_T *kxijy;
  emxArray_real_T *ddx;
  emxArray_real_T *tmp;
  emxArray_real_T *b_B;
  emxArray_creal_T *R;
  emxArray_real_T *b_C;
  emxArray_real_T *c_C;
  emxArray_real_T *d_C;
  emxArray_real_T *e_C;
  emxArray_int32_T *r2;
  emxArray_real_T *b_y;
  emxArray_real_T *c_y;
  emxArray_real_T *d_y;
  cell_wrap_0 reshapes[2];
  emxArray_real_T *e_y;
  emxArray_real_T *f_y;
  emxArray_real_T *g_y;
  emxArray_real_T *h_y;
  emxArray_real_T *b_dd;
  emxArray_real_T *b_abi;
  emxArray_real_T *f_C;
  emxArray_real_T *d_V;
  emxArray_real_T *b_Vfast;
  emxArray_real_T *b_xfast;
  emxArray_real_T *g_C;
  emxArray_real_T *e_V;
  emxArray_real_T *c_B;
  emxArray_int32_T *b_Ifast;
  emxArray_int32_T *c_Ifast;
  emxArray_real_T *b_dc;
  emxArray_real_T *c_abi;
  emxArray_int32_T *d_Ifast;
  emxArray_int32_T *e_Ifast;
  boolean_T exitg1;
  int b_ip;
  int K;
  int iter;
  double c_m;
  int Ifast_idx_0;
  int j;
  double h2;
  boolean_T guard1 = false;
  boolean_T b0;
  boolean_T exitg2;
  boolean_T b_guard1 = false;
  unsigned int a_idx_0;
  int k1;
  double b_a;
  double d;
  double ndbl;
  double apnd;
  double cdiff;
  emxInit_real_T(&qx, 2);
  emxInit_real_T(&a, 2);

  /*  */
  /*  [estB,cv] = MAVE(x, y, method, CV); */
  /*  INPUT */
  /*  x: nxp matrix of predictors */
  /*  y: nx1 of response; */
  /*  method: 'meanOPG'; 'meanMAVE'; 'csMAVE'; 'csOPG'; 'kSIR' */
  /*          by default, method = 'csOPG' */
  /*          (1) 'meanOPG' and 'meanMAVE' estimate dimension reduction space */
  /*              for conditional mean */
  /*          (2) 'csMAVE' and 'csOPG' estimate the central dimension reduction */
  /*              space */
  /*          (3) 'kSIR' is a kernel version of SIR (Li, 1991). It is fast, but */
  /*              with poor efficiency */
  /*  CV:  'yes'--to calculate the CV values, by default */
  /*       'no'--dont calculate CV values to make the calculation faster */
  /*  OUTPUT */
  /*  estB: pxpxp matrix */
  /*      estB(:,1:d, d) is the central space with dimension=d */
  /*          d = 1, 2, ..., p */
  /*  cv: Cross-Validation (CV) values for d = 1, 2, ..., p. */
  /*      the dimension with smallest CV is the estimated dimension. */
  /*  */
  /*  Reference,  Xia et al (2002, JRSSB); Xia (2007, AoS), Wang and Xia (2008, JASA) */
  /*  */
  /*  %Example */
  /*   n = 200; p = 10; */
  /*   B = zeros(p,2); */
  /*   B(1,1)=1; */
  /*   B(2,2)=1; */
  /*   x = randn(n,p); */
  /*   y = (x*B(:,1)) + (x*B(:,2)).*randn(n,1); */
  /*   [estB, cv] = MAVE(x,y,'csMAVE'); */
  /*   %Estimated directions */
  /*     B1 = estB(:,1:1,1)  % when dimension is fixed as 1 */
  /*     B2 = estB(:,1:2,2)  % when dimension is fixed as 2 */
  /*     B3 = estB(:,1:3,3)  % when dimension is fixed as 3 */
  /*     ... */
  /*   % Estimated dimension */
  /*      find(cv ==min(cv)) */
  p = x->size[1];
  n = x->size[0];
  //printf("n=%d p=%d\n",n,p);printf("x=\n");for(int i=0;i<x->size[0];++i){for(int j=0;j<x->size[1];++j) printf("%.3lf ",x->data[j*x->size[0]+i]);printf("\n");}
  b_mean(x, qx);
  repmat(qx, (double)x->size[0], a);
  i15 = x->size[0] * x->size[1];
  emxEnsureCapacity((emxArray__common *)x, i15, (int)sizeof(double));
  loop_ub = x->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    br = x->size[0];
    for (i16 = 0; i16 < br; i16++) {
      x->data[i16 + x->size[0] * i15] -= a->data[i16 + a->size[0] * i15];
    }
  }
  //printf("x=\n");for(int i=0;i<x->size[0];++i){for(int j=0;j<x->size[1];++j) printf("%.3lf ",x->data[j*x->size[0]+i]);printf("\n");}
  i15 = a->size[0] * a->size[1];
  a->size[0] = x->size[1];//as
  a->size[1] = x->size[0];
  emxEnsureCapacity((emxArray__common *)a, i15, (int)sizeof(double));
  loop_ub = x->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    br = x->size[1];
    for (i16 = 0; i16 < br; i16++) {
      a->data[i16 + a->size[0] * i15] = x->data[i15 + x->size[0] * i16];
    }
  }

  emxInit_real_T(&ss, 2);
  if ((a->size[1] == 1) || (x->size[0] == 1)) {
    i15 = ss->size[0] * ss->size[1];
    ss->size[0] = a->size[0];
    ss->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i15, (int)sizeof(double));
    loop_ub = a->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = x->size[1];
      for (i16 = 0; i16 < br; i16++) {
        ss->data[i15 + ss->size[0] * i16] = 0.0;
        nx = a->size[1];
        for (i17 = 0; i17 < nx; i17++) {
          ss->data[i15 + ss->size[0] * i16] += a->data[i15 + a->size[0] * i17] *
            x->data[i17 + x->size[0] * i16];
        }
      }
    }
  } else {
    k = a->size[1];
    vstride = a->size[0];
    nx = x->size[1];
    i15 = ss->size[0] * ss->size[1];
    ss->size[0] = vstride;
    ss->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)ss, i15, (int)sizeof(double));
    m = a->size[0];
    i15 = ss->size[0] * ss->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i15, (int)sizeof(double));
    loop_ub = ss->size[1];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = ss->size[0];
      for (i16 = 0; i16 < br; i16++) {
        ss->data[i16 + ss->size[0] * i15] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (x->size[1] == 0)) {
    } else {
      vstride = a->size[0] * (x->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i15 = nx + m;
        for (ic = nx; ic + 1 <= i15; ic++) {
          ss->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i15 = br + k;
        for (ib = br; ib + 1 <= i15; ib++) {
          if (x->data[ib] != 0.0) {
            ia = ar;
            i16 = nx + m;
            for (ic = nx; ic + 1 <= i16; ic++) {
              ia++;
              ss->data[ic] += x->data[ib] * a->data[ia - 1];
            }
          }

          ar += m;
        }

        br += k;
        nx += m;
      }
    }
  }

  emxInit_real_T(&b_ss, 2);
  pp = rt_powd_snf((double)n, 3.0);
  eye((double)p, a);
  i15 = b_ss->size[0] * b_ss->size[1];
  b_ss->size[0] = ss->size[0];
  b_ss->size[1] = ss->size[1];
  emxEnsureCapacity((emxArray__common *)b_ss, i15, (int)sizeof(double));
  loop_ub = ss->size[0] * ss->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    b_ss->data[i15] = ss->data[i15] / (double)n + a->data[i15] / pp;
  }

  emxInit_real_T(&V, 2);
  emxInit_creal_T1(&Vc, 2);
  emxInit_creal_T1(&Dc, 2);
  inv(b_ss, ss);
  eig(ss, Vc, Dc);
  i15 = V->size[0] * V->size[1];
  V->size[0] = Vc->size[0];
  V->size[1] = Vc->size[1];
  emxEnsureCapacity((emxArray__common *)V, i15, (int)sizeof(double));
  loop_ub = Vc->size[0] * Vc->size[1];
  emxFree_real_T(&b_ss);
  for (i15 = 0; i15 < loop_ub; i15++) {
    V->data[i15] = Vc->data[i15].re;
  }

  emxInit_real_T(&b, 2);
  i15 = b->size[0] * b->size[1];
  b->size[0] = Dc->size[0];
  b->size[1] = Dc->size[1];
  emxEnsureCapacity((emxArray__common *)b, i15, (int)sizeof(double));
  loop_ub = Dc->size[0] * Dc->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    b->data[i15] = Dc->data[i15].re;
  }

  c_sqrt(b);
  emxInit_real_T(&yj, 2);
  if ((V->size[1] == 1) || (b->size[0] == 1)) {
    i15 = yj->size[0] * yj->size[1];
    yj->size[0] = V->size[0];
    yj->size[1] = b->size[1];
    emxEnsureCapacity((emxArray__common *)yj, i15, (int)sizeof(double));
    loop_ub = V->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = b->size[1];
      for (i16 = 0; i16 < br; i16++) {
        yj->data[i15 + yj->size[0] * i16] = 0.0;
        nx = V->size[1];
        for (i17 = 0; i17 < nx; i17++) {
          yj->data[i15 + yj->size[0] * i16] += V->data[i15 + V->size[0] * i17] *
            b->data[i17 + b->size[0] * i16];
        }
      }
    }
  } else {
    k = V->size[1];
    vstride = V->size[0];
    nx = b->size[1];
    i15 = yj->size[0] * yj->size[1];
    yj->size[0] = vstride;
    yj->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)yj, i15, (int)sizeof(double));
    m = V->size[0];
    i15 = yj->size[0] * yj->size[1];
    emxEnsureCapacity((emxArray__common *)yj, i15, (int)sizeof(double));
    loop_ub = yj->size[1];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = yj->size[0];
      for (i16 = 0; i16 < br; i16++) {
        yj->data[i16 + yj->size[0] * i15] = 0.0;
      }
    }

    if ((V->size[0] == 0) || (b->size[1] == 0)) {
    } else {
      vstride = V->size[0] * (b->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i15 = nx + m;
        for (ic = nx; ic + 1 <= i15; ic++) {
          yj->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i15 = br + k;
        for (ib = br; ib + 1 <= i15; ib++) {
          if (b->data[ib] != 0.0) {
            ia = ar;
            i16 = nx + m;
            for (ic = nx; ic + 1 <= i16; ic++) {
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

  i15 = b->size[0] * b->size[1];
  b->size[0] = V->size[1];
  b->size[1] = V->size[0];
  emxEnsureCapacity((emxArray__common *)b, i15, (int)sizeof(double));
  loop_ub = V->size[0];
  for (i15 = 0; i15 < loop_ub; i15++) {
    br = V->size[1];
    for (i16 = 0; i16 < br; i16++) {
      b->data[i16 + b->size[0] * i15] = V->data[i15 + V->size[0] * i16];
    }
  }

  if ((yj->size[1] == 1) || (b->size[0] == 1)) {
    i15 = ss->size[0] * ss->size[1];
    ss->size[0] = yj->size[0];
    ss->size[1] = b->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i15, (int)sizeof(double));
    loop_ub = yj->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = b->size[1];
      for (i16 = 0; i16 < br; i16++) {
        ss->data[i15 + ss->size[0] * i16] = 0.0;
        nx = yj->size[1];
        for (i17 = 0; i17 < nx; i17++) {
          ss->data[i15 + ss->size[0] * i16] += yj->data[i15 + yj->size[0] * i17]
            * b->data[i17 + b->size[0] * i16];
        }
      }
    }
  } else {
    k = yj->size[1];
    vstride = yj->size[0];
    nx = b->size[1];
    i15 = ss->size[0] * ss->size[1];
    ss->size[0] = vstride;
    ss->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)ss, i15, (int)sizeof(double));
    m = yj->size[0];
    i15 = ss->size[0] * ss->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i15, (int)sizeof(double));
    loop_ub = ss->size[1];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = ss->size[0];
      for (i16 = 0; i16 < br; i16++) {
        ss->data[i16 + ss->size[0] * i15] = 0.0;
      }
    }

    if ((yj->size[0] == 0) || (b->size[1] == 0)) {
    } else {
      vstride = yj->size[0] * (b->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i15 = nx + m;
        for (ic = nx; ic + 1 <= i15; ic++) {
          ss->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i15 = br + k;
        for (ib = br; ib + 1 <= i15; ib++) {
          if (b->data[ib] != 0.0) {
            ia = ar;
            i16 = nx + m;
            for (ic = nx; ic + 1 <= i16; ic++) {
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

  i15 = a->size[0] * a->size[1];
  a->size[0] = x->size[0];
  a->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)a, i15, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    a->data[i15] = x->data[i15];
  }

  emxInit_real_T(&b_x, 2);
  if ((x->size[1] == 1) || (ss->size[0] == 1)) {
    i15 = b_x->size[0] * b_x->size[1];
    b_x->size[0] = x->size[0];
    b_x->size[1] = ss->size[1];
    emxEnsureCapacity((emxArray__common *)b_x, i15, (int)sizeof(double));
    loop_ub = x->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = ss->size[1];
      for (i16 = 0; i16 < br; i16++) {
        b_x->data[i15 + b_x->size[0] * i16] = 0.0;
        nx = x->size[1];
        for (i17 = 0; i17 < nx; i17++) {
          b_x->data[i15 + b_x->size[0] * i16] += x->data[i15 + x->size[0] * i17]
            * ss->data[i17 + ss->size[0] * i16];
        }
      }
    }

    i15 = x->size[0] * x->size[1];
    x->size[0] = b_x->size[0];
    x->size[1] = b_x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i15, (int)sizeof(double));
    loop_ub = b_x->size[1];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = b_x->size[0];
      for (i16 = 0; i16 < br; i16++) {
        x->data[i16 + x->size[0] * i15] = b_x->data[i16 + b_x->size[0] * i15];
      }
    }
  } else {
    k = x->size[1];
    vstride = x->size[0];
    nx = ss->size[1];
    i15 = x->size[0] * x->size[1];
    x->size[0] = vstride;
    x->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)x, i15, (int)sizeof(double));
    m = a->size[0];
    i15 = x->size[0] * x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i15, (int)sizeof(double));
    loop_ub = x->size[1];
    for (i15 = 0; i15 < loop_ub; i15++) {
      br = x->size[0];
      for (i16 = 0; i16 < br; i16++) {
        x->data[i16 + x->size[0] * i15] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (ss->size[1] == 0)) {
    } else {
      vstride = a->size[0] * (ss->size[1] - 1);
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        i15 = nx + m;
        for (ic = nx; ic + 1 <= i15; ic++) {
          x->data[ic] = 0.0;
        }

        nx += m;
      }

      br = 0;
      nx = 0;
      while ((m > 0) && (nx <= vstride)) {
        ar = 0;
        i15 = br + k;
        for (ib = br; ib + 1 <= i15; ib++) {
          if (ss->data[ib] != 0.0) {
            ia = ar;
            i16 = nx + m;
            for (ic = nx; ic + 1 <= i16; ic++) {
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
  empty_non_axis_sizes = b_strcmp(method);
  emxInit_real_T(&ky1, 2);
  emxInit_real_T(&ky2, 2);
  emxInit_real_T(&DD, 2);
  emxInit_real_T(&b_qx, 2);
  emxInit_real_T(&C, 2);
  emxInit_real_T(&U, 2);
  emxInit_real_T1(&s, 1);
  emxInit_real_T1(&dc, 1);
  emxInit_int32_T(&iidx, 1);
  emxInit_boolean_T1(&c_x, 1);
  emxInit_int32_T(&ii, 1);
  emxInit_real_T(&b_s, 2);
  emxInit_real_T(&b_Dc, 2);
  emxInit_boolean_T(&c_ss, 2);
  emxInit_real_T(&c_qx, 2);
  emxInit_real_T(&b_ky, 2);
  emxInit_real_T(&b_V, 2);
  emxInit_real_T(&c_V, 2);
  guard2 = false;
  if (empty_non_axis_sizes) {
    guard2 = true;
  } else {
    empty_non_axis_sizes = c_strcmp(method);
    if (empty_non_axis_sizes) {
      guard2 = true;
    } else {
      if (99 <= n) {
        b_m = 99;
      } else {
        b_m = n;
      }

      pp = std::floor(rt_powd_snf((double)n, 0.6));
      if ((b_m >= pp) || rtIsNaN(pp)) {
        pp = b_m;
      }

      /*  define 100 quantiles */
      if (rtIsNaN(pp)) {
        i15 = qx->size[0] * qx->size[1];
        qx->size[0] = 1;
        qx->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)qx, i15, (int)sizeof(double));
        qx->data[0] = rtNaN;
      } else if (pp < 1.0) {
        i15 = qx->size[0] * qx->size[1];
        qx->size[0] = 1;
        qx->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)qx, i15, (int)sizeof(double));
      } else if (rtIsInf(pp) && (1.0 == pp)) {
        i15 = qx->size[0] * qx->size[1];
        qx->size[0] = 1;
        qx->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)qx, i15, (int)sizeof(double));
        qx->data[0] = rtNaN;
      } else {
        i15 = qx->size[0] * qx->size[1];
        qx->size[0] = 1;
        qx->size[1] = (int)(pp - 1.0) + 1;
        emxEnsureCapacity((emxArray__common *)qx, i15, (int)sizeof(double));
        loop_ub = (int)(pp - 1.0);
        for (i15 = 0; i15 <= loop_ub; i15++) {
          qx->data[qx->size[0] * i15] = 1.0 + (double)i15;
        }
      }

      i15 = c_qx->size[0] * c_qx->size[1];
      c_qx->size[0] = 1;
      c_qx->size[1] = qx->size[1];
      emxEnsureCapacity((emxArray__common *)c_qx, i15, (int)sizeof(double));
      loop_ub = qx->size[0] * qx->size[1];
      for (i15 = 0; i15 < loop_ub; i15++) {
        c_qx->data[i15] = qx->data[i15] / (pp + 1.0);
      }

      quantile(y, c_qx, b_qx);
      if ((b_qx->size[0] == 0) || (b_qx->size[1] == 0)) {
        b_m = 0;
      } else {
        nx = b_qx->size[0];
        b_m = b_qx->size[1];
        if (nx >= b_m) {
          b_m = nx;
        }
      }

      nx = b_qx->size[0] * b_qx->size[1];
      i15 = qx->size[0] * qx->size[1];
      qx->size[0] = 1;
      qx->size[1] = b_m;
      emxEnsureCapacity((emxArray__common *)qx, i15, (int)sizeof(double));
      for (k = 0; k + 1 <= nx; k++) {
        qx->data[k] = b_qx->data[k];
      }

      /* qx = unique(qx); %Jan 13 Change for C */
      c_repmat(y, (double)qx->size[1], ss);
      repmat(qx, (double)n, yj);
      i15 = c_ss->size[0] * c_ss->size[1];
      c_ss->size[0] = ss->size[0];
      c_ss->size[1] = ss->size[1];
      emxEnsureCapacity((emxArray__common *)c_ss, i15, (int)sizeof(boolean_T));
      loop_ub = ss->size[0] * ss->size[1];
      for (i15 = 0; i15 < loop_ub; i15++) {
        c_ss->data[i15] = (ss->data[i15] - yj->data[i15] < 0.0);
      }

      b_abs(c_ss, ky);
      b_mean(ky, qx);
      repmat(qx, (double)n, a);
      i15 = ky->size[0] * ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i15, (int)sizeof(double));
      b_m = ky->size[0];
      vstride = ky->size[1];
      loop_ub = b_m * vstride;
      for (i15 = 0; i15 < loop_ub; i15++) {
        ky->data[i15] -= a->data[i15];
      }

      c_mean(ky, s);
      b_repmat(s, (double)ky->size[1], ky1);
      i15 = ky1->size[0] * ky1->size[1];
      ky1->size[0] = ky->size[0];
      ky1->size[1] = ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky1, i15, (int)sizeof(double));
      loop_ub = ky->size[0] * ky->size[1];
      for (i15 = 0; i15 < loop_ub; i15++) {
        ky1->data[i15] = ky->data[i15] - ky1->data[i15];
      }

      /*  LATEST CHANGED */
      b_mean(ky1, qx);
      repmat(qx, (double)n, a);
      i15 = ky1->size[0] * ky1->size[1];
      emxEnsureCapacity((emxArray__common *)ky1, i15, (int)sizeof(double));
      b_m = ky1->size[0];
      vstride = ky1->size[1];
      loop_ub = b_m * vstride;
      for (i15 = 0; i15 < loop_ub; i15++) {
        ky1->data[i15] -= a->data[i15];
      }

      /*  LATEST CHANGED */
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
      if (empty_non_axis_sizes || (!((ky->size[0] == 0) || (ky->size[1] == 0))))
      {
        b_m = ky->size[1];
      } else {
        b_m = 0;
      }

      if (empty_non_axis_sizes || (!((ky1->size[0] == 0) || (ky1->size[1] == 0))))
      {
        vstride = ky1->size[1];
      } else {
        vstride = 0;
      }

      i15 = b_ky->size[0] * b_ky->size[1];
      b_ky->size[0] = nx;
      b_ky->size[1] = b_m + vstride;
      emxEnsureCapacity((emxArray__common *)b_ky, i15, (int)sizeof(double));
      for (i15 = 0; i15 < b_m; i15++) {
        for (i16 = 0; i16 < nx; i16++) {
          b_ky->data[i16 + b_ky->size[0] * i15] = ky->data[i16 + nx * i15];
        }
      }

      for (i15 = 0; i15 < vstride; i15++) {
        for (i16 = 0; i16 < nx; i16++) {
          b_ky->data[i16 + b_ky->size[0] * (i15 + b_m)] = ky1->data[i16 + nx *
            i15];
        }
      }

      i15 = ky->size[0] * ky->size[1];
      ky->size[0] = b_ky->size[0];
      ky->size[1] = b_ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i15, (int)sizeof(double));
      loop_ub = b_ky->size[1];
      for (i15 = 0; i15 < loop_ub; i15++) {
        br = b_ky->size[0];
        for (i16 = 0; i16 < br; i16++) {
          ky->data[i16 + ky->size[0] * i15] = b_ky->data[i16 + b_ky->size[0] *
            i15];
        }
      }

      /*  LATEST CHANGED */
      /* clear yi ky1; %Jan 13 not supported */
      /* yi = []; ky1 = []; %%% latest changes by XIA */
      i15 = a->size[0] * a->size[1];
      a->size[0] = ky->size[1];
      a->size[1] = ky->size[0];
      emxEnsureCapacity((emxArray__common *)a, i15, (int)sizeof(double));
      loop_ub = ky->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        br = ky->size[1];
        for (i16 = 0; i16 < br; i16++) {
          a->data[i16 + a->size[0] * i15] = ky->data[i15 + ky->size[0] * i16];
        }
      }

      if ((a->size[1] == 1) || (ky->size[0] == 1)) {
        i15 = C->size[0] * C->size[1];
        C->size[0] = a->size[0];
        C->size[1] = ky->size[1];
        emxEnsureCapacity((emxArray__common *)C, i15, (int)sizeof(double));
        loop_ub = a->size[0];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = ky->size[1];
          for (i16 = 0; i16 < br; i16++) {
            C->data[i15 + C->size[0] * i16] = 0.0;
            nx = a->size[1];
            for (i17 = 0; i17 < nx; i17++) {
              C->data[i15 + C->size[0] * i16] += a->data[i15 + a->size[0] * i17]
                * ky->data[i17 + ky->size[0] * i16];
            }
          }
        }
      } else {
        k = a->size[1];
        vstride = a->size[0];
        nx = ky->size[1];
        i15 = C->size[0] * C->size[1];
        C->size[0] = vstride;
        C->size[1] = nx;
        emxEnsureCapacity((emxArray__common *)C, i15, (int)sizeof(double));
        m = a->size[0];
        i15 = C->size[0] * C->size[1];
        emxEnsureCapacity((emxArray__common *)C, i15, (int)sizeof(double));
        loop_ub = C->size[1];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = C->size[0];
          for (i16 = 0; i16 < br; i16++) {
            C->data[i16 + C->size[0] * i15] = 0.0;
          }
        }

        if ((a->size[0] == 0) || (ky->size[1] == 0)) {
        } else {
          vstride = a->size[0] * (ky->size[1] - 1);
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            i15 = nx + m;
            for (ic = nx; ic + 1 <= i15; ic++) {
              C->data[ic] = 0.0;
            }

            nx += m;
          }

          br = 0;
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            ar = 0;
            i15 = br + k;
            for (ib = br; ib + 1 <= i15; ib++) {
              if (ky->data[ib] != 0.0) {
                ia = ar;
                i16 = nx + m;
                for (ic = nx; ic + 1 <= i16; ic++) {
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

      /* Jan 7 */
      eig(C, Vc, Dc);

      /* Jan 7 */
      i15 = V->size[0] * V->size[1];
      V->size[0] = Vc->size[0];
      V->size[1] = Vc->size[1];
      emxEnsureCapacity((emxArray__common *)V, i15, (int)sizeof(double));
      loop_ub = Vc->size[0] * Vc->size[1];
      for (i15 = 0; i15 < loop_ub; i15++) {
        V->data[i15] = Vc->data[i15].re;
      }

      /* Jan 13 */
      /* clear C; %Jan 13 not supported %Jan 7 */
      i15 = b_Dc->size[0] * b_Dc->size[1];
      b_Dc->size[0] = Dc->size[0];
      b_Dc->size[1] = Dc->size[1];
      emxEnsureCapacity((emxArray__common *)b_Dc, i15, (int)sizeof(double));
      loop_ub = Dc->size[0] * Dc->size[1];
      for (i15 = 0; i15 < loop_ub; i15++) {
        b_Dc->data[i15] = Dc->data[i15].re;
      }

      diag(b_Dc, s);
      c_abs(s, dc);
      d_sort(dc, iidx);

      /* Jan 7 */
      b_m = 2;
      if (dc->size[0] != 1) {
        b_m = 1;
      }

      i15 = s->size[0];
      s->size[0] = dc->size[0];
      emxEnsureCapacity((emxArray__common *)s, i15, (int)sizeof(double));
      loop_ub = dc->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        s->data[i15] = dc->data[i15];
      }

      if (b_m <= 1) {
        i15 = dc->size[0];
      } else {
        i15 = 1;
      }

      if ((!(dc->size[0] == 0)) && (i15 > 1)) {
        vstride = 1;
        k = 1;
        while (k <= b_m - 1) {
          vstride *= dc->size[0];
          k = 2;
        }

        for (j = 0; j + 1 <= vstride; j++) {
          for (k = 1; k < i15; k++) {
            s->data[j + k * vstride] += s->data[j + (k - 1) * vstride];
          }
        }
      }

      pp = b_sum(dc);
      i15 = c_x->size[0];
      c_x->size[0] = s->size[0];
      emxEnsureCapacity((emxArray__common *)c_x, i15, (int)sizeof(boolean_T));
      loop_ub = s->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        c_x->data[i15] = (s->data[i15] / pp >= 0.99);
      }

      nx = c_x->size[0];
      vstride = 0;
      i15 = ii->size[0];
      ii->size[0] = c_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      b_m = 1;
      exitg2 = false;
      while ((!exitg2) && (b_m <= nx)) {
        b_guard1 = false;
        if (c_x->data[b_m - 1]) {
          vstride++;
          ii->data[vstride - 1] = b_m;
          if (vstride >= nx) {
            exitg2 = true;
          } else {
            b_guard1 = true;
          }
        } else {
          b_guard1 = true;
        }

        if (b_guard1) {
          b_m++;
        }
      }

      if (c_x->size[0] == 1) {
        if (vstride == 0) {
          i15 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
        }
      } else {
        i15 = ii->size[0];
        if (1 > vstride) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = vstride;
        }

        emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      }

      i15 = s->size[0];
      s->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)s, i15, (int)sizeof(double));
      loop_ub = ii->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        s->data[i15] = ii->data[i15];
      }

      nx = s->size[0];
      vstride = (int)s->data[0];
      if (s->size[0] > 1) {
        for (b_m = 1; b_m + 1 <= nx; b_m++) {
          if ((int)s->data[b_m] < vstride) {
            vstride = (int)s->data[b_m];
          }
        }
      }

      /* Jan 7 */
      b_m = V->size[0];
      i15 = b_V->size[0] * b_V->size[1];
      b_V->size[0] = b_m;
      b_V->size[1] = iidx->size[0];
      emxEnsureCapacity((emxArray__common *)b_V, i15, (int)sizeof(double));
      loop_ub = iidx->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        for (i16 = 0; i16 < b_m; i16++) {
          b_V->data[i16 + b_V->size[0] * i15] = V->data[i16 + V->size[0] *
            (iidx->data[i15] - 1)];
        }
      }

      i15 = V->size[0] * V->size[1];
      V->size[0] = b_V->size[0];
      V->size[1] = b_V->size[1];
      emxEnsureCapacity((emxArray__common *)V, i15, (int)sizeof(double));
      loop_ub = b_V->size[1];
      for (i15 = 0; i15 < loop_ub; i15++) {
        br = b_V->size[0];
        for (i16 = 0; i16 < br; i16++) {
          V->data[i16 + V->size[0] * i15] = b_V->data[i16 + b_V->size[0] * i15];
        }
      }

      /* Jan 7 */
      if (1 > vstride) {
        loop_ub = 0;
      } else {
        loop_ub = vstride;
      }

      /* Jan 7 */
      if (1 > vstride) {
        br = 0;
      } else {
        br = vstride;
      }

      b_m = V->size[0];
      i15 = c_V->size[0] * c_V->size[1];
      c_V->size[0] = b_m;
      c_V->size[1] = br;
      emxEnsureCapacity((emxArray__common *)c_V, i15, (int)sizeof(double));
      for (i15 = 0; i15 < br; i15++) {
        for (i16 = 0; i16 < b_m; i16++) {
          c_V->data[i16 + c_V->size[0] * i15] = V->data[i16 + V->size[0] * i15];
        }
      }

      i15 = V->size[0] * V->size[1];
      V->size[0] = c_V->size[0];
      V->size[1] = c_V->size[1];
      emxEnsureCapacity((emxArray__common *)V, i15, (int)sizeof(double));
      br = c_V->size[1];
      for (i15 = 0; i15 < br; i15++) {
        nx = c_V->size[0];
        for (i16 = 0; i16 < nx; i16++) {
          V->data[i16 + V->size[0] * i15] = c_V->data[i16 + c_V->size[0] * i15];
        }
      }

      /* Jan 7 */
      if ((ky->size[1] == 1) || (V->size[0] == 1)) {
        i15 = U->size[0] * U->size[1];
        U->size[0] = ky->size[0];
        U->size[1] = V->size[1];
        emxEnsureCapacity((emxArray__common *)U, i15, (int)sizeof(double));
        br = ky->size[0];
        for (i15 = 0; i15 < br; i15++) {
          nx = V->size[1];
          for (i16 = 0; i16 < nx; i16++) {
            U->data[i15 + U->size[0] * i16] = 0.0;
            vstride = ky->size[1];
            for (i17 = 0; i17 < vstride; i17++) {
              U->data[i15 + U->size[0] * i16] += ky->data[i15 + ky->size[0] *
                i17] * V->data[i17 + V->size[0] * i16];
            }
          }
        }
      } else {
        k = ky->size[1];
        vstride = ky->size[0];
        nx = V->size[1];
        i15 = U->size[0] * U->size[1];
        U->size[0] = vstride;
        U->size[1] = nx;
        emxEnsureCapacity((emxArray__common *)U, i15, (int)sizeof(double));
        m = ky->size[0];
        i15 = U->size[0] * U->size[1];
        emxEnsureCapacity((emxArray__common *)U, i15, (int)sizeof(double));
        br = U->size[1];
        for (i15 = 0; i15 < br; i15++) {
          nx = U->size[0];
          for (i16 = 0; i16 < nx; i16++) {
            U->data[i16 + U->size[0] * i15] = 0.0;
          }
        }

        if ((ky->size[0] == 0) || (V->size[1] == 0)) {
        } else {
          vstride = ky->size[0] * (V->size[1] - 1);
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            i15 = nx + m;
            for (ic = nx; ic + 1 <= i15; ic++) {
              U->data[ic] = 0.0;
            }

            nx += m;
          }

          br = 0;
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            ar = 0;
            i15 = br + k;
            for (ib = br; ib + 1 <= i15; ib++) {
              if (V->data[ib] != 0.0) {
                ia = ar;
                i16 = nx + m;
                for (ic = nx; ic + 1 <= i16; ic++) {
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

      /* Jan 7 */
      i15 = s->size[0];
      s->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)s, i15, (int)sizeof(double));
      for (i15 = 0; i15 < loop_ub; i15++) {
        s->data[i15] = dc->data[i15];
      }

      d_sqrt(s);

      /* Jan 7 */
      /* ky1 = bsxfun(@(x,c)x./c, U, s'); %Jan 7 */
      i15 = b_s->size[0] * b_s->size[1];
      b_s->size[0] = 1;
      b_s->size[1] = s->size[0];
      emxEnsureCapacity((emxArray__common *)b_s, i15, (int)sizeof(double));
      loop_ub = s->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        b_s->data[b_s->size[0] * i15] = s->data[i15];
      }

      repmat(b_s, (double)U->size[0], a);
      rdivide(U, a, ky1);

      /* Jan 13 */
      b_diag(s, a);
      i15 = b->size[0] * b->size[1];
      b->size[0] = V->size[1];
      b->size[1] = V->size[0];
      emxEnsureCapacity((emxArray__common *)b, i15, (int)sizeof(double));
      loop_ub = V->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        br = V->size[1];
        for (i16 = 0; i16 < br; i16++) {
          b->data[i16 + b->size[0] * i15] = V->data[i15 + V->size[0] * i16];
        }
      }

      if ((a->size[1] == 1) || (b->size[0] == 1)) {
        i15 = ky2->size[0] * ky2->size[1];
        ky2->size[0] = a->size[0];
        ky2->size[1] = b->size[1];
        emxEnsureCapacity((emxArray__common *)ky2, i15, (int)sizeof(double));
        loop_ub = a->size[0];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = b->size[1];
          for (i16 = 0; i16 < br; i16++) {
            ky2->data[i15 + ky2->size[0] * i16] = 0.0;
            nx = a->size[1];
            for (i17 = 0; i17 < nx; i17++) {
              ky2->data[i15 + ky2->size[0] * i16] += a->data[i15 + a->size[0] *
                i17] * b->data[i17 + b->size[0] * i16];
            }
          }
        }
      } else {
        k = a->size[1];
        vstride = a->size[0];
        nx = b->size[1];
        i15 = ky2->size[0] * ky2->size[1];
        ky2->size[0] = vstride;
        ky2->size[1] = nx;
        emxEnsureCapacity((emxArray__common *)ky2, i15, (int)sizeof(double));
        m = a->size[0];
        i15 = ky2->size[0] * ky2->size[1];
        emxEnsureCapacity((emxArray__common *)ky2, i15, (int)sizeof(double));
        loop_ub = ky2->size[1];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = ky2->size[0];
          for (i16 = 0; i16 < br; i16++) {
            ky2->data[i16 + ky2->size[0] * i15] = 0.0;
          }
        }

        if ((a->size[0] == 0) || (b->size[1] == 0)) {
        } else {
          vstride = a->size[0] * (b->size[1] - 1);
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            i15 = nx + m;
            for (ic = nx; ic + 1 <= i15; ic++) {
              ky2->data[ic] = 0.0;
            }

            nx += m;
          }

          br = 0;
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            ar = 0;
            i15 = br + k;
            for (ib = br; ib + 1 <= i15; ib++) {
              if (b->data[ib] != 0.0) {
                ia = ar;
                i16 = nx + m;
                for (ic = nx; ic + 1 <= i16; ic++) {
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

      /* Jan 7 */
      if ((ky1->size[1] == 1) || (ky2->size[0] == 1)) {
        i15 = ky->size[0] * ky->size[1];
        ky->size[0] = ky1->size[0];
        ky->size[1] = ky2->size[1];
        emxEnsureCapacity((emxArray__common *)ky, i15, (int)sizeof(double));
        loop_ub = ky1->size[0];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = ky2->size[1];
          for (i16 = 0; i16 < br; i16++) {
            ky->data[i15 + ky->size[0] * i16] = 0.0;
            nx = ky1->size[1];
            for (i17 = 0; i17 < nx; i17++) {
              ky->data[i15 + ky->size[0] * i16] += ky1->data[i15 + ky1->size[0] *
                i17] * ky2->data[i17 + ky2->size[0] * i16];
            }
          }
        }
      } else {
        k = ky1->size[1];
        vstride = ky1->size[0];
        nx = ky2->size[1];
        i15 = ky->size[0] * ky->size[1];
        ky->size[0] = vstride;
        ky->size[1] = nx;
        emxEnsureCapacity((emxArray__common *)ky, i15, (int)sizeof(double));
        m = ky1->size[0];
        i15 = ky->size[0] * ky->size[1];
        emxEnsureCapacity((emxArray__common *)ky, i15, (int)sizeof(double));
        loop_ub = ky->size[1];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = ky->size[0];
          for (i16 = 0; i16 < br; i16++) {
            ky->data[i16 + ky->size[0] * i15] = 0.0;
          }
        }

        if ((ky1->size[0] == 0) || (ky2->size[1] == 0)) {
        } else {
          vstride = ky1->size[0] * (ky2->size[1] - 1);
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            i15 = nx + m;
            for (ic = nx; ic + 1 <= i15; ic++) {
              ky->data[ic] = 0.0;
            }

            nx += m;
          }

          br = 0;
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            ar = 0;
            i15 = br + k;
            for (ib = br; ib + 1 <= i15; ib++) {
              if (ky2->data[ib] != 0.0) {
                ia = ar;
                i16 = nx + m;
                for (ic = nx; ic + 1 <= i16; ic++) {
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
      i15 = b->size[0] * b->size[1];
      b->size[0] = ky2->size[1];
      b->size[1] = ky2->size[0];
      emxEnsureCapacity((emxArray__common *)b, i15, (int)sizeof(double));
      loop_ub = ky2->size[0];
      for (i15 = 0; i15 < loop_ub; i15++) {
        br = ky2->size[1];
        for (i16 = 0; i16 < br; i16++) {
          b->data[i16 + b->size[0] * i15] = ky2->data[i15 + ky2->size[0] * i16];
        }
      }

      if ((ky2->size[1] == 1) || (b->size[0] == 1)) {
        i15 = DD->size[0] * DD->size[1];
        DD->size[0] = ky2->size[0];
        DD->size[1] = b->size[1];
        emxEnsureCapacity((emxArray__common *)DD, i15, (int)sizeof(double));
        loop_ub = ky2->size[0];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = b->size[1];
          for (i16 = 0; i16 < br; i16++) {
            DD->data[i15 + DD->size[0] * i16] = 0.0;
            nx = ky2->size[1];
            for (i17 = 0; i17 < nx; i17++) {
              DD->data[i15 + DD->size[0] * i16] += ky2->data[i15 + ky2->size[0] *
                i17] * b->data[i17 + b->size[0] * i16];
            }
          }
        }
      } else {
        k = ky2->size[1];
        vstride = ky2->size[0];
        nx = b->size[1];
        i15 = DD->size[0] * DD->size[1];
        DD->size[0] = vstride;
        DD->size[1] = nx;
        emxEnsureCapacity((emxArray__common *)DD, i15, (int)sizeof(double));
        m = ky2->size[0];
        i15 = DD->size[0] * DD->size[1];
        emxEnsureCapacity((emxArray__common *)DD, i15, (int)sizeof(double));
        loop_ub = DD->size[1];
        for (i15 = 0; i15 < loop_ub; i15++) {
          br = DD->size[0];
          for (i16 = 0; i16 < br; i16++) {
            DD->data[i16 + DD->size[0] * i15] = 0.0;
          }
        }

        if ((ky2->size[0] == 0) || (b->size[1] == 0)) {
        } else {
          vstride = ky2->size[0] * (b->size[1] - 1);
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            i15 = nx + m;
            for (ic = nx; ic + 1 <= i15; ic++) {
              DD->data[ic] = 0.0;
            }

            nx += m;
          }

          br = 0;
          nx = 0;
          while ((m > 0) && (nx <= vstride)) {
            ar = 0;
            i15 = br + k;
            for (ib = br; ib + 1 <= i15; ib++) {
              if (b->data[ib] != 0.0) {
                ia = ar;
                i16 = nx + m;
                for (ic = nx; ic + 1 <= i16; ic++) {
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

      /* U = []; V = []; C = []; %%% latest changes by XIA */
    }
  }

  if (guard2) {
    i15 = ky->size[0] * ky->size[1];
    ky->size[0] = y->size[0];
    ky->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)ky, i15, (int)sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (i15 = 0; i15 < loop_ub; i15++) {
      ky->data[i15] = y->data[i15];
    }

    i15 = ky1->size[0] * ky1->size[1];
    ky1->size[0] = y->size[0];
    ky1->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)ky1, i15, (int)sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (i15 = 0; i15 < loop_ub; i15++) {
      ky1->data[i15] = y->data[i15];
    }

    /* %% latest changes by XIA */
    i15 = ky2->size[0] * ky2->size[1];
    ky2->size[0] = 1;
    ky2->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)ky2, i15, (int)sizeof(double));
    ky2->data[0] = 1.0;

    /* %% latest changes by XIA */
    i15 = DD->size[0] * DD->size[1];
    DD->size[0] = 1;
    DD->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)DD, i15, (int)sizeof(double));
    DD->data[0] = 1.0;

    /* %% latest changes by XIA */
  }

  emxFree_real_T(&c_V);
  emxFree_real_T(&b_V);
  emxFree_real_T(&b_ky);
  emxFree_real_T(&c_qx);
  emxFree_boolean_T(&c_ss);
  emxFree_real_T(&b_Dc);
  emxFree_real_T(&b_s);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&c_x);
  emxFree_real_T(&b_qx);
  emxInit_real_T2(&BB, 3);
  i15 = BB->size[0] * BB->size[1] * BB->size[2];
  BB->size[0] = p;
  BB->size[1] = p;
  BB->size[2] = p;
  emxEnsureCapacity((emxArray__common *)BB, i15, (int)sizeof(double));
  loop_ub = p * p * p;
  for (i15 = 0; i15 < loop_ub; i15++) {
    BB->data[i15] = 0.0;
  }

  emxInit_real_T(&r1, 2);
  eye((double)p, r1);
  loop_ub = r1->size[1];
  for (i15 = 0; i15 < loop_ub; i15++) {
    br = r1->size[0];
    for (i16 = 0; i16 < br; i16++) {
      BB->data[(i16 + BB->size[0] * i15) + BB->size[0] * BB->size[1] * (p - 1)] =
        r1->data[i16 + r1->size[0] * i15];
    }
  }

  emxFree_real_T(&r1);
  emxInit_real_T(&onexi, 2);
  i15 = onexi->size[0] * onexi->size[1];
  onexi->size[0] = n;
  onexi->size[1] = p + 1;
  emxEnsureCapacity((emxArray__common *)onexi, i15, (int)sizeof(double));
  loop_ub = n * (p + 1);
  for (i15 = 0; i15 < loop_ub; i15++) {
    onexi->data[i15] = 1.0;
  }

  emxInit_real_T(&B, 2);
  eye((double)p, B);
  i15 = ss->size[0] * ss->size[1];
  ss->size[0] = p;
  ss->size[1] = p;
  emxEnsureCapacity((emxArray__common *)ss, i15, (int)sizeof(double));
  loop_ub = p * p;
  for (i15 = 0; i15 < loop_ub; i15++) {
    ss->data[i15] = 0.0;
  }

  /* Jan 13 */
  if (p > 1) {
    i15 = (int)((1.0 + (-1.0 - (double)p)) / -1.0);
    ip = 0;
    emxInit_real_T(&D, 2);
    emxInit_real_T(&Ifast, 2);
    emxInit_real_T(&Vfast, 2);
    emxInit_real_T(&xfast, 2);
    emxInit_real_T(&xij, 2);
    emxInit_real_T1(&dxij, 1);
    emxInit_real_T(&abi, 2);
    emxInit_real_T(&dd, 2);
    emxInit_real_T(&kxijy, 2);
    emxInit_real_T(&ddx, 2);
    emxInit_real_T(&tmp, 2);
    emxInit_real_T1(&b_B, 1);
    emxInit_creal_T(&R, 1);
    emxInit_real_T(&b_C, 2);
    emxInit_real_T(&c_C, 2);
    emxInit_real_T(&d_C, 2);
    emxInit_real_T1(&e_C, 1);
    emxInit_int32_T1(&r2, 2);
    emxInit_real_T(&b_y, 2);
    emxInit_real_T(&c_y, 2);
    emxInit_real_T(&d_y, 2);
    emxInitMatrix_cell_wrap_0(reshapes);
    emxInit_real_T(&e_y, 2);
    emxInit_real_T(&f_y, 2);
    emxInit_real_T(&g_y, 2);
    emxInit_real_T(&h_y, 2);
    emxInit_real_T(&b_dd, 2);
    emxInit_real_T(&b_abi, 2);
    emxInit_real_T(&f_C, 2);
    emxInit_real_T(&d_V, 2);
    emxInit_real_T(&b_Vfast, 2);
    emxInit_real_T(&b_xfast, 2);
    emxInit_real_T(&g_C, 2);
    emxInit_real_T1(&e_V, 1);
    emxInit_real_T(&c_B, 2);
    emxInit_int32_T1(&b_Ifast, 2);
    emxInit_int32_T1(&c_Ifast, 2);
    emxInit_real_T1(&b_dc, 1);
    emxInit_real_T(&c_abi, 2);
    emxInit_int32_T1(&d_Ifast, 2);
    emxInit_int32_T1(&e_Ifast, 2);
    exitg1 = false;
    while ((!exitg1) && (ip <= i15 - 1)) {
      b_ip = p - ip;
      b_m = (int)std::floor((double)n / 2.0);
      if (b_ip <= b_m) {
        m = b_ip;
      } else {
        m = b_m;
      }

      K = ((b_ip <= 10) * 5 + ((b_ip <= 20) * (b_ip > 10) << 1)) + (b_ip > 20);
      if (d_strcmp(method)) {
        K = 1;
      }

      for (iter = 0; iter < K; iter++) {
        /* Jan 13 Change for C */
        if (d_strcmp(method)) {
          i16 = V->size[0] * V->size[1];
          V->size[0] = y->size[0];
          V->size[1] = y->size[1];
          emxEnsureCapacity((emxArray__common *)V, i16, (int)sizeof(double));
          loop_ub = y->size[0] * y->size[1];
          for (i16 = 0; i16 < loop_ub; i16++) {
            V->data[i16] = y->data[i16];
          }

          m = 1;
        } else if ((x->size[1] == 1) || (B->size[0] == 1)) {
          i16 = V->size[0] * V->size[1];
          V->size[0] = x->size[0];
          V->size[1] = B->size[1];
          emxEnsureCapacity((emxArray__common *)V, i16, (int)sizeof(double));
          loop_ub = x->size[0];
          for (i16 = 0; i16 < loop_ub; i16++) {
            br = B->size[1];
            for (i17 = 0; i17 < br; i17++) {
              V->data[i16 + V->size[0] * i17] = 0.0;
              nx = x->size[1];
              for (b_m = 0; b_m < nx; b_m++) {
                V->data[i16 + V->size[0] * i17] += x->data[i16 + x->size[0] *
                  b_m] * B->data[b_m + B->size[0] * i17];
              }
            }
          }
        } else {
          k = x->size[1];
          vstride = x->size[0];
          nx = B->size[1];
          i16 = V->size[0] * V->size[1];
          V->size[0] = vstride;
          V->size[1] = nx;
          emxEnsureCapacity((emxArray__common *)V, i16, (int)sizeof(double));
          b_m = x->size[0];
          i16 = V->size[0] * V->size[1];
          emxEnsureCapacity((emxArray__common *)V, i16, (int)sizeof(double));
          loop_ub = V->size[1];
          for (i16 = 0; i16 < loop_ub; i16++) {
            br = V->size[0];
            for (i17 = 0; i17 < br; i17++) {
              V->data[i17 + V->size[0] * i16] = 0.0;
            }
          }

          if ((x->size[0] == 0) || (B->size[1] == 0)) {
          } else {
            vstride = x->size[0] * (B->size[1] - 1);
            nx = 0;
            while ((b_m > 0) && (nx <= vstride)) {
              i16 = nx + b_m;
              for (ic = nx; ic + 1 <= i16; ic++) {
                V->data[ic] = 0.0;
              }

              nx += b_m;
            }

            br = 0;
            nx = 0;
            while ((b_m > 0) && (nx <= vstride)) {
              ar = 0;
              i16 = br + k;
              for (ib = br; ib + 1 <= i16; ib++) {
                if (B->data[ib] != 0.0) {
                  ia = ar;
                  i17 = nx + b_m;
                  for (ic = nx; ic + 1 <= i17; ic++) {
                    ia++;
                    V->data[ic] += B->data[ib] * x->data[ia - 1];
                  }
                }

                ar += b_m;
              }

              br += k;
              nx += b_m;
            }
          }
        }

        if (100 <= n) {
          b_m = 100;
        } else {
          b_m = n;
        }

        pp = std::floor(std::sqrt((double)n));

        /*  Jan 13 */
        if ((b_m >= pp) || rtIsNaN(pp)) {
          c_m = b_m;
        } else {
          c_m = pp;
        }

        unifD(V, c_m, Ifast);

        /*  Jan 13 */
        i16 = b_Ifast->size[0] * b_Ifast->size[1];
        b_Ifast->size[0] = Ifast->size[0];
        b_Ifast->size[1] = Ifast->size[1];
        emxEnsureCapacity((emxArray__common *)b_Ifast, i16, (int)sizeof(int));
        loop_ub = Ifast->size[1];
        for (i16 = 0; i16 < loop_ub; i16++) {
          br = Ifast->size[0];
          for (i17 = 0; i17 < br; i17++) {
            b_Ifast->data[i17 + b_Ifast->size[0] * i16] = (int)Ifast->data[i17 +
              Ifast->size[0] * i16];
          }
        }

        Ifast_idx_0 = Ifast->size[0] * Ifast->size[1];
        loop_ub = V->size[1];
        i16 = Vfast->size[0] * Vfast->size[1];
        Vfast->size[0] = Ifast_idx_0;
        Vfast->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)Vfast, i16, (int)sizeof(double));
        for (i16 = 0; i16 < loop_ub; i16++) {
          for (i17 = 0; i17 < Ifast_idx_0; i17++) {
            Vfast->data[i17 + Vfast->size[0] * i16] = V->data[(b_Ifast->data[i17]
              + V->size[0] * i16) - 1];
          }
        }

        /*  Jan 13 */
        i16 = c_Ifast->size[0] * c_Ifast->size[1];
        c_Ifast->size[0] = Ifast->size[0];
        c_Ifast->size[1] = Ifast->size[1];
        emxEnsureCapacity((emxArray__common *)c_Ifast, i16, (int)sizeof(int));
        loop_ub = Ifast->size[1];
        for (i16 = 0; i16 < loop_ub; i16++) {
          br = Ifast->size[0];
          for (i17 = 0; i17 < br; i17++) {
            c_Ifast->data[i17 + c_Ifast->size[0] * i16] = (int)Ifast->data[i17 +
              Ifast->size[0] * i16];
          }
        }

        Ifast_idx_0 = Ifast->size[0] * Ifast->size[1];
        loop_ub = x->size[1];
        i16 = xfast->size[0] * xfast->size[1];
        xfast->size[0] = Ifast_idx_0;
        xfast->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)xfast, i16, (int)sizeof(double));
        for (i16 = 0; i16 < loop_ub; i16++) {
          for (i17 = 0; i17 < Ifast_idx_0; i17++) {
            xfast->data[i17 + xfast->size[0] * i16] = x->data[(c_Ifast->data[i17]
              + x->size[0] * i16) - 1];
          }
        }

        /*  Jan 13 */
        Ifast_idx_0 = Ifast->size[0] * Ifast->size[1];

        /*  Jan 13 */
        b_std(V, qx);
        pp = 1.2 * mean(qx);
        pp /= rt_powd_snf((double)n, 1.0 / ((double)m + 4.0));
        h2 = 2.0 * pp * pp;
        guard1 = false;
        if (e_strcmp(method)) {
          guard1 = true;
        } else {
          empty_non_axis_sizes = d_strcmp(method);
          b0 = c_strcmp(method);
          if (empty_non_axis_sizes || b0 || (b_ip > 5)) {
            guard1 = true;
          } else {
            i16 = dd->size[0] * dd->size[1];
            dd->size[0] = (int)((double)m * (double)p);
            dd->size[1] = (int)((double)m * (double)p);
            emxEnsureCapacity((emxArray__common *)dd, i16, (int)sizeof(double));
            loop_ub = (int)((double)m * (double)p) * (int)((double)m * (double)p);
            for (i16 = 0; i16 < loop_ub; i16++) {
              dd->data[i16] = 0.0;
            }

            i16 = dc->size[0];
            dc->size[0] = (int)((double)m * (double)p);
            emxEnsureCapacity((emxArray__common *)dc, i16, (int)sizeof(double));
            loop_ub = (int)((double)m * (double)p);
            for (i16 = 0; i16 < loop_ub; i16++) {
              dc->data[i16] = 0.0;
            }

            i16 = D->size[0] * D->size[1];
            D->size[0] = m;
            D->size[1] = m;
            emxEnsureCapacity((emxArray__common *)D, i16, (int)sizeof(double));
            loop_ub = m * m;
            for (i16 = 0; i16 < loop_ub; i16++) {
              D->data[i16] = 0.0;
            }

            /* Jan 13 Change for C */
            for (j = 0; j < Ifast_idx_0; j++) {
              loop_ub = x->size[1];
              i16 = b_xfast->size[0] * b_xfast->size[1];
              b_xfast->size[0] = 1;
              b_xfast->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_xfast, i16, (int)sizeof
                                (double));
              for (i16 = 0; i16 < loop_ub; i16++) {
                b_xfast->data[b_xfast->size[0] * i16] = xfast->data[j +
                  xfast->size[0] * i16];
              }

              repmat(b_xfast, (double)n, xij);
              i16 = xij->size[0] * xij->size[1];
              xij->size[0] = x->size[0];
              xij->size[1] = x->size[1];
              emxEnsureCapacity((emxArray__common *)xij, i16, (int)sizeof(double));
              loop_ub = x->size[0] * x->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                xij->data[i16] = x->data[i16] - xij->data[i16];
              }

              /*  NEW changes */
              loop_ub = V->size[1];
              i16 = b_Vfast->size[0] * b_Vfast->size[1];
              b_Vfast->size[0] = 1;
              b_Vfast->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_Vfast, i16, (int)sizeof
                                (double));
              for (i16 = 0; i16 < loop_ub; i16++) {
                b_Vfast->data[b_Vfast->size[0] * i16] = Vfast->data[j +
                  Vfast->size[0] * i16];
              }

              repmat(b_Vfast, (double)n, a);
              i16 = d_V->size[0] * d_V->size[1];
              d_V->size[0] = V->size[0];
              d_V->size[1] = V->size[1];
              emxEnsureCapacity((emxArray__common *)d_V, i16, (int)sizeof(double));
              loop_ub = V->size[0] * V->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                d_V->data[i16] = V->data[i16] - a->data[i16];
              }

              power(d_V, a);
              sum(a, dxij);

              /*  NEW changes */
              i16 = s->size[0];
              s->size[0] = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)s, i16, (int)sizeof(double));
              loop_ub = dxij->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                s->data[i16] = dxij->data[i16];
              }

              f_sort(s);
              if ((h2 >= s->data[(int)(2.0 * (double)m) - 1]) || rtIsNaN(s->
                   data[(int)(2.0 * (double)m) - 1])) {
                pp = h2;
              } else {
                pp = s->data[(int)(2.0 * (double)m) - 1];
              }

              i16 = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)dxij, i16, (int)sizeof
                                (double));
              loop_ub = dxij->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                dxij->data[i16] = -dxij->data[i16] / pp;
              }

              b_exp(dxij);
              b_repmat(dxij, (double)p + 1.0, C);
              if ((xij->size[1] == 1) || (B->size[0] == 1)) {
                i16 = d_y->size[0] * d_y->size[1];
                d_y->size[0] = xij->size[0];
                d_y->size[1] = B->size[1];
                emxEnsureCapacity((emxArray__common *)d_y, i16, (int)sizeof
                                  (double));
                loop_ub = xij->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = B->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    d_y->data[i16 + d_y->size[0] * i17] = 0.0;
                    nx = xij->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      d_y->data[i16 + d_y->size[0] * i17] += xij->data[i16 +
                        xij->size[0] * b_m] * B->data[b_m + B->size[0] * i17];
                    }
                  }
                }
              } else {
                k = xij->size[1];
                vstride = xij->size[0];
                nx = B->size[1];
                i16 = d_y->size[0] * d_y->size[1];
                d_y->size[0] = vstride;
                d_y->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)d_y, i16, (int)sizeof
                                  (double));
                b_m = xij->size[0];
                i16 = d_y->size[0] * d_y->size[1];
                emxEnsureCapacity((emxArray__common *)d_y, i16, (int)sizeof
                                  (double));
                loop_ub = d_y->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = d_y->size[0];
                  for (i17 = 0; i17 < br; i17++) {
                    d_y->data[i17 + d_y->size[0] * i16] = 0.0;
                  }
                }

                if ((xij->size[0] == 0) || (B->size[1] == 0)) {
                } else {
                  vstride = xij->size[0] * (B->size[1] - 1);
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    i16 = nx + b_m;
                    for (ic = nx; ic + 1 <= i16; ic++) {
                      d_y->data[ic] = 0.0;
                    }

                    nx += b_m;
                  }

                  br = 0;
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (B->data[ib] != 0.0) {
                        ia = ar;
                        i17 = nx + b_m;
                        for (ic = nx; ic + 1 <= i17; ic++) {
                          ia++;
                          d_y->data[ic] += B->data[ib] * xij->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    nx += b_m;
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
                b_m = d_y->size[1];
              } else {
                b_m = 0;
              }

              if (empty_non_axis_sizes || (!(n == 0))) {
                vstride = 1;
              } else {
                vstride = 0;
              }

              i16 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
              reshapes[1].f1->size[0] = nx;
              reshapes[1].f1->size[1] = vstride;
              emxEnsureCapacity((emxArray__common *)reshapes[1].f1, i16, (int)
                                sizeof(double));
              loop_ub = nx * vstride;
              for (i16 = 0; i16 < loop_ub; i16++) {
                reshapes[1].f1->data[i16] = 1.0;
              }

              i16 = onexi->size[0] * onexi->size[1];
              onexi->size[0] = nx;
              onexi->size[1] = b_m + reshapes[1].f1->size[1];
              emxEnsureCapacity((emxArray__common *)onexi, i16, (int)sizeof
                                (double));
              for (i16 = 0; i16 < b_m; i16++) {
                for (i17 = 0; i17 < nx; i17++) {
                  onexi->data[i17 + onexi->size[0] * i16] = d_y->data[i17 + nx *
                    i16];
                }
              }

              loop_ub = reshapes[1].f1->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = reshapes[1].f1->size[0];
                for (i17 = 0; i17 < br; i17++) {
                  onexi->data[i17 + onexi->size[0] * (i16 + b_m)] = reshapes[1].
                    f1->data[i17 + reshapes[1].f1->size[0] * i16];
                }
              }

              i16 = U->size[0] * U->size[1];
              U->size[0] = onexi->size[1];
              U->size[1] = onexi->size[0];
              emxEnsureCapacity((emxArray__common *)U, i16, (int)sizeof(double));
              loop_ub = onexi->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = onexi->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  U->data[i17 + U->size[0] * i16] = onexi->data[i16 +
                    onexi->size[0] * i17] * C->data[i16 + C->size[0] * i17];
                }
              }

              if ((U->size[1] == 1) || (onexi->size[0] == 1)) {
                i16 = d_C->size[0] * d_C->size[1];
                d_C->size[0] = U->size[0];
                d_C->size[1] = onexi->size[1];
                emxEnsureCapacity((emxArray__common *)d_C, i16, (int)sizeof
                                  (double));
                loop_ub = U->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = onexi->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    d_C->data[i16 + d_C->size[0] * i17] = 0.0;
                    nx = U->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      d_C->data[i16 + d_C->size[0] * i17] += U->data[i16 +
                        U->size[0] * b_m] * onexi->data[b_m + onexi->size[0] *
                        i17];
                    }
                  }
                }
              } else {
                k = U->size[1];
                vstride = U->size[0];
                nx = onexi->size[1];
                i16 = d_C->size[0] * d_C->size[1];
                d_C->size[0] = vstride;
                d_C->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)d_C, i16, (int)sizeof
                                  (double));
                b_m = U->size[0];
                i16 = d_C->size[0] * d_C->size[1];
                emxEnsureCapacity((emxArray__common *)d_C, i16, (int)sizeof
                                  (double));
                loop_ub = d_C->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = d_C->size[0];
                  for (i17 = 0; i17 < br; i17++) {
                    d_C->data[i17 + d_C->size[0] * i16] = 0.0;
                  }
                }

                if ((U->size[0] == 0) || (onexi->size[1] == 0)) {
                } else {
                  vstride = U->size[0] * (onexi->size[1] - 1);
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    i16 = nx + b_m;
                    for (ic = nx; ic + 1 <= i16; ic++) {
                      d_C->data[ic] = 0.0;
                    }

                    nx += b_m;
                  }

                  br = 0;
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (onexi->data[ib] != 0.0) {
                        ia = ar;
                        i17 = nx + b_m;
                        for (ic = nx; ic + 1 <= i17; ic++) {
                          ia++;
                          d_C->data[ic] += onexi->data[ib] * U->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    nx += b_m;
                  }
                }
              }

              eye((double)B->size[1] + 1.0, a);
              i16 = f_C->size[0] * f_C->size[1];
              f_C->size[0] = d_C->size[0];
              f_C->size[1] = d_C->size[1];
              emxEnsureCapacity((emxArray__common *)f_C, i16, (int)sizeof(double));
              loop_ub = d_C->size[0] * d_C->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                f_C->data[i16] = d_C->data[i16] + a->data[i16] / (double)n;
              }

              inv(f_C, ss);
              if ((U->size[1] == 1) || (ky1->size[0] == 1)) {
                i16 = e_y->size[0] * e_y->size[1];
                e_y->size[0] = U->size[0];
                e_y->size[1] = ky1->size[1];
                emxEnsureCapacity((emxArray__common *)e_y, i16, (int)sizeof
                                  (double));
                loop_ub = U->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = ky1->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    e_y->data[i16 + e_y->size[0] * i17] = 0.0;
                    nx = U->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      e_y->data[i16 + e_y->size[0] * i17] += U->data[i16 +
                        U->size[0] * b_m] * ky1->data[b_m + ky1->size[0] * i17];
                    }
                  }
                }
              } else {
                k = U->size[1];
                vstride = U->size[0];
                nx = ky1->size[1];
                i16 = e_y->size[0] * e_y->size[1];
                e_y->size[0] = vstride;
                e_y->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)e_y, i16, (int)sizeof
                                  (double));
                b_m = U->size[0];
                i16 = e_y->size[0] * e_y->size[1];
                emxEnsureCapacity((emxArray__common *)e_y, i16, (int)sizeof
                                  (double));
                loop_ub = e_y->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = e_y->size[0];
                  for (i17 = 0; i17 < br; i17++) {
                    e_y->data[i17 + e_y->size[0] * i16] = 0.0;
                  }
                }

                if ((U->size[0] == 0) || (ky1->size[1] == 0)) {
                } else {
                  vstride = U->size[0] * (ky1->size[1] - 1);
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    i16 = nx + b_m;
                    for (ic = nx; ic + 1 <= i16; ic++) {
                      e_y->data[ic] = 0.0;
                    }

                    nx += b_m;
                  }

                  br = 0;
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (ky1->data[ib] != 0.0) {
                        ia = ar;
                        i17 = nx + b_m;
                        for (ic = nx; ic + 1 <= i17; ic++) {
                          ia++;
                          e_y->data[ic] += ky1->data[ib] * U->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    nx += b_m;
                  }
                }
              }

              if ((ss->size[1] == 1) || (e_y->size[0] == 1)) {
                i16 = f_y->size[0] * f_y->size[1];
                f_y->size[0] = ss->size[0];
                f_y->size[1] = e_y->size[1];
                emxEnsureCapacity((emxArray__common *)f_y, i16, (int)sizeof
                                  (double));
                loop_ub = ss->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = e_y->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    f_y->data[i16 + f_y->size[0] * i17] = 0.0;
                    nx = ss->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      f_y->data[i16 + f_y->size[0] * i17] += ss->data[i16 +
                        ss->size[0] * b_m] * e_y->data[b_m + e_y->size[0] * i17];
                    }
                  }
                }
              } else {
                k = ss->size[1];
                vstride = ss->size[0];
                nx = e_y->size[1];
                i16 = f_y->size[0] * f_y->size[1];
                f_y->size[0] = vstride;
                f_y->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)f_y, i16, (int)sizeof
                                  (double));
                b_m = ss->size[0];
                i16 = f_y->size[0] * f_y->size[1];
                emxEnsureCapacity((emxArray__common *)f_y, i16, (int)sizeof
                                  (double));
                loop_ub = f_y->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = f_y->size[0];
                  for (i17 = 0; i17 < br; i17++) {
                    f_y->data[i17 + f_y->size[0] * i16] = 0.0;
                  }
                }

                if ((ss->size[0] == 0) || (e_y->size[1] == 0)) {
                } else {
                  vstride = ss->size[0] * (e_y->size[1] - 1);
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    i16 = nx + b_m;
                    for (ic = nx; ic + 1 <= i16; ic++) {
                      f_y->data[ic] = 0.0;
                    }

                    nx += b_m;
                  }

                  br = 0;
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (e_y->data[ib] != 0.0) {
                        ia = ar;
                        i17 = nx + b_m;
                        for (ic = nx; ic + 1 <= i17; ic++) {
                          ia++;
                          f_y->data[ic] += e_y->data[ib] * ss->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    nx += b_m;
                  }
                }
              }

              if ((f_y->size[1] == 1) || (ky2->size[0] == 1)) {
                i16 = abi->size[0] * abi->size[1];
                abi->size[0] = f_y->size[0];
                abi->size[1] = ky2->size[1];
                emxEnsureCapacity((emxArray__common *)abi, i16, (int)sizeof
                                  (double));
                loop_ub = f_y->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = ky2->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    abi->data[i16 + abi->size[0] * i17] = 0.0;
                    nx = f_y->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      abi->data[i16 + abi->size[0] * i17] += f_y->data[i16 +
                        f_y->size[0] * b_m] * ky2->data[b_m + ky2->size[0] * i17];
                    }
                  }
                }
              } else {
                k = f_y->size[1];
                vstride = f_y->size[0];
                nx = ky2->size[1];
                i16 = abi->size[0] * abi->size[1];
                abi->size[0] = vstride;
                abi->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)abi, i16, (int)sizeof
                                  (double));
                b_m = f_y->size[0];
                i16 = abi->size[0] * abi->size[1];
                emxEnsureCapacity((emxArray__common *)abi, i16, (int)sizeof
                                  (double));
                loop_ub = abi->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = abi->size[0];
                  for (i17 = 0; i17 < br; i17++) {
                    abi->data[i17 + abi->size[0] * i16] = 0.0;
                  }
                }

                if ((f_y->size[0] == 0) || (ky2->size[1] == 0)) {
                } else {
                  vstride = f_y->size[0] * (ky2->size[1] - 1);
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    i16 = nx + b_m;
                    for (ic = nx; ic + 1 <= i16; ic++) {
                      abi->data[ic] = 0.0;
                    }

                    nx += b_m;
                  }

                  br = 0;
                  nx = 0;
                  while ((b_m > 0) && (nx <= vstride)) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (ky2->data[ib] != 0.0) {
                        ia = ar;
                        i17 = nx + b_m;
                        for (ic = nx; ic + 1 <= i17; ic++) {
                          ia++;
                          abi->data[ic] += ky2->data[ib] * f_y->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    nx += b_m;
                  }
                }
              }

              /* Jan 7 */
              i16 = yj->size[0] * yj->size[1];
              yj->size[0] = xij->size[1];
              yj->size[1] = xij->size[0];
              emxEnsureCapacity((emxArray__common *)yj, i16, (int)sizeof(double));
              loop_ub = xij->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = xij->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  yj->data[i17 + yj->size[0] * i16] = xij->data[i16 + xij->size
                    [0] * i17] * C->data[i16 + C->size[0] * i17];
                }
              }

              loop_ub = abi->size[1];
              b_m = (int)(m + 1U);
              i16 = b_abi->size[0] * b_abi->size[1];
              b_abi->size[0] = 1;
              b_abi->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_abi, i16, (int)sizeof
                                (double));
              for (i16 = 0; i16 < loop_ub; i16++) {
                b_abi->data[b_abi->size[0] * i16] = abi->data[(b_m + abi->size[0]
                  * i16) - 1];
              }

              repmat(b_abi, (double)n, b);
              i16 = b->size[0] * b->size[1];
              b->size[0] = ky->size[0];
              b->size[1] = ky->size[1];
              emxEnsureCapacity((emxArray__common *)b, i16, (int)sizeof(double));
              loop_ub = ky->size[0] * ky->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                b->data[i16] = ky->data[i16] - b->data[i16];
              }

              if ((yj->size[1] == 1) || (b->size[0] == 1)) {
                i16 = kxijy->size[0] * kxijy->size[1];
                kxijy->size[0] = yj->size[0];
                kxijy->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)kxijy, i16, (int)sizeof
                                  (double));
                loop_ub = yj->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = b->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    kxijy->data[i16 + kxijy->size[0] * i17] = 0.0;
                    nx = yj->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      kxijy->data[i16 + kxijy->size[0] * i17] += yj->data[i16 +
                        yj->size[0] * b_m] * b->data[b_m + b->size[0] * i17];
                    }
                  }
                }
              } else {
                k = yj->size[1];
                vstride = yj->size[0];
                nx = b->size[1];
                i16 = kxijy->size[0] * kxijy->size[1];
                kxijy->size[0] = vstride;
                kxijy->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)kxijy, i16, (int)sizeof
                                  (double));
                b_m = yj->size[0];
                i16 = kxijy->size[0] * kxijy->size[1];
                emxEnsureCapacity((emxArray__common *)kxijy, i16, (int)sizeof
                                  (double));
                loop_ub = kxijy->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = kxijy->size[0];
                  for (i17 = 0; i17 < br; i17++) {
                    kxijy->data[i17 + kxijy->size[0] * i16] = 0.0;
                  }
                }

                if (b->size[1] != 0) {
                  vstride = yj->size[0] * (b->size[1] - 1);
                  for (nx = 0; nx <= vstride; nx += b_m) {
                    i16 = nx + b_m;
                    for (ic = nx; ic + 1 <= i16; ic++) {
                      kxijy->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (nx = 0; nx <= vstride; nx += b_m) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (b->data[ib] != 0.0) {
                        ia = ar;
                        i17 = nx + b_m;
                        for (ic = nx; ic + 1 <= i17; ic++) {
                          ia++;
                          kxijy->data[ic] += b->data[ib] * yj->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                  }
                }
              }

              if ((yj->size[1] == 1) || (xij->size[0] == 1)) {
                i16 = ddx->size[0] * ddx->size[1];
                ddx->size[0] = yj->size[0];
                ddx->size[1] = xij->size[1];
                emxEnsureCapacity((emxArray__common *)ddx, i16, (int)sizeof
                                  (double));
                loop_ub = yj->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = xij->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    ddx->data[i16 + ddx->size[0] * i17] = 0.0;
                    nx = yj->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      ddx->data[i16 + ddx->size[0] * i17] += yj->data[i16 +
                        yj->size[0] * b_m] * xij->data[b_m + xij->size[0] * i17];
                    }
                  }
                }
              } else {
                k = yj->size[1];
                vstride = yj->size[0];
                nx = xij->size[1];
                i16 = ddx->size[0] * ddx->size[1];
                ddx->size[0] = vstride;
                ddx->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)ddx, i16, (int)sizeof
                                  (double));
                b_m = yj->size[0];
                i16 = ddx->size[0] * ddx->size[1];
                emxEnsureCapacity((emxArray__common *)ddx, i16, (int)sizeof
                                  (double));
                loop_ub = ddx->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = ddx->size[0];
                  for (i17 = 0; i17 < br; i17++) {
                    ddx->data[i17 + ddx->size[0] * i16] = 0.0;
                  }
                }

                if (xij->size[1] != 0) {
                  vstride = yj->size[0] * (xij->size[1] - 1);
                  for (nx = 0; nx <= vstride; nx += b_m) {
                    i16 = nx + b_m;
                    for (ic = nx; ic + 1 <= i16; ic++) {
                      ddx->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (nx = 0; nx <= vstride; nx += b_m) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (xij->data[ib] != 0.0) {
                        ia = ar;
                        i17 = nx + b_m;
                        for (ic = nx; ic + 1 <= i17; ic++) {
                          ia++;
                          ddx->data[ic] += xij->data[ib] * yj->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                  }
                }
              }

              for (k1 = 0; k1 < m; k1++) {
                b_a = ((1.0 + (double)k1) - 1.0) * (double)p + 1.0;
                d = (1.0 + (double)k1) * (double)p;
                if (d < b_a) {
                  i16 = qx->size[0] * qx->size[1];
                  qx->size[0] = 1;
                  qx->size[1] = 0;
                  emxEnsureCapacity((emxArray__common *)qx, i16, (int)sizeof
                                    (double));
                } else if (b_a == b_a) {
                  i16 = qx->size[0] * qx->size[1];
                  qx->size[0] = 1;
                  qx->size[1] = (int)(d - b_a) + 1;
                  emxEnsureCapacity((emxArray__common *)qx, i16, (int)sizeof
                                    (double));
                  loop_ub = (int)(d - b_a);
                  for (i16 = 0; i16 <= loop_ub; i16++) {
                    qx->data[qx->size[0] * i16] = b_a + (double)i16;
                  }
                } else {
                  ndbl = std::floor((d - b_a) + 0.5);
                  apnd = b_a + ndbl;
                  cdiff = apnd - d;
                  if (b_a >= d) {
                    pp = b_a;
                  } else {
                    pp = d;
                  }

                  if (std::abs(cdiff) < 4.4408920985006262E-16 * pp) {
                    ndbl++;
                    apnd = d;
                  } else if (cdiff > 0.0) {
                    apnd = b_a + (ndbl - 1.0);
                  } else {
                    ndbl++;
                  }

                  if (ndbl >= 0.0) {
                    nx = (int)ndbl;
                  } else {
                    nx = 0;
                  }

                  i16 = qx->size[0] * qx->size[1];
                  qx->size[0] = 1;
                  qx->size[1] = nx;
                  emxEnsureCapacity((emxArray__common *)qx, i16, (int)sizeof
                                    (double));
                  if (nx > 0) {
                    qx->data[0] = b_a;
                    if (nx > 1) {
                      qx->data[nx - 1] = apnd;
                      vstride = (nx - 1) / 2;
                      for (k = 1; k < vstride; k++) {
                        qx->data[k] = b_a + (double)k;
                        qx->data[(nx - k) - 1] = apnd - (double)k;
                      }

                      if (vstride << 1 == nx - 1) {
                        qx->data[vstride] = (b_a + apnd) / 2.0;
                      } else {
                        qx->data[vstride] = b_a + (double)vstride;
                        qx->data[vstride + 1] = apnd - (double)vstride;
                      }
                    }
                  }
                }

                loop_ub = abi->size[1];
                i16 = s->size[0];
                s->size[0] = loop_ub;
                emxEnsureCapacity((emxArray__common *)s, i16, (int)sizeof(double));
                for (i16 = 0; i16 < loop_ub; i16++) {
                  s->data[i16] = abi->data[k1 + abi->size[0] * i16];
                }

                if ((kxijy->size[1] == 1) || (s->size[0] == 1)) {
                  i16 = e_C->size[0];
                  e_C->size[0] = kxijy->size[0];
                  emxEnsureCapacity((emxArray__common *)e_C, i16, (int)sizeof
                                    (double));
                  loop_ub = kxijy->size[0];
                  for (i16 = 0; i16 < loop_ub; i16++) {
                    e_C->data[i16] = 0.0;
                    br = kxijy->size[1];
                    for (i17 = 0; i17 < br; i17++) {
                      e_C->data[i16] += kxijy->data[i16 + kxijy->size[0] * i17] *
                        s->data[i17];
                    }
                  }
                } else {
                  k = kxijy->size[1];
                  a_idx_0 = (unsigned int)kxijy->size[0];
                  i16 = e_C->size[0];
                  e_C->size[0] = (int)a_idx_0;
                  emxEnsureCapacity((emxArray__common *)e_C, i16, (int)sizeof
                                    (double));
                  b_m = kxijy->size[0];
                  vstride = e_C->size[0];
                  i16 = e_C->size[0];
                  e_C->size[0] = vstride;
                  emxEnsureCapacity((emxArray__common *)e_C, i16, (int)sizeof
                                    (double));
                  for (i16 = 0; i16 < vstride; i16++) {
                    e_C->data[i16] = 0.0;
                  }

                  nx = 0;
                  while (nx <= 0) {
                    for (ic = 1; ic <= b_m; ic++) {
                      e_C->data[ic - 1] = 0.0;
                    }

                    nx = b_m;
                  }

                  br = 0;
                  nx = 0;
                  while (nx <= 0) {
                    ar = 0;
                    i16 = br + k;
                    for (ib = br; ib + 1 <= i16; ib++) {
                      if (s->data[ib] != 0.0) {
                        ia = ar;
                        for (ic = 0; ic + 1 <= b_m; ic++) {
                          ia++;
                          e_C->data[ic] += s->data[ib] * kxijy->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    nx = b_m;
                  }
                }

                i16 = r2->size[0] * r2->size[1];
                r2->size[0] = 1;
                r2->size[1] = qx->size[1];
                emxEnsureCapacity((emxArray__common *)r2, i16, (int)sizeof(int));
                loop_ub = qx->size[0] * qx->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  r2->data[i16] = (int)qx->data[i16];
                }

                i16 = b_dc->size[0];
                b_dc->size[0] = qx->size[1];
                emxEnsureCapacity((emxArray__common *)b_dc, i16, (int)sizeof
                                  (double));
                loop_ub = qx->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  b_dc->data[i16] = dc->data[(int)qx->data[qx->size[0] * i16] -
                    1] + e_C->data[i16];
                }

                loop_ub = r2->size[1];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  dc->data[r2->data[r2->size[0] * i16] - 1] = b_dc->data[(*(int
                    (*)[2])r2->size)[0] * i16];
                }
              }

              if (1 > m) {
                loop_ub = 0;
                br = 0;
              } else {
                loop_ub = m;
                br = m;
              }

              nx = abi->size[1];
              i16 = a->size[0] * a->size[1];
              a->size[0] = loop_ub;
              a->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)a, i16, (int)sizeof(double));
              for (i16 = 0; i16 < nx; i16++) {
                for (i17 = 0; i17 < loop_ub; i17++) {
                  a->data[i17 + a->size[0] * i16] = abi->data[i17 + abi->size[0]
                    * i16];
                }
              }

              nx = abi->size[1];
              i16 = b->size[0] * b->size[1];
              b->size[0] = nx;
              b->size[1] = br;
              emxEnsureCapacity((emxArray__common *)b, i16, (int)sizeof(double));
              for (i16 = 0; i16 < br; i16++) {
                for (i17 = 0; i17 < nx; i17++) {
                  b->data[i17 + b->size[0] * i16] = abi->data[i16 + abi->size[0]
                    * i17];
                }
              }

              i16 = abi->size[1];
              if ((i16 == 1) || (b->size[0] == 1)) {
                i16 = tmp->size[0] * tmp->size[1];
                tmp->size[0] = a->size[0];
                tmp->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)tmp, i16, (int)sizeof
                                  (double));
                loop_ub = a->size[0];
                for (i16 = 0; i16 < loop_ub; i16++) {
                  br = b->size[1];
                  for (i17 = 0; i17 < br; i17++) {
                    tmp->data[i16 + tmp->size[0] * i17] = 0.0;
                    nx = a->size[1];
                    for (b_m = 0; b_m < nx; b_m++) {
                      tmp->data[i16 + tmp->size[0] * i17] += a->data[i16 +
                        a->size[0] * b_m] * b->data[b_m + b->size[0] * i17];
                    }
                  }
                }
              } else {
                i16 = abi->size[1];
                nx = b->size[1];
                i17 = tmp->size[0] * tmp->size[1];
                tmp->size[0] = loop_ub;
                tmp->size[1] = nx;
                emxEnsureCapacity((emxArray__common *)tmp, i17, (int)sizeof
                                  (double));
                i17 = tmp->size[0] * tmp->size[1];
                emxEnsureCapacity((emxArray__common *)tmp, i17, (int)sizeof
                                  (double));
                br = tmp->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  nx = tmp->size[0];
                  for (b_m = 0; b_m < nx; b_m++) {
                    tmp->data[b_m + tmp->size[0] * i17] = 0.0;
                  }
                }

                if ((loop_ub == 0) || (b->size[1] == 0)) {
                } else {
                  vstride = loop_ub * (b->size[1] - 1);
                  nx = 0;
                  while ((loop_ub > 0) && (nx <= vstride)) {
                    i17 = nx + loop_ub;
                    for (ic = nx; ic + 1 <= i17; ic++) {
                      tmp->data[ic] = 0.0;
                    }

                    nx += loop_ub;
                  }

                  br = 0;
                  nx = 0;
                  while ((loop_ub > 0) && (nx <= vstride)) {
                    ar = 0;
                    i17 = br + i16;
                    for (ib = br; ib + 1 <= i17; ib++) {
                      if (b->data[ib] != 0.0) {
                        ia = ar;
                        b_m = nx + loop_ub;
                        for (ic = nx; ic + 1 <= b_m; ic++) {
                          ia++;
                          tmp->data[ic] += b->data[ib] * a->data[ia - 1];
                        }
                      }

                      ar += loop_ub;
                    }

                    br += i16;
                    nx += loop_ub;
                  }
                }
              }

              kron(tmp, ddx, a);
              i16 = dd->size[0] * dd->size[1];
              emxEnsureCapacity((emxArray__common *)dd, i16, (int)sizeof(double));
              vstride = dd->size[0];
              nx = dd->size[1];
              loop_ub = vstride * nx;
              for (i16 = 0; i16 < loop_ub; i16++) {
                dd->data[i16] += a->data[i16];
              }

              i16 = D->size[0] * D->size[1];
              emxEnsureCapacity((emxArray__common *)D, i16, (int)sizeof(double));
              vstride = D->size[0];
              nx = D->size[1];
              loop_ub = vstride * nx;
              for (i16 = 0; i16 < loop_ub; i16++) {
                D->data[i16] += tmp->data[i16];
              }
            }

            eye((double)dc->size[0], a);
            i16 = b_dd->size[0] * b_dd->size[1];
            b_dd->size[0] = dd->size[0];
            b_dd->size[1] = dd->size[1];
            emxEnsureCapacity((emxArray__common *)b_dd, i16, (int)sizeof(double));
            loop_ub = dd->size[0] * dd->size[1];
            for (i16 = 0; i16 < loop_ub; i16++) {
              b_dd->data[i16] = dd->data[i16] + a->data[i16] / (double)n;
            }

            inv(b_dd, a);
            if ((a->size[1] == 1) || (dc->size[0] == 1)) {
              i16 = b_B->size[0];
              b_B->size[0] = a->size[0];
              emxEnsureCapacity((emxArray__common *)b_B, i16, (int)sizeof(double));
              loop_ub = a->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                b_B->data[i16] = 0.0;
                br = a->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  b_B->data[i16] += a->data[i16 + a->size[0] * i17] * dc->
                    data[i17];
                }
              }
            } else {
              k = a->size[1];
              a_idx_0 = (unsigned int)a->size[0];
              i16 = b_B->size[0];
              b_B->size[0] = (int)a_idx_0;
              emxEnsureCapacity((emxArray__common *)b_B, i16, (int)sizeof(double));
              b_m = a->size[0];
              vstride = b_B->size[0];
              i16 = b_B->size[0];
              b_B->size[0] = vstride;
              emxEnsureCapacity((emxArray__common *)b_B, i16, (int)sizeof(double));
              for (i16 = 0; i16 < vstride; i16++) {
                b_B->data[i16] = 0.0;
              }

              if (a->size[0] != 0) {
                nx = 0;
                while ((b_m > 0) && (nx <= 0)) {
                  for (ic = 1; ic <= b_m; ic++) {
                    b_B->data[ic - 1] = 0.0;
                  }

                  nx = b_m;
                }

                br = 0;
                nx = 0;
                while ((b_m > 0) && (nx <= 0)) {
                  ar = 0;
                  i16 = br + k;
                  for (ib = br; ib + 1 <= i16; ib++) {
                    if (dc->data[ib] != 0.0) {
                      ia = ar;
                      for (ic = 0; ic + 1 <= b_m; ic++) {
                        ia++;
                        b_B->data[ic] += dc->data[ib] * a->data[ia - 1];
                      }
                    }

                    ar += b_m;
                  }

                  br += k;
                  nx = b_m;
                }
              }
            }

            i16 = ss->size[0] * ss->size[1];
            ss->size[0] = p;
            ss->size[1] = m;
            emxEnsureCapacity((emxArray__common *)ss, i16, (int)sizeof(double));
            for (k = 0; k + 1 <= b_B->size[0]; k++) {
              ss->data[k] = b_B->data[k];
            }

            if ((ss->size[1] == 1) || (D->size[0] == 1)) {
              i16 = b_y->size[0] * b_y->size[1];
              b_y->size[0] = ss->size[0];
              b_y->size[1] = D->size[1];
              emxEnsureCapacity((emxArray__common *)b_y, i16, (int)sizeof(double));
              loop_ub = ss->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = D->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  b_y->data[i16 + b_y->size[0] * i17] = 0.0;
                  nx = ss->size[1];
                  for (b_m = 0; b_m < nx; b_m++) {
                    b_y->data[i16 + b_y->size[0] * i17] += ss->data[i16 +
                      ss->size[0] * b_m] * D->data[b_m + D->size[0] * i17];
                  }
                }
              }
            } else {
              k = ss->size[1];
              vstride = ss->size[0];
              nx = D->size[1];
              i16 = b_y->size[0] * b_y->size[1];
              b_y->size[0] = vstride;
              b_y->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)b_y, i16, (int)sizeof(double));
              b_m = ss->size[0];
              i16 = b_y->size[0] * b_y->size[1];
              emxEnsureCapacity((emxArray__common *)b_y, i16, (int)sizeof(double));
              loop_ub = b_y->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = b_y->size[0];
                for (i17 = 0; i17 < br; i17++) {
                  b_y->data[i17 + b_y->size[0] * i16] = 0.0;
                }
              }

              if (D->size[1] != 0) {
                vstride = ss->size[0] * (D->size[1] - 1);
                for (nx = 0; nx <= vstride; nx += b_m) {
                  i16 = nx + b_m;
                  for (ic = nx; ic + 1 <= i16; ic++) {
                    b_y->data[ic] = 0.0;
                  }
                }

                br = 0;
                for (nx = 0; nx <= vstride; nx += b_m) {
                  ar = 0;
                  i16 = br + k;
                  for (ib = br; ib + 1 <= i16; ib++) {
                    if (D->data[ib] != 0.0) {
                      ia = ar;
                      i17 = nx + b_m;
                      for (ic = nx; ic + 1 <= i17; ic++) {
                        ia++;
                        b_y->data[ic] += D->data[ib] * ss->data[ia - 1];
                      }
                    }

                    ar += b_m;
                  }

                  br += k;
                }
              }
            }

            i16 = b->size[0] * b->size[1];
            b->size[0] = ss->size[1];
            b->size[1] = ss->size[0];
            emxEnsureCapacity((emxArray__common *)b, i16, (int)sizeof(double));
            loop_ub = ss->size[0];
            for (i16 = 0; i16 < loop_ub; i16++) {
              br = ss->size[1];
              for (i17 = 0; i17 < br; i17++) {
                b->data[i17 + b->size[0] * i16] = ss->data[i16 + ss->size[0] *
                  i17];
              }
            }

            if ((b_y->size[1] == 1) || (b->size[0] == 1)) {
              i16 = c_y->size[0] * c_y->size[1];
              c_y->size[0] = b_y->size[0];
              c_y->size[1] = b->size[1];
              emxEnsureCapacity((emxArray__common *)c_y, i16, (int)sizeof(double));
              loop_ub = b_y->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = b->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  c_y->data[i16 + c_y->size[0] * i17] = 0.0;
                  nx = b_y->size[1];
                  for (b_m = 0; b_m < nx; b_m++) {
                    c_y->data[i16 + c_y->size[0] * i17] += b_y->data[i16 +
                      b_y->size[0] * b_m] * b->data[b_m + b->size[0] * i17];
                  }
                }
              }
            } else {
              k = b_y->size[1];
              vstride = b_y->size[0];
              nx = b->size[1];
              i16 = c_y->size[0] * c_y->size[1];
              c_y->size[0] = vstride;
              c_y->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)c_y, i16, (int)sizeof(double));
              b_m = b_y->size[0];
              i16 = c_y->size[0] * c_y->size[1];
              emxEnsureCapacity((emxArray__common *)c_y, i16, (int)sizeof(double));
              loop_ub = c_y->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = c_y->size[0];
                for (i17 = 0; i17 < br; i17++) {
                  c_y->data[i17 + c_y->size[0] * i16] = 0.0;
                }
              }

              vstride = b_y->size[0] * (b->size[1] - 1);
              for (nx = 0; nx <= vstride; nx += b_m) {
                i16 = nx + b_m;
                for (ic = nx; ic + 1 <= i16; ic++) {
                  c_y->data[ic] = 0.0;
                }
              }

              br = 0;
              for (nx = 0; nx <= vstride; nx += b_m) {
                ar = 0;
                i16 = br + k;
                for (ib = br; ib + 1 <= i16; ib++) {
                  if (b->data[ib] != 0.0) {
                    ia = ar;
                    i17 = nx + b_m;
                    for (ic = nx; ic + 1 <= i17; ic++) {
                      ia++;
                      c_y->data[ic] += b->data[ib] * b_y->data[ia - 1];
                    }
                  }

                  ar += b_m;
                }

                br += k;
              }
            }

            eig(c_y, Vc, Dc);
            i16 = B->size[0] * B->size[1];
            B->size[0] = Vc->size[0];
            B->size[1] = Vc->size[1];
            emxEnsureCapacity((emxArray__common *)B, i16, (int)sizeof(double));
            loop_ub = Vc->size[0] * Vc->size[1];
            for (i16 = 0; i16 < loop_ub; i16++) {
              B->data[i16] = Vc->data[i16].re;
            }

            /* Jan 13 Change for C */
            c_diag(Dc, R);
            g_sort(R, iidx);
            i16 = s->size[0];
            s->size[0] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)s, i16, (int)sizeof(double));
            loop_ub = iidx->size[0];
            for (i16 = 0; i16 < loop_ub; i16++) {
              s->data[i16] = iidx->data[i16];
            }

            loop_ub = B->size[0];
            i16 = ss->size[0] * ss->size[1];
            ss->size[0] = loop_ub;
            ss->size[1] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)ss, i16, (int)sizeof(double));
            br = iidx->size[0];
            for (i16 = 0; i16 < br; i16++) {
              for (i17 = 0; i17 < loop_ub; i17++) {
                ss->data[i17 + ss->size[0] * i16] = B->data[i17 + B->size[0] *
                  (iidx->data[i16] - 1)];
              }
            }

            if (1 > b_ip) {
              loop_ub = 0;
            } else {
              loop_ub = b_ip;
            }

            br = B->size[0];
            i16 = B->size[0] * B->size[1];
            B->size[0] = br;
            B->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)B, i16, (int)sizeof(double));
            for (i16 = 0; i16 < loop_ub; i16++) {
              for (i17 = 0; i17 < br; i17++) {
                B->data[i17 + B->size[0] * i16] = ss->data[i17 + ss->size[0] *
                  i16];
              }
            }
          }
        }

        if (guard1) {
          i16 = yj->size[0] * yj->size[1];
          yj->size[0] = p;
          yj->size[1] = p;
          emxEnsureCapacity((emxArray__common *)yj, i16, (int)sizeof(double));
          loop_ub = p * p;
          for (i16 = 0; i16 < loop_ub; i16++) {
            yj->data[i16] = 0.0;
          }

          for (j = 0; j < Ifast_idx_0; j++) {
            i16 = xij->size[0] * xij->size[1];
            xij->size[0] = n;
            xij->size[1] = p;
            emxEnsureCapacity((emxArray__common *)xij, i16, (int)sizeof(double));
            for (vstride = 0; vstride < p; vstride++) {
              loop_ub = x->size[0] - 1;
              i16 = e_Ifast->size[0] * e_Ifast->size[1];
              e_Ifast->size[0] = Ifast->size[0];
              e_Ifast->size[1] = Ifast->size[1];
              emxEnsureCapacity((emxArray__common *)e_Ifast, i16, (int)sizeof
                                (int));
              br = Ifast->size[1];
              for (i16 = 0; i16 < br; i16++) {
                nx = Ifast->size[0];
                for (i17 = 0; i17 < nx; i17++) {
                  e_Ifast->data[i17 + e_Ifast->size[0] * i16] = (int)Ifast->
                    data[i17 + Ifast->size[0] * i16];
                }
              }

              pp = x->data[(e_Ifast->data[j] + x->size[0] * vstride) - 1];
              for (i16 = 0; i16 <= loop_ub; i16++) {
                xij->data[i16 + xij->size[0] * vstride] = x->data[i16 + x->size
                  [0] * vstride] - pp;
              }
            }

            i16 = dxij->size[0];
            dxij->size[0] = n;
            emxEnsureCapacity((emxArray__common *)dxij, i16, (int)sizeof(double));
            for (i16 = 0; i16 < n; i16++) {
              dxij->data[i16] = 0.0;
            }

            for (vstride = 0; vstride < m; vstride++) {
              loop_ub = V->size[0];
              i16 = d_Ifast->size[0] * d_Ifast->size[1];
              d_Ifast->size[0] = Ifast->size[0];
              d_Ifast->size[1] = Ifast->size[1];
              emxEnsureCapacity((emxArray__common *)d_Ifast, i16, (int)sizeof
                                (int));
              br = Ifast->size[1];
              for (i16 = 0; i16 < br; i16++) {
                nx = Ifast->size[0];
                for (i17 = 0; i17 < nx; i17++) {
                  d_Ifast->data[i17 + d_Ifast->size[0] * i16] = (int)Ifast->
                    data[i17 + Ifast->size[0] * i16];
                }
              }

              pp = V->data[(d_Ifast->data[j] + V->size[0] * vstride) - 1];
              i16 = e_V->size[0];
              e_V->size[0] = loop_ub;
              emxEnsureCapacity((emxArray__common *)e_V, i16, (int)sizeof(double));
              for (i16 = 0; i16 < loop_ub; i16++) {
                e_V->data[i16] = V->data[i16 + V->size[0] * vstride] - pp;
              }

              b_power(e_V, s);
              i16 = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)dxij, i16, (int)sizeof
                                (double));
              loop_ub = dxij->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                dxij->data[i16] += s->data[i16];
              }
            }

            i16 = s->size[0];
            s->size[0] = dxij->size[0];
            emxEnsureCapacity((emxArray__common *)s, i16, (int)sizeof(double));
            loop_ub = dxij->size[0];
            for (i16 = 0; i16 < loop_ub; i16++) {
              s->data[i16] = dxij->data[i16];
            }

            f_sort(s);
            if ((h2 >= s->data[(int)(2.0 * (double)m) - 1]) || rtIsNaN(s->data
                 [(int)(2.0 * (double)m) - 1])) {
              pp = h2;
            } else {
              pp = s->data[(int)(2.0 * (double)m) - 1];
            }

            i16 = dxij->size[0];
            emxEnsureCapacity((emxArray__common *)dxij, i16, (int)sizeof(double));
            loop_ub = dxij->size[0];
            for (i16 = 0; i16 < loop_ub; i16++) {
              dxij->data[i16] = -dxij->data[i16] / pp;
            }

            b_exp(dxij);
            loop_ub = xij->size[1];
            for (i16 = 0; i16 < loop_ub; i16++) {
              br = xij->size[0];
              for (i17 = 0; i17 < br; i17++) {
                onexi->data[i17 + onexi->size[0] * i16] = xij->data[i17 +
                  xij->size[0] * i16];
              }
            }

            b_repmat(dxij, (double)p + 1.0, a);
            i16 = U->size[0] * U->size[1];
            U->size[0] = onexi->size[1];
            U->size[1] = onexi->size[0];
            emxEnsureCapacity((emxArray__common *)U, i16, (int)sizeof(double));
            loop_ub = onexi->size[0];
            for (i16 = 0; i16 < loop_ub; i16++) {
              br = onexi->size[1];
              for (i17 = 0; i17 < br; i17++) {
                U->data[i17 + U->size[0] * i16] = onexi->data[i16 + onexi->size
                  [0] * i17] * a->data[i16 + a->size[0] * i17];
              }
            }

            if ((U->size[1] == 1) || (onexi->size[0] == 1)) {
              i16 = b_C->size[0] * b_C->size[1];
              b_C->size[0] = U->size[0];
              b_C->size[1] = onexi->size[1];
              emxEnsureCapacity((emxArray__common *)b_C, i16, (int)sizeof(double));
              loop_ub = U->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = onexi->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  b_C->data[i16 + b_C->size[0] * i17] = 0.0;
                  nx = U->size[1];
                  for (b_m = 0; b_m < nx; b_m++) {
                    b_C->data[i16 + b_C->size[0] * i17] += U->data[i16 + U->
                      size[0] * b_m] * onexi->data[b_m + onexi->size[0] * i17];
                  }
                }
              }
            } else {
              k = U->size[1];
              vstride = U->size[0];
              nx = onexi->size[1];
              i16 = b_C->size[0] * b_C->size[1];
              b_C->size[0] = vstride;
              b_C->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)b_C, i16, (int)sizeof(double));
              b_m = U->size[0];
              i16 = b_C->size[0] * b_C->size[1];
              emxEnsureCapacity((emxArray__common *)b_C, i16, (int)sizeof(double));
              loop_ub = b_C->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = b_C->size[0];
                for (i17 = 0; i17 < br; i17++) {
                  b_C->data[i17 + b_C->size[0] * i16] = 0.0;
                }
              }

              if ((U->size[0] == 0) || (onexi->size[1] == 0)) {
              } else {
                vstride = U->size[0] * (onexi->size[1] - 1);
                nx = 0;
                while ((b_m > 0) && (nx <= vstride)) {
                  i16 = nx + b_m;
                  for (ic = nx; ic + 1 <= i16; ic++) {
                    b_C->data[ic] = 0.0;
                  }

                  nx += b_m;
                }

                br = 0;
                nx = 0;
                while ((b_m > 0) && (nx <= vstride)) {
                  ar = 0;
                  i16 = br + k;
                  for (ib = br; ib + 1 <= i16; ib++) {
                    if (onexi->data[ib] != 0.0) {
                      ia = ar;
                      i17 = nx + b_m;
                      for (ic = nx; ic + 1 <= i17; ic++) {
                        ia++;
                        b_C->data[ic] += onexi->data[ib] * U->data[ia - 1];
                      }
                    }

                    ar += b_m;
                  }

                  br += k;
                  nx += b_m;
                }
              }
            }

            eye((double)p + 1.0, a);
            i16 = g_C->size[0] * g_C->size[1];
            g_C->size[0] = b_C->size[0];
            g_C->size[1] = b_C->size[1];
            emxEnsureCapacity((emxArray__common *)g_C, i16, (int)sizeof(double));
            loop_ub = b_C->size[0] * b_C->size[1];
            for (i16 = 0; i16 < loop_ub; i16++) {
              g_C->data[i16] = b_C->data[i16] + a->data[i16] / (double)n;
            }

            inv(g_C, ss);
            if ((U->size[1] == 1) || (ky1->size[0] == 1)) {
              i16 = g_y->size[0] * g_y->size[1];
              g_y->size[0] = U->size[0];
              g_y->size[1] = ky1->size[1];
              emxEnsureCapacity((emxArray__common *)g_y, i16, (int)sizeof(double));
              loop_ub = U->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = ky1->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  g_y->data[i16 + g_y->size[0] * i17] = 0.0;
                  nx = U->size[1];
                  for (b_m = 0; b_m < nx; b_m++) {
                    g_y->data[i16 + g_y->size[0] * i17] += U->data[i16 + U->
                      size[0] * b_m] * ky1->data[b_m + ky1->size[0] * i17];
                  }
                }
              }
            } else {
              k = U->size[1];
              vstride = U->size[0];
              nx = ky1->size[1];
              i16 = g_y->size[0] * g_y->size[1];
              g_y->size[0] = vstride;
              g_y->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)g_y, i16, (int)sizeof(double));
              b_m = U->size[0];
              i16 = g_y->size[0] * g_y->size[1];
              emxEnsureCapacity((emxArray__common *)g_y, i16, (int)sizeof(double));
              loop_ub = g_y->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = g_y->size[0];
                for (i17 = 0; i17 < br; i17++) {
                  g_y->data[i17 + g_y->size[0] * i16] = 0.0;
                }
              }

              if ((U->size[0] == 0) || (ky1->size[1] == 0)) {
              } else {
                vstride = U->size[0] * (ky1->size[1] - 1);
                nx = 0;
                while ((b_m > 0) && (nx <= vstride)) {
                  i16 = nx + b_m;
                  for (ic = nx; ic + 1 <= i16; ic++) {
                    g_y->data[ic] = 0.0;
                  }

                  nx += b_m;
                }

                br = 0;
                nx = 0;
                while ((b_m > 0) && (nx <= vstride)) {
                  ar = 0;
                  i16 = br + k;
                  for (ib = br; ib + 1 <= i16; ib++) {
                    if (ky1->data[ib] != 0.0) {
                      ia = ar;
                      i17 = nx + b_m;
                      for (ic = nx; ic + 1 <= i17; ic++) {
                        ia++;
                        g_y->data[ic] += ky1->data[ib] * U->data[ia - 1];
                      }
                    }

                    ar += b_m;
                  }

                  br += k;
                  nx += b_m;
                }
              }
            }

            if ((ss->size[1] == 1) || (g_y->size[0] == 1)) {
              i16 = abi->size[0] * abi->size[1];
              abi->size[0] = ss->size[0];
              abi->size[1] = g_y->size[1];
              emxEnsureCapacity((emxArray__common *)abi, i16, (int)sizeof(double));
              loop_ub = ss->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = g_y->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  abi->data[i16 + abi->size[0] * i17] = 0.0;
                  nx = ss->size[1];
                  for (b_m = 0; b_m < nx; b_m++) {
                    abi->data[i16 + abi->size[0] * i17] += ss->data[i16 +
                      ss->size[0] * b_m] * g_y->data[b_m + g_y->size[0] * i17];
                  }
                }
              }
            } else {
              k = ss->size[1];
              vstride = ss->size[0];
              nx = g_y->size[1];
              i16 = abi->size[0] * abi->size[1];
              abi->size[0] = vstride;
              abi->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)abi, i16, (int)sizeof(double));
              b_m = ss->size[0];
              i16 = abi->size[0] * abi->size[1];
              emxEnsureCapacity((emxArray__common *)abi, i16, (int)sizeof(double));
              loop_ub = abi->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = abi->size[0];
                for (i17 = 0; i17 < br; i17++) {
                  abi->data[i17 + abi->size[0] * i16] = 0.0;
                }
              }

              if ((ss->size[0] == 0) || (g_y->size[1] == 0)) {
              } else {
                vstride = ss->size[0] * (g_y->size[1] - 1);
                nx = 0;
                while ((b_m > 0) && (nx <= vstride)) {
                  i16 = nx + b_m;
                  for (ic = nx; ic + 1 <= i16; ic++) {
                    abi->data[ic] = 0.0;
                  }

                  nx += b_m;
                }

                br = 0;
                nx = 0;
                while ((b_m > 0) && (nx <= vstride)) {
                  ar = 0;
                  i16 = br + k;
                  for (ib = br; ib + 1 <= i16; ib++) {
                    if (g_y->data[ib] != 0.0) {
                      ia = ar;
                      i17 = nx + b_m;
                      for (ic = nx; ic + 1 <= i17; ic++) {
                        ia++;
                        abi->data[ic] += g_y->data[ib] * ss->data[ia - 1];
                      }
                    }

                    ar += b_m;
                  }

                  br += k;
                  nx += b_m;
                }
              }
            }

            /* *ky2;   %Jan 7 */
            i16 = abi->size[1];
            if ((i16 == 1) || (DD->size[0] == 1)) {
              loop_ub = abi->size[1];
              i16 = c_abi->size[0] * c_abi->size[1];
              c_abi->size[0] = p;
              c_abi->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)c_abi, i16, (int)sizeof
                                (double));
              for (i16 = 0; i16 < loop_ub; i16++) {
                for (i17 = 0; i17 < p; i17++) {
                  c_abi->data[i17 + c_abi->size[0] * i16] = abi->data[i17 +
                    abi->size[0] * i16];
                }
              }

              i16 = h_y->size[0] * h_y->size[1];
              h_y->size[0] = c_abi->size[0];
              h_y->size[1] = DD->size[1];
              emxEnsureCapacity((emxArray__common *)h_y, i16, (int)sizeof(double));
              loop_ub = c_abi->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = DD->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  h_y->data[i16 + h_y->size[0] * i17] = 0.0;
                  nx = c_abi->size[1];
                  for (b_m = 0; b_m < nx; b_m++) {
                    h_y->data[i16 + h_y->size[0] * i17] += c_abi->data[i16 +
                      c_abi->size[0] * b_m] * DD->data[b_m + DD->size[0] * i17];
                  }
                }
              }
            } else {
              i16 = abi->size[1];
              nx = DD->size[1];
              i17 = h_y->size[0] * h_y->size[1];
              h_y->size[0] = p;
              h_y->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)h_y, i17, (int)sizeof(double));
              i17 = h_y->size[0] * h_y->size[1];
              emxEnsureCapacity((emxArray__common *)h_y, i17, (int)sizeof(double));
              loop_ub = h_y->size[1];
              for (i17 = 0; i17 < loop_ub; i17++) {
                br = h_y->size[0];
                for (b_m = 0; b_m < br; b_m++) {
                  h_y->data[b_m + h_y->size[0] * i17] = 0.0;
                }
              }

              if (DD->size[1] != 0) {
                vstride = p * (DD->size[1] - 1);
                for (nx = 0; nx <= vstride; nx += p) {
                  i17 = nx + p;
                  for (ic = nx; ic + 1 <= i17; ic++) {
                    h_y->data[ic] = 0.0;
                  }
                }

                br = 0;
                for (nx = 0; nx <= vstride; nx += p) {
                  ar = 0;
                  i17 = br + i16;
                  for (ib = br; ib + 1 <= i17; ib++) {
                    if (DD->data[ib] != 0.0) {
                      ia = ar;
                      b_m = nx + p;
                      for (ic = nx; ic + 1 <= b_m; ic++) {
                        ia++;
                        h_y->data[ic] += DD->data[ib] * abi->data[(ia - 1) % p +
                          abi->size[0] * ((ia - 1) / p)];
                      }
                    }

                    ar += p;
                  }

                  br += i16;
                }
              }
            }

            loop_ub = abi->size[1];
            i16 = b->size[0] * b->size[1];
            b->size[0] = loop_ub;
            b->size[1] = p;
            emxEnsureCapacity((emxArray__common *)b, i16, (int)sizeof(double));
            for (i16 = 0; i16 < p; i16++) {
              for (i17 = 0; i17 < loop_ub; i17++) {
                b->data[i17 + b->size[0] * i16] = abi->data[i16 + abi->size[0] *
                  i17];
              }
            }

            if ((h_y->size[1] == 1) || (b->size[0] == 1)) {
              i16 = c_C->size[0] * c_C->size[1];
              c_C->size[0] = h_y->size[0];
              c_C->size[1] = b->size[1];
              emxEnsureCapacity((emxArray__common *)c_C, i16, (int)sizeof(double));
              loop_ub = h_y->size[0];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = b->size[1];
                for (i17 = 0; i17 < br; i17++) {
                  c_C->data[i16 + c_C->size[0] * i17] = 0.0;
                  nx = h_y->size[1];
                  for (b_m = 0; b_m < nx; b_m++) {
                    c_C->data[i16 + c_C->size[0] * i17] += h_y->data[i16 +
                      h_y->size[0] * b_m] * b->data[b_m + b->size[0] * i17];
                  }
                }
              }
            } else {
              k = h_y->size[1];
              vstride = h_y->size[0];
              nx = b->size[1];
              i16 = c_C->size[0] * c_C->size[1];
              c_C->size[0] = vstride;
              c_C->size[1] = nx;
              emxEnsureCapacity((emxArray__common *)c_C, i16, (int)sizeof(double));
              b_m = h_y->size[0];
              i16 = c_C->size[0] * c_C->size[1];
              emxEnsureCapacity((emxArray__common *)c_C, i16, (int)sizeof(double));
              loop_ub = c_C->size[1];
              for (i16 = 0; i16 < loop_ub; i16++) {
                br = c_C->size[0];
                for (i17 = 0; i17 < br; i17++) {
                  c_C->data[i17 + c_C->size[0] * i16] = 0.0;
                }
              }

              vstride = h_y->size[0] * (b->size[1] - 1);
              for (nx = 0; nx <= vstride; nx += b_m) {
                i16 = nx + b_m;
                for (ic = nx; ic + 1 <= i16; ic++) {
                  c_C->data[ic] = 0.0;
                }
              }

              br = 0;
              for (nx = 0; nx <= vstride; nx += b_m) {
                ar = 0;
                i16 = br + k;
                for (ib = br; ib + 1 <= i16; ib++) {
                  if (b->data[ib] != 0.0) {
                    ia = ar;
                    i17 = nx + b_m;
                    for (ic = nx; ic + 1 <= i17; ic++) {
                      ia++;
                      c_C->data[ic] += b->data[ib] * h_y->data[ia - 1];
                    }
                  }

                  ar += b_m;
                }

                br += k;
              }
            }

            i16 = yj->size[0] * yj->size[1];
            emxEnsureCapacity((emxArray__common *)yj, i16, (int)sizeof(double));
            vstride = yj->size[0];
            nx = yj->size[1];
            loop_ub = vstride * nx;
            for (i16 = 0; i16 < loop_ub; i16++) {
              yj->data[i16] += c_C->data[i16];
            }
          }

          eig(yj, Vc, Dc);
          i16 = B->size[0] * B->size[1];
          B->size[0] = Vc->size[0];
          B->size[1] = Vc->size[1];
          emxEnsureCapacity((emxArray__common *)B, i16, (int)sizeof(double));
          loop_ub = Vc->size[0] * Vc->size[1];
          for (i16 = 0; i16 < loop_ub; i16++) {
            B->data[i16] = Vc->data[i16].re;
          }

          /* Jan 13 Change for C */
          c_diag(Dc, R);
          g_sort(R, iidx);
          i16 = s->size[0];
          s->size[0] = iidx->size[0];
          emxEnsureCapacity((emxArray__common *)s, i16, (int)sizeof(double));
          loop_ub = iidx->size[0];
          for (i16 = 0; i16 < loop_ub; i16++) {
            s->data[i16] = iidx->data[i16];
          }

          loop_ub = B->size[0];
          i16 = ss->size[0] * ss->size[1];
          ss->size[0] = loop_ub;
          ss->size[1] = iidx->size[0];
          emxEnsureCapacity((emxArray__common *)ss, i16, (int)sizeof(double));
          br = iidx->size[0];
          for (i16 = 0; i16 < br; i16++) {
            for (i17 = 0; i17 < loop_ub; i17++) {
              ss->data[i17 + ss->size[0] * i16] = B->data[i17 + B->size[0] *
                (iidx->data[i16] - 1)];
            }
          }

          if (1 > b_ip) {
            loop_ub = 0;
          } else {
            loop_ub = b_ip;
          }

          br = B->size[0];
          i16 = B->size[0] * B->size[1];
          B->size[0] = br;
          B->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)B, i16, (int)sizeof(double));
          for (i16 = 0; i16 < loop_ub; i16++) {
            for (i17 = 0; i17 < br; i17++) {
              B->data[i17 + B->size[0] * i16] = ss->data[i17 + ss->size[0] * i16];
            }
          }
        }
      }

      loop_ub = B->size[1];
      for (i16 = 0; i16 < loop_ub; i16++) {
        br = B->size[0];
        for (i17 = 0; i17 < br; i17++) {
          BB->data[(i17 + BB->size[0] * i16) + BB->size[0] * BB->size[1] * (b_ip
            - 1)] = B->data[i17 + B->size[0] * i16];
        }
      }

      if (d_strcmp(method)) {
        i15 = (int)((1.0 + (-1.0 - ((double)p - 1.0))) / -1.0);
        for (b_m = 0; b_m < i15; b_m++) {
          vstride = (p - b_m) - 1;
          if (1 > vstride) {
            loop_ub = -1;
          } else {
            loop_ub = vstride - 1;
          }

          br = ss->size[0] - 1;
          for (i16 = 0; i16 <= loop_ub; i16++) {
            for (i17 = 0; i17 <= br; i17++) {
              BB->data[(i17 + BB->size[0] * i16) + BB->size[0] * BB->size[1] *
                (vstride - 1)] = ss->data[i17 + ss->size[0] * i16];
            }
          }
        }

        exitg1 = true;
      } else {
        if (1 > b_ip - 1) {
          loop_ub = 0;
        } else {
          loop_ub = b_ip - 1;
        }

        vstride = B->size[0];
        i16 = c_B->size[0] * c_B->size[1];
        c_B->size[0] = vstride;
        c_B->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)c_B, i16, (int)sizeof(double));
        for (i16 = 0; i16 < loop_ub; i16++) {
          for (i17 = 0; i17 < vstride; i17++) {
            c_B->data[i17 + c_B->size[0] * i16] = B->data[i17 + B->size[0] * i16];
          }
        }

        i16 = B->size[0] * B->size[1];
        B->size[0] = c_B->size[0];
        B->size[1] = c_B->size[1];
        emxEnsureCapacity((emxArray__common *)B, i16, (int)sizeof(double));
        loop_ub = c_B->size[1];
        for (i16 = 0; i16 < loop_ub; i16++) {
          br = c_B->size[0];
          for (i17 = 0; i17 < br; i17++) {
            B->data[i17 + B->size[0] * i16] = c_B->data[i17 + c_B->size[0] * i16];
          }
        }

        ip++;
      }
    }

    emxFree_int32_T(&e_Ifast);
    emxFree_int32_T(&d_Ifast);
    emxFree_real_T(&c_abi);
    emxFree_real_T(&b_dc);
    emxFree_int32_T(&c_Ifast);
    emxFree_int32_T(&b_Ifast);
    emxFree_real_T(&c_B);
    emxFree_real_T(&e_V);
    emxFree_real_T(&g_C);
    emxFree_real_T(&b_xfast);
    emxFree_real_T(&b_Vfast);
    emxFree_real_T(&d_V);
    emxFree_real_T(&f_C);
    emxFree_real_T(&b_abi);
    emxFree_real_T(&b_dd);
    emxFree_real_T(&h_y);
    emxFree_real_T(&g_y);
    emxFree_real_T(&f_y);
    emxFree_real_T(&e_y);
    emxFreeMatrix_cell_wrap_0(reshapes);
    emxFree_real_T(&d_y);
    emxFree_real_T(&c_y);
    emxFree_real_T(&b_y);
    emxFree_int32_T(&r2);
    emxFree_real_T(&e_C);
    emxFree_real_T(&d_C);
    emxFree_real_T(&c_C);
    emxFree_real_T(&b_C);
    emxFree_creal_T(&R);
    emxFree_real_T(&b_B);
    emxFree_real_T(&tmp);
    emxFree_real_T(&ddx);
    emxFree_real_T(&kxijy);
    emxFree_real_T(&dd);
    emxFree_real_T(&abi);
    emxFree_real_T(&dxij);
    emxFree_real_T(&xij);
    emxFree_real_T(&xfast);
    emxFree_real_T(&Vfast);
    emxFree_real_T(&Ifast);
    emxFree_real_T(&D);
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&qx);
  emxFree_creal_T(&Dc);
  emxFree_creal_T(&Vc);
  emxFree_real_T(&dc);
  emxFree_real_T(&B);
  emxFree_real_T(&onexi);
  emxFree_real_T(&s);
  emxFree_real_T(&U);
  emxFree_real_T(&C);
  emxFree_real_T(&yj);
  emxFree_real_T(&DD);
  emxFree_real_T(&ky2);
  emxFree_real_T(&ky1);
  emxFree_real_T(&V);
  emxFree_real_T(&ss);
  nx = BB->size[0] * BB->size[1] * BB->size[2];
  i15 = BB1D->size[0];
  BB1D->size[0] = (int)rt_powd_snf((double)p, 3.0);
  emxEnsureCapacity((emxArray__common *)BB1D, i15, (int)sizeof(double));
  for (k = 0; k + 1 <= nx; k++) {
    BB1D->data[k] = BB->data[k];
  }

  emxFree_real_T(&BB);
}

/* End of code generation (MAVEfast.cpp) */
