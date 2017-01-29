/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * CVfast.cpp
 *
 * Code generation for function 'CVfast'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "sum.h"
#include "power.h"
#include "CVfast_emxutil.h"
#include "repmat.h"
#include "mean.h"
#include "rdivide.h"
#include "combine_vector_elements.h"
#include "sort1.h"
#include "rand.h"
#include "std.h"
#include "CVfast_rtwutil.h"
#include "CVfast_data.h"

/* Function Definitions */
void CVfast(const emxArray_real_T *x, const emxArray_real_T *ky, const
            emxArray_real_T *BB1D, emxArray_real_T *cv)
{
  int n;
  int nD;
  int i0;
  int ar;
  int p;
  emxArray_real_T *nBB;
  emxArray_real_T *nx;
  emxArray_real_T *dx;
  emxArray_int32_T *Ia;
  emxArray_int32_T *Ib;
  emxArray_real_T *yA;
  emxArray_real_T *dxj;
  emxArray_real_T *sd;
  emxArray_real_T *hj;
  emxArray_real_T *ye;
  emxArray_real_T *unusedU1;
  emxArray_int32_T *iidx;
  emxArray_real_T *z;
  emxArray_real_T *y;
  emxArray_real_T *b_x;
  emxArray_real_T *r0;
  emxArray_real_T *b_nx;
  emxArray_real_T *c_nx;
  emxArray_real_T *b_z;
  emxArray_int32_T *c_z;
  double h;
  double b_y;
  int d_nx;
  int k;
  unsigned int uv0[2];
  int m;
  int br;
  int i1;
  int cr;
  int na;
  int ic;
  int j;
  int ib;
  unsigned int r;
  int ia;
  double extremum;

  /*  NEW function to select the points where gradients are calculated */
  /*  3-fold CV with random splitting of 100 times */
  n = x->size[0];
  nD = x->size[1];
  i0 = cv->size[0] * cv->size[1];
  cv->size[0] = 1;
  cv->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)cv, i0, (int)sizeof(double));
  ar = x->size[1];
  for (i0 = 0; i0 < ar; i0++) {
    cv->data[i0] = 0.0;
  }

  p = 0;
  emxInit_real_T(&nBB, 2);
  emxInit_real_T(&nx, 2);
  emxInit_real_T(&dx, 2);
  emxInit_int32_T(&Ia, 1);
  emxInit_int32_T(&Ib, 1);
  emxInit_real_T(&yA, 2);
  emxInit_real_T(&dxj, 2);
  emxInit_real_T(&sd, 2);
  emxInit_real_T(&hj, 2);
  emxInit_real_T(&ye, 2);
  emxInit_real_T1(&unusedU1, 1);
  emxInit_int32_T(&iidx, 1);
  emxInit_real_T(&z, 2);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&b_x, 2);
  emxInit_real_T1(&r0, 1);
  emxInit_real_T(&b_nx, 2);
  emxInit_real_T(&c_nx, 2);
  emxInit_real_T1(&b_z, 1);
  emxInit_int32_T(&c_z, 1);
  while (p <= nD - 1) {
    /* for matlab coder */
    h = (double)nD * (double)nD * ((1.0 + (double)p) - 1.0);
    b_y = (1.0 + (double)p) * (double)nD;
    i0 = z->size[0] * z->size[1];
    z->size[0] = 1;
    z->size[1] = (int)(b_y - 1.0) + 1;
    emxEnsureCapacity((emxArray__common *)z, i0, (int)sizeof(double));
    ar = (int)(b_y - 1.0);
    for (i0 = 0; i0 <= ar; i0++) {
      z->data[z->size[0] * i0] = h + (1.0 + (double)i0);
    }

    /* for matlab coder */
    i0 = c_z->size[0];
    c_z->size[0] = z->size[1];
    emxEnsureCapacity((emxArray__common *)c_z, i0, (int)sizeof(int));
    ar = z->size[1];
    for (i0 = 0; i0 < ar; i0++) {
      c_z->data[i0] = (int)z->data[z->size[0] * i0];
    }

    d_nx = c_z->size[0];
    i0 = nBB->size[0] * nBB->size[1];
    nBB->size[0] = nD;
    nBB->size[1] = p + 1;
    emxEnsureCapacity((emxArray__common *)nBB, i0, (int)sizeof(double));
    for (k = 0; k + 1 <= d_nx; k++) {
      nBB->data[k] = BB1D->data[(int)z->data[z->size[0] * k] - 1];
    }

    /* for matlab coder */
    if ((x->size[1] == 1) || (nBB->size[0] == 1)) {
      i0 = nx->size[0] * nx->size[1];
      nx->size[0] = x->size[0];
      nx->size[1] = nBB->size[1];
      emxEnsureCapacity((emxArray__common *)nx, i0, (int)sizeof(double));
      ar = x->size[0];
      for (i0 = 0; i0 < ar; i0++) {
        br = nBB->size[1];
        for (i1 = 0; i1 < br; i1++) {
          nx->data[i0 + nx->size[0] * i1] = 0.0;
          cr = x->size[1];
          for (d_nx = 0; d_nx < cr; d_nx++) {
            nx->data[i0 + nx->size[0] * i1] += x->data[i0 + x->size[0] * d_nx] *
              nBB->data[d_nx + nBB->size[0] * i1];
          }
        }
      }
    } else {
      k = x->size[1];
      uv0[0] = (unsigned int)x->size[0];
      uv0[1] = (unsigned int)nBB->size[1];
      i0 = nx->size[0] * nx->size[1];
      nx->size[0] = (int)uv0[0];
      nx->size[1] = (int)uv0[1];
      emxEnsureCapacity((emxArray__common *)nx, i0, (int)sizeof(double));
      m = x->size[0];
      i0 = nx->size[0] * nx->size[1];
      emxEnsureCapacity((emxArray__common *)nx, i0, (int)sizeof(double));
      ar = nx->size[1];
      for (i0 = 0; i0 < ar; i0++) {
        br = nx->size[0];
        for (i1 = 0; i1 < br; i1++) {
          nx->data[i1 + nx->size[0] * i0] = 0.0;
        }
      }

      if (x->size[0] != 0) {
        d_nx = x->size[0] * (nBB->size[1] - 1);
        cr = 0;
        while ((m > 0) && (cr <= d_nx)) {
          i0 = cr + m;
          for (ic = cr; ic + 1 <= i0; ic++) {
            nx->data[ic] = 0.0;
          }

          cr += m;
        }

        br = 0;
        cr = 0;
        while ((m > 0) && (cr <= d_nx)) {
          ar = 0;
          i0 = br + k;
          for (ib = br; ib + 1 <= i0; ib++) {
            if (nBB->data[ib] != 0.0) {
              ia = ar;
              i1 = cr + m;
              for (ic = cr; ic + 1 <= i1; ic++) {
                ia++;
                nx->data[ic] += nBB->data[ib] * x->data[ia - 1];
              }
            }

            ar += m;
          }

          br += k;
          cr += m;
        }
      }
    }

    /* for matlab coder */
    b_std(nx, z);
    h = mean(z);
    h /= rt_powd_snf((double)n, 1.0 / ((1.0 + (double)p) + 3.0));
    h *= 2.0 * h;
    na = (int)std::floor((double)n * 2.0 / 3.0);
    i0 = dx->size[0] * dx->size[1];
    dx->size[0] = n;
    dx->size[1] = n;
    emxEnsureCapacity((emxArray__common *)dx, i0, (int)sizeof(double));
    for (d_nx = 0; d_nx < n; d_nx++) {
      ar = nx->size[1];
      i0 = c_nx->size[0] * c_nx->size[1];
      c_nx->size[0] = 1;
      c_nx->size[1] = ar;
      emxEnsureCapacity((emxArray__common *)c_nx, i0, (int)sizeof(double));
      for (i0 = 0; i0 < ar; i0++) {
        c_nx->data[c_nx->size[0] * i0] = nx->data[d_nx + nx->size[0] * i0];
      }

      repmat(c_nx, (double)n, nBB);
      i0 = b_nx->size[0] * b_nx->size[1];
      b_nx->size[0] = nx->size[0];
      b_nx->size[1] = nx->size[1];
      emxEnsureCapacity((emxArray__common *)b_nx, i0, (int)sizeof(double));
      ar = nx->size[0] * nx->size[1];
      for (i0 = 0; i0 < ar; i0++) {
        b_nx->data[i0] = nx->data[i0] - nBB->data[i0];
      }

      power(b_nx, nBB);
      sum(nBB, r0);
      ar = r0->size[0];
      for (i0 = 0; i0 < ar; i0++) {
        dx->data[i0 + dx->size[0] * d_nx] = r0->data[i0];
      }
    }

    for (j = 0; j < 100; j++) {
      /* rand('state', j); not support; */
      r = (unsigned int)(1 + j);
      state[0] = (unsigned int)(1 + j);
      for (d_nx = 0; d_nx < 623; d_nx++) {
        r = (r ^ r >> 30U) * 1812433253U + (1 + d_nx);
        state[d_nx + 1] = r;
      }

      state[624] = 624U;
      b_rand((double)n, unusedU1);
      sort(unusedU1, iidx);
      i0 = unusedU1->size[0];
      unusedU1->size[0] = iidx->size[0];
      emxEnsureCapacity((emxArray__common *)unusedU1, i0, (int)sizeof(double));
      ar = iidx->size[0];
      for (i0 = 0; i0 < ar; i0++) {
        unusedU1->data[i0] = iidx->data[i0];
      }

      if (1 > na) {
        ar = 0;
      } else {
        ar = na;
      }

      i0 = Ia->size[0];
      Ia->size[0] = ar;
      emxEnsureCapacity((emxArray__common *)Ia, i0, (int)sizeof(int));
      for (i0 = 0; i0 < ar; i0++) {
        Ia->data[i0] = (int)unusedU1->data[i0];
      }

      if (na + 1 > n) {
        i0 = 1;
        i1 = 0;
      } else {
        i0 = na + 1;
        i1 = n;
      }

      d_nx = Ib->size[0];
      Ib->size[0] = (i1 - i0) + 1;
      emxEnsureCapacity((emxArray__common *)Ib, d_nx, (int)sizeof(int));
      br = (i1 - i0) + 1;
      for (i1 = 0; i1 < br; i1++) {
        Ib->data[i1] = (int)unusedU1->data[(i0 + i1) - 1];
      }

      br = ky->size[1];
      i0 = yA->size[0] * yA->size[1];
      yA->size[0] = Ia->size[0];
      yA->size[1] = br;
      emxEnsureCapacity((emxArray__common *)yA, i0, (int)sizeof(double));
      for (i0 = 0; i0 < br; i0++) {
        cr = Ia->size[0];
        for (i1 = 0; i1 < cr; i1++) {
          yA->data[i1 + yA->size[0] * i0] = ky->data[(Ia->data[i1] + ky->size[0]
            * i0) - 1];
        }
      }

      i0 = dxj->size[0] * dxj->size[1];
      dxj->size[0] = Ia->size[0];
      dxj->size[1] = Ib->size[0];
      emxEnsureCapacity((emxArray__common *)dxj, i0, (int)sizeof(double));
      br = Ib->size[0];
      for (i0 = 0; i0 < br; i0++) {
        cr = Ia->size[0];
        for (i1 = 0; i1 < cr; i1++) {
          dxj->data[i1 + dxj->size[0] * i0] = dx->data[(Ia->data[i1] + dx->size
            [0] * (Ib->data[i0] - 1)) - 1];
        }
      }

      i0 = sd->size[0] * sd->size[1];
      sd->size[0] = dxj->size[0];
      sd->size[1] = dxj->size[1];
      emxEnsureCapacity((emxArray__common *)sd, i0, (int)sizeof(double));
      br = dxj->size[0] * dxj->size[1];
      for (i0 = 0; i0 < br; i0++) {
        sd->data[i0] = dxj->data[i0];
      }

      b_sort(sd);
      b_y = 1.0 / h;
      br = sd->size[1];
      i0 = z->size[0] * z->size[1];
      z->size[0] = 1;
      z->size[1] = br;
      emxEnsureCapacity((emxArray__common *)z, i0, (int)sizeof(double));
      for (i0 = 0; i0 < br; i0++) {
        z->data[z->size[0] * i0] = 1.0 / sd->data[1 + sd->size[0] * i0];
      }

      d_nx = z->size[1];
      i0 = hj->size[0] * hj->size[1];
      hj->size[0] = 1;
      hj->size[1] = z->size[1];
      emxEnsureCapacity((emxArray__common *)hj, i0, (int)sizeof(double));
      for (k = 0; k + 1 <= d_nx; k++) {
        if ((b_y <= z->data[k]) || rtIsNaN(z->data[k])) {
          extremum = b_y;
        } else {
          extremum = z->data[k];
        }

        hj->data[k] = extremum;
      }

      repmat(hj, (double)na, nBB);
      i0 = sd->size[0] * sd->size[1];
      sd->size[0] = nBB->size[0];
      sd->size[1] = nBB->size[1];
      emxEnsureCapacity((emxArray__common *)sd, i0, (int)sizeof(double));
      br = nBB->size[0] * nBB->size[1];
      for (i0 = 0; i0 < br; i0++) {
        sd->data[i0] = nBB->data[i0];
      }

      i0 = dxj->size[0] * dxj->size[1];
      emxEnsureCapacity((emxArray__common *)dxj, i0, (int)sizeof(double));
      d_nx = dxj->size[0];
      cr = dxj->size[1];
      br = d_nx * cr;
      for (i0 = 0; i0 < br; i0++) {
        dxj->data[i0] *= sd->data[i0];
      }

      /* ker start */
      i0 = dxj->size[0] * dxj->size[1];
      emxEnsureCapacity((emxArray__common *)dxj, i0, (int)sizeof(double));
      d_nx = dxj->size[0];
      cr = dxj->size[1];
      br = d_nx * cr;
      for (i0 = 0; i0 < br; i0++) {
        dxj->data[i0] = (1.0 - dxj->data[i0]) * (double)(1.0 > dxj->data[i0]);
      }

      /* ye start */
      i0 = nBB->size[0] * nBB->size[1];
      nBB->size[0] = dxj->size[1];
      nBB->size[1] = dxj->size[0];
      emxEnsureCapacity((emxArray__common *)nBB, i0, (int)sizeof(double));
      br = dxj->size[0];
      for (i0 = 0; i0 < br; i0++) {
        cr = dxj->size[1];
        for (i1 = 0; i1 < cr; i1++) {
          nBB->data[i1 + nBB->size[0] * i0] = dxj->data[i0 + dxj->size[0] * i1];
        }
      }

      br = ky->size[1];
      i0 = ye->size[0] * ye->size[1];
      ye->size[0] = Ia->size[0];
      ye->size[1] = br;
      emxEnsureCapacity((emxArray__common *)ye, i0, (int)sizeof(double));
      for (i0 = 0; i0 < br; i0++) {
        cr = Ia->size[0];
        for (i1 = 0; i1 < cr; i1++) {
          ye->data[i1 + ye->size[0] * i0] = ky->data[(Ia->data[i1] + ky->size[0]
            * i0) - 1];
        }
      }

      if ((nBB->size[1] == 1) || (ar == 1)) {
        i0 = y->size[0] * y->size[1];
        y->size[0] = nBB->size[0];
        y->size[1] = yA->size[1];
        emxEnsureCapacity((emxArray__common *)y, i0, (int)sizeof(double));
        ar = nBB->size[0];
        for (i0 = 0; i0 < ar; i0++) {
          br = yA->size[1];
          for (i1 = 0; i1 < br; i1++) {
            y->data[i0 + y->size[0] * i1] = 0.0;
            cr = nBB->size[1];
            for (d_nx = 0; d_nx < cr; d_nx++) {
              y->data[i0 + y->size[0] * i1] += nBB->data[i0 + nBB->size[0] *
                d_nx] * yA->data[d_nx + yA->size[0] * i1];
            }
          }
        }
      } else {
        k = nBB->size[1];
        i0 = ky->size[1];
        uv0[0] = (unsigned int)nBB->size[0];
        i1 = y->size[0] * y->size[1];
        y->size[0] = (int)uv0[0];
        y->size[1] = i0;
        emxEnsureCapacity((emxArray__common *)y, i1, (int)sizeof(double));
        m = nBB->size[0];
        i0 = y->size[0] * y->size[1];
        emxEnsureCapacity((emxArray__common *)y, i0, (int)sizeof(double));
        ar = y->size[1];
        for (i0 = 0; i0 < ar; i0++) {
          br = y->size[0];
          for (i1 = 0; i1 < br; i1++) {
            y->data[i1 + y->size[0] * i0] = 0.0;
          }
        }

        if (nBB->size[0] == 0) {
        } else {
          i0 = ky->size[1];
          if (i0 == 0) {
          } else {
            i0 = ky->size[1] - 1;
            d_nx = nBB->size[0] * i0;
            cr = 0;
            while ((m > 0) && (cr <= d_nx)) {
              i0 = cr + m;
              for (ic = cr; ic + 1 <= i0; ic++) {
                y->data[ic] = 0.0;
              }

              cr += m;
            }

            br = 0;
            cr = 0;
            while ((m > 0) && (cr <= d_nx)) {
              ar = 0;
              i0 = br + k;
              for (ib = br; ib + 1 <= i0; ib++) {
                if (yA->data[ib] != 0.0) {
                  ia = ar;
                  i1 = cr + m;
                  for (ic = cr; ic + 1 <= i1; ic++) {
                    ia++;
                    y->data[ic] += ye->data[ib] * nBB->data[ia - 1];
                  }
                }

                ar += m;
              }

              br += k;
              cr += m;
            }
          }
        }
      }

      b_combine_vector_elements(dxj, z);
      i0 = b_z->size[0];
      b_z->size[0] = z->size[1];
      emxEnsureCapacity((emxArray__common *)b_z, i0, (int)sizeof(double));
      ar = z->size[1];
      for (i0 = 0; i0 < ar; i0++) {
        b_z->data[i0] = z->data[z->size[0] * i0] + 1.0E-6;
      }

      b_repmat(b_z, (double)ky->size[1], nBB);
      rdivide(y, nBB, ye);

      /* cv */
      ar = ky->size[1];
      i0 = b_x->size[0] * b_x->size[1];
      b_x->size[0] = Ib->size[0];
      b_x->size[1] = ar;
      emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(double));
      for (i0 = 0; i0 < ar; i0++) {
        br = Ib->size[0];
        for (i1 = 0; i1 < br; i1++) {
          b_x->data[i1 + b_x->size[0] * i0] = ky->data[(Ib->data[i1] + ky->size
            [0] * i0) - 1] - ye->data[i1 + ye->size[0] * i0];
        }
      }

      for (i0 = 0; i0 < 2; i0++) {
        uv0[i0] = (unsigned int)b_x->size[i0];
      }

      i0 = nBB->size[0] * nBB->size[1];
      nBB->size[0] = (int)uv0[0];
      nBB->size[1] = (int)uv0[1];
      emxEnsureCapacity((emxArray__common *)nBB, i0, (int)sizeof(double));
      d_nx = b_x->size[0] * b_x->size[1];
      for (k = 0; k + 1 <= d_nx; k++) {
        nBB->data[k] = std::abs(b_x->data[k]);
      }

      b_mean(nBB, z);
      cv->data[p] += mean(z);
    }

    p++;
  }

  emxFree_int32_T(&c_z);
  emxFree_real_T(&b_z);
  emxFree_real_T(&c_nx);
  emxFree_real_T(&b_nx);
  emxFree_real_T(&r0);
  emxFree_real_T(&b_x);
  emxFree_real_T(&y);
  emxFree_real_T(&z);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&unusedU1);
  emxFree_real_T(&ye);
  emxFree_real_T(&hj);
  emxFree_real_T(&sd);
  emxFree_real_T(&dxj);
  emxFree_real_T(&yA);
  emxFree_int32_T(&Ib);
  emxFree_int32_T(&Ia);
  emxFree_real_T(&dx);
  emxFree_real_T(&nx);
  emxFree_real_T(&nBB);
}

/* End of code generation (CVfast.cpp) */
