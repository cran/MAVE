/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * prctile.cpp
 *
 * Code generation for function 'prctile'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "prctile.h"
#include "MAVEfast_emxutil.h"

/* Function Declarations */
static void percentile_array(const emxArray_real_T *x, const emxArray_real_T *p,
  emxArray_real_T *pct);
static double rt_roundd_snf(double u);

/* Function Definitions */
static void percentile_array(const emxArray_real_T *x, const emxArray_real_T *p,
  emxArray_real_T *pct)
{
  unsigned int sz_idx_0;
  int j;
  emxArray_real_T *wk;
  emxArray_real_T *pctv;
  emxArray_int32_T *idx;
  emxArray_int32_T *iwork;
  int pEnd;
  int ix;
  int iy;
  int xi;
  int ixstart;
  int nj;
  int iystart;
  int k;
  int n;
  boolean_T b_p;
  int c_p;
  int q;
  int qEnd;
  double r;
  double i;
  int kEnd;
  sz_idx_0 = (unsigned int)p->size[1];
  j = pct->size[0] * pct->size[1];
  pct->size[0] = (int)sz_idx_0;
  pct->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)pct, j, (int)sizeof(double));
  emxInit_real_T1(&wk, 1);
  emxInit_real_T(&pctv, 2);
  emxInit_int32_T1(&idx, 1);
  emxInit_int32_T1(&iwork, 1);
  if ((x->size[0] == 0) || (x->size[1] == 0) || ((pct->size[0] == 0) ||
       (pct->size[1] == 0))) {
    j = pct->size[0] * pct->size[1];
    emxEnsureCapacity((emxArray__common *)pct, j, (int)sizeof(double));
    pEnd = pct->size[1];
    for (j = 0; j < pEnd; j++) {
      ixstart = pct->size[0];
      for (nj = 0; nj < ixstart; nj++) {
        pct->data[nj + pct->size[0] * j] = rtNaN;
      }
    }
  } else {
    j = wk->size[0];
    wk->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)wk, j, (int)sizeof(double));
    ix = -1;
    iy = -1;
    for (xi = 1; xi <= x->size[1]; xi++) {
      ixstart = ix + 1;
      iystart = iy + 1;
      ix++;
      wk->data[0] = x->data[ixstart];
      for (k = 2; k <= x->size[0]; k++) {
        ix++;
        wk->data[k - 1] = x->data[ix];
      }

      iy++;
      j = pctv->size[0] * pctv->size[1];
      pctv->size[0] = 1;
      pctv->size[1] = p->size[1];
      emxEnsureCapacity((emxArray__common *)pctv, j, (int)sizeof(double));
      if (p->size[1] == 0) {
        j = pctv->size[0] * pctv->size[1];
        pctv->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)pctv, j, (int)sizeof(double));
        pEnd = pctv->size[1];
        for (j = 0; j < pEnd; j++) {
          pctv->data[pctv->size[0] * j] = rtNaN;
        }
      } else {
        n = wk->size[0] + 1;
        sz_idx_0 = (unsigned int)wk->size[0];
        j = idx->size[0];
        idx->size[0] = (int)sz_idx_0;
        emxEnsureCapacity((emxArray__common *)idx, j, (int)sizeof(int));
        pEnd = (int)sz_idx_0;
        for (j = 0; j < pEnd; j++) {
          idx->data[j] = 0;
        }

        j = iwork->size[0];
        iwork->size[0] = (int)sz_idx_0;
        emxEnsureCapacity((emxArray__common *)iwork, j, (int)sizeof(int));
        for (k = 1; k <= n - 2; k += 2) {
          if ((wk->data[k - 1] <= wk->data[k]) || rtIsNaN(wk->data[k])) {
            b_p = true;
          } else {
            b_p = false;
          }

          if (b_p) {
            idx->data[k - 1] = k;
            idx->data[k] = k + 1;
          } else {
            idx->data[k - 1] = k + 1;
            idx->data[k] = k;
          }
        }

        if ((wk->size[0] & 1) != 0) {
          idx->data[wk->size[0] - 1] = wk->size[0];
        }

        ixstart = 2;
        while (ixstart < n - 1) {
          nj = ixstart << 1;
          j = 1;
          for (pEnd = 1 + ixstart; pEnd < n; pEnd = qEnd + ixstart) {
            c_p = j;
            q = pEnd - 1;
            qEnd = j + nj;
            if (qEnd > n) {
              qEnd = n;
            }

            k = 0;
            kEnd = qEnd - j;
            while (k + 1 <= kEnd) {
              if ((wk->data[idx->data[c_p - 1] - 1] <= wk->data[idx->data[q] - 1])
                  || rtIsNaN(wk->data[idx->data[q] - 1])) {
                b_p = true;
              } else {
                b_p = false;
              }

              if (b_p) {
                iwork->data[k] = idx->data[c_p - 1];
                c_p++;
                if (c_p == pEnd) {
                  while (q + 1 < qEnd) {
                    k++;
                    iwork->data[k] = idx->data[q];
                    q++;
                  }
                }
              } else {
                iwork->data[k] = idx->data[q];
                q++;
                if (q + 1 == qEnd) {
                  while (c_p < pEnd) {
                    k++;
                    iwork->data[k] = idx->data[c_p - 1];
                    c_p++;
                  }
                }
              }

              k++;
            }

            for (k = 0; k + 1 <= kEnd; k++) {
              idx->data[(j + k) - 1] = iwork->data[k];
            }

            j = qEnd;
          }

          ixstart = nj;
        }

        nj = wk->size[0];
        while ((nj > 0) && rtIsNaN(wk->data[idx->data[nj - 1] - 1])) {
          nj--;
        }

        if (nj < 1) {
          j = pctv->size[0] * pctv->size[1];
          pctv->size[0] = 1;
          emxEnsureCapacity((emxArray__common *)pctv, j, (int)sizeof(double));
          pEnd = pctv->size[1];
          for (j = 0; j < pEnd; j++) {
            pctv->data[pctv->size[0] * j] = rtNaN;
          }
        } else if (nj == 1) {
          j = pctv->size[0] * pctv->size[1];
          pctv->size[0] = 1;
          emxEnsureCapacity((emxArray__common *)pctv, j, (int)sizeof(double));
          pEnd = pctv->size[1];
          for (j = 0; j < pEnd; j++) {
            pctv->data[pctv->size[0] * j] = wk->data[idx->data[0] - 1];
          }
        } else {
          for (k = 0; k < p->size[1]; k++) {
            r = p->data[k] / 100.0 * (double)nj;
            i = rt_roundd_snf(r);
            if (i < 1.0) {
              pctv->data[k] = wk->data[idx->data[0] - 1];
            } else if (nj <= i) {
              pctv->data[k] = wk->data[idx->data[nj - 1] - 1];
            } else {
              r -= i;
              pctv->data[k] = (0.5 - r) * wk->data[idx->data[(int)i - 1] - 1] +
                (0.5 + r) * wk->data[idx->data[(int)(i + 1.0) - 1] - 1];
            }
          }
        }
      }

      pct->data[iystart] = pctv->data[0];
      for (k = 2; k <= p->size[1]; k++) {
        iy++;
        pct->data[iy] = pctv->data[k - 1];
      }
    }
  }

  emxFree_int32_T(&iwork);
  emxFree_int32_T(&idx);
  emxFree_real_T(&pctv);
  emxFree_real_T(&wk);
}

static double rt_roundd_snf(double u)
{
  double y;
  if (std::abs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

void prctile(const emxArray_real_T *x, const emxArray_real_T *p, emxArray_real_T
             *y)
{
  percentile_array(x, p, y);
}

/* End of code generation (prctile.cpp) */
