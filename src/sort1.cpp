/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sort1.cpp
 *
 * Code generation for function 'sort1'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "sort1.h"
#include "sortIdx.h"
#include "CVfast_emxutil.h"
#include "relop.h"

/* Function Declarations */
static void c_sort(emxArray_real_T *x, int dim);
static void e_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx);
static void h_sort(emxArray_creal_T *x, int dim, emxArray_int32_T *idx);

/* Function Definitions */
static void c_sort(emxArray_real_T *x, int dim)
{
  emxArray_real_T *vwork;
  int i14;
  int vstride;
  int k;
  int npages;
  int pagesize;
  int i;
  emxArray_int32_T *b_vwork;
  int pageoffset;
  int j;
  int idx0;
  emxInit_real_T1(&vwork, 1);
  i14 = x->size[dim - 1];
  vstride = x->size[dim - 1];
  k = vwork->size[0];
  vwork->size[0] = vstride;
  emxEnsureCapacity((emxArray__common *)vwork, k, (int)sizeof(double));
  vstride = 1;
  k = 1;
  while (k <= dim - 1) {
    vstride *= x->size[0];
    k = 2;
  }

  npages = 1;
  k = dim + 1;
  while (k < 3) {
    npages *= x->size[1];
    k = 3;
  }

  pagesize = x->size[dim - 1] * vstride;
  i = 1;
  emxInit_int32_T(&b_vwork, 1);
  while (i <= npages) {
    pageoffset = (i - 1) * pagesize;
    for (j = 0; j + 1 <= vstride; j++) {
      idx0 = pageoffset + j;
      for (k = 0; k + 1 <= i14; k++) {
        vwork->data[k] = x->data[idx0 + k * vstride];
      }

      sortIdx(vwork, b_vwork);
      for (k = 0; k + 1 <= i14; k++) {
        x->data[idx0 + k * vstride] = vwork->data[k];
      }
    }

    i++;
  }

  emxFree_int32_T(&b_vwork);
  emxFree_real_T(&vwork);
}

static void e_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx)
{
  int i22;
  emxArray_real_T *vwork;
  int vstride;
  int x_idx_0;
  int j;
  emxArray_int32_T *iidx;
  if (dim <= 1) {
    i22 = x->size[0];
  } else {
    i22 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  vstride = vwork->size[0];
  vwork->size[0] = i22;
  emxEnsureCapacity((emxArray__common *)vwork, vstride, (int)sizeof(double));
  x_idx_0 = x->size[0];
  vstride = idx->size[0];
  idx->size[0] = x_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, vstride, (int)sizeof(int));
  vstride = 1;
  x_idx_0 = 1;
  while (x_idx_0 <= dim - 1) {
    vstride *= x->size[0];
    x_idx_0 = 2;
  }

  j = 0;
  emxInit_int32_T(&iidx, 1);
  while (j + 1 <= vstride) {
    for (x_idx_0 = 0; x_idx_0 + 1 <= i22; x_idx_0++) {
      vwork->data[x_idx_0] = x->data[j + x_idx_0 * vstride];
    }

    b_sortIdx(vwork, iidx);
    for (x_idx_0 = 0; x_idx_0 + 1 <= i22; x_idx_0++) {
      x->data[j + x_idx_0 * vstride] = vwork->data[x_idx_0];
      idx->data[j + x_idx_0 * vstride] = iidx->data[x_idx_0];
    }

    j++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

static void h_sort(emxArray_creal_T *x, int dim, emxArray_int32_T *idx)
{
  int i24;
  emxArray_creal_T *vwork;
  int i;
  int x_idx_0;
  int vstride;
  int j;
  emxArray_int32_T *b_idx;
  emxArray_int32_T *iwork;
  emxArray_creal_T *xwork;
  int n;
  unsigned int unnamed_idx_0;
  creal_T b_vwork;
  creal_T c_vwork;
  boolean_T p;
  int i2;
  int b_j;
  int pEnd;
  int b_p;
  int q;
  int qEnd;
  int kEnd;
  if (dim <= 1) {
    i24 = x->size[0];
  } else {
    i24 = 1;
  }

  emxInit_creal_T(&vwork, 1);
  i = vwork->size[0];
  vwork->size[0] = i24;
  emxEnsureCapacity((emxArray__common *)vwork, i, (int)sizeof(creal_T));
  x_idx_0 = x->size[0];
  i = idx->size[0];
  idx->size[0] = x_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, i, (int)sizeof(int));
  vstride = 1;
  x_idx_0 = 1;
  while (x_idx_0 <= dim - 1) {
    vstride *= x->size[0];
    x_idx_0 = 2;
  }

  j = 0;
  emxInit_int32_T(&b_idx, 1);
  emxInit_int32_T(&iwork, 1);
  emxInit_creal_T(&xwork, 1);
  while (j + 1 <= vstride) {
    for (x_idx_0 = 0; x_idx_0 + 1 <= i24; x_idx_0++) {
      vwork->data[x_idx_0] = x->data[j + x_idx_0 * vstride];
    }

    n = vwork->size[0];
    unnamed_idx_0 = (unsigned int)vwork->size[0];
    i = b_idx->size[0];
    b_idx->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)b_idx, i, (int)sizeof(int));
    x_idx_0 = (int)unnamed_idx_0;
    for (i = 0; i < x_idx_0; i++) {
      b_idx->data[i] = 0;
    }

    if (vwork->size[0] != 0) {
      i = b_idx->size[0];
      b_idx->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)b_idx, i, (int)sizeof(int));
      x_idx_0 = (int)unnamed_idx_0;
      for (i = 0; i < x_idx_0; i++) {
        b_idx->data[i] = 0;
      }

      i = iwork->size[0];
      iwork->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)iwork, i, (int)sizeof(int));
      for (x_idx_0 = 1; x_idx_0 <= n - 1; x_idx_0 += 2) {
        b_vwork = vwork->data[x_idx_0 - 1];
        c_vwork = vwork->data[x_idx_0];
        if (relop(b_vwork, c_vwork) || (rtIsNaN(vwork->data[x_idx_0 - 1].re) ||
             rtIsNaN(vwork->data[x_idx_0 - 1].im))) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          b_idx->data[x_idx_0 - 1] = x_idx_0;
          b_idx->data[x_idx_0] = x_idx_0 + 1;
        } else {
          b_idx->data[x_idx_0 - 1] = x_idx_0 + 1;
          b_idx->data[x_idx_0] = x_idx_0;
        }
      }

      if ((vwork->size[0] & 1) != 0) {
        b_idx->data[vwork->size[0] - 1] = vwork->size[0];
      }

      i = 2;
      while (i < n) {
        i2 = i << 1;
        b_j = 1;
        for (pEnd = 1 + i; pEnd < n + 1; pEnd = qEnd + i) {
          b_p = b_j - 1;
          q = pEnd;
          qEnd = b_j + i2;
          if (qEnd > n + 1) {
            qEnd = n + 1;
          }

          x_idx_0 = 0;
          kEnd = qEnd - b_j;
          while (x_idx_0 + 1 <= kEnd) {
            b_vwork = vwork->data[b_idx->data[b_p] - 1];
            c_vwork = vwork->data[b_idx->data[q - 1] - 1];
            if (relop(b_vwork, c_vwork) || (rtIsNaN(vwork->data[b_idx->data[b_p]
                  - 1].re) || rtIsNaN(vwork->data[b_idx->data[b_p] - 1].im))) {
              p = true;
            } else {
              p = false;
            }

            if (p) {
              iwork->data[x_idx_0] = b_idx->data[b_p];
              b_p++;
              if (b_p + 1 == pEnd) {
                while (q < qEnd) {
                  x_idx_0++;
                  iwork->data[x_idx_0] = b_idx->data[q - 1];
                  q++;
                }
              }
            } else {
              iwork->data[x_idx_0] = b_idx->data[q - 1];
              q++;
              if (q == qEnd) {
                while (b_p + 1 < pEnd) {
                  x_idx_0++;
                  iwork->data[x_idx_0] = b_idx->data[b_p];
                  b_p++;
                }
              }
            }

            x_idx_0++;
          }

          for (x_idx_0 = 0; x_idx_0 + 1 <= kEnd; x_idx_0++) {
            b_idx->data[(b_j + x_idx_0) - 1] = iwork->data[x_idx_0];
          }

          b_j = qEnd;
        }

        i = i2;
      }

      unnamed_idx_0 = (unsigned int)vwork->size[0];
      i = xwork->size[0];
      xwork->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)xwork, i, (int)sizeof(creal_T));
      for (x_idx_0 = 0; x_idx_0 + 1 <= n; x_idx_0++) {
        xwork->data[x_idx_0] = vwork->data[x_idx_0];
      }

      for (x_idx_0 = 0; x_idx_0 + 1 <= n; x_idx_0++) {
        vwork->data[x_idx_0] = xwork->data[b_idx->data[x_idx_0] - 1];
      }
    }

    for (x_idx_0 = 0; x_idx_0 + 1 <= i24; x_idx_0++) {
      x->data[j + x_idx_0 * vstride] = vwork->data[x_idx_0];
      idx->data[j + x_idx_0 * vstride] = b_idx->data[x_idx_0];
    }

    j++;
  }

  emxFree_creal_T(&xwork);
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&b_idx);
  emxFree_creal_T(&vwork);
}

void b_sort(emxArray_real_T *x)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  c_sort(x, dim);
}

void d_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  e_sort(x, dim, idx);
}

void f_sort(emxArray_real_T *x)
{
  int dim;
  int i23;
  emxArray_real_T *vwork;
  int j;
  int vstride;
  int k;
  emxArray_int32_T *b_vwork;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  if (dim <= 1) {
    i23 = x->size[0];
  } else {
    i23 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  j = vwork->size[0];
  vwork->size[0] = i23;
  emxEnsureCapacity((emxArray__common *)vwork, j, (int)sizeof(double));
  vstride = 1;
  k = 1;
  while (k <= dim - 1) {
    vstride *= x->size[0];
    k = 2;
  }

  j = 0;
  emxInit_int32_T(&b_vwork, 1);
  while (j + 1 <= vstride) {
    for (k = 0; k + 1 <= i23; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    sortIdx(vwork, b_vwork);
    for (k = 0; k + 1 <= i23; k++) {
      x->data[j + k * vstride] = vwork->data[k];
    }

    j++;
  }

  emxFree_int32_T(&b_vwork);
  emxFree_real_T(&vwork);
}

void g_sort(emxArray_creal_T *x, emxArray_int32_T *idx)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  h_sort(x, dim, idx);
}

void sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int dim;
  int i13;
  emxArray_real_T *vwork;
  int j;
  int vstride;
  int k;
  emxArray_int32_T *iidx;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  if (dim <= 1) {
    i13 = x->size[0];
  } else {
    i13 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  j = vwork->size[0];
  vwork->size[0] = i13;
  emxEnsureCapacity((emxArray__common *)vwork, j, (int)sizeof(double));
  vstride = x->size[0];
  j = idx->size[0];
  idx->size[0] = vstride;
  emxEnsureCapacity((emxArray__common *)idx, j, (int)sizeof(int));
  vstride = 1;
  k = 1;
  while (k <= dim - 1) {
    vstride *= x->size[0];
    k = 2;
  }

  j = 0;
  emxInit_int32_T(&iidx, 1);
  while (j + 1 <= vstride) {
    for (k = 0; k + 1 <= i13; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    sortIdx(vwork, iidx);
    for (k = 0; k + 1 <= i13; k++) {
      x->data[j + k * vstride] = vwork->data[k];
      idx->data[j + k * vstride] = iidx->data[k];
    }

    j++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

/* End of code generation (sort1.cpp) */
