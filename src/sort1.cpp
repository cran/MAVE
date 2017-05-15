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
#include "MAVEfast.h"
#include "sort1.h"
#include "sortIdx.h"
#include "MAVEfast_emxutil.h"

/* Function Declarations */
static void b_sort(emxArray_real_T *x);
static void d_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx);

/* Function Definitions */
static void b_sort(emxArray_real_T *x)
{
  emxArray_real_T *vwork;
  int i16;
  int x_idx_0;
  int i17;
  emxArray_int32_T *b_vwork;
  emxInit_real_T1(&vwork, 1);
  i16 = x->size[1];
  x_idx_0 = x->size[1];
  i17 = vwork->size[0];
  vwork->size[0] = x_idx_0;
  emxEnsureCapacity((emxArray__common *)vwork, i17, (int)sizeof(double));
  for (x_idx_0 = 0; x_idx_0 + 1 <= i16; x_idx_0++) {
    vwork->data[x_idx_0] = x->data[x_idx_0];
  }

  emxInit_int32_T1(&b_vwork, 1);
  sortIdx(vwork, b_vwork);
  x_idx_0 = 0;
  emxFree_int32_T(&b_vwork);
  while (x_idx_0 + 1 <= i16) {
    x->data[x_idx_0] = vwork->data[x_idx_0];
    x_idx_0++;
  }

  emxFree_real_T(&vwork);
}

static void d_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx)
{
  int i18;
  emxArray_real_T *vwork;
  int vstride;
  int x_idx_0;
  int j;
  emxArray_int32_T *iidx;
  if (dim <= 1) {
    i18 = x->size[0];
  } else {
    i18 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  vstride = vwork->size[0];
  vwork->size[0] = i18;
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
  emxInit_int32_T1(&iidx, 1);
  while (j + 1 <= vstride) {
    for (x_idx_0 = 0; x_idx_0 + 1 <= i18; x_idx_0++) {
      vwork->data[x_idx_0] = x->data[j + x_idx_0 * vstride];
    }

    sortIdx(vwork, iidx);
    for (x_idx_0 = 0; x_idx_0 + 1 <= i18; x_idx_0++) {
      x->data[j + x_idx_0 * vstride] = vwork->data[x_idx_0];
      idx->data[j + x_idx_0 * vstride] = iidx->data[x_idx_0];
    }

    j++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

void c_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  d_sort(x, dim, idx);
}

void e_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int dim;
  int i19;
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
    i19 = x->size[0];
  } else {
    i19 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  j = vwork->size[0];
  vwork->size[0] = i19;
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
  emxInit_int32_T1(&iidx, 1);
  while (j + 1 <= vstride) {
    for (k = 0; k + 1 <= i19; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    b_sortIdx(vwork, iidx);
    for (k = 0; k + 1 <= i19; k++) {
      x->data[j + k * vstride] = vwork->data[k];
      idx->data[j + k * vstride] = iidx->data[k];
    }

    j++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

void f_sort(emxArray_real_T *x)
{
  int dim;
  int i20;
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
    i20 = x->size[0];
  } else {
    i20 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  j = vwork->size[0];
  vwork->size[0] = i20;
  emxEnsureCapacity((emxArray__common *)vwork, j, (int)sizeof(double));
  vstride = 1;
  k = 1;
  while (k <= dim - 1) {
    vstride *= x->size[0];
    k = 2;
  }

  j = 0;
  emxInit_int32_T1(&b_vwork, 1);
  while (j + 1 <= vstride) {
    for (k = 0; k + 1 <= i20; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    b_sortIdx(vwork, b_vwork);
    for (k = 0; k + 1 <= i20; k++) {
      x->data[j + k * vstride] = vwork->data[k];
    }

    j++;
  }

  emxFree_int32_T(&b_vwork);
  emxFree_real_T(&vwork);
}

void sort(emxArray_real_T *x)
{
  b_sort(x);
}

/* End of code generation (sort1.cpp) */
