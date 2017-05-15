/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * mean.cpp
 *
 * Code generation for function 'mean'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "mean.h"
#include "MAVEfast_emxutil.h"
#include "combine_vector_elements.h"

/* Function Definitions */
void b_mean(const emxArray_real_T *x, emxArray_real_T *y)
{
  int b_x;
  int i4;
  int loop_ub;
  combine_vector_elements(x, y);
  b_x = x->size[1];
  i4 = y->size[0];
  emxEnsureCapacity((emxArray__common *)y, i4, (int)sizeof(double));
  loop_ub = y->size[0];
  for (i4 = 0; i4 < loop_ub; i4++) {
    y->data[i4] /= (double)b_x;
  }
}

double c_mean(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[1]; k++) {
      y += x->data[k - 1];
    }
  }

  y /= (double)x->size[1];
  return y;
}

void mean(const emxArray_real_T *x, emxArray_real_T *y)
{
  int i;
  int vlen;
  int xoffset;
  double s;
  int k;
  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
    vlen = y->size[1];
    for (i = 0; i < vlen; i++) {
      y->data[y->size[0] * i] = 0.0;
    }
  } else {
    vlen = x->size[0];
    for (i = 0; i + 1 <= x->size[1]; i++) {
      xoffset = i * vlen;
      s = x->data[xoffset];
      for (k = 2; k <= vlen; k++) {
        s += x->data[(xoffset + k) - 1];
      }

      y->data[i] = s;
    }
  }

  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
  i = y->size[0];
  vlen = y->size[1];
  xoffset = x->size[0];
  vlen *= i;
  for (i = 0; i < vlen; i++) {
    y->data[i] /= (double)xoffset;
  }
}

/* End of code generation (mean.cpp) */
