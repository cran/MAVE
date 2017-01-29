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
#include "CVfast.h"
#include "MAVEfast.h"
#include "mean.h"
#include "CVfast_emxutil.h"
#include "combine_vector_elements.h"

/* Function Definitions */
void b_mean(const emxArray_real_T *x, emxArray_real_T *y)
{
  int b_y;
  int c_y;
  int b_x;
  b_combine_vector_elements(x, y);
  b_y = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)y, b_y, (int)sizeof(double));
  b_y = y->size[0];
  c_y = y->size[1];
  b_x = x->size[0];
  c_y *= b_y;
  for (b_y = 0; b_y < c_y; b_y++) {
    y->data[b_y] /= (double)b_x;
  }
}

void c_mean(const emxArray_real_T *x, emxArray_real_T *y)
{
  int b_x;
  int i11;
  int loop_ub;
  combine_vector_elements(x, y);
  b_x = x->size[1];
  i11 = y->size[0];
  emxEnsureCapacity((emxArray__common *)y, i11, (int)sizeof(double));
  loop_ub = y->size[0];
  for (i11 = 0; i11 < loop_ub; i11++) {
    y->data[i11] /= (double)b_x;
  }
}

double mean(const emxArray_real_T *x)
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

/* End of code generation (mean.cpp) */
