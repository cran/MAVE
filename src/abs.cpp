/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * abs.cpp
 *
 * Code generation for function 'abs'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "abs.h"
#include "CVfast_emxutil.h"

/* Function Definitions */
void b_abs(const emxArray_boolean_T *x, emxArray_real_T *y)
{
  unsigned int uv2[2];
  int n;
  int k;
  for (n = 0; n < 2; n++) {
    uv2[n] = (unsigned int)x->size[n];
  }

  n = y->size[0] * y->size[1];
  y->size[0] = (int)uv2[0];
  y->size[1] = (int)uv2[1];
  emxEnsureCapacity((emxArray__common *)y, n, (int)sizeof(double));
  n = x->size[0] * x->size[1];
  for (k = 0; k + 1 <= n; k++) {
    y->data[k] = x->data[k];
  }
}

void c_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int x_idx_0;
  int k;
  x_idx_0 = (unsigned int)x->size[0];
  k = y->size[0];
  y->size[0] = (int)x_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k + 1 <= x->size[0]; k++) {
    y->data[k] = std::abs(x->data[k]);
  }
}

/* End of code generation (abs.cpp) */
