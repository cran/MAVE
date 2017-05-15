/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * std.cpp
 *
 * Code generation for function 'std'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "std.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_std(const emxArray_real_T *varargin_1, emxArray_real_T *y)
{
  emxArray_real_T *x;
  int i7;
  int loop_ub;
  int n;
  double d;
  emxArray_real_T *b_y;
  int i;
  int ix;
  double c_y;
  double xbar;
  int k;
  double r;
  emxInit_real_T(&x, 2);
  i7 = x->size[0] * x->size[1];
  x->size[0] = varargin_1->size[0];
  x->size[1] = varargin_1->size[1];
  emxEnsureCapacity((emxArray__common *)x, i7, (int)sizeof(double));
  loop_ub = varargin_1->size[0] * varargin_1->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    x->data[i7] = varargin_1->data[i7];
  }

  n = varargin_1->size[0];
  if (varargin_1->size[0] > 1) {
    d = (double)varargin_1->size[0] - 1.0;
  } else {
    d = varargin_1->size[0];
  }

  emxInit_real_T(&b_y, 2);
  i7 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = varargin_1->size[1];
  emxEnsureCapacity((emxArray__common *)b_y, i7, (int)sizeof(double));
  if ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0)) {
    i7 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)b_y, i7, (int)sizeof(double));
    loop_ub = b_y->size[1];
    for (i7 = 0; i7 < loop_ub; i7++) {
      b_y->data[b_y->size[0] * i7] = 0.0;
    }
  } else {
    loop_ub = varargin_1->size[1];

//#pragma omp parallel for \
 num_threads(omp_get_max_threads()) \
 private(ix,c_y,xbar,k,r)

    for (i = 1; i <= loop_ub; i++) {
      ix = x->size[0];
      if (ix == 0) {
        c_y = rtNaN;
      } else {
        ix = 0;
        xbar = x->data[x->size[0] * (i - 1)];
        for (k = 2; k <= n; k++) {
          ix++;
          xbar += x->data[ix + x->size[0] * (i - 1)];
        }

        xbar /= (double)n;
        ix = 0;
        r = x->data[x->size[0] * (i - 1)] - xbar;
        c_y = r * r;
        for (k = 2; k <= n; k++) {
          ix++;
          r = x->data[ix + x->size[0] * (i - 1)] - xbar;
          c_y += r * r;
        }

        c_y /= d;
      }

      b_y->data[i - 1] = c_y;
    }
  }

  emxFree_real_T(&x);
  i7 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
  loop_ub = b_y->size[0] * b_y->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    y->data[i7] = b_y->data[i7];
  }

  for (loop_ub = 0; loop_ub + 1 <= b_y->size[1]; loop_ub++) {
    y->data[loop_ub] = std::sqrt(y->data[loop_ub]);
  }

  emxFree_real_T(&b_y);
}

/* End of code generation (std.cpp) */
