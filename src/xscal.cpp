/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xscal.cpp
 *
 * Code generation for function 'xscal'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xscal.h"

/* Function Definitions */
void xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i13;
  int k;
  i13 = (ix0 + n) - 1;
  for (k = ix0; k <= i13; k++) {
    x->data[k - 1] *= a;
  }
}

/* End of code generation (xscal.cpp) */
