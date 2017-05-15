/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * rdivide.cpp
 *
 * Code generation for function 'rdivide'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "rdivide.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void rdivide(const emxArray_real_T *x, const emxArray_real_T *y, emxArray_real_T
             *z)
{
  int i5;
  int loop_ub;
  i5 = z->size[0] * z->size[1];
  z->size[0] = x->size[0];
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i5, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i5 = 0; i5 < loop_ub; i5++) {
    z->data[i5] = x->data[i5] / y->data[i5];
  }
}

/* End of code generation (rdivide.cpp) */
