/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * power.cpp
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "power.h"
#include "CVfast_emxutil.h"

/* Function Definitions */
void b_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  emxArray_real_T *x;
  int ub_loop;
  int loop_ub;
  unsigned int a_idx_0;
  int k;
  emxInit_real_T1(&x, 1);
  ub_loop = x->size[0];
  x->size[0] = a->size[0];
  emxEnsureCapacity((emxArray__common *)x, ub_loop, (int)sizeof(double));
  loop_ub = a->size[0];
  for (ub_loop = 0; ub_loop < loop_ub; ub_loop++) {
    x->data[ub_loop] = a->data[ub_loop];
  }

  a_idx_0 = (unsigned int)a->size[0];
  ub_loop = y->size[0];
  y->size[0] = (int)a_idx_0;
  emxEnsureCapacity((emxArray__common *)y, ub_loop, (int)sizeof(double));
  ub_loop = a->size[0];

#pragma omp parallel for \
 num_threads(omp_get_max_threads())

  for (k = 1; k <= ub_loop; k++) {
    y->data[k - 1] = x->data[k - 1] * x->data[k - 1];
  }

  emxFree_real_T(&x);
}

void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  emxArray_real_T *x;
  int ub_loop;
  int loop_ub;
  unsigned int uv1[2];
  int k;
  emxInit_real_T(&x, 2);
  ub_loop = x->size[0] * x->size[1];
  x->size[0] = a->size[0];
  x->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)x, ub_loop, (int)sizeof(double));
  loop_ub = a->size[0] * a->size[1];
  for (ub_loop = 0; ub_loop < loop_ub; ub_loop++) {
    x->data[ub_loop] = a->data[ub_loop];
  }

  for (ub_loop = 0; ub_loop < 2; ub_loop++) {
    uv1[ub_loop] = (unsigned int)a->size[ub_loop];
  }

  ub_loop = y->size[0] * y->size[1];
  y->size[0] = (int)uv1[0];
  y->size[1] = (int)uv1[1];
  emxEnsureCapacity((emxArray__common *)y, ub_loop, (int)sizeof(double));
  ub_loop = a->size[0] * a->size[1];

#pragma omp parallel for \
 num_threads(omp_get_max_threads())

  for (k = 1; k <= ub_loop; k++) {
    y->data[k - 1] = x->data[k - 1] * x->data[k - 1];
  }

  emxFree_real_T(&x);
}

/* End of code generation (power.cpp) */
