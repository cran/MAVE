/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * combine_vector_elements.cpp
 *
 * Code generation for function 'combine_vector_elements'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "combine_vector_elements.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void combine_vector_elements(const emxArray_real_T *x, emxArray_real_T *y)
{
  int vstride;
  int j;
  double s;
  int k;
  vstride = y->size[0];
  y->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)y, vstride, (int)sizeof(double));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    j = y->size[0];
    vstride = y->size[0];
    y->size[0] = j;
    emxEnsureCapacity((emxArray__common *)y, vstride, (int)sizeof(double));
    for (vstride = 0; vstride < j; vstride++) {
      y->data[vstride] = 0.0;
    }
  } else {
    vstride = x->size[0];
    for (j = 0; j + 1 <= vstride; j++) {
      s = x->data[j];
      for (k = 2; k <= x->size[1]; k++) {
        s += x->data[j + (k - 1) * vstride];
      }

      y->data[j] = s;
    }
  }
}

/* End of code generation (combine_vector_elements.cpp) */
