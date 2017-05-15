/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * inv.cpp
 *
 * Code generation for function 'inv'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "inv.h"
#include "MAVEfast_emxutil.h"
#include "xtrsm.h"
#include "colon.h"
#include "xzgetrf.h"

/* Function Definitions */
void inv(const emxArray_real_T *x, emxArray_real_T *y)
{
  int n;
  int pipk;
  int c;
  emxArray_real_T *b_x;
  emxArray_int32_T *p;
  emxArray_int32_T *ipiv;
  int k;
  int i;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    pipk = y->size[0] * y->size[1];
    y->size[0] = x->size[0];
    y->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)y, pipk, (int)sizeof(double));
    c = x->size[0] * x->size[1];
    for (pipk = 0; pipk < c; pipk++) {
      y->data[pipk] = x->data[pipk];
    }
  } else {
    n = x->size[0];
    pipk = y->size[0] * y->size[1];
    y->size[0] = x->size[0];
    y->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)y, pipk, (int)sizeof(double));
    c = x->size[0] * x->size[1];
    for (pipk = 0; pipk < c; pipk++) {
      y->data[pipk] = 0.0;
    }

    emxInit_real_T(&b_x, 2);
    pipk = b_x->size[0] * b_x->size[1];
    b_x->size[0] = x->size[0];
    b_x->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)b_x, pipk, (int)sizeof(double));
    c = x->size[0] * x->size[1];
    for (pipk = 0; pipk < c; pipk++) {
      b_x->data[pipk] = x->data[pipk];
    }

    emxInit_int32_T(&p, 2);
    emxInit_int32_T(&ipiv, 2);
    xzgetrf(x->size[0], x->size[0], b_x, x->size[0], ipiv, &pipk);
    eml_signed_integer_colon(x->size[0], p);
    for (k = 0; k < ipiv->size[1]; k++) {
      if (ipiv->data[k] > 1 + k) {
        pipk = p->data[ipiv->data[k] - 1];
        p->data[ipiv->data[k] - 1] = p->data[k];
        p->data[k] = pipk;
      }
    }

    emxFree_int32_T(&ipiv);
    for (k = 0; k + 1 <= n; k++) {
      c = p->data[k] - 1;
      y->data[k + y->size[0] * (p->data[k] - 1)] = 1.0;
      for (pipk = k; pipk + 1 <= n; pipk++) {
        if (y->data[pipk + y->size[0] * c] != 0.0) {
          for (i = pipk + 1; i + 1 <= n; i++) {
            y->data[i + y->size[0] * c] -= y->data[pipk + y->size[0] * c] *
              b_x->data[i + b_x->size[0] * pipk];
          }
        }
      }
    }

    emxFree_int32_T(&p);
    xtrsm(x->size[0], x->size[0], b_x, x->size[0], y, x->size[0]);
    emxFree_real_T(&b_x);
  }
}

/* End of code generation (inv.cpp) */
