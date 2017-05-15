/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * diag.cpp
 *
 * Code generation for function 'diag'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "diag.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int unnamed_idx_0;
  int unnamed_idx_1;
  int i6;
  unnamed_idx_0 = v->size[0];
  unnamed_idx_1 = v->size[0];
  i6 = d->size[0] * d->size[1];
  d->size[0] = unnamed_idx_0;
  d->size[1] = unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)d, i6, (int)sizeof(double));
  unnamed_idx_0 *= unnamed_idx_1;
  for (i6 = 0; i6 < unnamed_idx_0; i6++) {
    d->data[i6] = 0.0;
  }

  for (unnamed_idx_0 = 0; unnamed_idx_0 + 1 <= v->size[0]; unnamed_idx_0++) {
    d->data[unnamed_idx_0 + d->size[0] * unnamed_idx_0] = v->data[unnamed_idx_0];
  }
}

void c_diag(const emxArray_creal_T *v, emxArray_creal_T *d)
{
  int j;
  int dlen;
  int stride;
  if ((v->size[0] == 1) && (v->size[1] == 1)) {
    j = d->size[0];
    d->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)d, j, (int)sizeof(creal_T));
    d->data[0] = v->data[0];
  } else {
    if (0 < v->size[1]) {
      if (v->size[0] <= v->size[1]) {
        dlen = v->size[0];
      } else {
        dlen = v->size[1];
      }

      stride = v->size[0] + 1;
    } else {
      dlen = 0;
      stride = 0;
    }

    j = d->size[0];
    d->size[0] = dlen;
    emxEnsureCapacity((emxArray__common *)d, j, (int)sizeof(creal_T));
    for (j = 0; j + 1 <= dlen; j++) {
      d->data[j] = v->data[j * stride];
    }
  }
}

void diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  int j;
  int dlen;
  int stride;
  if ((v->size[0] == 1) && (v->size[1] == 1)) {
    j = d->size[0];
    d->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)d, j, (int)sizeof(double));
    d->data[0] = v->data[0];
  } else {
    if (0 < v->size[1]) {
      if (v->size[0] <= v->size[1]) {
        dlen = v->size[0];
      } else {
        dlen = v->size[1];
      }

      stride = v->size[0] + 1;
    } else {
      dlen = 0;
      stride = 0;
    }

    j = d->size[0];
    d->size[0] = dlen;
    emxEnsureCapacity((emxArray__common *)d, j, (int)sizeof(double));
    for (j = 0; j + 1 <= dlen; j++) {
      d->data[j] = v->data[j * stride];
    }
  }
}

/* End of code generation (diag.cpp) */
