/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzggbal.cpp
 *
 * Code generation for function 'xzggbal'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "xzggbal.h"
#include "CVfast_emxutil.h"

/* Function Definitions */
void xzggbal(emxArray_creal_T *A, int *ilo, int *ihi, emxArray_int32_T *rscale)
{
  int jj;
  int loop_ub;
  emxArray_creal_T *b_A;
  int exitg2;
  int i;
  int j;
  boolean_T found;
  int ii;
  boolean_T exitg5;
  int nzcount;
  boolean_T exitg6;
  int exitg1;
  boolean_T c_A;
  boolean_T guard2 = false;
  boolean_T exitg3;
  double atmp_re;
  double atmp_im;
  boolean_T exitg4;
  boolean_T guard1 = false;
  jj = rscale->size[0];
  rscale->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)rscale, jj, (int)sizeof(int));
  loop_ub = A->size[0];
  for (jj = 0; jj < loop_ub; jj++) {
    rscale->data[jj] = 1;
  }

  *ilo = 1;
  *ihi = A->size[0];
  if (A->size[0] <= 1) {
    *ihi = 1;
  } else {
    emxInit_creal_T1(&b_A, 2);
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      found = false;
      ii = *ihi;
      exitg5 = false;
      while ((!exitg5) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = *ihi;
        jj = 1;
        exitg6 = false;
        while ((!exitg6) && (jj <= *ihi)) {
          c_A = ((A->data[(ii + A->size[0] * (jj - 1)) - 1].re != 0.0) ||
                 (A->data[(ii + A->size[0] * (jj - 1)) - 1].im != 0.0));
          guard2 = false;
          if (c_A || (ii == jj)) {
            if (nzcount == 0) {
              j = jj;
              nzcount = 1;
              guard2 = true;
            } else {
              nzcount = 2;
              exitg6 = true;
            }
          } else {
            guard2 = true;
          }

          if (guard2) {
            jj++;
          }
        }

        if (nzcount < 2) {
          found = true;
          exitg5 = true;
        } else {
          ii--;
        }
      }

      if (!found) {
        exitg2 = 2;
      } else {
        jj = b_A->size[0] * b_A->size[1];
        b_A->size[0] = A->size[0];
        b_A->size[1] = A->size[1];
        emxEnsureCapacity((emxArray__common *)b_A, jj, (int)sizeof(creal_T));
        loop_ub = A->size[0] * A->size[1];
        for (jj = 0; jj < loop_ub; jj++) {
          b_A->data[jj] = A->data[jj];
        }

        if (i != *ihi) {
          for (ii = 0; ii + 1 <= A->size[0]; ii++) {
            atmp_re = b_A->data[(i + b_A->size[0] * ii) - 1].re;
            atmp_im = b_A->data[(i + b_A->size[0] * ii) - 1].im;
            b_A->data[(i + b_A->size[0] * ii) - 1] = b_A->data[(*ihi + b_A->
              size[0] * ii) - 1];
            b_A->data[(*ihi + b_A->size[0] * ii) - 1].re = atmp_re;
            b_A->data[(*ihi + b_A->size[0] * ii) - 1].im = atmp_im;
          }
        }

        if (j != *ihi) {
          for (ii = 0; ii + 1 <= *ihi; ii++) {
            atmp_re = b_A->data[ii + b_A->size[0] * (j - 1)].re;
            atmp_im = b_A->data[ii + b_A->size[0] * (j - 1)].im;
            b_A->data[ii + b_A->size[0] * (j - 1)] = b_A->data[ii + b_A->size[0]
              * (*ihi - 1)];
            b_A->data[ii + b_A->size[0] * (*ihi - 1)].re = atmp_re;
            b_A->data[ii + b_A->size[0] * (*ihi - 1)].im = atmp_im;
          }
        }

        jj = A->size[0] * A->size[1];
        A->size[0] = b_A->size[0];
        A->size[1] = b_A->size[1];
        emxEnsureCapacity((emxArray__common *)A, jj, (int)sizeof(creal_T));
        loop_ub = b_A->size[1];
        for (jj = 0; jj < loop_ub; jj++) {
          ii = b_A->size[0];
          for (nzcount = 0; nzcount < ii; nzcount++) {
            A->data[nzcount + A->size[0] * jj] = b_A->data[nzcount + b_A->size[0]
              * jj];
          }
        }

        rscale->data[*ihi - 1] = j;
        (*ihi)--;
        if (*ihi == 1) {
          rscale->data[0] = 1;
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0);

    if (exitg2 == 1) {
    } else {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        found = false;
        jj = *ilo;
        exitg3 = false;
        while ((!exitg3) && (jj <= *ihi)) {
          nzcount = 0;
          i = *ihi;
          j = jj;
          ii = *ilo;
          exitg4 = false;
          while ((!exitg4) && (ii <= *ihi)) {
            c_A = ((A->data[(ii + A->size[0] * (jj - 1)) - 1].re != 0.0) ||
                   (A->data[(ii + A->size[0] * (jj - 1)) - 1].im != 0.0));
            guard1 = false;
            if (c_A || (ii == jj)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                guard1 = true;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              guard1 = true;
            }

            if (guard1) {
              ii++;
            }
          }

          if (nzcount < 2) {
            found = true;
            exitg3 = true;
          } else {
            jj++;
          }
        }

        if (!found) {
          exitg1 = 1;
        } else {
          jj = b_A->size[0] * b_A->size[1];
          b_A->size[0] = A->size[0];
          b_A->size[1] = A->size[1];
          emxEnsureCapacity((emxArray__common *)b_A, jj, (int)sizeof(creal_T));
          loop_ub = A->size[0] * A->size[1];
          for (jj = 0; jj < loop_ub; jj++) {
            b_A->data[jj] = A->data[jj];
          }

          if (i != *ilo) {
            for (ii = *ilo - 1; ii + 1 <= A->size[0]; ii++) {
              atmp_re = b_A->data[(i + b_A->size[0] * ii) - 1].re;
              atmp_im = b_A->data[(i + b_A->size[0] * ii) - 1].im;
              b_A->data[(i + b_A->size[0] * ii) - 1] = b_A->data[(*ilo +
                b_A->size[0] * ii) - 1];
              b_A->data[(*ilo + b_A->size[0] * ii) - 1].re = atmp_re;
              b_A->data[(*ilo + b_A->size[0] * ii) - 1].im = atmp_im;
            }
          }

          if (j != *ilo) {
            for (ii = 0; ii + 1 <= *ihi; ii++) {
              atmp_re = b_A->data[ii + b_A->size[0] * (j - 1)].re;
              atmp_im = b_A->data[ii + b_A->size[0] * (j - 1)].im;
              b_A->data[ii + b_A->size[0] * (j - 1)] = b_A->data[ii + b_A->size
                [0] * (*ilo - 1)];
              b_A->data[ii + b_A->size[0] * (*ilo - 1)].re = atmp_re;
              b_A->data[ii + b_A->size[0] * (*ilo - 1)].im = atmp_im;
            }
          }

          jj = A->size[0] * A->size[1];
          A->size[0] = b_A->size[0];
          A->size[1] = b_A->size[1];
          emxEnsureCapacity((emxArray__common *)A, jj, (int)sizeof(creal_T));
          loop_ub = b_A->size[1];
          for (jj = 0; jj < loop_ub; jj++) {
            ii = b_A->size[0];
            for (nzcount = 0; nzcount < ii; nzcount++) {
              A->data[nzcount + A->size[0] * jj] = b_A->data[nzcount + b_A->
                size[0] * jj];
            }
          }

          rscale->data[*ilo - 1] = j;
          (*ilo)++;
          if (*ilo == *ihi) {
            rscale->data[*ilo - 1] = *ilo;
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }

    emxFree_creal_T(&b_A);
  }
}

/* End of code generation (xzggbal.cpp) */
