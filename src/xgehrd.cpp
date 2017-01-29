/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xgehrd.cpp
 *
 * Code generation for function 'xgehrd'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "xgehrd.h"
#include "xzlarf.h"
#include "xscal.h"
#include "relop.h"
#include "xnrm2.h"
#include "CVfast_emxutil.h"
#include "CVfast_rtwutil.h"

/* Function Definitions */
void xgehrd(emxArray_real_T *a, emxArray_real_T *tau)
{
  int n;
  int ntau;
  emxArray_real_T *work;
  int i18;
  int i;
  int im1n;
  int in;
  int b_i;
  int c;
  double alpha1;
  double d2;
  double xnorm;
  int jy;
  int lastv;
  int knt;
  int lastc;
  boolean_T exitg2;
  int ix;
  int ia;
  int exitg1;
  n = a->size[0];
  if (a->size[0] < 1) {
    ntau = 0;
  } else {
    ntau = a->size[0] - 1;
  }

  emxInit_real_T1(&work, 1);
  i18 = tau->size[0];
  tau->size[0] = ntau;
  emxEnsureCapacity((emxArray__common *)tau, i18, (int)sizeof(double));
  ntau = a->size[0];
  i18 = work->size[0];
  work->size[0] = ntau;
  emxEnsureCapacity((emxArray__common *)work, i18, (int)sizeof(double));
  for (i18 = 0; i18 < ntau; i18++) {
    work->data[i18] = 0.0;
  }

  for (i = 0; i + 1 < n; i++) {
    im1n = i * n + 2;
    in = (i + 1) * n;
    if (i + 3 <= n) {
      b_i = i + 3;
    } else {
      b_i = n;
    }

    ntau = b_i + i * n;
    c = (n - i) - 2;
    alpha1 = a->data[(i + a->size[0] * i) + 1];
    d2 = 0.0;
    if (!(c + 1 <= 0)) {
      xnorm = xnrm2(c, a, ntau);
      if (xnorm != 0.0) {
        xnorm = rt_hypotd_snf(a->data[(i + a->size[0] * i) + 1], xnorm);
        if (a->data[(i + a->size[0] * i) + 1] >= 0.0) {
          xnorm = -xnorm;
        }

        if (std::abs(xnorm) < 1.0020841800044864E-292) {
          knt = 0;
          do {
            knt++;
            xscal(c, 9.9792015476736E+291, a, ntau);
            xnorm *= 9.9792015476736E+291;
            alpha1 *= 9.9792015476736E+291;
          } while (!(std::abs(xnorm) >= 1.0020841800044864E-292));

          xnorm = xnrm2(c, a, ntau);
          xnorm = rt_hypotd_snf(alpha1, xnorm);
          if (alpha1 >= 0.0) {
            xnorm = -xnorm;
          }

          d2 = (xnorm - alpha1) / xnorm;
          xscal(c, 1.0 / (alpha1 - xnorm), a, ntau);
          for (ntau = 1; ntau <= knt; ntau++) {
            xnorm *= 1.0020841800044864E-292;
          }

          alpha1 = xnorm;
        } else {
          d2 = (xnorm - a->data[(i + a->size[0] * i) + 1]) / xnorm;
          xscal(c, 1.0 / (a->data[(i + a->size[0] * i) + 1] - xnorm), a, ntau);
          alpha1 = xnorm;
        }
      }
    }

    tau->data[i] = d2;
    a->data[(i + a->size[0] * i) + 1] = 1.0;
    c = (n - i) - 3;
    jy = (i + im1n) - 1;
    if (tau->data[i] != 0.0) {
      lastv = c + 2;
      ntau = jy + c;
      while ((lastv > 0) && (a->data[ntau + 1] == 0.0)) {
        lastv--;
        ntau--;
      }

      lastc = n;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        ntau = in + lastc;
        ia = ntau;
        do {
          exitg1 = 0;
          if ((n > 0) && (ia <= ntau + (lastv - 1) * n)) {
            if (a->data[ia - 1] != 0.0) {
              exitg1 = 1;
            } else {
              ia += n;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      if (lastc != 0) {
        for (ntau = 1; ntau <= lastc; ntau++) {
          work->data[ntau - 1] = 0.0;
        }

        ix = jy;
        i18 = (in + n * (lastv - 1)) + 1;
        knt = in + 1;
        while ((n > 0) && (knt <= i18)) {
          ntau = 0;
          c = (knt + lastc) - 1;
          for (ia = knt; ia <= c; ia++) {
            work->data[ntau] += a->data[ia - 1] * a->data[ix];
            ntau++;
          }

          ix++;
          knt += n;
        }
      }

      if (!(-tau->data[i] == 0.0)) {
        ntau = in;
        for (knt = 1; knt <= lastv; knt++) {
          if (a->data[jy] != 0.0) {
            xnorm = a->data[jy] * -tau->data[i];
            ix = 0;
            i18 = lastc + ntau;
            for (c = ntau; c + 1 <= i18; c++) {
              a->data[c] += work->data[ix] * xnorm;
              ix++;
            }
          }

          jy++;
          ntau += n;
        }
      }
    }

    xzlarf((n - i) - 1, (n - i) - 1, i + im1n, tau->data[i], a, (i + in) + 2, n,
           work);
    a->data[(i + a->size[0] * i) + 1] = alpha1;
  }

  emxFree_real_T(&work);
}

/* End of code generation (xgehrd.cpp) */
