/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzggev.cpp
 *
 * Code generation for function 'xzggev'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "xzggev.h"
#include "CVfast_emxutil.h"
#include "xzlartg.h"
#include "xztgevc.h"
#include "xzhgeqz.h"
#include "xzggbal.h"
#include "relop.h"
#include "CVfast_rtwutil.h"

/* Function Definitions */
void xzggev(emxArray_creal_T *A, int *info, emxArray_creal_T *alpha1,
            emxArray_creal_T *beta1, emxArray_creal_T *V)
{
  int n;
  int jrow;
  int m;
  double anrm;
  boolean_T ilascl;
  boolean_T notdone;
  int k;
  boolean_T exitg1;
  double anrmto;
  double absxk;
  emxArray_int32_T *rscale;
  double ctoc;
  emxArray_int8_T *I;
  int ilo;
  int ihi;
  int b_n;
  double cfrom1;
  double cto1;
  double stemp_im;
  int jcol;
  creal_T b_A;
  creal_T c_A;
  int i;
  double c;
  creal_T tmp;
  int j;
  double stemp_re;
  *info = 0;
  n = A->size[0];
  jrow = alpha1->size[0];
  alpha1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)alpha1, jrow, (int)sizeof(creal_T));
  m = A->size[0];
  for (jrow = 0; jrow < m; jrow++) {
    alpha1->data[jrow].re = 0.0;
    alpha1->data[jrow].im = 0.0;
  }

  jrow = beta1->size[0];
  beta1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)beta1, jrow, (int)sizeof(creal_T));
  m = A->size[0];
  for (jrow = 0; jrow < m; jrow++) {
    beta1->data[jrow].re = 0.0;
    beta1->data[jrow].im = 0.0;
  }

  jrow = V->size[0] * V->size[1];
  V->size[0] = A->size[0];
  V->size[1] = A->size[0];
  emxEnsureCapacity((emxArray__common *)V, jrow, (int)sizeof(creal_T));
  m = A->size[0] * A->size[0];
  for (jrow = 0; jrow < m; jrow++) {
    V->data[jrow].re = 0.0;
    V->data[jrow].im = 0.0;
  }

  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    anrm = 0.0;
    ilascl = (A->size[0] == 0);
    notdone = (A->size[1] == 0);
    if (!(ilascl || notdone)) {
      k = 0;
      exitg1 = false;
      while ((!exitg1) && (k <= A->size[0] * A->size[1] - 1)) {
        absxk = rt_hypotd_snf(A->data[k].re, A->data[k].im);
        if (rtIsNaN(absxk)) {
          anrm = rtNaN;
          exitg1 = true;
        } else {
          if (absxk > anrm) {
            anrm = absxk;
          }

          k++;
        }
      }
    }

    if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
      jrow = alpha1->size[0];
      alpha1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)alpha1, jrow, (int)sizeof(creal_T));
      m = A->size[0];
      for (jrow = 0; jrow < m; jrow++) {
        alpha1->data[jrow].re = rtNaN;
        alpha1->data[jrow].im = 0.0;
      }

      jrow = beta1->size[0];
      beta1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)beta1, jrow, (int)sizeof(creal_T));
      m = A->size[0];
      for (jrow = 0; jrow < m; jrow++) {
        beta1->data[jrow].re = rtNaN;
        beta1->data[jrow].im = 0.0;
      }

      jrow = V->size[0] * V->size[1];
      V->size[0] = A->size[0];
      V->size[1] = A->size[0];
      emxEnsureCapacity((emxArray__common *)V, jrow, (int)sizeof(creal_T));
      m = A->size[0] * A->size[0];
      for (jrow = 0; jrow < m; jrow++) {
        V->data[jrow].re = rtNaN;
        V->data[jrow].im = 0.0;
      }
    } else {
      ilascl = false;
      anrmto = anrm;
      if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
        anrmto = 6.7178761075670888E-139;
        ilascl = true;
      } else {
        if (anrm > 1.4885657073574029E+138) {
          anrmto = 1.4885657073574029E+138;
          ilascl = true;
        }
      }

      if (ilascl) {
        absxk = anrm;
        ctoc = anrmto;
        notdone = true;
        while (notdone) {
          cfrom1 = absxk * 2.0041683600089728E-292;
          cto1 = ctoc / 4.9896007738368E+291;
          if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
            stemp_im = 2.0041683600089728E-292;
            absxk = cfrom1;
          } else if (cto1 > absxk) {
            stemp_im = 4.9896007738368E+291;
            ctoc = cto1;
          } else {
            stemp_im = ctoc / absxk;
            notdone = false;
          }

          jrow = A->size[0] * A->size[1];
          emxEnsureCapacity((emxArray__common *)A, jrow, (int)sizeof(creal_T));
          m = A->size[1];
          for (jrow = 0; jrow < m; jrow++) {
            k = A->size[0];
            for (jcol = 0; jcol < k; jcol++) {
              A->data[jcol + A->size[0] * jrow].re *= stemp_im;
              A->data[jcol + A->size[0] * jrow].im *= stemp_im;
            }
          }
        }
      }

      emxInit_int32_T(&rscale, 1);
      emxInit_int8_T(&I, 2);
      xzggbal(A, &ilo, &ihi, rscale);
      b_n = A->size[0];
      jrow = I->size[0] * I->size[1];
      I->size[0] = A->size[0];
      I->size[1] = A->size[0];
      emxEnsureCapacity((emxArray__common *)I, jrow, (int)sizeof(signed char));
      m = A->size[0] * A->size[0];
      for (jrow = 0; jrow < m; jrow++) {
        I->data[jrow] = 0;
      }

      if (A->size[0] > 0) {
        for (k = 0; k + 1 <= b_n; k++) {
          I->data[k + I->size[0] * k] = 1;
        }
      }

      jrow = V->size[0] * V->size[1];
      V->size[0] = I->size[0];
      V->size[1] = I->size[1];
      emxEnsureCapacity((emxArray__common *)V, jrow, (int)sizeof(creal_T));
      m = I->size[0] * I->size[1];
      for (jrow = 0; jrow < m; jrow++) {
        V->data[jrow].re = I->data[jrow];
        V->data[jrow].im = 0.0;
      }

      emxFree_int8_T(&I);
      if ((!(A->size[0] <= 1)) && (!(ihi < ilo + 2))) {
        for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
          for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
            b_A = A->data[(jrow + A->size[0] * jcol) - 1];
            c_A = A->data[jrow + A->size[0] * jcol];
            xzlartg(b_A, c_A, &c, &tmp, &A->data[(jrow + A->size[0] * jcol) - 1]);
            A->data[jrow + A->size[0] * jcol].re = 0.0;
            A->data[jrow + A->size[0] * jcol].im = 0.0;
            for (j = jcol + 1; j + 1 <= b_n; j++) {
              absxk = tmp.re * A->data[jrow + A->size[0] * j].re - tmp.im *
                A->data[jrow + A->size[0] * j].im;
              ctoc = tmp.re * A->data[jrow + A->size[0] * j].im + tmp.im *
                A->data[jrow + A->size[0] * j].re;
              stemp_re = c * A->data[(jrow + A->size[0] * j) - 1].re + absxk;
              stemp_im = c * A->data[(jrow + A->size[0] * j) - 1].im + ctoc;
              absxk = A->data[(jrow + A->size[0] * j) - 1].re;
              ctoc = A->data[(jrow + A->size[0] * j) - 1].im;
              cfrom1 = A->data[(jrow + A->size[0] * j) - 1].im;
              cto1 = A->data[(jrow + A->size[0] * j) - 1].re;
              A->data[jrow + A->size[0] * j].re = c * A->data[jrow + A->size[0] *
                j].re - (tmp.re * absxk + tmp.im * ctoc);
              A->data[jrow + A->size[0] * j].im = c * A->data[jrow + A->size[0] *
                j].im - (tmp.re * cfrom1 - tmp.im * cto1);
              A->data[(jrow + A->size[0] * j) - 1].re = stemp_re;
              A->data[(jrow + A->size[0] * j) - 1].im = stemp_im;
            }

            tmp.re = -tmp.re;
            tmp.im = -tmp.im;
            for (i = 0; i + 1 <= ihi; i++) {
              absxk = tmp.re * A->data[i + A->size[0] * (jrow - 1)].re - tmp.im *
                A->data[i + A->size[0] * (jrow - 1)].im;
              ctoc = tmp.re * A->data[i + A->size[0] * (jrow - 1)].im + tmp.im *
                A->data[i + A->size[0] * (jrow - 1)].re;
              stemp_re = c * A->data[i + A->size[0] * jrow].re + absxk;
              stemp_im = c * A->data[i + A->size[0] * jrow].im + ctoc;
              absxk = A->data[i + A->size[0] * jrow].re;
              ctoc = A->data[i + A->size[0] * jrow].im;
              cfrom1 = A->data[i + A->size[0] * jrow].im;
              cto1 = A->data[i + A->size[0] * jrow].re;
              A->data[i + A->size[0] * (jrow - 1)].re = c * A->data[i + A->size
                [0] * (jrow - 1)].re - (tmp.re * absxk + tmp.im * ctoc);
              A->data[i + A->size[0] * (jrow - 1)].im = c * A->data[i + A->size
                [0] * (jrow - 1)].im - (tmp.re * cfrom1 - tmp.im * cto1);
              A->data[i + A->size[0] * jrow].re = stemp_re;
              A->data[i + A->size[0] * jrow].im = stemp_im;
            }

            for (i = 0; i + 1 <= b_n; i++) {
              absxk = tmp.re * V->data[i + V->size[0] * (jrow - 1)].re - tmp.im *
                V->data[i + V->size[0] * (jrow - 1)].im;
              ctoc = tmp.re * V->data[i + V->size[0] * (jrow - 1)].im + tmp.im *
                V->data[i + V->size[0] * (jrow - 1)].re;
              stemp_re = c * V->data[i + V->size[0] * jrow].re + absxk;
              stemp_im = c * V->data[i + V->size[0] * jrow].im + ctoc;
              absxk = V->data[i + V->size[0] * jrow].re;
              ctoc = V->data[i + V->size[0] * jrow].im;
              cfrom1 = V->data[i + V->size[0] * jrow].im;
              cto1 = V->data[i + V->size[0] * jrow].re;
              V->data[i + V->size[0] * (jrow - 1)].re = c * V->data[i + V->size
                [0] * (jrow - 1)].re - (tmp.re * absxk + tmp.im * ctoc);
              V->data[i + V->size[0] * (jrow - 1)].im = c * V->data[i + V->size
                [0] * (jrow - 1)].im - (tmp.re * cfrom1 - tmp.im * cto1);
              V->data[i + V->size[0] * jrow].re = stemp_re;
              V->data[i + V->size[0] * jrow].im = stemp_im;
            }
          }
        }
      }

      xzhgeqz(A, ilo, ihi, V, &jcol, alpha1, beta1);
      *info = jcol;
      if (jcol == 0) {
        xztgevc(A, V);
        b_n = V->size[0];
        m = V->size[1];
        if (ilo > 1) {
          for (i = ilo - 2; i + 1 >= 1; i--) {
            k = rscale->data[i] - 1;
            if (rscale->data[i] != i + 1) {
              for (j = 0; j + 1 <= m; j++) {
                tmp = V->data[i + V->size[0] * j];
                V->data[i + V->size[0] * j] = V->data[k + V->size[0] * j];
                V->data[k + V->size[0] * j] = tmp;
              }
            }
          }
        }

        if (ihi < b_n) {
          while (ihi + 1 <= b_n) {
            k = rscale->data[ihi] - 1;
            if (rscale->data[ihi] != ihi + 1) {
              for (j = 0; j + 1 <= m; j++) {
                tmp = V->data[ihi + V->size[0] * j];
                V->data[ihi + V->size[0] * j] = V->data[k + V->size[0] * j];
                V->data[k + V->size[0] * j] = tmp;
              }
            }

            ihi++;
          }
        }

        for (jcol = 0; jcol < n; jcol++) {
          absxk = std::abs(V->data[V->size[0] * jcol].re) + std::abs(V->data
            [V->size[0] * jcol].im);
          if (n > 1) {
            for (k = 1; k - 1 <= n - 2; k++) {
              ctoc = std::abs(V->data[k + V->size[0] * jcol].re) + std::abs
                (V->data[k + V->size[0] * jcol].im);
              if (ctoc > absxk) {
                absxk = ctoc;
              }
            }
          }

          if (absxk >= 6.7178761075670888E-139) {
            absxk = 1.0 / absxk;
            for (k = 0; k < n; k++) {
              V->data[k + V->size[0] * jcol].re *= absxk;
              V->data[k + V->size[0] * jcol].im *= absxk;
            }
          }
        }

        if (ilascl) {
          notdone = true;
          while (notdone) {
            cfrom1 = anrmto * 2.0041683600089728E-292;
            cto1 = anrm / 4.9896007738368E+291;
            if ((cfrom1 > anrm) && (anrm != 0.0)) {
              stemp_im = 2.0041683600089728E-292;
              anrmto = cfrom1;
            } else if (cto1 > anrmto) {
              stemp_im = 4.9896007738368E+291;
              anrm = cto1;
            } else {
              stemp_im = anrm / anrmto;
              notdone = false;
            }

            jrow = alpha1->size[0];
            emxEnsureCapacity((emxArray__common *)alpha1, jrow, (int)sizeof
                              (creal_T));
            m = alpha1->size[0];
            for (jrow = 0; jrow < m; jrow++) {
              alpha1->data[jrow].re *= stemp_im;
              alpha1->data[jrow].im *= stemp_im;
            }
          }
        }
      }

      emxFree_int32_T(&rscale);
    }
  }
}

/* End of code generation (xzggev.cpp) */
