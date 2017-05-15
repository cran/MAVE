/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xtrsm.h
 *
 * Code generation for function 'xtrsm'
 *
 */

#ifndef XTRSM_H
#define XTRSM_H

/* Include files */
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
//#include "omp.h"
#include "MAVEfast_types.h"

/* Function Declarations */
extern void xtrsm(int m, int n, const emxArray_real_T *A, int lda,
                  emxArray_real_T *B, int ldb);

#endif

/* End of code generation (xtrsm.h) */
