/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzgetrf.h
 *
 * Code generation for function 'xzgetrf'
 *
 */

#ifndef XZGETRF_H
#define XZGETRF_H

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
extern void xzgetrf(int m, int n, emxArray_real_T *A, int lda, emxArray_int32_T *
                    ipiv, int *info);

#endif

/* End of code generation (xzgetrf.h) */
