/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * mean.h
 *
 * Code generation for function 'mean'
 *
 */

#ifndef MEAN_H
#define MEAN_H

/* Include files */
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "omp.h"
#include "CVfast_types.h"

/* Function Declarations */
extern void b_mean(const emxArray_real_T *x, emxArray_real_T *y);
extern void c_mean(const emxArray_real_T *x, emxArray_real_T *y);
extern double mean(const emxArray_real_T *x);

#endif

/* End of code generation (mean.h) */
