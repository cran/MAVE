/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sort1.h
 *
 * Code generation for function 'sort1'
 *
 */

#ifndef SORT1_H
#define SORT1_H

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
extern void b_sort(emxArray_real_T *x);
extern void d_sort(emxArray_real_T *x, emxArray_int32_T *idx);
extern void f_sort(emxArray_real_T *x);
extern void g_sort(emxArray_creal_T *x, emxArray_int32_T *idx);
extern void sort(emxArray_real_T *x, emxArray_int32_T *idx);

#endif

/* End of code generation (sort1.h) */
