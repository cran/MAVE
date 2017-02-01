/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * CVfast.h
 *
 * Code generation for function 'CVfast'
 *
 */

#ifndef CVFAST_H
#define CVFAST_H

/* Include files */
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
//#include "omp.h"
#include "CVfast_types.h"

/* Function Declarations */
extern void CVfast(const emxArray_real_T *x, const emxArray_real_T *ky, const
                   emxArray_real_T *BB1D, emxArray_real_T *cv);

#endif

/* End of code generation (CVfast.h) */
