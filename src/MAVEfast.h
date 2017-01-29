/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * MAVEfast.h
 *
 * Code generation for function 'MAVEfast'
 *
 */

#ifndef MAVEFAST_H
#define MAVEFAST_H

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
extern void MAVEfast(emxArray_real_T *x, const emxArray_real_T *y, const
                     emxArray_char_T *method, emxArray_real_T *BB1D,
                     emxArray_real_T *ky);

#endif

/* End of code generation (MAVEfast.h) */
