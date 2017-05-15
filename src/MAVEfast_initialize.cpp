/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * MAVEfast_initialize.cpp
 *
 * Code generation for function 'MAVEfast_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "MAVEfast_initialize.h"
#include "MAVEfast_data.h"

/* Function Definitions */
void MAVEfast_initialize()
{
  rt_InitInfAndNaN(8U);
//  omp_init_nest_lock(&emlrtNestLockGlobal);
}

/* End of code generation (MAVEfast_initialize.cpp) */
