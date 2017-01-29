/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * CVfast_initialize.cpp
 *
 * Code generation for function 'CVfast_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "CVfast_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "CVfast_data.h"

/* Function Definitions */
void CVfast_initialize()
{
  rt_InitInfAndNaN(8U);
  //omp_init_nest_lock(&emlrtNestLockGlobal);
  c_eml_rand_mt19937ar_stateful_i();
}

/* End of code generation (CVfast_initialize.cpp) */
