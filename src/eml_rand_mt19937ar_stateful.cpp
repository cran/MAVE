/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * eml_rand_mt19937ar_stateful.cpp
 *
 * Code generation for function 'eml_rand_mt19937ar_stateful'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CVfast.h"
#include "MAVEfast.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "CVfast_data.h"

/* Function Definitions */
void c_eml_rand_mt19937ar_stateful_i()
{
  unsigned int r;
  int mti;
  memset(&state[0], 0, 625U * sizeof(unsigned int));
  r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = (r ^ r >> 30U) * 1812433253U + (1 + mti);
    state[mti + 1] = r;
  }

  state[624] = 624U;
}

/* End of code generation (eml_rand_mt19937ar_stateful.cpp) */
