#ifdef __cplusplus
extern "C" {
#endif

// include guard
#ifndef QPALM_FORTRAN_H
#define QPALM_FORTRAN_H

// required packages
#include <stdint.h>
#include <qpalm.h>

#ifdef QPALM_FORTRAN_64BIT_INDICES
typedef int64_t f_int;
#else
typedef int32_t f_int;
#endif

#ifdef QPALM_FORTRAN_SINGLE_PRECISION
typedef float f_float;
#else
typedef double f_float;
#endif

void qpalm_fortran_c( f_int n,
                      f_int m,
                      f_int h_ne,
                      f_int H_ptr[],
                      f_int H_row[],
                      f_float H_val[],
                      f_float g[n],
                      f_float f,
                      f_int a_ne,
                      f_int A_ptr[],
                      f_int A_row[],
                      f_float A_val[],
                      f_float c_l[],
                      f_float c_u[],
                      QPALMSettings settings_c,
                      f_float x[],
                      f_float y[],
                      QPALMInfo *info_c );

// end include guard
#endif

#ifdef __cplusplus
} /* extern "C" */
#endif
