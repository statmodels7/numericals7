/* What R's src/nmath/nmath.h and dpq.h give toms708.c and bessel_k.c, for
 * the copies of them in n7_toms708.c and n7_bessel_k.c. The warnings R raises from these sources are compiled to
 * nothing here: a warning calls into the R API, and the point of the copy
 * is a distribution function a worker thread may call. The status codes the
 * algorithm reports are kept and read by the caller. Adapted from R
 * (R Core Team, GPL (>= 2)); see inst/COPYRIGHTS. */
#ifndef NUMERICALS7_N7_NMATH_H
#define NUMERICALS7_N7_NMATH_H

#include <math.h>
#include <float.h>
#include <limits.h>
#include <Rmath.h>
#include <R_ext/Boolean.h>

#ifndef IEEE_754
#define IEEE_754 1
#endif

#define ML_POSINF (1.0 / 0.0)
#define ML_NEGINF ((-1.0) / 0.0)
#define ML_NAN (0.0 / 0.0)
#define ISNAN(x) (isnan(x) != 0)
#define R_FINITE(x) isfinite(x)

#define ML_WARN_return_NAN { return ML_NAN; }
#define ML_WARNING(x, s) ((void) 0)
#define _(String) (String)
#define MATHLIB_WARNING(fmt, x) ((void) 0)
#define MATHLIB_WARNING2(fmt, x, x2) ((void) 0)
#define MATHLIB_WARNING3(fmt, x, x2, x3) ((void) 0)
#define MATHLIB_WARNING4(fmt, x, x2, x3, x4) ((void) 0)
#define MATHLIB_WARNING5(fmt, x, x2, x3, x4, x5) ((void) 0)
#define MATHLIB_WARNING6(fmt, x, x2, x3, x4, x5, x6) ((void) 0)
#define REprintf(...) ((void) 0)

#define attribute_hidden

#define R_D__0 (log_p ? ML_NEGINF : 0.)
#define R_D__1 (log_p ? 0. : 1.)
#define R_DT_0 (lower_tail ? R_D__0 : R_D__1)
#define R_DT_1 (lower_tail ? R_D__1 : R_D__0)
#define R_D_Cval(p) (lower_tail ? (0.5 - (p) + 0.5) : (p))
#define R_D_exp(x) (log_p ? (x) : exp(x))

/* R's d1mach() and i1mach() for IEEE double (src/nmath/d1mach.c,
 * i1mach.c), at the arguments toms708.c reads */
static inline double n7_d1mach(int i) {
  switch (i) {
  case 1: return DBL_MIN;
  case 2: return DBL_MAX;
  case 3: return 0.5 * DBL_EPSILON;
  case 4: return DBL_EPSILON;
  case 5: return M_LOG10_2;
  default: return 0.0;
  }
}
static inline int n7_i1mach(int i) {
  switch (i) {
  case 15: return DBL_MIN_EXP;
  case 16: return DBL_MAX_EXP;
  default: return 0;
  }
}
#define Rf_d1mach n7_d1mach
#define Rf_i1mach n7_i1mach

#ifdef __cplusplus
extern "C" {
#endif
void n7_bratio(double a, double b, double x, double y, double *w, double *w1,
               int *ierr, int log_p);
double n7_bessel_k_ex(double x, double alpha, double expo, double *bk);
#ifdef __cplusplus
}
#endif

#endif
