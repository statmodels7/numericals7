#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include "n7_par.h"
using namespace Rcpp;

// The Student t distribution function T_n(x), with real degrees of freedom,
// written so that it never reaches the R API: R's pt() goes through pbeta's
// bratio(), which may signal a precision warning, and a warning raised from
// a worker thread kills the process. n7_pt() is R's pt() (src/nmath/pt.c,
// R 4.6.0) and n7_pbeta() its pbeta() (src/nmath/pbeta.c), over
// n7_bratio(), the copy of R's TOMS 708 in n7_toms708.c; the code is R's
// line for line except that the warnings are gone and the status bratio()
// reports is dropped, as pbeta_raw() drops it after warning. The values are
// therefore R's. The R sources are GPL (>= 2); see inst/COPYRIGHTS.

extern "C" void n7_bratio(double a, double b, double x, double y, double *w,
                          double *w1, int *ierr, int log_p);

namespace {

double n7_pbeta_raw(double x, double a, double b, int lower_tail,
                    int log_p) {
  const double d0 = log_p ? R_NegInf : 0.0, d1 = log_p ? 0.0 : 1.0;
  const double dt0 = lower_tail ? d0 : d1, dt1 = lower_tail ? d1 : d0;
  if (x >= 1) return dt1;
  if (a == 0 || b == 0 || !std::isfinite(a) || !std::isfinite(b)) {
    if (a == 0 && b == 0) return log_p ? -M_LN2 : 0.5;
    if (a == 0 || a / b == 0) return dt1;
    if (b == 0 || b / a == 0) return dt0;
    if (x < 0.5) return dt0;
    return dt1;
  }
  if (x <= 0) return dt0;
  double x1 = 0.5 - x + 0.5, w, wc;
  int ierr;
  n7_bratio(a, b, x, x1, &w, &wc, &ierr, log_p);
  return lower_tail ? w : wc;
}

double n7_pbeta(double x, double a, double b, int lower_tail, int log_p) {
  if (std::isnan(x) || std::isnan(a) || std::isnan(b)) return x + a + b;
  if (a < 0 || b < 0) return R_NaN;
  return n7_pbeta_raw(x, a, b, lower_tail, log_p);
}

}  // namespace

extern "C" {

// T_n(x) as R's pt(x, n, lower_tail, log_p) returns it
double n7_pt(double x, double n, int lower_tail, int log_p) {
  double val, nx;
  if (std::isnan(x) || std::isnan(n)) return x + n;
  if (n <= 0.0) return R_NaN;
  const double d0 = log_p ? R_NegInf : 0.0, d1 = log_p ? 0.0 : 1.0;
  if (!std::isfinite(x)) {
    const double dt0 = lower_tail ? d0 : d1, dt1 = lower_tail ? d1 : d0;
    return (x < 0) ? dt0 : dt1;
  }
  if (!std::isfinite(n)) return R::pnorm5(x, 0.0, 1.0, lower_tail, log_p);
  nx = 1 + (x / n) * x;
  if (nx > 1e100) {
    double lval;
    lval = -0.5 * n * (2 * std::log(std::fabs(x)) - std::log(n)) -
      R::lbeta(0.5 * n, 0.5) - std::log(0.5 * n);
    val = log_p ? lval : std::exp(lval);
  } else {
    val = (n > x * x)
      ? n7_pbeta(x * x / (n + x * x), 0.5, n / 2., /*lower_tail*/0, log_p)
      : n7_pbeta(1. / nx, n / 2., 0.5, /*lower_tail*/1, log_p);
  }
  if (x <= 0.) lower_tail = !lower_tail;
  if (log_p) {
    if (lower_tail) return std::log1p(-0.5 * std::exp(val));
    return val - M_LN2;
  }
  val /= 2.;
  return lower_tail ? (0.5 - val + 0.5) : val;
}

}  // extern "C"

// [[Rcpp::export]]
NumericVector student_t_cdf_cpp(NumericVector q, NumericVector df,
                                bool lower_tail, bool log_p, int threads = 1) {
  const R_xlen_t nq = q.size(), nd = df.size();
  const R_xlen_t n = (nq == 0 || nd == 0) ? 0 : std::max(nq, nd);
  NumericVector out(n);
  const double *qp = q.begin(), *dp = df.begin();
  double* op = out.begin();
  int lt = lower_tail, lp = log_p;
  n7::par_for(n, threads, n7::kMinCostly, [&](std::size_t i) {
    op[i] = n7_pt(qp[i % nq], dp[i % nd], lt, lp);
  });
  return out;
}

// [[Rcpp::init]]
void n7_register_student_t(DllInfo* dll) {
  R_RegisterCCallable("numericals7", "n7_pt", (DL_FUNC) n7_pt);
}
