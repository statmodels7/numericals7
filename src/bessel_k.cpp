#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <cstdlib>
#include "n7_par.h"
using namespace Rcpp;

// The modified Bessel function of the third kind K_alpha(x), written so that
// it never reaches the R API: R's bessel_k() allocates its work array with
// R_alloc() and may signal a warning, and neither may happen on a worker
// thread. n7_bessel_k() is R's bessel_k() (src/nmath/bessel_k.c, R 4.6.0)
// over n7_bessel_k_ex(), the copy in n7_bessel_k.c, with the work array on
// the stack, or from malloc() beyond 64 orders; the values are therefore
// R's. expo is 1 for K_alpha(x) and 2 for exp(x) K_alpha(x), as in R. The
// R sources are GPL (>= 2); see inst/COPYRIGHTS.

extern "C" double n7_bessel_k_ex(double x, double alpha, double expo,
                                 double *bk);

extern "C" {

// K_alpha(x) as R's bessel_k(x, alpha, expo) returns it
double n7_bessel_k(double x, double alpha, double expo) {
  if (std::isnan(x) || std::isnan(alpha)) return x + alpha;
  if (x < 0) return R_NaN;
  const double a = std::fabs(alpha);
  if (a >= 2147483646.0) return R_NaN;
  const int nb = 1 + (int) std::floor(a);
  if (nb <= 64) {
    double bk[64];
    return n7_bessel_k_ex(x, alpha, expo, bk);
  }
  double *bk = (double *) std::malloc((std::size_t) nb * sizeof(double));
  if (bk == nullptr) return R_NaN;
  const double val = n7_bessel_k_ex(x, alpha, expo, bk);
  std::free(bk);
  return val;
}

}  // extern "C"

// [[Rcpp::export]]
NumericVector bessel_k_cpp(NumericVector x, NumericVector nu, bool scaled,
                           int threads = 1) {
  const R_xlen_t nx = x.size(), nn = nu.size();
  const R_xlen_t n = (nx == 0 || nn == 0) ? 0 : std::max(nx, nn);
  NumericVector out(n);
  const double *xp = x.begin(), *np = nu.begin();
  double* op = out.begin();
  const double expo = scaled ? 2.0 : 1.0;
  n7::par_for(n, threads, n7::kMinCostly, [&](std::size_t i) {
    op[i] = n7_bessel_k(xp[i % nx], np[i % nn], expo);
  });
  return out;
}

// [[Rcpp::init]]
void n7_register_bessel_k(DllInfo* dll) {
  R_RegisterCCallable("numericals7", "n7_bessel_k", (DL_FUNC) n7_bessel_k);
}
