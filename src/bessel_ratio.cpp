#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <limits>
#include "n7_par.h"
#include "n7_dd.h"
#include "bessel_ratio_coefs.h"
using namespace Rcpp;

// The ratio A(k) = I1(k)/I0(k), its derivatives in k, and the inverse map
// k(rho) with its derivatives, one function per order. Three regimes:
//
//  k < 0.5     the power series of A at 0, differentiated term by term. The
//              Riccati identity A' = 1 - A/k - A^2 cancels there: its terms
//              are of order k^-n at derivative n while the result is of
//              order one or k, so the identity in double loses everything
//              below about k = 1e-5 at the fourth order.
//  0.5 <= k < 30
//              A from its continued fraction r_j = k/(2j + k r_{j+1}),
//              truncated at j = floor(k) + 20: backwards in double for the
//              value, forwards by its convergents in double-double for the
//              derivatives, which then follow from the Riccati identity in
//              double-double. At k = 20 the identity in double loses 4e-10
//              at the fourth order (A' is of order k^-2 against terms of
//              order one); 32 digits absorb that.
//  k >= 30     the asymptotic series of A in 1/k (the quotient of those of
//              I1 and I0, 30 terms), differentiated term by term.
//
// Measured against 150-digit values (stabilita/bessel_ratio_measure.R) the
// three regimes reach the last bit of A and of every derivative; the
// thresholds are where the series still do. Nothing here touches the R API,
// so every kernel is threaded and every scalar function is also a C entry
// point for the packages above (registered at the end of the file).

namespace {

const int BR_NA = sizeof(BR_A) / sizeof(double);
const int BR_NQ = sizeof(BR_Q) / sizeof(double);
const double BR_SMALL = 0.5, BR_LARGE = 30.0;

// order M of sum_n BR_A[n] k^(2n+1), only that order
template <int M>
inline double br_small(double k) {
  const int nmin = M / 2;               // the first power 2n+1 >= M
  const double k2 = k * k;
  double acc = 0.0;
  for (int n = BR_NA - 1; n >= nmin; --n) {
    const int p = 2 * n + 1;
    double fac = 1.0;
    for (int j = 0; j < M; ++j) fac *= (double) (p - j);
    acc = acc * k2 + BR_A[n] * fac;
  }
  return (2 * nmin + 1 - M) ? acc * k : acc;
}

// order M of sum_n BR_Q[n] k^(-n), only that order:
// d^M k^-n = (-1)^M n (n+1) ... (n+M-1) k^(-n-M)
template <int M>
inline double br_large(double k) {
  const double u = 1.0 / k;
  double acc = 0.0;
  for (int n = BR_NQ - 1; n >= 0; --n) {
    double fac = 1.0;
    for (int j = 0; j < M; ++j) fac *= (double) (n + j);
    acc = acc * u + BR_Q[n] * fac;
  }
  double um = 1.0;
  for (int j = 0; j < M; ++j) um *= u;
  return (M % 2 ? -acc : acc) * um;
}

// 1 - A from the asymptotic series, without forming A
inline double br_large_complement(double k) {
  const double u = 1.0 / k;
  double acc = 0.0;
  for (int n = BR_NQ - 1; n >= 1; --n) acc = acc * u - BR_Q[n];
  return acc * u;
}

// A by the continued fraction, in double and in double-double
inline double br_cf(double k) {
  double r = 0.0;
  for (int j = (int) k + 20; j >= 1; --j) r = k / (2.0 * j + k * r);
  return r;
}
// In double-double the same truncated fraction is evaluated forwards, by
// the Wallis recurrences of its convergents P_j/Q_j: with x = k/2 it reads
// x/(1 + x^2/(2 + x^2/(3 + ...))), so P_j = j P_{j-1} + x^2 P_{j-2} and
// likewise Q_j. Every term is positive, so the recurrences are stable, and
// they cost multiplications and a single division, where the backward form
// costs one double-double division per step (four times the time).
inline n7::dd br_cf_dd(double k) {
  const double x = 0.5 * k;
  const n7::dd x2 = n7::dd_two_prod(x, x);
  n7::dd p0(0.0), p1(x), q0(1.0), q1(1.0);    // P_0, P_1, Q_0, Q_1
  const int n = (int) k + 20;
  for (int j = 2; j <= n; ++j) {
    const n7::dd p2 = (double) j * p1 + x2 * p0;
    const n7::dd q2 = (double) j * q1 + x2 * q0;
    p0 = p1; p1 = p2; q0 = q1; q1 = q2;
  }
  return p1 / q1;
}

// the Riccati identity differentiated, each order from the ones below it
inline n7::dd br_rd1(const n7::dd& A, const n7::dd& ik) {
  return 1.0 - A * ik - A * A;
}
inline n7::dd br_rd2(const n7::dd& A, const n7::dd& d1, const n7::dd& ik) {
  return -(d1 * ik) + A * ik * ik - 2.0 * A * d1;
}
inline n7::dd br_rd3(const n7::dd& A, const n7::dd& d1, const n7::dd& d2,
                     const n7::dd& ik) {
  const n7::dd ik2 = ik * ik;
  return -(d2 * ik) + 2.0 * d1 * ik2 - 2.0 * A * ik2 * ik - 2.0 * d1 * d1 -
    2.0 * A * d2;
}
inline n7::dd br_rd4(const n7::dd& A, const n7::dd& d1, const n7::dd& d2,
                     const n7::dd& d3, const n7::dd& ik) {
  const n7::dd ik2 = ik * ik;
  return -(d3 * ik) + 3.0 * d2 * ik2 - 6.0 * d1 * ik2 * ik +
    6.0 * A * ik2 * ik2 - 6.0 * d1 * d2 - 2.0 * A * d3;
}

inline bool br_bad(double k) { return ISNAN(k) || k < 0.0; }
inline double br_badval(double k) { return ISNAN(k) ? k : R_NaN; }

}  // namespace

// --- the scalar functions, one per order ----------------------------------

extern "C" {

double n7_bessel_ratio(double k) {
  if (br_bad(k)) return br_badval(k);
  if (k < BR_SMALL) return br_small<0>(k);
  if (k < BR_LARGE) return br_cf(k);
  return br_large<0>(k);
}

double n7_bessel_ratio_d1(double k) {
  if (br_bad(k)) return br_badval(k);
  if (k < BR_SMALL) return br_small<1>(k);
  if (k < BR_LARGE) {
    const n7::dd A = br_cf_dd(k), ik = n7::dd(1.0) / n7::dd(k);
    return br_rd1(A, ik).hi;
  }
  return br_large<1>(k);
}

double n7_bessel_ratio_d2(double k) {
  if (br_bad(k)) return br_badval(k);
  if (k < BR_SMALL) return br_small<2>(k);
  if (k < BR_LARGE) {
    const n7::dd A = br_cf_dd(k), ik = n7::dd(1.0) / n7::dd(k);
    const n7::dd d1 = br_rd1(A, ik);
    return br_rd2(A, d1, ik).hi;
  }
  return br_large<2>(k);
}

double n7_bessel_ratio_d3(double k) {
  if (br_bad(k)) return br_badval(k);
  if (k < BR_SMALL) return br_small<3>(k);
  if (k < BR_LARGE) {
    const n7::dd A = br_cf_dd(k), ik = n7::dd(1.0) / n7::dd(k);
    const n7::dd d1 = br_rd1(A, ik), d2 = br_rd2(A, d1, ik);
    return br_rd3(A, d1, d2, ik).hi;
  }
  return br_large<3>(k);
}

double n7_bessel_ratio_d4(double k) {
  if (br_bad(k)) return br_badval(k);
  if (k < BR_SMALL) return br_small<4>(k);
  if (k < BR_LARGE) {
    const n7::dd A = br_cf_dd(k), ik = n7::dd(1.0) / n7::dd(k);
    const n7::dd d1 = br_rd1(A, ik), d2 = br_rd2(A, d1, ik),
      d3 = br_rd3(A, d1, d2, ik);
    return br_rd4(A, d1, d2, d3, ik).hi;
  }
  return br_large<4>(k);
}

// A and its derivatives to order n (0 <= n <= 4) into out[0..n], for a
// caller that uses every one of them: one continued fraction serves them all
void n7_bessel_ratio_upto(double k, int n, double* out) {
  if (br_bad(k)) {
    for (int m = 0; m <= n; ++m) out[m] = br_badval(k);
    return;
  }
  if (k < BR_SMALL || k >= BR_LARGE) {
    const bool s = k < BR_SMALL;
    out[0] = s ? br_small<0>(k) : br_large<0>(k);
    if (n >= 1) out[1] = s ? br_small<1>(k) : br_large<1>(k);
    if (n >= 2) out[2] = s ? br_small<2>(k) : br_large<2>(k);
    if (n >= 3) out[3] = s ? br_small<3>(k) : br_large<3>(k);
    if (n >= 4) out[4] = s ? br_small<4>(k) : br_large<4>(k);
    return;
  }
  const n7::dd A = br_cf_dd(k), ik = n7::dd(1.0) / n7::dd(k);
  out[0] = A.hi;
  if (n < 1) return;
  const n7::dd d1 = br_rd1(A, ik);
  out[1] = d1.hi;
  if (n < 2) return;
  const n7::dd d2 = br_rd2(A, d1, ik);
  out[2] = d2.hi;
  if (n < 3) return;
  const n7::dd d3 = br_rd3(A, d1, d2, ik);
  out[3] = d3.hi;
  if (n < 4) return;
  out[4] = br_rd4(A, d1, d2, d3, ik).hi;
}

// k = A^-1(rho) by Newton's method. A is increasing and concave, so after
// the first step every iterate lies left of the root and rises to it; a step
// below 2 rho is replaced by it (A(k) < k/2 puts the root above). The
// residual A(k) - rho is formed as (1 - rho) - (1 - A(k)) when rho >= 1/2,
// where 1 - rho is exact and 1 - A comes without cancellation, so the
// iteration resolves k to its last bits up to the conditioning of rho
// itself, dk/k = -d(1 - rho)/(1 - rho) for k large.
double n7_bessel_ratio_inverse(double rho) {
  if (ISNAN(rho)) return rho;
  if (!(rho > 0.0 && rho < 1.0)) return NA_REAL;
  const double r = rho;
  double g = r < 0.53 ? 2.0 * r + r * r * r + 5.0 * std::pow(r, 5) / 6.0
    : (r < 0.85 ? -0.4 + 1.39 * r + 0.43 / (1.0 - r)
       : 1.0 / (r * r * r - 4.0 * r * r + 3.0 * r));
  const double lo = 2.0 * r;
  if (!(g >= lo)) g = lo;
  const bool upper = r >= 0.5;
  const double c = 1.0 - r;
  for (int it = 0; it < 200; ++it) {
    const double kk = g;
    double res, d1;
    if (kk < BR_SMALL) {
      res = br_small<0>(kk) - r;
      d1 = br_small<1>(kk);
    } else if (kk < BR_LARGE) {
      const n7::dd A = br_cf_dd(kk), ik = n7::dd(1.0) / n7::dd(kk);
      res = upper ? (c - (1.0 - A)).hi : (A + (-r)).hi;
      d1 = br_rd1(A, ik).hi;
    } else {
      res = upper ? c - br_large_complement(kk) : br_large<0>(kk) - r;
      d1 = br_large<1>(kk);
    }
    double kn = kk - res / d1;
    if (!std::isfinite(kn) || kn < lo) kn = lo;
    g = kn;
    if (std::fabs(kn - kk) <= 4.0 * std::numeric_limits<double>::epsilon() * kn ||
        (it > 0 && kn <= kk)) break;
  }
  return g;
}

// The derivatives of the inverse in rho, at rho = A(k), by the inverse
// function rule; A' > 0 keeps every denominator away from zero
double n7_bessel_ratio_inverse_d1(double k) {
  if (br_bad(k)) return br_badval(k);
  return 1.0 / n7_bessel_ratio_d1(k);
}

double n7_bessel_ratio_inverse_d2(double k) {
  if (br_bad(k)) return br_badval(k);
  double a[3];
  n7_bessel_ratio_upto(k, 2, a);
  return -a[2] / (a[1] * a[1] * a[1]);
}

double n7_bessel_ratio_inverse_d3(double k) {
  if (br_bad(k)) return br_badval(k);
  double a[4];
  n7_bessel_ratio_upto(k, 3, a);
  const double p1 = a[1], p2 = a[2], p3 = a[3];
  const double p12 = p1 * p1;
  return (3.0 * p2 * p2 - p1 * p3) / (p12 * p12 * p1);
}

double n7_bessel_ratio_inverse_d4(double k) {
  if (br_bad(k)) return br_badval(k);
  double a[5];
  n7_bessel_ratio_upto(k, 4, a);
  const double p1 = a[1], p2 = a[2], p3 = a[3], p4 = a[4];
  const double p12 = p1 * p1;
  return (-15.0 * p2 * p2 * p2 + 10.0 * p1 * p2 * p3 - p12 * p4) /
    (p12 * p12 * p12 * p1);
}

}  // extern "C"

// --- vectorized kernels ----------------------------------------------------

namespace {
template <typename F>
NumericVector br_map(NumericVector x, int threads, F f) {
  const std::size_t n = x.size();
  NumericVector out(n);
  const double* xp = x.begin();
  double* op = out.begin();
  n7::par_for(n, threads, n7::kMinCostly,
              [&](std::size_t i) { op[i] = f(xp[i]); });
  return out;
}
}  // namespace

// [[Rcpp::export]]
NumericVector bessel_ratio_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_d1_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_d1);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_d2_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_d2);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_d3_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_d3);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_d4_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_d4);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_inverse_cpp(NumericVector rho, int threads = 1) {
  return br_map(rho, threads, n7_bessel_ratio_inverse);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_inverse_d1_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_inverse_d1);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_inverse_d2_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_inverse_d2);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_inverse_d3_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_inverse_d3);
}
// [[Rcpp::export]]
NumericVector bessel_ratio_inverse_d4_cpp(NumericVector kappa, int threads = 1) {
  return br_map(kappa, threads, n7_bessel_ratio_inverse_d4);
}

// The scalar functions as C entry points: a consumer resolves them once with
// R_GetCCallable("numericals7", ...) and calls plain function pointers in
// its loops, threaded or not.
// [[Rcpp::init]]
void n7_register_bessel_ratio(DllInfo* dll) {
  R_RegisterCCallable("numericals7", "n7_bessel_ratio", (DL_FUNC) n7_bessel_ratio);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_d1", (DL_FUNC) n7_bessel_ratio_d1);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_d2", (DL_FUNC) n7_bessel_ratio_d2);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_d3", (DL_FUNC) n7_bessel_ratio_d3);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_d4", (DL_FUNC) n7_bessel_ratio_d4);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_upto", (DL_FUNC) n7_bessel_ratio_upto);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_inverse",
                      (DL_FUNC) n7_bessel_ratio_inverse);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_inverse_d1",
                      (DL_FUNC) n7_bessel_ratio_inverse_d1);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_inverse_d2",
                      (DL_FUNC) n7_bessel_ratio_inverse_d2);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_inverse_d3",
                      (DL_FUNC) n7_bessel_ratio_inverse_d3);
  R_RegisterCCallable("numericals7", "n7_bessel_ratio_inverse_d4",
                      (DL_FUNC) n7_bessel_ratio_inverse_d4);
}
