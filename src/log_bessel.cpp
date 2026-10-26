// The compiled kernels behind log_bessel_i() and log_bessel_k(): the same
// branches and formulas as the R twins kept beside them, one scalar loop
// per element. The u_k polynomial coefficients are injected once from the
// R table at load, so the two routes share the table and nothing else.
#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "n7_par.h"
#include <cmath>
using namespace Rcpp;

static const int UK_MAX = 13;
// u_k polynomial coefficients, filled from R once via lb_set_uk_cpp(). The
// exponents of u_k are k, k + 2, ..., 3k, so only the coefficients are kept
// and lb_u_sum() evaluates u_k(t) = t^k P_k(t^2) by Horner; the setter
// checks the exponents, so a table of another shape fails at load.
static std::vector<std::vector<double>> UK_C;

// [[Rcpp::export]]
void lb_set_uk_cpp(List uk) {
  UK_C.clear();
  for (int k = 0; k < uk.size(); ++k) {
    List el = uk[k];
    NumericVector e = el["e"], c = el["c"];
    if (e.size() != c.size()) stop("u_%d: exponents and coefficients differ in length", k + 1);
    for (R_xlen_t j = 0; j < e.size(); ++j) {
      if (e[j] != (double) (k + 1 + 2 * j)) stop("u_%d: exponent %d is not %d", k + 1,
                                                 (int) j + 1, (int) (k + 1 + 2 * j));
    }
    UK_C.push_back(std::vector<double>(c.begin(), c.end()));
  }
}

static inline int branch(double x, double v) {
  double lx = std::log(x), lv = std::log(v);
  if ((x > 1400 && v < 3.05) ||
      ((0.6229 * lx - 3.2318) > lv && v > 3.1 && v * v < 0.002 * x)) return 1;
  if ((x > 30 && v < 15.3919) ||
      ((0.5113 * lx + 0.7939) > lv && x > 59.6925 && v * v < 2 * x)) return 2;
  if ((x > 274.2377 && v > 0.3) || (v > 163.6993)) return 3;
  if ((x > 84.4153 && v > 0.46) || (v > 56.9971)) return 4;
  if ((x > 35.9074 && v > 0.6) || (v > 20.1534)) return 5;
  if ((x > 19.6931 && v > 0.7) || (v > 12.6964)) return 6;
  return 0;
}

static double lb_mu(double x, double v, int K, double sign) {
  double mu = 4 * v * v, s = 1.0, ck = 1.0;
  for (int k = 1; k <= K; ++k) {
    ck *= sign * (mu - (2.0 * k - 1) * (2.0 * k - 1)) / (k * 8.0 * x);
    s += ck;
  }
  return std::log(std::fabs(s));
}

static double lb_u_sum(double t, double v, int K, double sign) {
  double s = 1.0, vk = 1.0, t2 = t * t, tk = 1.0;
  for (int k = 1; k <= K; ++k) {
    vk *= sign / v;
    tk *= t;
    const std::vector<double>& C = UK_C[k - 1];
    double p = C.back();
    for (int j = (int) C.size() - 2; j >= 0; --j) p = p * t2 + C[j];
    s += vk * (tk * p);
  }
  return std::log(std::fabs(s));
}

// log I by its ascending series, sum_k t_k with t_0 = 1 and
// t_k = t_{k-1} q / (k (k + v)), q = x^2 / 4, scaled by (x/2)^v / Gamma(v + 1).
// The sum is 0F1(; v + 1; q), which decreases in v, so on this branch
// (x <= 30) it is at most I_0(30) < 1e12 and needs no log scale. Past the
// peak (q < k (k + v)) every term is smaller than the one before, so the
// first term that leaves the sum unchanged ends the loop with the bits that
// all 90 terms would give.
static double lbi_series(double x, double v) {
  const int nmax = 90;
  double q = 0.25 * x * x, t = 1.0, s = 1.0;
  for (int k = 1; k <= nmax; ++k) {
    double den = k * (k + v);
    t *= q / den;
    if (q < den && s + t == s) break;
    s += t;
  }
  return v * std::log(x / 2) - std::lgamma(v + 1.0) + std::log(s);
}

static double lbi_one(double x, double v);

static double lbk_integral(double x, double v) {
  const int N = 600;
  const double n8 = 8.0;
  double beta = 2 * n8 / (2 * v + 1);
  double la[2 * N], best = -1e308;
  int m = 0;
  for (int k = 1; k < N; ++k) {
    double u = (double) k / N, w = (k % 2) ? 4.0 : 2.0;
    double ub = std::pow(u, beta);
    la[m] = std::log(w) + std::log(beta) - ub +
      (v - 0.5) * std::log(2 * x + ub) + (n8 - 1) * std::log(u);
    if (la[m] > best) best = la[m];
    ++m;
    la[m] = std::log(w) - 1 / u - (2 * v + 1) * std::log(u) +
      (v - 0.5) * std::log(2 * x * u + 1);
    if (la[m] > best) best = la[m];
    ++m;
  }
  double g1 = std::log(beta) - 1 + (v - 0.5) * std::log(2 * x + 1);
  double h1 = -1 + (v - 0.5) * std::log(2 * x + 1);
  double s = std::exp(g1 - best) + std::exp(h1 - best);
  for (int j = 0; j < m; ++j) s += std::exp(la[j] - best);
  double lint = best + std::log(s) - std::log(3.0 * N);
  return 0.5 * std::log(M_PI) - std::lgamma(v + 0.5) - v * std::log(2 * x) -
    x + lint;
}

static inline void u_parts(double x, double v, double& t, double& eta,
                           double& lr) {
  double xp = x / v, r = std::sqrt(1 + xp * xp);
  t = 1 / r;
  eta = r + std::log(xp / (1 + r));
  lr = std::log1p(xp * xp);
}

static double lbi_one(double x, double v) {
  if (x == 0) return v == 0 ? 0.0 : R_NegInf;
  int b = branch(x, v);
  double t, eta, lr;
  switch (b) {
  case 1: return x - 0.5 * std::log(2 * M_PI * x) + lb_mu(x, v, 3, -1);
  case 2: return x - 0.5 * std::log(2 * M_PI * x) + lb_mu(x, v, 20, -1);
  case 3: case 4: case 5: case 6: {
    int K = (b == 3) ? 4 : (b == 4) ? 6 : (b == 5) ? 9 : 13;
    u_parts(x, v, t, eta, lr);
    return -0.5 * std::log(2 * M_PI * v) + v * eta - 0.25 * lr +
      lb_u_sum(t, v, K, +1);
  }
  default: return lbi_series(x, v);
  }
}

static double lbk_one(double x, double v) {
  if (x == 0) return R_PosInf;
  v = std::fabs(v);
  int b = branch(x, v);
  double t, eta, lr;
  switch (b) {
  case 1: return 0.5 * (std::log(M_PI) - std::log(2 * x)) - x + lb_mu(x, v, 3, +1);
  case 2: return 0.5 * (std::log(M_PI) - std::log(2 * x)) - x + lb_mu(x, v, 20, +1);
  case 3: case 4: case 5: case 6: {
    int K = (b == 3) ? 4 : (b == 4) ? 6 : (b == 5) ? 9 : 13;
    u_parts(x, v, t, eta, lr);
    return 0.5 * (std::log(M_PI) - std::log(2 * v)) - v * eta - 0.25 * lr +
      lb_u_sum(t, v, K, -1);
  }
  default: {
    // hybrid: R's scaled besselK where it does not overflow
    if (v * std::log(2 / std::max(x, 1e-300)) +
        std::lgamma(std::max(v, 0.5)) <= 690) {
      return std::log(R::bessel_k(x, v, 2)) - x;
    }
    return lbk_integral(x, v);
  }
  }
}

// Every branch of lbi_one() is this file's own arithmetic -- a series or a
// uniform asymptotic expansion -- with no call into Rmath, which is what
// makes the body admissible inside a worker. Since 0.24.0 every branch costs
// about 0.2 microseconds an element (3 to 4 ms over 20000 points), and the
// threaded region pays from about 500 elements (measured: 0.5x at 128,
// 1.2x at 512, 2.6x at 2048 on 16 threads). The kernel still takes the
// shared threshold kMinCostly (128), which was set for the earlier cost.
// [[Rcpp::export]]
NumericVector log_bessel_i_cpp(NumericVector x, NumericVector nu,
                               int threads = 1) {
  R_xlen_t n = std::max(x.size(), nu.size());
  NumericVector out(n);
  const R_xlen_t nx = x.size(), nv = nu.size();
  const double* xp = x.begin();
  const double* vp = nu.begin();
  double* op = out.begin();
  n7::par_for((std::size_t) n, threads, n7::kMinCostly, [&](std::size_t i) {
    double xi = xp[i % nx], vi = vp[i % nv];
    op[i] = (ISNAN(xi) || ISNAN(vi) || xi < 0 || vi < 0)
      ? NA_REAL : lbi_one(xi, vi);
  });
  return out;
}

// log I_nu(x) for one argument, as log_bessel_i() returns it, for the
// compiled code of other packages: R_GetCCallable("numericals7",
// "n7_log_bessel_i"). The u_k table it reads is filled when the package is
// loaded, so a call from a worker thread reads it only.
extern "C" double n7_log_bessel_i(double x, double nu) {
  return (ISNAN(x) || ISNAN(nu) || x < 0 || nu < 0) ? NA_REAL : lbi_one(x, nu);
}

// the entry point resolved as other packages resolve it, for the tests
// [[Rcpp::export]]
NumericVector n7_log_bessel_i_probe(NumericVector x, NumericVector nu) {
  typedef double (*Fn)(double, double);
  Fn f = (Fn) R_GetCCallable("numericals7", "n7_log_bessel_i");
  R_xlen_t n = std::max(x.size(), nu.size());
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) out[i] = f(x[i % x.size()], nu[i % nu.size()]);
  return out;
}

// [[Rcpp::init]]
void n7_register_log_bessel(DllInfo* dll) {
  R_RegisterCCallable("numericals7", "n7_log_bessel_i", (DL_FUNC) n7_log_bessel_i);
}

// This one is NOT threaded and takes no count: lbk_one()'s hybrid branch
// calls R::bessel_k, which can raise a warning, and a warning from a worker
// thread kills the process (n7_par.h). It costs 4 to 7 ms over 20000 points,
// against log I's 3 to 4, so what the restriction costs is small, and
// stating it is better than an argument that
// would be read by nobody.
// [[Rcpp::export]]
NumericVector log_bessel_k_cpp(NumericVector x, NumericVector nu) {
  R_xlen_t n = std::max(x.size(), nu.size());
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = x[i % x.size()], vi = nu[i % nu.size()];
    out[i] = (ISNAN(xi) || ISNAN(vi) || xi < 0) ? NA_REAL : lbk_one(xi, vi);
  }
  return out;
}
