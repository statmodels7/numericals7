// The compiled kernels behind log_bessel_i() and log_bessel_k(): the same
// branches and formulas as the R twins kept beside them, one scalar loop
// per element. The u_k polynomial coefficients are injected once from the
// R table at load, so the two routes share the table and nothing else.
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static const int UK_MAX = 13;
// u_k polynomial coefficients, filled from R once via lb_set_uk_cpp()
static std::vector<std::vector<double>> UK_E, UK_C;

// [[Rcpp::export]]
void lb_set_uk_cpp(List uk) {
  UK_E.clear(); UK_C.clear();
  for (int k = 0; k < uk.size(); ++k) {
    List el = uk[k];
    NumericVector e = el["e"], c = el["c"];
    UK_E.push_back(std::vector<double>(e.begin(), e.end()));
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
  double s = 1.0, vk = 1.0;
  for (int k = 1; k <= K; ++k) {
    vk *= sign / v;
    double uk = 0.0;
    const std::vector<double>& E = UK_E[k - 1];
    const std::vector<double>& C = UK_C[k - 1];
    for (size_t j = 0; j < E.size(); ++j) uk += C[j] * std::pow(t, E[j]);
    s += vk * uk;
  }
  return std::log(std::fabs(s));
}

static double lbi_series(double x, double v) {
  const int nmax = 90;
  double l2 = 2 * std::log(x) - std::log(4.0);
  double best = -1e308, la[nmax + 1];
  for (int k = 0; k <= nmax; ++k) {
    la[k] = k * l2 - std::lgamma(k + 1.0) - std::lgamma(k + v + 1.0);
    if (la[k] > best) best = la[k];
  }
  double s = 0.0;
  for (int k = 0; k <= nmax; ++k) s += std::exp(la[k] - best);
  return v * std::log(x / 2) + best + std::log(s);
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

// [[Rcpp::export]]
NumericVector log_bessel_i_cpp(NumericVector x, NumericVector nu) {
  R_xlen_t n = std::max(x.size(), nu.size());
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = x[i % x.size()], vi = nu[i % nu.size()];
    out[i] = (ISNAN(xi) || ISNAN(vi) || xi < 0 || vi < 0)
      ? NA_REAL : lbi_one(xi, vi);
  }
  return out;
}

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
