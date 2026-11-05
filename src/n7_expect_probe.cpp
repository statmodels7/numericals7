// Probes of inst/include/numericals7/n7_expect.h for its twin tests: each
// runs the compiled routine with an R function as the integrand, term or
// target, so that the tests compare it with quad_vec(), series_vec(),
// uniroot() and mean() on the same function.
#include <Rcpp.h>
#include "../inst/include/numericals7/n7_expect.h"

// [[Rcpp::export]]
double n7_mean_probe(Rcpp::NumericVector x) {
    return n7::r_mean(x.begin(), x.size());
}

// [[Rcpp::export]]
Rcpp::NumericVector n7_fmod_probe(Rcpp::NumericVector x, double y) {
    Rcpp::NumericVector out(x.size());
    for (R_xlen_t i = 0; i < x.size(); ++i) out[i] = n7::r_fmod(x[i], y);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector n7_series_steps_probe() {
    const std::vector<double>& s = n7::series_steps64();
    return Rcpp::NumericVector(s.begin(), s.end());
}

// one row per (lower, upper); f(x, i) is called with the fresh nodes of
// that row as a matrix of 15 columns, one row per panel, and i the row index
// [[Rcpp::export]]
Rcpp::NumericVector n7_quad_rows_probe(Rcpp::Function f, Rcpp::NumericVector lower,
                                       Rcpp::NumericVector upper, double atol,
                                       double rtol, int max_panels) {
    const int n = lower.size();
    Rcpp::NumericVector out(n);
    for (int r = 0; r < n; ++r) {
        auto fx = [&](const double* x, int nx, double* val) {
            const int np = nx / 15;
            Rcpp::NumericMatrix X(np, 15);
            for (int p = 0; p < np; ++p)
                for (int j = 0; j < 15; ++j) X(p, j) = x[p * 15 + j];
            Rcpp::NumericVector F = f(X, Rcpp::IntegerVector(np, r + 1));
            for (int p = 0; p < np; ++p)
                for (int j = 0; j < 15; ++j) val[p * 15 + j] = F[p + np * j];
        };
        double v;
        int st = n7::quad_row(fx, lower[r], upper[r], atol, rtol, max_panels, &v);
        out[r] = st == 0 ? v : NA_REAL;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector n7_series_rows_probe(Rcpp::Function term, int n, int from,
                                         double atol, double rtol, int max_terms) {
    Rcpp::NumericVector out(n);
    for (int r = 0; r < n; ++r) {
        auto tm = [&](const int* k, int nk, double* val) {
            Rcpp::IntegerVector K(k, k + nk);
            Rcpp::NumericVector V = term(K, Rcpp::IntegerVector(nk, r + 1));
            for (int j = 0; j < nk; ++j) val[j] = V[j];
        };
        double v;
        int st = n7::series_row(tm, from, atol, rtol, max_terms, &v);
        out[r] = st == 0 ? v : NA_REAL;
    }
    return out;
}

// [[Rcpp::export]]
double n7_uniroot_probe(Rcpp::Function f, double lower, double upper,
                        double tol, int maxiter) {
    auto g = [&](double x) { return Rcpp::as<double>(f(x)); };
    const double fl = g(lower), fu = g(upper);
    double root;
    int st = n7::uniroot(g, lower, upper, fl, fu, tol, maxiter, &root);
    if (st != 0) Rcpp::stop("uniroot code %d", st);
    return root;
}
