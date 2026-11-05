// The numerical routines that the R route of an expected information runs,
// transcribed for compiled callers so that they return the same doubles:
// R's mean() of a double vector, R's %% on doubles, R's Brent zero finder
// with the wrapper uniroot() puts around it, one row of quad_vec() and one
// row of series_vec(). Header-only; read with LinkingTo: numericals7.
//
// Sources: R 4.6.0 (tag R-4-6-0 of wch/r-source): summary.c real_mean(),
// arithmetic.c myfmod(), stats/src/zeroin.c R_zeroin2(), stats/src/optimize.c
// fcn2(), the R function uniroot(); numericals7 R/quadrature.R quad_vec() and
// R/series.R series_vec(). Every quantity R sums in LDOUBLE is summed here in
// long double, every one R sums in double in double, in R's order.

#ifndef NUMERICALS7_N7_EXPECT_H
#define NUMERICALS7_N7_EXPECT_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include <algorithm>
#include <climits>

namespace n7 {

// ------------------------------------------------------------------ mean
// summary.c real_mean() of R 4.6.0
inline double r_mean(const double* x, R_xlen_t n) {
    long double s = 0.0;
    for (R_xlen_t k = 0; k < n; k++) s += x[k];
    const bool finite_s = R_FINITE((double) s) != 0;
    if (finite_s) {
        s /= n;
    } else {
        s = 0.;
        for (R_xlen_t k = 0; k < n; k++) s += x[k] / n;
    }
    if (finite_s && R_FINITE((double) s)) {
        long double t = 0.0;
        for (R_xlen_t k = 0; k < n; k++) t += (x[k] - s);
        s += t / n;
    } else if (R_FINITE((double) s)) {
        long double t = 0.0;
        for (R_xlen_t k = 0; k < n; k++) t += (x[k] - s) / n;
        s += t;
    }
    return (double) s;
}

// ------------------------------------------------------------- x %% y
// arithmetic.c myfmod() of R 4.6.0 (the warning on a complete loss of
// accuracy is not reproduced; it does not change the value)
inline double r_fmod(double x1, double x2) {
    const double c_eps = DBL_EPSILON;
    if (x2 == 0.0) return R_NaN;
    if (std::fabs(x2) * c_eps > 1 && R_FINITE(x1) && std::fabs(x1) <= std::fabs(x2)) {
        return (std::fabs(x1) == std::fabs(x2)) ? 0 :
            ((x1 < 0 && x2 > 0) || (x2 < 0 && x1 > 0)) ? x1 + x2 : x1;
    }
    double q = x1 / x2;
    long double tmp = (long double) x1 - std::floor(q) * (long double) x2;
    return (double) (tmp - floorl(tmp / x2) * x2);
}

// x^y as R's vectorized ^ computes it (arithmetic.h R_POW)
inline double r_pow(double x, double y) {
    return y == 2.0 ? x * x : R_pow(x, y);
}

// R's pmax(a, b) of two doubles: NaN if either is NaN
inline double r_pmax(double a, double b) {
    if (ISNAN(a)) return a;
    if (ISNAN(b)) return b;
    return a >= b ? a : b;
}

// ---------------------------------------------------------------- Brent
// stats/src/zeroin.c R_zeroin2() of R 4.6.0. f is called through fcn(),
// which applies fcn2()'s replacement of a non-finite value.
template <class F>
inline double zeroin2_body(double ax, double bx, double fa, double fb, F& fcn,
                           double* Tol, int* Maxit) {
    double a, b, c, fc;
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = *Tol;

    if (fa == 0.0) {
        *Tol = 0.0;
        *Maxit = 0;
        return a;
    }
    if (fb == 0.0) {
        *Tol = 0.0;
        *Maxit = 0;
        return b;
    }

    while (maxit--) {
        double prev_step = b - a;
        double tol_act;
        double p;
        double q;
        double new_step;

        if (std::fabs(fc) < std::fabs(fb)) {
            a = b;  b = c;  c = a;
            fa = fb;  fb = fc;  fc = fa;
        }
        tol_act = 2 * DBL_EPSILON * std::fabs(b) + tol / 2;
        new_step = (c - b) / 2;

        if (std::fabs(new_step) <= tol_act || fb == (double) 0) {
            *Maxit -= maxit;
            *Tol = std::fabs(c - b);
            return b;
        }

        if (std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
            double t1, cb, t2;
            cb = c - b;
            if (a == c) {
                t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                q = fa / fc;  t1 = fb / fc;  t2 = fb / fa;
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if (p > (double) 0)
                q = -q;
            else
                p = -p;

            if (p < (0.75 * cb * q - std::fabs(tol_act * q) / 2)
                && p < std::fabs(prev_step * q / 2))
                new_step = p / q;
        }

        if (std::fabs(new_step) < tol_act) {
            if (new_step > (double) 0)
                new_step = tol_act;
            else
                new_step = -tol_act;
        }
        a = b;  fa = fb;
        b += new_step;  fb = fcn(b);
        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a;  fc = fa;
        }
    }
    *Tol = std::fabs(c - b);
    *Maxit = -1;
    return b;
}

// uniroot(f, c(lower, upper), f.lower = fl, f.upper = fu, tol = tol,
// maxiter = maxiter)$root with extendInt = "no" and check.conv = FALSE.
// Returns 0 and the root in *root, or a nonzero code where uniroot()
// signals an error: 1 lower >= upper, 2 an NA end value, 3 no sign change.
// A run that reaches maxiter returns its last point with code 0, as uniroot()
// does (with a warning).
template <class F>
inline int uniroot(F&& f, double lower, double upper, double fl, double fu,
                   double tol, int maxiter, double* root) {
    if (!(lower < upper)) return 1;
    if (ISNAN(fl) || ISNAN(fu)) return 2;
    auto trunc = [](double x) {
        return std::max(std::min(x, DBL_MAX), -DBL_MAX);
    };
    const double flo = trunc(fl), fup = trunc(fu);
    // sign(f.lower) * sign(f.upper) <= 0
    auto sgn = [](double x) { return (double) ((x > 0) - (x < 0)); };
    if (!(sgn(fl) * sgn(fu) <= 0)) return 3;
    auto fcn = [&](double x) -> double {
        double v = f(x);
        if (!R_FINITE(v)) {
            if (v == R_NegInf) return -DBL_MAX;
            return DBL_MAX;
        }
        return v;
    };
    double t = tol;
    int it = maxiter;
    *root = zeroin2_body(lower, upper, flo, fup, fcn, &t, &it);
    return 0;
}

// ------------------------------------------------------------ quad_vec
// The 15-point Gauss-Kronrod pair as R holds it (gauss_kronrod15()), by bit
// pattern, nodes ascending.
struct GK15 {
    double node[15], kw[15], gw[15];
    GK15() {
        static const double xh[8] = {
            0x1.fba009d4d09b1p-1, 0x1.e5f178e7c6229p-1, 0x1.bacf827b9bb3ep-1,
            0x1.7ba9f9be3a1d6p-1, 0x1.2c13a049dfa24p-1, 0x1.9f95df119fd62p-2,
            0x1.a98b2892e0c77p-3, 0x0p+0};
        static const double wk[8] = {
            0x1.77c5b67d5747p-6, 0x1.026cdaa7b61c4p-4, 0x1.ad384a34814c6p-4,
            0x1.200ed0f46e8c1p-3, 0x1.5a1f266e47d5cp-3, 0x1.85d6861c80eb1p-3,
            0x1.a2adbcbec9cd8p-3, 0x1.ad04f9087090fp-3};
        static const double wg[8] = {
            0, 0x1.092f69f826d57p-3, 0, 0x1.1e6b1713d8644p-2, 0,
            0x1.86fe74ee32b3dp-2, 0, 0x1.abfd7e03c2fa6p-2};
        for (int j = 0; j < 7; ++j) {
            node[j] = -xh[j];        node[14 - j] = xh[j];
            kw[j] = wk[j];           kw[14 - j] = wk[j];
            gw[j] = wg[j];           gw[14 - j] = wg[j];
        }
        node[7] = xh[7]; kw[7] = wk[7]; gw[7] = wg[7];
    }
};

inline const GK15& gk15() {
    static const GK15 r;
    return r;
}

// One row of quad_vec(f, lower, upper, atol, rtol, max_depth, max_panels)
// with the GK15 rule. f(const double* x, int nx, double* out) evaluates the
// integrand elementwise at nx points (15 per fresh panel, panel by panel).
// Returns 0 and the integral in *value where the row converges, 1 where
// quad_vec() returns NA for it.
template <class F>
inline int quad_row(F&& f, double lower, double upper, double atol, double rtol,
                    int max_panels, double* value, int max_depth = 48) {
    const GK15& R = gk15();
    const int m = 15;
    int type;
    if (R_FINITE(lower) && R_FINITE(upper)) type = 0;
    else if (R_FINITE(lower)) type = 1;
    else if (R_FINITE(upper)) type = 2;
    else type = 3;
    const double t0 = type == 0 ? lower : (type == 3 ? -1 : 0);
    const double t1 = type == 0 ? upper : 1;

    struct Panel { double a, b, K, err; int depth; bool fresh; };
    std::vector<Panel> pl;
    pl.push_back({t0, t1, NA_REAL, R_PosInf, 0, true});

    std::vector<double> X, J, Fv;
    int iter = 0;
    const int iter_max = 40 + 4 * max_depth;
    while (iter < iter_max) {
        ++iter;
        // evaluate the fresh panels, in list order
        int nf = 0;
        for (const Panel& P : pl) if (P.fresh) ++nf;
        if (nf) {
            X.assign((std::size_t) nf * m, 0.0);
            J.assign((std::size_t) nf * m, 1.0);
            int r = 0;
            for (const Panel& P : pl) {
                if (!P.fresh) continue;
                const double half = (P.b - P.a) / 2;
                for (int j = 0; j < m; ++j) {
                    const double T = (P.a + half) * 1 + half * R.node[j];
                    double x = T, jac = 1;
                    if (type == 1) {
                        x = lower + T / (1 - T);
                        jac = 1 / ((1 - T) * (1 - T));
                    } else if (type == 2) {
                        x = upper - T / (1 - T);
                        jac = 1 / ((1 - T) * (1 - T));
                    } else if (type == 3) {
                        x = T / (1 - T * T);
                        jac = (1 + T * T) / ((1 - T * T) * (1 - T * T));
                    }
                    X[(std::size_t) r * m + j] = x;
                    J[(std::size_t) r * m + j] = jac;
                }
                ++r;
            }
            Fv.assign((std::size_t) nf * m, 0.0);
            f(X.data(), nf * m, Fv.data());
            r = 0;
            for (Panel& P : pl) {
                if (!P.fresh) continue;
                const double half = (P.b - P.a) / 2;
                double sk = 0.0, sg = 0.0;
                for (int j = 0; j < m; ++j) {
                    const double fm = Fv[(std::size_t) r * m + j] * J[(std::size_t) r * m + j];
                    sk = sk + fm * R.kw[j];
                    sg = sg + fm * R.gw[j];
                }
                P.K = half * sk;
                const double G = half * sg;
                P.err = std::fabs(P.K - G);
                if (!R_FINITE(P.err)) P.err = R_PosInf;
                P.fresh = false;
                ++r;
            }
        }
        double tot = 0.0, te = 0.0;
        for (const Panel& P : pl) { tot = tot + P.K; te = te + P.err; }
        const double budget = r_pmax(atol, rtol * std::fabs(tot));
        if (R_FINITE(tot) && te <= budget) {
            *value = tot;
            return 0;
        }
        double irr = 0.0, mx = R_NegInf;
        int nsp = 0;
        const int cnt = (int) pl.size();
        for (const Panel& P : pl) {
            const bool sp = P.depth < max_depth;
            irr = irr + (sp ? 0 : P.err);
            if (sp) { ++nsp; if (P.err > mx) mx = P.err; }
        }
        if (nsp == 0 || irr > budget || cnt >= max_panels) return 1;
        std::vector<Panel> keep, left, right;
        for (const Panel& P : pl) {
            const bool sp = P.depth < max_depth;
            if (sp && P.err >= 0.5 * mx) {
                const double mid = (P.a + P.b) / 2;
                left.push_back({P.a, mid, NA_REAL, R_PosInf, P.depth + 1, true});
                right.push_back({mid, P.b, NA_REAL, R_PosInf, P.depth + 1, true});
            } else {
                keep.push_back(P);
            }
        }
        pl.clear();
        pl.insert(pl.end(), keep.begin(), keep.end());
        pl.insert(pl.end(), left.begin(), left.end());
        pl.insert(pl.end(), right.begin(), right.end());
    }
    return 1;
}

// ----------------------------------------------------------- series_vec
// The probe distances of series_vec() for a block of 64:
// c(0, unique(floor(64 * (1 + 1/128)^(0:nj)))).
inline const std::vector<double>& series_steps64() {
    static std::vector<double> s;
    if (s.empty()) {
        const int block = 64;
        const double nj = std::ceil(std::log(R_pow(2.0, 31.0) / block) / std::log1p(1.0 / 128));
        s.push_back(0);
        double last = R_NaN;
        // the powers are nondecreasing, so unique() drops repeats in a row
        for (int k = 0; k <= (int) nj; ++k) {
            const double v = std::floor(block * r_pow(1 + 1.0 / 128, (double) k));
            if (!(v == last)) s.push_back(v);
            last = v;
        }
    }
    return s;
}

// One row of series_vec(term, n, from, atol, rtol, max_terms, block = 64).
// term(const int* k, int nk, double* out) evaluates the term elementwise.
// Returns 0 and the sum in *value, 1 where series_vec() returns NA for the
// row, 2 where an index would pass INT_MAX (R's integer overflow, handed
// back to the caller).
template <class T>
inline int series_row(T&& term, int from, double atol, double rtol, int max_terms,
                      double* value) {
    const int block = 64, qb = 16;
    const std::vector<double>& steps = series_steps64();
    double acc = 0.0;
    bool seen = false;
    long long kst = from;
    int used = 0;
    std::vector<int> K(block);
    std::vector<double> M(block);
    for (;;) {
        if (kst + block - 1 > INT_MAX) return 2;
        for (int j = 0; j < block; ++j) K[j] = (int) (kst + j);
        term(K.data(), block, M.data());
        long double c = 0.0, s3 = 0.0, s4 = 0.0, sa = 0.0;
        for (int j = 0; j < block; ++j) c += M[j];
        const double contrib = (double) c;
        acc = acc + contrib;
        const double budget = r_pmax(atol, rtol * std::fabs(acc));
        for (int j = block - 2 * qb; j < block - qb; ++j) s3 += std::fabs(M[j]);
        for (int j = block - qb; j < block; ++j) s4 += std::fabs(M[j]);
        for (int j = 0; j < block; ++j) sa += std::fabs(M[j]);
        const double a3 = (double) s3, a4 = (double) s4;
        seen = seen || ((double) sa > 0);
        const double q = a4 / a3;
        const double tail = (a4 == 0) ? 0 : 2 * a4 * q / (1 - q);
        bool done = (a4 < a3 || a4 == 0) && tail <= budget &&
                    std::fabs(contrib) <= budget && seen;
        // R ends an NA expression FALSE; a NaN comparison is false in C, and
        // the expression has no negation, so the two agree
        used += block;
        kst += block;
        if (!done && !seen) {
            // the probe on the geometric grid
            int first = -1;
            std::vector<int> pk;
            std::vector<int> idx;
            for (std::size_t s = 0; s < steps.size(); ++s) {
                const double P = (double) kst + steps[s];
                if (P > INT_MAX) continue;
                pk.push_back((int) P);
                idx.push_back((int) s);
            }
            std::vector<double> pv(pk.size());
            if (!pk.empty()) term(pk.data(), (int) pk.size(), pv.data());
            for (std::size_t s = 0; s < pk.size(); ++s) {
                if (ISNAN(pv[s]) || pv[s] != 0) { first = idx[s]; break; }
            }
            if (first < 0) {
                done = true;
            } else if (first > 0) {
                kst = (long long) ((double) (kst) + steps[first - 1]);
            }
        }
        if (done) { *value = acc; return 0; }
        if (used >= max_terms) return 1;
    }
}

}  // namespace n7

#endif
