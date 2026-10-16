#include <Rcpp.h>
#include <algorithm>
#include <map>
#include <vector>
#include "n7_par.h"

// The tensor-product stencil of fd_tensor(). Two halves: the PLAN, which
// depends only on the integer structure of the request (which components,
// which multi-indices each analytic quantity supplies, the accuracy) and is
// cached by the R caller; and the COMBINATION, the weighted sum over the node
// values, which runs on vectors of length n. The node values themselves are
// computed in R between the two, because they are calls to R functions.

namespace {

// Central weights for the derivative of order m on the nodes -R..R, by
// Fornberg's recursion in extended precision. The rule is then made exactly
// symmetric (w(-k) = (-1)^m w(k)), so that the centre of an odd order is an
// exact zero and is dropped, as the Vandermonde solve of fd_weights() does.
std::vector<long double> central_weights(int m, int R) {
  const int N = 2 * R + 1;
  std::vector<long double> x(N);
  for (int i = 0; i < N; ++i) x[i] = (long double)(i - R);
  std::vector<std::vector<long double>> c(N, std::vector<long double>(m + 1, 0.0L));
  long double c1 = 1.0L, c4 = x[0];
  c[0][0] = 1.0L;
  for (int i = 1; i < N; ++i) {
    const int mn = std::min(i, m);
    long double c2 = 1.0L;
    const long double c5 = c4;
    c4 = x[i];
    for (int j = 0; j < i; ++j) {
      const long double c3 = x[i] - x[j];
      c2 *= c3;
      if (j == i - 1) {
        for (int k = mn; k >= 1; --k)
          c[i][k] = c1 * (k * c[i - 1][k - 1] - c5 * c[i - 1][k]) / c2;
        c[i][0] = -c1 * c5 * c[i - 1][0] / c2;
      }
      for (int k = mn; k >= 1; --k)
        c[j][k] = (c4 * c[j][k] - k * c[j][k - 1]) / c3;
      c[j][0] = c4 * c[j][0] / c3;
    }
    c1 = c2;
  }
  std::vector<long double> w(N);
  const long double sgn = (m % 2 == 0) ? 1.0L : -1.0L;
  for (int k = 1; k <= R; ++k) {
    const long double a = 0.5L * (c[R + k][m] + sgn * c[R - k][m]);
    w[R + k] = a;
    w[R - k] = sgn * a;
  }
  w[R] = (m % 2 == 0) ? c[R][m] : 0.0L;
  return w;
}

// The reach of fd_offsets(): the central rule of order m at accuracy a uses
// the nodes -R..R.
int stencil_reach(int m, int a) {
  return std::max(1, (m + a + 1) / 2 - 1);
}

struct Stencil {
  std::vector<int> off;
  std::vector<double> w;
};

const Stencil& stencil_of(std::map<int, Stencil>& cache, int m, int a) {
  auto it = cache.find(m);
  if (it != cache.end()) return it->second;
  const int R = stencil_reach(m, a);
  std::vector<long double> w = central_weights(m, R);
  Stencil s;
  for (int k = 0; k < (int)w.size(); ++k) {
    if (w[k] != 0.0L) {
      s.off.push_back(k - R);
      s.w.push_back((double)w[k]);
    }
  }
  return cache.emplace(m, s).first->second;
}

// Is candidate (rem_a, beta_a) to be preferred over (rem_b, beta_b) within one
// analytic quantity? The flatter residual first, compared as the residual
// sorted in decreasing order: this is what taking the most repeated
// coordinates into the analytic part first produces. Then the beta that is
// lexicographically larger, which is the tie rule of that greedy choice.
bool flatter(const std::vector<int>& rem_a, const std::vector<int>& beta_a,
             const std::vector<int>& rem_b, const std::vector<int>& beta_b) {
  std::vector<int> sa(rem_a), sb(rem_b);
  std::sort(sa.begin(), sa.end(), std::greater<int>());
  std::sort(sb.begin(), sb.end(), std::greater<int>());
  if (sa != sb) return sa < sb;
  return beta_a > beta_b;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List fd_tensor_plan_cpp(Rcpp::IntegerMatrix comp, Rcpp::List known,
                              int accuracy) {
  const int m = comp.nrow();
  const int p = comp.ncol();
  const int K = known.size();
  std::vector<Rcpp::IntegerMatrix> kn;
  for (int s = 0; s < K; ++s) kn.push_back(Rcpp::as<Rcpp::IntegerMatrix>(known[s]));

  Rcpp::IntegerVector comp_src(m), comp_col(m), comp_r(m);
  Rcpp::IntegerMatrix comp_rem(m, p);

  std::map<int, Stencil> st;
  std::map<std::vector<int>, int> node_id;
  std::vector<std::vector<int>> node_key;
  std::vector<int> t_comp, t_node;
  std::vector<double> t_w;
  std::map<int, std::vector<int>> reach;  // residual order -> reach per direction

  for (int c = 0; c < m; ++c) {
    std::vector<int> alpha(p);
    for (int j = 0; j < p; ++j) alpha[j] = comp(c, j);

    int bs = -1, bq = -1, bres = 0;
    std::vector<int> brem, bbeta;
    for (int s = 0; s < K; ++s) {
      const Rcpp::IntegerMatrix& B = kn[s];
      for (int q = 0; q < B.nrow(); ++q) {
        std::vector<int> beta(p), rem(p);
        bool ok = true;
        int res = 0;
        for (int j = 0; j < p; ++j) {
          beta[j] = B(q, j);
          rem[j] = alpha[j] - beta[j];
          if (rem[j] < 0) { ok = false; break; }
          res += rem[j];
        }
        if (!ok) continue;
        bool take;
        if (bs < 0) take = true;
        else if (res != bres) take = res < bres;
        else if (s != bs) {
          // more analytic derivatives first, then the order of declaration
          if (B.nrow() != kn[bs].nrow()) take = B.nrow() > kn[bs].nrow();
          else take = s < bs;
        } else take = flatter(rem, beta, brem, bbeta);
        if (take) { bs = s; bq = q; bres = res; brem = rem; bbeta = beta; }
      }
    }
    if (bs < 0) Rcpp::stop("No declared quantity lies below component %d.", c + 1);
    comp_src[c] = bs;
    comp_col[c] = bq;
    comp_r[c] = bres;
    for (int j = 0; j < p; ++j) comp_rem(c, j) = brem[j];

    std::vector<int>& rch = reach[bres];
    if (rch.empty()) rch.assign(p, 0);
    std::vector<int> dirs;
    std::vector<const Stencil*> sts;
    for (int j = 0; j < p; ++j) {
      if (brem[j] > 0) {
        dirs.push_back(j);
        sts.push_back(&stencil_of(st, brem[j], accuracy));
        rch[j] = std::max(rch[j], stencil_reach(brem[j], accuracy));
      }
    }
    // the tensor product, the first direction varying fastest
    std::vector<int> pos(dirs.size(), 0);
    while (true) {
      std::vector<int> key(p + 2, 0);
      key[0] = bs;
      double w = 1.0;
      bool zero = true;
      for (std::size_t d = 0; d < dirs.size(); ++d) {
        const int o = sts[d]->off[pos[d]];
        key[2 + dirs[d]] = o;
        w *= sts[d]->w[pos[d]];
        if (o != 0) zero = false;
      }
      // the centre is the same point whatever the step, so it is shared
      // across residual orders; any other node belongs to the steps of one
      key[1] = zero ? -1 : bres;
      auto it = node_id.find(key);
      int id;
      if (it == node_id.end()) {
        id = (int)node_key.size();
        node_id.emplace(key, id);
        node_key.push_back(key);
      } else id = it->second;
      t_comp.push_back(c);
      t_node.push_back(id);
      t_w.push_back(w);
      std::size_t d = 0;
      while (d < dirs.size()) {
        if (++pos[d] < (int)sts[d]->off.size()) break;
        pos[d] = 0;
        ++d;
      }
      if (d == dirs.size()) break;
    }
  }

  const int N = (int)node_key.size();
  Rcpp::IntegerVector node_src(N), node_r(N);
  Rcpp::IntegerMatrix node_off(N, p);
  for (int i = 0; i < N; ++i) {
    node_src[i] = node_key[i][0];
    node_r[i] = node_key[i][1];
    for (int j = 0; j < p; ++j) node_off(i, j) = node_key[i][2 + j];
  }
  const int nr = (int)reach.size();
  Rcpp::IntegerVector r_val(nr);
  Rcpp::IntegerMatrix r_reach(nr, p);
  int k = 0;
  for (auto& e : reach) {
    r_val[k] = e.first;
    for (int j = 0; j < p; ++j) r_reach(k, j) = e.second[j];
    ++k;
  }
  // terms were appended component by component, so they are already grouped
  Rcpp::IntegerVector t_start(m + 1);
  for (std::size_t t = 0; t < t_comp.size(); ++t) t_start[t_comp[t] + 1]++;
  for (int c = 0; c < m; ++c) t_start[c + 1] += t_start[c];

  return Rcpp::List::create(
    Rcpp::Named("comp_src") = comp_src, Rcpp::Named("comp_col") = comp_col,
    Rcpp::Named("comp_r") = comp_r, Rcpp::Named("comp_rem") = comp_rem,
    Rcpp::Named("node_src") = node_src, Rcpp::Named("node_r") = node_r,
    Rcpp::Named("node_off") = node_off,
    Rcpp::Named("term_start") = t_start,
    Rcpp::Named("term_node") = Rcpp::wrap(t_node),
    Rcpp::Named("term_w") = Rcpp::wrap(t_w),
    Rcpp::Named("r_val") = r_val, Rcpp::Named("r_reach") = r_reach);
}

// The weighted sums. vals holds one vector per (node, column) the terms read,
// of length 1 or n; term_val maps each term to its entry. Element i of every
// component is written by one thread, and its sum runs over the terms in the
// plan's order, so the result does not depend on the thread count.
// [[Rcpp::export]]
Rcpp::List fd_tensor_combine_cpp(Rcpp::List vals, Rcpp::IntegerVector term_start,
                                 Rcpp::IntegerVector term_val,
                                 Rcpp::NumericVector term_w, Rcpp::List dens,
                                 int n, int threads) {
  const int m = term_start.size() - 1;
  const int nv = vals.size();
  std::vector<const double*> vp(nv);
  std::vector<int> vlen(nv);
  for (int v = 0; v < nv; ++v) {
    Rcpp::NumericVector x = vals[v];
    vp[v] = x.begin();
    vlen[v] = x.size();
  }
  std::vector<const double*> dp(m);
  std::vector<int> dlen(m);
  for (int c = 0; c < m; ++c) {
    Rcpp::NumericVector x = dens[c];
    dp[c] = x.begin();
    dlen[c] = x.size();
  }
  Rcpp::List out(m);
  std::vector<double*> op(m);
  for (int c = 0; c < m; ++c) {
    Rcpp::NumericVector o(n);
    op[c] = o.begin();
    out[c] = o;
  }
  const int* ts = term_start.begin();
  const int* tv = term_val.begin();
  const double* tw = term_w.begin();

  // The work is cut into blocks of elements, and within a block each term is
  // one pass over contiguous memory, so the inner loop streams two arrays
  // rather than gathering from every node vector at each element. The terms
  // of an element are still summed in the plan's order.
  constexpr std::size_t kBlock = 512;
  // a few multiply-adds per element and term: below this many blocks the
  // pool costs more than it saves
  constexpr std::size_t kMinBlocks = 32;
  const std::size_t nn = (std::size_t)n;
  const std::size_t nb = (nn + kBlock - 1) / kBlock;
  n7::par_for(nb, threads, kMinBlocks, [&](std::size_t b) {
    const std::size_t lo = b * kBlock;
    const std::size_t hi = std::min(nn, lo + kBlock);
    for (int c = 0; c < m; ++c) {
      double* o = op[c];
      for (std::size_t i = lo; i < hi; ++i) o[i] = 0.0;
      for (int t = ts[c]; t < ts[c + 1]; ++t) {
        const int v = tv[t];
        const double w = tw[t];
        const double* x = vp[v];
        if (vlen[v] == 1) {
          const double a = w * x[0];
          for (std::size_t i = lo; i < hi; ++i) o[i] += a;
        } else {
          for (std::size_t i = lo; i < hi; ++i) o[i] += w * x[i];
        }
      }
      const double* d = dp[c];
      if (dlen[c] == 1) {
        for (std::size_t i = lo; i < hi; ++i) o[i] /= d[0];
      } else {
        for (std::size_t i = lo; i < hi; ++i) o[i] /= d[i];
      }
    }
  });
  return out;
}
