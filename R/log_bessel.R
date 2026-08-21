#' @include quadrature.R bessel_uk_table.R
NULL

# The logarithm of the modified Bessel functions, after Plesner, Sorensen &
# Hauberg (2024): every intermediate quantity is carried on the log scale,
# so the results are finite wherever the functions' own logs are, including
# the regions where the unscaled functions overflow and the exponentially
# scaled ones underflow to an exact zero.

# ---- the switching table (Table 1 / Algorithm 1 of the paper) --------------
# 1 mu3, 2 mu20, 3 U4, 4 U6, 5 U9, 6 U13, 0 fallback. Two guards are ours,
# beyond the paper's table, both measured on the Wronskian identity:
# the mu expansions truncate, and at the paper's own large-order boundary
# they leave ~1e-9 of error where the U_K route sits at ~6e-12. mu20's
# large-order clause therefore demands v^2 < 2x (mu/8x < 1/4), and mu3's,
# whose error grows like (v^2/2x)^4, demands v^2 < 0.002 x; what they lose
# falls to clauses that are supersets there.
.lb_branch <- function(x, v) {
  b <- integer(length(x))
  lx <- log(x)
  lv <- log(v)
  m3 <- (x > 1400 & v < 3.05) |
    ((0.6229 * lx - 3.2318) > lv & v > 3.1 & v * v < 0.002 * x)
  m20 <- (x > 30 & v < 15.3919) |
    ((0.5113 * lx + 0.7939) > lv & x > 59.6925 & v * v < 2 * x)
  u4 <- (x > 274.2377 & v > 0.3) | (v > 163.6993)
  u6 <- (x > 84.4153 & v > 0.46) | (v > 56.9971)
  u9 <- (x > 35.9074 & v > 0.6) | (v > 20.1534)
  u13 <- (x > 19.6931 & v > 0.7) | (v > 12.6964)
  b[u13] <- 6L; b[u9] <- 5L; b[u6] <- 4L; b[u4] <- 3L
  b[m20] <- 2L; b[m3] <- 1L
  b
}

# ---- the series for log I (moderate x and v) --------------------------------
# log a_k = k log(x^2/4) - lgamma(k+1) - lgamma(k+v+1), summed by the
# log-sum-exp anchored at the largest term. The number of relevant terms
# grows like 9.2 sqrt(x) (Plesner et al.), so the truncation adapts to the
# chunk instead of always paying the worst case.
.lbi_series <- function(x, v) {
  nmax <- min(120L, as.integer(ceiling(9.2 * sqrt(max(x))) + 20L))
  k <- 0:nmax
  la <- outer(2 * log(x) - log(4), k) -
    rep(lgamma(k + 1), each = length(x)) -
    lgamma(outer(v, k, "+") + 1)
  m <- apply(la, 1L, max)
  v * log(x / 2) + m + log(rowSums(exp(la - m)))
}

# ---- the mu_K expansion (large argument) ------------------------------------
# sign is -1 for I and +1 for K; the absolute value guards the truncated
# series against a numerically negative sum
.lb_mu <- function(x, v, K, sign) {
  mu <- 4 * v * v
  s <- rep(1, length(x))
  ck <- rep(1, length(x))
  for (k in seq_len(K)) {
    ck <- sign * ck * (mu - (2 * k - 1)^2) / (k * 8 * x)
    s <- s + ck
  }
  log(abs(s))
}

# ---- the U_K uniform expansion (large order) --------------------------------
# sign is +1 for I and -1 for K
.lb_u_sum <- function(t, v, K, sign) {
  s <- rep(1, length(t))
  vk <- rep(1, length(t))
  for (k in seq_len(K)) {
    vk <- vk * (sign / v)
    ek <- .bessel_uk[[k]]$e
    ck <- .bessel_uk[[k]]$c
    uk <- rep(0, length(t))
    for (j in seq_along(ek)) uk <- uk + ck[j] * t^ek[j]
    s <- s + vk * uk
  }
  log(abs(s))
}

.lb_u_parts <- function(x, v) {
  xp <- x / v
  r <- sqrt(1 + xp * xp)
  list(t = 1 / r, eta = r + log(xp / (1 + r)), lr = log1p(xp * xp))
}

.lbi_u <- function(x, v, K) {
  p <- .lb_u_parts(x, v)
  -0.5 * log(2 * pi * v) + v * p$eta - 0.25 * p$lr + .lb_u_sum(p$t, v, K, +1)
}

.lbk_u <- function(x, v, K) {
  p <- .lb_u_parts(x, v)
  0.5 * (log(pi) - log(2 * v)) - v * p$eta - 0.25 * p$lr +
    .lb_u_sum(p$t, v, K, -1)
}

# ---- log K at moderate inputs ------------------------------------------------
# R's exponentially scaled besselK is machine-precision exact in the
# fallback region wherever it does not overflow, and one call beats a
# 600-node Simpson rule; the Rothwell integral serves the corner where the
# scaled value overflows, and doubles as the independent reference in the
# tests.
.lbk_moderate <- function(x, v) {
  out <- rep(NA_real_, length(x))
  big <- v * log(2 / pmax(x, 1e-300)) + lgamma(pmax(v, 0.5)) > 690
  if (any(!big)) {
    out[!big] <- log(besselK(x[!big], v[!big], expon.scaled = TRUE)) - x[!big]
  }
  if (any(big)) out[big] <- .lbk_integral(x[big], v[big])
  out
}

# The integral representation of Rothwell (2006), eq. (26)-(27), with the
# integrand evaluated on the log scale and summed by the anchored
# log-sum-exp over the composite Simpson rule.
.lbk_integral <- function(x, v, nsimp = 600L) {
  n <- 8
  beta <- 2 * n / (2 * v + 1)
  u <- seq_len(nsimp - 1L) / nsimp
  lw <- log(rep(c(4, 2), length.out = nsimp - 1L))

  Ub <- outer(rep(1, length(x)), u)
  B <- matrix(beta, length(x), nsimp - 1L)
  lg <- log(B) - Ub^B + (v - 0.5) * log(2 * x + Ub^B) + (n - 1) * log(Ub)
  lh <- -1 / Ub - (2 * v + 1) * log(Ub) + (v - 0.5) * log(2 * x * Ub + 1)

  lg1 <- log(beta) - 1 + (v - 0.5) * log(2 * x + 1)
  lh1 <- -1 + (v - 0.5) * log(2 * x + 1)

  la <- cbind(lg + rep(lw, each = length(x)), lg1,
              lh + rep(lw, each = length(x)), lh1)
  m <- apply(la, 1L, max)
  lint <- m + log(rowSums(exp(la - m))) - log(3 * nsimp)

  0.5 * log(pi) - lgamma(v + 0.5) - v * log(2 * x) - x + lint
}

.lb_dispatch <- function(x, v, funs) {
  out <- rep(NA_real_, length(x))
  b <- .lb_branch(x, v)
  for (code in unique(b)) {
    i <- b == code
    out[i] <- funs[[code + 1L]](x[i], v[i])
  }
  out
}

#' Logarithm of the Modified Bessel Function of the First Kind
#'
#' @description
#' \eqn{\log I_\nu(x)} for \eqn{x \ge 0} and \eqn{\nu \ge 0}, computed with
#' every intermediate quantity on the log scale, so the result is finite and
#' accurate wherever \eqn{\log I_\nu(x)} itself is representable -- including
#' where the unscaled function overflows (from about \eqn{x = 700}) and where
#' the exponentially scaled one underflows or loses its precision (large
#' order with a small argument, and arguments beyond about \eqn{10^5}).
#'
#' @details
#' The algorithm follows Plesner, Sørensen and Hauberg (2024): the truncated
#' power series with its terms carried as logarithms and summed by the
#' log-sum-exp anchored at the largest term, the large-argument asymptotic
#' expansion, and the large-order uniform asymptotic expansion in the
#' \eqn{u_k(t)} polynomials, selected by the input-range table of that paper.
#' Two switching guards are tightened relative to the paper, measured on the
#' Wronskian identity \eqn{I_\nu K_{\nu+1} + I_{\nu+1} K_\nu = 1/x}: the
#' large-argument expansions are used only where their truncation error is
#' below \eqn{10^{-11}}, the uniform expansion taking over otherwise.
#'
#' @param x A numeric vector of non-negative arguments, recycled against
#'   \code{nu}.
#' @param nu A numeric vector of non-negative orders.
#' @param threads How many threads the kernel may use, a plain integer, as
#'   \code{\link{thread_count}} reads it off an \code{\link{n_threads}}
#'   policy. Every branch of the kernel is this package's own arithmetic, so
#'   element \eqn{i} is computed and written by one thread and the result is
#'   bit-identical at any count; below an internal threshold the sequential
#'   path is taken whatever the count says. \code{\link{log_bessel_k}} takes
#'   no such argument: its hybrid branch calls R's own scaled \code{besselK},
#'   which can raise a warning, and a warning from a worker thread ends the
#'   session.
#'
#' @return A numeric vector of \eqn{\log I_\nu(x)} values; \code{-Inf} at
#'   \eqn{x = 0} with \eqn{\nu > 0}, \code{0} at \eqn{x = 0} with
#'   \eqn{\nu = 0}, \code{NA} for arguments outside the domain.
#'
#' @references
#' Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
#' computation of the logarithm of modified Bessel functions on GPUs.
#' \emph{Proceedings of the 38th ACM International Conference on
#' Supercomputing (ICS '24)}. arXiv:2409.08729.
#'
#' Olver, F. W. J., et al. (2024). \emph{NIST Digital Library of
#' Mathematical Functions}, chapter 10. https://dlmf.nist.gov/.
#'
#' @seealso \code{\link{log_bessel_k}}, \code{\link{log_bessel_i_derivs}},
#'   \code{\link{bessel_i_ratio}}
#'
#' @examples
#' log_bessel_i(c(1, 100, 1e6), 2)      # far past the unscaled overflow
#' log_bessel_i(0.001, 500)             # far past the scaled underflow
#'
#' @export
log_bessel_i <- function(x, nu, threads = 1L) {
  log_bessel_i_cpp(as.numeric(x), as.numeric(nu), as.integer(threads))
}

#' The R Twin of the Compiled log I Kernel
#'
#' @description
#' The same branches and formulas as the compiled kernel behind
#' \code{\link{log_bessel_i}}, in vectorized R. The compiled route measured
#' 1.1x faster on a mixed workload of one million points; the twin is kept
#' as the independent reference the tests compare against.
#'
#' @param x A non-negative numeric vector.
#' @param nu A non-negative numeric vector.
#' @return A numeric vector.
#' @keywords internal
.log_bessel_i_r <- function(x, nu) {
  nn <- max(length(x), length(nu))
  x <- rep_len(x, nn)
  v <- rep_len(nu, nn)
  out <- rep(NA_real_, nn)
  bad <- is.na(x) | is.na(v) | x < 0 | v < 0
  zero <- !bad & x == 0
  out[zero] <- ifelse(v[zero] == 0, 0, -Inf)
  todo <- which(!bad & !zero)
  if (length(todo)) {
    out[todo] <- .lb_dispatch(x[todo], v[todo], list(
      .lbi_series,
      function(x, v) .lbi_mu_full(x, v, 3L),
      function(x, v) .lbi_mu_full(x, v, 20L),
      function(x, v) .lbi_u(x, v, 4L),
      function(x, v) .lbi_u(x, v, 6L),
      function(x, v) .lbi_u(x, v, 9L),
      function(x, v) .lbi_u(x, v, 13L)
    ))
  }
  out
}

.lbi_mu_full <- function(x, v, K) {
  x - 0.5 * log(2 * pi * x) + .lb_mu(x, v, K, -1)
}

.lbk_mu_full <- function(x, v, K) {
  0.5 * (log(pi) - log(2 * x)) - x + .lb_mu(x, v, K, +1)
}

#' Logarithm of the Modified Bessel Function of the Second Kind
#'
#' @description
#' \eqn{\log K_\nu(x)} for \eqn{x > 0} and any real order, computed with
#' every intermediate quantity on the log scale; \eqn{K} is even in its
#' order, so \eqn{\nu} enters as \eqn{|\nu|}.
#'
#' @details
#' The large-argument and large-order branches follow Plesner, Sørensen and
#' Hauberg (2024) exactly as in \code{\link{log_bessel_i}}. At moderate
#' inputs the exponentially scaled \code{\link[base]{besselK}} is
#' machine-precision exact wherever it does not overflow and one call beats
#' a quadrature, so it serves that region; the corner where the scaled value
#' itself overflows (a small argument with the order near the switching
#' boundary) goes through the integral representation of Rothwell (2006),
#' evaluated on the log scale over a composite Simpson rule.
#'
#' @param x A numeric vector of positive arguments, recycled against
#'   \code{nu}.
#' @param nu A numeric vector of orders, any sign.
#'
#' @return A numeric vector of \eqn{\log K_\nu(x)} values; \code{Inf} at
#'   \eqn{x = 0}, \code{NA} for arguments outside the domain.
#'
#' @references
#' Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
#' computation of the logarithm of modified Bessel functions on GPUs.
#' \emph{Proceedings of the 38th ACM International Conference on
#' Supercomputing (ICS '24)}. arXiv:2409.08729.
#'
#' Rothwell, E. J. (2006). Computation of the logarithm of Bessel functions
#' of complex argument and fractional order. \emph{Communications in
#' Numerical Methods in Engineering} 24(3), 237--249.
#'
#' @seealso \code{\link{log_bessel_i}}, \code{\link{log_bessel_k_derivs}}
#'
#' @examples
#' log_bessel_k(c(0.01, 1, 1000), 2.5)
#' log_bessel_k(1, 500)                 # the unscaled function overflows here
#'
#' @export
log_bessel_k <- function(x, nu) {
  log_bessel_k_cpp(as.numeric(x), as.numeric(nu))
}

#' The R Twin of the Compiled log K Kernel
#'
#' @description
#' The same branches and formulas as the compiled kernel behind
#' \code{\link{log_bessel_k}}, in vectorized R. The compiled route measured
#' 2.9x faster on the mixed workload; the twin is kept as the independent
#' reference the tests compare against.
#'
#' @param x A non-negative numeric vector.
#' @param nu A numeric vector; the function is even in the order.
#' @return A numeric vector.
#' @keywords internal
.log_bessel_k_r <- function(x, nu) {
  nn <- max(length(x), length(nu))
  x <- rep_len(x, nn)
  v <- abs(rep_len(nu, nn))
  out <- rep(NA_real_, nn)
  bad <- is.na(x) | is.na(v) | x < 0
  zero <- !bad & x == 0
  out[zero] <- Inf
  todo <- which(!bad & !zero)
  if (length(todo)) {
    out[todo] <- .lb_dispatch(x[todo], v[todo], list(
      .lbk_moderate,
      function(x, v) .lbk_mu_full(x, v, 3L),
      function(x, v) .lbk_mu_full(x, v, 20L),
      function(x, v) .lbk_u(x, v, 4L),
      function(x, v) .lbk_u(x, v, 6L),
      function(x, v) .lbk_u(x, v, 9L),
      function(x, v) .lbk_u(x, v, 13L)
    ))
  }
  out
}

# ---- derivatives in the argument ---------------------------------------------
# The first derivative is the ratio identity with the ratio formed as
# exp(difference of logs); the higher orders follow from the modified Bessel
# equation, which both I and K satisfy: with L = log B and r = L',
#   L'' = 1 + (v/x)^2 - r/x - r^2,
# and each further order differentiates that expression.
.lb_deriv_chain <- function(x, v, r1) {
  L2 <- 1 + (v / x)^2 - r1 / x - r1 * r1
  L3 <- -2 * v^2 / x^3 + r1 / x^2 - L2 / x - 2 * r1 * L2
  L4 <- 6 * v^2 / x^4 - 2 * r1 / x^3 + 2 * L2 / x^2 - L3 / x -
    2 * L2 * L2 - 2 * r1 * L3
  list(d1 = r1, d2 = L2, d3 = L3, d4 = L4)
}

#' Derivatives of the Logarithm of the Modified Bessel Function I
#'
#' @description
#' \eqn{\log I_\nu(x)} together with its first four derivatives with respect
#' to the \emph{argument}. The first derivative is the ratio identity
#' \eqn{(\log I_\nu)' = \nu/x + I_{\nu+1}/I_\nu} with the ratio formed as
#' the exponential of a difference of logarithms, which is finite wherever
#' the logs are; the higher orders follow from the modified Bessel equation
#' and cost no further Bessel evaluations. Derivatives with respect to the
#' order have no elementary form and are not provided.
#'
#' @param x A numeric vector of positive arguments, recycled against
#'   \code{nu}.
#' @param nu A numeric vector of non-negative orders.
#'
#' @return A list with the value \code{l} and the derivatives \code{d1} to
#'   \code{d4} in the argument.
#'
#' @references
#' Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
#' computation of the logarithm of modified Bessel functions on GPUs.
#' \emph{Proceedings of the 38th ACM International Conference on
#' Supercomputing (ICS '24)}. arXiv:2409.08729.
#'
#' @seealso \code{\link{log_bessel_i}}, \code{\link{log_bessel_k_derivs}}
#'
#' @examples
#' log_bessel_i_derivs(2, 0.5)$d1   # equals coth(2) - 1/(2*2) exactly
#'
#' @export
log_bessel_i_derivs <- function(x, nu) {
  nn <- max(length(x), length(nu))
  x <- rep_len(x, nn)
  v <- rep_len(nu, nn)
  l0 <- log_bessel_i(x, v)
  r1 <- v / x + exp(log_bessel_i(x, v + 1) - l0)
  c(list(l = l0), .lb_deriv_chain(x, v, r1))
}

#' Derivatives of the Logarithm of the Modified Bessel Function K
#'
#' @description
#' \eqn{\log K_\nu(x)} together with its first four derivatives with respect
#' to the \emph{argument}, from the ratio identity
#' \eqn{(\log K_\nu)' = \nu/x - K_{\nu+1}/K_\nu} and the modified Bessel
#' equation, exactly as in \code{\link{log_bessel_i_derivs}}.
#'
#' @param x A numeric vector of positive arguments, recycled against
#'   \code{nu}.
#' @param nu A numeric vector of orders, any sign.
#'
#' @return A list with the value \code{l} and the derivatives \code{d1} to
#'   \code{d4} in the argument.
#'
#' @references
#' Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
#' computation of the logarithm of modified Bessel functions on GPUs.
#' \emph{Proceedings of the 38th ACM International Conference on
#' Supercomputing (ICS '24)}. arXiv:2409.08729.
#'
#' @seealso \code{\link{log_bessel_k}}, \code{\link{log_bessel_i_derivs}}
#'
#' @examples
#' # d/dx log K_{1/2}(x) is exactly -1/(2x) - 1
#' log_bessel_k_derivs(2, 0.5)$d1
#'
#' @export
log_bessel_k_derivs <- function(x, nu) {
  nn <- max(length(x), length(nu))
  x <- rep_len(x, nn)
  v <- abs(rep_len(nu, nn))
  l0 <- log_bessel_k(x, v)
  r1 <- v / x - exp(log_bessel_k(x, v + 1) - l0)
  c(list(l = l0), .lb_deriv_chain(x, v, r1))
}
