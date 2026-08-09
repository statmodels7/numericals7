#' @include quadrature.R
NULL

# Special functions the toolkit's distributions are written in, each carrying
# the overflow discipline learned on it: the Mills ratio on the log scale,
# Owen's T through one batched quadrature, the Bessel ratio through the
# exponentially scaled Bessel functions.

#' The Mills Ratio and Its Derivative
#'
#' @description
#' Returns \eqn{R(t) = \phi(t)/\Phi(t)} and \eqn{R'(t) = -R(t)\{t + R(t)\}},
#' the two quantities every derivative of a skew normal log-density is built
#' from.
#'
#' @details
#' The ratio is formed on the log scale. Written directly it is \eqn{0/0} for
#' \eqn{t} below about \eqn{-38}, where both the density and the distribution
#' function underflow, while the ratio itself is finite there and close to
#' \eqn{-t}. The identity for \eqn{R'} follows from differentiating the
#' quotient and using \eqn{\phi'(t) = -t\phi(t)}.
#'
#' @param t A numeric vector.
#'
#' @return A list with \code{r} and \code{dr}.
#'
#' @examples
#' mills_ratio(c(-40, 0, 3))$r
#'
#' @seealso \code{\link{owen_t}}, \code{\link{bessel_i_ratio}}, \code{\link{log_bessel_i}}, \code{\link{log_bessel_k}}
#' @export
mills_ratio <- function(t) {
  r <- exp(stats::dnorm(t, log = TRUE) - stats::pnorm(t, log.p = TRUE))
  list(r = r, dr = -r * (t + r))
}

#' Owen's T Function
#'
#' @description
#' Computes \eqn{T(h, a) = \dfrac{1}{2\pi}\displaystyle\int_0^{a}
#'   \dfrac{e^{-h^2(1 + x^2)/2}}{1 + x^2}\,dx}, which is what the skew normal
#' distribution function is written in.
#'
#' @details
#' The integrand is bounded and smooth over a finite range, so quadrature
#' evaluates it to near machine precision; every element of the input goes
#' into one batched call of \code{\link{quad_vec}}, one row per element. Two
#' identities keep the extremes exact rather than quadrature-bound:
#' \eqn{T(h, a) = -T(h, -a)}, and \eqn{T(h, \infty) = \tfrac{1}{2}\Phi(-|h|)}.
#'
#' @param h A numeric vector.
#' @param a A numeric vector, recycled against \code{h}.
#'
#' @return A numeric vector.
#'
#' @references
#' Owen, D. B. (1956). Tables for computing bivariate normal probabilities.
#' \emph{Annals of Mathematical Statistics} 27, 1075-1090.
#'
#' @examples
#' owen_t(0, 1) - atan(1) / (2 * pi)
#'
#' @seealso \code{\link{mills_ratio}}, \code{\link{bessel_i_ratio}}, \code{\link{log_bessel_i}}, \code{\link{log_bessel_k}}
#' @export
owen_t <- function(h, a) {
  n <- max(length(h), length(a))
  h <- rep_len(h, n)
  a <- rep_len(a, n)

  out <- numeric(n)
  sgn <- sign(a)
  aa <- abs(a)

  inf_a <- is.finite(h) & is.infinite(aa)
  out[inf_a] <- sgn[inf_a] * 0.5 * stats::pnorm(-abs(h[inf_a]))

  todo <- which(is.finite(h) & is.finite(aa) & aa > 0)
  if (length(todo)) {
    hv <- h[todo]
    f <- function(x, i) {
      xv <- as.numeric(x)
      hh <- rep(hv[i], times = ncol(x))
      exp(-hh^2 * (1 + xv^2) / 2) / (1 + xv^2)
    }
    v <- quad_vec(f, 0, aa[todo], atol = 1e-15, rtol = 1e-12)
    out[todo] <- sgn[todo] * v / (2 * pi)
  }
  out
}

#' The Ratio of Modified Bessel Functions
#'
#' @description
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}, a strictly increasing bijection
#' from \eqn{(0, \infty)} onto \eqn{(0, 1)}. For a von Mises distribution it
#' is the mean resultant length, the expected cosine of the deviation from
#' the mean direction.
#'
#' @details
#' Both Bessel functions are taken exponentially scaled, so the factor
#' \eqn{e^{\kappa}} they share cancels in the ratio where the unscaled
#' functions would overflow, from about \eqn{\kappa = 700}. The scaled
#' functions themselves underflow to an exact zero between \eqn{10^5} and
#' \eqn{10^6}, so past \eqn{\kappa = 10^4} the ratio is taken from its
#' asymptotic expansion
#' \eqn{1 - 1/(2\kappa) - 1/(8\kappa^2) - 1/(8\kappa^3)}, whose next term is
#' already below the resolution of a double at the switch point; the result
#' is therefore finite and accurate for an argument of any size.
#'
#' @param kappa A positive numeric vector.
#'
#' @return A numeric vector in \eqn{(0, 1)}.
#'
#' @seealso \code{\link{bessel_i_ratio_derivs}}, \code{\link{bessel_i_ratio_inverse}}
#'
#' @examples
#' bessel_i_ratio(c(0.5, 2, 1000))
#'
#' @export
bessel_i_ratio <- function(kappa) {
  out <- rep(NA_real_, length(kappa))
  small <- !is.na(kappa) & kappa < 1e4
  large <- !is.na(kappa) & kappa >= 1e4
  if (any(small)) {
    k <- kappa[small]
    out[small] <- besselI(k, 1, expon.scaled = TRUE) /
      besselI(k, 0, expon.scaled = TRUE)
  }
  if (any(large)) {
    k <- kappa[large]
    out[large] <- 1 - 1 / (2 * k) - 1 / (8 * k^2) - 1 / (8 * k^3)
  }
  out
}

#' Derivatives of the Bessel Ratio
#'
#' @description
#' \eqn{A(\kappa)} and its first four derivatives, obtained by
#' differentiating \eqn{A' = 1 - A/\kappa - A^2} repeatedly.
#'
#' @details
#' Each order is written in the orders below it, so the whole table costs the
#' two Bessel evaluations of \code{\link{bessel_i_ratio}} and nothing more.
#' The first identity follows from \eqn{I_0' = I_1} and
#' \eqn{I_1' = I_0 - I_1/\kappa}; the alternative, a Bessel function of
#' higher order per derivative, costs more and is less accurate at large
#' \eqn{\kappa}, where the functions themselves overflow and only their ratio
#' does not. \eqn{A'} is the variance of a cosine and therefore positive.
#'
#' @param kappa A positive numeric vector.
#'
#' @return A named list with \code{A} and its derivatives \code{d1} to
#'   \code{d4}.
#'
#' @seealso \code{\link{bessel_i_ratio}}, \code{\link{bessel_i_ratio_inverse}}
#'
#' @examples
#' bessel_i_ratio_derivs(2)$d1
#'
#' @export
bessel_i_ratio_derivs <- function(kappa) {
  k <- kappa
  A <- bessel_i_ratio(k)
  d1 <- 1 - A / k - A * A
  d2 <- -d1 / k + A / k^2 - 2 * A * d1
  d3 <- -d2 / k + 2 * d1 / k^2 - 2 * A / k^3 - 2 * d1^2 - 2 * A * d2
  d4 <- -d3 / k + 3 * d2 / k^2 - 6 * d1 / k^3 + 6 * A / k^4 -
    6 * d1 * d2 - 2 * A * d3
  list(A = A, d1 = d1, d2 = d2, d3 = d3, d4 = d4)
}

#' The Inverse of the Bessel Ratio
#'
#' @description
#' \eqn{\kappa = A^{-1}(\rho)} by root finding, together with the four
#' derivatives of the inverse in \eqn{\rho}.
#'
#' @details
#' \eqn{A} has no elementary inverse, so \eqn{\kappa} is found by root
#' finding from a bracket built around the standard series approximation and
#' widened until it straddles the root; the asymptotic branch of
#' \code{\link{bessel_i_ratio}} keeps the function evaluable over the whole
#' bracket, however concentrated. The derivatives come from the inverse
#' function rule on \code{\link{bessel_i_ratio_derivs}}:
#' \deqn{\kappa' = \dfrac{1}{A'}, \qquad
#'       \kappa'' = -\dfrac{A''}{(A')^3}, \qquad
#'       \kappa''' = \dfrac{3(A'')^2 - A'A'''}{(A')^5},}
#' and the fourth in the same pattern; \eqn{A' > 0} keeps every denominator
#' away from zero in the interior. A \code{rho} outside \eqn{(0, 1)} returns
#' \code{NA}.
#'
#' @param rho A numeric vector in \eqn{(0, 1)}.
#'
#' @return A named list with \code{kappa} and its derivatives \code{d1} to
#'   \code{d4} in \code{rho}.
#'
#' @seealso \code{\link{bessel_i_ratio}}, \code{\link{bessel_i_ratio_derivs}}
#'
#' @examples
#' bessel_i_ratio(bessel_i_ratio_inverse(0.7)$kappa)
#'
#' @export
bessel_i_ratio_inverse <- function(rho) {
  one <- function(r) {
    if (!is.finite(r) || r <= 0 || r >= 1) return(NA_real_)
    g <- if (r < 0.53) {
      2 * r + r^3 + 5 * r^5 / 6
    } else if (r < 0.85) {
      -0.4 + 1.39 * r + 0.43 / (1 - r)
    } else {
      1 / (r^3 - 4 * r^2 + 3 * r)
    }
    g <- min(max(g, 1e-8), 1e12)
    f <- function(k) bessel_i_ratio(k) - r
    lo <- g / 2
    hi <- g * 2
    it <- 0L
    while (f(lo) > 0 && lo > 1e-10 && it < 60L) {
      lo <- lo / 2
      it <- it + 1L
    }
    it <- 0L
    while (f(hi) < 0 && hi < 1e13 && it < 60L) {
      hi <- hi * 2
      it <- it + 1L
    }
    stats::uniroot(f, c(lo, hi), tol = .Machine$double.eps^0.75)$root
  }
  k <- vapply(rho, one, numeric(1))
  a <- bessel_i_ratio_derivs(k)
  p1 <- a$d1
  list(
    kappa = k,
    d1 = 1 / p1,
    d2 = -a$d2 / p1^3,
    d3 = (3 * a$d2^2 - p1 * a$d3) / p1^5,
    d4 = (-15 * a$d2^3 + 10 * p1 * a$d2 * a$d3 - p1^2 * a$d4) / p1^7
  )
}
