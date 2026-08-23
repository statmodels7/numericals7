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
#' @param t A numeric vector of any values, the whole real line included. No
#'   argument is out of range and none is special-cased.
#'
#' @return A list of two numeric vectors, each the length of `t`:
#'   \describe{
#'     \item{`r`}{the ratio \eqn{R(t) = \phi(t)/\Phi(t)}, positive and
#'       decreasing, asymptotic to \eqn{-t} as \eqn{t \to -\infty}.}
#'     \item{`dr`}{its derivative \eqn{R'(t)}, which lies in \eqn{(-1, 0)}.}
#'   }
#'
#' @examples
#' mills_ratio(c(-5, 0, 3))$r
#'
#' # The point of the log-scale form. Written directly the ratio is 0/0 well
#' # inside the range a skew normal reaches, while the true value is finite
#' # and close to -t.
#' dnorm(-40) / pnorm(-40)
#' mills_ratio(-400)$r + (-400)
#'
#' # The derivative is the stated identity, exactly.
#' m <- mills_ratio(c(-5, 0, 3))
#' max(abs(m$dr - (-m$r * (c(-5, 0, 3) + m$r))))
#'
#' @seealso [owen_t()], [bessel_i_ratio()], [log_bessel_i()], [log_bessel_k()]
#' @export
mills_ratio <- function(t) {
  r <- exp(stats::dnorm(t, log = TRUE) - stats::pnorm(t, log.p = TRUE))
  list(r = r, dr = -r * (t + r))
}

#' Owen's T Function
#'
#' @description
#' Computes \eqn{T(h, a) = \dfrac{1}{2\pi}\displaystyle\int_0^{a}
#'   \dfrac{e^{-h^2(1 + x^2)/2}}{1 + x^2}\,\mathrm{d}x}, the function the skew
#' normal distribution function is written in. \eqn{T(h, a)} is the probability
#' that a standard bivariate normal pair falls in the wedge below the line of
#' slope \eqn{a} beyond \eqn{h}, so it is bounded by \eqn{1/4} and is odd in
#' \eqn{a}.
#'
#' @details
#' The integrand is bounded and smooth over a finite range, so quadrature
#' evaluates it to near machine precision. Every element of the input goes into
#' one batched call of [quad_vec()], one row per element, so a whole vector of
#' skew normal probabilities costs a single quadrature.
#'
#' Two identities are applied in closed form, so the extremes are exact where
#' quadrature would merely be accurate: \eqn{T(h, a) = -T(h, -a)} handles a
#' negative second argument, and \eqn{T(h, \infty) = \tfrac{1}{2}\Phi(-|h|)}
#' handles an infinite one.
#'
#' @param h A numeric vector, the offset. Any finite value.
#' @param a A numeric vector of slopes, recycled against `h`. May be negative or
#'   infinite; both are taken by identity.
#'
#' @return A numeric vector of the recycled length of `h` and `a`, bounded in
#'   \eqn{[-1/4, 1/4]}.
#'
#' @references
#' Owen, D. B. (1956). Tables for computing bivariate normal probabilities.
#' *Annals of Mathematical Statistics* 27, 1075-1090.
#'
#' @examples
#' # At h = 0 the integral is elementary: T(0, a) = atan(a) / (2 pi).
#' a <- c(0.5, 1, 4)
#' max(abs(owen_t(0, a) - atan(a) / (2 * pi)))
#'
#' # Odd in the second argument, and the infinite case is a normal tail.
#' owen_t(1, 2) + owen_t(1, -2)
#' owen_t(1.3, Inf) - pnorm(-1.3) / 2
#'
#' @seealso [mills_ratio()], [bessel_i_ratio()], [log_bessel_i()], [log_bessel_k()]
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
#' Computes \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)}, a strictly increasing
#' bijection from \eqn{(0, \infty)} onto \eqn{(0, 1)}. For a von Mises
#' distribution it is the mean resultant length, the expected cosine of the
#' deviation from the mean direction, so it is the map between a concentration
#' and the moment a method of moments estimates.
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
#' @param kappa A numeric vector of concentrations, positive and of any size.
#'   Zero returns 0, the limit. A negative value returns `NaN`.
#'
#' @return A numeric vector the length of `kappa`, in \eqn{(0, 1)} and
#'   increasing in its argument.
#'
#' @seealso [bessel_i_ratio_derivs()] for its derivatives,
#'   [bessel_i_ratio_inverse()] for the map back, [bessel_i_ratios()] for the
#'   sequence of higher orders.
#'
#' @examples
#' bessel_i_ratio(c(0.5, 2, 1000))
#'
#' # It agrees with the scaled Bessel functions where those still evaluate.
#' k <- c(0.5, 2, 1e3, 1e4)
#' max(abs(bessel_i_ratio(k) - besselI(k, 1, TRUE) / besselI(k, 0, TRUE)))
#'
#' # Past that the scaled functions underflow to zero and their ratio is NaN,
#' # while the asymptotic branch carries the answer to any concentration.
#' suppressWarnings(besselI(1e6, 1, TRUE) / besselI(1e6, 0, TRUE))
#' bessel_i_ratio(c(1e6, 1e12))
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


#' The Sequence of Modified Bessel Ratios
#'
#' @description
#' Computes \eqn{I_j(\kappa)/I_0(\kappa)} for \eqn{j = 1, \dots, m} by Miller's
#' backward recurrence, vectorized over \eqn{\kappa}. A series in these ratios
#' is what a von Mises distribution function costs, so getting all \eqn{m} of
#' them for the price of one matters.
#'
#' @details
#' # Why the recurrence runs backwards
#'
#' The three-term recurrence \eqn{I_{j-1} - I_{j+1} = (2j/\kappa) I_j} has two
#' solutions, one growing and one decaying. Run upwards it should follow the
#' decaying one and instead follows rounding error into the growing one, so it
#' is unstable. Run downwards the roles swap and it is stable, which is
#' Miller's algorithm. The ratios \eqn{r_j = I_j/I_{j-1}} satisfy
#' \eqn{r_j = 1/(2j/\kappa + r_{j+1})}, started from \eqn{r_{n_0+1} = 0} at
#' an index far enough above both \eqn{m} and \eqn{\kappa}; the answer is
#' their running product, and the normalization by \eqn{I_0} is free because
#' the product starts there.
#'
#' # Cost
#'
#' The loop runs over the series index, never over the data, so a vector of
#' \eqn{\kappa} costs the same number of vectorized steps as a single value.
#' A series over these ratios therefore costs less than a quadrature per
#' observation, which is why the von Mises distribution function stopped being
#' one.
#'
#' [bessel_i_ratio()] is the first of them and carries an asymptotic branch past
#' \eqn{\kappa = 10^4}, where the scaled Bessel functions underflow. There is no
#' such branch here, and none is wanted: the recurrence needs a starting index
#' above \eqn{\kappa}, so its cost grows with the concentration, and a caller
#' that far out is already past the point where a series in these ratios
#' converges in any useful number of terms.
#'
#' @param kappa A numeric vector of concentrations, positive.
#' @param m How many ratios to return, a positive whole number. It sets the
#'   number of columns and, with `kappa`, the starting index of the recurrence.
#'
#' @return A numeric matrix of `length(kappa)` rows and `m` columns. Entry
#'   \eqn{(i, j)} is \eqn{I_j(\kappa_i)/I_0(\kappa_i)}, decreasing along a row.
#'
#' @seealso [bessel_i_ratio()] for the first ratio alone, [log_bessel_i()] for
#'   the functions themselves.
#'
#' @examples
#' # Four ratios at two concentrations, one row each.
#' r <- bessel_i_ratios(c(1, 5), 4)
#' round(r, 6)
#'
#' # They agree with the scaled Bessel functions to the last bit.
#' r[2L, ] - besselI(5, 1:4, TRUE) / besselI(5, 0, TRUE)
#'
#' # And decrease along a row: a higher order is a smaller ratio.
#' all(diff(r[2L, ]) < 0)
#'
#' @export
bessel_i_ratios <- function(kappa, m) {
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 1L) {
    stop("'m' must be a single positive integer.", call. = FALSE)
  }
  kappa <- as.numeric(kappa)
  if (!length(kappa)) return(matrix(numeric(0), 0L, m))
  if (any(!is.na(kappa) & kappa <= 0)) {
    stop("'kappa' must be positive.", call. = FALSE)
  }
  kmax <- suppressWarnings(max(kappa[is.finite(kappa)], 0))
  # far enough above both the order asked for and the argument: below either
  # the downward recurrence has not yet forgotten its starting value
  n0 <- m + max(30L, ceiling(sqrt(40 * m)), ceiling(kmax))
  r <- numeric(length(kappa))
  out <- matrix(NA_real_, length(kappa), m)
  ok <- is.finite(kappa) & kappa > 0
  k <- kappa[ok]
  rj <- numeric(length(k))
  keep <- matrix(0, length(k), m)
  for (j in seq.int(n0, 1L)) {
    rj <- 1 / (2 * j / k + rj)
    if (j <= m) keep[, j] <- rj
  }
  if (m > 1L) keep <- t(apply(keep, 1L, cumprod))
  out[ok, ] <- keep
  out
}

#' Derivatives of the Bessel Ratio
#'
#' @description
#' Computes \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} and its first four
#' derivatives, by differentiating the identity \eqn{A' = 1 - A/\kappa - A^2}
#' repeatedly. A von Mises family needs all five to reach fourth-order
#' derivatives in its concentration.
#'
#' @details
#' Each order is written in the orders below it, so the whole table costs the
#' two Bessel evaluations of [bessel_i_ratio()] and nothing more.
#' The first identity follows from \eqn{I_0' = I_1} and
#' \eqn{I_1' = I_0 - I_1/\kappa}; the alternative, a Bessel function of
#' higher order per derivative, costs more and is less accurate at large
#' \eqn{\kappa}, where the functions themselves overflow and only their ratio
#' does not. \eqn{A'} is the variance of a cosine and therefore positive.
#'
#' @param kappa A numeric vector of concentrations, positive.
#'
#' @return A named list of five numeric vectors, each the length of `kappa`:
#'   `A`, the ratio itself, and `d1` to `d4`, its derivatives in \eqn{\kappa}.
#'   `d1` is strictly positive, being a variance.
#'
#' @seealso [bessel_i_ratio()] for the value alone,
#'   [bessel_i_ratio_inverse()] for the derivatives of the inverse map.
#'
#' @examples
#' str(bessel_i_ratio_derivs(2))
#'
#' # The first derivative is the variance of a cosine, so it is positive at
#' # every concentration and vanishes as the distribution concentrates.
#' bessel_i_ratio_derivs(c(0.1, 1, 100))$d1
#'
#' # It satisfies the identity the whole table is built from.
#' d <- bessel_i_ratio_derivs(2)
#' d$d1 - (1 - d$A / 2 - d$A^2)
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
#' Computes \eqn{\kappa = A^{-1}(\rho)} by root finding, together with the four
#' derivatives of the inverse in \eqn{\rho}. This is the map a von Mises method
#' of moments runs: it turns an observed mean resultant length back into the
#' concentration that produced it.
#'
#' @details
#' \eqn{A} has no elementary inverse, so \eqn{\kappa} is found by root
#' finding from a bracket built around the standard series approximation and
#' widened until it straddles the root; the asymptotic branch of
#' [bessel_i_ratio()] keeps the function evaluable over the whole
#' bracket, however concentrated. The derivatives come from the inverse
#' function rule on [bessel_i_ratio_derivs()]:
#' \deqn{\kappa' = \dfrac{1}{A'}, \qquad
#'       \kappa'' = -\dfrac{A''}{(A')^3}, \qquad
#'       \kappa''' = \dfrac{3(A'')^2 - A'A'''}{(A')^5},}
#' and the fourth in the same pattern; \eqn{A' > 0} keeps every denominator
#' away from zero in the interior. A `rho` outside \eqn{(0, 1)} returns
#' `NA`.
#'
#' @param rho A numeric vector of mean resultant lengths, strictly inside
#'   \eqn{(0, 1)}. Anything outside, the endpoints included, returns `NA`
#'   without a warning, the inverse having no finite value there.
#'
#' @return A named list of five numeric vectors, each the length of `rho`:
#'   `kappa`, the concentration, and `d1` to `d4`, the derivatives of the
#'   inverse in \eqn{\rho}. `NA` wherever `rho` left \eqn{(0, 1)}.
#'
#' @seealso [bessel_i_ratio()] for the forward map,
#'   [bessel_i_ratio_derivs()] for the derivatives it inverts.
#'
#' @examples
#' # The round trip closes to machine precision across the range.
#' rho <- c(0.1, 0.5, 0.99)
#' bessel_i_ratio(bessel_i_ratio_inverse(rho)$kappa) - rho
#'
#' # The first derivative is the reciprocal of A', the inverse function rule.
#' inv <- bessel_i_ratio_inverse(0.7)
#' inv$d1 - 1 / bessel_i_ratio_derivs(inv$kappa)$d1
#'
#' # Outside the open unit interval there is no concentration to return.
#' bessel_i_ratio_inverse(c(0, 0.5, 1))$kappa
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
