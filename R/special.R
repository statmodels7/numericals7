#' @include quadrature.R
NULL

# Special functions the toolkit's distributions are written in, each carrying
# the overflow discipline learned on it: the Mills ratio on the log scale,
# Owen's T through one batched quadrature, the Bessel ratio through its
# series and continued fraction (src/bessel_ratio.cpp).

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
#' The ratio is evaluated in compiled code without the Bessel functions
#' themselves, which overflow from about \eqn{\kappa = 700} and, exponentially
#' scaled, underflow between \eqn{10^5} and \eqn{10^6}. Below
#' \eqn{\kappa = 0.5} it is the power series of \eqn{A} at zero; from there to
#' \eqn{\kappa = 30} it is the continued fraction
#' \eqn{A = \kappa/(2 + \kappa^2/(4 + \kappa^2/(6 + \cdots)))}, evaluated
#' backwards from the index \eqn{\lfloor\kappa\rfloor + 20}; above, it is the
#' asymptotic series of \eqn{A} in \eqn{1/\kappa} with 30 terms, the quotient
#' of the asymptotic series of \eqn{I_1} and \eqn{I_0}. The result is finite
#' and accurate to the last bits for an argument of any size.
#'
#' @param kappa A numeric vector of concentrations, non-negative and of any
#'   size. Zero returns 0, the limit, and `Inf` returns 1. A negative value
#'   returns `NaN`.
#' @param threads The number of threads, a positive whole number.
#'
#' @return A numeric vector the length of `kappa`, in \eqn{[0, 1]} and
#'   increasing in its argument.
#'
#' @seealso [bessel_i_ratio_d1()] for its derivatives,
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
#' # while the asymptotic series carries the answer to any concentration.
#' suppressWarnings(besselI(1e6, 1, TRUE) / besselI(1e6, 0, TRUE))
#' bessel_i_ratio(c(1e6, 1e12))
#'
#' @export
bessel_i_ratio <- function(kappa, threads = 1L) {
  bessel_ratio_cpp(as.numeric(kappa), as.integer(threads))
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
#' [bessel_i_ratio()] is the first of them and carries an asymptotic series
#' from \eqn{\kappa = 30}. There is no such branch here, and none is
#' wanted: the recurrence needs a starting index
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
#' Compute the derivatives of \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} in
#' \eqn{\kappa}, one function per order: `bessel_i_ratio_d1()` returns
#' \eqn{A'}, `bessel_i_ratio_d2()` \eqn{A''}, and so on to the fourth. Each
#' computes its own order and the orders below it that the computation
#' needs, and nothing above.
#'
#' @details
#' The derivatives follow from the identity \eqn{A' = 1 - A/\kappa - A^2},
#' a consequence of \eqn{I_0' = I_1} and \eqn{I_1' = I_0 - I_1/\kappa},
#' differentiated repeatedly:
#' \deqn{A'' = -\frac{A'}{\kappa} + \frac{A}{\kappa^2} - 2AA', \qquad
#'       A''' = -\frac{A''}{\kappa} + \frac{2A'}{\kappa^2} -
#'       \frac{2A}{\kappa^3} - 2(A')^2 - 2AA'',}
#' and the fourth in the same pattern. In double precision this identity
#' cancels at both ends of the range: at small \eqn{\kappa} its terms are of
#' order \eqn{\kappa^{-n}} at derivative \eqn{n} while the result is of order
#' one or \eqn{\kappa}, and at large \eqn{\kappa} the derivative \eqn{A'} is
#' of order \eqn{\kappa^{-2}} against terms of order one. The functions
#' therefore use three regimes, as [bessel_i_ratio()] does: below
#' \eqn{\kappa = 0.5} the power series of \eqn{A} at zero differentiated term
#' by term; from 0.5 to 30 the continued fraction for \eqn{A} and the
#' identity above, both in double-double arithmetic (about 32 digits); from 30
#' the asymptotic series in \eqn{1/\kappa} differentiated term by term. Each
#' order is accurate to the last bits over the whole range.
#'
#' \eqn{A'} is the variance of \eqn{\cos(\Theta - \mu)} under a von Mises
#' distribution and is therefore positive; the higher derivatives are its
#' cumulants of order three to five.
#'
#' @inheritParams bessel_i_ratio
#'
#' @return A numeric vector the length of `kappa`. `bessel_i_ratio_d1()` is
#'   strictly positive at every finite concentration and tends to 1/2 at zero.
#'
#' @seealso [bessel_i_ratio()] for the value,
#'   [bessel_i_ratio_inverse_d1()] for the derivatives of the inverse map.
#'
#' @examples
#' bessel_i_ratio_d1(c(0.1, 1, 100))
#'
#' # The first derivative satisfies the identity the others are built from.
#' k <- 2
#' bessel_i_ratio_d1(k) - (1 - bessel_i_ratio(k) / k - bessel_i_ratio(k)^2)
#'
#' # At a small concentration the series keeps every digit: A''' tends to
#' # -3/8 and A'''' to 5 kappa / 4.
#' c(bessel_i_ratio_d3(1e-6), bessel_i_ratio_d4(1e-6) / 1e-6)
#'
#' @name bessel_i_ratio_d1
NULL

#' @rdname bessel_i_ratio_d1
#' @export
bessel_i_ratio_d1 <- function(kappa, threads = 1L) {
  bessel_ratio_d1_cpp(as.numeric(kappa), as.integer(threads))
}

#' @rdname bessel_i_ratio_d1
#' @export
bessel_i_ratio_d2 <- function(kappa, threads = 1L) {
  bessel_ratio_d2_cpp(as.numeric(kappa), as.integer(threads))
}

#' @rdname bessel_i_ratio_d1
#' @export
bessel_i_ratio_d3 <- function(kappa, threads = 1L) {
  bessel_ratio_d3_cpp(as.numeric(kappa), as.integer(threads))
}

#' @rdname bessel_i_ratio_d1
#' @export
bessel_i_ratio_d4 <- function(kappa, threads = 1L) {
  bessel_ratio_d4_cpp(as.numeric(kappa), as.integer(threads))
}

#' The Inverse of the Bessel Ratio
#'
#' @description
#' Computes \eqn{\kappa = A^{-1}(\rho)}, the concentration whose ratio
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} equals \eqn{\rho}. This is the map
#' a von Mises method of moments runs: it turns an observed mean resultant
#' length back into the concentration that produced it.
#'
#' @details
#' \eqn{A} has no elementary inverse, so \eqn{\kappa} is found by Newton's
#' method from the standard series approximation, in compiled code. \eqn{A}
#' is increasing and concave, so after the first step every iterate lies on
#' the left of the root and rises to it; a step that would fall below
#' \eqn{2\rho} is replaced by it, which is still on the left since
#' \eqn{A(\kappa) < \kappa/2}. The iteration ends where a step is no larger
#' than the spacing of the doubles at the iterate.
#'
#' For \eqn{\rho \ge 1/2} the residual \eqn{A(\kappa) - \rho} is formed as
#' \eqn{(1 - \rho) - (1 - A(\kappa))}, where \eqn{1 - \rho} is exact and
#' \eqn{1 - A} is computed without forming \eqn{A}. Near \eqn{\rho = 1} the
#' inverse behaves as \eqn{\kappa \approx 1/(2(1 - \rho))}, so its relative
#' error is that of \eqn{1 - \rho}: the result is the concentration whose
#' ratio is the double `rho`, to the last bits, at any concentration.
#'
#' @param rho A numeric vector of mean resultant lengths, strictly inside
#'   \eqn{(0, 1)}. Anything outside, the endpoints included, returns `NA`
#'   without a warning, the inverse having no finite value there.
#' @inheritParams bessel_i_ratio
#'
#' @return A numeric vector of concentrations the length of `rho`, `NA`
#'   wherever `rho` left \eqn{(0, 1)}.
#'
#' @seealso [bessel_i_ratio()] for the forward map,
#'   [bessel_i_ratio_inverse_d1()] for the derivatives of the inverse.
#'
#' @examples
#' # The round trip closes to machine precision across the range.
#' rho <- c(0.1, 0.5, 0.99)
#' bessel_i_ratio(bessel_i_ratio_inverse(rho)) - rho
#'
#' # Outside the open unit interval there is no concentration to return.
#' bessel_i_ratio_inverse(c(0, 0.5, 1))
#'
#' @export
bessel_i_ratio_inverse <- function(rho, threads = 1L) {
  bessel_ratio_inverse_cpp(as.numeric(rho), as.integer(threads))
}

#' Derivatives of the Inverse Bessel Ratio
#'
#' @description
#' Compute the derivatives of the inverse map \eqn{\kappa(\rho) =
#' A^{-1}(\rho)} in \eqn{\rho}, one function per order, evaluated at
#' \eqn{\rho = A(\kappa)}. They take the concentration rather than \eqn{\rho},
#' so a caller that has already inverted \eqn{\rho} does not invert it again.
#'
#' @details
#' The derivatives come from the inverse function rule on the derivatives of
#' \eqn{A} ([bessel_i_ratio_d1()] and the following orders):
#' \deqn{\kappa' = \frac{1}{A'}, \qquad
#'       \kappa'' = -\frac{A''}{(A')^3}, \qquad
#'       \kappa''' = \frac{3(A'')^2 - A'A'''}{(A')^5},}
#' \deqn{\kappa'''' = \frac{-15(A'')^3 + 10A'A''A''' - (A')^2A''''}{(A')^7}.}
#' The derivative of order \eqn{n} needs \eqn{A'} to \eqn{A^{(n)}}, which are
#' computed together at the cost of one evaluation of the continued fraction.
#' \eqn{A' > 0} keeps every denominator away from zero.
#'
#' @inheritParams bessel_i_ratio
#'
#' @return A numeric vector the length of `kappa`: the derivative of the
#'   inverse map at \eqn{\rho = A(\kappa)}.
#'
#' @seealso [bessel_i_ratio_inverse()] for the inverse itself.
#'
#' @examples
#' k <- bessel_i_ratio_inverse(0.7)
#'
#' # The first derivative is the reciprocal of A', the inverse function rule.
#' bessel_i_ratio_inverse_d1(k) - 1 / bessel_i_ratio_d1(k)
#'
#' # The second against a central difference of the first in rho.
#' h <- 1e-5
#' c(bessel_i_ratio_inverse_d2(k),
#'   (bessel_i_ratio_inverse_d1(bessel_i_ratio_inverse(0.7 + h)) -
#'      bessel_i_ratio_inverse_d1(bessel_i_ratio_inverse(0.7 - h))) / (2 * h))
#'
#' @name bessel_i_ratio_inverse_d1
NULL

#' @rdname bessel_i_ratio_inverse_d1
#' @export
bessel_i_ratio_inverse_d1 <- function(kappa, threads = 1L) {
  bessel_ratio_inverse_d1_cpp(as.numeric(kappa), as.integer(threads))
}

#' @rdname bessel_i_ratio_inverse_d1
#' @export
bessel_i_ratio_inverse_d2 <- function(kappa, threads = 1L) {
  bessel_ratio_inverse_d2_cpp(as.numeric(kappa), as.integer(threads))
}

#' @rdname bessel_i_ratio_inverse_d1
#' @export
bessel_i_ratio_inverse_d3 <- function(kappa, threads = 1L) {
  bessel_ratio_inverse_d3_cpp(as.numeric(kappa), as.integer(threads))
}

#' @rdname bessel_i_ratio_inverse_d1
#' @export
bessel_i_ratio_inverse_d4 <- function(kappa, threads = 1L) {
  bessel_ratio_inverse_d4_cpp(as.numeric(kappa), as.integer(threads))
}
