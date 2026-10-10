#' @include quadrature.R
NULL

# Special functions the toolkit's distributions are written in, each carrying
# the overflow discipline learned on it: the Mills ratio on the log scale,
# Owen's T through one batched quadrature, the Bessel ratio through its
# series and continued fraction (src/bessel_ratio.cpp).

#' The Mills Ratio and Its Derivative
#'
#' @description
#' Returns \eqn{R(t) = \phi(t)/\Phi(t)}, often called the inverse Mills ratio,
#' and its derivative \eqn{R'(t) = -R(t)\{t + R(t)\}}. Every derivative of a
#' skew normal log-density is built from these two quantities.
#'
#' @details
#' The ratio is formed on the log scale. Written directly it is \eqn{0/0} for
#' \eqn{t} below about \eqn{-38}, where both the density and the distribution
#' function underflow, while the ratio itself is finite there and close to
#' \eqn{-t}. The identity for \eqn{R'} follows from differentiating the
#' quotient and using \eqn{\phi'(t) = -t\phi(t)}.
#'
#' @param t A numeric vector. Large finite values of either sign are handled
#'   on the log scale and are not special-cased. An infinite value gives `NaN`
#'   in `dr`, and `-Inf` gives `NaN` in `r` as well.
#'
#' @return A list of two numeric vectors, each the length of `t`:
#'   \describe{
#'     \item{`r`}{the ratio \eqn{R(t) = \phi(t)/\Phi(t)}, positive and
#'       decreasing, asymptotic to \eqn{-t} as \eqn{t \to -\infty}. In double
#'       precision it underflows to zero above about \eqn{t = 38.6}.}
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
#'   \dfrac{e^{-h^2(1 + x^2)/2}}{1 + x^2}\,\mathrm{d}x}, the function in which
#' the skew normal distribution function is written. \eqn{T(h, a)} is the probability
#' that a pair of independent standard normal variables falls in the wedge
#' below the line of
#' slope \eqn{a} beyond \eqn{h}, so it is bounded by \eqn{1/4} and is odd in
#' \eqn{a}.
#'
#' @details
#' # The integral
#'
#' \eqn{T} is even in \eqn{h} and odd in \eqn{a}, so the computation runs on
#' \eqn{\lvert h\rvert} and \eqn{\lvert a\rvert}. For \eqn{\lvert a\rvert \le 1}
#' the factor \eqn{e^{-h^2/2}} is taken out of the integral,
#'
#' \deqn{T(h, a) = \frac{e^{-h^2/2}}{2\pi} \int_0^{a}
#'   \frac{e^{-h^2 x^2/2}}{1 + x^2}\,\mathrm{d}x,}
#'
#' so the integrand equals one at zero whatever \eqn{h} is, and [quad_vec()]
#' evaluates the integral to a relative tolerance of \eqn{10^{-13}} with no
#' absolute tolerance, which keeps the relative accuracy of a value far below
#' one. Every such element goes into one batched call, one row per element, so
#' a whole vector of skew normal probabilities costs a single quadrature.
#'
#' # A slope above one
#'
#' For \eqn{a > 1} the integrand is concentrated within a distance of order
#' \eqn{1/h} of zero, which a quadrature over \eqn{[0, a]} can miss. The
#' reflection identity (Owen, 1956)
#'
#' \deqn{T(h, a) = \tfrac{1}{2}\bigl\{Q(h) + Q(ah)\bigr\} - Q(h)\,Q(ah)
#'   - T(ah, 1/a), \qquad h \ge 0,}
#'
#' with \eqn{Q = 1 - \Phi} computed as an upper tail, replaces it by an integral
#' with a slope below one. The result is at least a quarter of the first
#' terms, so the subtraction loses at most two bits.
#'
#' # Closed forms
#'
#' \eqn{T(h, 0) = 0}, \eqn{T(h, \infty) = \tfrac{1}{2}Q(\lvert h\rvert)} and
#' \eqn{T(\pm\infty, a) = 0} are set directly, and a missing argument gives
#' `NA`.
#'
#' Against values computed to 50 digits on a grid of \eqn{h} from 0 to 37 and
#' \eqn{a} from \eqn{10^{-6}} to \eqn{10^6}, the largest relative error is of
#' order \eqn{10^{-15}}.
#'
#' @param h A numeric vector, the offset, of any sign and size.
#' @param a A numeric vector of slopes, recycled against `h`, of any sign and
#'   size, `Inf` included.
#'
#' @return A numeric vector of the recycled length of `h` and `a`, bounded in
#'   \eqn{[-1/4, 1/4]}, with `NA` where either argument is missing.
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
#' # A steep slope, where the integrand sits within 1/h of zero: the skew
#' # normal distribution function at alpha = 1000 against its limit
#' # 2 Phi(z) - 1.
#' z <- 2.5
#' (pnorm(z) - 2 * owen_t(z, 1000)) - (2 * pnorm(z) - 1)
#'
#' @seealso [mills_ratio()], [bessel_i_ratio()], [log_bessel_i()], [log_bessel_k()]
#' @export
owen_t <- function(h, a) {
  n <- max(length(h), length(a))
  h <- abs(rep_len(as.numeric(h), n))
  a <- rep_len(as.numeric(a), n)

  out <- rep(NA_real_, n)
  ok <- !is.na(h) & !is.na(a)
  sgn <- sign(a)
  aa <- abs(a)

  out[ok & aa == 0] <- 0
  out[ok & is.infinite(h) & aa > 0] <- 0
  inf_a <- ok & is.finite(h) & is.infinite(aa)
  out[inf_a] <- sgn[inf_a] * 0.5 * stats::pnorm(h[inf_a], lower.tail = FALSE)

  fin <- ok & is.finite(h) & is.finite(aa) & aa > 0
  small <- which(fin & aa <= 1)
  big <- which(fin & aa > 1)
  # a slope above one is reflected onto T(a h, 1/a), so every integral runs
  # over a slope in (0, 1], where the scaled integrand starts at one and its
  # peak of width 1/h cannot fall between the nodes of the first panel
  hh <- c(h[small], aa[big] * h[big])
  sl <- c(aa[small], 1 / aa[big])
  if (length(hh)) {
    f <- function(x, i) {
      xv <- as.numeric(x)
      hv <- rep(hh[i], times = ncol(x))
      exp(-hv^2 * xv^2 / 2) / (1 + xv^2)
    }
    # no absolute budget: the integral is at least 0.3 min(a, 1/h), and an
    # absolute floor would cost the relative accuracy of a tiny T
    v <- exp(-hh^2 / 2) * quad_vec(f, 0, sl, atol = 0, rtol = 1e-13) / (2 * pi)
    ns <- length(small)
    out[small] <- sgn[small] * v[seq_len(ns)]
    if (length(big)) {
      q1 <- stats::pnorm(h[big], lower.tail = FALSE)
      q2 <- stats::pnorm(aa[big] * h[big], lower.tail = FALSE)
      out[big] <- sgn[big] * ((q1 + q2) / 2 - q1 * q2 - v[ns + seq_along(big)])
    }
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
#' and the moment that a method of moments estimates.
#'
#' @details
#' The ratio is evaluated in compiled code without the Bessel functions
#' themselves, which overflow from about \eqn{\kappa = 709} and, exponentially
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
#' backward recurrence, vectorized over \eqn{\kappa}. The von Mises
#' distribution function is a series in these ratios, and the recurrence
#' returns all \eqn{m} of them in one pass.
#'
#' @details
#' # The backward recurrence
#'
#' The three-term recurrence \eqn{I_{j-1} - I_{j+1} = (2j/\kappa) I_j} has two
#' solutions, one growing and one decaying. Run upwards it should follow the
#' decaying one and instead follows rounding error into the growing one, so it
#' is unstable. Run downwards the roles swap and it is stable, which is
#' Miller's algorithm. The ratios \eqn{r_j = I_j/I_{j-1}} satisfy
#' \eqn{r_j = 1/(2j/\kappa + r_{j+1})}, started from \eqn{r_{n_0+1} = 0} at
#' an index far enough above both \eqn{m} and \eqn{\kappa}; the result is
#' their running product, and the normalization by \eqn{I_0} requires no extra
#' work because the product starts there.
#'
#' # Cost
#'
#' The recurrence loop runs over the series index and is vectorized over the
#' data, so the number of steps depends on \eqn{m} and on the largest
#' \eqn{\kappa}, and not on the length of the vector. For \eqn{m > 1} the
#' running products are then formed row by row. A series over these ratios
#' therefore costs less than one quadrature per observation, and the von Mises
#' distribution function is evaluated this way.
#'
#' [bessel_i_ratio()] is the first of these ratios and switches to an
#' asymptotic series from \eqn{\kappa = 30}. No such branch is implemented
#' here: the recurrence needs a starting index above \eqn{\kappa}, so its cost
#' grows with the concentration, and at such concentrations a series in these
#' ratios does not converge in a useful number of terms.
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
#' order is accurate to the last bits wherever its value is a normal double;
#' at very large concentrations the value becomes subnormal and carries fewer
#' significant digits.
#'
#' \eqn{A'} is the variance of \eqn{\cos(\Theta - \mu)} under a von Mises
#' distribution and is therefore positive; the higher derivatives are its
#' cumulants of order three to five.
#'
#' @param kappa A numeric vector of concentrations, non-negative and of any
#'   size. A negative value returns `NaN`.
#' @inheritParams bessel_i_ratio
#'
#' @return A numeric vector the length of `kappa`. `bessel_i_ratio_d1()` is
#'   positive, equal to 1/2 at zero, and decays like \eqn{1/(2\kappa^2)}, so
#'   in double precision it underflows to zero beyond about
#'   \eqn{\kappa = 10^{154}}.
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
#' \eqn{A(\kappa) = I_1(\kappa)/I_0(\kappa)} equals \eqn{\rho}. It is the map
#' used by a von Mises method of moments, which converts an observed mean
#' resultant length into the concentration that produced it.
#'
#' @details
#' \eqn{A} has no elementary inverse, so \eqn{\kappa} is found by Newton's
#' method in compiled code, started from a piecewise approximation to the
#' inverse (a short series for \eqn{\rho < 0.53} and rational forms above).
#' \eqn{A} is increasing and concave, so after the first step every iterate
#' lies on the left of the root and rises to it; a step that would fall below
#' \eqn{2\rho} is replaced by \eqn{2\rho}, which is still on the left since
#' \eqn{A(\kappa) < \kappa/2}. The iteration ends when a step is no larger
#' than four times the machine epsilon times the iterate, or when an iterate
#' fails to rise above the previous one, and in any case after 200 iterations.
#'
#' For \eqn{\rho \ge 1/2} the residual \eqn{A(\kappa) - \rho} is formed as
#' \eqn{(1 - \rho) - (1 - A(\kappa))}, where \eqn{1 - \rho} is exact and
#' \eqn{1 - A} keeps its relative accuracy: it is computed from \eqn{A} in
#' double-double arithmetic below \eqn{\kappa = 30}, and from its own
#' asymptotic series above, without forming \eqn{A}. Near \eqn{\rho = 1} the
#' inverse behaves as \eqn{\kappa \approx 1/(2(1 - \rho))}, so its relative
#' error is that of \eqn{1 - \rho}: the result is the concentration whose
#' ratio is the double `rho`, to the last bits, at any concentration.
#'
#' @param rho A numeric vector of mean resultant lengths, strictly inside
#'   \eqn{(0, 1)}. Anything outside, the endpoints included, returns `NA`
#'   without a warning (at \eqn{\rho = 0} the inverse is the limit 0, and at
#'   \eqn{\rho = 1} it diverges).
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
#' # Outside the open unit interval, the endpoints included, the result is NA.
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
#' computed together in one evaluation. Since \eqn{A'} decays like
#' \eqn{1/(2\kappa^2)}, its powers in the denominators would underflow at
#' large concentrations, so the formulas are evaluated on \eqn{w = 1/A'} and
#' the ratios \eqn{r_j = A^{(j)}/A'},
#'
#' \deqn{\kappa'' = -r_2 w^2, \qquad \kappa''' = (3r_2^2 - r_3)\,w^3, \qquad
#'       \kappa'''' = (-15r_2^3 + 10r_2r_3 - r_4)\,w^4,}
#'
#' with the ratios taken from the asymptotic series without its powers of
#' \eqn{1/\kappa} from \eqn{\kappa = 30}. Each derivative is finite and
#' accurate until the derivative itself overflows (near
#' \eqn{\kappa = 10^{60}} for the fourth).
#'
#' @param kappa A numeric vector of concentrations, non-negative. A negative
#'   value returns `NaN`.
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
