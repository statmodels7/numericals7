#' @include enumerate.R
NULL

# The one stencil library. Before this file the toolkit carried three
# independent finite-difference implementations: linkfunctions7's four
# hard-coded central stencils, basis7's Vandermonde weights, and
# distributions7's five-point fd5_* family. They were all instances of one
# construction -- solve a Vandermonde system for the weights, apply ONE
# stencil, never compose lower-order differences -- and now they are.

#' Finite-Difference Weights for Any Stencil
#'
#' @description
#' Solves for the weights that combine function values at the given offsets into
#' an estimate of a derivative. The weights are the unique solution of the
#' Vandermonde system that makes the stencil exact on polynomials, so a stencil
#' on `n` nodes reproduces the derivative of every polynomial of degree `n - 1`
#' or less without error. They are returned for a **unit step**: divide by
#' `h^order` before using them at spacing `h`, or call [fd_derivative()], which
#' assembles the whole expression.
#'
#' @details
#' # The system solved
#'
#' Write \eqn{s_1, \dots, s_n} for the offsets and \eqn{d} for the order. The
#' weights \eqn{w} solve the \eqn{n} moment conditions
#'
#' \deqn{\sum_{j=1}^{n} w_j\, s_j^{\,m} \;=\; d!\;[\,m = d\,],
#'       \qquad m = 0, 1, \dots, n-1,}
#'
#' one equation per polynomial degree the stencil must reproduce. The matrix is
#' Vandermonde in the offsets, so it is nonsingular exactly when they are
#' distinct and the weights exist and are unique.
#'
#' Two properties of the answer are cheap to check against a returned vector.
#' The weights sum to zero for every \eqn{d \ge 1}, this being the \eqn{m = 0}
#' condition. On a symmetric stencil they are antisymmetric for odd \eqn{d} and
#' symmetric for even \eqn{d}.
#'
#' # Applying them
#'
#' At a step \eqn{h} the derivative estimate is
#'
#' \deqn{f^{(d)}(x) \;\approx\; h^{-d} \sum_{j=1}^{n} w_j\, f(x + s_j h).}
#'
#' The factor \eqn{h^{-d}} belongs to the caller. The weights carry no step, so
#' one solve serves every step size. [fd_step()] supplies a step balanced for
#' the order and [fd_derivative()] puts the three together.
#'
#' # Accuracy
#'
#' Because the stencil is exact to degree \eqn{n-1}, the leading error is the
#' first term it cannot reproduce, giving an estimate accurate to
#' \eqn{O(h^{\,n-d})} with constant
#'
#' \deqn{\frac{1}{n!}\sum_{j=1}^{n} w_j\, s_j^{\,n}.}
#'
#' Five nodes therefore give a fourth-order first derivative and a second-order
#' fourth derivative. [fd_offsets()] sizes a stencil from the order and the
#' accuracy asked of it.
#'
#' # One stencil, never nested
#'
#' Reaching a high order by composing low-order differences multiplies the error
#' of each stage into the next, and a fourth derivative built from four nested
#' first differences is noise. Every numerical fallback in the toolkit takes one
#' stencil of the order it wants, and these weights are what it takes.
#'
#' @param offsets A numeric vector of distinct offsets, in units of the step.
#'   They need not be sorted, symmetric, or whole numbers: `0:2` gives a
#'   one-sided stencil for use at a boundary, and `c(-1, 0, 3)` an uneven one.
#'   At least `order + 1` of them are needed. Duplicated offsets throw an error,
#'   the Vandermonde system being singular.
#' @param order The derivative order, **strictly smaller than
#'   `length(offsets)`**; an order at or above the node count throws an error
#'   naming both numbers. Order `0` is legal and returns interpolation weights
#'   at the origin. A negative or fractional order is not checked and has no
#'   defined meaning: a negative one warns and returns `NaN`, a fractional one
#'   returns weights for the truncated order scaled by `factorial(order)`.
#'
#' @return A numeric vector of weights for a unit step, one per offset and in
#'   the order the offsets were given. For `order >= 1` the entries sum to zero.
#'
#' @seealso [fd_offsets()] for the offsets to pass in, [fd_step()] for the step
#'   to divide by, and [fd_derivative()] for all three assembled into a
#'   derivative.
#'
#' @examples
#' # Three nodes give the classical central differences.
#' fd_weights(c(-1, 0, 1), 1)                 # -1/2, 0, 1/2
#' fd_weights(c(-1, 0, 1), 2)                 #    1,  -2,  1
#'
#' # Five nodes, first derivative: the familiar (1, -8, 0, 8, -1)/12.
#' fd_weights(-2:2, 1) * 12
#'
#' # The weights are for a unit step, so the caller divides by h^order.
#' h <- 1e-3
#' w <- fd_weights(c(-1, 0, 1), 2)
#' sum(w * exp(1 + c(-1, 0, 1) * h)) / h^2    # exp'' at 1
#' exp(1)
#'
#' # Exactness is what they are solved for. On five nodes a first derivative
#' # reproduces every polynomial up to degree four and fails at degree five.
#' s <- -2:2
#' w <- fd_weights(s, 1)
#' vapply(0:5, function(m) sum(w * s^m), numeric(1))
#'
#' # That first failure is the error constant: -4/5! = -1/30.
#' sum(w * s^5) / factorial(5)
#'
#' # Which predicts the error to three figures at a step of 0.05.
#' sum(w * sin(0.7 + s * 0.05)) / 0.05 - cos(0.7)
#' -1 / 30 * 0.05^4 * cos(0.7)
#'
#' # Offsets need not straddle the point; this one-sided stencil is what a
#' # fallback uses where a symmetric one would leave the domain.
#' fd_weights(0:2, 1)
#'
#' @export
fd_weights <- function(offsets, order) {
  n <- length(offsets)
  if (order >= n) {
    stop(sprintf(
      "A stencil on %d node(s) has no derivative of order %d to offer: the order must be smaller than the number of nodes.",
      n, order
    ), call. = FALSE)
  }
  if (anyDuplicated(offsets)) {
    stop("The offsets must be distinct.", call. = FALSE)
  }
  # Row i is the offsets raised to the power i - 1: the exponent indexes the
  # rows, the offsets the columns.
  A <- outer(seq_len(n) - 1L, offsets, function(power, s) s^power)
  rhs <- numeric(n)
  rhs[order + 1L] <- factorial(order)
  solve(A, rhs)
}


#' Stencil Offsets for a Derivative Order
#'
#' @description
#' The symmetric offsets used away from a boundary, and the one-sided ones
#' used where a symmetric stencil would not fit, sized from the derivative
#' order and the accuracy asked of it.
#'
#' @details
#' The reach is \eqn{r = \lceil (d + a)/2 \rceil - 1} for order \eqn{d} and
#' accuracy \eqn{a}, giving \eqn{2r + 1} nodes: at the default accuracy of
#' two this is the three-point stencil for the first and second derivatives
#' and the five-point one for the third and fourth, and at accuracy four it
#' is the five-point stencils for the first and second -- every stencil the
#' toolkit's packages had written out by hand, from one formula.
#'
#' @param order The derivative order.
#' @param accuracy The order of the error term, a positive integer. Central
#'   stencils gain a free order on symmetry, so an odd request costs the same
#'   nodes as the even one above it.
#'
#' @return A list with the \code{reach} and the three offset vectors
#'   \code{central}, \code{forward} and \code{backward}.
#'
#' @seealso \code{\link{fd_weights}}, \code{\link{fd_derivative}}
#'
#' @examples
#' fd_offsets(1)               # three points
#' fd_offsets(1, accuracy = 4) # five points
#' fd_offsets(4)$central       # the five-point fourth-derivative stencil
#'
#' @export
fd_offsets <- function(order, accuracy = 2L) {
  if (accuracy < 1L) stop("'accuracy' must be a positive integer.", call. = FALSE)
  r <- max(1L, as.integer(ceiling((order + accuracy) / 2)) - 1L)
  list(
    reach = r,
    central = seq.int(-r, r),
    forward = seq.int(0L, 2L * r),
    backward = seq.int(-2L * r, 0L)
  )
}


#' A Step Size for One Stencil
#'
#' @description
#' The magnitude-scaled step that balances truncation against rounding for a
#' single stencil of the given order and accuracy, clamped so the whole
#' stencil stays inside the domain when bounds are given.
#'
#' @details
#' The step is \eqn{\varepsilon^{1/(d+a)}\max(1, |x|)}: truncation grows like
#' \eqn{h^{a}} and rounding like \eqn{\varepsilon h^{-d}}, and this is where
#' they balance. Near a finite bound the step is shrunk so that the farthest
#' node of the stencil stays strictly inside, since a node outside the domain
#' does not make a derivative inaccurate, it makes it \code{NaN}.
#'
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order.
#' @param accuracy The accuracy the stencil will be built at.
#' @param bounds An optional length-two numeric vector of domain bounds.
#'
#' @return A numeric vector of steps, the same length as \code{x}.
#'
#' @seealso \code{\link{fd_derivative}}
#'
#' @examples
#' fd_step(c(0.5, 1000), 2)
#' # near a boundary the stencil is kept inside
#' fd_step(0.01, 2, bounds = c(0, Inf))
#'
#' @export
fd_step <- function(x, order, accuracy = 2L, bounds = NULL) {
  reach <- fd_offsets(order, accuracy)$reach
  h <- .Machine$double.eps^(1 / (order + accuracy)) * pmax(1, abs(x))
  if (!is.null(bounds)) {
    if (is.finite(bounds[1])) h <- pmin(h, 0.49 * (x - bounds[1]) / reach)
    if (is.finite(bounds[2])) h <- pmin(h, 0.49 * (bounds[2] - x) / reach)
  }
  h
}


#' One Stencil, Applied
#'
#' @description
#' The \code{order}-th derivative of \code{f} at \code{x} from a single
#' finite-difference stencil: the weighted sum of function values at the
#' stencil's nodes, divided by \eqn{h^{d}}.
#'
#' @details
#' This is the applicator every fallback in the toolkit speaks through, and
#' it enforces the one rule they share: one stencil of the requested order,
#' never a composition of lower-order differences. What it deliberately does
#' \strong{not} choose is the policy around it -- which order to fall back
#' from, when a reference can be trusted, what to do at a domain boundary
#' beyond keeping the nodes inside. Those belong to the callers, who know
#' what they are differentiating.
#'
#' \code{f} must be vectorized in its argument; \code{x} and \code{h} may be
#' vectors, and the stencil is applied elementwise.
#'
#' @param f A vectorized function of one numeric argument.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order, 1 to 4.
#' @param h A numeric vector of steps; \code{\link{fd_step}} when missing.
#' @param accuracy The order of the error term. Two gives the classical
#'   compact stencils; four gives the five-point first and second
#'   derivatives.
#' @param side \code{"central"} away from boundaries, \code{"forward"} or
#'   \code{"backward"} where a symmetric stencil would leave the domain.
#'
#' @return A numeric vector of the same length as \code{x}.
#'
#' @seealso \code{\link{fd_weights}}, \code{\link{fd_step}}
#'
#' @examples
#' # the third derivative of exp is exp
#' fd_derivative(exp, 1, order = 3)
#' exp(1)
#'
#' # a five-point first derivative is fourth-order accurate
#' fd_derivative(sin, 0.7, order = 1, accuracy = 4) - cos(0.7)
#'
#' # one-sided at a boundary
#' fd_derivative(sqrt, 1e-4, order = 1, side = "forward",
#'               h = fd_step(1e-4, 1, bounds = c(0, Inf)))
#'
#' @export
fd_derivative <- function(f, x, order, h = NULL,
                          accuracy = 2L,
                          side = c("central", "forward", "backward")) {
  side <- match.arg(side)
  off <- fd_offsets(order, accuracy)
  s <- off[[side]]
  w <- fd_weights(s, order)
  if (is.null(h)) h <- fd_step(x, order, accuracy)
  acc <- 0
  for (k in seq_along(s)) {
    if (w[k] == 0) next
    acc <- acc + w[k] * f(x + s[k] * h)
  }
  acc / h^order
}
