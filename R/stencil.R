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
#' The weights that combine function values at the given offsets into an
#' approximation of the derivative of the given order, obtained by solving
#' the Vandermonde system that makes the stencil exact on polynomials.
#'
#' @details
#' A stencil on \eqn{n} nodes is exact for polynomials up to degree
#' \eqn{n - 1}, so its error is \eqn{O(h^{n - d})} for the \eqn{d}-th
#' derivative -- fourth-order accurate for a first derivative on five nodes,
#' second-order for a fourth derivative on the same five.
#'
#' Building the weights this way, rather than composing lower-order
#' differences, is what keeps a high order usable: each numerical
#' differentiation multiplies the error of the one before it, so a fourth
#' derivative reached by four nested first differences is noise. One stencil,
#' never nested.
#'
#' @param offsets A numeric vector of distinct stencil offsets, in units of
#'   the step.
#' @param order The derivative order, smaller than \code{length(offsets)}.
#'
#' @return A numeric vector of weights, the same length as \code{offsets}.
#'
#' @seealso \code{\link{fd_offsets}}, \code{\link{fd_derivative}}
#'
#' @examples
#' fd_weights(c(-1, 0, 1), 1)             # the central first difference
#' fd_weights(c(-1, 0, 1), 2)             # 1, -2, 1
#' fd_weights(-2:2, 1) * 12               # the five-point first derivative
#' fd_weights(0:2, 1)                     # one-sided, for a boundary
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
