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
#' @param order The derivative order: a single non-negative whole number,
#'   **strictly smaller than `length(offsets)`**. Order `0` is legal and returns
#'   interpolation weights at the origin. Anything else throws: an order at or
#'   above the node count with a message naming both numbers, and a negative,
#'   fractional, missing or non-scalar order with a message naming the
#'   requirement.
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
  # Checked before the count, because a non-finite order makes the comparison
  # below NA and the error would come from the `if` rather than from here. A
  # fractional order used to fall through and return the weights of the
  # truncated order scaled by factorial(order): a plausible-looking stencil
  # that solves no moment condition.
  if (length(order) != 1L || !is.finite(order) ||
      order < 0 || order != trunc(order)) {
    stop("'order' must be a single non-negative whole number.", call. = FALSE)
  }
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
#' Sizes a stencil from the derivative order and the accuracy asked of it, and
#' returns the offsets to evaluate at: the symmetric ones used away from a
#' boundary, and the one-sided ones used where a symmetric stencil would leave
#' the domain. Pass them to [fd_weights()] for the weights and to
#' [fd_derivative()] to apply the whole thing.
#'
#' @details
#' # The reach
#'
#' For order \eqn{d} and accuracy \eqn{a} the half-width is
#'
#' \deqn{r = \Bigl\lceil \tfrac{d + a}{2} \Bigr\rceil - 1,}
#'
#' floored at one, giving \eqn{2r + 1} nodes. At the default accuracy of two
#' that is the three-point stencil for the first and second derivatives and the
#' five-point one for the third and fourth. At accuracy four it is the
#' five-point stencils for the first and second. Those are every stencil the
#' toolkit's packages had written out by hand, from one formula.
#'
#' # An odd accuracy is rounded, and which way depends on the order
#'
#' A central stencil is symmetric, so the odd powers cancel from its error
#' expansion and the accuracy it delivers is always even. An odd request is
#' therefore served by an even neighbor, and the parity of \eqn{d + a} decides
#' which one. Measured on \eqn{\exp} by halving the step:
#'
#' \tabular{lrrrr}{
#'   \strong{order}    \tab 1 \tab 2 \tab 3 \tab 4 \cr
#'   accuracy 2 \tab 2 \tab 2 \tab 2 \tab 2 \cr
#'   accuracy 3 \tab 2 \tab 4 \tab 2 \tab 4 \cr
#'   accuracy 4 \tab 4 \tab 4 \tab 4 \tab 4
#' }
#'
#' At an odd order the request rounds down and costs nothing extra; at an even
#' order it rounds up and buys two more nodes. Ask for an even accuracy and the
#' question does not arise.
#'
#' @param order The derivative order \eqn{d}. Not validated here, though
#'   [fd_weights()] rejects anything but a non-negative whole number when the
#'   offsets reach it.
#' @param accuracy The order of the error term, a positive integer, `2` by
#'   default. Zero or below throws. See above for what an odd value does.
#'
#' @return A list of four components:
#'   \describe{
#'     \item{`reach`}{integer, the half-width \eqn{r}, at least 1.}
#'     \item{`central`}{integer vector `-r:r`, the symmetric stencil.}
#'     \item{`forward`}{integer vector `0:(2r)`, for the lower boundary.}
#'     \item{`backward`}{integer vector `(-2r):0`, for the upper one.}
#'   }
#'   All three offset vectors have the same length, \eqn{2r + 1}, so the three
#'   sides cost the same number of evaluations.
#'
#' @seealso [fd_weights()] for the weights at these offsets, [fd_step()] for the
#'   step to pair with them, [fd_derivative()] for all three assembled.
#'
#' @examples
#' # Three points for a first or second derivative, five for a third or fourth.
#' fd_offsets(1)$central
#' fd_offsets(4)$central
#'
#' # Accuracy four buys two more nodes for a first derivative.
#' fd_offsets(1, accuracy = 4)$central
#'
#' # The one-sided sets are the same size, so a boundary costs no more.
#' str(fd_offsets(2))
#'
#' # Accuracy must be positive.
#' try(fd_offsets(2, accuracy = 0))
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
#' Returns the step at which a single stencil of the given order and accuracy
#' balances truncation error against rounding error, scaled by the magnitude of
#' the point, and shrunk near a finite bound so that the whole stencil stays
#' inside the domain.
#'
#' @details
#' # Where the balance falls
#'
#' \deqn{h = \varepsilon^{1/(d + a)} \max(1, \lvert x \rvert).}
#'
#' Truncation grows like \eqn{h^{a}} and rounding like \eqn{\varepsilon h^{-d}},
#' and this is where the two meet. The factor \eqn{\max(1, \lvert x \rvert)}
#' makes the step relative for a large argument and absolute for a small one, so
#' a point near zero does not get a step below the resolution of its own
#' neighborhood.
#'
#' At the default accuracy the exponent is \eqn{1/4} for a second derivative and
#' \eqn{1/6} for a fourth, giving steps of about `1.2e-4` and `2.5e-3` at
#' \eqn{x = 1}. A high order wants a *large* step, since rounding is what
#' dominates there.
#'
#' # Staying inside the domain
#'
#' Given `bounds`, the step is shrunk so the farthest node of the stencil stays
#' strictly inside: a node outside the domain does not make a derivative
#' inaccurate, it makes it `NaN`. The margin is 0.49 of the distance to the
#' bound, divided by the reach.
#'
#' # The one case to guard
#'
#' A point sitting exactly on a finite bound gets a step of **zero**, and a
#' stencil divided by \eqn{h^{d}} is then `NaN`. The evaluation point is the
#' caller's, so this is not checked here; either keep the point off the bound or
#' use a one-sided stencil with a step of your own.
#'
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order \eqn{d}.
#' @param accuracy The accuracy the stencil will be built at, `2` by default.
#'   Matches [fd_offsets()]'s argument of the same name.
#' @param bounds An optional numeric vector of two domain bounds, either of
#'   which may be infinite. `NULL`, the default, applies no clamp.
#'
#' @return A numeric vector of steps the same length as `x`, positive except at
#'   a point sitting on a finite bound, where it is zero.
#'
#' @seealso [fd_derivative()], which calls this when no step is given,
#'   [fd_offsets()] for the reach the clamp divides by.
#'
#' @examples
#' # Magnitude-scaled: absolute near zero, relative far from it.
#' fd_step(c(0.5, 1000), 2)
#'
#' # A higher order wants a larger step, rounding being what dominates.
#' c(order2 = fd_step(1, 2), order4 = fd_step(1, 4))
#'
#' # Near a boundary the stencil is kept inside the domain.
#' fd_step(0.01, 2, bounds = c(0, Inf))
#'
#' # On the boundary the step is zero, which no stencil can use.
#' fd_step(0, 2, bounds = c(0, Inf))
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
#' Estimates the `order`-th derivative of `f` at `x` from a single
#' finite-difference stencil: the weighted sum of function values at the
#' stencil's nodes, divided by \eqn{h^{d}}. It picks the offsets with
#' [fd_offsets()], the weights with [fd_weights()] and, unless given one, the
#' step with [fd_step()].
#'
#' @details
#' # One stencil, never nested
#'
#' This is the applicator every numerical fallback in the toolkit speaks
#' through, and it enforces the rule they share: one stencil of the order
#' requested, never a composition of lower-order differences. Each numerical
#' differentiation multiplies the error of the one before it, so a fourth
#' derivative reached by four nested first differences is noise.
#'
#' # What it deliberately leaves to the caller
#'
#' The policy around the stencil. Which order to fall back from, when a
#' reference can be trusted, and what to do at a domain boundary beyond keeping
#' the nodes inside are all decisions that need to know what is being
#' differentiated, and this function does not.
#'
#' # Vectorization
#'
#' `f` must be vectorized in its argument. `x` and `h` may be vectors, and the
#' stencil is applied elementwise, so a whole vector of points costs
#' \eqn{2r + 1} calls to `f` and no more.
#'
#' @param f A vectorized function of one numeric argument.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order. One to four is what the toolkit uses and
#'   what the step rule is tuned for. A higher order is accepted and sized by
#'   [fd_offsets()], but nothing is claimed about its accuracy.
#' @param h A numeric vector of steps, recycled against `x`. `NULL`, the
#'   default, takes [fd_step()] at the same order and accuracy, with no bounds;
#'   pass a step explicitly when the domain is bounded.
#' @param accuracy The order of the error term, `2` by default. Two gives the
#'   classical compact stencils, four the five-point first and second
#'   derivatives.
#' @param side `"central"` away from boundaries, `"forward"` or `"backward"`
#'   where a symmetric stencil would leave the domain. All three use
#'   \eqn{2r + 1} nodes, so a one-sided estimate costs the same and is one order
#'   less accurate.
#'
#' @return A numeric vector of the same length as `x`.
#'
#' @seealso [fd_weights()], [fd_offsets()] and [fd_step()], the three pieces
#'   this assembles.
#'
#' @examples
#' # The third derivative of exp is exp.
#' fd_derivative(exp, 1, order = 3) - exp(1)
#'
#' # Accuracy four is worth about two decades here.
#' c(acc2 = fd_derivative(sin, 0.7, 1) - cos(0.7),
#'   acc4 = fd_derivative(sin, 0.7, 1, accuracy = 4) - cos(0.7))
#'
#' # At a boundary, one-sided with a step that keeps the nodes inside. The
#' # central stencil would reach below zero, where sqrt is not defined.
#' fd_derivative(sqrt, 1e-4, order = 1, side = "forward",
#'               h = fd_step(1e-4, 1, bounds = c(0, Inf)))
#' 0.5 / sqrt(1e-4)
#'
#' # Vectorized in the point: one call per node, not one per point.
#' fd_derivative(sin, c(0, pi / 4, pi / 2), order = 1) - cos(c(0, pi / 4, pi / 2))
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
