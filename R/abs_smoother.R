#' @title S7 Class for Smoothers of the Absolute Value
#'
#' @description
#' A smooth replacement \eqn{s(u)} for \eqn{\lvert u \rvert}, carrying its
#' derivatives in \eqn{u} up to order five as functions. Every non-smooth
#' primitive the toolkit uses is generated from the absolute value:
#'
#' \deqn{\operatorname{sign}(u) = \frac{d\lvert u\rvert}{du}, \qquad
#'   \mathbb{1}(u \geq 0) = \frac{1 + \operatorname{sign}(u)}{2}, \qquad
#'   (u)_{+} = \frac{u + \lvert u\rvert}{2},}
#'
#' so one contract serves them all. The smooth sign is \eqn{s'}, the smooth
#' step \eqn{(1 + s'(u))/2} and the smooth hinge \eqn{(u + s(u))/2} follow by
#' composition, and a consumer that replaces \eqn{\lvert u\rvert} by \eqn{s(u)}
#' has replaced every one of them consistently.
#'
#' @details
#' # What a smoother buys
#'
#' A break-point term smoothed this way becomes an ordinary nonlinear term
#' whose design block is the true Jacobian, so a random or penalized
#' development of its break-points becomes fittable. A kinked penalty
#' smoothed this way becomes a proper separable one, at the price of the exact
#' zeros the kink produced.
#'
#' # The derivatives are functions
#'
#' They are functions and not expressions. A piecewise smoother, the quintic
#' here, has branches that [stats::deriv()] does not read, and with the
#' derivatives written out in the constructor those branches are ordinary code.
#' Each takes `(u, width)` and vectorizes in both, so a per-group width is one
#' value per observation.
#'
#' # The width
#'
#' `width` is the transition scale: `h` for a smoother whose parameter is a
#' length, the bent-cable reading of a transition of width \eqn{h}, and `c` for
#' the hyperbolic, whose parameter is a squared length. `NULL`, the default of
#' all three constructors, asks the consumer to resolve it from the data at
#' build through [smoother_width()]; a break-point term hands it the median
#' spacing of its covariate. `per_group` asks for one width per group where a
#' grouping is available, the validity window of a Laplace approximation being
#' per-subject.
#'
#' # The scale correction
#'
#' `tau_correction` is a property of the mollifier and of no model. The probit
#' smoother satisfies an exact convolution identity, smoothing with width
#' \eqn{h} being the same as convolving the break-point with \eqn{N(0, h^2)},
#' so an apparent scale \eqn{\tau} of a random break-point composes as
#' \eqn{\tau^2_{\mathrm{true}} = \tau^2 - h^2} and the smoother declares that
#' correction. Smoothers with no such identity declare `NULL`, and a consumer
#' then reports the apparent scale alone.
#'
#' @param smoother_name A single non-empty string naming the smoother, used by
#'   `print()` and by a consumer reporting which smoother a term carries.
#' @param width The transition width, a single positive number, or `NULL` to be
#'   resolved at build.
#' @param width_name What the width parameter is called, a single non-empty
#'   string: `"h"` for a length and `"c"` for a squared one. `"h"` by default.
#' @param per_group `TRUE` when the width is resolved per group, `FALSE` (the
#'   default) for one width throughout.
#' @param s A list of exactly six functions of `(u, width)`: \eqn{s} and its
#'   derivatives in \eqn{u} of orders one to five, in that order. Each must
#'   take at least two arguments; the validator checks the count and the arity.
#' @param width_from_spacing `NULL` (the identity) or a function carrying a
#'   spacing in covariate units onto the width parameter's own scale, which is
#'   the square for the hyperbolic.
#' @param tau_correction `NULL`, or a function `(tau, width)` returning the
#'   corrected scale of a random break-point.
#' @param exact_radius `NULL`, or a function of the width returning the radius
#'   beyond which \eqn{s(u) = \lvert u\rvert} exactly. Only the quintic has
#'   one.
#'
#' @return An S7 object of class `abs_smoother` carrying the eight properties
#'   above.
#'
#' @references
#' Bacon, D. W. and Watts, D. G. (1971). Estimating the transition between
#' two intersecting straight lines. *Biometrika*, **58**(3), 525--534.
#'
#' Tishler, A. and Zang, I. (1981). A new maximum likelihood algorithm for
#' piecewise regression. *Journal of the American Statistical
#' Association*, **76**(376), 980--987.
#'
#' Seo, M. H. and Linton, O. (2007). A smoothed least squares estimator
#' for threshold regression models. *Journal of Econometrics*,
#' **141**(2), 704--735.
#'
#' @seealso [smooth_probit()], [smooth_hyperbolic()] and [smooth_quintic()]
#'   for the three that ship, [smoother_deriv()] to evaluate one,
#'   [smoother_width()] and [smoother_width_floor()] for the width,
#'   [check_abs_smoother()] to verify one of your own.
#'
#' @examples
#' sm <- smooth_probit(h = 0.3)
#' S7::S7_inherits(sm, abs_smoother)
#' sm
#'
#' # The smooth sign, step and hinge all follow from the same object.
#' u <- c(-1, -0.1, 0, 0.1, 1)
#' smoother_deriv(sm, u, order = 1)                    # smooth sign
#' (1 + smoother_deriv(sm, u, order = 1)) / 2          # smooth step
#' (u + smoother_deriv(sm, u, order = 0)) / 2          # smooth hinge
#'
#' # Three widths out they agree with the sharp versions to five decimals,
#' # and inside the transition they are what replaces them.
#' far <- c(-1.2, -0.9, 0.9, 1.2)
#' round(smoother_deriv(sm, far, order = 1) - sign(far), 5)
#' round((far + smoother_deriv(sm, far, order = 0)) / 2 - pmax(far, 0), 5)
#'
#' @export
abs_smoother <- S7::new_class(
  name = "abs_smoother",
  properties = list(
    smoother_name = S7::class_character,
    width = S7::class_any,
    width_name = S7::new_property(S7::class_character, default = "h"),
    per_group = S7::new_property(S7::class_logical, default = FALSE),
    s = S7::class_list,
    width_from_spacing = S7::class_any,
    tau_correction = S7::class_any,
    exact_radius = S7::class_any
  ),
  validator = function(self) {
    if (length(self@smoother_name) != 1L || is.na(self@smoother_name) ||
        !nzchar(self@smoother_name)) {
      return("@smoother_name must be a single non-empty string")
    }
    if (!is.null(self@width) &&
        (!is.numeric(self@width) || length(self@width) != 1L ||
         !is.finite(self@width) || self@width <= 0)) {
      return("@width must be NULL or a single positive number")
    }
    if (length(self@width_name) != 1L || !nzchar(self@width_name)) {
      return("@width_name must be a single non-empty string")
    }
    if (length(self@per_group) != 1L || is.na(self@per_group)) {
      return("@per_group must be TRUE or FALSE")
    }
    if (length(self@s) != 6L ||
        !all(vapply(self@s, is.function, logical(1)))) {
      return(paste("@s must be a list of six functions of (u, width): s and",
                   "its derivatives of orders one to five"))
    }
    if (any(vapply(self@s, function(f) length(formals(f)) < 2L,
                   logical(1)))) {
      return("every function in @s must take (u, width)")
    }
    for (nm in c("width_from_spacing", "tau_correction", "exact_radius")) {
      v <- S7::prop(self, nm)
      if (!is.null(v) && !is.function(v)) {
        return(sprintf("@%s must be NULL or a function", nm))
      }
    }
    NULL
  }
)

#' The Probit Smoother of the Absolute Value
#'
#' @description
#' \eqn{s(u) = \mathbb{E}\lvert u + hZ\rvert} for \eqn{Z} standard normal, the
#' one to reach for unless you have a reason not to. Its excess over
#' \eqn{\lvert u\rvert} has gaussian tails, so the smoothing bias is confined
#' to a window of width \eqn{h} around the kink, and it is the only one of the
#' three that declares a scale correction for a random break-point.
#'
#' @details
#' # The closed form
#'
#' \deqn{s(u) = u\,\bigl(2\Phi(u/h) - 1\bigr) + 2h\,\phi(u/h), \qquad
#'   s'(u) = 2\Phi(u/h) - 1, \qquad s''(u) = \frac{2}{h}\,\phi(u/h),}
#'
#' with the third to fifth derivatives following from
#' \eqn{\phi'(z) = -z\phi(z)}. At the kink \eqn{s(0) = 2h\phi(0) \approx
#' 0.798h}, which is the intrinsic scale [check_abs_smoother()] places its grid
#' by.
#'
#' # Why the tails matter
#'
#' The excess \eqn{s(u) - \lvert u\rvert} decays like \eqn{\phi(u/h)}, so it
#' is gone within a few widths. Measured at \eqn{h = 0.3}: `6.7e-05` at
#' \eqn{u = 1} and exactly `0` at \eqn{u = 3}, against `4.4e-02` and `1.5e-02`
#' for [smooth_hyperbolic()] at the same transition width. A model smoothed
#' this way is the sharp model to within rounding a short distance from the
#' break-point.
#'
#' # The convolution identity
#'
#' Smoothing the step with width \eqn{h} is exactly convolving the break-point
#' with \eqn{N(0, h^2)}, so a random break-point of true scale \eqn{\tau}
#' appears at \eqn{\sqrt{\tau^2 + h^2}}, and the smoother declares
#' \eqn{\tau_{\mathrm{true}} = \sqrt{\tau^2 - h^2}}, floored at zero. No other
#' smoother here has such an identity.
#'
#' @param h The transition width, a single positive number, or `NULL` (the
#'   default) to be resolved at build from the covariate's spacing through
#'   [smoother_width()].
#' @param per_group `TRUE` to resolve the width once per group, `FALSE` (the
#'   default) for one width throughout. The validity window of a Laplace
#'   approximation is per-subject, so a hierarchical break-point wants `TRUE`.
#'
#' @return An [abs_smoother()] named `"probit"`, with `width_name` `"h"`, no
#'   `width_from_spacing` (the spacing is the width), a `tau_correction`, and
#'   no `exact_radius`.
#'
#' @examples
#' sm <- smooth_probit(h = 0.3)
#'
#' # At the kink the value is 2 h phi(0), not zero.
#' smoother_deriv(sm, 0, order = 0)
#' 2 * 0.3 * stats::dnorm(0)
#'
#' # And the excess over |u| is gone within a few widths.
#' u <- c(0, 0.3, 1, 3)
#' smoother_deriv(sm, u, order = 0) - abs(u)
#'
#' # The scale correction it declares, applied to an apparent tau.
#' sm@tau_correction(0.5, 0.231)
#' sqrt(0.5^2 - 0.231^2)
#'
#' @seealso [abs_smoother()] for the contract, [smooth_hyperbolic()] and
#'   [smooth_quintic()] for the alternatives, [smoother_width()] for
#'   resolving `h` from data, [check_abs_smoother()] to verify it.
#' @export
smooth_probit <- function(h = NULL, per_group = FALSE) {
  abs_smoother(
    smoother_name = "probit", width = h, width_name = "h",
    per_group = per_group,
    s = list(
      function(u, width) {
        z <- u / width
        u * (2 * stats::pnorm(z) - 1) + 2 * width * stats::dnorm(z)
      },
      function(u, width) 2 * stats::pnorm(u / width) - 1,
      function(u, width) 2 * stats::dnorm(u / width) / width,
      function(u, width) {
        z <- u / width
        -2 * z * stats::dnorm(z) / width^2
      },
      function(u, width) {
        z <- u / width
        2 * (z^2 - 1) * stats::dnorm(z) / width^3
      },
      function(u, width) {
        z <- u / width
        2 * (3 * z - z^3) * stats::dnorm(z) / width^4
      }
    ),
    width_from_spacing = NULL,
    tau_correction = function(tau, width) sqrt(pmax(tau^2 - width^2, 0)),
    exact_radius = NULL
  )
}

#' The Hyperbolic Smoother of the Absolute Value
#'
#' @description
#' \eqn{s(u) = \sqrt{u^2 + c}}, the simplest smoother of the three and the one
#' with the worst tails. Its parameter \eqn{c} is a **squared** length, so the
#' transition width in covariate units is \eqn{\sqrt{c}}.
#'
#' @details
#' # The closed form
#'
#' \deqn{s(u) = \sqrt{u^2 + c}, \qquad s'(u) = \frac{u}{\sqrt{u^2 + c}},
#'   \qquad s''(u) = \frac{c}{(u^2 + c)^{3/2}},}
#'
#' every order a rational function of \eqn{u} over a half-integer power, so
#' nothing here needs a special function or a branch.
#'
#' # The cost of the tails
#'
#' The excess over \eqn{\lvert u\rvert} decays only as \eqn{c/(2\lvert
#' u\rvert)}, so the smoothing bias spreads well away from the kink. Measured
#' at a transition width of `0.3`, so \eqn{c = 0.09}: the excess is `4.4e-02`
#' at \eqn{u = 1} and `1.5e-02` at \eqn{u = 3}, against `6.7e-05` and exactly
#' `0` for [smooth_probit()] at the same width. A fit smoothed this way is
#' perturbed everywhere, not only near the break-point.
#'
#' No convolution identity relates \eqn{c} to the scale of a random
#' break-point, so `tau_correction` is `NULL` and a consumer reports the
#' apparent scale alone.
#'
#' @param c The squared transition width, a single positive number, or `NULL`
#'   (the default) to be resolved at build as the **square** of the covariate's
#'   spacing, through [smoother_width()].
#'
#' @return An [abs_smoother()] named `"hyperbolic"`, with `width_name` `"c"`, a
#'   `width_from_spacing` that squares, and neither a `tau_correction` nor an
#'   `exact_radius`.
#'
#' @examples
#' # c is a squared length, so a transition of 0.3 is c = 0.09.
#' sm <- smooth_hyperbolic(c = 0.09)
#' smoother_deriv(sm, 0, order = 0)
#' sqrt(0.09)
#'
#' # And resolving from a spacing squares it for you.
#' smoother_width(smooth_hyperbolic(), 0.3)
#'
#' # The excess over |u| is still there three units out.
#' u <- c(0, 0.3, 1, 3)
#' smoother_deriv(sm, u, order = 0) - abs(u)
#' smoother_deriv(smooth_probit(h = 0.3), u, order = 0) - abs(u)
#'
#' @seealso [abs_smoother()] for the contract, [smooth_probit()] for the
#'   recommended default, [smooth_quintic()] for exactness outside the
#'   transition, [check_abs_smoother()] to verify it.
#' @export
smooth_hyperbolic <- function(c = NULL) {
  abs_smoother(
    smoother_name = "hyperbolic", width = c, width_name = "c",
    per_group = FALSE,
    s = list(
      function(u, width) sqrt(u^2 + width),
      function(u, width) u / sqrt(u^2 + width),
      function(u, width) width / (u^2 + width)^1.5,
      function(u, width) -3 * width * u / (u^2 + width)^2.5,
      function(u, width) 3 * width * (4 * u^2 - width) / (u^2 + width)^3.5,
      function(u, width) {
        15 * width * u * (3 * width - 4 * u^2) / (u^2 + width)^4.5
      }
    ),
    width_from_spacing = function(gap) gap^2,
    tau_correction = NULL,
    exact_radius = NULL
  )
}

#' The Quintic Smoother of the Absolute Value
#'
#' @description
#' \eqn{\lvert u\rvert} replaced inside \eqn{[-h, h]} by the even polynomial
#' matching it to third order at the seam, and **exact** outside. It is the
#' only smoother here whose bias is zero beyond a finite radius, and the only
#' one that is not analytic.
#'
#' @details
#' # The closed form
#'
#' With \eqn{t = u/h},
#'
#' \deqn{s(u) = \frac{h}{16}\,\bigl(5 + 15t^2 - 5t^4 + t^6\bigr)
#'   \quad \text{for } \lvert t\rvert < 1, \qquad
#'   s(u) = \lvert u\rvert \quad \text{otherwise},}
#'
#' so \eqn{s''(u) = 15(1 - t^2)^2/(8h)} vanishes to second order at \eqn{\pm h}
#' and the function is \eqn{C^3}, with a jump in the fourth derivative at the
#' seam. At the kink \eqn{s(0) = 5h/16}.
#'
#' # Exact outside the transition
#'
#' The smoothing bias is exactly zero outside \eqn{[-h, h]}, which is the
#' cleanest fixed-width choice: measured at \eqn{h = 0.3} the excess over
#' \eqn{\lvert u\rvert} is `0` at \eqn{u = 1} and at \eqn{u = 3}, where
#' [smooth_probit()] leaves `6.7e-05` and [smooth_hyperbolic()] `4.4e-02`. The
#' smoother declares an `exact_radius`, so a consumer can say where its answer
#' is the sharp model's.
#'
#' The branches are why [abs_smoother()] carries the derivatives as functions:
#' [stats::deriv()] does not read a clamp, and written out in the constructor
#' the branches are ordinary code.
#'
#' There is no convolution identity here, so `tau_correction` is `NULL`.
#'
#' @param h The transition half-width, a single positive number, or `NULL` (the
#'   default) to be resolved at build from the covariate's spacing through
#'   [smoother_width()].
#'
#' @return An [abs_smoother()] named `"quintic"`, with `width_name` `"h"`, no
#'   `width_from_spacing`, no `tau_correction`, and an `exact_radius` equal to
#'   the width.
#'
#' @examples
#' sm <- smooth_quintic(h = 0.5)
#'
#' # Outside the transition it is the absolute value, exactly.
#' smoother_deriv(sm, 1, order = 0)
#' smoother_deriv(sm, c(0.6, 1, 4), order = 0) - c(0.6, 1, 4)
#'
#' # At the kink the value is 5h/16.
#' smoother_deriv(sm, 0, order = 0)
#' 5 * 0.5 / 16
#'
#' # The second derivative vanishes to second order at the seam, which is
#' # what makes it C^3 there.
#' smoother_deriv(sm, c(0.4, 0.49, 0.5, 0.51), order = 2)
#'
#' @seealso [abs_smoother()] for the contract, [smooth_probit()] for the
#'   recommended default, [smooth_hyperbolic()] for the simplest form,
#'   [check_abs_smoother()] to verify it.
#' @export
smooth_quintic <- function(h = NULL) {
  inside <- function(u, width) abs(u) < width
  abs_smoother(
    smoother_name = "quintic", width = h, width_name = "h",
    per_group = FALSE,
    s = list(
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width),
               width * (5 + 15 * t^2 - 5 * t^4 + t^6) / 16, abs(u))
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), (15 * t - 10 * t^3 + 3 * t^5) / 8, sign(u))
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), 15 * (1 - t^2)^2 / (8 * width), 0)
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), -15 * t * (1 - t^2) / (2 * width^2), 0)
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), 15 * (3 * t^2 - 1) / (2 * width^3), 0)
      },
      function(u, width) {
        t <- u / width
        ifelse(inside(u, width), 45 * t / width^4, 0)
      }
    ),
    width_from_spacing = NULL,
    tau_correction = NULL,
    exact_radius = function(width) width
  )
}

#' Evaluate a Smoother or One of Its Derivatives
#'
#' @description
#' Returns \eqn{d^k s/du^k} at the given points, for \eqn{k} from 0 to 5. The
#' width is the smoother's own unless one is supplied, which is how a consumer
#' that resolved the width at build evaluates without mutating the object.
#'
#' @details
#' A vector `width` is one value per point, which is how a per-group width
#' arrives: a break-point term with `per_group = TRUE` resolves one width per
#' subject and passes the whole vector.
#'
#' The derivatives are the ones the smoother carries, evaluated directly. None
#' is differenced, so order five is as accurate as order zero.
#'
#' @param smoother An [abs_smoother()]. Anything else is rejected.
#' @param u A numeric vector of points.
#' @param width The width, a single positive number or one per point, or `NULL`
#'   (the default) for the smoother's own. A smoother built with `width = NULL`
#'   and given none here is rejected, with the two ways to supply one named.
#' @param order The derivative order, a whole number from 0 to 5. `0L` by
#'   default, which is \eqn{s} itself. Anything outside that range is rejected.
#'
#' @return A numeric vector as long as `u`.
#'
#' @examples
#' # The smooth sign: odd, bounded by one, and passing through zero.
#' smoother_deriv(smooth_hyperbolic(c = 1), 0:3, order = 1)
#'
#' # A per-point width, as a per-group smoother supplies.
#' smoother_deriv(smooth_probit(), c(-1, 0, 1), width = c(0.1, 0.5, 1),
#'                order = 0)
#'
#' # An unresolved width is an error rather than a guess.
#' try(smoother_deriv(smooth_probit(), 1))
#'
#' @seealso [abs_smoother()] for the contract, [smoother_width()] for
#'   resolving a width, [check_abs_smoother()] for verifying the derivatives.
#' @export
smoother_deriv <- function(smoother, u, width = NULL, order = 0L) {
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 0L || order > 5L) {
    stop("'order' must be an integer from 0 to 5.", call. = FALSE)
  }
  if (is.null(width)) width <- smoother@width
  if (is.null(width)) {
    stop(paste("the smoother's width is unresolved; supply 'width' or",
               "construct the smoother with one."), call. = FALSE)
  }
  smoother@s[[order + 1L]](as.numeric(u), width)
}

#' Resolve a Smoother's Width from a Spacing
#'
#' @description
#' Returns the width the smoother carries where it carries one, and otherwise
#' the given spacing carried onto the width parameter's own scale. This is what
#' a consumer calls at build: a break-point term hands it the median spacing of
#' its covariate, the smallest transition the data can tell from a step.
#'
#' @details
#' A smoother constructed with an explicit width keeps it, and `spacing` is
#' then not even looked at. A smoother constructed with `NULL` takes the
#' spacing through its own `width_from_spacing`, which is the identity for
#' [smooth_probit()] and [smooth_quintic()], whose parameter is a length, and
#' the square for [smooth_hyperbolic()], whose parameter is a squared one.
#'
#' The result is nothing about whether the width is large enough for the
#' arithmetic; [smoother_width_floor()] answers that separately, and a consumer
#' takes the larger of the two.
#'
#' @param smoother An [abs_smoother()]. Anything else is rejected.
#' @param spacing A spacing in covariate units, one value or one per group.
#'   Every entry must be positive and none may be missing; a smoother that
#'   already carries a width never reaches the check.
#'
#' @return The width, on the smoother's own scale. The same length as `spacing`
#'   when resolved from it, and a single number when the smoother carries one.
#'
#' @examples
#' # A length-parametrized smoother takes the spacing as it stands.
#' smoother_width(smooth_probit(), 0.3)
#'
#' # The hyperbolic's parameter is a squared length, so it is squared.
#' smoother_width(smooth_hyperbolic(), 0.3)
#'
#' # A smoother that carries a width keeps it, whatever the spacing.
#' smoother_width(smooth_probit(h = 0.5), 0.3)
#'
#' # One width per group.
#' smoother_width(smooth_probit(per_group = TRUE), c(0.2, 0.4, 0.35))
#'
#' @seealso [smoother_width_floor()] for the lower bound the arithmetic
#'   imposes, [abs_smoother()] for what a width means.
#' @export
smoother_width <- function(smoother, spacing) {
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  if (!is.null(smoother@width)) return(smoother@width)
  spacing <- as.numeric(spacing)
  if (!length(spacing) || anyNA(spacing) || any(spacing <= 0)) {
    stop("'spacing' must be positive.", call. = FALSE)
  }
  f <- smoother@width_from_spacing
  if (is.null(f)) spacing else f(spacing)
}

#' The Smallest Width a Consumer May Use
#'
#' @description
#' Returns the floor of the width, derived from the expression that binds. A
#' consumer resolving a width from data takes the larger of this and
#' [smoother_width()]'s answer.
#'
#' @details
#' # Where the bound comes from
#'
#' The derivatives of a smoother scale as \eqn{s^{(k)} \sim h^{1-k}}, so the
#' Jacobian column of a smoothed break-point carries \eqn{s''(0)/2 \sim 1/h}
#' against covariate columns of order the range \eqn{D}. Holding the design's
#' condition number below \eqn{\epsilon^{-1/2}}, which leaves half the digits
#' of a double to a QR of the design, gives
#'
#' \deqn{h \geq \sqrt{\epsilon}\,D,}
#'
#' carried onto the width parameter's own scale for a smoother parametrized by
#' a squared length. At a covariate range of 10 that is `1.49e-07` for
#' [smooth_probit()] and its square, `2.22e-14`, for [smooth_hyperbolic()].
#'
#' The bound is derived, which is the toolkit's rule for a guard constant, and
#' it is the same argument the break-point schedule's own scaling floor rests
#' on.
#'
#' @param smoother An [abs_smoother()]. Anything else is rejected.
#' @param scale The scale of the covariate, its range, a single positive
#'   number. Anything else is rejected.
#'
#' @return A single number, on the smoother's own width scale.
#'
#' @examples
#' # A length-parametrized smoother: sqrt(eps) times the range.
#' smoother_width_floor(smooth_probit(), scale = 10)
#' sqrt(.Machine$double.eps) * 10
#'
#' # The hyperbolic's floor is the square of that, its parameter being one.
#' smoother_width_floor(smooth_hyperbolic(), scale = 10)
#'
#' # A consumer takes the larger of the data's spacing and the floor.
#' max(smoother_width(smooth_probit(), 0.3),
#'     smoother_width_floor(smooth_probit(), scale = 10))
#'
#' @seealso [smoother_width()] for the width the data suggests,
#'   [abs_smoother()] for what a width means.
#' @export
smoother_width_floor <- function(smoother, scale) {
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  scale <- as.numeric(scale)
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("'scale' must be a single positive number.", call. = FALSE)
  }
  h <- sqrt(.Machine$double.eps) * scale
  f <- smoother@width_from_spacing
  if (is.null(f)) h else f(h)
}

#' @title Print a Smoother of the Absolute Value
#' @name print.abs_smoother
#'
#' @description
#' Writes one line naming the smoother and its width, and adds a line for each
#' of the two optional properties it declares: a scale correction for a random
#' break-point, and exactness outside the transition.
#'
#' @details
#' The width line reads `h resolved at build` where the smoother carries none,
#' and names the width parameter by its own name, so a hyperbolic smoother
#' shows `c` where a probit shows `h`. `per group` is appended when the width
#' is to be resolved per group.
#'
#' @param x An [abs_smoother()] object.
#' @param ... Unused, accepted for consistency with [print()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' smooth_probit()
#' smooth_probit(h = 0.2, per_group = TRUE)
#' smooth_hyperbolic()
#' smooth_quintic(h = 0.4)
#'
#' @keywords internal
S7::method(print, abs_smoother) <- function(x, ...) {
  cat(sprintf("<abs_smoother> %s: %s%s\n", x@smoother_name,
              if (is.null(x@width)) {
                sprintf("%s resolved at build", x@width_name)
              } else {
                sprintf("%s = %g", x@width_name, x@width)
              },
              if (isTRUE(x@per_group)) ", per group" else ""))
  if (!is.null(x@tau_correction)) {
    cat("  declares a scale correction for a random break-point\n")
  }
  if (!is.null(x@exact_radius)) {
    cat("  exact outside the transition\n")
  }
  invisible(x)
}

#' Check a Smoother of the Absolute Value Numerically
#'
#' @description
#' Compares a smoother against the properties every smoother of
#' \eqn{\lvert u\rvert} must have, and each derivative order against one
#' numerical differentiation of the analytic order below it. Returns one row
#' per check with the worst error and a pass or fail. Write a smoother of your
#' own and this says whether its five derivatives are right.
#'
#' @details
#' # The rows
#'
#' Ten checks, and two more for a smoother declaring a scale correction:
#'
#' | row | what it asserts |
#' |---|---|
#' | `s is even` | \eqn{s(-u) = s(u)}, relative to the size of \eqn{s} |
#' | `s' is odd` | \eqn{s'(-u) = -s'(u)} |
#' | `\|s'\| bounded by one` | the smooth sign never leaves \eqn{[-1, 1]} |
#' | `s convex` | \eqn{s'' \ge 0} on the grid |
#' | `matches \|u\| in the tails` | the excess ten intrinsic widths out is at most a fifth of the excess at the kink |
#' | `order k vs numDeriv on order k-1` | five rows, \eqn{k} from 1 to 5 |
#' | `tau_correction bounded by tau` | the correction never increases a scale |
#' | `tau_correction is the identity at width zero` | it vanishes as the smoothing does |
#'
#' Order \eqn{k} goes against \pkg{numDeriv} applied **once** to the analytic
#' order \eqn{k-1}, so a wrong derivative is caught against an independent
#' route while the reference never degenerates into a difference of
#' differences. Over the three shipped smoothers at four widths each, every row
#' passes and the worst error is `1.3e-10`.
#'
#' The tail row is deliberately loose, at a fifth of the excess at the kink,
#' because [smooth_hyperbolic()]'s polynomial tails would fail a tight one and
#' are a legitimate choice: its excess is \eqn{c/(2\lvert u\rvert)}, a
#' twentieth of the kink's value ten widths out.
#'
#' # Where the grid goes
#'
#' The grid spans the transition and the tails, placed by the smoother's own
#' intrinsic scale \eqn{s(0)}, which is of order the transition width for any
#' smoother of the absolute value, so no reading of what the width parameter
#' means is needed. Every point sits off zero and off the seam: a piecewise
#' smoother has measure-zero points where a one-sided derivative is read, and a
#' difference straddling one compares nothing.
#'
#' # What it catches
#'
#' This is not decoration. Writing [smooth_quintic()], a factor-2 error in its
#' own \eqn{s'''} was caught here before anything shipped: the order-3 row read
#' a relative `0.5` and the order-4 row `1.0`, the wrong third derivative
#' becoming the reference for the fourth. The examples below reproduce it.
#'
#' @param smoother An [abs_smoother()]. Anything else is rejected.
#' @param width The width to check at, a single positive number, or `NULL` (the
#'   default) for the smoother's own. A smoother carrying none is checked at
#'   `0.5`, a scale chosen to exercise the formulas and saying nothing about
#'   any data.
#' @param tol The relative error above which a row is reported as `FAILED`. A
#'   single positive number, `1e-6` by default. The worst error the three
#'   shipped smoothers reach is `1.3e-10`, four orders under it, so the default
#'   separates a correct smoother from one whose derivative is wrong in its
#'   fourth digit.
#' @param verbose `TRUE`, the default, prints a header naming the smoother and
#'   the width, then the table without row names, and returns invisibly.
#'   `FALSE` returns the table.
#'
#' @return A data frame with one row per check and three columns: `check`
#'   (character), `max_error` (numeric) and `status` (character, `"OK"` or
#'   `"FAILED"`). Ten rows, or twelve when the smoother declares a
#'   `tau_correction`. Returned invisibly when `verbose` is `TRUE`.
#'
#' @section Errors:
#' \pkg{numDeriv} is in `Suggests` and every derivative row needs it, so the
#' function stops when it is not installed.
#'
#' @examples
#' res <- check_abs_smoother(smooth_probit(h = 0.3))
#' all(res$status == "OK")
#'
#' # Ten rows for a smoother with no scale correction, twelve with one.
#' nrow(check_abs_smoother(smooth_quintic(h = 0.3), verbose = FALSE))
#' nrow(check_abs_smoother(smooth_probit(h = 0.3), verbose = FALSE))
#'
#' # The error this caught while the quintic was being written: its own
#' # third derivative with a 4 in the denominator where the algebra gives 2.
#' broken <- smooth_quintic(h = 0.3)
#' fns <- broken@s
#' fns[[4]] <- function(u, width) {
#'   t <- u / width
#'   ifelse(abs(u) < width, -15 * t * (1 - t^2) / (4 * width^2), 0)
#' }
#' broken@s <- fns
#' bad <- check_abs_smoother(broken, verbose = FALSE)
#' bad[bad$status != "OK", c("check", "max_error")]
#'
#' @seealso [abs_smoother()] for the contract being checked,
#'   [smoother_deriv()] for the derivatives it reads,
#'   `penalties7::check_penalty()` for the sibling on a penalty.
#' @export
check_abs_smoother <- function(smoother, width = NULL, tol = 1e-6,
                               verbose = TRUE) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("check_abs_smoother() needs the numDeriv package.", call. = FALSE)
  }
  if (!S7::S7_inherits(smoother, abs_smoother)) {
    stop("'smoother' must be an abs_smoother.", call. = FALSE)
  }
  if (is.null(width)) width <- smoother@width
  if (is.null(width)) width <- 0.5

  # the intrinsic scale: s(0) is of order the transition width for any
  # smoother of the absolute value, so the grid needs no reading of what
  # the width parameter means
  w0 <- max(smoother_deriv(smoother, 0, width, 0L), sqrt(width) * 1e-3)
  # off the seam and off zero: a piecewise smoother has measure-zero points
  # where a one-sided derivative is read, and a difference straddling one
  # compares nothing
  grid <- w0 * c(-6.3, -2.7, -1.13, -0.41, 0.17, 0.59, 1.21, 3.1, 7.7)

  checks <- list()
  add <- function(name, gap) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, max_error = gap,
      status = if (is.finite(gap) && gap < tol) "OK" else "FAILED"
    )
  }

  s0 <- smoother_deriv(smoother, grid, width, 0L)
  add("s is even", max(abs(s0 - smoother_deriv(smoother, -grid, width, 0L))) /
        max(1, max(abs(s0))))
  s1 <- smoother_deriv(smoother, grid, width, 1L)
  add("s' is odd", max(abs(s1 + smoother_deriv(smoother, -grid, width, 1L))))
  add("|s'| bounded by one", max(0, max(abs(s1)) - 1))
  add("s convex", max(0, -min(smoother_deriv(smoother, grid, width, 2L))))
  # the excess over |u| leaves the transition: at ten intrinsic widths it is
  # at most a fifth of what it is at the kink, which the hyperbolic's
  # polynomial tails still satisfy (excess c / (4|u|), a twentieth there)
  e0 <- smoother_deriv(smoother, 0, width, 0L)
  e10 <- smoother_deriv(smoother, 10 * w0, width, 0L) - 10 * w0
  add("matches |u| in the tails", max(0, e10 - e0 / 5))

  for (k in 1:5) {
    fk <- vapply(grid, function(u) {
      numDeriv::grad(function(v) smoother_deriv(smoother, v, width, k - 1L),
                     u)
    }, numeric(1))
    ak <- smoother_deriv(smoother, grid, width, k)
    add(sprintf("order %d vs numDeriv on order %d", k, k - 1L),
        max(abs(ak - fk)) / max(1, max(abs(fk))))
  }

  tc <- smoother@tau_correction
  if (!is.null(tc)) {
    tau <- 3 * w0
    add("tau_correction bounded by tau", max(0, tc(tau, width) - tau))
    add("tau_correction is the identity at width zero",
        abs(tc(tau, width * 1e-8) - tau) / tau)
  }

  out <- do.call(rbind, checks)
  rownames(out) <- NULL
  if (verbose) {
    cat(sprintf("check_abs_smoother: %s at %s = %g\n",
                smoother@smoother_name, smoother@width_name, width))
    print(out, row.names = FALSE)
    return(invisible(out))
  }
  out
}
