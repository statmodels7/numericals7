# S7 Class for Smoothers of the Absolute Value

A smooth replacement \\s(u)\\ for \\\lvert u \rvert\\, carrying its
derivatives in \\u\\ up to order five as functions. Every non-smooth
primitive the toolkit uses is generated from the absolute value:

\$\$\operatorname{sign}(u) = \frac{d\lvert u\rvert}{du}, \qquad
\mathbb{1}(u \geq 0) = \frac{1 + \operatorname{sign}(u)}{2}, \qquad
(u)\_{+} = \frac{u + \lvert u\rvert}{2},\$\$

so one contract serves them all. The smooth sign is \\s'\\, the smooth
step \\(1 + s'(u))/2\\ and the smooth hinge \\(u + s(u))/2\\ follow by
composition, and a consumer that replaces \\\lvert u\rvert\\ by \\s(u)\\
has replaced every one of them consistently.

## Usage

``` r
abs_smoother(
  smoother_name = character(0),
  width = NULL,
  width_name = "h",
  per_group = FALSE,
  s = list(),
  width_from_spacing = NULL,
  tau_correction = NULL,
  exact_radius = NULL
)
```

## Arguments

- smoother_name:

  A single non-empty string naming the smoother, used by
  [`print()`](https://rdrr.io/r/base/print.html) and by a consumer
  reporting which smoother a term carries.

- width:

  The transition width, a single positive number, or `NULL` to be
  resolved at build.

- width_name:

  What the width parameter is called, a single non-empty string: `"h"`
  for a length and `"c"` for a squared one. `"h"` by default.

- per_group:

  `TRUE` when the width is resolved per group, `FALSE` (the default) for
  one width throughout.

- s:

  A list of exactly six functions of `(u, width)`: \\s\\ and its
  derivatives in \\u\\ of orders one to five, in that order. Each must
  take at least two arguments; the validator checks the count and the
  arity.

- width_from_spacing:

  `NULL` (the identity) or a function carrying a spacing in covariate
  units onto the width parameter's own scale, which is the square for
  the hyperbolic.

- tau_correction:

  `NULL`, or a function `(tau, width)` returning the corrected scale of
  a random break-point.

- exact_radius:

  `NULL`, or a function of the width returning the radius beyond which
  \\s(u) = \lvert u\rvert\\ exactly. Only the quintic has one.

## Value

An S7 object of class `abs_smoother` carrying the eight properties
above.

## What a smoother buys

A break-point term smoothed this way becomes an ordinary nonlinear term
whose design block is the true Jacobian, so a random or penalized
development of its break-points becomes fittable. A kinked penalty
smoothed this way becomes a proper separable one, at the price of the
exact zeros the kink produced.

## The derivatives are functions

They are functions and not expressions. A piecewise smoother, the
quintic here, has branches that
[`stats::deriv()`](https://rdrr.io/r/stats/deriv.html) does not read,
and with the derivatives written out in the constructor those branches
are ordinary code. Each takes `(u, width)` and vectorizes in both, so a
per-group width is one value per observation.

## The width

`width` is the transition scale: `h` for a smoother whose parameter is a
length, the bent-cable reading of a transition of width \\h\\, and `c`
for the hyperbolic, whose parameter is a squared length. `NULL`, the
default of all three constructors, asks the consumer to resolve it from
the data at build through
[`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md);
a break-point term hands it the median spacing of its covariate.
`per_group` asks for one width per group where a grouping is available,
the validity window of a Laplace approximation being per-subject.

## The scale correction

`tau_correction` is a property of the mollifier and of no model. The
probit smoother satisfies an exact convolution identity, smoothing with
width \\h\\ being the same as convolving the break-point with \\N(0,
h^2)\\, so an apparent scale \\\tau\\ of a random break-point composes
as \\\tau^2\_{\mathrm{true}} = \tau^2 - h^2\\ and the smoother declares
that correction. Smoothers with no such identity declare `NULL`, and a
consumer then reports the apparent scale alone.

## References

Bacon, D. W. and Watts, D. G. (1971). Estimating the transition between
two intersecting straight lines. *Biometrika*, **58**(3), 525–534.

Tishler, A. and Zang, I. (1981). A new maximum likelihood algorithm for
piecewise regression. *Journal of the American Statistical Association*,
**76**(376), 980–987.

Seo, M. H. and Linton, O. (2007). A smoothed least squares estimator for
threshold regression models. *Journal of Econometrics*, **141**(2),
704–735.

## See also

[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md),
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)
and
[`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md)
for the three that ship,
[`smoother_deriv()`](https://statmodels7.github.io/numericals7/reference/smoother_deriv.md)
to evaluate one,
[`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)
and
[`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
for the width,
[`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
to verify one of your own.

## Examples

``` r
sm <- smooth_probit(h = 0.3)
S7::S7_inherits(sm, abs_smoother)
#> [1] TRUE
sm
#> <abs_smoother> probit: h = 0.3
#>   declares a scale correction for a random break-point

# The smooth sign, step and hinge all follow from the same object.
u <- c(-1, -0.1, 0, 0.1, 1)
smoother_deriv(sm, u, order = 1)                    # smooth sign
#> [1] -0.9991419 -0.2611173  0.0000000  0.2611173  0.9991419
(1 + smoother_deriv(sm, u, order = 1)) / 2          # smooth step
#> [1] 0.0004290603 0.3694413402 0.5000000000 0.6305586598 0.9995709397
(u + smoother_deriv(sm, u, order = 0)) / 2          # smooth hinge
#> [1] 3.362337e-05 7.627083e-02 1.196827e-01 1.762708e-01 1.000034e+00

# Three widths out they agree with the sharp versions to five decimals,
# and inside the transition they are what replaces them.
far <- c(-1.2, -0.9, 0.9, 1.2)
round(smoother_deriv(sm, far, order = 1) - sign(far), 5)
#> [1]  0.00006  0.00270 -0.00270 -0.00006
round((far + smoother_deriv(sm, far, order = 0)) / 2 - pmax(far, 0), 5)
#> [1] 0.00000 0.00011 0.00011 0.00000
```
