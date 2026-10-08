# The Quintic Smoother of the Absolute Value

\\\lvert u\rvert\\ replaced inside \\\[-h, h\]\\ by the even polynomial
matching it to third order at the seam, and **exact** outside. It is the
only smoother here whose bias is zero beyond a finite radius, and the
only one that is not analytic.

## Usage

``` r
smooth_quintic(h = NULL)
```

## Arguments

- h:

  The transition half-width, a single positive number, or `NULL` (the
  default) to be resolved at build from the covariate's spacing through
  [`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md),
  as \\5/(2\log 2)\\ times the spacing.

## Value

An
[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
named `"quintic"`, with `width_name` `"h"`, a `width_from_spacing`
multiplying the spacing by \\5/(2\log 2)\\, no `tau_correction`, and an
`exact_radius` equal to the width.

## The closed form

With \\t = u/h\\,

\$\$s(u) = \frac{h}{16}\\\bigl(5 + 15t^2 - 5t^4 + t^6\bigr) \quad
\text{for } \lvert t\rvert \< 1, \qquad s(u) = \lvert u\rvert \quad
\text{otherwise},\$\$

so \\s''(u) = 15(1 - t^2)^2/(8h)\\ vanishes to second order at \\\pm h\\
and the function is \\C^3\\, with a jump in the fourth derivative at the
seam. At the kink \\s(0) = 5h/16\\.

## Exact outside the transition

The smoothing bias is exactly zero outside \\\[-h, h\]\\, which is the
cleanest fixed-width choice: measured at \\h = 0.3\\ the excess over
\\\lvert u\rvert\\ is `0` at \\u = 1\\ and at \\u = 3\\, where
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
leaves `6.7e-05` and
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)
`4.4e-02`. The smoother declares an `exact_radius`, so a consumer can
say where its answer is the sharp model's.

The branches are why
[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
carries the derivatives as functions:
[`stats::deriv()`](https://rdrr.io/r/stats/deriv.html) does not read a
clamp, and written out in the constructor the branches are ordinary
code.

There is no convolution identity here, so `tau_correction` is `NULL`.

## The width resolved from data

Exactness outside the transition has a cost for a break-point. The
derivative of \\s\\ in the break-point is \\s''\\, which is zero outside
\\\[-h, h\]\\, so where no observation falls within \\h\\ of the fitted
break-point the break-point has no curvature at all and is not
identified. At a width of one median spacing this happens in 25 to 55
per cent of fits of a discontinuous break-point on a uniform covariate,
because a quarter of the gaps are wider than \\2h\\. The width resolved
from a spacing \\g\\ is therefore

\$\$h = \frac{5}{2\log 2}\\g \approx 3.61\\g,\$\$

which puts five observations inside \\\[-h, h\]\\ on average for a
uniform covariate, whose median spacing is \\\log 2\\ times its mean
spacing. Measured on the same fits, this removes the unidentified
break-point and brings its error to that of
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md),
at a fit error on a true sharp step 10 to 15 per cent above the
probit's. It does not guarantee an observation inside the width for a
covariate with large gaps, since the largest of \\n\\ gaps is about
\\\log_2 n\\ median spacings, so a consumer that passes the largest gap
to
[`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)
also has the width raised to at least \\0.55\\ times it. A width
supplied as `h` is used as it stands.

## Under an outer criterion

The fourth derivative jumps at \\\pm h\\ and the fifth is a point mass
there. An exact outer gradient of a marginal criterion reads the fourth
derivative of \\s\\ through a break-point term and an exact outer
Hessian reads the fifth, so over the hyperparameters the criterion is
smooth only between the values at which an observation crosses \\\psi
\pm h\\, and the Hessian leaves out the point masses. Measured on a
random break-point under REML, the jump of the gradient at such a
crossing is below the resolution of a step of \\10^{-4}\\ in the
hyperparameter, and the exact gradient agrees with a difference of the
criterion to \\2.7 \times 10^{-7}\\.
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md),
which is analytic, has neither effect and is the choice where a
criterion is estimated over a smoothed break-point.

## See also

[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
for the contract,
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
for the recommended default,
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)
for the simplest form,
[`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
to verify it.

## Examples

``` r
sm <- smooth_quintic(h = 0.5)

# Outside the transition it is the absolute value, exactly.
smoother_deriv(sm, 1, order = 0)
#> [1] 1
smoother_deriv(sm, c(0.6, 1, 4), order = 0) - c(0.6, 1, 4)
#> [1] 0 0 0

# At the kink the value is 5h/16.
smoother_deriv(sm, 0, order = 0)
#> [1] 0.15625
5 * 0.5 / 16
#> [1] 0.15625

# The second derivative vanishes to second order at the seam, which is
# what makes it C^3 there.
smoother_deriv(sm, c(0.4, 0.49, 0.5, 0.51), order = 2)
#> [1] 0.4860000 0.0058806 0.0000000 0.0000000
```
