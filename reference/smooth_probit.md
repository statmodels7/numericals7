# The Probit Smoother of the Absolute Value

\\s(u) = \mathbb{E}\lvert u + hZ\rvert\\ for \\Z\\ standard normal, the
one to reach for unless you have a reason not to. Its excess over
\\\lvert u\rvert\\ has gaussian tails, so the smoothing bias is confined
to a window of width \\h\\ around the kink, and it is the only one of
the three that declares a scale correction for a random break-point.

## Usage

``` r
smooth_probit(h = NULL, per_group = FALSE)
```

## Arguments

- h:

  The transition width, a single positive number, or `NULL` (the
  default) to be resolved at build from the covariate's spacing through
  [`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md).

- per_group:

  `TRUE` to resolve the width once per group, `FALSE` (the default) for
  one width throughout. The validity window of a Laplace approximation
  is per-subject, so a hierarchical break-point wants `TRUE`.

## Value

An
[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
named `"probit"`, with `width_name` `"h"`, no `width_from_spacing` (the
spacing is the width), a `tau_correction`, and no `exact_radius`.

## The closed form

\$\$s(u) = u\\\bigl(2\Phi(u/h) - 1\bigr) + 2h\\\phi(u/h), \qquad s'(u) =
2\Phi(u/h) - 1, \qquad s''(u) = \frac{2}{h}\\\phi(u/h),\$\$

with the third to fifth derivatives following from \\\phi'(z) =
-z\phi(z)\\. At the kink \\s(0) = 2h\phi(0) \approx 0.798h\\, which is
the intrinsic scale
[`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
places its grid by.

## Why the tails matter

The excess \\s(u) - \lvert u\rvert\\ decays like \\\phi(u/h)\\, so it is
gone within a few widths. Measured at \\h = 0.3\\: `6.7e-05` at \\u =
1\\ and exactly `0` at \\u = 3\\, against `4.4e-02` and `1.5e-02` for
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)
at the same transition width. A model smoothed this way is the sharp
model to within rounding a short distance from the break-point.

## The convolution identity

Smoothing the step with width \\h\\ is exactly convolving the
break-point with \\N(0, h^2)\\, so a random break-point of true scale
\\\tau\\ appears at \\\sqrt{\tau^2 + h^2}\\, and the smoother declares
\\\tau\_{\mathrm{true}} = \sqrt{\tau^2 - h^2}\\, floored at zero. No
other smoother here has such an identity.

## See also

[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
for the contract,
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)
and
[`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md)
for the alternatives,
[`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)
for resolving `h` from data,
[`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
to verify it.

## Examples

``` r
sm <- smooth_probit(h = 0.3)

# At the kink the value is 2 h phi(0), not zero.
smoother_deriv(sm, 0, order = 0)
#> [1] 0.2393654
2 * 0.3 * stats::dnorm(0)
#> [1] 0.2393654

# And the excess over |u| is gone within a few widths.
u <- c(0, 0.3, 1, 3)
smoother_deriv(sm, u, order = 0) - abs(u)
#> [1] 2.393654e-01 4.998928e-02 6.724673e-05 0.000000e+00

# The scale correction it declares, applied to an apparent tau.
sm@tau_correction(0.5, 0.231)
#> [1] 0.44344
sqrt(0.5^2 - 0.231^2)
#> [1] 0.44344
```
