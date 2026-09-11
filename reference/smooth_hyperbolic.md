# The Hyperbolic Smoother of the Absolute Value

\\s(u) = \sqrt{u^2 + c}\\, the simplest smoother of the three and the
one with the worst tails. Its parameter \\c\\ is a **squared** length,
so the transition width in covariate units is \\\sqrt{c}\\.

## Usage

``` r
smooth_hyperbolic(c = NULL)
```

## Arguments

- c:

  The squared transition width, a single positive number, or `NULL` (the
  default) to be resolved at build as the **square** of the covariate's
  spacing, through
  [`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md).

## Value

An
[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
named `"hyperbolic"`, with `width_name` `"c"`, a `width_from_spacing`
that squares, and neither a `tau_correction` nor an `exact_radius`.

## The closed form

\$\$s(u) = \sqrt{u^2 + c}, \qquad s'(u) = \frac{u}{\sqrt{u^2 + c}},
\qquad s''(u) = \frac{c}{(u^2 + c)^{3/2}},\$\$

every order a rational function of \\u\\ over a half-integer power, so
nothing here needs a special function or a branch.

## The cost of the tails

The excess over \\\lvert u\rvert\\ decays only as \\c/(2\lvert
u\rvert)\\, so the smoothing bias spreads well away from the kink.
Measured at a transition width of `0.3`, so \\c = 0.09\\: the excess is
`4.4e-02` at \\u = 1\\ and `1.5e-02` at \\u = 3\\, against `6.7e-05` and
exactly `0` for
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
at the same width. A fit smoothed this way is perturbed everywhere, not
only near the break-point.

No convolution identity relates \\c\\ to the scale of a random
break-point, so `tau_correction` is `NULL` and a consumer reports the
apparent scale alone.

## See also

[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
for the contract,
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
for the recommended default,
[`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md)
for exactness outside the transition,
[`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
to verify it.

## Examples

``` r
# c is a squared length, so a transition of 0.3 is c = 0.09.
sm <- smooth_hyperbolic(c = 0.09)
smoother_deriv(sm, 0, order = 0)
#> [1] 0.3
sqrt(0.09)
#> [1] 0.3

# And resolving from a spacing squares it for you.
smoother_width(smooth_hyperbolic(), 0.3)
#> [1] 0.09

# The excess over |u| is still there three units out.
u <- c(0, 0.3, 1, 3)
smoother_deriv(sm, u, order = 0) - abs(u)
#> [1] 0.30000000 0.12426407 0.04403065 0.01496269
smoother_deriv(smooth_probit(h = 0.3), u, order = 0) - abs(u)
#> [1] 2.393654e-01 4.998928e-02 6.724673e-05 0.000000e+00
```
