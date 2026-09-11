# The Smallest Width a Consumer May Use

Returns the floor of the width, derived from the expression that binds.
A consumer resolving a width from data takes the larger of this and
[`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)'s
answer.

## Usage

``` r
smoother_width_floor(smoother, scale)
```

## Arguments

- smoother:

  An
  [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md).
  Anything else is rejected.

- scale:

  The scale of the covariate, its range, a single positive number.
  Anything else is rejected.

## Value

A single number, on the smoother's own width scale.

## Where the bound comes from

The derivatives of a smoother scale as \\s^{(k)} \sim h^{1-k}\\, so the
Jacobian column of a smoothed break-point carries \\s''(0)/2 \sim 1/h\\
against covariate columns of order the range \\D\\. Holding the design's
condition number below \\\epsilon^{-1/2}\\, which leaves half the digits
of a double to a QR of the design, gives

\$\$h \geq \sqrt{\epsilon}\\D,\$\$

carried onto the width parameter's own scale for a smoother parametrized
by a squared length. At a covariate range of 10 that is `1.49e-07` for
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
and its square, `2.22e-14`, for
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md).

The bound is derived, which is the toolkit's rule for a guard constant,
and it is the same argument the break-point schedule's own scaling floor
rests on.

## See also

[`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)
for the width the data suggests,
[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
for what a width means.

## Examples

``` r
# A length-parametrized smoother: sqrt(eps) times the range.
smoother_width_floor(smooth_probit(), scale = 10)
#> [1] 1.490116e-07
sqrt(.Machine$double.eps) * 10
#> [1] 1.490116e-07

# The hyperbolic's floor is the square of that, its parameter being one.
smoother_width_floor(smooth_hyperbolic(), scale = 10)
#> [1] 2.220446e-14

# A consumer takes the larger of the data's spacing and the floor.
max(smoother_width(smooth_probit(), 0.3),
    smoother_width_floor(smooth_probit(), scale = 10))
#> [1] 0.3
```
