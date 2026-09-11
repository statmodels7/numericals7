# Resolve a Smoother's Width from a Spacing

Returns the width the smoother carries where it carries one, and
otherwise the given spacing carried onto the width parameter's own
scale. This is what a consumer calls at build: a break-point term hands
it the median spacing of its covariate, the smallest transition the data
can tell from a step.

## Usage

``` r
smoother_width(smoother, spacing)
```

## Arguments

- smoother:

  An
  [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md).
  Anything else is rejected.

- spacing:

  A spacing in covariate units, one value or one per group. Every entry
  must be positive and none may be missing; a smoother that already
  carries a width never reaches the check.

## Value

The width, on the smoother's own scale. The same length as `spacing`
when resolved from it, and a single number when the smoother carries
one.

## Details

A smoother constructed with an explicit width keeps it, and `spacing` is
then not even looked at. A smoother constructed with `NULL` takes the
spacing through its own `width_from_spacing`, which is the identity for
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
and
[`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md),
whose parameter is a length, and the square for
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md),
whose parameter is a squared one.

The result is nothing about whether the width is large enough for the
arithmetic;
[`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
answers that separately, and a consumer takes the larger of the two.

## See also

[`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
for the lower bound the arithmetic imposes,
[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
for what a width means.

## Examples

``` r
# A length-parametrized smoother takes the spacing as it stands.
smoother_width(smooth_probit(), 0.3)
#> [1] 0.3

# The hyperbolic's parameter is a squared length, so it is squared.
smoother_width(smooth_hyperbolic(), 0.3)
#> [1] 0.09

# A smoother that carries a width keeps it, whatever the spacing.
smoother_width(smooth_probit(h = 0.5), 0.3)
#> [1] 0.5

# One width per group.
smoother_width(smooth_probit(per_group = TRUE), c(0.2, 0.4, 0.35))
#> [1] 0.20 0.40 0.35
```
