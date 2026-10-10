# Resolve a Smoother's Width from a Spacing

Returns the width carried by the smoother when it has one, and otherwise
the given spacing carried onto the width parameter's own scale. A
consumer calls this function at build: a break-point term passes it the
median spacing of its covariate, which is the smallest transition that
the data can distinguish from a step.

## Usage

``` r
smoother_width(smoother, spacing, max_gap = NULL)
```

## Arguments

- smoother:

  An
  [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md).
  Anything else is rejected.

- spacing:

  A spacing in covariate units, one value or one per group. Every entry
  must be positive and none may be missing; the check is not made when
  the smoother already carries a width.

- max_gap:

  `NULL` (the default), or the largest gap between consecutive distinct
  values of the covariate over the range a break-point may take, one
  value or one per entry of `spacing`. Read only by a smoother with an
  `exact_radius`.

## Value

The width, on the smoother's own scale. The same length as `spacing`
when resolved from it, and a single number when the smoother carries
one.

## Details

A smoother constructed with an explicit width keeps it, and `spacing` is
then ignored. A smoother constructed with `NULL` takes the spacing
through its own `width_from_spacing`: the identity for
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md),
the square for
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md),
whose parameter is a squared length, and \\5/(2\log 2)\\ times the
spacing for
[`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md),
for the reason given on its page.

A smoother that declares an `exact_radius` is equal to \\\lvert
u\rvert\\ beyond that radius, so a break-point smoothed with it has no
curvature where no observation falls within the radius, and a gap in the
covariate wider than twice the radius can hold such a break-point. Where
`max_gap` is given, the width is therefore raised until the radius is at
least \\0.55\\ times it, the radius being taken to grow in proportion to
the width, which holds for
[`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md),
whose radius is its width. The median spacing does not bound the largest
gap: among \\n\\ uniform points the largest gap is about \\(\log n +
\gamma)/\log 2\\ median spacings, with \\\gamma\\ Euler's constant (9.5
at \\n = 400\\), while the quintic's width is 3.61 median spacings. A
smoother without an `exact_radius` ignores `max_gap`.

The result does not indicate whether the width is large enough for the
arithmetic;
[`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
returns that bound separately, and a consumer takes the larger of the
two.

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

# The quintic is raised to 0.55 times a large gap; the probit is not.
smoother_width(smooth_quintic(), 0.01, max_gap = 0.1)
#> [1] 0.055
smoother_width(smooth_probit(), 0.01, max_gap = 0.1)
#> [1] 0.01
```
