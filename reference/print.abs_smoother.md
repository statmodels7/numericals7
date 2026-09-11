# Print a Smoother of the Absolute Value

Writes one line naming the smoother and its width, and adds a line for
each of the two optional properties it declares: a scale correction for
a random break-point, and exactness outside the transition.

## Arguments

- x:

  An
  [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
  object.

- ...:

  Unused, accepted for consistency with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly.

## Details

The width line reads `h resolved at build` where the smoother carries
none, and names the width parameter by its own name, so a hyperbolic
smoother shows `c` where a probit shows `h`. `per group` is appended
when the width is to be resolved per group.

## Examples

``` r
smooth_probit()
#> <abs_smoother> probit: h resolved at build
#>   declares a scale correction for a random break-point
smooth_probit(h = 0.2, per_group = TRUE)
#> <abs_smoother> probit: h = 0.2, per group
#>   declares a scale correction for a random break-point
smooth_hyperbolic()
#> <abs_smoother> hyperbolic: c resolved at build
smooth_quintic(h = 0.4)
#> <abs_smoother> quintic: h = 0.4
#>   exact outside the transition
```
