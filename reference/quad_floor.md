# The Smallest Relative Budget a Quadrature Rule Can Meet

Returns \\\lvert\sum w_k - \sum w_g\rvert / 2\\ plus one unit in the
last place: the relative error estimate an embedded pair gives on a
constant integrand, on every panel, together with the rounding of the
two sums.
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
refuses a relative budget below it when no absolute budget is given.

## Usage

``` r
quad_floor(rule)
```

## Arguments

- rule:

  An embedded pair: a list of `nodes`, `wk` and `wg` of equal length, as
  [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  returns.

## Value

A single positive number.

## Details

On a constant integrand both rules are exact up to their weight sums, so
the difference of the two estimates, which an adaptive routine reads as
the error, is the difference of the sums scaled by the panel's
half-width. Relative to the integral it is the same on every panel, and
bisection does not lower it. For
[`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
the two sums are 2 and the value is about \\2.2 \times 10^{-16}\\; with
the fifteen-decimal constants this package carried before 0.14.0 it was
\\3.7 \times 10^{-15}\\.

## See also

[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md),
[`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
