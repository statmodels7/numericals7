# Evaluate a Smoother or One of Its Derivatives

Returns \\d^k s/du^k\\ at the given points, for \\k\\ from 0 to 5. The
width is the smoother's own unless one is supplied, which is how a
consumer that resolved the width at build evaluates without mutating the
object.

## Usage

``` r
smoother_deriv(smoother, u, width = NULL, order = 0L)
```

## Arguments

- smoother:

  An
  [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md).
  Anything else is rejected.

- u:

  A numeric vector of points.

- width:

  The width, a single positive number or one per point, or `NULL` (the
  default) for the smoother's own. A smoother built with `width = NULL`
  and given none here is rejected, with the two ways to supply one
  named.

- order:

  The derivative order, a whole number from 0 to 5. `0L` by default,
  which is \\s\\ itself. Anything outside that range is rejected.

## Value

A numeric vector as long as `u`.

## Details

A vector `width` is one value per point, which is how a per-group width
arrives: a break-point term with `per_group = TRUE` resolves one width
per subject and passes the whole vector.

The derivatives are the ones the smoother carries, evaluated directly.
None is differenced, so order five is as accurate as order zero.

## See also

[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
for the contract,
[`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)
for resolving a width,
[`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
for verifying the derivatives.

## Examples

``` r
# The smooth sign: odd, bounded by one, and passing through zero.
smoother_deriv(smooth_hyperbolic(c = 1), 0:3, order = 1)
#> [1] 0.0000000 0.7071068 0.8944272 0.9486833

# A per-point width, as a per-group smoother supplies.
smoother_deriv(smooth_probit(), c(-1, 0, 1), width = c(0.1, 0.5, 1),
               order = 0)
#> [1] 1.0000000 0.3989423 1.1666309

# An unresolved width is an error rather than a guess.
try(smoother_deriv(smooth_probit(), 1))
#> Error : the smoother's width is unresolved; supply 'width' or construct the smoother with one.
```
