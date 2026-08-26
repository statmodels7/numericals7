# One Stencil, Applied

Estimates the `order`-th derivative of `f` at `x` from a single
finite-difference stencil: the weighted sum of function values at the
stencil's nodes, divided by \\h^{d}\\. It picks the offsets with
[`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md),
the weights with
[`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
and, unless given one, the step with
[`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md).

## Usage

``` r
fd_derivative(
  f,
  x,
  order,
  h = NULL,
  accuracy = 2L,
  side = c("central", "forward", "backward")
)
```

## Arguments

- f:

  A vectorized function of one numeric argument.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order. One to four is what the toolkit uses and what
  the step rule is tuned for. A higher order is accepted and sized by
  [`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md),
  but nothing is claimed about its accuracy.

- h:

  A numeric vector of steps, recycled against `x`. `NULL`, the default,
  takes
  [`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
  at the same order and accuracy, with no bounds; pass a step explicitly
  when the domain is bounded.

- accuracy:

  The order of the error term, `2` by default. Two gives the classical
  compact stencils, four the five-point first and second derivatives.

- side:

  `"central"` away from boundaries, `"forward"` or `"backward"` where a
  symmetric stencil would leave the domain. All three use \\2r + 1\\
  nodes, so a one-sided estimate costs the same and is one order less
  accurate.

## Value

A numeric vector of the same length as `x`.

## One stencil, never nested

This is the applicator every numerical fallback in the toolkit speaks
through, and it enforces the rule they share: one stencil of the order
requested, never a composition of lower-order differences. Each
numerical differentiation multiplies the error of the one before it, so
a fourth derivative reached by four nested first differences is noise.

## What it deliberately leaves to the caller

The policy around the stencil. Which order to fall back from, when a
reference can be trusted, and what to do at a domain boundary beyond
keeping the nodes inside are all decisions that need to know what is
being differentiated, and this function does not.

## Vectorization

`f` must be vectorized in its argument. `x` and `h` may be vectors, and
the stencil is applied elementwise, so a whole vector of points costs
\\2r + 1\\ calls to `f` and no more.

## See also

[`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md),
[`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
and
[`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md),
the three pieces this assembles.

## Examples

``` r
# The third derivative of exp is exp.
fd_derivative(exp, 1, order = 3) - exp(1)
#> [1] -1.201531e-07

# Accuracy four is worth about two decades here.
c(acc2 = fd_derivative(sin, 0.7, 1) - cos(0.7),
  acc4 = fd_derivative(sin, 0.7, 1, accuracy = 4) - cos(0.7))
#>         acc2         acc4 
#> 3.798628e-12 1.665335e-14 

# At a boundary, one-sided with a step that keeps the nodes inside. The
# central stencil would reach below zero, where sqrt is not defined.
fd_derivative(sqrt, 1e-4, order = 1, side = "forward",
              h = fd_step(1e-4, 1, bounds = c(0, Inf)))
#> [1] 49.9589
0.5 / sqrt(1e-4)
#> [1] 50

# Vectorized in the point: one call per node, not one per point.
fd_derivative(sin, c(0, pi / 4, pi / 2), order = 1) - cos(c(0, pi / 4, pi / 2))
#> [1] -6.111445e-12 -7.690515e-12 -6.123234e-17
```
