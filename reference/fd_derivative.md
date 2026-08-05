# One Stencil, Applied

The `order`-th derivative of `f` at `x` from a single finite-difference
stencil: the weighted sum of function values at the stencil's nodes,
divided by \\h^{d}\\.

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

  The derivative order, 1 to 4.

- h:

  A numeric vector of steps;
  [`fd_step`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
  when missing.

- accuracy:

  The order of the error term. Two gives the classical compact stencils;
  four gives the five-point first and second derivatives.

- side:

  `"central"` away from boundaries, `"forward"` or `"backward"` where a
  symmetric stencil would leave the domain.

## Value

A numeric vector of the same length as `x`.

## Details

This is the applicator every fallback in the toolkit speaks through, and
it enforces the one rule they share: one stencil of the requested order,
never a composition of lower-order differences. What it deliberately
does **not** choose is the policy around it – which order to fall back
from, when a reference can be trusted, what to do at a domain boundary
beyond keeping the nodes inside. Those belong to the callers, who know
what they are differentiating.

`f` must be vectorized in its argument; `x` and `h` may be vectors, and
the stencil is applied elementwise.

## See also

[`fd_weights`](https://statmodels7.github.io/numericals7/reference/fd_weights.md),
[`fd_step`](https://statmodels7.github.io/numericals7/reference/fd_step.md)

## Examples

``` r
# the third derivative of exp is exp
fd_derivative(exp, 1, order = 3)
#> [1] 2.718282
exp(1)
#> [1] 2.718282

# a five-point first derivative is fourth-order accurate
fd_derivative(sin, 0.7, order = 1, accuracy = 4) - cos(0.7)
#> [1] 1.665335e-14

# one-sided at a boundary
fd_derivative(sqrt, 1e-4, order = 1, side = "forward",
              h = fd_step(1e-4, 1, bounds = c(0, Inf)))
#> [1] 49.9589
```
