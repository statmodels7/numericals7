# Finite-Difference Weights for Any Stencil

The weights that combine function values at the given offsets into an
approximation of the derivative of the given order, obtained by solving
the Vandermonde system that makes the stencil exact on polynomials.

## Usage

``` r
fd_weights(offsets, order)
```

## Arguments

- offsets:

  A numeric vector of distinct stencil offsets, in units of the step.

- order:

  The derivative order, smaller than `length(offsets)`.

## Value

A numeric vector of weights, the same length as `offsets`.

## Details

A stencil on \\n\\ nodes is exact for polynomials up to degree \\n -
1\\, so its error is \\O(h^{n - d})\\ for the \\d\\-th derivative –
fourth-order accurate for a first derivative on five nodes, second-order
for a fourth derivative on the same five.

Building the weights this way, rather than composing lower-order
differences, is what keeps a high order usable: each numerical
differentiation multiplies the error of the one before it, so a fourth
derivative reached by four nested first differences is noise. One
stencil, never nested.

## See also

[`fd_offsets`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md),
[`fd_derivative`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)

## Examples

``` r
fd_weights(c(-1, 0, 1), 1)             # the central first difference
#> [1] -0.5  0.0  0.5
fd_weights(c(-1, 0, 1), 2)             # 1, -2, 1
#> [1]  1 -2  1
fd_weights(-2:2, 1) * 12               # the five-point first derivative
#> [1]  1 -8  0  8 -1
fd_weights(0:2, 1)                     # one-sided, for a boundary
#> [1] -1.5  2.0 -0.5
```
