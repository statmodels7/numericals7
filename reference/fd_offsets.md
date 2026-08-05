# Stencil Offsets for a Derivative Order

The symmetric offsets used away from a boundary, and the one-sided ones
used where a symmetric stencil would not fit, sized from the derivative
order and the accuracy asked of it.

## Usage

``` r
fd_offsets(order, accuracy = 2L)
```

## Arguments

- order:

  The derivative order.

- accuracy:

  The order of the error term, a positive integer. Central stencils gain
  a free order on symmetry, so an odd request costs the same nodes as
  the even one above it.

## Value

A list with the `reach` and the three offset vectors `central`,
`forward` and `backward`.

## Details

The reach is \\r = \lceil (d + a)/2 \rceil - 1\\ for order \\d\\ and
accuracy \\a\\, giving \\2r + 1\\ nodes: at the default accuracy of two
this is the three-point stencil for the first and second derivatives and
the five-point one for the third and fourth, and at accuracy four it is
the five-point stencils for the first and second – every stencil the
toolkit's packages had written out by hand, from one formula.

## See also

[`fd_weights`](https://statmodels7.github.io/numericals7/reference/fd_weights.md),
[`fd_derivative`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)

## Examples

``` r
fd_offsets(1)               # three points
#> $reach
#> [1] 1
#> 
#> $central
#> [1] -1  0  1
#> 
#> $forward
#> [1] 0 1 2
#> 
#> $backward
#> [1] -2 -1  0
#> 
fd_offsets(1, accuracy = 4) # five points
#> $reach
#> [1] 2
#> 
#> $central
#> [1] -2 -1  0  1  2
#> 
#> $forward
#> [1] 0 1 2 3 4
#> 
#> $backward
#> [1] -4 -3 -2 -1  0
#> 
fd_offsets(4)$central       # the five-point fourth-derivative stencil
#> [1] -2 -1  0  1  2
```
