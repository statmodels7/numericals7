# A Step Size for One Stencil

The magnitude-scaled step that balances truncation against rounding for
a single stencil of the given order and accuracy, clamped so the whole
stencil stays inside the domain when bounds are given.

## Usage

``` r
fd_step(x, order, accuracy = 2L, bounds = NULL)
```

## Arguments

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order.

- accuracy:

  The accuracy the stencil will be built at.

- bounds:

  An optional length-two numeric vector of domain bounds.

## Value

A numeric vector of steps, the same length as `x`.

## Details

The step is \\\varepsilon^{1/(d+a)}\max(1, \|x\|)\\: truncation grows
like \\h^{a}\\ and rounding like \\\varepsilon h^{-d}\\, and this is
where they balance. Near a finite bound the step is shrunk so that the
farthest node of the stencil stays strictly inside, since a node outside
the domain does not make a derivative inaccurate, it makes it `NaN`.

## See also

[`fd_derivative`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)

## Examples

``` r
fd_step(c(0.5, 1000), 2)
#> [1] 0.0001220703 0.1220703125
# near a boundary the stencil is kept inside
fd_step(0.01, 2, bounds = c(0, Inf))
#> [1] 0.0001220703
```
