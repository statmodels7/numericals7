# An Analytic Quantity Declared to a Tensor Stencil

Declares a function that returns analytic partial derivatives of the
function differentiated by
[`fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.md),
together with the multi-indices of the derivatives it returns.
[`fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.md)
differences these values instead of the function itself wherever that
leaves a lower order to the stencil.

## Usage

``` r
fd_known(fun, order = NULL, indices = NULL)
```

## Arguments

- fun:

  A function of one argument, the named list of coordinates.

- order:

  A single non-negative whole number, the total order of every
  derivative `fun` returns, or `NULL`.

- indices:

  An integer matrix of counts with one row per returned derivative and
  one column per coordinate, or `NULL`. Exactly one of `order` and
  `indices` is given.

## Value

An object of class `fd_known`, with the properties `fun`, `order` and
`indices`.

## Details

A multi-index is a vector of counts, one per coordinate of `x`: the row
`c(2, 0, 1)` over the coordinates `(a, b, c)` is \\\partial^3 f /
\partial a^2 \partial c\\. The rows of `indices` are the multi-indices
of the values `fun` returns, in the order in which it returns them. With
`order = k` the rows are every multi-index of total order \\k\\, in the
order of
[`tuple_indices(p, k)`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
over the \\p\\ coordinates;
[`fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.md)
expands the shortcut from the length of its `x`.

`fun` is called with a named list of the same shape as `x` and returns
either a list of numeric vectors, one per row of `indices`, or a numeric
matrix with one column per row, or a single vector when there is one
row. Each vector has length 1 or \\n\\. The names of the returned list
are not read.

## See also

[`fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.md),
which reads these declarations.

## Examples

``` r
f <- function(x) x$a^3 * exp(x$b)
# the gradient, declared by its order
g <- fd_known(function(x) list(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b)),
              order = 1)
# the same quantity, declared by its multi-indices
g2 <- fd_known(function(x) cbind(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b)),
               indices = rbind(c(1, 0), c(0, 1)))
g2@indices
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```
