# Partial Derivatives by One Tensor-Product Stencil per Component

Computes partial derivatives of a function of several coordinates. Each
component is a single linear combination of values of the highest
analytic quantity below it (the function itself, or a derivative
declared through
[`fd_known()`](https://statmodels7.github.io/numericals7/reference/fd_known.md)):
a tensor product of central one-dimensional stencils, one in each
direction that still has to be differentiated, and never a difference of
a difference.

## Usage

``` r
fd_tensor(
  f,
  x,
  order = NULL,
  components = NULL,
  known = list(),
  bounds = NULL,
  accuracy = 4L,
  h_rel = NULL
)
```

## Arguments

- f:

  A function of one argument, the named list of coordinates, returning a
  numeric vector of length 1 or \\n\\. When every coordinate has length
  1, \\n\\ is the length of the values (a response held by `f`, for
  instance).

- x:

  A named list of \\p\\ numeric coordinates, each of length 1 or a
  common length \\n\\. The derivatives are taken at these values,
  elementwise along the vectors of length \\n\\.

- order:

  A single positive whole number: every component of this total order is
  computed, in the order of
  [`tuple_indices(p, order)`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md).
  Exactly one of `order` and `components` is given.

- components:

  An integer matrix of counts, one row per component and one column per
  coordinate.

- known:

  A list of
  [`fd_known()`](https://statmodels7.github.io/numericals7/reference/fd_known.md)
  objects, or a single one.

- bounds:

  `NULL`, or a list with one entry per coordinate (matched by name when
  named), each `NULL` or a numeric vector `c(lower, upper)`.

- accuracy:

  The order of accuracy of every one-dimensional stencil, a positive
  even whole number. The default, 4, makes the truncation error of order
  \\h^4\\; 2 uses fewer nodes in each direction, and the number of
  evaluations falls by a factor of two to four.

- h_rel:

  `NULL` for the step above, or a single positive number that replaces
  \\\varepsilon^{1/(r + a)}\\.

## Value

A named list of numeric vectors of length \\n\\, one per component,
named by the coordinates of each multi-index joined with `_` (`"a_a_b"`
for \\\partial^3 f / \partial a^2 \partial b\\). The attribute
`"fd_source"` is a data frame with one row per component and the columns
`component`, `source` (0 for `f`, \\k\\ for `known[[k]]`), `index` (the
row of that quantity's multi-indices) and `residual` (the order left to
the stencil).

## The source of each component

A component is a multi-index \\\alpha\\, a vector of counts over the
coordinates. Every declared multi-index \\\beta \le \alpha\\ can serve
as its source, with residual \\\alpha - \beta\\ and residual order \\r =
\|\alpha\| - \|\beta\|\\; the function itself is the source with \\\beta
= 0\\. The source with the smallest \\r\\ is used. Among sources of
equal \\r\\, the quantity that declares more derivatives is preferred,
and then the one declared first in `known`. Within one quantity the
analytic part takes the most repeated coordinates of \\\alpha\\ first,
so that the residual is spread over the directions as evenly as
possible; among equally even residuals the earlier coordinates go to the
analytic part. The attribute `"fd_source"` of the result records the
choice.

## The stencil

With residual \\\rho = \alpha - \beta\\ and source \\g\\, the component
is \$\$\sum\_{k} \Bigl(\prod\_{j:\\\rho_j \> 0} w^{(\rho_j)}\_{k_j} /
h_j^{\rho_j}\Bigr)\\ g(x_1 + k_1 h_1, \dots, x_p + k_p h_p),\$\$ where
\\w^{(m)}\\ are the central weights of order \\m\\ on the nodes of
[`fd_offsets(m, accuracy)`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md),
computed by Fornberg's recursion and exactly symmetric, and the sum runs
over the nodes with a nonzero weight. A node shared by several
components is evaluated once.

## The step

A component of residual order \\r\\ uses, in direction \\j\\, the step
\\h_j = \varepsilon^{1/(r + a)}\max(1, \|x_j\|)\\, with \\a\\ the
accuracy: rounding grows as \\\varepsilon / h^r\\ and the truncation as
\\h^a\\, and the two balance there. `h_rel` replaces
\\\varepsilon^{1/(r + a)}\\. Where `bounds` gives a finite bound, the
step is reduced to \\0.49\\d / R_j\\, with \\d\\ the distance to the
bound and \\R_j\\ the largest offset that the stencils of that residual
order use in direction \\j\\, so every node lies inside the domain.

## The computation

The plan (the source of each component, the distinct nodes and the
sparse list of weights) depends only on the multi-indices and the
accuracy. It is computed in C++ and cached for the session. The node
values are computed in R, one call of the source function per node, and
the weighted sums are formed in C++ over the vectors of length \\n\\.

## See also

[`fd_known()`](https://statmodels7.github.io/numericals7/reference/fd_known.md)
for the declarations;
[`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
and
[`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
for the one-dimensional rules;
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
for one coordinate.

## Examples

``` r
f <- function(x) x$a^3 * exp(x$b)
x <- list(a = c(0.5, 1, 2), b = 0.3)
d3 <- fd_tensor(f, x, order = 3)
# the exact a_a_b component is 6 a exp(b)
rbind(d3$a_a_b, 6 * x$a * exp(x$b))
#>          [,1]     [,2]     [,3]
#> [1,] 4.049576 8.099153 16.19831
#> [2,] 4.049576 8.099153 16.19831

# with the analytic gradient, each third derivative is a second difference
g <- fd_known(function(x) list(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b)),
              order = 1)
d3g <- fd_tensor(f, x, order = 3, known = g)
attr(d3g, "fd_source")
#>   component source index residual
#> 1     a_a_a      1     1        2
#> 2     a_a_b      1     1        2
#> 3     a_b_b      1     2        2
#> 4     b_b_b      1     2        2
rbind(d3g$a_a_b, 6 * x$a * exp(x$b))
#>          [,1]     [,2]     [,3]
#> [1,] 4.049576 8.099153 16.19831
#> [2,] 4.049576 8.099153 16.19831

# one component only, at second-order accuracy
fd_tensor(f, x, components = rbind(c(2, 1)), accuracy = 2)
#> $a_a_b
#> [1]  4.049577  8.099153 16.198307
#> 
#> attr(,"fd_source")
#>   component source index residual
#> 1     a_a_b      0     1        3
```
