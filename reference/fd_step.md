# A Step Size for One Stencil

Returns the step at which a single stencil of the given order and
accuracy balances truncation error against rounding error, scaled by the
magnitude of the point, and shrunk near a finite bound so that the whole
stencil stays inside the domain.

## Usage

``` r
fd_step(x, order, accuracy = 2L, bounds = NULL)
```

## Arguments

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order \\d\\.

- accuracy:

  The accuracy the stencil will be built at, `2` by default. Matches
  [`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)'s
  argument of the same name.

- bounds:

  An optional numeric vector of two domain bounds, either of which may
  be infinite. `NULL`, the default, applies no clamp.

## Value

A numeric vector of steps the same length as `x`, positive except at a
point sitting on a finite bound, where it is zero.

## Where the balance falls

\$\$h = \varepsilon^{1/(d + a)} \max(1, \lvert x \rvert).\$\$

Truncation grows like \\h^{a}\\ and rounding like \\\varepsilon
h^{-d}\\, and this is where the two meet. The factor \\\max(1, \lvert x
\rvert)\\ makes the step relative for a large argument and absolute for
a small one, so a point near zero does not get a step below the
resolution of its own neighborhood.

At the default accuracy the exponent is \\1/4\\ for a second derivative
and \\1/6\\ for a fourth, giving steps of about `1.2e-4` and `2.5e-3` at
\\x = 1\\. A high order wants a *large* step, since rounding is what
dominates there.

## Staying inside the domain

Given `bounds`, the step is shrunk so the farthest node of the stencil
stays strictly inside: a node outside the domain does not make a
derivative inaccurate, it makes it `NaN`. The margin is 0.49 of the
distance to the bound, divided by the reach.

## The one case to guard

A point sitting exactly on a finite bound gets a step of **zero**, and a
stencil divided by \\h^{d}\\ is then `NaN`. The evaluation point is the
caller's, so this is not checked here; either keep the point off the
bound or use a one-sided stencil with a step of your own.

## See also

[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md),
which calls this when no step is given,
[`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
for the reach the clamp divides by.

## Examples

``` r
# Magnitude-scaled: absolute near zero, relative far from it.
fd_step(c(0.5, 1000), 2)
#> [1] 0.0001220703 0.1220703125

# A higher order wants a larger step, rounding being what dominates.
c(order2 = fd_step(1, 2), order4 = fd_step(1, 4))
#>       order2       order4 
#> 0.0001220703 0.0024607833 

# Near a boundary the stencil is kept inside the domain.
fd_step(0.01, 2, bounds = c(0, Inf))
#> [1] 0.0001220703

# On the boundary the step is zero, which no stencil can use.
fd_step(0, 2, bounds = c(0, Inf))
#> [1] 0
```
