# numericals7

Every package of the [statmodels7](https://statmodels7.github.io)
toolkit needs the same numerical machinery: finite-difference stencils,
quadrature, special functions and the combinatorics behind a chain rule.
Before this package existed several packages carried their own copies
(the set-partition enumeration was written twice, finite differences
three times). [numericals7](https://statmodels7.github.io/numericals7/)
contains that machinery once, at the bottom of the toolkit, and every
package above it uses this copy.

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/numericals7")
```

## The enumerations behind a chain rule

A derivative of order four over several variables is a sum over
combinatorial objects. The enumerations are defined once, so every
package that uses them indexes the derivatives in the same order.

``` r

# the multi-indexes that key a derivative list (diagonal first at order 2)
tuple_indices(2, 2)
#> [[1]]
#> [1] 1 1
#> 
#> [[2]]
#> [1] 2 2
#> 
#> [[3]]
#> [1] 1 2

# the set partitions a Faa di Bruno formula sums over: 1, 2, 5, 15
lengths(lapply(1:4, set_partitions))
#> [1]  1  2  5 15

# the weak compositions of n into k parts -- the support of a multinomial
compositions(3, 2)
#>      [,1] [,2]
#> [1,]    0    3
#> [2,]    1    2
#> [3,]    2    1
#> [4,]    3    0
```

## One stencil library

[`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
builds the finite-difference weights for any derivative order at any
offsets from one Vandermonde system,
[`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
chooses the offsets,
[`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
the step, and
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
assembles the three. The toolkit applies one stencil to the highest
analytic order available, never a chain of first differences. The
packages that call these functions choose the order from which to
differentiate, and the functions here build and apply the stencil.

``` r

# the classical five-point first-derivative weights, and a derivative
fd_weights(-2:2, 1L) * 12
#> [1]  1 -8  0  8 -1
fd_derivative(exp, 1, order = 3L, accuracy = 4L) - exp(1)
#> [1] 3.076484e-09
```

[`fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.md)
extends the same rule to partial derivatives of any order in several
variables. It differentiates the highest quantity that the caller
declares with
[`fd_known()`](https://statmodels7.github.io/numericals7/reference/fd_known.md),
so each component takes one stencil of the order that remains, and the
attribute `"fd_source"` records which quantity each component came from.

``` r

f <- function(x) x$a^3 * exp(x$b)
x <- list(a = c(0.5, 1, 2), b = 0.3)
# The third derivative in a, a and b is computed from the function alone.
fd_tensor(f, x, order = 3)$a_a_b
#> [1]  4.049576  8.099153 16.198306
# With the analytic gradient declared, the same component is one second difference.
g <- fd_known(function(x) list(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b)),
              order = 1)
fd_tensor(f, x, order = 3, known = g)$a_a_b
#> [1]  4.049576  8.099153 16.198306
6 * x$a * exp(x$b)   # the exact value
#> [1]  4.049576  8.099153 16.198306
```

## Quadrature and series, vectorized over the parameters

A regression model gives every observation its own parameter value, so
the integrals that a fit needs come in families: the same integrand at
hundreds of parameter rows.
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
integrates all of them in one adaptive run, with a single call to the
integrand per refinement pass: the integrand receives a matrix of points
and a row index, so a thousand parameter values cost matrix evaluations
instead of a thousand adaptive runs.
[`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
does the same for infinite sums. A row that does not reach the requested
accuracy returns `NA` with a warning.

``` r

shp <- c(0.5, 2, 5, 50)
f <- function(x, i) x * dgamma(x, shape = shp[i], rate = 1)
quad_vec(f, 0, rep(Inf, 4))          # the means, one call
#> [1]  0.5  2.0  5.0 50.0

lam <- c(0.5, 4, 300)
series_vec(function(k, i) dpois(k, lam[i]), n = 3)   # masses sum to one
#> [1] 1 1 1
```

## The logarithm of the modified Bessel functions

Bessel functions recur across distributions: the von Mises carries
$`I_0`$, and the Poisson-inverse Gaussian carries $`K_{y-1/2}`$.
$`I_\nu(x)`$ overflows from about $`x = 709`$ and $`K_\nu(x)`$
underflows to zero from about $`x = 705`$, while the exponentially
scaled forms underflow past $`10^5`$ and at large orders. Following
Plesner, Sørensen and Hauberg (ICS 2024, arXiv:2409.08729),
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
and
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
work on the log scale, so the result is finite and accurate wherever the
logarithm of the function is representable. The `_derivs` variants add
four derivatives in the argument from the ratio identity and the Bessel
equation, at the cost of one more log-Bessel evaluation.

``` r

log_bessel_i(1e7, 2)              # I overflows near 709; its logarithm is finite
#> [1] 9999991
log_bessel_i(0.001, 5000)         # the scaled form loses this order entirely
#> [1] -75595.66
log_bessel_k_derivs(2, 0.5)$d1    # -1/(2x) - 1 at nu = 1/2, to rounding
#> [1] -1.25
```

## Special functions

Each of these functions was first written for a distribution, and each
is evaluated on the scale where it stays representable. The Mills ratio
is formed on the log scale, because the direct quotient loses precision
from about `t = -38` and is 0/0 from about `t = -38.6`. Owen’s T batches
every element into one
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
call. The Bessel ratio is evaluated from its power series, its continued
fraction and its asymptotic series, without the Bessel functions
themselves, so it is finite at any argument, and its derivatives and
those of its inverse are computed one order per function.

``` r

mills_ratio(-40)$r        # finite, close to 40
#> [1] 40.02497
owen_t(0, 1) * 2 * pi     # atan(1) = pi/4
#> [1] 0.7853982
bessel_i_ratio(1000)      # the unscaled ratio overflows past 700
#> [1] 0.9994999
bessel_i_ratio(bessel_i_ratio_inverse(0.7))
#> [1] 0.7
```

## Functions for worker threads

[`student_t_cdf()`](https://statmodels7.github.io/numericals7/reference/student_t_cdf.md)
and
[`bessel_k()`](https://statmodels7.github.io/numericals7/reference/bessel_k.md)
return the values of [`stats::pt()`](https://rdrr.io/r/stats/TDist.html)
and [`base::besselK()`](https://rdrr.io/r/base/Bessel.html) from
compiled code that never calls the R API, so the compiled kernels of the
toolkit can call them from a worker thread. Each one is also available
to C code through `R_GetCCallable()`. The `threads` argument of a fit
such as `distributions7::fit_distrib()` or `statmodels7::statmod()` sets
the number of threads inside one R process, which
[`local_threads()`](https://statmodels7.github.io/numericals7/reference/local_threads.md)
applies. The default of one thread takes the sequential code path, and
the compiled kernels return the same values at any number of threads
([`?local_threads`](https://statmodels7.github.io/numericals7/reference/local_threads.md)
describes the one exception). The number of separate R processes is the
`cores` argument of the functions that distribute independent fits.

``` r

# The values are those of pt() and besselK().
student_t_cdf(c(-2, 0, 3), df = 2.5) - pt(c(-2, 0, 3), df = 2.5)
#> [1] 0 0 0
bessel_k(c(0.5, 10), nu = 1.5) - besselK(c(0.5, 10), 1.5)
#> [1] 0 0
```

## Smoothers of the absolute value

An `abs_smoother` replaces $`\lvert u \rvert`$ with a smooth $`s(u)`$
and carries the derivatives of $`s`$ to fifth order as functions. With
it, a term with a break-point becomes an ordinary differentiable model.
[`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
is the recommended default: its excess over $`\lvert u \rvert`$ has
Gaussian tails, which fall to rounding error within a few widths, and
smoothing with width $`h`$ is a convolution of the break-point with
$`N(0, h^2)`$, for which the smoother supplies a closed-form scale
correction.

``` r

sm <- smooth_probit(h = 0.2)
u  <- c(-1, -0.1, 0, 0.1, 1)
rbind(abs = abs(u), smooth = smoother_deriv(sm, u, order = 0))
#>        [,1]      [,2]      [,3]      [,4] [,5]
#> abs       1 0.1000000 0.0000000 0.1000000    1
#> smooth    1 0.1791186 0.1595769 0.1791186    1
```

The width cannot be made arbitrarily small.
[`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
derives the smallest width that the design can carry, from the
observation that the Jacobian column near the kink is of order $`1/h`$
against columns of order the range of the data.

``` r

smoother_width_floor(sm, scale = 10)   # a covariate spanning ten units
#> [1] 1.490116e-07
```

[`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
validates a user-written smoother, comparing each order with one
numerical differentiation of the analytic order below it.
