# Numerical kernels

Almost every R package that computes a derivative it has no formula for
writes its own finite difference, and almost every one that sums a
series writes its own stopping rule. They are written as internal
helpers, so nothing outside can reuse them, and each one is written to
the accuracy its author needed on the day. A stencil is a fixed piece of
arithmetic; so is a quadrature rule, so is the enumeration of the set
partitions of four elements. This package holds one copy of each and the
packages above it consume them, so an enumeration cannot disagree with
itself between two of them.

There are no classes here, deliberately. What follows is what the four
groups do and where each of them stops.

## Derivatives with no closed form

A finite-difference stencil is the solution of a small Vandermonde
system: given offsets, the weights that reproduce a derivative of a
given order exactly on polynomials of as high a degree as the offsets
allow.
[`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
solves it and returns the weights.

``` r

fd_weights(-1:1, order = 1)     # the textbook central difference
#> [1] -0.5  0.0  0.5
fd_weights(-2:2, order = 1)     # five points, accurate two orders further
#> [1]  0.08333333 -0.66666667  0.00000000  0.66666667 -0.08333333
```

[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
applies one to a function. `accuracy` is the order of the truncation
error, and buying more of it costs evaluations of the integrand and
nothing else:

``` r

f  <- function(x) exp(-x^2 / 2)
d1 <- function(x) -x * exp(-x^2 / 2)          # the closed form, for comparison

sapply(c(2, 4, 6), function(a)
  abs(fd_derivative(f, 0.7, order = 1, accuracy = a) - d1(0.7)))
#> [1] 7.842615e-12 1.512124e-13 8.881784e-15
```

The step is not a constant. Rounding grows with the order of the
difference, truncation falls with it, and the step that balances the two
is $`\varepsilon^{1/(k+2)}`$ scaled by the point:

``` r

data.frame(order = 1:4,
           step  = fd_step(1, 1:4),
           rule  = .Machine$double.eps^(1 / (1:4 + 2)))
#>   order         step         rule
#> 1     1 6.055454e-06 6.055454e-06
#> 2     2 1.220703e-04 1.220703e-04
#> 3     3 7.400960e-04 7.400960e-04
#> 4     4 2.460783e-03 2.460783e-03
```

so a fourth derivative is differenced at a step four hundred times a
first derivative’s, and asking for one at a first derivative’s step
would return noise.

### One stencil, never a chain

The rule the whole toolkit follows is that a derivative of order $`k`$
is reached by **one** stencil applied to the highest analytic order
available, never by differencing a difference. Two nested central
differences in the same variable multiply their rounding, and by the
third order the answer is noise: the identity function’s exactly zero
third derivative comes back of order one. Where a caller has an analytic
first derivative, the fourth derivative is one stencil of order three on
it, and the accuracy is several digits better than four stencils of
order one.

## Integrals and sums over a parameter vector

A modeling routine rarely wants one integral. It wants the same integral
at every parameter setting in a vector, and computing them one at a time
in a loop pays R’s call overhead once per setting.
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
takes a vector of endpoints and returns a vector of integrals,
evaluating the integrand once for every node of every panel of every
row.

The integrand receives a matrix of evaluation points and an integer
vector saying which parameter setting each row belongs to, so the
indexing is ordinary recycling:

``` r

a <- c(0.5, 1, 2, 5, 20)
got <- quad_vec(function(x, i) x^(a[i] - 1) * exp(-x),
                lower = rep(0, length(a)), upper = rep(Inf, length(a)))
rbind(quad_vec = got, gamma = gamma(a))
#>              [,1] [,2] [,3] [,4]         [,5]
#> quad_vec 1.772454    1    1   24 1.216451e+17
#> gamma    1.772454    1    1   24 1.216451e+17
max(abs(got / gamma(a) - 1))
#> [1] 4.344783e-09
```

Note that the number of rows comes from the **endpoints**, so scalar
endpoints give one integral however long the parameter vector is.
Infinite endpoints are mapped to finite ones by a rational transform
whose Jacobian multiplies the integrand, and rows of different kinds may
share one call.

Adaptivity is judged per row on the **sum** of that row’s panel errors,
as QUADPACK judges it. A budget allocated panel by panel cannot meet an
integrable endpoint singularity at any depth, the error there being
concentrated in the innermost panel however deep the bisection goes.

[`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
is the same idea for a sum, and its stopping rule needs three conditions
rather than two. A block sum that is small and a last term that is small
are both true long **before** the mode of a hump-shaped term:

``` r

# The first 64 terms of E[Y^2] for a Poisson mean of 300.
sum((0:63)^2 * dpois(0:63, 300))   # the block sum: tiny
#> [1] 1.479217e-58
63^2 * dpois(63, 300)              # the last term: tiny
#> [1] 1.179611e-58
300 + 300^2                        # what the sum actually is
#> [1] 90300
```

The third condition is that the terms are not growing across the block,
and it is what tells a premature block from a finished tail:

``` r

lam <- c(0.5, 3, 40, 300)
got <- series_vec(function(k, i) k^2 * dpois(k, lam[i]), n = length(lam))
max(abs(got / (lam + lam^2) - 1))
#> [1] 2.220446e-16
```

### A failure is reported as one

A row whose panels still exceed their budget at `max_depth` returns
`NA`, with one warning naming every such row:

``` r

suppressWarnings(
  quad_vec(function(x, i) 1 / x, lower = 0, upper = 1, max_depth = 6)
)
#> [1] NA
```

The integral diverges, so there is no accuracy to reach. An `NA` says
the budget was not met; a plausible number would say nothing and be
believed.

## Functions that leave the doubles

Several quantities a likelihood needs are finite and perfectly ordinary
while the expression normally used to compute them is not. Each of these
is written on the scale where it stays representable.

The inverse Mills ratio $`R(t) = \phi(t)/\Phi(t)`$ is asymptotic to
$`-t`$, so it is about 100 at $`t = -100`$. Both of its factors
underflow well before that:

``` r

t <- c(-5, -38, -100, -400)
m <- mills_ratio(t)
data.frame(t, naive = dnorm(t) / pnorm(t), mills = m$r, asymptote = -t)
#>      t     naive      mills asymptote
#> 1   -5  5.186504   5.186504         5
#> 2  -38 38.026280  38.026279        38
#> 3 -100       NaN 100.009998       100
#> 4 -400       NaN 400.002500       400
```

It returns the derivative alongside, which lies in $`(-1, 0)`$ and which
a caller would otherwise compute from `r` a second time. The modified
Bessel function of the first kind overflows at an argument no larger
than a fitted concentration:

``` r

x <- c(10, 500, 800, 5000)
data.frame(x, naive = besselI(x, 0), log_bessel = log_bessel_i(x, 0))
#>      x         naive  log_bessel
#> 1   10  2.815717e+03    7.942972
#> 2  500 2.504809e+215  495.974008
#> 3  800           Inf  795.738912
#> 4 5000           Inf 4994.822490
```

and the ratio $`I_1/I_0`$, which a von Mises score needs, is finite
everywhere while the ratio of the two scaled functions is not:

``` r

k <- c(1, 1e4, 1e5, 1e6)
data.frame(
  kappa = k,
  naive = besselI(k, 1, expon.scaled = TRUE) /
          besselI(k, 0, expon.scaled = TRUE),
  ratio = bessel_i_ratio(k)
)
#>   kappa    naive     ratio
#> 1 1e+00 0.446390 0.4463900
#> 2 1e+04 0.999950 0.9999500
#> 3 1e+05 0.999995 0.9999950
#> 4 1e+06      NaN 0.9999995
```

Owen’s T is here because a skew-normal distribution function is one
bounded one-dimensional quadrature of it, which beats integrating the
density over a semi-infinite interval:

``` r

h <- 0.7; a <- 1.3
owen_t(h, a)
#> [1] 0.1034353
integrate(function(x) exp(-h^2 * (1 + x^2) / 2) / (1 + x^2), 0, a)$value / (2 * pi)
#> [1] 0.1034353
```

## The enumerations

A higher-order chain rule is a sum over set partitions, a symmetric
derivative array is indexed by non-decreasing tuples, and the support of
a multinomial is the weak compositions. All three are fixed
combinatorial objects, and the reason they live here is that three
packages needed them and three copies could disagree.

``` r

length(set_partitions(4))          # the Bell number B4
#> [1] 15
choose(2 + 3 - 1, 3)               # tuples of order 3 on 2 parameters
#> [1] 4
length(tuple_indices(2, order = 3))
#> [1] 4
```

`set_partitions(n)` returns each partition as a list of blocks, from the
one block containing everything to the $`n`$ singletons.
`tuple_indices(d, order)` returns non-decreasing index vectors, which is
exactly the set of distinct components of a symmetric derivative array:

``` r

sapply(tuple_indices(2, order = 3), paste, collapse = "")
#> [1] "111" "112" "122" "222"
```

`compositions(n, k)` returns the ways of writing $`n`$ as an ordered sum
of $`k`$ non-negative parts, one per row:

``` r

compositions(3, 2)
#>      [,1] [,2]
#> [1,]    0    3
#> [2,]    1    2
#> [3,]    2    1
#> [4,]    3    0
rowSums(compositions(3, 2))
#> [1] 3 3 3 3
```

## Threads

[`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
carries the thread policy of the whole toolkit. It is an object passed
as an argument from `statmod()` and `fit_distrib()` down to the kernels,
with no global state, and it lives at the root of the dependency graph
because the packages holding the kernels cannot import the one at the
top.

``` r

th <- n_threads(threads = 4, workers = 2)
c(threads = thread_count(th), workers = worker_count(th))
#> threads workers 
#>       4       2
```

The guarantee on its page is that a kernel’s result does not depend on
the count bit for bit:

``` r

x <- seq(0.5, 50, length.out = 2000)
identical(log_bessel_i(x, 0, threads = 1), log_bessel_i(x, 0, threads = 4))
#> [1] TRUE
```

[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
deliberately takes no count. Its hybrid branch calls R’s own `besselK`,
which can warn, and a warning raised from a worker thread reads R’s
stack bounds against a foreign stack pointer and kills the process.
Stating the restriction beats offering an argument that is unsafe on one
branch.

## Summary

- **Stencils.**
  [`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
  and
  [`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
  give the arithmetic,
  [`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
  applies it,
  [`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
  gives the step that balances rounding against truncation. One stencil
  on the highest analytic order available, never a difference of a
  difference.
- **Batched engines.**
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  and
  [`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
  are vectorized over the parameter, take their row count from the
  endpoints, judge convergence on the sum of a row’s errors, and return
  `NA` for a row whose budget was not met.
- **Special functions.**
  [`mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.md),
  [`owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.md),
  [`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md),
  [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
  and
  [`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md),
  each written on the scale where it stays representable, and each
  returning its derivatives where a caller would otherwise recompute
  them.
- **Enumerations.**
  [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md),
  [`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
  and
  [`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md),
  one copy each for the whole toolkit.
- **Threads.**
  [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  is passed as an argument, and a kernel’s answer is bit-identical at
  any count.
