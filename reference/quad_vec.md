# Integrate One Function at Many Parameter Values

\\\int\_{a_i}^{b_i} f(x; \theta_i)\\dx\\ for every row \\i\\ at once, by
matrix evaluation: the nodes of every panel of every row go into `f` in
a single call per refinement pass.

## Usage

``` r
quad_vec(
  f,
  lower,
  upper,
  atol = 1e-10,
  rtol = 1e-08,
  max_depth = 48L,
  rule = gauss_kronrod15()
)
```

## Arguments

- f:

  The integrand, as described above.

- lower, upper:

  Numeric vectors of endpoints, recycled to a common length; either may
  be infinite.

- atol, rtol:

  The absolute and relative error budgets per row.

- max_depth:

  The maximum number of bisections a panel may undergo. The default
  reaches the integrable endpoint singularities of the mild kind a
  density with shape below one has; a harsher one is rejected.

- rule:

  The embedded quadrature pair, by default
  [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md).

## Value

A numeric vector of integrals, one per row, with `NA` where the
requested accuracy was not reached.

## Details

**The integrand contract.** `f(x, i)` receives a numeric matrix `x` of
evaluation points and an integer vector `i`, one entry per row of `x`,
saying which parameter set that row belongs to. It returns the values
elementwise, as a matrix like `x` or as a vector in column-major order.
A caller holding parameter vectors of length \\n\\ writes, for a gamma
mean,


    f <- function(x, i) x * dgamma(x, shape = shp[i], rate = rt[i])

and the elementwise recycling does the rest: `shp[i]` has one entry per
row and recycles down each column of `x`.

**Batched adaptivity.** Each panel carries the error estimate of its
Gauss-Kronrod pair, and a row converges when the *sum* of its panel
errors fits the budget \\\max(\mathrm{atol}, \mathrm{rtol}\\\lvert I_i
\rvert)\\. Until then the row's worst panels are bisected, and the
splits from all rows join the next single evaluation, so one hard row
refines its own panels without serializing the others. Judging the sum,
rather than giving each panel a share of the budget proportional to its
length, is what lets an integrable endpoint singularity converge: near
such a point the error is concentrated however deep the bisection goes,
and a per-length share would demand of the innermost panel an accuracy
no depth can reach.

**Infinite endpoints** are mapped to finite ones with the rational
transforms \\x = a + t/(1-t)\\, \\x = b - t/(1-t)\\ and \\x =
t/(1-t^2)\\, whose Jacobians multiply the integrand; rows of different
kinds may share a call.

**Rejection over plausibility.** A row whose panels still exceed their
budget at `max_depth` returns `NA` with a warning naming it. An `NA`
names a failure; a plausible number would hide one.

## See also

[`series_vec`](https://statmodels7.github.io/numericals7/reference/series_vec.md),
[`gauss_kronrod15`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)

## Examples

``` r
# thirty gamma densities integrate to one, in one call
shp <- seq(0.5, 15, length.out = 30)
f <- function(x, i) dgamma(x, shape = shp[i], rate = 1)
quad_vec(f, lower = 0, upper = rep(Inf, 30))
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# and their means, against the closed form
g <- function(x, i) x * dgamma(x, shape = shp[i], rate = 1)
range(quad_vec(g, 0, rep(Inf, 30)) - shp)
#> [1] -1.889955e-11  2.317291e-10
```
