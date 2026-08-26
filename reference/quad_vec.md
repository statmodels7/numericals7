# Integrate One Function at Many Parameter Values

Computes \\\int\_{a_i}^{b_i} f(x; \theta_i)\\\mathrm{d}x\\ for every row
\\i\\ at once. The nodes of every panel of every row reach `f` in a
single call per refinement pass, so the parameter index is a matrix
dimension and not a loop. The toolkit's integrals are almost always this
shape, one integrand at many parameter values, and a scalar integrator
called in a loop pays its overhead once per value.

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

  The integrand, obeying the contract above.

- lower, upper:

  Numeric vectors of endpoints, recycled to a common length, either of
  which may be infinite. Every `lower` must be strictly below its
  `upper`, or the call throws.

- atol, rtol:

  The absolute and relative error budgets per row, defaulting to `1e-10`
  and `1e-8`. A row is judged against the larger of the two, so `atol`
  governs an integral near zero and `rtol` a large one.

- max_depth:

  The greatest number of bisections one panel may undergo, `48` by
  default. It is the lever for an endpoint singularity, and the measured
  reach is narrower than "integrable" suggests: at the default a gamma
  density of shape 0.5 converges and one of shape 0.45 does not, while
  `max_depth = 200` reaches shape 0.2 and still not 0.1. A row past the
  budget returns `NA`.

- rule:

  The embedded quadrature pair,
  [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  by default. Any list of `nodes`, `wk` and `wg` of equal length serves.

## Value

A numeric vector of integrals, one per row, of the recycled length of
`lower` and `upper`. `NA` in any row that did not reach the requested
accuracy, with a warning naming those rows.

## The integrand contract

`f(x, i)` receives a numeric matrix `x` of evaluation points and an
integer vector `i` with one entry per row of `x`, saying which parameter
set that row belongs to. It returns the values elementwise, either as a
matrix shaped like `x` or as a vector in column-major order.

Elementwise recycling does the indexing. For a gamma mean at parameter
vectors of length \\n\\:

    f <- function(x, i) x * dgamma(x, shape = shp[i], rate = rt[i])

`shp[i]` has one entry per row of `x` and recycles down each column.

## Batched adaptivity

Each panel carries the error estimate of its Gauss-Kronrod pair. A row
converges when the *sum* of its panel errors fits the budget
\\\max(\mathrm{atol}, \mathrm{rtol}\\\lvert I_i \rvert)\\. Until then
the row's worst panels are bisected, and the splits from every row join
the next single evaluation, so one hard row refines its own panels
without serializing the others.

The budget is judged on the sum for a reason worth knowing, because the
obvious alternative fails. Giving each panel a share proportional to its
length cannot integrate an endpoint singularity at all: near such a
point the error stays concentrated in the innermost panel however deep
the bisection goes, so a per-length share demands of that panel an
accuracy no depth reaches. Judging the sum lets the smooth panels carry
the row.

## Infinite endpoints

Mapped to finite ones by the rational transforms \\x = a + t/(1-t)\\,
\\x = b - t/(1-t)\\ and \\x = t/(1-t^2)\\, whose Jacobians multiply the
integrand. Rows of different kinds may share one call, so a vector of
endpoints mixing finite and infinite costs nothing extra.

## A failure is reported as one

A row whose panels still exceed their budget at `max_depth` returns
`NA`, and one warning names every such row. An `NA` says the accuracy
was not reached; a plausible number would say nothing and be believed.

## See also

[`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
for the discrete counterpart,
[`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
for the default rule.

## Examples

``` r
# Thirty gamma densities integrate to one, in one call rather than thirty.
shp <- seq(0.5, 15, length.out = 30)
f <- function(x, i) dgamma(x, shape = shp[i], rate = 1)
range(quad_vec(f, lower = 0, upper = rep(Inf, 30)) - 1)
#> [1] -4.344782e-09  4.634579e-10

# Their means, against the closed form.
g <- function(x, i) x * dgamma(x, shape = shp[i], rate = 1)
range(quad_vec(g, 0, rep(Inf, 30)) - shp)
#> [1] -1.889955e-11  2.317291e-10

# A shape below one puts an integrable singularity at the origin. The
# sum-judged budget reaches shape 0.5 at the default depth, and a harsher
# one needs a deeper budget rather than a looser tolerance.
quad_vec(function(x, i) dgamma(x, shape = 0.5, rate = 1), 0, Inf)
#> [1] 1
quad_vec(function(x, i) dgamma(x, shape = 0.2, rate = 1), 0, Inf,
         max_depth = 200L)
#> [1] 1

# Endpoints of different kinds share a call.
quad_vec(function(x, i) dnorm(x), c(-Inf, -1, 0), c(0, 1, Inf))
#> [1] 0.5000000 0.6826895 0.5000000

# A divergent integral is refused, not estimated.
suppressWarnings(quad_vec(function(x, i) 1 / x, 0, 1))
#> [1] NA
```
