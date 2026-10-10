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
  max_panels = 4096L,
  rule = gauss_kronrod15()
)
```

## Arguments

- f:

  The integrand, a function of `x` and `i` as described above.

- lower, upper:

  Numeric vectors of endpoints, recycled to a common length, either of
  which may be infinite. Every `lower` must be strictly below its
  `upper`, or the call signals an error.

- atol, rtol:

  The absolute and relative error budgets per row, defaulting to `1e-10`
  and `1e-8`. A row is judged against the larger of the two, so `atol`
  governs an integral near zero and `rtol` a large one. With `atol = 0`
  an `rtol` below the floor of the rule signals an error; with a
  positive `atol` any `rtol` is accepted, `rtol = 0` included.

- max_depth:

  The greatest number of bisections one panel may undergo, `48` by
  default. It governs the integration of an endpoint singularity, and
  not every integrable singularity is reached: at the default the
  density of a gamma distribution with shape 0.5 is integrated and the
  one with shape 0.45 returns `NA`, while with `max_depth = 200` shape
  0.2 is integrated and shape 0.1 returns `NA`. A row past the budget
  returns `NA`.

- max_panels:

  The greatest number of panels one row may hold, `4096` by default. A
  row still over its budget once it holds that many returns `NA`. It
  bounds the memory of a row whose error estimate does not fall with
  bisection, and at the default it admits `cos(16000 x)` on the unit
  interval, which needs between 2048 and 4096 panels. The count is per
  row, so a row's result does not depend on the other rows of the call.

- rule:

  The embedded quadrature pair,
  [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  by default: a list of `nodes` on \\\[-1, 1\]\\ and the weights `wk`
  and `wg` of the two rules at those nodes, three numeric vectors of
  equal length.

## Value

A numeric vector of integrals, one per row, of the recycled length of
`lower` and `upper`. `NA` in any row that did not reach the requested
accuracy within `max_depth` and `max_panels`, with a warning that lists
the first eight of those rows.

## The integrand

`f(x, i)` receives a numeric matrix `x` of evaluation points and an
integer vector `i` with one entry per row of `x`, which gives the
parameter set to which that row belongs. It returns the values
elementwise, either as a matrix shaped like `x` or as a vector in
column-major order.

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

The budget is judged on the sum because the alternative of giving each
panel a share proportional to its length cannot integrate an endpoint
singularity at all: near such a point the error stays concentrated in
the innermost panel however deep the bisection goes, so a per-length
share demands of that panel an accuracy that no depth reaches. With the
sum, the smooth panels absorb the error of the panel at the singularity.

## Infinite endpoints

Infinite endpoints are mapped to finite ones by the rational transforms
\\x = a + t/(1-t)\\, \\x = b - t/(1-t)\\ and \\x = t/(1-t^2)\\, whose
Jacobians multiply the integrand. Rows of different kinds may share one
call, so a vector of endpoints mixing finite and infinite ones is
integrated in the same passes.

## Rows that do not converge

A row whose panels still exceed their budget at `max_depth`, or once it
holds `max_panels` panels, returns `NA`, and one warning lists the first
eight such rows (followed by an ellipsis if there are more). The value
`NA` indicates that the requested accuracy was not reached.

`max_depth` ends a row whose error is concentrated at a point, where
bisection narrows one panel after another, and `max_panels` ends a row
whose error is spread over its whole interval, from rounding or from an
oscillation faster than the panels resolve. In the second case every
panel carries a similar share of the error, so most panels are split at
every pass, and the count grows geometrically long before any panel
reaches `max_depth`.

## The floor of the rule

On a constant integrand the two rules of the pair differ by half the
difference of their weight sums on every panel, so no refinement brings
a row's relative error estimate below that value plus the rounding of
the sums. For
[`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
it is \\2.2 \times 10^{-16}\\ where the two computed sums agree to the
last bit, as on x86_64, and \\4.4 \times 10^{-16}\\ where they differ by
one unit in the last place of 2, as on the arm64 build of R for macOS. A
relative budget below it with `atol = 0` could never be met, and the
call signals an error naming the floor before the integrand is
evaluated.

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
#> [1] -4.344779e-09  4.634610e-10

# Their means, against the closed form.
g <- function(x, i) x * dgamma(x, shape = shp[i], rate = 1)
range(quad_vec(g, 0, rep(Inf, 30)) - shp)
#> [1] -1.889511e-11  2.317305e-10

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

# A relative budget below the floor of the rule cannot be met, and with no
# absolute budget beside it the call is refused before anything is evaluated.
try(quad_vec(function(x, i) matrix(1, nrow(x), ncol(x)), 0, 1,
             atol = 0, rtol = 1e-17))
#> Error : quad_vec: rtol = 1e-17 is below the floor of the quadrature rule, 2.22e-16,
#>   and atol = 0 leaves no absolute budget to stop at instead.
#>   Use rtol >= 2.22e-16, or give a positive atol.

# An oscillation far faster than the panels resolve stops at max_panels.
suppressWarnings(quad_vec(function(x, i) sin(1e9 * x), 0, 1,
                          max_panels = 256L))
#> [1] NA
```
