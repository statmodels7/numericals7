# Sum One Series at Many Parameter Values

Computes \\\sum\_{k \ge k_0} t(k; \theta_i)\\ for every row \\i\\ at
once. Terms are evaluated in blocks as one matrix, and a row retires as
soon as it converges, so a slow row does not keep the finished ones
paying. This is the discrete counterpart of
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
and it exists for the same reason: the toolkit's sums are one series at
many parameter values.

## Usage

``` r
series_vec(
  term,
  n,
  from = 0L,
  atol = 1e-12,
  rtol = 1e-10,
  max_terms = 100000L,
  block = 64L
)
```

## Arguments

- term:

  The term function, obeying the contract above.

- n:

  The number of parameter rows, a positive whole number. It fixes the
  length of the answer and the range `i` takes.

- from:

  The first summation index, `0` by default. Pass `1` for a series
  indexed from one.

- atol, rtol:

  The absolute and relative budgets per row, defaulting to `1e-12` and
  `1e-10`. A row is judged against the larger of the two, so `atol`
  governs a sum near zero and `rtol` a large one. Both are tighter than
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)'s,
  a term being cheaper than a panel.

- max_terms:

  The greatest number of terms one row may consume, `100000` by default.

- block:

  How many terms are evaluated per pass, `64` by default. It is also the
  window the tail guard judges growth over, so a very small block
  weakens the guard.

## Value

A numeric vector of sums, one per row, of length `n`. `NA` in any row
that did not converge within `max_terms`, with a warning naming those
rows.

## The term contract

`term(k, i)` receives two integer vectors of equal length and returns
the terms elementwise. `k` is the summation index and `i` says which
parameter set each term belongs to. A Poisson mass at a rate vector
`lam` is

    term <- function(k, i) dpois(k, lam[i])

## Convergence, row by row

A row retires when three conditions hold together. Its last block
contributed less than \\\max(\mathrm{atol}, \mathrm{rtol}\\\lvert S_i
\rvert)\\; the final term of that block is itself below the same budget;
and the terms are not growing across the block.

The last two are the tail guard, and the third is the one that earns its
place. A hump-shaped term can open with a block that sums to nearly
nothing and ends on a term smaller still, while the whole series lies
beyond it. `dpois(0:63, 300)` sums to 3.8e-62 and ends at 3.0e-62, and
every series there is still to come. That the terms are rising is the
signal that the block sits before the mode.

Terms are assumed eventually decreasing in magnitude, which every series
in the toolkit satisfies.

## A failure is reported as one

A row still unconverged after `max_terms` terms returns `NA`, and one
warning names every such row.

## See also

[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md),
the continuous counterpart.

## Examples

``` r
# Four geometric series against the closed form.
r <- c(0.1, 0.5, 0.9, 0.99)
series_vec(function(k, i) r[i]^k, n = 4) - 1 / (1 - r)
#> [1]  0.000000e+00  0.000000e+00 -2.131628e-14 -8.779850e-09

# Poisson masses sum to one at every rate, in one call.
lam <- c(0.5, 4, 60)
series_vec(function(k, i) dpois(k, lam[i]), n = 3)
#> [1] 1 1 1

# The tail guard at work. The first 64 terms of a Poisson at rate 300 sum
# to almost nothing and are still rising, so the row does not retire there.
sum(dpois(0:63, 300))
#> [1] 3.75793e-62
series_vec(function(k, i) dpois(k, 300), n = 1)
#> [1] 1

# A series indexed from one, against its closed form. Convergence has to be
# geometric or better: the terms of sum 1/k^2 decay too slowly to retire a
# row within max_terms, and that row comes back NA.
series_vec(function(k, i) 1 / factorial(k), n = 1, from = 1L) - (exp(1) - 1)
#> [1] 2.220446e-16

# A divergent series is refused, not estimated.
suppressWarnings(series_vec(function(k, i) 1 / (k + 1), n = 1,
                            max_terms = 500L))
#> [1] NA
```
