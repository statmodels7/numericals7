# Sum One Series at Many Parameter Values

\\\sum\_{k \ge k_0} t(k; \theta_i)\\ for every row \\i\\ at once, in
blocks of terms evaluated as one matrix, with rows retiring as they
converge so that a slow row does not keep the finished ones paying.

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

  The term function, as described above.

- n:

  The number of parameter rows.

- from:

  The first summation index.

- atol, rtol:

  The absolute and relative budgets per row.

- max_terms:

  The largest number of terms a row may consume.

- block:

  How many terms are evaluated per pass.

## Value

A numeric vector of sums, one per row, with `NA` where the series did
not converge.

## Details

**The term contract.** `term(k, i)` receives two integer vectors of
equal length and returns the terms elementwise: `k` is the summation
index, `i` says which parameter set. A Poisson mass, for a rate vector
`lam`, is `function(k, i) dpois(k, lam[i])`.

**Convergence, per row.** A row retires when its last block contributed
less than \\\max(\mathrm{atol}, \mathrm{rtol}\\\lvert S_i \rvert)\\, the
final term of the block is itself below that budget, *and* the terms are
not growing across the block. The last two conditions are the tail
guard: a hump-shaped term – a Poisson mass at a large rate, say – can
open with a block that sums to nearly nothing and ends on a term smaller
still, while the whole series lies beyond it; that its terms are rising
is what says the block came before the mode rather than after it. The
terms are assumed eventually decreasing in magnitude, which every series
in the toolkit satisfies.

**Refusal over plausibility.** A row not converged after `max_terms`
terms returns `NA` with a warning naming it.

## See also

[`quad_vec`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)

## Examples

``` r
# geometric series against the closed form
r <- c(0.1, 0.5, 0.9, 0.99)
series_vec(function(k, i) r[i]^k, n = 4)
#> [1]   1.111111   2.000000  10.000000 100.000000
1 / (1 - r)
#> [1]   1.111111   2.000000  10.000000 100.000000

# Poisson masses sum to one for every rate at once
lam <- c(0.5, 4, 60)
series_vec(function(k, i) dpois(k, lam[i]), n = 3)
#> [1] 1 1 1
```
