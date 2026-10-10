# Sum One Series at Many Parameter Values

Computes \\\sum\_{k \ge k_0} t(k; \theta_i)\\ for every row \\i\\ at
once. Terms are evaluated in blocks as one matrix, and a row retires as
soon as it converges, so rows that have converged are not evaluated
again while slower rows continue. This is the discrete counterpart of
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md),
written for the same reason: the toolkit's sums are one series at many
parameter values.

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

  The term function, a function of `k` and `i` as described above.

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

  The number of terms after which a row that has not converged returns
  `NA`, `100000` by default. Terms are evaluated in whole blocks, so a
  row may consume up to the first multiple of `block` at or above this
  number; the terms skipped after a probe and the probes themselves are
  not counted.

- block:

  How many terms are evaluated per pass, a whole number of at least 4,
  `64` by default. Its last two quarters give the decay ratio of the
  tail estimate, so a very small block makes the estimate noisy.

## Value

A numeric vector of sums, one per row, of length `n`. `NA` in any row
that did not converge within `max_terms`, with a warning that lists the
first eight of those rows.

## The term function

`term(k, i)` receives two integer vectors of equal length and returns
the terms elementwise. `k` is the summation index and `i` gives the
parameter set to which each term belongs. A Poisson mass at a rate
vector `lam` is

    term <- function(k, i) dpois(k, lam[i])

## Convergence, row by row

Write \\b = \max(\mathrm{atol}, \mathrm{rtol}\\\lvert S_i \rvert)\\ for
the budget of row \\i\\, and \\A_3\\ and \\A_4\\ for the sums of the
absolute terms over the third and the fourth quarter of the last block.
The ratio \\q = A_4/A_3\\ measures the decay of the terms at the end of
the block, and the neglected tail is estimated by continuing that decay
geometrically,

\$\$R = 2\\A_4\\\frac{q}{1 - q}.\$\$

A row retires when \\q \< 1\\ (or \\A_4 = 0\\), \\R \le b\\, and the
last block contributed at most \\b\\.

The estimate bounds the tail of a series whose terms decay geometrically
or faster. For terms that decay like \\k^{-p}\\, geometric continuation
underestimates the tail by a factor that tends to \\(p - 1)/p\\, and the
factor 2 covers every \\p \ge 2\\; a slower decay, \\1/k^2\\ included,
does not reach its budget within the default `max_terms` and returns
`NA`. The sum of \\1/k^3\\ from one retires after about 75000 terms with
a relative error below the default `rtol`.

Rising terms give \\q \> 1\\, so a block that lies before the mode of a
hump-shaped term does not retire its row, however small its terms are:
`dpois(0:63, 300)` sums to 3.8e-62 and is still rising.

## Rows of zeros

A block whose terms are all zero retires a row only after the row has
seen a nonzero term, so the zeros past a finite support end the sum. A
row that has seen only zeros is probed ahead on a geometric grid of
relative spacing 1/128, up to the index \\2^{31}\\, in one call to
`term`. If every probe is zero the row is taken as identically zero and
returns 0, as a structurally null component of an expected Hessian does;
otherwise the summation resumes from the last zero probe, so leading
terms that underflow to zero (`dpois(k, 1e4)` for \\k\\ below about
8000) do not end the sum. The resumption assumes that the terms rise
before the first nonzero one, and a head narrower than the spacing of
the grid (a Poisson mass at a rate above about \\10^8\\) falls between
the probes and gives 0.

## Rows that do not converge

A row still unconverged after `max_terms` terms returns `NA`, and one
warning lists the first eight such rows (followed by an ellipsis if
there are more).

## See also

[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md),
the continuous counterpart.

## Examples

``` r
# Four geometric series against the closed form.
r <- c(0.1, 0.5, 0.9, 0.99)
series_vec(function(k, i) r[i]^k, n = 4) - 1 / (1 - r)
#> [1]  0.000000e+00  0.000000e+00 -2.131628e-14 -4.614648e-09

# Poisson masses sum to one at every rate, in one call.
lam <- c(0.5, 4, 60)
series_vec(function(k, i) dpois(k, lam[i]), n = 3)
#> [1] 1 1 1

# The first 64 terms of a Poisson at rate 300 sum to almost nothing and are
# still rising, so the row does not retire there.
sum(dpois(0:63, 300))
#> [1] 3.75793e-62
series_vec(function(k, i) dpois(k, 300), n = 1)
#> [1] 1

# At rate 1e4 the first 8000 or so terms underflow to zero, and the sum
# continues past them.
series_vec(function(k, i) dpois(k, 1e4), n = 1)
#> [1] 1

# Series indexed from one, against their closed forms: a factorial decay
# and a polynomial one.
series_vec(function(k, i) 1 / factorial(k), n = 1, from = 1L) - (exp(1) - 1)
#> [1] 2.220446e-16
series_vec(function(k, i) 1 / k^4, n = 1, from = 1L) - pi^4 / 90
#> [1] -7.228129e-11

# A divergent series is refused, not estimated.
suppressWarnings(series_vec(function(k, i) 1 / (k + 1), n = 1,
                            max_terms = 500L))
#> [1] NA
```
