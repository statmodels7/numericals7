# The Student t Distribution Function

Computes \\T\_\nu(q) = P(T \le q)\\ for a Student \\t\\ variable with
\\\nu \> 0\\ degrees of freedom, real-valued, or the upper tail, on the
natural or the log scale. The values are those of
[`stats::pt()`](https://rdrr.io/r/stats/TDist.html), from a compiled
function that never reaches the R API, so that it may be called from a
worker thread. Its C entry point is `n7_pt`, resolved with
`R_GetCCallable("numericals7", "n7_pt")` and taking the arguments
`(x, n, lower_tail, log_p)` that R's C function
[`pt()`](https://rdrr.io/r/stats/TDist.html) takes.

## Usage

``` r
student_t_cdf(q, df, lower.tail = TRUE, log.p = FALSE, threads = 1L)
```

## Arguments

- q:

  A numeric vector of quantiles.

- df:

  A numeric vector of positive degrees of freedom, recycled against `q`.
  `Inf` gives the standard normal.

- lower.tail:

  If `TRUE`, the default, \\P(T \le q)\\; otherwise \\P(T \> q)\\.

- log.p:

  If `TRUE`, the logarithm of the probability.

- threads:

  A single positive integer, how many threads the computation may use.
  Defaults to `1L`.

## Value

A numeric vector of the recycled length of `q` and `df`: the
probability, or its logarithm. A non-positive `df` gives `NaN`.

## Details

The computation is R's own: [`pt()`](https://rdrr.io/r/stats/TDist.html)
from R's mathematical library and the incomplete beta ratio of Algorithm
708 of the ACM Transactions on Mathematical Software (Didonato and
Morris, 1992) as R carries it, copied with the warnings R raises from
them removed. A warning calls into the R API, and a call into the R API
from a worker thread terminates the process.

## References

Didonato, A. R. and Morris, A. H. (1992). Algorithm 708: Significant
digit computation of the incomplete beta function ratios. *ACM
Transactions on Mathematical Software* 18, 360-373.

## See also

[`stats::pt()`](https://rdrr.io/r/stats/TDist.html), whose values these
are, and
[`mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.md).

## Examples

``` r
q <- c(-40, -2, 0, 1.5, 30)
identical(student_t_cdf(q, df = 3.5), pt(q, df = 3.5))
#> [1] TRUE

# the far left tail, on the log scale
student_t_cdf(-1e6, df = 2.5, log.p = TRUE)
#> [1] -34.8682
```
