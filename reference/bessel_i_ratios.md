# The Sequence of Modified Bessel Ratios

\\I_j(\kappa)/I_0(\kappa)\\ for \\j = 1, \dots, m\\, by the backward
recurrence, vectorized over \\\kappa\\.

## Usage

``` r
bessel_i_ratios(kappa, m)
```

## Arguments

- kappa:

  A positive numeric vector.

- m:

  How many ratios to return.

## Value

A `length(kappa)` by `m` matrix.

## Details

The three-term recurrence \\I\_{j-1} - I\_{j+1} = (2j/\kappa) I_j\\ is
unstable run upwards – the recessive solution it should follow is
swamped by the dominant one – and stable run downwards, which is
Miller's algorithm. The ratios \\r_j = I_j/I\_{j-1}\\ satisfy \\r_j =
1/(2j/\kappa + r\_{j+1})\\, started from \\r\_{n_0+1} = 0\\ at an index
far enough above both \\m\\ and \\\kappa\\; the answer is their running
product, and the normalization by \\I_0\\ is free because the product
starts there.

The loop runs over the series index and not over the data, so a vector
of \\\kappa\\ costs the same number of vectorized steps as a single
value. That is what makes a series over these ratios cheaper than a
quadrature per observation.

[`bessel_i_ratio`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
is the first of them, \\I_1/I_0\\, and carries an asymptotic branch for
an argument past \\10^4\\ where the scaled Bessel functions underflow.
There is no such branch here: the recurrence needs a starting index
above \\\kappa\\, so the cost grows with it, and a caller working at a
concentration that large is past the point where a series over these
ratios converges in any useful number of terms.

## See also

[`bessel_i_ratio`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md),
[`log_bessel_i`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)

## Examples

``` r
r <- bessel_i_ratios(c(1, 5), 4)
r[2L, ]
#> [1] 0.8933831 0.6426467 0.3792657 0.1875279
besselI(5, 1:4, expon.scaled = TRUE) / besselI(5, 0, expon.scaled = TRUE)
#> [1] 0.8933831 0.6426467 0.3792657 0.1875279
```
