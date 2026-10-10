# The Inverse of the Bessel Ratio

Computes \\\kappa = A^{-1}(\rho)\\, the concentration whose ratio
\\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\ equals \\\rho\\. It is the map
used by a von Mises method of moments, which converts an observed mean
resultant length into the concentration that produced it.

## Usage

``` r
bessel_i_ratio_inverse(rho, threads = 1L)
```

## Arguments

- rho:

  A numeric vector of mean resultant lengths, strictly inside \\(0,
  1)\\. Anything outside, the endpoints included, returns `NA` without a
  warning (at \\\rho = 0\\ the inverse is the limit 0, and at \\\rho =
  1\\ it diverges).

- threads:

  The number of threads, a positive whole number.

## Value

A numeric vector of concentrations the length of `rho`, `NA` wherever
`rho` left \\(0, 1)\\.

## Details

\\A\\ has no elementary inverse, so \\\kappa\\ is found by Newton's
method in compiled code, started from a piecewise approximation to the
inverse (a short series for \\\rho \< 0.53\\ and rational forms above).
\\A\\ is increasing and concave, so after the first step every iterate
lies on the left of the root and rises to it; a step that would fall
below \\2\rho\\ is replaced by \\2\rho\\, which is still on the left
since \\A(\kappa) \< \kappa/2\\. The iteration ends when a step is no
larger than four times the machine epsilon times the iterate, or when an
iterate fails to rise above the previous one, and in any case after 200
iterations.

For \\\rho \ge 1/2\\ the residual \\A(\kappa) - \rho\\ is formed as
\\(1 - \rho) - (1 - A(\kappa))\\, where \\1 - \rho\\ is exact and \\1 -
A\\ keeps its relative accuracy: it is computed from \\A\\ in
double-double arithmetic below \\\kappa = 30\\, and from its own
asymptotic series above, without forming \\A\\. Near \\\rho = 1\\ the
inverse behaves as \\\kappa \approx 1/(2(1 - \rho))\\, so its relative
error is that of \\1 - \rho\\: the result is the concentration whose
ratio is the double `rho`, to the last bits, at any concentration.

## See also

[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
for the forward map,
[`bessel_i_ratio_inverse_d1()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse_d1.md)
for the derivatives of the inverse.

## Examples

``` r
# The round trip closes to machine precision across the range.
rho <- c(0.1, 0.5, 0.99)
bessel_i_ratio(bessel_i_ratio_inverse(rho)) - rho
#> [1] 0 0 0

# Outside the open unit interval, the endpoints included, the result is NA.
bessel_i_ratio_inverse(c(0, 0.5, 1))
#> [1]      NA 1.15932      NA
```
