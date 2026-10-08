# The Inverse of the Bessel Ratio

Computes \\\kappa = A^{-1}(\rho)\\, the concentration whose ratio
\\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\ equals \\\rho\\. This is the map
a von Mises method of moments runs: it turns an observed mean resultant
length back into the concentration that produced it.

## Usage

``` r
bessel_i_ratio_inverse(rho, threads = 1L)
```

## Arguments

- rho:

  A numeric vector of mean resultant lengths, strictly inside \\(0,
  1)\\. Anything outside, the endpoints included, returns `NA` without a
  warning, the inverse having no finite value there.

- threads:

  The number of threads, a positive whole number.

## Value

A numeric vector of concentrations the length of `rho`, `NA` wherever
`rho` left \\(0, 1)\\.

## Details

\\A\\ has no elementary inverse, so \\\kappa\\ is found by Newton's
method from the standard series approximation, in compiled code. \\A\\
is increasing and concave, so after the first step every iterate lies on
the left of the root and rises to it; a step that would fall below
\\2\rho\\ is replaced by it, which is still on the left since
\\A(\kappa) \< \kappa/2\\. The iteration ends where a step is no larger
than the spacing of the doubles at the iterate.

For \\\rho \ge 1/2\\ the residual \\A(\kappa) - \rho\\ is formed as
\\(1 - \rho) - (1 - A(\kappa))\\, where \\1 - \rho\\ is exact and \\1 -
A\\ is computed without forming \\A\\. Near \\\rho = 1\\ the inverse
behaves as \\\kappa \approx 1/(2(1 - \rho))\\, so its relative error is
that of \\1 - \rho\\: the result is the concentration whose ratio is the
double `rho`, to the last bits, at any concentration.

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

# Outside the open unit interval there is no concentration to return.
bessel_i_ratio_inverse(c(0, 0.5, 1))
#> [1]      NA 1.15932      NA
```
