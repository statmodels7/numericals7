# The Inverse of the Bessel Ratio

Computes \\\kappa = A^{-1}(\rho)\\ by Newton's method, together with the
four derivatives of the inverse in \\\rho\\. This is the map a von Mises
method of moments runs: it turns an observed mean resultant length back
into the concentration that produced it.

## Usage

``` r
bessel_i_ratio_inverse(rho)
```

## Arguments

- rho:

  A numeric vector of mean resultant lengths, strictly inside \\(0,
  1)\\. Anything outside, the endpoints included, returns `NA` without a
  warning, the inverse having no finite value there.

## Value

A named list of five numeric vectors, each the length of `rho`: `kappa`,
the concentration, and `d1` to `d4`, the derivatives of the inverse in
\\\rho\\. `NA` wherever `rho` left \\(0, 1)\\.

## Details

\\A\\ has no elementary inverse, so \\\kappa\\ is found by Newton's
method, vectorized over `rho`, from the standard series approximation.
\\A\\ is increasing and concave, so after the first step every iterate
lies on the left of the root and rises to it; a step that would fall
below \\2\rho\\ is replaced by it, which is still on the left since
\\A(\kappa) \< \kappa/2\\. The iteration ends where a step is no larger
than the spacing of the doubles at the iterate. Near \\\rho = 1\\ the
inverse is ill conditioned, \\d\kappa = 2\kappa^2 d\rho\\, so one unit
in the last place of `rho` moves \\\kappa\\ by a relative \\10^{-3}\\ at
\\\kappa = 5 \times 10^{12}\\: the answer there is one of the
concentrations the forward map sends to `rho`. The derivatives come from
the inverse function rule on
[`bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md):
\$\$\kappa' = \dfrac{1}{A'}, \qquad \kappa'' = -\dfrac{A''}{(A')^3},
\qquad \kappa''' = \dfrac{3(A'')^2 - A'A'''}{(A')^5},\$\$ and the fourth
in the same pattern; \\A' \> 0\\ keeps every denominator away from zero
in the interior. A `rho` outside \\(0, 1)\\ returns `NA`.

## See also

[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
for the forward map,
[`bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md)
for the derivatives it inverts.

## Examples

``` r
# The round trip closes to machine precision across the range.
rho <- c(0.1, 0.5, 0.99)
bessel_i_ratio(bessel_i_ratio_inverse(rho)$kappa) - rho
#> [1] -1.387779e-17 -1.110223e-16 -1.110223e-16

# The first derivative is the reciprocal of A', the inverse function rule.
inv <- bessel_i_ratio_inverse(0.7)
inv$d1 - 1 / bessel_i_ratio_derivs(inv$kappa)$d1
#> [1] 0

# Outside the open unit interval there is no concentration to return.
bessel_i_ratio_inverse(c(0, 0.5, 1))$kappa
#> [1]      NA 1.15932      NA
```
