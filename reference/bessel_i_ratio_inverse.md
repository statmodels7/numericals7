# The Inverse of the Bessel Ratio

\\\kappa = A^{-1}(\rho)\\ by root finding, together with the four
derivatives of the inverse in \\\rho\\.

## Usage

``` r
bessel_i_ratio_inverse(rho)
```

## Arguments

- rho:

  A numeric vector in \\(0, 1)\\.

## Value

A named list with `kappa` and its derivatives `d1` to `d4` in `rho`.

## Details

\\A\\ has no elementary inverse, so \\\kappa\\ is found by root finding
from a bracket built around the standard series approximation and
widened until it straddles the root; the asymptotic branch of
[`bessel_i_ratio`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
keeps the function evaluable over the whole bracket, however
concentrated. The derivatives come from the inverse function rule on
[`bessel_i_ratio_derivs`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md):
\$\$\kappa' = \dfrac{1}{A'}, \qquad \kappa'' = -\dfrac{A''}{(A')^3},
\qquad \kappa''' = \dfrac{3(A'')^2 - A'A'''}{(A')^5},\$\$ and the fourth
in the same pattern; \\A' \> 0\\ keeps every denominator away from zero
in the interior. A `rho` outside \\(0, 1)\\ returns `NA`.

## See also

[`bessel_i_ratio`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md),
[`bessel_i_ratio_derivs`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md)

## Examples

``` r
bessel_i_ratio(bessel_i_ratio_inverse(0.7)$kappa)
#> [1] 0.7
```
