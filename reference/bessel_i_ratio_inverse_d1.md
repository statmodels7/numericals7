# Derivatives of the Inverse Bessel Ratio

Compute the derivatives of the inverse map \\\kappa(\rho) =
A^{-1}(\rho)\\ in \\\rho\\, one function per order, evaluated at \\\rho
= A(\kappa)\\. They take the concentration rather than \\\rho\\, so a
caller that has already inverted \\\rho\\ does not invert it again.

## Usage

``` r
bessel_i_ratio_inverse_d1(kappa, threads = 1L)

bessel_i_ratio_inverse_d2(kappa, threads = 1L)

bessel_i_ratio_inverse_d3(kappa, threads = 1L)

bessel_i_ratio_inverse_d4(kappa, threads = 1L)
```

## Arguments

- kappa:

  A numeric vector of concentrations, non-negative. A negative value
  returns `NaN`.

- threads:

  The number of threads, a positive whole number.

## Value

A numeric vector the length of `kappa`: the derivative of the inverse
map at \\\rho = A(\kappa)\\.

## Details

The derivatives come from the inverse function rule on the derivatives
of \\A\\
([`bessel_i_ratio_d1()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_d1.md)
and the following orders): \$\$\kappa' = \frac{1}{A'}, \qquad \kappa'' =
-\frac{A''}{(A')^3}, \qquad \kappa''' = \frac{3(A'')^2 -
A'A'''}{(A')^5},\$\$ \$\$\kappa'''' = \frac{-15(A'')^3 + 10A'A''A''' -
(A')^2A''''}{(A')^7}.\$\$ The derivative of order \\n\\ needs \\A'\\ to
\\A^{(n)}\\, which are computed together in one evaluation. Since \\A'\\
decays like \\1/(2\kappa^2)\\, its powers in the denominators would
underflow at large concentrations, so the formulas are evaluated on \\w
= 1/A'\\ and the ratios \\r_j = A^{(j)}/A'\\,

\$\$\kappa'' = -r_2 w^2, \qquad \kappa''' = (3r_2^2 - r_3)\\w^3, \qquad
\kappa'''' = (-15r_2^3 + 10r_2r_3 - r_4)\\w^4,\$\$

with the ratios taken from the asymptotic series without its powers of
\\1/\kappa\\ from \\\kappa = 30\\. Each derivative is finite and
accurate until the derivative itself overflows (near \\\kappa =
10^{60}\\ for the fourth).

## See also

[`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
for the inverse itself.

## Examples

``` r
k <- bessel_i_ratio_inverse(0.7)

# The first derivative is the reciprocal of A', the inverse function rule.
bessel_i_ratio_inverse_d1(k) - 1 / bessel_i_ratio_d1(k)
#> [1] 0

# The second against a central difference of the first in rho.
h <- 1e-5
c(bessel_i_ratio_inverse_d2(k),
  (bessel_i_ratio_inverse_d1(bessel_i_ratio_inverse(0.7 + h)) -
     bessel_i_ratio_inverse_d1(bessel_i_ratio_inverse(0.7 - h))) / (2 * h))
#> [1] 31.61031 31.61031
```
