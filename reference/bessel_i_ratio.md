# The Ratio of Modified Bessel Functions

\\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\, a strictly increasing bijection
from \\(0, \infty)\\ onto \\(0, 1)\\. For a von Mises distribution it is
the mean resultant length, the expected cosine of the deviation from the
mean direction.

## Usage

``` r
bessel_i_ratio(kappa)
```

## Arguments

- kappa:

  A positive numeric vector.

## Value

A numeric vector in \\(0, 1)\\.

## Details

Both Bessel functions are taken exponentially scaled, so the factor
\\e^{\kappa}\\ they share cancels in the ratio where the unscaled
functions would overflow, from about \\\kappa = 700\\. The scaled
functions themselves underflow to an exact zero between \\10^5\\ and
\\10^6\\, so past \\\kappa = 10^4\\ the ratio is taken from its
asymptotic expansion \\1 - 1/(2\kappa) - 1/(8\kappa^2) -
1/(8\kappa^3)\\, whose next term is already below the resolution of a
double at the switch point; the result is therefore finite and accurate
for an argument of any size.

## See also

[`bessel_i_ratio_derivs`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md),
[`bessel_i_ratio_inverse`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)

## Examples

``` r
bessel_i_ratio(c(0.5, 2, 1000))
#> [1] 0.2424996 0.6977747 0.9994999
```
