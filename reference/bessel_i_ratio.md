# The Ratio of Modified Bessel Functions

Computes \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\, a strictly increasing
bijection from \\(0, \infty)\\ onto \\(0, 1)\\. For a von Mises
distribution it is the mean resultant length, the expected cosine of the
deviation from the mean direction, so it is the map between a
concentration and the moment a method of moments estimates.

## Usage

``` r
bessel_i_ratio(kappa)
```

## Arguments

- kappa:

  A numeric vector of concentrations, positive and of any size. Zero
  returns 0, the limit. A negative value returns `NaN`.

## Value

A numeric vector the length of `kappa`, in \\(0, 1)\\ and increasing in
its argument.

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

[`bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md)
for its derivatives,
[`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
for the map back,
[`bessel_i_ratios()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.md)
for the sequence of higher orders.

## Examples

``` r
bessel_i_ratio(c(0.5, 2, 1000))
#> [1] 0.2424996 0.6977747 0.9994999

# It agrees with the scaled Bessel functions where those still evaluate.
k <- c(0.5, 2, 1e3, 1e4)
max(abs(bessel_i_ratio(k) - besselI(k, 1, TRUE) / besselI(k, 0, TRUE)))
#> [1] 1.110223e-16

# Past that the scaled functions underflow to zero and their ratio is NaN,
# while the asymptotic branch carries the answer to any concentration.
suppressWarnings(besselI(1e6, 1, TRUE) / besselI(1e6, 0, TRUE))
#> [1] NaN
bessel_i_ratio(c(1e6, 1e12))
#> [1] 0.9999995 1.0000000
```
