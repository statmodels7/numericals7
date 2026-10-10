# The Ratio of Modified Bessel Functions

Computes \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\, a strictly increasing
bijection from \\(0, \infty)\\ onto \\(0, 1)\\. For a von Mises
distribution it is the mean resultant length, the expected cosine of the
deviation from the mean direction, so it is the map between a
concentration and the moment that a method of moments estimates.

## Usage

``` r
bessel_i_ratio(kappa, threads = 1L)
```

## Arguments

- kappa:

  A numeric vector of concentrations, non-negative and of any size. Zero
  returns 0, the limit, and `Inf` returns 1. A negative value returns
  `NaN`.

- threads:

  The number of threads, a positive whole number.

## Value

A numeric vector the length of `kappa`, in \\\[0, 1\]\\ and increasing
in its argument.

## Details

The ratio is evaluated in compiled code without the Bessel functions
themselves, which overflow from about \\\kappa = 709\\ and,
exponentially scaled, underflow between \\10^5\\ and \\10^6\\. Below
\\\kappa = 0.5\\ it is the power series of \\A\\ at zero; from there to
\\\kappa = 30\\ it is the continued fraction \\A = \kappa/(2 +
\kappa^2/(4 + \kappa^2/(6 + \cdots)))\\, evaluated backwards from the
index \\\lfloor\kappa\rfloor + 20\\; above, it is the asymptotic series
of \\A\\ in \\1/\kappa\\ with 30 terms, the quotient of the asymptotic
series of \\I_1\\ and \\I_0\\. The result is finite and accurate to the
last bits for an argument of any size.

## See also

[`bessel_i_ratio_d1()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_d1.md)
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
# while the asymptotic series carries the answer to any concentration.
suppressWarnings(besselI(1e6, 1, TRUE) / besselI(1e6, 0, TRUE))
#> [1] NaN
bessel_i_ratio(c(1e6, 1e12))
#> [1] 0.9999995 1.0000000
```
