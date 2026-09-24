# Derivatives of the Bessel Ratio

Computes \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\ and its first four
derivatives, by differentiating the identity \\A' = 1 - A/\kappa - A^2\\
repeatedly. A von Mises family needs all five to reach fourth-order
derivatives in its concentration.

## Usage

``` r
bessel_i_ratio_derivs(kappa)
```

## Arguments

- kappa:

  A numeric vector of concentrations, positive.

## Value

A named list of five numeric vectors, each the length of `kappa`: `A`,
the ratio itself, and `d1` to `d4`, its derivatives in \\\kappa\\. `d1`
is strictly positive, being a variance.

## Details

Each order is written in the orders below it, so the whole table costs
the two Bessel evaluations of
[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
and nothing more.

**At a large concentration the recursion cancels**, \\A'\\ being \\1 -
A/\kappa - A^2\\, three terms of order one whose sum is of order
\\\kappa^{-2}\\, and each derivative above it losing a further factor.
Measured against the asymptotic series, the third derivative is out by
3.7e-06 at \\\kappa = 300\\, 3.0e-04 at \\10^3\\ and 0.61 at \\10^4\\.
From \\\kappa = 20\\ the four derivatives are therefore taken from the
series of \\A = I_1/I_0\\ in \\1/\kappa\\, the quotient of the two
asymptotic series of \\I_1\\ and \\I_0\\, differentiated term by term
with 21 terms. The crossover is where the two routes agree best, 5e-13
to 4e-11 over the four orders; the value \\A\\ itself stays
[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)'s.
The first identity follows from \\I_0' = I_1\\ and \\I_1' = I_0 -
I_1/\kappa\\; the alternative, a Bessel function of higher order per
derivative, costs more and is less accurate at large \\\kappa\\, where
the functions themselves overflow and only their ratio does not. \\A'\\
is the variance of a cosine and therefore positive.

## See also

[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
for the value alone,
[`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
for the derivatives of the inverse map.

## Examples

``` r
str(bessel_i_ratio_derivs(2))
#> List of 5
#>  $ A : num 0.698
#>  $ d1: num 0.164
#>  $ d2: num -0.137
#>  $ d3: num 0.113
#>  $ d4: num -0.0437

# The first derivative is the variance of a cosine, so it is positive at
# every concentration and vanishes as the distribution concentrates.
bessel_i_ratio_derivs(c(0.1, 1, 100))$d1
#> [1] 4.981302e-01 3.543460e-01 5.025383e-05

# It satisfies the identity the whole table is built from.
d <- bessel_i_ratio_derivs(2)
d$d1 - (1 - d$A / 2 - d$A^2)
#> [1] 0
```
