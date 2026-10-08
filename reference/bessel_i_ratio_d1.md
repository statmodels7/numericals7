# Derivatives of the Bessel Ratio

Compute the derivatives of \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\ in
\\\kappa\\, one function per order: `bessel_i_ratio_d1()` returns
\\A'\\, `bessel_i_ratio_d2()` \\A''\\, and so on to the fourth. Each
computes its own order and the orders below it that the computation
needs, and nothing above.

## Usage

``` r
bessel_i_ratio_d1(kappa, threads = 1L)

bessel_i_ratio_d2(kappa, threads = 1L)

bessel_i_ratio_d3(kappa, threads = 1L)

bessel_i_ratio_d4(kappa, threads = 1L)
```

## Arguments

- kappa:

  A numeric vector of concentrations, non-negative and of any size. Zero
  returns 0, the limit, and `Inf` returns 1. A negative value returns
  `NaN`.

- threads:

  The number of threads, a positive whole number.

## Value

A numeric vector the length of `kappa`. `bessel_i_ratio_d1()` is
strictly positive at every finite concentration and tends to 1/2 at
zero.

## Details

The derivatives follow from the identity \\A' = 1 - A/\kappa - A^2\\, a
consequence of \\I_0' = I_1\\ and \\I_1' = I_0 - I_1/\kappa\\,
differentiated repeatedly: \$\$A'' = -\frac{A'}{\kappa} +
\frac{A}{\kappa^2} - 2AA', \qquad A''' = -\frac{A''}{\kappa} +
\frac{2A'}{\kappa^2} - \frac{2A}{\kappa^3} - 2(A')^2 - 2AA'',\$\$ and
the fourth in the same pattern. In double precision this identity
cancels at both ends of the range: at small \\\kappa\\ its terms are of
order \\\kappa^{-n}\\ at derivative \\n\\ while the result is of order
one or \\\kappa\\, and at large \\\kappa\\ the derivative \\A'\\ is of
order \\\kappa^{-2}\\ against terms of order one. The functions
therefore use three regimes, as
[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
does: below \\\kappa = 0.5\\ the power series of \\A\\ at zero
differentiated term by term; from 0.5 to 30 the continued fraction for
\\A\\ and the identity above, both in double-double arithmetic (about 32
digits); from 30 the asymptotic series in \\1/\kappa\\ differentiated
term by term. Each order is accurate to the last bits over the whole
range.

\\A'\\ is the variance of \\\cos(\Theta - \mu)\\ under a von Mises
distribution and is therefore positive; the higher derivatives are its
cumulants of order three to five.

## See also

[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
for the value,
[`bessel_i_ratio_inverse_d1()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse_d1.md)
for the derivatives of the inverse map.

## Examples

``` r
bessel_i_ratio_d1(c(0.1, 1, 100))
#> [1] 4.981302e-01 3.543460e-01 5.025383e-05

# The first derivative satisfies the identity the others are built from.
k <- 2
bessel_i_ratio_d1(k) - (1 - bessel_i_ratio(k) / k - bessel_i_ratio(k)^2)
#> [1] 0

# At a small concentration the series keeps every digit: A''' tends to
# -3/8 and A'''' to 5 kappa / 4.
c(bessel_i_ratio_d3(1e-6), bessel_i_ratio_d4(1e-6) / 1e-6)
#> [1] -0.375  1.250
```
