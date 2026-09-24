# The Bessel Ratio's Derivatives From Its Asymptotic Series

The first four derivatives of \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\
from \\A \sim \sum\_{n=0}^{20} q_n \kappa^{-n}\\, the quotient of the
asymptotic series of \\I_1\\ and \\I_0\\, differentiated term by term:
\\d^m \kappa^{-n}/d\kappa^m = (-1)^m n(n+1)\cdots(n+m-1)\kappa^{-n-m}\\.

## Usage

``` r
bessel_ratio_series_derivs(kappa)
```

## Arguments

- kappa:

  A numeric vector of concentrations, at least 20.

## Value

A matrix with one row per concentration and four columns, the first to
fourth derivatives.

## Details

The coefficients are dyadic rationals, \\q_n 2^{2n+1}\\ an integer
through \\n = 13\\: 1, -1/2, -1/8, -1/8, -25/128, ... They are the ones
\\\sum_k (-1)^k a_k(\nu)\kappa^{-k}\\ gives for \\\nu = 1\\ divided by
\\\nu = 0\\, with \\a_k(\nu) = \prod\_{j=1}^k (4\nu^2 - (2j-1)^2) /
(k!\\8^k)\\.

## See also

[`bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md)
