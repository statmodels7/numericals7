# Derivatives of the Bessel Ratio

\\A(\kappa)\\ and its first four derivatives, obtained by
differentiating \\A' = 1 - A/\kappa - A^2\\ repeatedly.

## Usage

``` r
bessel_i_ratio_derivs(kappa)
```

## Arguments

- kappa:

  A positive numeric vector.

## Value

A named list with `A` and its derivatives `d1` to `d4`.

## Details

Each order is written in the orders below it, so the whole table costs
the two Bessel evaluations of
[`bessel_i_ratio`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
and nothing more. The first identity follows from \\I_0' = I_1\\ and
\\I_1' = I_0 - I_1/\kappa\\; the alternative, a Bessel function of
higher order per derivative, costs more and is less accurate at large
\\\kappa\\, where the functions themselves overflow and only their ratio
does not. \\A'\\ is the variance of a cosine and therefore positive.

## See also

[`bessel_i_ratio`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md),
[`bessel_i_ratio_inverse`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)

## Examples

``` r
bessel_i_ratio_derivs(2)$d1
#> [1] 0.1642232
```
