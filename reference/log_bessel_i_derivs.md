# Derivatives of the Logarithm of the Modified Bessel Function I

\\\log I\_\nu(x)\\ together with its first four derivatives with respect
to the *argument*. The first derivative is the ratio identity \\(\log
I\_\nu)' = \nu/x + I\_{\nu+1}/I\_\nu\\ with the ratio formed as the
exponential of a difference of logarithms, which is finite wherever the
logs are; the higher orders follow from the modified Bessel equation and
cost no further Bessel evaluations. Derivatives with respect to the
order have no elementary form and are not provided.

## Usage

``` r
log_bessel_i_derivs(x, nu)
```

## Arguments

- x:

  A numeric vector of positive arguments, recycled against `nu`.

- nu:

  A numeric vector of non-negative orders.

## Value

A list with the value `l` and the derivatives `d1` to `d4` in the
argument.

## References

Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
computation of the logarithm of modified Bessel functions on GPUs.
*Proceedings of the 38th ACM International Conference on Supercomputing
(ICS '24)*. arXiv:2409.08729.

## See also

[`log_bessel_i`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
[`log_bessel_k_derivs`](https://statmodels7.github.io/numericals7/reference/log_bessel_k_derivs.md)

## Examples

``` r
log_bessel_i_derivs(2, 0.5)$d1   # equals coth(2) - 1/(2*2) exactly
#> [1] 0.7873147
```
