# Derivatives of the Logarithm of the Modified Bessel Function K

\\\log K\_\nu(x)\\ together with its first four derivatives with respect
to the *argument*, from the ratio identity \\(\log K\_\nu)' = \nu/x -
K\_{\nu+1}/K\_\nu\\ and the modified Bessel equation, exactly as in
[`log_bessel_i_derivs`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md).

## Usage

``` r
log_bessel_k_derivs(x, nu)
```

## Arguments

- x:

  A numeric vector of positive arguments, recycled against `nu`.

- nu:

  A numeric vector of orders, any sign.

## Value

A list with the value `l` and the derivatives `d1` to `d4` in the
argument.

## References

Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
computation of the logarithm of modified Bessel functions on GPUs.
*Proceedings of the 38th ACM International Conference on Supercomputing
(ICS '24)*. arXiv:2409.08729.

## See also

[`log_bessel_k`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md),
[`log_bessel_i_derivs`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md)

## Examples

``` r
# d/dx log K_{1/2}(x) is exactly -1/(2x) - 1
log_bessel_k_derivs(2, 0.5)$d1
#> [1] -1.25
```
