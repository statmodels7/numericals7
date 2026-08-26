# Derivatives of the Logarithm of the Modified Bessel Function I

Computes \\\log I\_\nu(x)\\ together with its first four derivatives
with respect to the *argument*. The first is the ratio identity \\(\log
I\_\nu)' = \nu/x + I\_{\nu+1}/I\_\nu\\, with the ratio formed as the
exponential of a difference of logarithms and therefore finite wherever
the logarithms are. The higher orders follow from the modified Bessel
equation and cost no further Bessel evaluations, so the whole table is
the price of two.

## Usage

``` r
log_bessel_i_derivs(x, nu)
```

## Arguments

- x:

  A numeric vector of arguments, positive, recycled against `nu`.

- nu:

  A numeric vector of orders, non-negative.

## Value

A named list of five numeric vectors, each of the recycled length of `x`
and `nu`: `l`, the value \\\log I\_\nu(x)\\, and `d1` to `d4`, its
derivatives in the argument.

## Details

Derivatives with respect to the *order* have no elementary form and are
not provided. A caller needing one differences
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
with a single stencil from
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md).

## References

Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
computation of the logarithm of modified Bessel functions on GPUs.
*Proceedings of the 38th ACM International Conference on Supercomputing
(ICS '24)*. arXiv:2409.08729.

## See also

[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
[`log_bessel_k_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k_derivs.md)

## Examples

``` r
str(log_bessel_i_derivs(2, 0.5))
#> List of 5
#>  $ l : num 0.716
#>  $ d1: num 0.787
#>  $ d2: num 0.049
#>  $ d3: num 0.0327
#>  $ d4: num -0.151

# At half-integer order the Bessel functions are elementary, so the first
# derivative has a closed form to check against.
log_bessel_i_derivs(2, 0.5)$d1 - (1 / tanh(2) - 1 / (2 * 2))
#> [1] -1.110223e-16

# Elsewhere, against one stencil on the value itself.
log_bessel_i_derivs(3, 2)$d1 -
  fd_derivative(function(z) log_bessel_i(z, 2), 3, 1, accuracy = 4)
#> [1] 5.262457e-14
```
