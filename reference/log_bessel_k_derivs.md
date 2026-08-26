# Derivatives of the Logarithm of the Modified Bessel Function K

Computes \\\log K\_\nu(x)\\ together with its first four derivatives
with respect to the *argument*, from the ratio identity \\(\log K\_\nu)'
= \nu/x - K\_{\nu+1}/K\_\nu\\ and the modified Bessel equation, exactly
as in
[`log_bessel_i_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md).
The sign is the one difference: \\K\\ decreases in its argument where
\\I\\ grows.

## Usage

``` r
log_bessel_k_derivs(x, nu)
```

## Arguments

- x:

  A numeric vector of arguments, positive, recycled against `nu`.

- nu:

  A numeric vector of orders, of any sign, entering as its absolute
  value.

## Value

A named list of five numeric vectors, each of the recycled length of `x`
and `nu`: `l`, the value, and `d1` to `d4`, its derivatives in the
argument. `d1` is negative at every order, `K` decreasing in its
argument.

## Details

Derivatives with respect to the *order* have no elementary form and are
not provided.

## References

Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
computation of the logarithm of modified Bessel functions on GPUs.
*Proceedings of the 38th ACM International Conference on Supercomputing
(ICS '24)*. arXiv:2409.08729.

## See also

[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md),
[`log_bessel_i_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md)

## Examples

``` r
str(log_bessel_k_derivs(2, 0.5))
#> List of 5
#>  $ l : num -2.12
#>  $ d1: num -1.25
#>  $ d2: num 0.125
#>  $ d3: num -0.125
#>  $ d4: num 0.187

# At half-integer order the closed form is elementary: the derivative of
# log K_{1/2}(x) is -1/(2x) - 1.
log_bessel_k_derivs(2, 0.5)$d1 - (-1 / (2 * 2) - 1)
#> [1] -4.440892e-16

# Elsewhere, against one stencil on the value itself.
log_bessel_k_derivs(3, 2)$d1 -
  fd_derivative(function(z) log_bessel_k(z, 2), 3, 1, accuracy = 4)
#> [1] -2.442491e-15
```
