# Logarithm of the Modified Bessel Function of the Second Kind

\\\log K\_\nu(x)\\ for \\x \> 0\\ and any real order, computed with
every intermediate quantity on the log scale; \\K\\ is even in its
order, so \\\nu\\ enters as \\\|\nu\|\\.

## Usage

``` r
log_bessel_k(x, nu)
```

## Arguments

- x:

  A numeric vector of positive arguments, recycled against `nu`.

- nu:

  A numeric vector of orders, any sign.

## Value

A numeric vector of \\\log K\_\nu(x)\\ values; `Inf` at \\x = 0\\, `NA`
for arguments outside the domain.

## Details

The large-argument and large-order branches follow Plesner, Sørensen and
Hauberg (2024) exactly as in
[`log_bessel_i`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md).
At moderate inputs the exponentially scaled
[`besselK`](https://rdrr.io/r/base/Bessel.html) is machine-precision
exact wherever it does not overflow and one call beats a quadrature, so
it serves that region; the corner where the scaled value itself
overflows (a small argument with the order near the switching boundary)
goes through the integral representation of Rothwell (2006), evaluated
on the log scale over a composite Simpson rule.

## References

Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
computation of the logarithm of modified Bessel functions on GPUs.
*Proceedings of the 38th ACM International Conference on Supercomputing
(ICS '24)*. arXiv:2409.08729.

Rothwell, E. J. (2006). Computation of the logarithm of Bessel functions
of complex argument and fractional order. *Communications in Numerical
Methods in Engineering* 24(3), 237–249.

## See also

[`log_bessel_i`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
[`log_bessel_k_derivs`](https://statmodels7.github.io/numericals7/reference/log_bessel_k_derivs.md)

## Examples

``` r
log_bessel_k(c(0.01, 1, 1000), 2.5)
#> [1]    12.837312     1.171702 -1003.225088
log_bessel_k(1, 500)                 # the unscaled function overflows here
#> [1] 2950.996
```
