# Logarithm of the Modified Bessel Function of the Second Kind

Computes \\\log K\_\nu(x)\\ for \\x \> 0\\ and any real order, carrying
every intermediate quantity on the log scale, so the result is finite
and accurate wherever \\\log K\_\nu(x)\\ itself is representable. \\K\\
is even in its order, so \\\nu\\ enters as \\\|\nu\|\\.

## Usage

``` r
log_bessel_k(x, nu)
```

## Arguments

- x:

  A numeric vector of arguments, recycled against `nu`. Positive; zero
  gives `Inf`, \\K\\ diverging there, and a negative value gives `NA`.

- nu:

  A numeric vector of orders, of any sign, entering as \\\|\nu\|\\.
  Unlike
  [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
  this takes no `threads` argument: its hybrid branch calls R's own
  scaled `besselK`, which can raise a warning, and a warning from a
  worker thread ends the session.

## Value

A numeric vector of \\\log K\_\nu(x)\\, of the recycled length of `x`
and `nu`. `Inf` at \\x = 0\\, and `NA` where `x` is negative or either
argument is missing.

## Details

The large-argument and large-order branches follow Plesner, Sørensen and
Hauberg (2024) exactly as in
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md).
At moderate inputs the exponentially scaled
[`base::besselK()`](https://rdrr.io/r/base/Bessel.html) is
machine-precision exact wherever it does not overflow and one call beats
a quadrature, so it serves that region; the corner where the scaled
value itself overflows (a small argument with the order near the
switching boundary) goes through the integral representation of Rothwell
(2006), evaluated on the log scale over a composite Simpson rule.

## References

Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
computation of the logarithm of modified Bessel functions on GPUs.
*Proceedings of the 38th ACM International Conference on Supercomputing
(ICS '24)*. arXiv:2409.08729.

Rothwell, E. J. (2006). Computation of the logarithm of Bessel functions
of complex argument and fractional order. *Communications in Numerical
Methods in Engineering* 24(3), 237–249.

## See also

[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
[`log_bessel_k_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k_derivs.md)

## Examples

``` r
# It agrees with R's own function wherever that one still evaluates.
x <- c(0.01, 1, 30)
max(abs(log_bessel_k(x, 2.5) - log(besselK(x, 2.5))))
#> [1] 0

# K underflows to zero at a large argument and overflows at a large order,
# and its logarithm is finite at both.
log(besselK(1000, 2.5))
#> [1] -Inf
log_bessel_k(1000, 2.5)
#> [1] -1003.225
log_bessel_k(1, 500)
#> [1] 2950.996

# Even in the order, exactly.
log_bessel_k(1, 2.5) - log_bessel_k(1, -2.5)
#> [1] 0
```
