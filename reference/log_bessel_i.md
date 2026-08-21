# Logarithm of the Modified Bessel Function of the First Kind

\\\log I\_\nu(x)\\ for \\x \ge 0\\ and \\\nu \ge 0\\, computed with
every intermediate quantity on the log scale, so the result is finite
and accurate wherever \\\log I\_\nu(x)\\ itself is representable –
including where the unscaled function overflows (from about \\x = 700\\)
and where the exponentially scaled one underflows or loses its precision
(large order with a small argument, and arguments beyond about
\\10^5\\).

## Usage

``` r
log_bessel_i(x, nu, threads = 1L)
```

## Arguments

- x:

  A numeric vector of non-negative arguments, recycled against `nu`.

- nu:

  A numeric vector of non-negative orders.

- threads:

  How many threads the kernel may use, a plain integer, as
  [`thread_count`](https://statmodels7.github.io/numericals7/reference/thread_count.md)
  reads it off an
  [`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  policy. Every branch of the kernel is this package's own arithmetic,
  so element \\i\\ is computed and written by one thread and the result
  is bit-identical at any count; below an internal threshold the
  sequential path is taken whatever the count says.
  [`log_bessel_k`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
  takes no such argument: its hybrid branch calls R's own scaled
  `besselK`, which can raise a warning, and a warning from a worker
  thread ends the session.

## Value

A numeric vector of \\\log I\_\nu(x)\\ values; `-Inf` at \\x = 0\\ with
\\\nu \> 0\\, `0` at \\x = 0\\ with \\\nu = 0\\, `NA` for arguments
outside the domain.

## Details

The algorithm follows Plesner, Sørensen and Hauberg (2024): the
truncated power series with its terms carried as logarithms and summed
by the log-sum-exp anchored at the largest term, the large-argument
asymptotic expansion, and the large-order uniform asymptotic expansion
in the \\u_k(t)\\ polynomials, selected by the input-range table of that
paper. Two switching guards are tightened relative to the paper,
measured on the Wronskian identity \\I\_\nu K\_{\nu+1} + I\_{\nu+1}
K\_\nu = 1/x\\: the large-argument expansions are used only where their
truncation error is below \\10^{-11}\\, the uniform expansion taking
over otherwise.

## References

Plesner, A., Sørensen, H. H. B., and Hauberg, S. (2024). Accurate
computation of the logarithm of modified Bessel functions on GPUs.
*Proceedings of the 38th ACM International Conference on Supercomputing
(ICS '24)*. arXiv:2409.08729.

Olver, F. W. J., et al. (2024). *NIST Digital Library of Mathematical
Functions*, chapter 10. https://dlmf.nist.gov/.

## See also

[`log_bessel_k`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md),
[`log_bessel_i_derivs`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md),
[`bessel_i_ratio`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)

## Examples

``` r
log_bessel_i(c(1, 100, 1e6), 2)      # far past the unscaled overflow
#> [1]     -1.996957     96.759632 999992.173304
log_bessel_i(0.001, 500)             # far past the scaled underflow
#> [1] -6411.782
```
