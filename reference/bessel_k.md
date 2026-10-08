# The Modified Bessel Function of the Third Kind

Computes \\K\_\nu(x)\\, or \\e^x K\_\nu(x)\\ when `expon.scaled = TRUE`,
for real orders \\\nu\\. The values are those of
[`base::besselK()`](https://rdrr.io/r/base/Bessel.html), from a compiled
function that never reaches the R API, so that it may be called from a
worker thread. Its C entry point is `n7_bessel_k`, resolved with
`R_GetCCallable("numericals7", "n7_bessel_k")` and taking the arguments
`(x, alpha, expo)` that R's C function `bessel_k()` takes, `expo` being
1 for the function and 2 for the scaled function.

## Usage

``` r
bessel_k(x, nu, expon.scaled = FALSE, threads = 1L)
```

## Arguments

- x:

  A numeric vector of non-negative arguments.

- nu:

  A numeric vector of orders, recycled against `x`; the order enters as
  \\\|\nu\|\\.

- expon.scaled:

  If `TRUE`, \\e^x K\_\nu(x)\\; otherwise, the default, \\K\_\nu(x)\\.

- threads:

  A single positive integer, how many threads the computation may use.
  Defaults to `1L`.

## Value

A numeric vector of the recycled length of `x` and `nu`. A negative `x`
gives `NaN`.

## Details

The computation is R's own: the routine RKBESL of W. J. Cody and L.
Stoltz, after J. B. Campbell's implementation of Temme's algorithm, as R
carries it, copied with the warnings R raises from it removed and with
its work array on the stack rather than in R's memory pool. A warning
calls into the R API, and R's allocator is not safe to call from a
worker thread.

## References

Campbell, J. B. (1980). On Temme's algorithm for the modified Bessel
function of the third kind. *ACM Transactions on Mathematical Software*
6, 581-586.

## See also

[`base::besselK()`](https://rdrr.io/r/base/Bessel.html), whose values
these are, and
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
for the logarithm where the function itself overflows or underflows.

## Examples

``` r
x <- c(0.01, 0.5, 3, 40)
identical(bessel_k(x, nu = 1.5), besselK(x, nu = 1.5))
#> [1] TRUE

# the scaled function at a large argument, where K itself underflows
bessel_k(1000, nu = 2, expon.scaled = TRUE)
#> [1] 0.03970762
```
