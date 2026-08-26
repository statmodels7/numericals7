# The R Twin of the Compiled log K Kernel

Computes \\\log K\_\nu(x)\\ in vectorized R, through the same branches
and the same formulas as the compiled kernel behind
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md).
It exists as the independent reference the tests compare that kernel
against, so a change to either side that is not a change to both shows
up as a disagreement. Not called on any production path;
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
is.

## Usage

``` r
.log_bessel_k_r(x, nu)
```

## Arguments

- x:

  A numeric vector of arguments, non-negative.

- nu:

  A numeric vector of orders, of any sign.

## Value

A numeric vector of \\\log K\_\nu(x)\\, of length
`max(length(x), length(nu))`. `Inf` at \\x = 0\\, and `NA` where `x` is
negative or either argument is missing.

## Details

`x` and `nu` are recycled against each other to the longer length, and
the order enters as \\\|\nu\|\\, \\K\\ being even in it. The branches
are the large-argument and large-order expansions, R's own scaled
`besselK` in the moderate region, and the Rothwell integral in the
corner where that scaled value overflows.

The compiled route measured 2.9x faster on a mixed workload spanning
every branch, so the twin costs little to keep.

## See also

[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md),
the compiled kernel this mirrors, and
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
for the first-kind counterpart.
