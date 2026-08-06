# The R Twin of the Compiled log I Kernel

The same branches and formulas as the compiled kernel behind
[`log_bessel_i`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
in vectorized R. The compiled route measured 1.1x faster on a mixed
workload of one million points; the twin is kept as the independent
reference the tests compare against.

## Usage

``` r
.log_bessel_i_r(x, nu)
```

## Arguments

- x:

  A non-negative numeric vector.

- nu:

  A non-negative numeric vector.

## Value

A numeric vector.
