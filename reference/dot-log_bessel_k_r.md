# The R Twin of the Compiled log K Kernel

The same branches and formulas as the compiled kernel behind
[`log_bessel_k`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md),
in vectorized R. The compiled route measured 2.9x faster on the mixed
workload; the twin is kept as the independent reference the tests
compare against.

## Usage

``` r
.log_bessel_k_r(x, nu)
```

## Arguments

- x:

  A non-negative numeric vector.

- nu:

  A numeric vector; the function is even in the order.

## Value

A numeric vector.
