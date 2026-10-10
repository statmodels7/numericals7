# The R Twin of the Compiled log I Kernel

Computes \\\log I\_\nu(x)\\ in vectorized R, through the same seven
branches and the same formulas as the compiled kernel behind
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md).
It serves as the independent reference against which the tests compare
that kernel, so a change to one side that is not made to the other shows
up as a disagreement. No production code calls it; the production route
is
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md).

## Usage

``` r
.log_bessel_i_r(x, nu)
```

## Arguments

- x:

  A numeric vector of arguments, non-negative.

- nu:

  A numeric vector of orders, non-negative.

## Value

A numeric vector of \\\log I\_\nu(x)\\, of length
`max(length(x), length(nu))`. `NA` where either argument is `NA` or
negative; `0` at `x = 0, nu = 0`, since \\I_0(0) = 1\\; and `-Inf` at
`x = 0` for any `nu > 0`. An argument outside the domain does not signal
an error.

## Details

`x` and `nu` are recycled against each other to the longer length.
Branch selection is that of `.lb_branch()`: the ascending series for a
small argument, the large-argument expansion at two truncation depths (3
and 20 terms), and the large-order uniform asymptotic expansion at four
(4, 6, 9 and 13 polynomials), chosen so that every branch is used where
its own error is smallest.

## See also

[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
the compiled kernel that this function mirrors, and
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
for the second-kind counterpart.
