# The Gauss-Kronrod 7-15 Pair

The 15-point Kronrod extension of the 7-point Gauss rule on \\\[-1,
1\]\\: one set of nodes, two sets of weights, so that a single
evaluation of the integrand yields both an estimate and, from the
difference of the two rules, an error for it.

## Usage

``` r
gauss_kronrod15()
```

## Value

A list with `nodes`, the Kronrod weights `wk` and the embedded Gauss
weights `wg`, each of length 15.

## Details

The Gauss weights are zero at the eight Kronrod-only nodes, which is
what lets both rules be formed from one matrix of function values. The
constants are the classical ones; the tests pin them by their defining
property rather than by their digits, checking that the 7-point rule
integrates polynomials to degree 13 exactly and the 15-point one to
degree 22.

## See also

[`quad_vec`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)

## Examples

``` r
r <- gauss_kronrod15()
sum(r$wk)   # weights of a rule on [-1, 1] sum to its length
#> [1] 2
sum(r$wg)
#> [1] 2
```
