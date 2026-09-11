# The Gauss-Kronrod 7-15 Pair

Returns the nodes and the two sets of weights of the 15-point Kronrod
extension of the 7-point Gauss rule on \\\[-1, 1\]\\. One set of nodes
carries both rules, so a single evaluation of the integrand gives an
estimate and, from the difference of the two rules, an error for it. An
adaptive routine reads that error to decide where to refine, and pays
nothing for it beyond the estimate.

## Usage

``` r
gauss_kronrod15()
```

## Value

A list of three numeric vectors, each of length 15:

- `nodes`:

  the abscissae on \\\[-1, 1\]\\, ascending and symmetric about zero.

- `wk`:

  the Kronrod weights, all positive.

- `wg`:

  the embedded Gauss weights, zero at the eight Kronrod-only nodes.

## Details

The eight Kronrod-only nodes carry a Gauss weight of zero, so both rules
are formed from one matrix of function values by two weighted sums.

The constants are QUADPACK's, carried at the precision it prints them.
An adaptive routine reads the difference of the two rules, and on a
constant integrand that difference is \\\lvert\sum w_k - \sum w_g\rvert
/ 2\\ of the integral on every panel. Transcribed to fifteen decimals,
as they were before version 0.14.0, the Kronrod weights summed to \\2 -
6.0 \times 10^{-15}\\, so every error estimate carried a floor of \\3.4
\times 10^{-15}\\ of the integral that no bisection lowers. At full
precision both sums are 2 within the rounding of the sum. The tests pin
the constants by their defining property, that the 7-point rule
integrates polynomials of degree 13 exactly and the 15-point one degree
22, and pin the two sums to a few units in the last place.

## References

Kronrod, A. S. (1965). *Nodes and Weights of Quadrature Formulas*.
Consultants Bureau, New York.

Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W. and Kahaner, D.
K. (1983). *QUADPACK: A Subroutine Package for Automatic Integration*.
Springer.

## See also

[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md),
which integrates with this pair by default.

## Examples

``` r
r <- gauss_kronrod15()

# The weights of a rule sum to the length of its interval.
c(kronrod = sum(r$wk), gauss = sum(r$wg))
#> kronrod   gauss 
#>       2       2 

# The Gauss rule lives on eight of the fifteen nodes, so both rules come
# from one set of function values.
sum(r$wg == 0)
#> [1] 8

# Its defining property: exact on polynomials to degree 13.
p <- function(x, d) x^d
vapply(c(13, 14), function(d) {
  exact <- if (d %% 2 == 0) 2 / (d + 1) else 0
  sum(r$wg * p(r$nodes, d)) - exact
}, numeric(1))
#> [1]  0.0000000000 -0.0001854659
```
