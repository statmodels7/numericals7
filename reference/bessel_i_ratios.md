# The Sequence of Modified Bessel Ratios

Computes \\I_j(\kappa)/I_0(\kappa)\\ for \\j = 1, \dots, m\\ by Miller's
backward recurrence, vectorized over \\\kappa\\. The von Mises
distribution function is a series in these ratios, and the recurrence
returns all \\m\\ of them in one pass.

## Usage

``` r
bessel_i_ratios(kappa, m)
```

## Arguments

- kappa:

  A numeric vector of concentrations, positive.

- m:

  How many ratios to return, a positive whole number. It sets the number
  of columns and, with `kappa`, the starting index of the recurrence.

## Value

A numeric matrix of `length(kappa)` rows and `m` columns. Entry \\(i,
j)\\ is \\I_j(\kappa_i)/I_0(\kappa_i)\\, decreasing along a row.

## The backward recurrence

The three-term recurrence \\I\_{j-1} - I\_{j+1} = (2j/\kappa) I_j\\ has
two solutions, one growing and one decaying. Run upwards it should
follow the decaying one and instead follows rounding error into the
growing one, so it is unstable. Run downwards the roles swap and it is
stable, which is Miller's algorithm. The ratios \\r_j = I_j/I\_{j-1}\\
satisfy \\r_j = 1/(2j/\kappa + r\_{j+1})\\, started from \\r\_{n_0+1} =
0\\ at an index far enough above both \\m\\ and \\\kappa\\; the result
is their running product, and the normalization by \\I_0\\ requires no
extra work because the product starts there.

## Cost

The recurrence loop runs over the series index and is vectorized over
the data, so the number of steps depends on \\m\\ and on the largest
\\\kappa\\, and not on the length of the vector. For \\m \> 1\\ the
running products are then formed row by row. A series over these ratios
therefore costs less than one quadrature per observation, and the von
Mises distribution function is evaluated this way.

[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
is the first of these ratios and switches to an asymptotic series from
\\\kappa = 30\\. No such branch is implemented here: the recurrence
needs a starting index above \\\kappa\\, so its cost grows with the
concentration, and at such concentrations a series in these ratios does
not converge in a useful number of terms.

## See also

[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
for the first ratio alone,
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
for the functions themselves.

## Examples

``` r
# Four ratios at two concentrations, one row each.
r <- bessel_i_ratios(c(1, 5), 4)
round(r, 6)
#>          [,1]     [,2]     [,3]     [,4]
#> [1,] 0.446390 0.107220 0.017510 0.002162
#> [2,] 0.893383 0.642647 0.379266 0.187528

# They agree with the scaled Bessel functions to the last bit.
r[2L, ] - besselI(5, 1:4, TRUE) / besselI(5, 0, TRUE)
#> [1] 0 0 0 0

# And decrease along a row: a higher order is a smaller ratio.
all(diff(r[2L, ]) < 0)
#> [1] TRUE
```
