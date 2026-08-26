# The Sequence of Modified Bessel Ratios

Computes \\I_j(\kappa)/I_0(\kappa)\\ for \\j = 1, \dots, m\\ by Miller's
backward recurrence, vectorized over \\\kappa\\. A series in these
ratios is what a von Mises distribution function costs, so getting all
\\m\\ of them for the price of one matters.

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

## Why the recurrence runs backwards

The three-term recurrence \\I\_{j-1} - I\_{j+1} = (2j/\kappa) I_j\\ has
two solutions, one growing and one decaying. Run upwards it should
follow the decaying one and instead follows rounding error into the
growing one, so it is unstable. Run downwards the roles swap and it is
stable, which is Miller's algorithm. The ratios \\r_j = I_j/I\_{j-1}\\
satisfy \\r_j = 1/(2j/\kappa + r\_{j+1})\\, started from \\r\_{n_0+1} =
0\\ at an index far enough above both \\m\\ and \\\kappa\\; the answer
is their running product, and the normalization by \\I_0\\ is free
because the product starts there.

## Cost

The loop runs over the series index, never over the data, so a vector of
\\\kappa\\ costs the same number of vectorized steps as a single value.
A series over these ratios therefore costs less than a quadrature per
observation, which is why the von Mises distribution function stopped
being one.

[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
is the first of them and carries an asymptotic branch past \\\kappa =
10^4\\, where the scaled Bessel functions underflow. There is no such
branch here, and none is wanted: the recurrence needs a starting index
above \\\kappa\\, so its cost grows with the concentration, and a caller
that far out is already past the point where a series in these ratios
converges in any useful number of terms.

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
