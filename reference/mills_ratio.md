# The Mills Ratio and Its Derivative

Returns \\R(t) = \phi(t)/\Phi(t)\\ and \\R'(t) = -R(t)\\t + R(t)\\\\,
the two quantities every derivative of a skew normal log-density is
built from.

## Usage

``` r
mills_ratio(t)
```

## Arguments

- t:

  A numeric vector.

## Value

A list with `r` and `dr`.

## Details

The ratio is formed on the log scale. Written directly it is \\0/0\\ for
\\t\\ below about \\-38\\, where both the density and the distribution
function underflow, while the ratio itself is finite there and close to
\\-t\\. The identity for \\R'\\ follows from differentiating the
quotient and using \\\phi'(t) = -t\phi(t)\\.

## Examples

``` r
mills_ratio(c(-40, 0, 3))$r
#> [1] 40.024968847  0.797884561  0.004437839
```
