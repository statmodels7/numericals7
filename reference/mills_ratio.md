# The Mills Ratio and Its Derivative

Returns \\R(t) = \phi(t)/\Phi(t)\\, often called the inverse Mills
ratio, and its derivative \\R'(t) = -R(t)\\t + R(t)\\\\. Every
derivative of a skew normal log-density is built from these two
quantities.

## Usage

``` r
mills_ratio(t)
```

## Arguments

- t:

  A numeric vector. Large finite values of either sign are handled on
  the log scale and are not special-cased. An infinite value gives `NaN`
  in `dr`, and `-Inf` gives `NaN` in `r` as well.

## Value

A list of two numeric vectors, each the length of `t`:

- `r`:

  the ratio \\R(t) = \phi(t)/\Phi(t)\\, positive and decreasing,
  asymptotic to \\-t\\ as \\t \to -\infty\\. In double precision it
  underflows to zero above about \\t = 38.6\\.

- `dr`:

  its derivative \\R'(t)\\, which lies in \\(-1, 0)\\.

## Details

The ratio is formed on the log scale. Written directly it is \\0/0\\ for
\\t\\ below about \\-38\\, where both the density and the distribution
function underflow, while the ratio itself is finite there and close to
\\-t\\. The identity for \\R'\\ follows from differentiating the
quotient and using \\\phi'(t) = -t\phi(t)\\.

## See also

[`owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.md),
[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md),
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)

## Examples

``` r
mills_ratio(c(-5, 0, 3))$r
#> [1] 5.186503967 0.797884561 0.004437839

# The point of the log-scale form. Written directly the ratio is 0/0 well
# inside the range a skew normal reaches, while the true value is finite
# and close to -t.
dnorm(-40) / pnorm(-40)
#> [1] NaN
mills_ratio(-400)$r + (-400)
#> [1] 0.002499971

# The derivative is the stated identity, exactly.
m <- mills_ratio(c(-5, 0, 3))
max(abs(m$dr - (-m$r * (c(-5, 0, 3) + m$r))))
#> [1] 0
```
