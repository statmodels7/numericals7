# Owen's T Function

Computes \\T(h, a) = \dfrac{1}{2\pi}\displaystyle\int_0^{a}
\dfrac{e^{-h^2(1 + x^2)/2}}{1 + x^2}\\\mathrm{d}x\\, the function in
which the skew normal distribution function is written. \\T(h, a)\\ is
the probability that a pair of independent standard normal variables
falls in the wedge below the line of slope \\a\\ beyond \\h\\, so it is
bounded by \\1/4\\ and is odd in \\a\\.

## Usage

``` r
owen_t(h, a)
```

## Arguments

- h:

  A numeric vector, the offset, of any sign and size.

- a:

  A numeric vector of slopes, recycled against `h`, of any sign and
  size, `Inf` included.

## Value

A numeric vector of the recycled length of `h` and `a`, bounded in
\\\[-1/4, 1/4\]\\, with `NA` where either argument is missing.

## The integral

\\T\\ is even in \\h\\ and odd in \\a\\, so the computation runs on
\\\lvert h\rvert\\ and \\\lvert a\rvert\\. For \\\lvert a\rvert \le 1\\
the factor \\e^{-h^2/2}\\ is taken out of the integral,

\$\$T(h, a) = \frac{e^{-h^2/2}}{2\pi} \int_0^{a} \frac{e^{-h^2
x^2/2}}{1 + x^2}\\\mathrm{d}x,\$\$

so the integrand equals one at zero whatever \\h\\ is, and
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
evaluates the integral to a relative tolerance of \\10^{-13}\\ with no
absolute tolerance, which keeps the relative accuracy of a value far
below one. Every such element goes into one batched call, one row per
element, so a whole vector of skew normal probabilities costs a single
quadrature.

## A slope above one

For \\a \> 1\\ the integrand is concentrated within a distance of order
\\1/h\\ of zero, which a quadrature over \\\[0, a\]\\ can miss. The
reflection identity (Owen, 1956)

\$\$T(h, a) = \tfrac{1}{2}\bigl\\Q(h) + Q(ah)\bigr\\ - Q(h)\\Q(ah) -
T(ah, 1/a), \qquad h \ge 0,\$\$

with \\Q = 1 - \Phi\\ computed as an upper tail, replaces it by an
integral with a slope below one. The result is at least a quarter of the
first terms, so the subtraction loses at most two bits.

## Closed forms

\\T(h, 0) = 0\\, \\T(h, \infty) = \tfrac{1}{2}Q(\lvert h\rvert)\\ and
\\T(\pm\infty, a) = 0\\ are set directly, and a missing argument gives
`NA`.

Against values computed to 50 digits on a grid of \\h\\ from 0 to 37 and
\\a\\ from \\10^{-6}\\ to \\10^6\\, the largest relative error is of
order \\10^{-15}\\.

## References

Owen, D. B. (1956). Tables for computing bivariate normal probabilities.
*Annals of Mathematical Statistics* 27, 1075-1090.

## See also

[`mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.md),
[`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md),
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md),
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)

## Examples

``` r
# At h = 0 the integral is elementary: T(0, a) = atan(a) / (2 pi).
a <- c(0.5, 1, 4)
max(abs(owen_t(0, a) - atan(a) / (2 * pi)))
#> [1] 2.775558e-17

# Odd in the second argument, and the infinite case is a normal tail.
owen_t(1, 2) + owen_t(1, -2)
#> [1] 0
owen_t(1.3, Inf) - pnorm(-1.3) / 2
#> [1] 0

# A steep slope, where the integrand sits within 1/h of zero: the skew
# normal distribution function at alpha = 1000 against its limit
# 2 Phi(z) - 1.
z <- 2.5
(pnorm(z) - 2 * owen_t(z, 1000)) - (2 * pnorm(z) - 1)
#> [1] 0
```
