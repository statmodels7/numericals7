# Owen's T Function

Computes \\T(h, a) = \dfrac{1}{2\pi}\displaystyle\int_0^{a}
\dfrac{e^{-h^2(1 + x^2)/2}}{1 + x^2}\\\mathrm{d}x\\, the function the
skew normal distribution function is written in. \\T(h, a)\\ is the
probability that a standard bivariate normal pair falls in the wedge
below the line of slope \\a\\ beyond \\h\\, so it is bounded by \\1/4\\
and is odd in \\a\\.

## Usage

``` r
owen_t(h, a)
```

## Arguments

- h:

  A numeric vector, the offset. Any finite value.

- a:

  A numeric vector of slopes, recycled against `h`. May be negative or
  infinite; both are taken by identity.

## Value

A numeric vector of the recycled length of `h` and `a`, bounded in
\\\[-1/4, 1/4\]\\.

## Details

The integrand is bounded and smooth over a finite range, so quadrature
evaluates it to near machine precision. Every element of the input goes
into one batched call of
[`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md),
one row per element, so a whole vector of skew normal probabilities
costs a single quadrature.

Two identities are applied in closed form, so the extremes are exact
where quadrature would merely be accurate: \\T(h, a) = -T(h, -a)\\
handles a negative second argument, and \\T(h, \infty) =
\tfrac{1}{2}\Phi(-\|h\|)\\ handles an infinite one.

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
#> [1] 6.106227e-16

# Odd in the second argument, and the infinite case is a normal tail.
owen_t(1, 2) + owen_t(1, -2)
#> [1] 0
owen_t(1.3, Inf) - pnorm(-1.3) / 2
#> [1] 0
```
