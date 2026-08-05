# Owen's T Function

Computes \\T(h, a) = \dfrac{1}{2\pi}\displaystyle\int_0^{a}
\dfrac{e^{-h^2(1 + x^2)/2}}{1 + x^2}\\dx\\, which is what the skew
normal distribution function is written in.

## Usage

``` r
owen_t(h, a)
```

## Arguments

- h:

  A numeric vector.

- a:

  A numeric vector, recycled against `h`.

## Value

A numeric vector.

## Details

The integrand is bounded and smooth over a finite range, so quadrature
evaluates it to near machine precision; every element of the input goes
into one batched call of
[`quad_vec`](https://statmodels7.github.io/numericals7/reference/quad_vec.md),
one row per element. Two identities keep the extremes exact rather than
quadrature-bound: \\T(h, a) = -T(h, -a)\\, and \\T(h, \infty) =
\tfrac{1}{2}\Phi(-\|h\|)\\.

## References

Owen, D. B. (1956). Tables for computing bivariate normal probabilities.
*Annals of Mathematical Statistics* 27, 1075-1090.

## Examples

``` r
owen_t(0, 1) - atan(1) / (2 * pi)
#> [1] -3.469447e-16
```
