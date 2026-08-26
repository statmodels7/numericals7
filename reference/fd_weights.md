# Finite-Difference Weights for Any Stencil

Solves for the weights that combine function values at the given offsets
into an estimate of a derivative. The weights are the unique solution of
the Vandermonde system that makes the stencil exact on polynomials, so a
stencil on `n` nodes reproduces the derivative of every polynomial of
degree `n - 1` or less without error. They are returned for a **unit
step**: divide by `h^order` before using them at spacing `h`, or call
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md),
which assembles the whole expression.

## Usage

``` r
fd_weights(offsets, order)
```

## Arguments

- offsets:

  A numeric vector of distinct offsets, in units of the step. They need
  not be sorted, symmetric, or whole numbers: `0:2` gives a one-sided
  stencil for use at a boundary, and `c(-1, 0, 3)` an uneven one. At
  least `order + 1` of them are needed. Duplicated offsets throw an
  error, the Vandermonde system being singular.

- order:

  The derivative order: a single non-negative whole number, **strictly
  smaller than `length(offsets)`**. Order `0` is legal and returns
  interpolation weights at the origin. Anything else throws: an order at
  or above the node count with a message naming both numbers, and a
  negative, fractional, missing or non-scalar order with a message
  naming the requirement.

## Value

A numeric vector of weights for a unit step, one per offset and in the
order the offsets were given. For `order >= 1` the entries sum to zero.

## The system solved

Write \\s_1, \dots, s_n\\ for the offsets and \\d\\ for the order. The
weights \\w\\ solve the \\n\\ moment conditions

\$\$\sum\_{j=1}^{n} w_j\\ s_j^{\\m} \\=\\ d!\\\[\\m = d\\\], \qquad m =
0, 1, \dots, n-1,\$\$

one equation per polynomial degree the stencil must reproduce. The
matrix is Vandermonde in the offsets, so it is nonsingular exactly when
they are distinct and the weights exist and are unique.

Two properties of the answer are cheap to check against a returned
vector. The weights sum to zero for every \\d \ge 1\\, this being the
\\m = 0\\ condition. On a symmetric stencil they are antisymmetric for
odd \\d\\ and symmetric for even \\d\\.

## Applying them

At a step \\h\\ the derivative estimate is

\$\$f^{(d)}(x) \\\approx\\ h^{-d} \sum\_{j=1}^{n} w_j\\ f(x + s_j
h).\$\$

The factor \\h^{-d}\\ belongs to the caller. The weights carry no step,
so one solve serves every step size.
[`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
supplies a step balanced for the order and
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
puts the three together.

## Accuracy

Because the stencil is exact to degree \\n-1\\, the leading error is the
first term it cannot reproduce, giving an estimate accurate to
\\O(h^{\\n-d})\\ with constant

\$\$\frac{1}{n!}\sum\_{j=1}^{n} w_j\\ s_j^{\\n}.\$\$

Five nodes therefore give a fourth-order first derivative and a
second-order fourth derivative.
[`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
sizes a stencil from the order and the accuracy asked of it.

## One stencil, never nested

Reaching a high order by composing low-order differences multiplies the
error of each stage into the next, and a fourth derivative built from
four nested first differences is noise. Every numerical fallback in the
toolkit takes one stencil of the order it wants, and these weights are
what it takes.

## See also

[`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
for the offsets to pass in,
[`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
for the step to divide by, and
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
for all three assembled into a derivative.

## Examples

``` r
# Three nodes give the classical central differences.
fd_weights(c(-1, 0, 1), 1)                 # -1/2, 0, 1/2
#> [1] -0.5  0.0  0.5
fd_weights(c(-1, 0, 1), 2)                 #    1,  -2,  1
#> [1]  1 -2  1

# Five nodes, first derivative: the familiar (1, -8, 0, 8, -1)/12.
fd_weights(-2:2, 1) * 12
#> [1]  1 -8  0  8 -1

# The weights are for a unit step, so the caller divides by h^order.
h <- 1e-3
w <- fd_weights(c(-1, 0, 1), 2)
sum(w * exp(1 + c(-1, 0, 1) * h)) / h^2    # exp'' at 1
#> [1] 2.718282
exp(1)
#> [1] 2.718282

# Exactness is what they are solved for. On five nodes a first derivative
# reproduces every polynomial up to degree four and fails at degree five.
s <- -2:2
w <- fd_weights(s, 1)
vapply(0:5, function(m) sum(w * s^m), numeric(1))
#> [1]  0  1  0  0  0 -4

# That first failure is the error constant: -4/5! = -1/30.
sum(w * s^5) / factorial(5)
#> [1] -0.03333333

# Which predicts the error to three figures at a step of 0.05.
sum(w * sin(0.7 + s * 0.05)) / 0.05 - cos(0.7)
#> [1] -1.592947e-07
-1 / 30 * 0.05^4 * cos(0.7)
#> [1] -1.593421e-07

# Offsets need not straddle the point; this one-sided stencil is what a
# fallback uses where a symmetric one would leave the domain.
fd_weights(0:2, 1)
#> [1] -1.5  2.0 -0.5
```
