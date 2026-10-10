# Stencil Offsets for a Derivative Order

Sizes a stencil from the derivative order and the requested accuracy,
and returns the offsets at which to evaluate: the symmetric ones used
away from a boundary, and the one-sided ones used where a symmetric
stencil would leave the domain. Pass them to
[`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
for the weights and to
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
to apply the whole thing.

## Usage

``` r
fd_offsets(order, accuracy = 2L)
```

## Arguments

- order:

  The derivative order \\d\\. Not validated here, though
  [`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
  rejects anything but a non-negative whole number when the offsets
  reach it.

- accuracy:

  The order of the error term, a positive integer, `2` by default. A
  value below one signals an error. The section on odd accuracy
  describes the treatment of an odd value.

## Value

A list of four components:

- `reach`:

  integer, the half-width \\r\\, at least 1.

- `central`:

  integer vector `-r:r`, the symmetric stencil.

- `forward`:

  integer vector `0:(2r)`, for the lower boundary.

- `backward`:

  integer vector `(-2r):0`, for the upper one.

All three offset vectors have the same length, \\2r + 1\\. At an odd
derivative order the central stencil has a zero weight at the origin,
which
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
does not evaluate, so the central side then uses \\2r\\ evaluations and
each one-sided side \\2r + 1\\.

## The reach

For order \\d\\ and accuracy \\a\\ the half-width is

\$\$r = \Bigl\lceil \tfrac{d + a}{2} \Bigr\rceil - 1,\$\$

floored at one, giving \\2r + 1\\ nodes. At the default accuracy of two
that is the three-point stencil for the first and second derivatives and
the five-point one for the third and fourth. At accuracy four it is the
five-point stencils for the first and second.

## Odd accuracy

A central stencil is symmetric, so the odd powers cancel from its error
expansion and the accuracy of the central stencil is always even. An odd
accuracy is therefore replaced by an even neighbor, and the parity of
\\d + a\\ decides which one. The observed orders on \\\exp\\, obtained
by halving the step, are:

|            |     |     |     |     |
|------------|-----|-----|-----|-----|
| **order**  | 1   | 2   | 3   | 4   |
| accuracy 2 | 2   | 2   | 2   | 2   |
| accuracy 3 | 2   | 4   | 2   | 4   |
| accuracy 4 | 4   | 4   | 4   | 4   |

At an odd order the accuracy is rounded down and the stencil keeps its
size; at an even order it is rounded up and the stencil gains two nodes.
An even accuracy is used as given.

## See also

[`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
for the weights at these offsets,
[`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
for the step to pair with them,
[`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
for all three assembled.

## Examples

``` r
# Three points for a first or second derivative, five for a third or fourth.
fd_offsets(1)$central
#> [1] -1  0  1
fd_offsets(4)$central
#> [1] -2 -1  0  1  2

# Accuracy four buys two more nodes for a first derivative.
fd_offsets(1, accuracy = 4)$central
#> [1] -2 -1  0  1  2

# The one-sided sets are the same size, so a boundary costs no more.
str(fd_offsets(2))
#> List of 4
#>  $ reach   : int 1
#>  $ central : int [1:3] -1 0 1
#>  $ forward : int [1:3] 0 1 2
#>  $ backward: int [1:3] -2 -1 0

# Accuracy must be positive.
try(fd_offsets(2, accuracy = 0))
#> Error : 'accuracy' must be a positive integer.
```
