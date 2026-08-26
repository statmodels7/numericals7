# Stencil Offsets for a Derivative Order

Sizes a stencil from the derivative order and the accuracy asked of it,
and returns the offsets to evaluate at: the symmetric ones used away
from a boundary, and the one-sided ones used where a symmetric stencil
would leave the domain. Pass them to
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

  The order of the error term, a positive integer, `2` by default. Zero
  or below throws. See above for what an odd value does.

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

All three offset vectors have the same length, \\2r + 1\\, so the three
sides cost the same number of evaluations.

## The reach

For order \\d\\ and accuracy \\a\\ the half-width is

\$\$r = \Bigl\lceil \tfrac{d + a}{2} \Bigr\rceil - 1,\$\$

floored at one, giving \\2r + 1\\ nodes. At the default accuracy of two
that is the three-point stencil for the first and second derivatives and
the five-point one for the third and fourth. At accuracy four it is the
five-point stencils for the first and second. Those are every stencil
the toolkit's packages had written out by hand, from one formula.

## An odd accuracy is rounded, and which way depends on the order

A central stencil is symmetric, so the odd powers cancel from its error
expansion and the accuracy it delivers is always even. An odd request is
therefore served by an even neighbor, and the parity of \\d + a\\
decides which one. Measured on \\\exp\\ by halving the step:

|            |     |     |     |     |
|------------|-----|-----|-----|-----|
| **order**  | 1   | 2   | 3   | 4   |
| accuracy 2 | 2   | 2   | 2   | 2   |
| accuracy 3 | 2   | 4   | 2   | 4   |
| accuracy 4 | 4   | 4   | 4   | 4   |

At an odd order the request rounds down and costs nothing extra; at an
even order it rounds up and buys two more nodes. Ask for an even
accuracy and the question does not arise.

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
