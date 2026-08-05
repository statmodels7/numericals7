# Package index

## Jets

A jet is a value carried together with every partial derivative up to
fourth order. Seed the variables, write the map in ordinary R, and the
arithmetic operators and mathematical functions propagate every
derivative exactly: no chain rule is transcribed and none can be
mistranscribed.

- [`jet_layout()`](https://statmodels7.github.io/numericals7/reference/jet_layout.md)
  : The Bookkeeping a Jet Needs
- [`jet_var()`](https://statmodels7.github.io/numericals7/reference/jet_var.md)
  : A Jet for One Variable
- [`jet_const()`](https://statmodels7.github.io/numericals7/reference/jet_const.md)
  : A Constant Jet

## Enumerations

The combinatorial objects a higher-order chain rule sums over, each in
one copy so that an enumeration cannot disagree with itself.

- [`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
  : The Index Tuples of a Given Width
- [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
  : The Set Partitions of the First n Integers
- [`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md)
  : Every Way to Write an Integer as an Ordered Sum

## Finite differences

One stencil library: weights from the Vandermonde system, offsets sized
from order and accuracy, a balanced step, and an applicator that never
composes lower-order differences.

- [`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
  : Finite-Difference Weights for Any Stencil
- [`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
  : Stencil Offsets for a Derivative Order
- [`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
  : A Step Size for One Stencil
- [`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
  : One Stencil, Applied
