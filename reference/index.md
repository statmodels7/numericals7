# Package index

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

## Quadrature and series

Integrals and sums vectorized over the parameters: one matrix evaluation
for many parameter values, adaptivity batched by row, and NA with a
named row where the requested accuracy cannot be reached.

- [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  : Integrate One Function at Many Parameter Values
- [`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
  : Sum One Series at Many Parameter Values
- [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  : The Gauss-Kronrod 7-15 Pair

## Special functions

- [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
  : Logarithm of the Modified Bessel Function of the First Kind
- [`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
  : Logarithm of the Modified Bessel Function of the Second Kind
- [`log_bessel_i_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md)
  : Derivatives of the Logarithm of the Modified Bessel Function I
- [`log_bessel_k_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k_derivs.md)
  : Derivatives of the Logarithm of the Modified Bessel Function K
- [`mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.md)
  : The Mills Ratio and Its Derivative
- [`owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.md)
  : Owen's T Function
- [`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
  : The Ratio of Modified Bessel Functions
- [`bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md)
  : Derivatives of the Bessel Ratio
- [`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
  : The Inverse of the Bessel Ratio
