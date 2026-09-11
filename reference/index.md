# Package index

## Enumerations

The combinatorial objects a higher-order chain rule sums over, each in
one copy so that an enumeration cannot disagree with itself.

- [`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
  : Multi-Indices for Derivatives of a Given Order
- [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
  : Set Partitions of the First n Integers
- [`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md)
  : Weak Compositions of an Integer

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

## Threads

The thread policy of the toolkit, one object passed as an argument from
the fit entry points down to the compiled kernels. The result does not
depend on the count, bit for bit.

- [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  [`print(`*`<n_threads>`*`)`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  : How Many Threads and Worker Processes a Fit May Use
- [`thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.md)
  : Read the Thread Count Off a Policy
- [`worker_count()`](https://statmodels7.github.io/numericals7/reference/worker_count.md)
  : Read the Worker Count Off a Policy
- [`local_threads()`](https://statmodels7.github.io/numericals7/reference/local_threads.md)
  : Apply a Thread Policy for the Duration of One Fit

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
- [`bessel_i_ratios()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratios.md)
  : The Sequence of Modified Bessel Ratios
- [`bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md)
  : Derivatives of the Bessel Ratio
- [`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
  : The Inverse of the Bessel Ratio

## Smoothers of the absolute value

A smooth replacement for the absolute value with its derivatives up to
order five, from which the smooth sign, step and hinge follow by
composition. A break-point term smoothed this way has a true Jacobian
block, and the validator checks a smoother written outside the package.

- [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
  : S7 Class for Smoothers of the Absolute Value
- [`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md)
  : The Probit Smoother of the Absolute Value
- [`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)
  : The Hyperbolic Smoother of the Absolute Value
- [`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md)
  : The Quintic Smoother of the Absolute Value
- [`smoother_deriv()`](https://statmodels7.github.io/numericals7/reference/smoother_deriv.md)
  : Evaluate a Smoother or One of Its Derivatives
- [`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)
  : Resolve a Smoother's Width from a Spacing
- [`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
  : The Smallest Width a Consumer May Use
- [`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md)
  : Check a Smoother of the Absolute Value Numerically
