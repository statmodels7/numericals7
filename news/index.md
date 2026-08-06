# Changelog

## numericals7 0.6.0

- The jets are removed. Every production consumer now carries its
  derivatives as written closed forms: the reparametrized families of
  distributions7 declare their map partials explicitly, and
  parameters7’s autoregressive family propagates its derivative arrays
  through the Levinson-Durbin recursion in compiled code, the product
  rule written out per order. Measured on the Poisson-inverse Gaussian
  kernels, the jet route’s fixed composition overhead was 2x to 36x the
  cost of the written forms.

## numericals7 0.5.0

- The logarithm of the modified Bessel functions, after Plesner,
  Sorensen and Hauberg (ICS 2024, arXiv:2409.08729):
  [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
  and
  [`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
  carry every intermediate quantity on the log scale – the power series
  through the log-sum-exp anchored at its largest term, the
  large-argument and large-order asymptotic expansions selected by the
  paper’s input-range table – and are finite and accurate wherever the
  logarithm itself is representable, where the unscaled functions
  overflow from x = 700 and the scaled ones underflow past 1e5 or lose
  the order. Two switching guards are tightened relative to the paper,
  measured on the Wronskian identity.
  [`log_bessel_i_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md)
  and
  [`log_bessel_k_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k_derivs.md)
  add the first four derivatives in the argument from the ratio identity
  and the Bessel equation. Measured against an Rcpp transcription of the
  same algorithm: the compiled version buys 0.9x to 2.7x, so the package
  stays pure R.

## numericals7 0.4.0

- Special functions, each carrying the overflow discipline learned on
  it:
  [`mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.md)
  on the log scale, finite where density and distribution function both
  underflow;
  [`owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.md)
  through one batched
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  call, with the closed identities at `a = 0` and `a = Inf`;
  [`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
  through the exponentially scaled Bessel functions, with its four
  derivatives from the recurrence
  ([`bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.md))
  and its inverse by root finding with the inverse-function-rule
  derivatives
  ([`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)).

## numericals7 0.3.0

- Quadrature and series vectorized over the parameters.
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  integrates one function at many parameter values by matrix evaluation
  – the nodes of every panel of every row go into the integrand in a
  single call per refinement pass – with the Gauss-Kronrod 7-15 pair
  supplying an error estimate from the same function values, adaptivity
  batched by row so that one hard row cannot serialize the others,
  rational maps for infinite endpoints, and refusal over plausibility: a
  row that cannot reach the requested accuracy returns NA with a warning
  naming it.
  [`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
  does the same for series, in blocks with a per-row convergence mask
  and a tail guard that sees past a block straddling the mode of a
  hump-shaped term.
  [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  exposes the rule, pinned in the tests by its defining property:
  degree-13 exactness for the embedded Gauss rule, degree-22 for the
  Kronrod extension, and a weight corrupted by 5% fails both the moment
  conditions and the gamma normalization.

- Measured motivation: the von Mises variance through the per-theta
  fallback costs 25 ms per parameter value and scales linearly. One
  matrix pass over hundreds of rows is what regression models, where
  parameters vary by observation, need.

## numericals7 0.2.0

- The stencil library, unifying the three finite-difference
  implementations the toolkit carried:
  [`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
  solves the Vandermonde system for any offsets and order (the
  construction basis7 had),
  [`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
  sizes a stencil from the order and the accuracy asked of it,
  [`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
  balances truncation against rounding and keeps every node inside a
  bounded domain, and
  [`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
  applies one stencil – never a composition of lower-order differences,
  since each numerical differentiation multiplies the error of the one
  before it. At the default accuracy these reproduce linkfunctions7’s
  four central stencils and basis7’s shapes exactly; at accuracy four
  they reproduce distributions7’s five-point `fd5_first` and
  `fd5_second`. The policy around a stencil – which order to fall back
  from, what a trustworthy reference is – deliberately stays with the
  callers.

## numericals7 0.1.0

- First release: the numerical layer of the statmodels7 toolkit,
  collecting what its packages had been rewriting privately. The census
  that motivated it found
  [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
  written twice (distributions7 and parameters7, independently) and
  finite-difference machinery written three times.

- Jets, moved here from parameters7. A jet carries a value together with
  every partial derivative to fourth order and propagates them exactly
  through sums, products and a vocabulary of smooth functions (`exp`,
  `log`, `log1p`, `expm1`, an arbitrary power, `sqrt`, `gamma`,
  `lgamma`, `digamma`, `trigamma`, `sin`, `cos`). `Ops` and `Math`
  dispatch on the class, so a map written as `mu / gamma(1 + 1 / sigma)`
  differentiates itself. Comparison operators and the non-smooth
  functions are refused: a branch taken on a jet would keep one side’s
  derivatives and report them as the whole expression’s.

- The enumerations a higher-order chain rule rests on, in one copy each:
  [`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
  (the multi-indexes of a derivative, diagonal-first at order two
  because that ordering is a contract),
  [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
  (the Bell recursion, blocks indexing positions so that a repeated
  variable carries its multiplicity), and
  [`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md)
  (the weak compositions, which are the support of a multinomial).
