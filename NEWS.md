# numericals7 0.3.0

* Quadrature and series vectorized over the parameters. `quad_vec()`
  integrates one function at many parameter values by matrix evaluation --
  the nodes of every panel of every row go into the integrand in a single
  call per refinement pass -- with the Gauss-Kronrod 7-15 pair supplying an
  error estimate from the same function values, adaptivity batched by row so
  that one hard row cannot serialize the others, rational maps for infinite
  endpoints, and refusal over plausibility: a row that cannot reach the
  requested accuracy returns NA with a warning naming it. `series_vec()`
  does the same for series, in blocks with a per-row convergence mask and a
  tail guard that sees past a block straddling the mode of a hump-shaped
  term. `gauss_kronrod15()` exposes the rule, pinned in the tests by its
  defining property: degree-13 exactness for the embedded Gauss rule,
  degree-22 for the Kronrod extension, and a weight corrupted by 5% fails
  both the moment conditions and the gamma normalization.

* Measured motivation: the von Mises variance through the per-theta fallback
  costs 25 ms per parameter value and scales linearly. One matrix pass over
  hundreds of rows is what regression models, where parameters vary by
  observation, need.

# numericals7 0.2.0

* The stencil library, unifying the three finite-difference implementations
  the toolkit carried: `fd_weights()` solves the Vandermonde system for any
  offsets and order (the construction basis7 had), `fd_offsets()` sizes a
  stencil from the order and the accuracy asked of it, `fd_step()` balances
  truncation against rounding and keeps every node inside a bounded domain,
  and `fd_derivative()` applies one stencil -- never a composition of
  lower-order differences, since each numerical differentiation multiplies
  the error of the one before it. At the default accuracy these reproduce
  linkfunctions7's four central stencils and basis7's shapes exactly; at
  accuracy four they reproduce distributions7's five-point `fd5_first` and
  `fd5_second`. The policy around a stencil -- which order to fall back from,
  what a trustworthy reference is -- deliberately stays with the callers.

# numericals7 0.1.0

* First release: the numerical layer of the statmodels7 toolkit, collecting
  what its packages had been rewriting privately. The census that motivated
  it found `set_partitions()` written twice (distributions7 and parameters7,
  independently) and finite-difference machinery written three times.

* Jets, moved here from parameters7. A jet carries a value together with
  every partial derivative to fourth order and propagates them exactly
  through sums, products and a vocabulary of smooth functions (`exp`, `log`,
  `log1p`, `expm1`, an arbitrary power, `sqrt`, `gamma`, `lgamma`, `digamma`,
  `trigamma`, `sin`, `cos`). `Ops` and `Math` dispatch on the class, so a map
  written as `mu / gamma(1 + 1 / sigma)` differentiates itself. Comparison
  operators and the non-smooth functions are refused: a branch taken on a jet
  would keep one side's derivatives and report them as the whole
  expression's.

* The enumerations a higher-order chain rule rests on, in one copy each:
  `tuple_indices()` (the multi-indexes of a derivative, diagonal-first at
  order two because that ordering is a contract), `set_partitions()` (the
  Bell recursion, blocks indexing positions so that a repeated variable
  carries its multiplicity), and `compositions()` (the weak compositions,
  which are the support of a multinomial).
