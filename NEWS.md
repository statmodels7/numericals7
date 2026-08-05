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
