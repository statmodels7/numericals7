# numericals7 0.9.2

* A second qualification on the `n_threads()` page: a kernel that calls
  into the platform's own math routines per element inherits that
  runtime's per-thread behavior in the last bit (one runtime was measured
  returning one-ulp differences between the main thread and a worker,
  deterministically), so bit-identity across counts is promised for the
  arithmetic a kernel computes itself.

# numericals7 0.9.1

* The `n_threads()` page qualifies its guarantee where a threaded kernel
  replaces a BLAS expression: the count never changes what a kernel
  computes, and the replacement itself is bit-exact against the reference
  BLAS R ships and within the rounding of one dot product against an
  optimized one, whose accumulation order is its own.

# numericals7 0.9.0

* `n_threads()` gains `workers`, the number of R processes the independent
  fits of a cross-validation's folds may use, read by `worker_count()`.
  The same object carries both levels of the toolkit's parallelism, which
  is what the one-argument constructor existed for; the two do not nest (a
  fit inside a worker is sequential by construction), and the result does
  not depend on either count, bit for bit.

# numericals7 0.8.0

* The thread policy of the toolkit lives here, at the root, being the one
  package below every compiled kernel. `n_threads(threads = 1)` constructs
  it, `thread_count()` reads the count off it, and `local_threads()` applies
  it to RcppParallel's process-level setting for one calling frame and
  restores the previous state on exit. The object is passed as an argument
  from the fit entry points (`statmodels7::statmod(threads =)`,
  `distributions7::fit_distrib(threads =)`) down to the kernels: no global
  state, and no package reads a setting that lives in another. The result
  of a fit does not depend on the count, bit for bit, because every
  parallel region in the toolkit decomposes its work over the elements of
  its output and never splits a reduction.

# numericals7 0.7.0

* The log-Bessel kernels are compiled. `log_bessel_i()` and
  `log_bessel_k()` now run scalar C++ loops over the same branches and
  formulas; the vectorized R implementations stay as internal twins
  (`.log_bessel_i_r`, `.log_bessel_k_r`) that a test compares against the
  compiled route on every branch. Measured on a mixed workload of one
  million points: log K 2.9x faster, log I 1.1x. The u_k polynomial table
  is injected into the kernels once at load, so the two routes share the
  table and nothing else.

# numericals7 0.6.0

* The jets are removed. Every production consumer now carries its
  derivatives as written closed forms: the reparametrized families of
  distributions7 declare their map partials explicitly, and parameters7's
  autoregressive family propagates its derivative arrays through the
  Levinson-Durbin recursion in compiled code, the product rule written out
  per order. Measured on the Poisson-inverse Gaussian kernels, the jet
  route's fixed composition overhead was 2x to 36x the cost of the written
  forms.

# numericals7 0.5.0

* The logarithm of the modified Bessel functions, after Plesner, Sorensen
  and Hauberg (ICS 2024, arXiv:2409.08729): `log_bessel_i()` and
  `log_bessel_k()` carry every intermediate quantity on the log scale --
  the power series through the log-sum-exp anchored at its largest term,
  the large-argument and large-order asymptotic expansions selected by the
  paper's input-range table -- and are finite and accurate wherever the
  logarithm itself is representable, where the unscaled functions overflow
  from x = 700 and the scaled ones underflow past 1e5 or lose the order.
  Two switching guards are tightened relative to the paper, measured on
  the Wronskian identity. `log_bessel_i_derivs()` and
  `log_bessel_k_derivs()` add the first four derivatives in the argument
  from the ratio identity and the Bessel equation. Measured against an
  Rcpp transcription of the same algorithm: the compiled version buys
  0.9x to 2.7x, so the package stays pure R.

# numericals7 0.4.0

* Special functions, each carrying the overflow discipline learned on it:
  `mills_ratio()` on the log scale, finite where density and distribution
  function both underflow; `owen_t()` through one batched `quad_vec()` call,
  with the closed identities at `a = 0` and `a = Inf`; `bessel_i_ratio()`
  through the exponentially scaled Bessel functions, with its four
  derivatives from the recurrence (`bessel_i_ratio_derivs()`) and its
  inverse by root finding with the inverse-function-rule derivatives
  (`bessel_i_ratio_inverse()`).

# numericals7 0.3.0

* Quadrature and series vectorized over the parameters. `quad_vec()`
  integrates one function at many parameter values by matrix evaluation --
  the nodes of every panel of every row go into the integrand in a single
  call per refinement pass -- with the Gauss-Kronrod 7-15 pair supplying an
  error estimate from the same function values, adaptivity batched by row so
  that one hard row cannot serialize the others, rational maps for infinite
  endpoints, and rejection over plausibility: a row that cannot reach the
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
  operators and the non-smooth functions are rejected: a branch taken on a jet
  would keep one side's derivatives and report them as the whole
  expression's.

* The enumerations a higher-order chain rule rests on, in one copy each:
  `tuple_indices()` (the multi-indexes of a derivative, diagonal-first at
  order two because that ordering is a contract), `set_partitions()` (the
  Bell recursion, blocks indexing positions so that a repeated variable
  carries its multiplicity), and `compositions()` (the weak compositions,
  which are the support of a multinomial).
