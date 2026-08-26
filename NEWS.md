# numericals7 0.12.0

* `set_partitions()` returns integer blocks whatever storage `n` arrives
  in. It coerced nothing, so the blocks inherited the caller's mode and
  `set_partitions(4)` gave doubles where `set_partitions(4L)` gave
  integers, and where both sibling enumerations, `tuple_indices()` and
  `compositions()`, give integers however they are called. The consequence
  was that `identical(sort(unlist(p)), 1:4)` was `FALSE` for every one of
  the fifteen partitions of four, and is `TRUE` for all fifteen now.
  Measured across the four call sites in `distributions7` and
  `parameters7` that consume the enumeration, which are the wrapper
  derivatives at orders three and four, the Bartlett expected Hessian, the
  reparametrized chain rule and `sum_struct()`'s log-determinant
  expansion: every computed value is unchanged bit for bit. The coercion
  is applied where `n` enters a block and not to `n` at the top, so that a
  zero, negative or fractional argument still recurses until the stack
  overflows, which is what the page documents.

# numericals7 0.11.0

* `fd_weights()` rejects an `order` that is not a single non-negative whole
  number. A fractional order used to fall through the existing checks and
  return the weights of the truncated order scaled by `factorial(order)` — a
  vector shaped like a stencil that satisfies no moment condition, and no
  warning. A negative one returned all `NaN` behind a warning about `NaN`s.
  `fd_derivative()` is covered by the same check, always going through
  `fd_weights()`. The orders that mean something are unaffected, order zero
  included.

* `bessel_i_ratios(kappa, m)` gives \eqn{I_j(\kappa)/I_0(\kappa)} for
  \eqn{j = 1, \dots, m} by Miller's backward recurrence, vectorized over
  the argument: the loop runs over the series index and not over the data, so
  a vector of concentrations costs the same number of steps as one. That is
  what makes a series over these ratios cheaper than a quadrature per
  observation, and it is why `distributions7`'s von Mises distribution
  function stopped being one. Checked against R's own `besselI` from
  \eqn{\kappa = 0.01} to 500 and orders to 200, agreeing to 1e-15 wherever
  the reference itself has not underflowed.

* `bessel_i_ratio()` is the first of them and keeps its asymptotic branch for
  an argument past \eqn{10^4}, where the scaled Bessel functions underflow.
  The sequence has none: the recurrence needs a starting index above the
  argument, so its cost grows with it, and a caller working there is past the
  point where a series over these ratios converges in any useful number of
  terms.

# numericals7 0.10.0

* `log_bessel_i()` takes a `threads` count and runs its elementwise loop
  over that many. Every branch of the kernel is this package's own
  arithmetic -- a series or a uniform asymptotic expansion, with no call
  into Rmath -- so element `i` is computed and written by one thread and the
  result is bit-identical at any count. Measured at 20000 points: 167 ms
  sequential against 52 at eight threads in the dearest branch (small
  argument), 3.6x, and 5.9x where the asymptotic branch runs. It is the
  package's first parallel kernel; until now nothing here was threaded
  although the policy object lives here.

* `log_bessel_k()` deliberately takes no count: its hybrid branch calls R's
  own scaled `besselK`, which can raise a warning, and a warning from a
  worker thread ends the session. It is also the cheap one of the pair, 8 to
  10 ms over the same 20000 points, so what the restriction costs is small.
  An argument that would be read by nobody is worse than one that is absent.

* `src/n7_par.h` carries the driver, in the shape the toolkit's other two
  have: the worker's loop is noinline, the sequential branch runs through
  the worker, the calling thread's floating-point environment is installed
  before the chunk, and the count is passed to `parallelFor()` rather than
  left to `RCPP_PARALLEL_NUM_THREADS`. `LinkingTo` gains RcppParallel and
  the namespace imports it, without which the package's own DLL does not
  find TBB at load.

# numericals7 0.9.3

* The guarantee on the `n_threads()` page states what now holds rather than
  what was true of the drivers at 0.9.2. The qualification for a kernel
  reading the platform's math routines is gone: those differences were
  measured again on 2026-08-21 and were neither deterministic nor
  unbindable, and a worker that installs the calling thread's
  floating-point environment reproduces the sequential value exactly. What
  remains qualified is the one case where a thread count changes which
  implementation runs -- a threaded kernel in place of a BLAS call.

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
