# Changelog

## numericals7 0.22.0

- [`owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.md)
  keeps its relative accuracy at steep slopes and at large `h`. It
  integrated the unscaled integrand over \[0, a\] with an absolute
  tolerance of 1e-15, so for a slope well above one the first panel’s
  nodes could miss the peak of width 1/h near zero (`owen_t(2.5, 1000)`
  returned 7e-28 for 3.1e-3, and the skew normal distribution function
  at alpha = 1000 was off by 6e-3), and a value far below one lost its
  digits to the absolute tolerance (`owen_t(10, 1)` was off by 2e-8
  relative). A slope above one is now reflected onto T(ah, 1/a) through
  Owen’s identity, with the normal tails computed as upper tails, and
  the integral over a slope in (0, 1\] takes the factor exp(-h^2/2) out
  of the integrand and runs with a relative tolerance of 1e-13 and no
  absolute one. Against values computed to 50 digits on 168 points with
  h from 0 to 37 and a from 1e-6 to 1e6, the largest relative error is
  6.5e-16, where 60 points were off by more than 1e-12 before; on 20000
  skew-normal-like arguments the time goes from 0.35 s to 0.09 s and the
  values move by at most 1.1e-16. A missing argument now gives `NA`,
  where it gave 0.

- [`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
  estimates the neglected tail. A row retires when its last block fits
  the budget and the geometric continuation, doubled, of the decay
  between the last two quarters of that block fits it too; rising terms
  give a decay ratio above one and keep the row alive, as the old growth
  test did. Before, a row retired on a small last block and a small last
  term, which bounded the tail only for geometric decay: the sum of
  1/k^3 retired at a relative error of 6e-9 and 0.999^k at 1.4e-9, both
  above the default `rtol` of 1e-10. A block of zeros now retires a row
  only after the row has seen a nonzero term. A row that has seen only
  zeros is probed on a geometric grid of relative spacing 1/128 up to
  the index 2^31 (about 2200 terms in one call): if every probe is zero
  the row returns 0, and otherwise the sum resumes from the last zero
  probe. A head that underflows to zero therefore no longer ends the
  sum: the Poisson masses at rate 1e4 summed to 0 and the Poisson mean
  at rate 2000 was 0, both on the first block, and the masses now sum to
  one within 4e-11 up to rate 1e6. Where the old rule was right the
  result is unchanged bit for bit (2000 Poisson rows of E\[Y^2\]).
  `block` must be a whole number of at least 4.

- [`bessel_i_ratio_inverse_d2()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse_d1.md)
  to `_d4()` stay accurate at large kappa. The inverse function rule
  divided by A’^(2n - 1), and A’ ~ 1/(2 kappa^2) made that power
  subnormal long before the derivative overflows: the fourth derivative
  was Inf or NaN from kappa ~ 1e22, the third from ~1e32 (and 0.4 per
  cent off just below), the second from ~1e52. The rule is now evaluated
  on w = 1/A’ and the ratios A^(j)/A’, taken from the asymptotic series
  without its powers of 1/kappa from kappa = 30, with the powers of w
  applied one at a time. Against mpmath at 120 digits on kappa from 1e-6
  to 1e110 every derivative is within 3.8e-15 relative until its true
  value overflows; the time is unchanged (0.043 s for 1e5 values of the
  fourth).

- [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  accepts an embedded rule of any length, where the 15 nodes of
  [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  were written into the evaluation; a rule whose three vectors differ in
  length signals an error.

- Documentation: the prose of the help pages, the README, the vignette,
  the NEWS file and `_pkgdown.yml` is rewritten in a plainer register,
  and factual errors are corrected. Among them: the observed order of a
  one-sided stencil (one order below the central one only at even
  orders), the number of evaluations of a central stencil at odd orders,
  the warning of
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  and
  [`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
  (it lists the first eight rows that fail), the degree of exactness of
  the Kronrod rule (23), the stopping rule and the starting value of
  [`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md),
  the inherited description of `kappa` on the pages of the derivatives
  of the Bessel ratio, the branch counts of the log-Bessel twins, the
  toolkit’s derivative orders on the pages of
  [`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
  and
  [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md),
  and the Bessel thresholds in the README. One error message is reworded
  ([`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
  for an order at or above the number of nodes).

## numericals7 0.21.1

- The test of
  [`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
  against 150-digit preimages writes its probabilities in hexadecimal.
  On arm64 macOS the decimal literal of 1 - 2^-40 parsed to a
  neighboring double, and near one a unit in the last place moves the
  preimage by 1.2e-4 relative; the function itself was correct on every
  platform.

## numericals7 0.21.0

- [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)’s
  compiled kernel is registered as the C entry point `n7_log_bessel_i`,
  taking the argument and the order, for the compiled code of other
  packages. It never calls the R API.

## numericals7 0.20.0

- [`bessel_k()`](https://statmodels7.github.io/numericals7/reference/bessel_k.md)
  computes the modified Bessel function of the third kind K_nu(x), or
  its exponentially scaled form, with the values of
  [`base::besselK()`](https://rdrr.io/r/base/Bessel.html). The code is
  R’s
  [`bessel_k()`](https://statmodels7.github.io/numericals7/reference/bessel_k.md)
  (the routine RKBESL) copied with the warnings removed and with its
  work array on the stack (allocated with `malloc()` beyond 64 orders),
  so that it never calls the R API and may run on a worker thread. It is
  threaded through `threads` and registered as the C entry point
  `n7_bessel_k`, with the arguments of R’s C function
  [`bessel_k()`](https://statmodels7.github.io/numericals7/reference/bessel_k.md).

## numericals7 0.19.0

- [`student_t_cdf()`](https://statmodels7.github.io/numericals7/reference/student_t_cdf.md)
  computes the Student t distribution function with real degrees of
  freedom, on the natural or the log scale. Its values are those of
  [`stats::pt()`](https://rdrr.io/r/stats/TDist.html): the code is R’s
  [`pt()`](https://rdrr.io/r/stats/TDist.html) and its incomplete beta
  ratio (ACM TOMS Algorithm 708), copied with the warnings removed, so
  that it never calls the R API and may run on a worker thread. It is
  threaded through `threads` and registered as the C entry point
  `n7_pt`, with the arguments of R’s C function
  [`pt()`](https://rdrr.io/r/stats/TDist.html).

- The license is GPL (\>= 2), which the copied R sources require. The
  copyright holders of those sources are listed in `inst/COPYRIGHTS` and
  in `Authors@R`.

## numericals7 0.18.0

- The derivatives of the Bessel ratio A(kappa) = I1(kappa)/I0(kappa) are
  one function per order,
  [`bessel_i_ratio_d1()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_d1.md)
  to
  [`bessel_i_ratio_d4()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_d1.md),
  and `bessel_i_ratio_derivs()`, which always computed all four, is
  removed. Each order computes only itself and the orders below it that
  its formula needs. The inverse follows the same pattern:
  [`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
  returns the concentration alone, and
  [`bessel_i_ratio_inverse_d1()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse_d1.md)
  to
  [`bessel_i_ratio_inverse_d4()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse_d1.md)
  give the derivatives of the inverse in rho, taking the concentration
  so that a caller that has inverted rho does not invert it again. All
  are compiled, threaded through `threads`, and also registered as C
  entry points (`n7_bessel_ratio`, `n7_bessel_ratio_d1`, …,
  `n7_bessel_ratio_upto`, `n7_bessel_ratio_inverse`,
  `n7_bessel_ratio_inverse_d1`, …) for the compiled code of other
  packages.

- The ratio and its derivatives are computed in three regimes. The
  Riccati identity A’ = 1 - A/kappa - A^2, differentiated in double
  precision, cancelled at both ends of the range: at a small kappa its
  terms are of order kappa^-n at order n while the result is of order
  one or kappa. Against 150-digit values the fourth derivative was out
  by 1.1e-4 relative at kappa = 1e-3 and had no correct digit below
  1e-5, and the second by 8e-11 at 1e-3; between 5 and 30 the fourth
  lost up to 1e-8 (4e-10 at kappa = 20). Now the power series of A at
  zero is used below kappa = 0.5, the continued fraction and the
  identity in double-double arithmetic from 0.5 to 30, and the
  asymptotic series in 1/kappa with 30 terms above (21 terms from 20
  before). Over 72 concentrations from 1e-8 to 1e6 the largest relative
  error of A and of every derivative is 6.0e-16. The value no longer
  calls [`besselI()`](https://rdrr.io/r/base/Bessel.html).

- [`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
  resolves the concentration near rho = 1. Newton’s residual is formed
  as (1 - rho) - (1 - A(kappa)) for rho \>= 1/2, where 1 - rho is exact
  and 1 - A keeps its relative accuracy (double-double arithmetic below
  kappa = 30, its own asymptotic series above). The relative error of
  kappa against 150-digit preimages of the same doubles was up to 4e-12
  at kappa = 3e3 and is now within 2.8e-16 over 72 concentrations from
  1e-8 to 1e6; the four derivatives of the inverse are within 3.5e-15.

- Time per value, at kappa drawn from an exponential of mean 3: the
  value 0.04 us, any one derivative 0.4 to 0.5 us (the table of four
  took 0.45 us), the inverse 1.6 us (2.8 before).

## numericals7 0.17.0

- [`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md)
  takes `max_gap`, the largest gap between consecutive distinct values
  of the covariate over the range a break-point may take. A smoother
  that declares an `exact_radius`, which is
  [`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md),
  then has its width raised until the radius is at least 0.55 times that
  gap. The quintic is exactly \|u\| beyond its radius, so a break-point
  in a gap wider than twice the radius has no curvature and is not
  identified, and the median spacing of 0.16.0 does not bound the
  largest gap: among n uniform points the largest gap is about (log n +
  gamma)/log 2 median spacings, 9.5 at n = 400, while the quintic’s
  width is 3.61 median spacings. On 400 uniform points, a `jump()` whose
  fitted break-point lay in a gap of 0.0155 against 2h = 0.0129 had its
  column exactly zero. Other smoothers ignore `max_gap`, and a width
  supplied as `h` is used as it stands.

- [`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md)’s
  page states the consequence of its C^3 continuity under an outer
  criterion: the fourth derivative jumps at +-h and the fifth is a point
  mass there, so the criterion over the hyperparameters is smooth only
  between the values at which an observation crosses psi +- h, and an
  exact outer Hessian leaves out the point masses. On a random
  break-point under REML, the jump of the gradient at a crossing was
  below the resolution of a step of 1e-4 in the hyperparameter.

## numericals7 0.16.0

- [`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md)
  resolves its width at 5/(2 log 2) spacings, about 3.61, instead of
  one. The quintic is exact outside \[-h, h\], so its second derivative,
  which is the curvature of a smoothed break-point, is zero there, and a
  fitted break-point with no observation within h is not identified. At
  one median spacing this happened in 25 to 55 per cent of fits of
  `jump()` and `jseg()` on a uniform covariate, and the break-point’s
  error did not shrink with n (0.037 at n = 200 and at 1000). At the new
  width, which puts five observations inside \[-h, h\] on average for a
  uniform covariate, none of 80 such fits left it unidentified and the
  error is that of the probit smoother (0.0035 against 0.004), with a
  fit error on a true sharp step 10 to 15 per cent above the probit’s.
  `seg()` was not affected. A width supplied as `h` is used as it
  stands, and
  [`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
  is multiplied by the same factor, which makes it more conservative.

## numericals7 0.15.0

- `bessel_i_ratio_derivs()` keeps its digits at a large concentration.
  The four derivatives of A(kappa) = I1(kappa)/I0(kappa) came from
  differentiating A’ = 1 - A/kappa - A^2 repeatedly, and at a large
  kappa that identity is three terms of order one summing to order
  kappa^-2, each higher order losing a further factor. Against the
  asymptotic series of A differentiated term by term, the third
  derivative was out by 3.7e-06 at kappa = 300, 3.0e-04 at 1e3 and 0.61
  at 1e4, and the second by 9.1e-05 at 1e4. From kappa = 20 the four
  derivatives are that series, the quotient of the asymptotic series of
  I1 and I0 in 1/kappa with 21 terms, whose coefficients are dyadic
  rationals; the internal `bessel_ratio_series_derivs()` evaluates it.
  The crossover is where the two routes agree best, 5e-13 to 4e-11 over
  the four orders. The value A itself is unchanged, and so is every
  derivative below kappa = 20. Every von Mises derivative that reads
  these changes above kappa = 20 by the recursion’s error there: 1e-10
  or less at kappa = 100, the whole value past 1e4.

- [`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)
  is Newton’s method vectorized over `rho`, where it was one
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) per element. A is
  increasing and concave, so after the first step every iterate lies on
  the left of the root and rises to it, and a step falling below 2 rho,
  which A(kappa) \< kappa/2 puts on the same side, is replaced by 2 rho.
  Over 4000 values the time goes from 0.38 s to 0.012 s, and a
  `vonmises2_distrib()` smooth at n = 4000, which spent 95 per cent of
  its fit in this function, goes from 54.7 s to 5.0 s with the
  coefficients within 6.7e-16 and the smoothing parameter within 3.7e-16
  relative. The residual of the forward map is at most 4.4e-16, against
  2.2e-13 for [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) at its
  tolerance of eps^0.75. A `rho` below about 5e-11, which raised an
  error because the bracket stopped at 1e-10 above a root of 2 rho, is
  now accepted. Near rho = 1 the two disagree by up to 2e-4 relative, at
  kappa = 5e12, which is within the effect of one unit in the last place
  of `rho` there (1.1e-3), and both send `rho` back to itself exactly.

## numericals7 0.14.1

- The test of the rule’s floor asserts a bound and not the last bit. It
  required `quad_floor(gauss_kronrod15())` to lie below two units in the
  last place, and the floor is 2.2e-16 on x86_64 and 4.4e-16 on the
  arm64 macOS runner, where the computed Kronrod and Gauss sums differ
  by one unit in the last place of 2. It now checks a bound of 1e-15,
  between that floor and the 3.7e-15 of the fifteen-decimal table, and
  that the arm64 value passes it. The pages of
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  and
  [`quad_floor()`](https://statmodels7.github.io/numericals7/reference/quad_floor.md)
  give both values of the floor. No computed value changes.

## numericals7 0.14.0

- [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  carries QUADPACK’s constants at full precision. Transcribed to fifteen
  decimals, the Kronrod weights summed to 2 - 6.0e-15. On a constant
  integrand the difference of the two rules is half the difference of
  their weight sums, so every panel of every integral carried an error
  estimate of 3.4e-15 of the integral that bisection does not lower. At
  full precision both sums are 2 within the rounding of the sum, and the
  Legendre moments that the Kronrod rule must annihilate, to degree 22,
  vanish within 1.2e-16 where the old table left 2.5e-15. Every integral
  computed through the rule changes in the fifteenth digit.

- [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  rejects a relative budget below the rule’s floor when no absolute
  budget is given. The floor is the relative error estimate of the rule
  on a constant integrand, `|sum(wk) - sum(wg)| / 2`, plus the machine
  epsilon: 2.2e-16 for the default rule on x86_64 and 4.4e-16 on the
  arm64 build of R for macOS, and 3.7e-15 for the fifteen-decimal table.
  With `atol = 0` such a row could never converge; the call signals an
  error naming the floor before the integrand is evaluated. A positive
  `atol` can still end the refinement, so `rtol = 0` beside one stays
  valid, which is the condition QUADPACK applies.

- `quad_vec(max_panels = 4096)`, the greatest number of panels one row
  may hold. A row still over its budget once it holds that many returns
  `NA` with the warning, which names both limits. It bounds the case
  that `max_depth` does not: an error spread over the whole row, from
  rounding or from an oscillation faster than the panels resolve, where
  every panel carries a similar share of it, most panels are split at
  every pass and the count grows geometrically long before any panel
  reaches `max_depth`. Without the cap, a constant integrand with a
  budget below the floor exhausted 30 GiB at depth 28. With the cap, a
  row of `sin(1e9 x)` on the unit interval stops after 54 passes with a
  peak 41 MB above the session’s baseline, and 100 such rows in one call
  at 360 MB. The default admits `cos(16000 x)` on the unit interval,
  which needs between 2048 and 4096 panels, and `abs(sin(50 x))` on \[0,
  10\], which needs between 1024 and 2048. The cap is counted per row,
  so a row’s result does not depend on the other rows of the call.

## numericals7 0.13.0

- The smoothers of the absolute value move here from `penalties7`:
  [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md),
  the three shipped smoothers
  ([`smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.md),
  [`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)
  and
  [`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md)),
  [`smoother_deriv()`](https://statmodels7.github.io/numericals7/reference/smoother_deriv.md),
  [`smoother_width()`](https://statmodels7.github.io/numericals7/reference/smoother_width.md),
  [`smoother_width_floor()`](https://statmodels7.github.io/numericals7/reference/smoother_width_floor.md)
  and
  [`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md).
  `penalties7` 0.23.0 no longer exports them. The smoothed families
  planned for `distributions7` need the class, and `distributions7`
  cannot import `penalties7`, while every package of the toolkit can
  import this one. The move is not a re-export, so the class has one
  home and one name.

- The code is unchanged. The R file differs from its last version in
  `penalties7` in two places of its documentation: the `@include` of a
  `penalties7` file is removed, and the link to `check_penalty()` is
  plain text. The test file is copied unchanged and passes, within a
  suite of 593 expectations in 67 blocks with none failing or skipped.
  The class of an object is now
  [`numericals7::abs_smoother`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
  instead of `penalties7::abs_smoother`.

- This is the package’s first S7 class, so `S7` joins `Imports` and the
  namespace imports it, as in `distributions7`, `linkfunctions7` and
  `optimizers7`; none of the 57 names exported by S7 collides with
  `base`, `stats`, `utils` or this package. The print method sits on a
  base generic, and
  [`S7::methods_register()`](https://rconsortium.github.io/S7/reference/methods_register.html)
  in `.onLoad()` registers it for an installed package; this was checked
  on the installed copy, because loading from source registers the
  method in any case.

- A set of 529 quantities computed with the smoothers in `penalties7`
  0.22.1 is reproduced
  [`identical()`](https://rdrr.io/r/base/identical.html) with them here:
  the six derivatives on a 17-point grid at three widths for each
  smoother, the width resolution and its floor, the tables and printed
  text of
  [`check_abs_smoother()`](https://statmodels7.github.io/numericals7/reference/check_abs_smoother.md),
  the blocks and the first and second block derivatives that
  `modelterms7` builds for `seg()`, `jump()` and `jseg()` under each
  smoother, and ten fits with their log-likelihoods, coefficients,
  variances, certificates and summary notes. A change of one unit in the
  last place or one character injected at 18 of the quantities is
  detected each time at the quantity it touches. The 168 comparisons
  that pin the rest of the toolkit’s fits are identical too.

## numericals7 0.12.0

- [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
  returns integer blocks whatever the storage of `n`. Before, the blocks
  inherited the caller’s mode: `set_partitions(4)` gave doubles where
  `set_partitions(4L)` gave integers, while
  [`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
  and
  [`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md)
  give integers however they are called, and
  `identical(sort(unlist(p)), 1:4)` was `FALSE` for all fifteen
  partitions of four. Across the four call sites in `distributions7` and
  `parameters7` that use the enumeration (the wrapper derivatives at
  orders three and four, the Bartlett expected Hessian, the
  reparametrized chain rule and `sum_struct()`’s log-determinant
  expansion) every computed value is unchanged bit for bit. The coercion
  is applied where `n` enters a block and not to `n` at the top, so a
  zero, negative or fractional argument still recurses until R signals a
  stack overflow error, as the page documents.

## numericals7 0.11.0

- [`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
  rejects an `order` that is not a single non-negative whole number. A
  fractional order used to pass the existing checks and return the
  weights of the truncated order scaled by `factorial(order)`, which
  satisfy no moment condition, with no warning. A negative one returned
  all `NaN` with a warning about `NaN`s.
  [`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
  is covered by the same check, since it always calls
  [`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md).
  Valid orders, order zero included, are unaffected.

- `bessel_i_ratios(kappa, m)` gives I_j(kappa)/I_0(kappa) for j = 1, …,
  m by Miller’s backward recurrence, vectorized over the argument: the
  recurrence loop runs over the series index and is vectorized over the
  data. A series over these ratios is therefore cheaper than a
  quadrature per observation, and `distributions7`’s von Mises
  distribution function is now evaluated this way. Checked against R’s
  own `besselI` from kappa = 0.01 to 500 and orders to 200, agreeing to
  1e-15 wherever the reference itself has not underflowed.

- [`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
  is the first of these ratios and keeps its asymptotic branch for an
  argument past 10^4. The sequence has no such branch: the recurrence
  needs a starting index above the argument, so its cost grows with it,
  and at such arguments a series over these ratios does not converge in
  a useful number of terms.

## numericals7 0.10.0

- [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
  takes a `threads` count and runs its elementwise loop over that many
  threads. Every branch of the kernel is this package’s own arithmetic
  (a series or a uniform asymptotic expansion, with no call into Rmath),
  so element `i` is computed and written by one thread and the result is
  bit-identical at any count. At 20000 points: 167 ms sequential against
  52 ms at eight threads in the most expensive branch (small argument),
  3.2x, and 5.9x where the asymptotic branch runs. It is the package’s
  first parallel kernel.

- [`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
  takes no count: its hybrid branch calls R’s own scaled `besselK`,
  which can raise a warning, and a warning from a worker thread ends the
  session. It is also the cheaper of the pair, 8 to 10 ms over the same
  20000 points.

- `src/n7_par.h` carries the driver, in the shape of the toolkit’s other
  two: the worker’s loop is noinline, the sequential branch runs through
  the worker, the calling thread’s floating-point environment is
  installed before the chunk, and the count is passed to `parallelFor()`
  instead of being left to `RCPP_PARALLEL_NUM_THREADS`. `LinkingTo`
  gains RcppParallel and the namespace imports it, without which the
  package’s own DLL does not find TBB at load.

## numericals7 0.9.3

- The guarantee on the
  [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  page states what holds now. The qualification for a kernel reading the
  platform’s math routines is removed: those differences, measured again
  on 2026-08-21, were neither deterministic nor unbindable, and a worker
  that installs the calling thread’s floating-point environment
  reproduces the sequential value exactly. The remaining qualification
  is the one case where a thread count changes which implementation
  runs, a threaded kernel in place of a BLAS call.

## numericals7 0.9.2

- A second qualification on the
  [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  page: a kernel that calls into the platform’s own math routines per
  element inherits that runtime’s per-thread behavior in the last bit
  (one runtime returned one-ulp differences between the main thread and
  a worker, deterministically), so bit-identity across counts is
  promised for the arithmetic a kernel computes itself.

## numericals7 0.9.1

- The
  [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  page qualifies its guarantee where a threaded kernel replaces a BLAS
  expression: the count does not change what a kernel computes, and the
  replacement is bit-exact against the reference BLAS shipped with R and
  within the rounding of one dot product against an optimized one, whose
  accumulation order is its own.

## numericals7 0.9.0

- [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
  gains `workers`, the number of R processes that the independent fits
  of a cross-validation’s folds may use, read by
  [`worker_count()`](https://statmodels7.github.io/numericals7/reference/worker_count.md).
  The same object carries both levels of the toolkit’s parallelism; the
  two do not nest (a fit inside a worker is sequential by construction),
  and the result does not depend on either count, bit for bit.

## numericals7 0.8.0

- The thread policy of the toolkit is defined here, at the root, since
  this is the one package below every compiled kernel.
  `n_threads(threads = 1)` constructs it,
  [`thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.md)
  reads the count, and
  [`local_threads()`](https://statmodels7.github.io/numericals7/reference/local_threads.md)
  applies it to RcppParallel’s process-level setting for one calling
  frame and restores the previous state on exit. The object is passed as
  an argument from the fit entry points
  (`statmodels7::statmod(threads =)`,
  `distributions7::fit_distrib(threads =)`) down to the kernels, so no
  package reads a setting that lives in another. The result of a fit
  does not depend on the count, bit for bit, because every parallel
  region in the toolkit decomposes its work over the elements of its
  output and never splits a reduction.

## numericals7 0.7.0

- The log-Bessel kernels are compiled.
  [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
  and
  [`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
  run scalar C++ loops over the same branches and formulas; the
  vectorized R implementations stay as internal twins
  (`.log_bessel_i_r`, `.log_bessel_k_r`) that a test compares with the
  compiled route on every branch. On a mixed workload of one million
  points: log K 2.9x faster, log I 1.1x. The u_k polynomial table is
  passed to the kernels once at load, so the two routes share the table.

## numericals7 0.6.0

- The jets are removed. Every production consumer carries its
  derivatives as written closed forms: the reparametrized families of
  distributions7 declare their map partials explicitly, and
  parameters7’s autoregressive family propagates its derivative arrays
  through the Levinson-Durbin recursion in compiled code, with the
  product rule written out per order. On the Poisson-inverse Gaussian
  kernels, the fixed composition overhead of the jet route was 2x to 36x
  the cost of the written forms.

## numericals7 0.5.0

- The logarithm of the modified Bessel functions, after Plesner,
  Sorensen and Hauberg (ICS 2024, arXiv:2409.08729):
  [`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
  and
  [`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
  work on the log scale (the power series through the log-sum-exp
  anchored at its largest term, the large-argument and large-order
  asymptotic expansions selected by the paper’s input-range table) and
  are finite and accurate wherever the logarithm itself is
  representable, while the unscaled I overflows from about x = 709, the
  unscaled K underflows from about x = 705, and the scaled I underflows
  past 1e5 or loses large orders. Two switching guards are tightened
  relative to the paper, measured on the Wronskian identity.
  [`log_bessel_i_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i_derivs.md)
  and
  [`log_bessel_k_derivs()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k_derivs.md)
  add the first four derivatives in the argument from the ratio identity
  and the Bessel equation. Against an Rcpp transcription of the same
  algorithm the compiled version was 0.9x to 2.7x as fast, so the
  package stayed pure R.

## numericals7 0.4.0

- Special functions:
  [`mills_ratio()`](https://statmodels7.github.io/numericals7/reference/mills_ratio.md)
  on the log scale, finite where the density and the distribution
  function both underflow;
  [`owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.md)
  through one batched
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  call, with the closed identities at `a = 0` and `a = Inf`;
  [`bessel_i_ratio()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio.md)
  through the exponentially scaled Bessel functions, with its four
  derivatives from the recurrence (`bessel_i_ratio_derivs()`) and its
  inverse by root finding with the inverse-function-rule derivatives
  ([`bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.md)).

## numericals7 0.3.0

- Quadrature and series vectorized over the parameters.
  [`quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.md)
  integrates one function at many parameter values by matrix evaluation
  (the nodes of every panel of every row go into the integrand in a
  single call per refinement pass), with the Gauss-Kronrod 7-15 pair
  supplying an error estimate from the same function values, adaptivity
  batched by row so that one hard row does not serialize the others, and
  rational maps for infinite endpoints. A row that does not reach the
  requested accuracy returns NA with a warning.
  [`series_vec()`](https://statmodels7.github.io/numericals7/reference/series_vec.md)
  does the same for series, in blocks with a per-row convergence mask
  and a tail guard that detects a block lying before the mode of a
  hump-shaped term.
  [`gauss_kronrod15()`](https://statmodels7.github.io/numericals7/reference/gauss_kronrod15.md)
  exposes the rule, pinned in the tests by its defining property:
  exactness to degree 13 for the embedded Gauss rule and to degree 23
  for the Kronrod extension (checked on the even degrees to 22), and a
  weight corrupted by 5% fails both the moment conditions and the gamma
  normalization.

- The von Mises variance through the per-theta fallback cost 25 ms per
  parameter value and scaled linearly, while regression models, where
  the parameters vary by observation, need hundreds of rows in one pass.

## numericals7 0.2.0

- The stencil library, replacing the three finite-difference
  implementations of the toolkit:
  [`fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.md)
  solves the Vandermonde system for any offsets and order (the
  construction basis7 had),
  [`fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.md)
  sizes a stencil from the order and the requested accuracy,
  [`fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.md)
  balances truncation against rounding and keeps every node inside a
  bounded domain, and
  [`fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.md)
  applies one stencil and never a composition of lower-order
  differences, since each numerical differentiation multiplies the error
  of the one before it. At the default accuracy these reproduce
  linkfunctions7’s four central stencils and basis7’s shapes exactly; at
  accuracy four they reproduce distributions7’s five-point `fd5_first`
  and `fd5_second`. The policy around a stencil (the order from which to
  fall back, and when a reference is reliable) stays with the callers.

## numericals7 0.1.0

- First release: the numerical layer of the statmodels7 toolkit,
  collecting what its packages had written separately.
  [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
  was written twice (distributions7 and parameters7, independently) and
  the finite-difference machinery three times.

- Jets, moved here from parameters7. A jet carries a value together with
  every partial derivative to fourth order and propagates them exactly
  through sums, products and a set of smooth functions (`exp`, `log`,
  `log1p`, `expm1`, an arbitrary power, `sqrt`, `gamma`, `lgamma`,
  `digamma`, `trigamma`, `sin`, `cos`). `Ops` and `Math` dispatch on the
  class, so a map written as `mu / gamma(1 + 1 / sigma)` is
  differentiated automatically. Comparison operators and the non-smooth
  functions are rejected, because a branch taken on a jet would keep one
  side’s derivatives and report them as those of the whole expression.

- The enumerations on which a higher-order chain rule rests, in one copy
  each:
  [`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
  (the multi-indices of a derivative, diagonal first at order two
  because that ordering is part of the interface),
  [`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
  (the Bell recursion, with blocks indexing positions so that a repeated
  variable carries its multiplicity), and
  [`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md)
  (the weak compositions, which are the support of a multinomial).
