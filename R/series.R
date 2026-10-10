#' @include quadrature.R
NULL

#' Sum One Series at Many Parameter Values
#'
#' @description
#' Computes \eqn{\sum_{k \ge k_0} t(k; \theta_i)} for every row \eqn{i} at once.
#' Terms are evaluated in blocks as one matrix, and a row retires as soon as it
#' converges, so rows that have converged are not evaluated again while slower
#' rows continue. This is the discrete counterpart of [quad_vec()], written for
#' the same reason: the toolkit's sums are one series at many parameter values.
#'
#' @details
#' # The term function
#'
#' `term(k, i)` receives two integer vectors of equal length and returns the
#' terms elementwise. `k` is the summation index and `i` gives the parameter
#' set to which each term belongs. A Poisson mass at a rate vector `lam` is
#'
#' ```
#' term <- function(k, i) dpois(k, lam[i])
#' ```
#'
#' # Convergence, row by row
#'
#' Write \eqn{b = \max(\mathrm{atol}, \mathrm{rtol}\,\lvert S_i \rvert)} for
#' the budget of row \eqn{i}, and \eqn{A_3} and \eqn{A_4} for the sums of the
#' absolute terms over the third and the fourth quarter of the last block. The
#' ratio \eqn{q = A_4/A_3} measures the decay of the terms at the end of the
#' block, and the neglected tail is estimated by continuing that decay
#' geometrically,
#'
#' \deqn{R = 2\,A_4\,\frac{q}{1 - q}.}
#'
#' A row retires when \eqn{q < 1} (or \eqn{A_4 = 0}), \eqn{R \le b}, and the
#' last block contributed at most \eqn{b}.
#'
#' The estimate bounds the tail of a series whose terms decay geometrically or
#' faster. For terms that decay like \eqn{k^{-p}}, geometric continuation
#' underestimates the tail by a factor that tends to \eqn{(p - 1)/p}, and the
#' factor 2 covers every \eqn{p \ge 2}; a slower decay, \eqn{1/k^2} included,
#' does not reach its budget within the default `max_terms` and returns `NA`.
#' The sum of \eqn{1/k^3} from one retires after about 75000 terms with a
#' relative error below the default `rtol`.
#'
#' Rising terms give \eqn{q > 1}, so a block that lies before the mode of a
#' hump-shaped term does not retire its row, however small its terms are:
#' `dpois(0:63, 300)` sums to 3.8e-62 and is still rising.
#'
#' # Rows of zeros
#'
#' A block whose terms are all zero retires a row only after the row has seen
#' a nonzero term, so the zeros past a finite support end the sum. A row that
#' has seen only zeros is probed ahead on a geometric grid of relative spacing
#' 1/128, up to the index \eqn{2^{31}}, in one call to `term`. If every probe
#' is zero the row is taken as identically zero and returns 0, as a
#' structurally null component of an expected Hessian does; otherwise the
#' summation resumes from the last zero probe, so leading terms that underflow
#' to zero (`dpois(k, 1e4)` for \eqn{k} below about 8000) do not end the sum.
#' The resumption assumes that the terms rise before the first nonzero one,
#' and a head narrower than the spacing of the grid (a Poisson mass at a rate
#' above about \eqn{10^8}) falls between the probes and gives 0.
#'
#' # Rows that do not converge
#'
#' A row still unconverged after `max_terms` terms returns `NA`, and one
#' warning lists the first eight such rows (followed by an ellipsis if there
#' are more).
#'
#' @param term The term function, a function of `k` and `i` as described above.
#' @param n The number of parameter rows, a positive whole number. It fixes the
#'   length of the answer and the range `i` takes.
#' @param from The first summation index, `0` by default. Pass `1` for a series
#'   indexed from one.
#' @param atol,rtol The absolute and relative budgets per row, defaulting to
#'   `1e-12` and `1e-10`. A row is judged against the larger of the two, so
#'   `atol` governs a sum near zero and `rtol` a large one. Both are tighter
#'   than [quad_vec()]'s, a term being cheaper than a panel.
#' @param max_terms The number of terms after which a row that has not
#'   converged returns `NA`, `100000` by default. Terms are evaluated in whole
#'   blocks, so a row may consume up to the first multiple of `block` at or
#'   above this number; the terms skipped after a probe and the probes
#'   themselves are not counted.
#' @param block How many terms are evaluated per pass, a whole number of at
#'   least 4, `64` by default. Its last two quarters give the decay ratio of
#'   the tail estimate, so a very small block makes the estimate noisy.
#'
#' @return A numeric vector of sums, one per row, of length `n`. `NA` in any row
#'   that did not converge within `max_terms`, with a warning that lists the
#'   first eight of those rows.
#'
#' @seealso [quad_vec()], the continuous counterpart.
#'
#' @examples
#' # Four geometric series against the closed form.
#' r <- c(0.1, 0.5, 0.9, 0.99)
#' series_vec(function(k, i) r[i]^k, n = 4) - 1 / (1 - r)
#'
#' # Poisson masses sum to one at every rate, in one call.
#' lam <- c(0.5, 4, 60)
#' series_vec(function(k, i) dpois(k, lam[i]), n = 3)
#'
#' # The first 64 terms of a Poisson at rate 300 sum to almost nothing and are
#' # still rising, so the row does not retire there.
#' sum(dpois(0:63, 300))
#' series_vec(function(k, i) dpois(k, 300), n = 1)
#'
#' # At rate 1e4 the first 8000 or so terms underflow to zero, and the sum
#' # continues past them.
#' series_vec(function(k, i) dpois(k, 1e4), n = 1)
#'
#' # Series indexed from one, against their closed forms: a factorial decay
#' # and a polynomial one.
#' series_vec(function(k, i) 1 / factorial(k), n = 1, from = 1L) - (exp(1) - 1)
#' series_vec(function(k, i) 1 / k^4, n = 1, from = 1L) - pi^4 / 90
#'
#' # A divergent series is refused, not estimated.
#' suppressWarnings(series_vec(function(k, i) 1 / (k + 1), n = 1,
#'                             max_terms = 500L))
#'
#' @export
series_vec <- function(term, n, from = 0L, atol = 1e-12, rtol = 1e-10,
                       max_terms = 100000L, block = 64L) {
  if (!is.numeric(block) || length(block) != 1L || is.na(block) ||
      block < 4 || block != trunc(block)) {
    stop("'block' must be a single whole number of at least 4.", call. = FALSE)
  }
  block <- as.integer(block)
  qb <- block %/% 4L
  q3 <- block - 2L * qb + seq_len(qb)
  q4 <- block - qb + seq_len(qb)
  acc <- numeric(n)
  seen <- logical(n)
  kst <- rep(as.integer(from), n)
  used <- integer(n)
  active <- seq_len(n)
  failed <- integer(0)
  offs <- seq.int(0L, block - 1L)
  # probe distances for a row that has seen only zeros: zero, then a
  # geometric grid of ratio 1 + 1/128 from one block to 2^31
  nj <- ceiling(log(2^31 / block) / log1p(1 / 128))
  steps <- c(0, unique(floor(block * (1 + 1 / 128)^(0:nj))))

  while (length(active)) {
    na <- length(active)
    K <- outer(kst[active], offs, "+")
    vals <- term(as.vector(K), rep(active, times = block))
    M <- matrix(as.numeric(vals), na, block)

    contrib <- rowSums(M)
    acc[active] <- acc[active] + contrib
    budget <- pmax(atol, rtol * abs(acc[active]))
    # the tail is estimated by continuing geometrically the decay between the
    # last two quarters of the block, doubled: exact for geometric terms,
    # conservative for faster ones, and covering k^-p for p >= 2, where the
    # geometric continuation alone underestimates by (p - 1)/p. Rising terms
    # give a ratio above one, which keeps a hump-shaped term alive before its
    # mode. A block of zeros says nothing about the decay, so it retires a row
    # only once that row has seen a nonzero term.
    A <- abs(M)
    a3 <- rowSums(A[, q3, drop = FALSE])
    a4 <- rowSums(A[, q4, drop = FALSE])
    seen[active] <- seen[active] | rowSums(A) > 0
    q <- a4 / a3
    tail <- ifelse(a4 == 0, 0, 2 * a4 * q / (1 - q))
    done <- (a4 < a3 | a4 == 0) & tail <= budget & abs(contrib) <= budget &
      seen[active]
    done[is.na(done)] <- FALSE
    used[active] <- used[active] + block
    kst[active] <- kst[active] + block

    # A row that has seen only zeros either has an underflowed head
    # (dpois(k, 1e4) is zero below k ~ 8000) or is zero throughout (a
    # structurally null component of an expected Hessian). It is probed on a
    # geometric grid of relative spacing 1/128: if every probe up to 2^31 is
    # zero the row is taken as identically zero and returns 0, and otherwise
    # it resumes from the last zero probe, the head being assumed to rise
    # before its first nonzero term.
    z <- which(!done & !seen[active])
    if (length(z)) {
      rz <- active[z]
      P <- outer(as.numeric(kst[rz]), steps, "+")
      P[P > .Machine$integer.max] <- NA
      okp <- !is.na(P)
      pv <- rep(NA_real_, length(P))
      pv[okp] <- as.numeric(term(as.integer(P[okp]),
                                 rep(rz, times = length(steps))[okp]))
      hit <- matrix(okp & (is.na(pv) | pv != 0), length(rz), length(steps))
      first <- max.col(hit, ties.method = "first")
      found <- rowSums(hit) > 0
      done[z[!found]] <- TRUE
      jump <- found & first > 1L
      kst[rz[jump]] <- as.integer(P[cbind(which(jump), first[jump] - 1L)])
    }

    active <- active[!done]
    out <- active[used[active] >= max_terms]
    failed <- c(failed, out)
    active <- setdiff(active, out)
  }

  active <- sort(failed)
  if (length(active)) {
    acc[active] <- NA_real_
    shown <- paste(utils::head(active, 8L), collapse = ", ")
    if (length(active) > 8L) shown <- paste0(shown, ", ...")
    warning(sprintf(
      "series_vec: %d row(s) did not converge within %d terms and return NA: rows %s.",
      length(active), max_terms, shown
    ), call. = FALSE)
  }
  acc
}
