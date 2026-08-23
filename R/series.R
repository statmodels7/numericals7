#' @include quadrature.R
NULL

#' Sum One Series at Many Parameter Values
#'
#' @description
#' Computes \eqn{\sum_{k \ge k_0} t(k; \theta_i)} for every row \eqn{i} at once.
#' Terms are evaluated in blocks as one matrix, and a row retires as soon as it
#' converges, so a slow row does not keep the finished ones paying. This is the
#' discrete counterpart of [quad_vec()] and it exists for the same reason: the
#' toolkit's sums are one series at many parameter values.
#'
#' @details
#' # The term contract
#'
#' `term(k, i)` receives two integer vectors of equal length and returns the
#' terms elementwise. `k` is the summation index and `i` says which parameter
#' set each term belongs to. A Poisson mass at a rate vector `lam` is
#'
#' ```
#' term <- function(k, i) dpois(k, lam[i])
#' ```
#'
#' # Convergence, row by row
#'
#' A row retires when three conditions hold together. Its last block
#' contributed less than \eqn{\max(\mathrm{atol}, \mathrm{rtol}\,\lvert S_i
#' \rvert)}; the final term of that block is itself below the same budget; and
#' the terms are not growing across the block.
#'
#' The last two are the tail guard, and the third is the one that earns its
#' place. A hump-shaped term can open with a block that sums to nearly nothing
#' and ends on a term smaller still, while the whole series lies beyond it.
#' `dpois(0:63, 300)` sums to 3.8e-62 and ends at 3.0e-62, and every series
#' there is still to come. That the terms are rising is the signal that the
#' block sits before the mode.
#'
#' Terms are assumed eventually decreasing in magnitude, which every series in
#' the toolkit satisfies.
#'
#' # A failure is reported as one
#'
#' A row still unconverged after `max_terms` terms returns `NA`, and one warning
#' names every such row.
#'
#' @param term The term function, obeying the contract above.
#' @param n The number of parameter rows, a positive whole number. It fixes the
#'   length of the answer and the range `i` takes.
#' @param from The first summation index, `0` by default. Pass `1` for a series
#'   indexed from one.
#' @param atol,rtol The absolute and relative budgets per row, defaulting to
#'   `1e-12` and `1e-10`. A row is judged against the larger of the two, so
#'   `atol` governs a sum near zero and `rtol` a large one. Both are tighter
#'   than [quad_vec()]'s, a term being cheaper than a panel.
#' @param max_terms The greatest number of terms one row may consume, `100000`
#'   by default.
#' @param block How many terms are evaluated per pass, `64` by default. It is
#'   also the window the tail guard judges growth over, so a very small block
#'   weakens the guard.
#'
#' @return A numeric vector of sums, one per row, of length `n`. `NA` in any row
#'   that did not converge within `max_terms`, with a warning naming those rows.
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
#' # The tail guard at work. The first 64 terms of a Poisson at rate 300 sum
#' # to almost nothing and are still rising, so the row does not retire there.
#' sum(dpois(0:63, 300))
#' series_vec(function(k, i) dpois(k, 300), n = 1)
#'
#' # A series indexed from one, against its closed form. Convergence has to be
#' # geometric or better: the terms of sum 1/k^2 decay too slowly to retire a
#' # row within max_terms, and that row comes back NA.
#' series_vec(function(k, i) 1 / factorial(k), n = 1, from = 1L) - (exp(1) - 1)
#'
#' # A divergent series is refused, not estimated.
#' suppressWarnings(series_vec(function(k, i) 1 / (k + 1), n = 1,
#'                             max_terms = 500L))
#'
#' @export
series_vec <- function(term, n, from = 0L, atol = 1e-12, rtol = 1e-10,
                       max_terms = 100000L, block = 64L) {
  acc <- numeric(n)
  active <- seq_len(n)
  k0 <- as.integer(from)
  used <- 0L

  while (length(active) && used < max_terms) {
    kb <- seq.int(k0, k0 + block - 1L)
    K <- rep(kb, each = length(active))
    I <- rep(active, times = block)
    vals <- term(K, I)
    M <- matrix(as.numeric(vals), length(active), block)

    contrib <- rowSums(M)
    acc[active] <- acc[active] + contrib
    budget <- pmax(atol, rtol * abs(acc[active]))
    # three conditions to retire: the block contributed little, its last term
    # is itself small, and the terms are not growing across the block. The
    # third is what keeps a hump-shaped term alive before its mode -- for a
    # Poisson mass at rate 300 the block over k = 0..63 sums to nearly
    # nothing AND ends on a term of about 1e-52, so the first two conditions
    # alone would retire a series that has not begun.
    done <- abs(contrib) <= budget &
      abs(M[, block]) <= budget &
      abs(M[, block]) <= abs(M[, 1L])
    done[is.na(done)] <- FALSE

    active <- active[!done]
    k0 <- k0 + block
    used <- used + block
  }

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
