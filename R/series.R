#' @include quadrature.R
NULL

#' Sum One Series at Many Parameter Values
#'
#' @description
#' \eqn{\sum_{k \ge k_0} t(k; \theta_i)} for every row \eqn{i} at once, in
#' blocks of terms evaluated as one matrix, with rows retiring as they
#' converge so that a slow row does not keep the finished ones paying.
#'
#' @details
#' \strong{The term contract.} \code{term(k, i)} receives two integer
#' vectors of equal length and returns the terms elementwise: \code{k} is
#' the summation index, \code{i} says which parameter set. A Poisson mass,
#' for a rate vector \code{lam}, is
#' \code{function(k, i) dpois(k, lam[i])}.
#'
#' \strong{Convergence, per row.} A row retires when its last block
#' contributed less than \eqn{\max(\mathrm{atol},
#' \mathrm{rtol}\,\lvert S_i \rvert)}, the final term of the block is itself
#' below that budget, \emph{and} the terms are not growing across the block.
#' The last two conditions are the tail guard: a hump-shaped term -- a
#' Poisson mass at a large rate, say -- can open with a block that sums to
#' nearly nothing and ends on a term smaller still, while the whole series
#' lies beyond it; that its terms are rising is what says the block came
#' before the mode rather than after it. The terms are assumed eventually
#' decreasing in magnitude, which every series in the toolkit satisfies.
#'
#' \strong{Rejection over plausibility.} A row not converged after
#' \code{max_terms} terms returns \code{NA} with a warning naming it.
#'
#' @param term The term function, as described above.
#' @param n The number of parameter rows.
#' @param from The first summation index.
#' @param atol,rtol The absolute and relative budgets per row.
#' @param max_terms The largest number of terms a row may consume.
#' @param block How many terms are evaluated per pass.
#'
#' @return A numeric vector of sums, one per row, with \code{NA} where the
#'   series did not converge.
#'
#' @seealso \code{\link{quad_vec}}
#'
#' @examples
#' # geometric series against the closed form
#' r <- c(0.1, 0.5, 0.9, 0.99)
#' series_vec(function(k, i) r[i]^k, n = 4)
#' 1 / (1 - r)
#'
#' # Poisson masses sum to one for every rate at once
#' lam <- c(0.5, 4, 60)
#' series_vec(function(k, i) dpois(k, lam[i]), n = 3)
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
