#' @include stencil.R
NULL

# Quadrature vectorized over the parameters. The toolkit's integrals are
# almost always the same integrand at many parameter values -- a moment per
# observation, a normalizing constant per row of a model matrix -- and
# calling a scalar integrator in a loop pays its overhead once per value.
# Here the parameter dimension is a matrix dimension: one call to the
# integrand evaluates every node of every panel of every parameter row.

#' The Gauss-Kronrod 7-15 Pair
#'
#' @description
#' Returns the nodes and the two sets of weights of the 15-point Kronrod
#' extension of the 7-point Gauss rule on \eqn{[-1, 1]}. One set of nodes
#' carries both rules, so a single evaluation of the integrand gives an
#' estimate and, from the difference of the two rules, an error for it. An
#' adaptive routine reads that error to decide where to refine, and pays
#' nothing for it beyond the estimate.
#'
#' @details
#' The eight Kronrod-only nodes carry a Gauss weight of zero, so both rules are
#' formed from one matrix of function values by two weighted sums.
#'
#' The constants are the classical ones. The tests pin them by their defining
#' property, checking that the 7-point rule integrates polynomials of degree 13
#' exactly and the 15-point one degree 22, which catches a transcription error
#' that comparing digits against a table would only catch if the table were
#' right.
#'
#' @return A list of three numeric vectors, each of length 15:
#'   \describe{
#'     \item{`nodes`}{the abscissae on \eqn{[-1, 1]}, ascending and symmetric
#'       about zero.}
#'     \item{`wk`}{the Kronrod weights, all positive.}
#'     \item{`wg`}{the embedded Gauss weights, zero at the eight
#'       Kronrod-only nodes.}
#'   }
#'
#' @references
#' Kronrod, A. S. (1965). *Nodes and Weights of Quadrature Formulas*.
#' Consultants Bureau, New York.
#'
#' Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W. and Kahaner, D. K.
#' (1983). *QUADPACK: A Subroutine Package for Automatic Integration*.
#' Springer.
#'
#' @seealso [quad_vec()], which integrates with this pair by default.
#'
#' @examples
#' r <- gauss_kronrod15()
#'
#' # The weights of a rule sum to the length of its interval.
#' c(kronrod = sum(r$wk), gauss = sum(r$wg))
#'
#' # The Gauss rule lives on eight of the fifteen nodes, so both rules come
#' # from one set of function values.
#' sum(r$wg == 0)
#'
#' # Its defining property: exact on polynomials to degree 13.
#' p <- function(x, d) x^d
#' vapply(c(13, 14), function(d) {
#'   exact <- if (d %% 2 == 0) 2 / (d + 1) else 0
#'   sum(r$wg * p(r$nodes, d)) - exact
#' }, numeric(1))
#'
#' @export
gauss_kronrod15 <- function() {
  xh <- c(0.991455371120813, 0.949107912342759, 0.864864423359769,
          0.741531185599394, 0.586087235467691, 0.405845151377397,
          0.207784955007898, 0)
  wkh <- c(0.022935322010529, 0.063092092629979, 0.104790010322250,
           0.140653259715525, 0.169004726639267, 0.190350578064785,
           0.204432940075298, 0.209482141084728)
  wgh <- c(0, 0.129484966168870, 0, 0.279705391489277,
           0, 0.381830050505119, 0, 0.417959183673469)
  list(
    nodes = c(-xh[1:7], xh[8], rev(xh[1:7])),
    wk = c(wkh[1:7], wkh[8], rev(wkh[1:7])),
    wg = c(wgh[1:7], wgh[8], rev(wgh[1:7]))
  )
}


#' Integrate One Function at Many Parameter Values
#'
#' @description
#' Computes \eqn{\int_{a_i}^{b_i} f(x; \theta_i)\,\mathrm{d}x} for every row
#' \eqn{i} at once. The nodes of every panel of every row reach `f` in a single
#' call per refinement pass, so the parameter index is a matrix dimension and
#' not a loop. The toolkit's integrals are almost always this shape, one
#' integrand at many parameter values, and a scalar integrator called in a loop
#' pays its overhead once per value.
#'
#' @details
#' # The integrand contract
#'
#' `f(x, i)` receives a numeric matrix `x` of evaluation points and an integer
#' vector `i` with one entry per row of `x`, saying which parameter set that row
#' belongs to. It returns the values elementwise, either as a matrix shaped like
#' `x` or as a vector in column-major order.
#'
#' Elementwise recycling does the indexing. For a gamma mean at parameter
#' vectors of length \eqn{n}:
#'
#' ```
#' f <- function(x, i) x * dgamma(x, shape = shp[i], rate = rt[i])
#' ```
#'
#' `shp[i]` has one entry per row of `x` and recycles down each column.
#'
#' # Batched adaptivity
#'
#' Each panel carries the error estimate of its Gauss-Kronrod pair. A row
#' converges when the *sum* of its panel errors fits the budget
#' \eqn{\max(\mathrm{atol}, \mathrm{rtol}\,\lvert I_i \rvert)}. Until then the
#' row's worst panels are bisected, and the splits from every row join the next
#' single evaluation, so one hard row refines its own panels without
#' serializing the others.
#'
#' The budget is judged on the sum for a reason worth knowing, because the
#' obvious alternative fails. Giving each panel a share proportional to its
#' length cannot integrate an endpoint singularity at all: near such a point the
#' error stays concentrated in the innermost panel however deep the bisection
#' goes, so a per-length share demands of that panel an accuracy no depth
#' reaches. Judging the sum lets the smooth panels carry the row.
#'
#' # Infinite endpoints
#'
#' Mapped to finite ones by the rational transforms \eqn{x = a + t/(1-t)},
#' \eqn{x = b - t/(1-t)} and \eqn{x = t/(1-t^2)}, whose Jacobians multiply the
#' integrand. Rows of different kinds may share one call, so a vector of
#' endpoints mixing finite and infinite costs nothing extra.
#'
#' # A failure is reported as one
#'
#' A row whose panels still exceed their budget at `max_depth` returns `NA`, and
#' one warning names every such row. An `NA` says the accuracy was not reached;
#' a plausible number would say nothing and be believed.
#'
#' @param f The integrand, obeying the contract above.
#' @param lower,upper Numeric vectors of endpoints, recycled to a common length,
#'   either of which may be infinite. Every `lower` must be strictly below its
#'   `upper`, or the call throws.
#' @param atol,rtol The absolute and relative error budgets per row, defaulting
#'   to `1e-10` and `1e-8`. A row is judged against the larger of the two, so
#'   `atol` governs an integral near zero and `rtol` a large one.
#' @param max_depth The greatest number of bisections one panel may undergo,
#'   `48` by default. It is the lever for an endpoint singularity, and the
#'   measured reach is narrower than "integrable" suggests: at the default a
#'   gamma density of shape 0.5 converges and one of shape 0.45 does not, while
#'   `max_depth = 200` reaches shape 0.2 and still not 0.1. A row past the
#'   budget returns `NA`.
#' @param rule The embedded quadrature pair, [gauss_kronrod15()] by default.
#'   Any list of `nodes`, `wk` and `wg` of equal length serves.
#'
#' @return A numeric vector of integrals, one per row, of the recycled length of
#'   `lower` and `upper`. `NA` in any row that did not reach the requested
#'   accuracy, with a warning naming those rows.
#'
#' @seealso [series_vec()] for the discrete counterpart, [gauss_kronrod15()] for
#'   the default rule.
#'
#' @examples
#' # Thirty gamma densities integrate to one, in one call rather than thirty.
#' shp <- seq(0.5, 15, length.out = 30)
#' f <- function(x, i) dgamma(x, shape = shp[i], rate = 1)
#' range(quad_vec(f, lower = 0, upper = rep(Inf, 30)) - 1)
#'
#' # Their means, against the closed form.
#' g <- function(x, i) x * dgamma(x, shape = shp[i], rate = 1)
#' range(quad_vec(g, 0, rep(Inf, 30)) - shp)
#'
#' # A shape below one puts an integrable singularity at the origin. The
#' # sum-judged budget reaches shape 0.5 at the default depth, and a harsher
#' # one needs a deeper budget rather than a looser tolerance.
#' quad_vec(function(x, i) dgamma(x, shape = 0.5, rate = 1), 0, Inf)
#' quad_vec(function(x, i) dgamma(x, shape = 0.2, rate = 1), 0, Inf,
#'          max_depth = 200L)
#'
#' # Endpoints of different kinds share a call.
#' quad_vec(function(x, i) dnorm(x), c(-Inf, -1, 0), c(0, 1, Inf))
#'
#' # A divergent integral is refused, not estimated.
#' suppressWarnings(quad_vec(function(x, i) 1 / x, 0, 1))
#'
#' @export
quad_vec <- function(f, lower, upper, atol = 1e-10, rtol = 1e-8,
                     max_depth = 48L, rule = gauss_kronrod15()) {
  n <- max(length(lower), length(upper))
  lower <- rep_len(lower, n)
  upper <- rep_len(upper, n)
  if (any(lower >= upper)) {
    stop("Every 'lower' must be smaller than its 'upper'.", call. = FALSE)
  }

  # 0 finite, 1 upper-infinite, 2 lower-infinite, 3 doubly infinite
  type <- ifelse(is.finite(lower) & is.finite(upper), 0L,
          ifelse(is.finite(lower), 1L,
          ifelse(is.finite(upper), 2L, 3L)))

  # panels live in the TRANSFORMED coordinate
  t0 <- ifelse(type == 0L, lower, ifelse(type == 3L, -1, 0))
  t1 <- ifelse(type == 0L, upper, 1)
  L0 <- t1 - t0

  # evaluate the raw integrand times the Jacobian of the row's transform
  eval_panels <- function(orig, a, b) {
    half <- (b - a) / 2
    Tm <- outer(a + half, rep(1, 15L)) + outer(half, rule$nodes)
    ty <- type[orig]
    X <- Tm
    J <- matrix(1, nrow(Tm), 15L)
    if (any(ty == 1L)) {
      r <- ty == 1L
      X[r, ] <- lower[orig[r]] + Tm[r, , drop = FALSE] / (1 - Tm[r, , drop = FALSE])
      J[r, ] <- 1 / (1 - Tm[r, , drop = FALSE])^2
    }
    if (any(ty == 2L)) {
      r <- ty == 2L
      X[r, ] <- upper[orig[r]] - Tm[r, , drop = FALSE] / (1 - Tm[r, , drop = FALSE])
      J[r, ] <- 1 / (1 - Tm[r, , drop = FALSE])^2
    }
    if (any(ty == 3L)) {
      r <- ty == 3L
      tt <- Tm[r, , drop = FALSE]
      X[r, ] <- tt / (1 - tt^2)
      J[r, ] <- (1 + tt^2) / (1 - tt^2)^2
    }
    Fv <- f(X, orig)
    Fm <- matrix(as.numeric(Fv), nrow(Tm), 15L) * J
    K <- half * as.numeric(Fm %*% rule$wk)
    G <- half * as.numeric(Fm %*% rule$wg)
    list(K = K, err = abs(K - G))
  }

  # Panels persist with their cached estimates until their row converges or
  # fails; only newly created panels are evaluated, so each pass costs one
  # call to f on the panels the previous pass split.
  result <- rep(NA_real_, n)
  failed <- rep(FALSE, n)
  live <- rep(TRUE, n)
  orig <- seq_len(n)
  a <- t0
  b <- t1
  depth <- rep(0L, n)
  K <- rep(NA_real_, n)
  err <- rep(Inf, n)
  fresh <- rep(TRUE, n)

  iter <- 0L
  while (any(live) && length(orig) && iter < 40L + 4L * max_depth) {
    iter <- iter + 1L
    if (any(fresh)) {
      ev <- eval_panels(orig[fresh], a[fresh], b[fresh])
      K[fresh] <- ev$K
      err[fresh] <- ev$err
      err[fresh][!is.finite(err[fresh])] <- Inf
      fresh <- rep(FALSE, length(orig))
    }

    # Convergence is judged on the SUM of a row's panel errors, as a global
    # adaptive integrator judges it. Allocating the budget by panel length
    # instead forces a uniform error density, which a row with an integrable
    # endpoint singularity cannot meet at any depth: that was the first
    # version, and it returned NA for a beta(0.8, 1.2) normalization that
    # stats::integrate handles routinely.
    tot <- rowsum(K, orig, reorder = FALSE)
    te <- rowsum(err, orig, reorder = FALSE)
    rows <- as.integer(rownames(tot))
    budget <- pmax(atol, rtol * abs(tot[, 1L]))
    # which() drops the NA a non-finite integrand puts into its row's budget;
    # such a row keeps its Inf error and fails at max_depth instead
    done_rows <- rows[which(te[, 1L] <= budget)]

    if (length(done_rows)) {
      result[done_rows] <- tot[match(done_rows, rows), 1L]
      live[done_rows] <- FALSE
    }

    still <- setdiff(rows, done_rows)
    if (!length(still)) {
      sel <- !(orig %in% done_rows)
      orig <- orig[sel]; a <- a[sel]; b <- b[sel]
      depth <- depth[sel]; K <- K[sel]; err <- err[sel]; fresh <- fresh[sel]
      next
    }

    # Split the worst panels of each unconverged row. A row fails when
    # nothing is left to split, and ALSO when the panels already at
    # max_depth carry, by themselves, more error than the whole budget:
    # that error is irreducible, so no amount of work elsewhere can save
    # the row. Without the second clause the loop keeps splitting the
    # smooth panels' rounding noise -- errors all of one tiny magnitude,
    # so half the list clears the threshold at every pass and the panel
    # count doubles until memory runs out. All of it is computed by row
    # aggregates rather than a loop over rows, whose which() scan cost
    # O(rows x panels) per pass.
    sp <- depth < max_depth
    g <- match(orig, rows)
    irr <- rowsum(err * !sp, orig, reorder = FALSE)[, 1L]
    nsp <- rowsum(sp + 0L, orig, reorder = FALSE)[, 1L]
    is_still <- rows %in% still
    fail_rows <- rows[is_still & (nsp == 0L | irr > budget)]
    if (length(fail_rows)) {
      failed[fail_rows] <- TRUE
      live[fail_rows] <- FALSE
      result[fail_rows] <- NA_real_
    }
    row_ok <- is_still & !(rows %in% fail_rows)
    mx <- vapply(split(ifelse(sp, err, -Inf), g), max, numeric(1))
    do_split <- sp & row_ok[g] & err >= 0.5 * mx[g]

    keep <- live[orig] & !do_split
    picked <- do_split & live[orig]
    so <- orig[picked]
    sa <- a[picked]
    sb <- b[picked]
    sd <- depth[picked] + 1L
    sm <- (sa + sb) / 2

    orig <- c(orig[keep], so, so)
    a <- c(a[keep], sa, sm)
    b <- c(b[keep], sm, sb)
    depth <- c(depth[keep], sd, sd)
    K <- c(K[keep], rep(NA_real_, 2L * length(so)))
    err <- c(err[keep], rep(Inf, 2L * length(so)))
    fresh <- c(rep(FALSE, sum(keep)), rep(TRUE, 2L * length(so)))
  }

  # rows still live when the iteration guard fires are failures too
  if (any(live)) {
    failed[live] <- TRUE
    result[live] <- NA_real_
  }
  banked <- result

  if (any(failed)) {
    banked[failed] <- NA_real_
    which_bad <- which(failed)
    shown <- paste(utils::head(which_bad, 8L), collapse = ", ")
    if (length(which_bad) > 8L) shown <- paste0(shown, ", ...")
    warning(sprintf(
      "quad_vec: %d row(s) did not reach the requested accuracy at max_depth = %d and return NA: rows %s.",
      length(which_bad), max_depth, shown
    ), call. = FALSE)
  }
  banked
}
