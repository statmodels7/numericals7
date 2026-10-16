#' @include stencil.R
NULL

# The tensor-product stencil. distributions7 wrote it twice (the derivatives
# in the parameters, and the mixed derivatives in the response and the
# parameters) and parameters7 a third time; the construction is the same in
# all three: each component is ONE product of central 1-D stencils applied to
# the highest analytic quantity below it, never a difference of a difference.
# The plan is integer structure and is computed once per request shape; the
# node values are calls to R functions; the weighted sums run in C++.

#' An Analytic Quantity Declared to a Tensor Stencil
#'
#' @description
#' Declares a function that returns analytic partial derivatives of the
#' function differentiated by [fd_tensor()], together with the multi-indices
#' of the derivatives that it returns. [fd_tensor()] differences these values
#' instead of the function itself wherever that leaves a lower order to the
#' stencil.
#'
#' @details
#' A multi-index is a vector of counts, one per coordinate of `x`: the row
#' `c(2, 0, 1)` over the coordinates `(a, b, c)` is
#' \eqn{\partial^3 f / \partial a^2 \partial c}. The rows of `indices` are the
#' multi-indices of the values that `fun` returns, in the order in which it
#' returns them. With `order = k` the rows are every multi-index of total
#' order \eqn{k}, in the order of [`tuple_indices(p, k)`][tuple_indices] over
#' the \eqn{p} coordinates; [fd_tensor()] expands the shortcut from the
#' length of its `x`.
#'
#' `fun` is called with a named list of the same shape as `x` and returns
#' either a list of numeric vectors, one per row of `indices`, or a numeric
#' matrix with one column per row, or a single vector when there is one row.
#' Each vector has length 1 or \eqn{n}. The names of the returned list are not
#' read.
#'
#' @param fun A function of one argument, the named list of coordinates.
#' @param order A single non-negative whole number, the total order of every
#'   derivative that `fun` returns, or `NULL`.
#' @param indices An integer matrix of counts with one row per returned
#'   derivative and one column per coordinate, or `NULL`. Exactly one of
#'   `order` and `indices` is given.
#'
#' @return An object of class `fd_known`, with the properties `fun`, `order`
#'   and `indices`.
#'
#' @seealso [fd_tensor()], which reads these declarations.
#'
#' @examples
#' f <- function(x) x$a^3 * exp(x$b)
#' # the gradient, declared by its order
#' g <- fd_known(function(x) list(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b)),
#'               order = 1)
#' # the same quantity, declared by its multi-indices
#' g2 <- fd_known(function(x) cbind(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b)),
#'                indices = rbind(c(1, 0), c(0, 1)))
#' g2@indices
#' @export
fd_known <- S7::new_class(
  name = "fd_known",
  properties = list(
    fun = S7::class_function,
    order = S7::class_any,
    indices = S7::class_any
  ),
  constructor = function(fun, order = NULL, indices = NULL) {
    if (is.null(order) == is.null(indices)) {
      stop("Give exactly one of 'order' and 'indices'.", call. = FALSE)
    }
    if (!is.null(order)) {
      if (!is_whole(order) || order < 0) {
        stop("'order' must be a single non-negative whole number.",
             call. = FALSE)
      }
      order <- as.integer(order)
    } else {
      indices <- as_counts(indices, NULL, "indices")
    }
    S7::new_object(S7::S7_object(), fun = fun, order = order,
                   indices = indices)
  }
)

#' Partial Derivatives by One Tensor-Product Stencil per Component
#'
#' @description
#' Computes partial derivatives of a function of several coordinates. Each
#' component is a single linear combination of values of the highest analytic
#' quantity below it (the function itself, or a derivative declared through
#' [fd_known()]): a tensor product of central one-dimensional stencils, one in
#' each direction that still has to be differentiated, and never a difference
#' of a difference.
#'
#' @details
#' # The source of each component
#'
#' A component is a multi-index \eqn{\alpha}, a vector of counts over the
#' coordinates. Every declared multi-index \eqn{\beta \le \alpha} can serve as
#' its source, with residual \eqn{\alpha - \beta} and residual order
#' \eqn{r = |\alpha| - |\beta|}; the function itself is the source with
#' \eqn{\beta = 0}. The source with the smallest \eqn{r} is used. Among
#' sources of equal \eqn{r}, the quantity that declares more derivatives is
#' preferred, and then the one declared first in `known`. Within one quantity
#' the analytic part takes the most repeated coordinates of \eqn{\alpha}
#' first, so that the residual is spread over the directions as evenly as
#' possible; among equally even residuals the earlier coordinates go to the
#' analytic part. The attribute `"fd_source"` of the result records the
#' choice.
#'
#' # The stencil
#'
#' With residual \eqn{\rho = \alpha - \beta} and source \eqn{g}, the component
#' is
#' \deqn{\sum_{k} \Bigl(\prod_{j:\,\rho_j > 0} w^{(\rho_j)}_{k_j} /
#'   h_j^{\rho_j}\Bigr)\, g(x_1 + k_1 h_1, \dots, x_p + k_p h_p),}
#' where \eqn{w^{(m)}} are the central weights of order \eqn{m} on the nodes
#' of [`fd_offsets(m, accuracy)`][fd_offsets], computed by Fornberg's
#' recursion and exactly symmetric, and the sum runs over the nodes with a
#' nonzero weight. A node shared by several components is evaluated once.
#'
#' # The step
#'
#' A component of residual order \eqn{r} uses, in direction \eqn{j}, the step
#' \eqn{h_j = \varepsilon^{1/(r + a)}\max(1, |x_j|)}, with \eqn{a} the
#' accuracy: rounding grows as \eqn{\varepsilon / h^r} and the truncation as
#' \eqn{h^a}, and the two balance there. `h_rel` replaces
#' \eqn{\varepsilon^{1/(r + a)}}. Where `bounds` gives a finite bound, the
#' step is reduced to \eqn{0.49\,d / R_j}, with \eqn{d} the distance to the
#' bound and \eqn{R_j} the largest offset that the stencils of that residual
#' order use in direction \eqn{j}, so every node lies inside the domain.
#'
#' # The computation
#'
#' The plan (the source of each component, the distinct nodes and the sparse
#' list of weights) depends only on the multi-indices and the accuracy. It is
#' computed in C++ and cached for the session. The node values are computed in
#' R, one call of the source function per node, and the weighted sums are
#' formed in C++ over the vectors of length \eqn{n}.
#'
#' @param f A function of one argument, the named list of coordinates,
#'   returning a numeric vector of length 1 or \eqn{n}. When every coordinate
#'   has length 1, \eqn{n} is the length of the values (a response held by
#'   `f`, for instance).
#' @param x A named list of \eqn{p} numeric coordinates, each of length 1 or a
#'   common length \eqn{n}. The derivatives are taken at these values,
#'   elementwise along the vectors of length \eqn{n}.
#' @param order A single positive whole number: every component of this total
#'   order is computed, in the order of [`tuple_indices(p, order)`][tuple_indices].
#'   Exactly one of `order` and `components` is given.
#' @param components An integer matrix of counts, one row per component and
#'   one column per coordinate.
#' @param known A list of [fd_known()] objects, or a single one.
#' @param bounds `NULL`, or a list with one entry per coordinate (matched by
#'   name when named), each `NULL` or a numeric vector `c(lower, upper)`.
#' @param accuracy The order of accuracy of every one-dimensional stencil, a
#'   positive even whole number. The default, 4, makes the truncation error
#'   of order \eqn{h^4}; 2 uses fewer nodes in each direction, and the number
#'   of evaluations falls by a factor of two to four.
#' @param h_rel `NULL` for the step above, or a single positive number that
#'   replaces \eqn{\varepsilon^{1/(r + a)}}.
#'
#' @return A named list of numeric vectors of length \eqn{n}, one per
#'   component, named by the coordinates of each multi-index joined with `_`
#'   (`"a_a_b"` for \eqn{\partial^3 f / \partial a^2 \partial b}). The
#'   attribute `"fd_source"` is a data frame with one row per component and
#'   the columns `component`, `source` (0 for `f`, \eqn{k} for `known[[k]]`),
#'   `index` (the row of that quantity's multi-indices) and `residual` (the
#'   order left to the stencil).
#'
#' @seealso [fd_known()] for the declarations; [fd_weights()] and
#'   [fd_offsets()] for the one-dimensional rules; [fd_derivative()] for one
#'   coordinate.
#'
#' @examples
#' f <- function(x) x$a^3 * exp(x$b)
#' x <- list(a = c(0.5, 1, 2), b = 0.3)
#' d3 <- fd_tensor(f, x, order = 3)
#' # the exact a_a_b component is 6 a exp(b)
#' rbind(d3$a_a_b, 6 * x$a * exp(x$b))
#'
#' # with the analytic gradient, each third derivative is a second difference
#' g <- fd_known(function(x) list(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b)),
#'               order = 1)
#' d3g <- fd_tensor(f, x, order = 3, known = g)
#' attr(d3g, "fd_source")
#' rbind(d3g$a_a_b, 6 * x$a * exp(x$b))
#'
#' # one component only, at second-order accuracy
#' fd_tensor(f, x, components = rbind(c(2, 1)), accuracy = 2)
#' @export
fd_tensor <- function(f, x, order = NULL, components = NULL, known = list(),
                      bounds = NULL, accuracy = 4L, h_rel = NULL) {
  if (!is.function(f)) stop("'f' must be a function.", call. = FALSE)
  nm <- names(x)
  if (!is.list(x) || !length(x) || is.null(nm) || anyNA(nm) ||
      any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("'x' must be a non-empty list of coordinates with distinct names.",
         call. = FALSE)
  }
  if (!all(vapply(x, is.numeric, logical(1)))) {
    stop("Every coordinate of 'x' must be numeric.", call. = FALSE)
  }
  p <- length(x)
  len <- lengths(x)
  n <- max(len)
  if (any(len != 1L & len != n)) {
    stop("Every coordinate of 'x' must have length 1 or the common length n.",
         call. = FALSE)
  }
  if (!is_whole(accuracy) || accuracy < 2 || accuracy %% 2 != 0) {
    stop("'accuracy' must be a positive even whole number.", call. = FALSE)
  }
  accuracy <- as.integer(accuracy)
  if (!is.null(h_rel) && (!is.numeric(h_rel) || length(h_rel) != 1L ||
                          !is.finite(h_rel) || h_rel <= 0)) {
    stop("'h_rel' must be NULL or a single positive number.", call. = FALSE)
  }

  if (is.null(order) == is.null(components)) {
    stop("Give exactly one of 'order' and 'components'.", call. = FALSE)
  }
  if (!is.null(order)) {
    if (!is_whole(order) || order < 1) {
      stop("'order' must be a single positive whole number.", call. = FALSE)
    }
    comp <- counts_of_order(p, as.integer(order))
    comp_key <- paste0("o", order)
  } else {
    comp <- as_counts(components, p, "components")
    comp_key <- paste(comp, collapse = ",")
  }
  comp_names <- apply(comp, 1L, function(a) paste(rep(nm, a), collapse = "_"))

  if (S7::S7_inherits(known, fd_known)) known <- list(known)
  if (!is.list(known) ||
      !all(vapply(known, S7::S7_inherits, logical(1), class = fd_known))) {
    stop("'known' must be a list of fd_known() objects.", call. = FALSE)
  }
  known_key <- vapply(known, function(k) {
    if (!is.null(k@order)) paste0("o", k@order) else
      paste(c(nrow(k@indices), k@indices), collapse = ",")
  }, character(1))
  kidx <- c(list(matrix(0L, 1L, p)), lapply(known, function(k) {
    if (!is.null(k@order)) return(counts_of_order(p, k@order))
    if (ncol(k@indices) != p) {
      stop(sprintf("The indices of an fd_known() object have %d columns for %d coordinates.",
                   ncol(k@indices), p), call. = FALSE)
    }
    k@indices
  }))
  funs <- c(list(f), lapply(known, function(k) k@fun))

  bnd <- vector("list", p)
  if (!is.null(bounds)) {
    if (!is.list(bounds)) stop("'bounds' must be NULL or a list.", call. = FALSE)
    if (!is.null(names(bounds))) {
      bnd <- bounds[nm]
    } else {
      if (length(bounds) != p) {
        stop("An unnamed 'bounds' must have one entry per coordinate.",
             call. = FALSE)
      }
      bnd <- bounds
    }
    for (b in bnd) {
      if (!is.null(b) && (!is.numeric(b) || length(b) != 2L || anyNA(b) ||
                          b[1] >= b[2])) {
        stop("Each entry of 'bounds' must be NULL or c(lower, upper) with lower < upper.",
             call. = FALSE)
      }
    }
  }

  plan <- tensor_plan(paste(c(accuracy, p, comp_key, known_key),
                            collapse = "|"), comp, kidx, accuracy)

  # one step per residual order and direction, reduced to keep the widest
  # stencil of that order inside the bounds
  steps <- lapply(seq_along(plan$r_val), function(k) {
    hr <- if (is.null(h_rel)) {
      .Machine$double.eps^(1 / (plan$r_val[k] + accuracy))
    } else h_rel
    lapply(seq_len(p), function(j) {
      reach <- plan$r_reach[k, j]
      if (reach == 0L) return(NULL)
      xj <- x[[j]]
      h <- hr * pmax(1, abs(xj))
      b <- bnd[[j]]
      if (!is.null(b)) {
        if (is.finite(b[1])) h <- pmin(h, 0.49 * (xj - b[1]) / reach)
        if (is.finite(b[2])) h <- pmin(h, 0.49 * (b[2] - xj) / reach)
      }
      h
    })
  })

  # the (node, column) pairs the terms read, each evaluated once
  m <- nrow(comp)
  t_comp <- rep.int(seq_len(m), diff(plan$term_start))
  t_node <- plan$term_node + 1L
  t_col <- plan$comp_col[t_comp] + 1L
  width <- max(vapply(kidx, nrow, integer(1)))
  key <- (t_node - 1L) * width + t_col
  ukey <- unique(key)
  t_val <- match(key, ukey)
  u_node <- (ukey - 1L) %/% width + 1L
  u_col <- (ukey - 1L) %% width + 1L

  vals <- vector("list", length(ukey))
  for (i in unique(u_node)) {
    src <- plan$node_src[i] + 1L
    off <- plan$node_off[i, ]
    xi <- x
    if (any(off != 0L)) {
      st <- steps[[match(plan$node_r[i], plan$r_val)]]
      for (j in which(off != 0L)) xi[[j]] <- x[[j]] + off[j] * st[[j]]
    }
    v <- as_columns(funs[[src]](xi), nrow(kidx[[src]]), src)
    for (u in which(u_node == i)) vals[[u]] <- v[[u_col[u]]]
  }
  # the values may be longer than the coordinates, when the vector the
  # derivatives run along is held by f rather than passed in x
  vlen <- lengths(vals)
  n_out <- max(n, vlen)
  if (any(vlen != 1L & vlen != n_out)) {
    stop("The values of 'f' and of the known quantities must have length 1 or a common length.",
         call. = FALSE)
  }

  dens <- lapply(seq_len(m), function(c) {
    st <- steps[[match(plan$comp_r[c], plan$r_val)]]
    den <- 1
    for (j in which(plan$comp_rem[c, ] > 0L)) {
      den <- den * st[[j]]^plan$comp_rem[c, j]
    }
    as.double(den)
  })

  out <- fd_tensor_combine_cpp(vals, plan$term_start, t_val - 1L,
                               plan$term_w, dens, n_out, 1L)
  names(out) <- comp_names
  attr(out, "fd_source") <- structure(
    list(component = unname(comp_names), source = plan$comp_src,
         index = plan$comp_col + 1L, residual = plan$comp_r),
    class = "data.frame", row.names = c(NA_integer_, -m))
  out
}

# The plans of the session, keyed by the integer structure of the request.
.fd_tensor_plans <- new.env(parent = emptyenv())

#' Helpers of the Tensor Stencil
#'
#' @description
#' These are internal functions of [fd_tensor()]. `tensor_plan()` returns the
#' plan of a
#' request from the session's cache, computing it in C++ on the first call
#' with that key. `is_whole()` tests for a single finite whole number.
#' `counts_of_order()` returns the multi-indices of one total order as rows of
#' counts, in the order of [tuple_indices()]. `as_counts()` validates a matrix
#' of multi-indices and stores it as integers. `as_columns()` turns the value
#' of a source function into a list of double vectors, one per declared
#' multi-index.
#'
#' @param key A string identifying the request: the accuracy, the number of
#'   coordinates, the components and the declarations.
#' @param comp,m An integer matrix of counts, one row per multi-index.
#' @param kidx A list of such matrices, of which the first is the function's
#'   own zero row.
#' @param accuracy The order of accuracy of the stencils.
#' @param x A value to test.
#' @param p,k The number of coordinates and a total order; in `as_columns()`,
#'   `k` is the number of declared multi-indices.
#' @param what The argument name used in an error message.
#' @param v The value returned by a source function.
#' @param src The position of the source: 1 for `f`, \eqn{j + 1} for
#'   `known[[j]]`.
#'
#' @return `tensor_plan()` a list of the plan's vectors and matrices;
#'   `is_whole()` a single logical; `counts_of_order()` and `as_counts()` an
#'   integer matrix; `as_columns()` a list of double vectors.
#' @keywords internal
#' @name fd_tensor_helpers
NULL

#' @rdname fd_tensor_helpers
tensor_plan <- function(key, comp, kidx, accuracy) {
  plan <- .fd_tensor_plans[[key]]
  if (is.null(plan)) {
    plan <- fd_tensor_plan_cpp(comp, kidx, accuracy)
    assign(key, plan, envir = .fd_tensor_plans)
  }
  plan
}

#' @rdname fd_tensor_helpers
is_whole <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == trunc(x)
}

#' @rdname fd_tensor_helpers
counts_of_order <- function(p, k) {
  if (k == 0L) return(matrix(0L, 1L, p))
  tup <- tuple_indices(p, k)
  matrix(unlist(lapply(tup, tabulate, nbins = p)), ncol = p, byrow = TRUE)
}

#' @rdname fd_tensor_helpers
as_counts <- function(m, p, what) {
  if (is.numeric(m) && is.null(dim(m))) m <- matrix(m, nrow = 1L)
  if (!is.numeric(m) || !is.matrix(m) || !nrow(m) || anyNA(m) ||
      any(m < 0) || any(m != trunc(m))) {
    stop(sprintf("'%s' must be a matrix of non-negative whole counts.", what),
         call. = FALSE)
  }
  if (!is.null(p) && ncol(m) != p) {
    stop(sprintf("'%s' must have one column per coordinate (%d).", what, p),
         call. = FALSE)
  }
  storage.mode(m) <- "integer"
  dimnames(m) <- NULL
  m
}

#' @rdname fd_tensor_helpers
as_columns <- function(v, k, src) {
  what <- function() {
    if (src == 1L) "'f'" else sprintf("The function of known[[%d]]", src - 1L)
  }
  if (is.matrix(v)) {
    v <- lapply(seq_len(ncol(v)), function(j) v[, j])
  } else if (!is.list(v)) {
    v <- list(v)
  }
  if (length(v) != k) {
    stop(sprintf("%s returned %d values for %d declared derivatives.",
                 what(), length(v), k), call. = FALSE)
  }
  lapply(v, function(e) {
    if (!is.numeric(e) || !length(e)) {
      stop(sprintf("%s must return non-empty numeric vectors.", what()),
           call. = FALSE)
    }
    as.double(e)
  })
}
