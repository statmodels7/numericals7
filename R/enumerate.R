# The combinatorial enumerations a higher-order chain rule rests on. Each of
# these existed somewhere in the toolkit before this package did -- two of
# them twice -- and the whole point of collecting them is that an enumeration
# living in one copy cannot disagree with itself.

#' The Index Tuples of a Given Width
#'
#' @description
#' Every multi-index of a derivative of the given order over `d`
#' variables: the tuples \eqn{(i_1 \le \dots \le i_k)} that key a list of
#' partial derivatives.
#'
#' @details
#' The ordering is part of the contract. At order two the diagonal comes
#' first and the off-diagonal pairs follow, because that is the order a
#' Hessian consumer indexes by; at orders three and four the enumeration is
#' the lexicographic one over non-decreasing tuples. Anything that holds
#' derivatives over \eqn{d} variables -- a constrained parameter, say --
#' keys its components by this enumeration, so nothing has to be reordered
#' when they meet.
#'
#' @param d The number of variables.
#' @param order The derivative order, 1 to 4.
#'
#' @return A list of integer vectors of length `order`.
#'
#' @seealso [set_partitions()]
#'
#' @examples
#' tuple_indices(2, 2)
#' lengths(tuple_indices(3, 4))[1:3]
#' # there are choose(d + k - 1, k) tuples of width k
#' length(tuple_indices(3, 4)) == choose(3 + 4 - 1, 4)
#'
#' @export
tuple_indices <- function(d, order = 2L) {
  if (!order %in% 1:4) stop("'order' must be 1, 2, 3 or 4.", call. = FALSE)
  if (d == 0L) return(list())
  if (order == 1L) return(lapply(seq_len(d), as.integer))
  if (order == 2L) {
    # Diagonal first: the ordering is part of the contract, because it is the
    # one a Hessian consumer indexes by.
    out <- lapply(seq_len(d), function(k) c(k, k))
    if (d > 1L) {
      for (k in seq_len(d - 1L)) {
        for (l in seq.int(k + 1L, d)) out[[length(out) + 1L]] <- c(k, l)
      }
    }
    return(lapply(out, as.integer))
  }
  grid <- utils::combn(seq_len(d + order - 1L), order)
  lapply(seq_len(ncol(grid)), function(j) {
    as.integer(grid[, j] - seq_len(order) + 1L)
  })
}


#' The Set Partitions of the First n Integers
#'
#' @description
#' Every way of splitting `1:n` into disjoint non-empty blocks, which is
#' what a chain rule of order \eqn{n} sums over.
#'
#' @details
#' Built by the standard recursion: the partitions of `1:n` are obtained
#' from those of `1:(n-1)` by placing `n` into each existing block
#' in turn and then into a block of its own. There are 1, 2, 5 and 15 of them
#' for \eqn{n = 1, \dots, 4}, the Bell numbers.
#'
#' The blocks index **positions** rather than variables, which is what
#' makes a repeated variable count with the right multiplicity without any
#' bookkeeping of its own -- the same subset-of-positions device with
#' subsets of positions.
#'
#' @param n A positive integer, at most four here.
#'
#' @return A list of partitions, each a list of integer vectors.
#'
#' @seealso [tuple_indices()]
#'
#' @examples
#' set_partitions(3)
#' lengths(lapply(1:4, set_partitions))
#'
#' @export
set_partitions <- function(n) {
  if (n == 1L) return(list(list(1L)))
  prev <- set_partitions(n - 1L)
  out <- list()
  for (p in prev) {
    for (k in seq_along(p)) {
      q <- p
      q[[k]] <- c(q[[k]], n)
      out[[length(out) + 1L]] <- q
    }
    out[[length(out) + 1L]] <- c(p, list(n))
  }
  out
}


#' Every Way to Write an Integer as an Ordered Sum
#'
#' @description
#' The weak compositions of `n` into `k` parts: every vector of
#' `k` non-negative integers summing to `n`, one per row.
#'
#' @details
#' The set enumerated is
#'
#' \deqn{\mathcal{C}(n, k) = \Bigl\{(c_1, \dots, c_k) \in \mathbb{N}_0^{k} :
#'   \textstyle\sum_{j=1}^{k} c_j = n\Bigr\},
#'   \qquad \lvert\mathcal{C}(n, k)\rvert = \binom{n + k - 1}{k - 1},}
#'
#' the count being the number of ways to place \eqn{k - 1} dividers among
#' \eqn{n} units.
#'
#' Built by recursion on the number of parts, which is what keeps the result
#' in a fixed order and avoids generating and filtering a full grid. There
#' are `choose(n + k - 1, k - 1)` of them, so the enumeration is only
#' practical for a moderate size: at `n = 20` and `k = 5` it is
#' 10626 rows, and at `k = 10` it is 10015005.
#'
#' A discrete distribution on a fixed total -- a multinomial -- has exactly
#' this set as its support, which is what makes its expectations exact sums.
#'
#' @param n The total, a non-negative integer.
#' @param k The number of parts, a positive integer.
#'
#' @return An integer matrix with `k` columns.
#'
#' @seealso [set_partitions()]
#'
#' @examples
#' compositions(3, 2)
#' nrow(compositions(5, 3)) == choose(5 + 3 - 1, 3 - 1)
#'
#' @export
compositions <- function(n, k) {
  n <- as.integer(n)
  k <- as.integer(k)
  if (k == 1L) return(matrix(n, nrow = 1L, ncol = 1L))
  out <- lapply(0:n, function(first) {
    rest <- compositions(n - first, k - 1L)
    cbind(first, rest, deparse.level = 0L)
  })
  do.call(rbind, out)
}
