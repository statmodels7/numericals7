# The combinatorial enumerations a higher-order chain rule rests on. Each of
# these existed somewhere in the toolkit before this package did -- two of
# them twice -- and the whole point of collecting them is that an enumeration
# living in one copy cannot disagree with itself.

#' Multi-Indices for Derivatives of a Given Order
#'
#' @description
#' Enumerates the multi-indices that key a list of partial derivatives: every
#' non-decreasing tuple \eqn{(i_1 \le \cdots \le i_k)} drawn from `1:d`, where
#' \eqn{k} is the derivative order and \eqn{d} the number of variables. There
#' is one tuple per distinct partial derivative, because a mixed partial does
#' not depend on the order the variables are differentiated in, so the count is
#' \eqn{\binom{d + k - 1}{k}}.
#'
#' @details
#' # The ordering is part of the interface
#'
#' At order two the diagonal comes first, \eqn{(1,1), (2,2), \dots, (d,d)}, and
#' the off-diagonal pairs follow in lexicographic order. That is the order a
#' Hessian consumer indexes by. At orders three and four the enumeration is
#' plain lexicographic over non-decreasing tuples.
#'
#' Every object in the toolkit holding derivatives over \eqn{d} variables keys
#' its components by this enumeration, so two of them meet without either being
#' reordered. Treat the order as fixed.
#'
#' # Counts
#'
#' \deqn{\lvert T(d, k)\rvert = \binom{d + k - 1}{k},}
#'
#' the number of multisets of size \eqn{k} drawn from \eqn{d} symbols. At
#' \eqn{d = 3} the four orders give 3, 6, 10 and 15 tuples.
#'
#' @param d The number of variables, a non-negative whole number. Zero gives an
#'   empty list. A negative value throws, from `seq_len()`, with a message
#'   about coercion to a non-negative integer.
#' @param order The derivative order \eqn{k}, one of `1`, `2`, `3` or `4`.
#'   Anything else throws: the enumeration is written out only that far,
#'   because a fourth derivative is as high as the toolkit carries.
#'
#' @return A list of \eqn{\binom{d + k - 1}{k}} integer vectors, each of length
#'   `order`, each non-decreasing. An empty list when `d` is zero.
#'
#' @seealso [set_partitions()], the other enumeration a higher-order chain rule
#'   needs, and [compositions()] for the ordered sums.
#'
#' @examples
#' # Second-order tuples of three variables: the diagonal first, then the
#' # off-diagonal pairs. A Hessian consumer indexes in this order.
#' tuple_indices(3, 2)
#'
#' # One tuple per distinct partial derivative, so the count is the number of
#' # multisets of size k drawn from d symbols.
#' vapply(1:4, function(k) length(tuple_indices(3, k)), integer(1))
#' choose(3 + 1:4 - 1, 1:4)
#'
#' # Every tuple is non-decreasing, so no derivative is enumerated twice.
#' all(vapply(tuple_indices(4, 3), function(i) !is.unsorted(i), logical(1)))
#'
#' # With no variables there are no derivatives.
#' tuple_indices(0, 2)
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


#' Set Partitions of the First n Integers
#'
#' @description
#' Enumerates every way of splitting `1:n` into disjoint non-empty blocks. A
#' chain rule of order \eqn{n} contributes one term per partition, so this is
#' the index set of Faà di Bruno's formula and of the Bell-polynomial
#' identities built on it. The number of partitions is the Bell number
#' \eqn{B_n}: 1, 2, 5, 15, 52 and 203 for \eqn{n} from one to six.
#'
#' @details
#' # How they are built
#'
#' By the standard recursion. The partitions of `1:n` come from those of
#' `1:(n-1)` by putting \eqn{n} into each existing block in turn, and then into
#' a block of its own.
#'
#' The cost is therefore \eqn{B_n}, which grows faster than any exponential:
#' \eqn{B_8} is 4140 and \eqn{B_{10}} is 115975. Four is as high as the
#' toolkit's derivatives go, where the sum has fifteen terms.
#'
#' # The blocks index positions
#'
#' A block holds positions within a multi-index, not variables. A variable
#' repeated in that multi-index therefore counts with its correct multiplicity
#' and needs no bookkeeping of its own: differentiating three times in one
#' variable and once in another sums over the same fifteen partitions as four
#' distinct variables do, and what differs is which derivative each block
#' names. [tuple_indices()] supplies the multi-index the positions point into.
#'
#' @param n A positive whole number. The recursion has no base case below one,
#'   so zero and negative values recurse until the stack overflows instead of
#'   throwing.
#'
#' @return A list of \eqn{B_n} partitions. Each partition is a list of integer
#'   vectors, its blocks, which between them contain each of `1:n` exactly
#'   once. The storage is integer whether `n` is given as `4` or as `4L`, as
#'   it is for [tuple_indices()] and [compositions()], so the three
#'   enumerations agree and a block may be compared against `1:n` with
#'   `identical()`.
#'
#' @references
#' Constantine, G. M. and Savits, T. H. (1996). A multivariate Faà di Bruno
#' formula with applications. *Transactions of the American Mathematical
#' Society* **348**, 503-520.
#'
#' @seealso [tuple_indices()] for the multi-indices the blocks index into, and
#'   [compositions()] for the ordered sums.
#'
#' @examples
#' # The five partitions of {1, 2, 3}.
#' set_partitions(3)
#'
#' # The counts are the Bell numbers.
#' lengths(lapply(1:6, set_partitions))
#'
#' # Every partition covers 1:n exactly once, the blocks being disjoint. The
#' # comparison is `identical()` because the blocks are integer however `n` was
#' # given, which is the convention all three enumerations follow.
#' all(vapply(set_partitions(4), function(p) identical(sort(unlist(p)), 1:4),
#'            logical(1)))
#'
#' # A fourth-order chain rule sums fifteen terms, one per partition, and the
#' # number of blocks in each is the number of factors in that term.
#' table(lengths(set_partitions(4)))
#'
#' @export
set_partitions <- function(n) {
  if (n == 1L) return(list(list(1L)))
  prev <- set_partitions(n - 1L)
  # The two places n enters a block, coerced so that a block's storage does not
  # follow the caller's: `tuple_indices()` and `compositions()` both return
  # integers whatever they are given, and an enumeration kept in one copy so it
  # cannot disagree with itself should not disagree about its return type. The
  # coercion is here rather than on `n` at the top so that it changes the
  # storage and nothing else -- coercing at the top would turn a fractional
  # argument from the documented stack overflow into a silent truncation.
  ni <- as.integer(n)
  out <- list()
  for (p in prev) {
    for (k in seq_along(p)) {
      q <- p
      q[[k]] <- c(q[[k]], ni)
      out[[length(out) + 1L]] <- q
    }
    out[[length(out) + 1L]] <- c(p, list(ni))
  }
  out
}


#' Weak Compositions of an Integer
#'
#' @description
#' Enumerates the weak compositions of \eqn{n} into \eqn{k} parts: every vector
#' of \eqn{k} non-negative whole numbers summing to \eqn{n}, returned one per
#' row of an integer matrix. *Weak* means a part may be zero, so the set is the
#' whole lattice simplex of counts and includes its faces. There are
#' \eqn{\binom{n + k - 1}{k - 1}} rows.
#'
#' @details
#' # The set
#'
#' \deqn{\mathcal{C}(n, k) = \Bigl\{(c_1, \dots, c_k) \in \mathbb{N}_0^{k} :
#'   \textstyle\sum_{j=1}^{k} c_j = n\Bigr\},
#'   \qquad \lvert\mathcal{C}(n, k)\rvert = \binom{n + k - 1}{k - 1}.}
#'
#' The count is the number of ways to place \eqn{k - 1} dividers among \eqn{n}
#' units.
#'
#' # Order and cost
#'
#' Rows come out in lexicographic order, the first part ascending from zero to
#' \eqn{n}. The recursion runs on the number of parts, so the rows are produced
#' in that order directly and no full grid is built and filtered.
#'
#' The whole matrix is materialized and the count grows quickly, which bounds
#' the practical size: \eqn{n = 20} with \eqn{k = 5} is 10626 rows, and the
#' same \eqn{n} with \eqn{k = 10} is 10015005.
#'
#' # Where the set turns up
#'
#' It is exactly the support of a multinomial with total \eqn{n} over \eqn{k}
#' categories, and of any other distribution on a fixed total. An expectation
#' under such a law is therefore a finite sum over these rows, evaluated
#' exactly, where a continuous family would need a quadrature.
#'
#' @param n The total, a non-negative whole number. Neither argument is
#'   validated. A negative `n` returns a matrix that answers no question, since
#'   the recursion walks `0:n` and that sequence runs downwards.
#' @param k The number of parts, a positive whole number. Zero recurses until
#'   the stack overflows.
#'
#' @return An integer matrix with `k` columns and
#'   \eqn{\binom{n + k - 1}{k - 1}} rows, every row summing to `n`.
#'
#' @seealso [set_partitions()] and [tuple_indices()], the other two
#'   enumerations collected here.
#'
#' @examples
#' # The four ways to split 3 into two ordered non-negative parts.
#' compositions(3, 2)
#'
#' # Every row sums to the total, and the row count is the divider formula.
#' cs <- compositions(5, 3)
#' all(rowSums(cs) == 5)
#' nrow(cs) == choose(5 + 3 - 1, 3 - 1)
#'
#' # Zero parts are allowed, so the faces of the simplex are included. Of the
#' # fifteen rows, twelve touch a face and three are strictly positive.
#' cs4 <- compositions(4, 3)
#' c(rows = nrow(cs4),
#'   with_a_zero = sum(apply(cs4, 1, function(r) any(r == 0))),
#'   positive = choose(4 - 1, 3 - 1))
#'
#' # One part, or a total of zero, leaves a single composition.
#' compositions(3, 1)
#' compositions(0, 3)
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
