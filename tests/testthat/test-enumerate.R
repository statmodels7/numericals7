# The enumerations. Each existed elsewhere in the toolkit before this package
# did, so the tests carry the assertions that came with them.

test_that("tuple_indices counts and orders as the contract says", {
  # choose(d + k - 1, k) tuples of width k over d variables
  for (d in 1:4) {
    for (o in 1:4) {
      tt <- tuple_indices(d, o)
      expect_length(tt, choose(d + o - 1, o))
      expect_true(all(lengths(tt) == o), label = paste(d, o))
      # non-decreasing within each tuple
      expect_true(all(vapply(tt, function(t) !is.unsorted(t), logical(1))))
      expect_identical(anyDuplicated(tt), 0L)
    }
  }

  # the order-2 ordering is a CONTRACT: diagonal first, then off-diagonal.
  # A Hessian consumer indexes by it, and pairing it with a lexicographic
  # enumeration mislabels components -- the mu_sigma/(sigma, sigma) trap.
  t2 <- tuple_indices(3, 2)
  expect_identical(t2[1:3], list(c(1L, 1L), c(2L, 2L), c(3L, 3L)))
  expect_identical(t2[4:6], list(c(1L, 2L), c(1L, 3L), c(2L, 3L)))

  expect_identical(tuple_indices(0, 2), list())
  expect_error(tuple_indices(2, 5), "1, 2, 3 or 4")
})


test_that("set_partitions has the Bell numbers and disjoint blocks", {
  expect_identical(lengths(lapply(1:4, set_partitions)), c(1L, 2L, 5L, 15L))
  for (n in 1:4) {
    for (p in set_partitions(n)) {
      expect_setequal(sort(unlist(p)), seq_len(n))
      expect_identical(length(unlist(p)), n)
    }
  }
})


test_that("the three enumerations agree on storage, whatever they are given", {
  # An enumeration kept in one copy so it cannot disagree with itself should
  # not disagree about its return type either. The blocks are integer whether
  # n arrives as a double or as an integer, so the comparison a caller reaches
  # for -- identical() against 1:n -- holds. Before this was pinned,
  # set_partitions(4) returned doubles and set_partitions(4L) integers.
  expect_identical(set_partitions(4), set_partitions(4L))
  for (n in 1:5) {
    p <- set_partitions(n)
    # `info` because an expectation that prints only TRUE or FALSE cannot be
    # read on a platform one does not have.
    expect_true(all(vapply(p, function(q) all(vapply(q, is.integer, TRUE)),
                           TRUE)),
                info = paste("a block is not integer at n =", n))
    expect_true(all(vapply(p, function(q) identical(sort(unlist(q)),
                                                    seq_len(n)), TRUE)),
                info = paste("a partition does not cover 1:n at n =", n))
  }

  expect_identical(tuple_indices(3, 2), tuple_indices(3L, 2L))
  expect_identical(compositions(3, 2), compositions(3L, 2L))
  expect_identical(typeof(unlist(tuple_indices(3, 4))), "integer")
  expect_identical(typeof(compositions(3, 2)), "integer")
})


test_that("compositions enumerates every way to split an integer", {
  expect_identical(compositions(3, 2),
                   matrix(c(0L, 1L, 2L, 3L, 3L, 2L, 1L, 0L), ncol = 2))
  for (n in 0:5) {
    for (k in 1:4) {
      m <- compositions(n, k)
      expect_identical(nrow(m), as.integer(choose(n + k - 1, k - 1)))
      expect_true(all(rowSums(m) == n))
      expect_true(all(m >= 0))
      expect_identical(anyDuplicated(m), 0L)
    }
  }
  # one part: the whole total in one cell
  expect_identical(compositions(4, 1), matrix(4L, 1, 1))
})
