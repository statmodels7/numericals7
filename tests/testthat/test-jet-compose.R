# The jet vocabulary. A jet carries a value together with every partial
# derivative to fourth order, so a quantity built from it is analytic at every
# order with no formula transcribed.
#
# The reference is a HAND-WRITTEN exact derivative, not a numerical one. This
# is the package's own rule and it is not a preference: nesting
# numDeriv::grad four deep to reach a fourth derivative returns 41703.95 for
# exp(x^3) at x = 0.9, where the exact value is 771.13. The reference would be
# wrong by a factor of fifty-four, and it is the reference that fails, not the
# jets.

jet_seed <- function(x) {
  lay <- jet_layout(length(x))
  list(lay = lay,
       v = lapply(seq_along(x), function(k) {
         jet_var(k, list(x[k], 1, 0, 0, 0), lay)
       }))
}

jet_at <- function(j, lay, tup) {
  slot <- get(paste(sort(tup), collapse = ","), envir = lay$pos)
  j$d[[slot[1L]]][slot[2L]]
}


test_that("set_partitions has the Bell numbers and disjoint blocks", {
  expect_length(set_partitions(1L), 1L)
  expect_length(set_partitions(2L), 2L)
  expect_length(set_partitions(3L), 5L)
  expect_length(set_partitions(4L), 15L)
  for (n in 1:4) {
    for (p in set_partitions(n)) {
      expect_setequal(sort(unlist(p)), seq_len(n))
      # disjoint: the blocks together use each position exactly once
      expect_identical(length(unlist(p)), n)
    }
  }
})


test_that("the scalar vocabulary is exact against written derivatives", {
  lay <- jet_layout(1L)
  x <- 0.7
  a <- jet_var(1L, list(x, 1, 0, 0, 0), lay)

  e <- jet_exp(a, lay)
  expect_equal(e$v, exp(x))
  for (o in 1:4) expect_equal(e$d[[o]][1], exp(x))

  l <- jet_log(a, lay)
  expect_equal(l$v, log(x))
  expect_equal(l$d[[1]][1], 1 / x)
  expect_equal(l$d[[2]][1], -1 / x^2)
  expect_equal(l$d[[3]][1], 2 / x^3)
  expect_equal(l$d[[4]][1], -6 / x^4)

  p <- jet_pow(a, 2.5, lay)
  expect_equal(p$d[[1]][1], 2.5 * x^1.5)
  expect_equal(p$d[[2]][1], 2.5 * 1.5 * x^0.5)
  expect_equal(p$d[[3]][1], 2.5 * 1.5 * 0.5 * x^-0.5)
  expect_equal(p$d[[4]][1], 2.5 * 1.5 * 0.5 * -0.5 * x^-1.5)

  expect_equal(jet_sqrt(a, lay)$v, sqrt(x))
  expect_equal(jet_inv(a, lay)$v, 1 / x)

  g <- jet_lgamma(a, lay)
  expect_equal(g$v, lgamma(x))
  expect_equal(g$d[[1]][1], digamma(x))
  expect_equal(g$d[[2]][1], trigamma(x))
  expect_equal(g$d[[3]][1], psigamma(x, 2L))
  expect_equal(g$d[[4]][1], psigamma(x, 3L))

  dg <- jet_digamma(a, lay)
  expect_equal(dg$v, digamma(x))
  expect_equal(dg$d[[1]][1], trigamma(x))
  expect_equal(dg$d[[2]][1], psigamma(x, 2L))

  # gamma is exp of lgamma, so it must reproduce the direct derivatives
  gm <- jet_gamma(a, lay)
  expect_equal(gm$v, gamma(x))
  expect_equal(gm$d[[1]][1], gamma(x) * digamma(x))
  expect_equal(gm$d[[2]][1], gamma(x) * (digamma(x)^2 + trigamma(x)))
})


test_that("a composition is exact at every order", {
  # exp(x^3): the chain rule is genuinely non-trivial at every order and the
  # four derivatives are short enough to write out.
  lay <- jet_layout(1L)
  x <- 0.9
  a <- jet_var(1L, list(x, 1, 0, 0, 0), lay)
  j <- jet_exp(jet_pow(a, 3, lay), lay)
  e <- exp(x^3)

  expect_equal(j$v, e)
  expect_equal(j$d[[1]][1], 3 * x^2 * e)
  expect_equal(j$d[[2]][1], (9 * x^4 + 6 * x) * e)
  expect_equal(j$d[[3]][1], (27 * x^6 + 54 * x^3 + 6) * e)
  expect_equal(j$d[[4]][1], (81 * x^8 + 324 * x^5 + 180 * x^2) * e)
})


test_that("every mixed component of a product is exact", {
  # x1^2 x2^3 has all nine non-zero partial derivatives writable, and they are
  # what a chain rule gets wrong when it counts a repeated index once instead
  # of with its multiplicity.
  s <- jet_seed(c(1.3, 0.8))
  X <- 1.3
  Y <- 0.8
  j <- jet_mul(jet_pow(s$v[[1]], 2, s$lay), jet_pow(s$v[[2]], 3, s$lay), s$lay)

  expect_equal(j$v, X^2 * Y^3)
  expect_equal(jet_at(j, s$lay, 1L), 2 * X * Y^3)
  expect_equal(jet_at(j, s$lay, 2L), 3 * X^2 * Y^2)
  expect_equal(jet_at(j, s$lay, c(1L, 1L)), 2 * Y^3)
  expect_equal(jet_at(j, s$lay, c(1L, 2L)), 6 * X * Y^2)
  expect_equal(jet_at(j, s$lay, c(2L, 2L)), 6 * X^2 * Y)
  expect_equal(jet_at(j, s$lay, c(1L, 1L, 2L)), 6 * Y^2)
  expect_equal(jet_at(j, s$lay, c(1L, 2L, 2L)), 12 * X * Y)
  expect_equal(jet_at(j, s$lay, c(2L, 2L, 2L)), 6 * X^2)
  expect_equal(jet_at(j, s$lay, c(1L, 1L, 2L, 2L)), 12 * Y)
  # and the components that must vanish do
  expect_equal(jet_at(j, s$lay, c(1L, 1L, 1L)), 0)
  expect_equal(jet_at(j, s$lay, c(2L, 2L, 2L, 2L)), 0)
})


test_that("a transcendental of several variables is exact", {
  # log(x1 + 2 x2) composed with a product, so that no component separates.
  s <- jet_seed(c(1.3, 0.8))
  X <- 1.3
  Y <- 0.8
  u <- X + 2 * Y
  j <- jet_log(jet_add(s$v[[1]], jet_cmul(s$v[[2]], 2)), s$lay)

  expect_equal(j$v, log(u))
  expect_equal(jet_at(j, s$lay, 1L), 1 / u)
  expect_equal(jet_at(j, s$lay, 2L), 2 / u)
  expect_equal(jet_at(j, s$lay, c(1L, 2L)), -2 / u^2)
  expect_equal(jet_at(j, s$lay, c(2L, 2L)), -4 / u^2)
  expect_equal(jet_at(j, s$lay, c(2L, 2L, 2L)), 16 / u^3)
  expect_equal(jet_at(j, s$lay, c(1L, 2L, 2L, 2L)), -48 / u^4)
})


test_that("jet_compose agrees with jet_mul on a square", {
  # Two routes to one object, sharing the value and nothing else: jet_mul
  # enumerates subsets of positions, jet_compose set partitions of them.
  lay <- jet_layout(2L)
  a <- jet_var(1L, list(1.4, 1, 0, 0, 0), lay)
  b <- jet_var(2L, list(0.6, 1, 0, 0, 0), lay)
  s <- jet_add(a, jet_cmul(b, 3))
  by_mul <- jet_mul(s, s, lay)
  by_pow <- jet_pow(s, 2, lay)
  expect_equal(by_pow$v, by_mul$v)
  for (o in 1:4) expect_equal(by_pow$d[[o]], by_mul$d[[o]])

  # and a cube, where the two enumerations disagree more if either is wrong
  cube_mul <- jet_mul(by_mul, s, lay)
  cube_pow <- jet_pow(s, 3, lay)
  for (o in 1:4) expect_equal(cube_pow$d[[o]], cube_mul$d[[o]])
})


test_that("jet_div is the product with the reciprocal", {
  s <- jet_seed(c(2.1, 0.7))
  q <- jet_div(s$v[[1]], s$v[[2]], s$lay)
  X <- 2.1
  Y <- 0.7
  expect_equal(q$v, X / Y)
  expect_equal(jet_at(q, s$lay, 1L), 1 / Y)
  expect_equal(jet_at(q, s$lay, 2L), -X / Y^2)
  expect_equal(jet_at(q, s$lay, c(2L, 2L)), 2 * X / Y^3)
  expect_equal(jet_at(q, s$lay, c(1L, 2L, 2L)), 2 / Y^3)
  expect_equal(jet_at(q, s$lay, c(2L, 2L, 2L, 2L)), 24 * X / Y^5)
})


test_that("one numerical differentiation of an exact gradient agrees", {
  # The one legal use of a numerical reference: ONE Richardson pass applied to
  # a gradient written by hand, never a pass applied to another pass.
  skip_if_not_installed("numDeriv")
  f <- function(z) sqrt(z[1] * z[2]) * log(z[1] + 2 * z[2])
  grad_exact <- function(z) {
    r <- sqrt(z[1] * z[2])
    l <- log(z[1] + 2 * z[2])
    u <- z[1] + 2 * z[2]
    c(z[2] / (2 * r) * l + r / u, z[1] / (2 * r) * l + 2 * r / u)
  }
  x <- c(1.3, 0.8)

  build <- function(v, lay) {
    jet_mul(jet_sqrt(jet_mul(v[[1]], v[[2]], lay), lay),
            jet_log(jet_add(v[[1]], jet_cmul(v[[2]], 2)), lay), lay)
  }
  s <- jet_seed(x)
  j <- build(s$v, s$lay)

  expect_equal(j$v, f(x))
  expect_equal(c(jet_at(j, s$lay, 1L), jet_at(j, s$lay, 2L)), grad_exact(x))

  # order two against one differentiation of the exact gradient
  H <- numDeriv::jacobian(grad_exact, x)
  expect_equal(jet_at(j, s$lay, c(1L, 1L)), H[1, 1], tolerance = 1e-7)
  expect_equal(jet_at(j, s$lay, c(1L, 2L)), H[1, 2], tolerance = 1e-7)
  expect_equal(jet_at(j, s$lay, c(2L, 2L)), H[2, 2], tolerance = 1e-7)
  expect_equal(H[1, 2], H[2, 1], tolerance = 1e-8)
})


test_that("a map written in ordinary R differentiates itself", {
  # The point of the operators: the expression below says nothing about
  # derivatives, and carries all thirty-five of them.
  s <- jet_seed(c(1.3, 0.8))
  x <- s$v[[1]]
  y <- s$v[[2]]
  X <- 1.3
  Y <- 0.8

  j <- x^2 * y^3
  expect_equal(j$v, X^2 * Y^3)
  expect_equal(jet_at(j, s$lay, 1L), 2 * X * Y^3)
  expect_equal(jet_at(j, s$lay, c(1L, 2L)), 6 * X * Y^2)
  expect_equal(jet_at(j, s$lay, c(1L, 1L, 2L, 2L)), 12 * Y)

  # a number on either side is a constant
  expect_equal((3 * x)$v, 3 * X)
  expect_equal((x * 3)$v, 3 * X)
  expect_equal((x + 3)$v, X + 3)
  expect_equal((3 - x)$v, 3 - X)
  expect_equal((x / 4)$v, X / 4)
  expect_equal((-x)$v, -X)
  expect_equal(jet_at(3 - x, s$lay, 1L), -1)
  expect_equal(jet_at(2 / x, s$lay, 1L), -2 / X^2)

  # and the operators agree with the functions they stand for
  expect_equal((x * y)$d, jet_mul(x, y, s$lay)$d)
  expect_equal((x / y)$d, jet_div(x, y, s$lay)$d)
  expect_equal((x^2.5)$d, jet_pow(x, 2.5, s$lay)$d)
  expect_equal(exp(x)$d, jet_exp(x, s$lay)$d)
  expect_equal(log(x)$d, jet_log(x, s$lay)$d)
  expect_equal(sqrt(x)$d, jet_sqrt(x, s$lay)$d)
  expect_equal(gamma(x)$d, jet_gamma(x, s$lay)$d)
  expect_equal(lgamma(x)$d, jet_lgamma(x, s$lay)$d)
  expect_equal(digamma(x)$d, jet_digamma(x, s$lay)$d)
})


test_that("the Weibull mean map differentiates itself", {
  # The map that reparametrising a Weibull to its mean needs, written as a
  # reader would write it. Its first derivative in sigma is checked against the
  # formula, which involves digamma through the gamma function.
  s <- jet_seed(c(4, 1.7))
  m <- s$v[[1]]
  sg <- s$v[[2]]
  M <- 4
  S <- 1.7

  scale <- m / gamma(1 + 1 / sg)
  g <- gamma(1 + 1 / S)
  expect_equal(scale$v, M / g)

  # d scale / d m = 1/g
  expect_equal(jet_at(scale, s$lay, 1L), 1 / g)
  # d scale / d sigma = M * digamma(1+1/S) / (S^2 g), by the chain rule
  expect_equal(jet_at(scale, s$lay, 2L), M * digamma(1 + 1 / S) / (S^2 * g))
  # and the pure second derivative in m vanishes, the map being linear in it
  expect_equal(jet_at(scale, s$lay, c(1L, 1L)), 0)
  # while the mixed one does not
  expect_equal(jet_at(scale, s$lay, c(1L, 2L)),
               digamma(1 + 1 / S) / (S^2 * g))
})


test_that("a jet refuses what it cannot carry", {
  s <- jet_seed(c(1.5, 2))
  x <- s$v[[1]]
  expect_error(x > 1, "not defined for a jet")
  expect_error(x == 1, "not defined for a jet")
  expect_error(abs(x), "not defined for a jet")
  expect_error(floor(x), "not defined for a jet")
  expect_error(log(x, base = 2), "no further arguments")
  expect_error(!x, "not defined for a jet")
})


test_that("a jet prints its value and its gradient", {
  s <- jet_seed(c(2, 3))
  out <- capture.output(print(s$v[[1]] * s$v[[2]]))
  expect_match(out[1], "<jet> value 6")
  expect_match(out[1], "2 variables")
  expect_match(out[2], "gradient")
})
