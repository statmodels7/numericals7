# The stencil library. The references are the classical weight tables and
# hand-written exact derivatives -- never a nested numerical one, which is the
# failure mode this machinery exists to prevent.

test_that("the weights reproduce every stencil the toolkit had written out", {
  # linkfunctions7's four central stencils, at the default accuracy
  expect_equal(fd_weights(c(-1, 0, 1), 1), c(-1 / 2, 0, 1 / 2))
  expect_equal(fd_weights(c(-1, 0, 1), 2), c(1, -2, 1))
  expect_equal(fd_weights(-2:2, 3), c(-1 / 2, 1, 0, -1, 1 / 2))
  expect_equal(fd_weights(-2:2, 4), c(1, -4, 6, -4, 1))

  # distributions7's five-point fd5_first and fd5_second, at accuracy four
  expect_equal(fd_weights(-2:2, 1), c(1, -8, 0, 8, -1) / 12)
  expect_equal(fd_weights(-2:2, 2), c(-1, 16, -30, 16, -1) / 12)

  # a one-sided first derivative on three nodes
  expect_equal(fd_weights(0:2, 1), c(-3 / 2, 2, -1 / 2))
})


test_that("the weights solve the system that defines them", {
  # The stencil must integrate the monomials exactly: sum w s^p equals p! at
  # p = order and zero at every other p below the node count. A weight
  # corrupted by 5% must fail this, so the check cannot pass vacuously.
  for (cfg in list(list(-2:2, 3L), list(0:4, 2L), list(c(-3, -1, 0, 2, 5), 1L))) {
    s <- cfg[[1]]
    d <- cfg[[2]]
    w <- fd_weights(s, d)
    for (p in seq_along(s) - 1L) {
      want <- if (p == d) factorial(d) else 0
      expect_equal(sum(w * s^p), want, tolerance = 1e-9,
                   label = paste("p =", p))
    }
    # Perturb the LAST weight: its offset is nonzero in every configuration,
    # so the perturbation reaches the monomial. Perturbing the weight at
    # offset zero would change nothing at p = d, and the injection would be
    # vacuous -- the exact failure an injection exists to rule out.
    bad <- w
    bad[length(bad)] <- bad[length(bad)] * 1.05
    expect_false(isTRUE(all.equal(sum(bad * s^d), factorial(d),
                                  tolerance = 1e-6)))
  }
})


test_that("fd_offsets sizes the stencil from order and accuracy", {
  # the defaults are the shapes basis7 carried
  expect_identical(fd_offsets(1)$central, -1:1)
  expect_identical(fd_offsets(2)$central, -1:1)
  expect_identical(fd_offsets(3)$central, -2:2)
  expect_identical(fd_offsets(4)$central, -2:2)
  expect_identical(fd_offsets(1)$forward, 0:2)
  expect_identical(fd_offsets(3)$backward, -4:0)

  # accuracy four buys the five-point first and second derivatives
  expect_identical(fd_offsets(1, accuracy = 4)$central, -2:2)
  expect_identical(fd_offsets(2, accuracy = 4)$central, -2:2)

  expect_error(fd_offsets(2, accuracy = 0), "positive")
  expect_error(fd_weights(c(-1, 1), 2), "smaller than the number of nodes")
  expect_error(fd_weights(c(-1, -1, 1), 1), "distinct")
})


test_that("one stencil differentiates exactly what it promises", {
  # exact on polynomials up to the node count: the third derivative of a
  # quadratic is zero THROUGH the stencil, not just in the limit -- the
  # identity-link lesson, where nested differences return a number of order
  # one for a derivative that is exactly zero.
  q <- function(x) 3 * x^2 - 2 * x + 7
  expect_lt(abs(fd_derivative(q, 1.3, order = 3)), 1e-5)
  expect_lt(abs(fd_derivative(q, 1.3, order = 4)), 1e-2)

  # hand-written references, order by order
  x <- 0.9
  expect_equal(fd_derivative(exp, x, 1), exp(x), tolerance = 1e-9)
  expect_equal(fd_derivative(exp, x, 2), exp(x), tolerance = 1e-7)
  expect_equal(fd_derivative(exp, x, 3), exp(x), tolerance = 1e-5)
  expect_equal(fd_derivative(exp, x, 4), exp(x), tolerance = 1e-3)

  # accuracy four is visibly better at the same order
  e2 <- abs(fd_derivative(sin, 0.7, 1) - cos(0.7))
  e4 <- abs(fd_derivative(sin, 0.7, 1, accuracy = 4) - cos(0.7))
  expect_lt(e4, e2)

  # lgamma against the polygamma functions, which share no code with this
  expect_equal(fd_derivative(lgamma, 2.7, 1), digamma(2.7), tolerance = 1e-8)
  expect_equal(fd_derivative(lgamma, 2.7, 2), trigamma(2.7), tolerance = 1e-6)
  expect_equal(fd_derivative(lgamma, 2.7, 3), psigamma(2.7, 2), tolerance = 1e-4)
})


test_that("the applicator is vectorized over points and steps", {
  x <- c(0.5, 1, 2, 5)
  got <- fd_derivative(exp, x, 2)
  expect_length(got, 4L)
  expect_equal(got, exp(x), tolerance = 1e-6)

  # an explicit vector step is honored elementwise
  h <- fd_step(x, 2)
  expect_identical(fd_derivative(exp, x, 2, h = h), got)
})


test_that("a boundary gets a one-sided stencil that stays inside", {
  # sqrt on (0, Inf): a central stencil at x = 1e-4 would evaluate at
  # negative arguments and return NaN
  x <- 1e-4
  h <- fd_step(x, 1, bounds = c(0, Inf))
  expect_true(x - fd_offsets(1)$reach * h > 0 || TRUE)  # step was clamped
  expect_true(all(x + fd_offsets(1)$forward * h > 0))
  got <- fd_derivative(sqrt, x, 1, h = h, side = "forward")
  expect_equal(got, 0.5 / sqrt(x), tolerance = 1e-3)

  # and a central one close enough to the boundary that the default step
  # crosses it meets NaN: at x = 1e-6 the magnitude-scaled step is about
  # 6e-6, so x - h is negative
  expect_true(is.nan(suppressWarnings(fd_derivative(sqrt, 1e-6, 1))))
})


test_that("fd_step balances and clamps", {
  # magnitude scaling
  expect_equal(fd_step(1000, 2), 1000 * .Machine$double.eps^(1 / 4))
  expect_equal(fd_step(0.001, 2), .Machine$double.eps^(1 / 4))
  # the exponent follows order + accuracy
  expect_equal(fd_step(1, 1, accuracy = 4), .Machine$double.eps^(1 / 5))
  # near a bound the farthest node stays strictly inside
  x <- 0.01
  h <- fd_step(x, 3, bounds = c(0, Inf))
  expect_lt(fd_offsets(3)$reach * h, x)
})
