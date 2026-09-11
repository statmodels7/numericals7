# The smoothers of the absolute value: the three shipped instances pass the
# validator, the closed derivatives agree with one numerical differentiation
# of the order below, the declared identities hold, and an injected defect is
# caught -- the check_link bargain, restated for smoothers.

test_that("the three shipped smoothers pass check_abs_smoother", {
  skip_if_not_installed("numDeriv")
  for (sm in list(smooth_probit(h = 0.3), smooth_hyperbolic(c = 0.09),
                  smooth_quintic(h = 0.3))) {
    res <- check_abs_smoother(sm, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = sm@smoother_name)
  }
})

test_that("the probit smoother is E|u + hZ| and declares the convolution
          correction", {
  sm <- smooth_probit(h = 0.4)
  # against the expectation computed by quadrature, which shares no code
  # with the closed form
  for (u in c(-1.1, -0.2, 0, 0.37, 2.5)) {
    ref <- stats::integrate(function(z) abs(u + 0.4 * z) * stats::dnorm(z),
                            -Inf, Inf, rel.tol = 1e-10)$value
    expect_equal(smoother_deriv(sm, u, order = 0L), ref, tolerance = 1e-8)
  }
  # tau^2_apparent = tau^2_true + h^2, exactly
  expect_equal(sm@tau_correction(sqrt(0.5^2 + 0.4^2), 0.4), 0.5,
               tolerance = 1e-12)
  # and it does not go imaginary below the floor
  expect_identical(sm@tau_correction(0.1, 0.4), 0)
})

test_that("the quintic is exact outside the transition and C^3 at the seam", {
  h <- 0.5
  sm <- smooth_quintic(h = h)
  u <- c(-3, -0.51, 0.51, 1.7, 12)
  expect_identical(smoother_deriv(sm, u, order = 0L), abs(u))
  expect_identical(smoother_deriv(sm, u, order = 1L), sign(u))
  expect_identical(smoother_deriv(sm, u, order = 2L), rep(0, 5L))
  # orders 0..3 continuous at the seam, order 4 jumps: the documented C^3
  eps <- 1e-9
  for (k in 0:3) {
    expect_equal(smoother_deriv(sm, h - eps, order = k),
                 smoother_deriv(sm, h + eps, order = k), tolerance = 1e-6)
  }
  expect_gt(abs(smoother_deriv(sm, h - eps, order = 4L) -
                  smoother_deriv(sm, h + eps, order = 4L)), 1)
  expect_identical(sm@exact_radius(h), h)
})

test_that("the hyperbolic smoother is sqrt(u^2 + c) with no correction", {
  sm <- smooth_hyperbolic(c = 0.04)
  u <- c(-2, -0.1, 0, 0.3, 5)
  expect_equal(smoother_deriv(sm, u, order = 0L), sqrt(u^2 + 0.04))
  expect_null(sm@tau_correction)
  # its parameter is a squared length: the width resolved from a spacing is
  # the spacing's square (asked of a smoother holding no width of its own;
  # one that does keeps it, which the next test asserts)
  expect_equal(smoother_width(smooth_hyperbolic(), 0.3), 0.09)
})

test_that("a width the smoother holds wins over the spacing", {
  expect_equal(smoother_width(smooth_probit(h = 0.5), 0.3), 0.5)
  expect_equal(smoother_width(smooth_probit(), 0.3), 0.3)
  expect_equal(smoother_width(smooth_probit(), c(0.3, 0.7)), c(0.3, 0.7))
  expect_error(smoother_width(smooth_probit(), -1), "positive")
})

test_that("the width floor is sqrt(eps) times the scale, on the width's own
          scale", {
  D <- 10
  expect_equal(smoother_width_floor(smooth_probit(), D),
               sqrt(.Machine$double.eps) * D)
  expect_equal(smoother_width_floor(smooth_hyperbolic(), D),
               .Machine$double.eps * D^2)
})

test_that("a per-observation width vectorizes through every order", {
  sm <- smooth_probit()
  u <- c(-0.5, 0, 0.5, 1)
  w <- c(0.1, 0.2, 0.3, 0.4)
  for (k in 0:5) {
    one <- vapply(seq_along(u),
                  function(i) smoother_deriv(sm, u[i], w[i], k), numeric(1))
    expect_equal(smoother_deriv(sm, u, w, k), one)
  }
})

test_that("the validator catches a wrong derivative and a missing order", {
  skip_if_not_installed("numDeriv")
  sm <- smooth_probit(h = 0.3)
  # a second derivative 5% out fails its own order and the one above,
  # whose reference it is
  broken <- sm
  S7::prop(broken, "s")[[3L]] <-
    function(u, width) 2.1 * stats::dnorm(u / width) / width
  res <- check_abs_smoother(broken, verbose = FALSE)
  expect_true(any(res$status == "FAILED"))
  # an unbuildable list of derivatives is refused at construction
  expect_error(
    abs_smoother(smoother_name = "short", width = 0.3,
                 s = list(function(u, width) sqrt(u^2 + width))),
    "six functions")
  # and the unbroken smoother still passes, so the check is not trivially red
  expect_true(all(check_abs_smoother(sm, verbose = FALSE)$status == "OK"))
})

test_that("an unresolved width is reported rather than guessed", {
  expect_error(smoother_deriv(smooth_probit(), 1), "unresolved")
})

test_that("a user-written smoother passes through the same validator", {
  skip_if_not_installed("numDeriv")
  # the logistic smoother s(u) = 2h log(2 cosh(u / (2h)))... written more
  # simply as h * (z + 2 * log1p(exp(-z))) with z = u/h scaled so that
  # s'' = 1/(h (1 + cosh(z))): even, convex, |s'| < 1
  s0 <- function(u, width) {
    z <- abs(u) / width
    width * (z + 2 * log1p(exp(-z)))
  }
  s1 <- function(u, width) {
    z <- u / width
    tanh(z / 2)
  }
  s2 <- function(u, width) {
    z <- u / width
    1 / (2 * width * cosh(z / 2)^2)
  }
  s3 <- function(u, width) {
    z <- u / width
    -sinh(z / 2) / (2 * width^2 * cosh(z / 2)^3)
  }
  s4 <- function(u, width) {
    z <- u / width
    (cosh(z) - 2) / (4 * width^3 * cosh(z / 2)^4)
  }
  s5 <- function(u, width) {
    z <- u / width
    sinh(z / 2) * (10 - 2 * cosh(z)) / (8 * width^4 * cosh(z / 2)^5)
  }
  sm <- abs_smoother(smoother_name = "logistic", width = 0.3,
                     s = list(s0, s1, s2, s3, s4, s5))
  res <- check_abs_smoother(sm, verbose = FALSE)
  expect_true(all(res$status == "OK"))
})
