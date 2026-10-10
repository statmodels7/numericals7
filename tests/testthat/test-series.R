# The vectorized series. References are closed forms and exact identities.

test_that("the geometric series matches its closed form for every ratio", {
  r <- c(0.1, 0.5, 0.9, 0.99)
  got <- series_vec(function(k, i) r[i]^k, n = 4)
  expect_equal(got, 1 / (1 - r), tolerance = 1e-9)
})


test_that("probability masses sum to one across very different rates", {
  lam <- c(0.05, 0.5, 4, 60, 300)
  got <- series_vec(function(k, i) dpois(k, lam[i]), n = 5)
  expect_equal(got, rep(1, 5), tolerance = 1e-10)

  # and the mean, a series whose early terms vanish
  m <- series_vec(function(k, i) k * dpois(k, lam[i]), n = 5)
  expect_equal(m, lam, tolerance = 1e-8)

  # negative binomial masses too, with the summation starting at zero
  mu <- c(2, 15)
  th <- c(0.7, 3)
  nb <- series_vec(function(k, i) dnbinom(k, size = th[i], mu = mu[i]), n = 2)
  expect_equal(nb, rep(1, 2), tolerance = 1e-10)
})


test_that("the tail guard sees past a block that sums to little", {
  # A hump-shaped term: for lambda = 300 the mass below k = 64 is essentially
  # zero, so the first block contributes nothing while the series has not
  # begun. Retiring on the block sum alone would return 0 for a sum that is 1.
  got <- series_vec(function(k, i) dpois(k, 300), n = 1, block = 64L)
  expect_equal(got, 1, tolerance = 1e-10)
})


test_that("a divergent row is refused and named, the others survive", {
  trm <- function(k, i) ifelse(i == 2L, 1 / (k + 1), 0.5^k)
  expect_warning(
    got <- series_vec(trm, n = 2, max_terms = 2000L),
    "rows 2"
  )
  expect_equal(got[1], 2, tolerance = 1e-9)
  expect_true(is.na(got[2]))
})


test_that("a nonzero starting index is honored", {
  # sum_{k=2}^inf r^k = r^2 / (1 - r)
  r <- c(0.3, 0.8)
  got <- series_vec(function(k, i) r[i]^k, n = 2, from = 2L)
  expect_equal(got, r^2 / (1 - r), tolerance = 1e-10)
})


test_that("convergence at a block boundary is not special", {
  # with block = 8 the geometric series retires on different passes per row
  r <- c(0.01, 0.6, 0.95)
  got <- series_vec(function(k, i) r[i]^k, n = 3, block = 8L)
  expect_equal(got, 1 / (1 - r), tolerance = 1e-9)
})


test_that("leading terms that underflow to zero do not end the sum", {
  # dpois(k, 1e4) is exactly zero below k ~ 8000 and k * dpois(k, 2000)
  # below k ~ 1400; the old rule retired both rows on their first block and
  # returned 0
  expect_equal(series_vec(function(k, i) dpois(k, 1e4), n = 1), 1,
               tolerance = 1e-10)
  expect_equal(series_vec(function(k, i) k * dpois(k, 2000), n = 1), 2000,
               tolerance = 1e-10)
  # the zeros past a finite support end it at once
  expect_equal(series_vec(function(k, i) dbinom(k, 10, 0.3), n = 1), 1,
               tolerance = 1e-14)
})


test_that("a polynomial tail is summed within the budget or returned NA", {
  # the tail estimate doubles the geometric continuation, which covers k^-p
  # for p >= 2; before, 1/k^3 retired at a relative error of 6e-9
  z3 <- 1.2020569031595942854
  expect_lt(abs(series_vec(function(k, i) 1 / k^3, n = 1, from = 1L) / z3 - 1),
            1e-10)
  expect_lt(abs(series_vec(function(k, i) 1 / k^4, n = 1, from = 1L) /
                  (pi^4 / 90) - 1), 1e-10)
  expect_lt(abs(series_vec(function(k, i) 0.999^k, n = 1) / 1000 - 1), 1e-10)
  # 1/k^2 does not reach its budget within the default max_terms
  expect_warning(v <- series_vec(function(k, i) 1 / k^2, n = 1, from = 1L),
                 "rows 1")
  expect_true(is.na(v))
})


test_that("a row of zeros returns 0 and block is validated", {
  # a structurally null series (a zero component of an expected Hessian)
  # is probed ahead and returns 0 without a warning
  expect_identical(series_vec(function(k, i) 0 * k, n = 1), 0)
  lam <- c(0.5, 1e4, 0)
  expect_equal(series_vec(function(k, i) dpois(k, lam[i]) * (lam[i] > 0), n = 3),
               c(1, 1, 0), tolerance = 1e-10)
  expect_equal(series_vec(function(k, i) dpois(k, 1e6), n = 1), 1,
               tolerance = 1e-10)
  expect_error(series_vec(function(k, i) 0.5^k, n = 1, block = 3L), "block")
  expect_error(series_vec(function(k, i) 0.5^k, n = 1, block = 6.5), "block")
})
