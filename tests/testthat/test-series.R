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
