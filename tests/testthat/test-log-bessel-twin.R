# The compiled kernels against their R twins, across every branch of the
# input-range table. The two share the u_k coefficient table and nothing
# else -- different loop shapes, different accumulation orders -- so
# agreement here checks the transcription both ways.

test_that("the compiled log-Bessel kernels match the R twins on every branch", {
  set.seed(11)
  n <- 20000
  x <- exp(runif(n, log(1e-4), log(1e6)))
  v <- exp(runif(n, log(1e-3), log(1e4)))
  i_c <- log_bessel_i(x, v)
  i_r <- numericals7:::.log_bessel_i_r(x, v)
  expect_lt(max(abs(i_c - i_r) / pmax(1, abs(i_r))), 1e-13)
  k_c <- log_bessel_k(x, v)
  k_r <- numericals7:::.log_bessel_k_r(x, v)
  expect_lt(max(abs(k_c - k_r) / pmax(1, abs(k_r))), 1e-13)

  # the edge cases the loops treat before any branch
  expect_identical(log_bessel_i(0, c(0, 2)), c(0, -Inf))
  expect_identical(log_bessel_k(0, 1), Inf)
  expect_true(is.na(log_bessel_i(NA, 1)))
  expect_true(is.na(log_bessel_i(-1, 1)))
  expect_true(is.na(log_bessel_k(-1, 1)))
  # K is even in the order
  expect_equal(log_bessel_k(3, -2.5), log_bessel_k(3, 2.5))
})
