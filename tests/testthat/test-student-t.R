# student_t_cdf() is R's pt() compiled without its warnings, so the check is
# identity with pt() over a grid of degrees of freedom and quantiles that
# reaches both tails, plus a comparison with values computed at 60 digits
# (mpmath's regularized incomplete beta), which pt() itself meets to 1e-13.

test_that("the values are pt()'s, bit for bit", {
  nu <- c(0.1, 0.3, 0.5, 1, 1.7, 2.5, 5, 10, 30, 100, 1e3, 1e4, 1e5, 4e5, 1e6)
  mag <- c(0, 10^seq(-8, 3, by = 0.25))
  x <- sort(unique(c(-mag, mag, -1e6, 1e6, -1e60, 1e60)))
  g <- expand.grid(x = x, nu = nu)
  for (lt in c(TRUE, FALSE)) for (lg in c(TRUE, FALSE)) {
    expect_identical(student_t_cdf(g$x, g$nu, lower.tail = lt, log.p = lg),
                     pt(g$x, g$nu, lower.tail = lt, log.p = lg))
  }
})

test_that("the log-probability meets 60-digit values to 1e-13", {
  # log P(T <= x) at 60 digits
  ref <- data.frame(
    x = c(-31.6227766, -1.77827941, 1.77827941, -31.6227766, -3, 5,
          -31.6227766, -1.77827941, 0.25),
    nu = c(2.5, 2.5, 2.5, 1e3, 1e5, 1e5, 1e6, 1e6, 1e6),
    lp = c(-8.96654337279219152660907863924, -2.34619648309069266666066246029,
           -0.100630156084857319511187612098, -350.601078367702585197118262721,
           -6.60748000329726318009686329422, -2.87135125155369746107976039145e-7,
           -504.123480671311460173780271706, -3.27864865249355387686840093104,
           -0.512984118297196467425126851237))
  got <- student_t_cdf(ref$x, ref$nu, log.p = TRUE)
  expect_lt(max(abs(got - ref$lp) / abs(ref$lp)), 1e-13)
})

test_that("the edges and the thread count", {
  expect_identical(student_t_cdf(c(-Inf, Inf), 3), c(0, 1))
  expect_identical(student_t_cdf(c(-Inf, Inf), 3, log.p = TRUE), c(-Inf, 0))
  expect_identical(student_t_cdf(1.3, Inf), pnorm(1.3))
  expect_true(is.nan(student_t_cdf(1, 0)))
  expect_true(is.nan(student_t_cdf(1, -2)))
  set.seed(1)
  q <- rt(5000, 3) * 2
  nu <- exp(runif(5000, -1, 8))
  expect_identical(student_t_cdf(q, nu, log.p = TRUE, threads = 4L),
                   student_t_cdf(q, nu, log.p = TRUE))
})
