# The log-Bessel functions. References: base R's scaled functions where they
# are healthy, half-integer closed forms, small- and large-argument
# asymptotics, the Wronskian identity across regimes, and the u_k recurrence.

test_that("log I and log K agree with base R where base R is healthy", {
  set.seed(1)
  x <- exp(runif(2000, log(0.05), log(600)))
  v <- exp(runif(2000, log(1e-2), log(80)))
  ref_i <- log(besselI(x, v, expon.scaled = TRUE)) + x
  ok <- is.finite(ref_i)
  expect_gt(sum(ok), 1900)
  expect_lt(max(abs(log_bessel_i(x[ok], v[ok]) - ref_i[ok]) /
                  pmax(1, abs(ref_i[ok]))), 1e-12)

  ref_k <- log(besselK(x, v, expon.scaled = TRUE)) - x
  ok <- is.finite(ref_k)
  expect_lt(max(abs(log_bessel_k(x[ok], v[ok]) - ref_k[ok]) /
                  pmax(1, abs(ref_k[ok]))), 1e-12)
})


test_that("the half-integer closed forms are reproduced", {
  x <- c(0.01, 0.1, 1, 10, 100, 1000, 1e5)
  ck12 <- 0.5 * log(pi / (2 * x)) - x
  expect_lt(max(abs(log_bessel_k(x, 0.5) - ck12) / pmax(1, abs(ck12))), 1e-14)
  ck32 <- ck12 + log1p(1 / x)
  expect_lt(max(abs(log_bessel_k(x, 1.5) - ck32) / pmax(1, abs(ck32))), 1e-14)
  # I_{1/2}(x) = sqrt(2/(pi x)) sinh(x)
  ci12 <- 0.5 * log(2 / (pi * x)) + x + log(-expm1(-2 * x)) - log(2)
  expect_lt(max(abs(log_bessel_i(x, 0.5) - ci12) / pmax(1, abs(ci12))), 1e-13)
})


test_that("the regimes where base R fails are finite and correct", {
  # small argument, large order: base R's scaled besselI loses its precision
  # and returns 0; the series limit v log(x/2) - lgamma(v+1) is exact there
  for (v in c(150, 500, 5000)) {
    ref <- v * log(5e-4) - lgamma(v + 1) + log1p(1e-6 / (4 * (v + 1)))
    expect_equal(log_bessel_i(1e-3, v), ref, tolerance = 1e-13)
  }
  # large argument, past the scaled underflow of base R (about 1e5-1e6)
  x <- 1e7
  ref <- x - 0.5 * log(2 * pi * x) - (4 * 0.3^2 - 1) / (8 * x)
  expect_equal(log_bessel_i(x, 0.3), ref, tolerance = 1e-13)
  # K at an order where the unscaled function overflows
  expect_true(is.finite(log_bessel_k(1, 500)))
  ref <- 500 * log(2 / 1) + lgamma(500) - log(2)   # K_v(x) ~ Gamma(v)/2 (2/x)^v
  expect_equal(log_bessel_k(1, 500), ref, tolerance = 1e-3)
})


test_that("the Wronskian identity holds across every regime", {
  # I_v K_{v+1} + I_{v+1} K_v = 1/x. At v >> x the two log terms reach 1e5
  # in magnitude and cancel down to -log x, so ~4e-11 of absolute error on
  # the sum is machine precision of the individual logs, not a defect.
  set.seed(2)
  x <- exp(runif(5000, log(1e-4), log(1e6)))
  v <- exp(runif(5000, log(1e-3), log(1e4)))
  a <- log_bessel_i(x, v) + log_bessel_k(x, v + 1)
  b <- log_bessel_i(x, v + 1) + log_bessel_k(x, v)
  m <- pmax(a, b)
  lhs <- m + log(exp(a - m) + exp(b - m))
  expect_lt(max(abs(lhs + log(x)) / pmax(1, abs(log(x)))), 1e-10)
})


test_that("the Rothwell integral agrees with base R in the moderate region", {
  # the integral serves the overflow corner, so its own accuracy is checked
  # against the healthy route it replaces elsewhere
  x <- c(0.5, 2, 8, 15)
  v <- c(0.3, 1.7, 6, 11)
  ref <- log(besselK(x, v, expon.scaled = TRUE)) - x
  expect_lt(max(abs(.lbk_integral(x, v) - ref) / pmax(1, abs(ref))), 1e-9)
})


test_that("the u_k table satisfies its recurrence, coefficient by coefficient", {
  # u_{k+1} = t^2(1-t^2)/2 u_k' + (1/8) int_0^t (1-5s^2) u_k ds. Evaluating
  # the high-k polynomials invites catastrophic cancellation (u_13 carries
  # coefficients of 1e12 that cancel to order one), so the recurrence is
  # checked on the COEFFICIENTS, where every operation is a division by a
  # small integer: a corrupted entry cannot pass, and no tolerance argument
  # is needed beyond machine precision.
  step_coefs <- function(e, cc) {
    out <- new.env()
    put <- function(ee, val) {
      key <- as.character(ee)
      assign(key, val + mget(key, out, ifnotfound = 0)[[1]], out)
    }
    for (j in seq_along(e)) {
      # t^2(1-t^2)/2 d/dt: c e/2 t^{e+1} - c e/2 t^{e+3}
      put(e[j] + 1, cc[j] * e[j] / 2)
      put(e[j] + 3, -cc[j] * e[j] / 2)
      # (1/8) int (1-5s^2): c/(8(e+1)) t^{e+1} - 5c/(8(e+3)) t^{e+3}
      put(e[j] + 1, cc[j] / (8 * (e[j] + 1)))
      put(e[j] + 3, -5 * cc[j] / (8 * (e[j] + 3)))
    }
    ee <- sort(as.integer(ls(out)))
    list(e = ee, c = vapply(as.character(ee), get, numeric(1), envir = out))
  }
  for (k in 1:12) {
    got <- step_coefs(.bessel_uk[[k]]$e, .bessel_uk[[k]]$c)
    keep <- abs(got$c) > 1e-8 * max(abs(got$c))
    expect_equal(got$e[keep], .bessel_uk[[k + 1]]$e,
                 label = paste("exponents of u", k + 1))
    expect_equal(unname(got$c[keep]), unname(.bessel_uk[[k + 1]]$c),
                 tolerance = 1e-13,
                 label = paste("coefficients of u", k + 1))
  }
  # and u_1 itself against DLMF: (3t - 5t^3)/24
  expect_equal(.bessel_uk[[1]]$c, c(3, -5) / 24, tolerance = 1e-15)
})


test_that("derivatives match closed forms and one numerical pass", {
  x <- c(0.5, 2, 40)
  dk <- log_bessel_k_derivs(x, 0.5)
  expect_equal(dk$d1, -0.5 / x - 1, tolerance = 1e-14)
  expect_equal(dk$d2, 0.5 / x^2, tolerance = 1e-13)
  expect_equal(dk$d3, -1 / x^3, tolerance = 1e-12)
  expect_equal(dk$d4, 3 / x^4, tolerance = 1e-11)

  # d/dx log I_{1/2} = coth(x) - 1/(2x)
  di <- log_bessel_i_derivs(x, 0.5)
  expect_equal(di$d1, 1 / tanh(x) - 0.5 / x, tolerance = 1e-13)

  skip_if_not_installed("numDeriv")
  for (p in list(c(3, 2), c(50, 0.4), c(8, 30), c(2000, 100), c(5, 800))) {
    di <- log_bessel_i_derivs(p[1], p[2])
    nd <- numDeriv::grad(function(z) log_bessel_i(z, p[2]), p[1])
    expect_equal(di$d1, nd, tolerance = 1e-6, label = sprintf("dI at %g,%g", p[1], p[2]))
    d2 <- numDeriv::grad(function(z) log_bessel_i_derivs(z, p[2])$d1, p[1])
    expect_equal(di$d2, d2, tolerance = 1e-5, label = sprintf("d2I at %g,%g", p[1], p[2]))
    dk <- log_bessel_k_derivs(p[1], p[2])
    nd <- numDeriv::grad(function(z) log_bessel_k(z, p[2]), p[1])
    expect_equal(dk$d1, nd, tolerance = 1e-6, label = sprintf("dK at %g,%g", p[1], p[2]))
  }
})


test_that("edges and recycling behave", {
  expect_identical(log_bessel_i(0, 0), 0)
  expect_identical(log_bessel_i(0, 2), -Inf)
  expect_identical(log_bessel_k(0, 1), Inf)
  expect_true(is.na(log_bessel_i(-1, 1)))
  expect_true(is.na(log_bessel_i(1, -1)))
  expect_true(is.na(log_bessel_i(NA, 1)))
  # K is even in its order
  expect_equal(log_bessel_k(3, -2.5), log_bessel_k(3, 2.5))
  # recycling
  expect_length(log_bessel_i(1:5, 2), 5L)
  expect_length(log_bessel_k(2, 1:7), 7L)
})
