# The special functions. References are closed identities, asymptotic forms,
# and one Richardson pass on an analytic quantity -- never a nested difference.

test_that("the Mills ratio survives the deep tail and matches its derivative", {
  m <- mills_ratio(c(-60, -40, -10, 0, 10))
  expect_true(all(is.finite(m$r)))
  expect_true(all(is.finite(m$dr)))

  # in the deep tail R(t) ~ -t - 1/t + O(t^-3)
  expect_equal(mills_ratio(-60)$r, 60 + 1 / 60, tolerance = 1e-4)

  # at zero, R(0) = phi(0)/Phi(0) = 2 phi(0) = sqrt(2/pi)
  expect_equal(mills_ratio(0)$r, sqrt(2 / pi), tolerance = 1e-14)

  # dr against a central difference of r (one stencil on an analytic value)
  h <- 1e-6
  num <- (mills_ratio(-2 + h)$r - mills_ratio(-2 - h)$r) / (2 * h)
  expect_equal(mills_ratio(-2)$dr, num, tolerance = 1e-7)
})


test_that("Owen's T matches its closed identities", {
  # a = 0, antisymmetry, and the infinite-a identity
  expect_equal(owen_t(c(0.5, -1.2), 0), c(0, 0))
  expect_equal(owen_t(1.3, Inf), 0.5 * pnorm(-1.3))
  expect_equal(owen_t(1.3, -2), -owen_t(1.3, 2))

  # T(0, a) = arctan(a) / (2 pi)
  a <- c(0.3, 1, 5)
  expect_equal(owen_t(0, a), atan(a) / (2 * pi), tolerance = 1e-12)

  # T(h, 1) = Phi(h) Phi(-h) / 2
  hs <- c(-2, -0.7, 0, 1.1, 3)
  expect_equal(owen_t(hs, 1), pnorm(hs) * pnorm(-hs) / 2, tolerance = 1e-12)

  # a non-finite h contributes zero
  expect_equal(owen_t(c(Inf, -Inf), 1), c(0, 0))

  # the batched evaluation agrees with one scalar quadrature per element
  hh <- c(-1.5, 0.2, 0.9, 2.4)
  aa <- c(0.4, 1.7, 3, 0.8)
  ref <- vapply(seq_along(hh), function(j) {
    stats::integrate(function(x) exp(-hh[j]^2 * (1 + x^2) / 2) / (1 + x^2),
                     0, aa[j], rel.tol = 1e-12)$value / (2 * pi)
  }, numeric(1))
  expect_equal(owen_t(hh, aa), ref, tolerance = 1e-11)
})


test_that("the Bessel ratio is finite at any argument and correct where naive is", {
  k <- c(0.01, 0.5, 2, 50, 600)
  a <- bessel_i_ratio(k)
  expect_true(all(a > 0 & a < 1))
  expect_true(all(diff(a) > 0))

  # against the unscaled ratio where that one does not overflow
  naive <- besselI(k, 1) / besselI(k, 0)
  expect_equal(a, naive, tolerance = 1e-12)

  # the two branches agree where both are exact: below the switch the
  # implementation is the Bessel ratio, and the test-side asymptotic series
  # must match it to machine precision
  for (kk in c(5e3, 9.99e3)) {
    expect_equal(bessel_i_ratio(kk),
                 1 - 1 / (2 * kk) - 1 / (8 * kk^2) - 1 / (8 * kk^3),
                 tolerance = 1e-15)
  }

  # far past the scaled Bessel underflow the ratio stays finite and ordered
  big <- bessel_i_ratio(c(1e5, 1e6, 1e9))
  expect_true(all(is.finite(big) & big > 0 & big < 1))
  expect_true(all(diff(big) > 0))
})


test_that("the ratio's derivatives match one numerical pass each", {
  skip_if_not_installed("numDeriv")
  for (k in c(0.3, 2, 20)) {
    a <- bessel_i_ratio_derivs(k)
    expect_equal(a$d1, numDeriv::grad(bessel_i_ratio, k), tolerance = 1e-8)
    expect_equal(a$d2, numDeriv::grad(function(z) bessel_i_ratio_derivs(z)$d1, k),
                 tolerance = 1e-7)
    expect_equal(a$d3, numDeriv::grad(function(z) bessel_i_ratio_derivs(z)$d2, k),
                 tolerance = 1e-6)
    expect_equal(a$d4, numDeriv::grad(function(z) bessel_i_ratio_derivs(z)$d3, k),
                 tolerance = 1e-5)
  }
  # d1 is the variance of a cosine, so positive across the range
  expect_true(all(bessel_i_ratio_derivs(c(0.01, 1, 100))$d1 > 0))
})


test_that("the inverse round-trips and refuses the boundary", {
  rho <- c(1e-6, 0.1, 0.53, 0.7, 0.85, 0.99, 0.999999)
  k <- bessel_i_ratio_inverse(rho)$kappa
  expect_equal(bessel_i_ratio(k), rho, tolerance = 1e-10)

  expect_true(is.na(bessel_i_ratio_inverse(0)$kappa))
  expect_true(is.na(bessel_i_ratio_inverse(1)$kappa))

  # the inverse function rule against the analytic forward derivative
  inv <- bessel_i_ratio_inverse(0.7)
  expect_equal(inv$d1, 1 / bessel_i_ratio_derivs(inv$kappa)$d1, tolerance = 1e-12)
})


test_that("the inverse's higher derivatives match one numerical pass each", {
  skip_if_not_installed("numDeriv")
  for (r in c(0.2, 0.6, 0.9)) {
    kd <- bessel_i_ratio_inverse(r)
    expect_equal(kd$d1,
                 numDeriv::grad(function(z) bessel_i_ratio_inverse(z)$kappa, r),
                 tolerance = 1e-7)
    expect_equal(kd$d2,
                 numDeriv::grad(function(z) bessel_i_ratio_inverse(z)$d1, r),
                 tolerance = 1e-6)
    expect_equal(kd$d3,
                 numDeriv::grad(function(z) bessel_i_ratio_inverse(z)$d2, r),
                 tolerance = 1e-5)
    expect_equal(kd$d4,
                 numDeriv::grad(function(z) bessel_i_ratio_inverse(z)$d3, r),
                 tolerance = 1e-4)
  }
})

test_that("the Bessel ratios agree with an independent evaluation", {
  # R's own besselI, which shares no arithmetic with the backward recurrence
  for (kap in c(0.01, 0.1, 1, 5, 20, 100, 500)) {
    for (m in c(5L, 40L)) {
      got <- bessel_i_ratios(kap, m)[1L, ]
      ref <- besselI(kap, seq_len(m), expon.scaled = TRUE) /
        besselI(kap, 0, expon.scaled = TRUE)
      ok <- is.finite(ref) & ref > 0
      expect_true(any(ok))
      expect_equal(got[ok], ref[ok], tolerance = 1e-12,
                   info = sprintf("kappa %g, m %d", kap, m))
    }
  }
  # the first of them is bessel_i_ratio()
  k <- c(0.3, 2, 40)
  expect_equal(bessel_i_ratios(k, 1L)[, 1L], bessel_i_ratio(k),
               tolerance = 1e-12)
})

test_that("the recurrence's loop runs over the order, not over the data", {
  # a vector of arguments gives exactly what one at a time does, which is
  # what says the ratios are not being recomputed per element
  k <- c(0.5, 3, 17)
  a <- bessel_i_ratios(k, 30L)
  b <- do.call(rbind, lapply(k, function(z) bessel_i_ratios(z, 30L)[1L, ]))
  expect_identical(a, b)
  expect_identical(dim(a), c(3L, 30L))
})

test_that("bessel_i_ratios rejects what it cannot answer", {
  expect_error(bessel_i_ratios(1, 0), "positive integer")
  expect_error(bessel_i_ratios(-1, 3), "must be positive")
  expect_identical(dim(bessel_i_ratios(numeric(0), 4L)), c(0L, 4L))
  expect_true(all(is.na(bessel_i_ratios(NA_real_, 3L))))
})
