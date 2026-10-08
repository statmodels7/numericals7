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

test_that("the Bessel ratio and its derivatives match 150-digit values", {
  # columns: kappa, A, A', A'', A''', A''''; computed with mpmath at 150
  # digits from I1/I0 and the Riccati identity (stabilita/
  # bessel_ratio_reference.py). The rows sit in each regime and on both
  # sides of the switches at 0.5 and 30.
  ref <- rbind(
    c(1e-08, 5.0e-9, 4.9999999999999998e-1, -3.7499999999999999e-9, -3.7499999999999994e-1, 1.2499999999999999e-8),
    c(0.001, 4.9999993750001043e-4, 4.9999981250005208e-1, -3.7499979166674187e-4, -3.7499937500037598e-1, 1.2499984960946852e-3),
    c(0.3, 1.4833742694087526e-1, 4.8353791796564295e-1, -1.0705296836782769e-1, -3.2168516985297446e-1, 3.3657739069179652e-1),
    c(0.49, 2.37929523394175e-1, 4.5781908599094114e-1, -1.6122190505394156e-1, -2.4461853007867078e-1, 4.5935026960722952e-1),
    c(0.51, 2.4705333764769773e-1, 4.5454633924399987e-1, -1.6602189051206942e-1, -2.3535514502028998e-1, 4.6681155264285959e-1),
    c(5.0, 8.9338313704408522e-1, 2.31899430364522e-2, -1.0337671241085641e-2, 7.0240549085431692e-3, -6.2939350554317247e-3),
    c(19.0, 9.7331803639515719e-1, 1.424734954446128e-3, -1.5225466075414019e-4, 2.4423784403210615e-5, -5.2280061689460967e-6),
    c(29.9, 9.8313283326580569e-1, 5.6913539000685184e-4, -3.8417150398757939e-5, 3.8907845283066981e-6, -5.2553861312592275e-7),
    c(30.1, 9.8324589715371098e-1, 5.6152908076906912e-4, -3.7649387110457785e-5, 3.787428243384775e-6, -5.0814139990910202e-7),
    c(1000.0, 9.9949987487480428e-1, 5.0025037578328756e-7, -1.0007515039184817e-9, 3.0030075235231662e-12, -1.2015045164748185e-14),
    c(100000000.0, 9.9999999499999999e-1, 5.0000000250000004e-17, -1.0000000075000001e-24, 3.0000000300000007e-32, -1.2000000150000005e-39))
  k <- ref[, 1]
  got <- cbind(bessel_i_ratio(k), bessel_i_ratio_d1(k), bessel_i_ratio_d2(k),
               bessel_i_ratio_d3(k), bessel_i_ratio_d4(k))
  rel <- abs(got / ref[, -1] - 1)
  expect_lt(max(rel), 4e-15)
})

test_that("the derivatives match one numerical pass of the order below", {
  skip_if_not_installed("numDeriv")
  f <- list(bessel_i_ratio, bessel_i_ratio_d1, bessel_i_ratio_d2,
            bessel_i_ratio_d3, bessel_i_ratio_d4)
  for (k in c(0.2, 2, 20, 60)) {
    for (m in 1:4) {
      expect_equal(f[[m + 1]](k), numDeriv::grad(f[[m]], k), tolerance = 1e-8,
                   info = sprintf("kappa %g, order %d", k, m))
    }
  }
})

test_that("the Bessel ratio is finite at any argument and correct where naive is", {
  k <- c(0.01, 0.5, 2, 50, 600)
  a <- bessel_i_ratio(k)
  expect_true(all(a > 0 & a < 1))
  expect_true(all(diff(a) > 0))
  # against the unscaled ratio where that one does not overflow
  expect_equal(a, besselI(k, 1) / besselI(k, 0), tolerance = 1e-14)
  # far past the scaled Bessel underflow the ratio stays finite and ordered
  big <- bessel_i_ratio(c(1e5, 1e6, 1e9))
  expect_true(all(is.finite(big) & big > 0 & big < 1))
  expect_true(all(diff(big) > 0))
  # A' is the variance of a cosine, positive across the range
  expect_true(all(bessel_i_ratio_d1(c(0, 1e-6, 1, 100, 1e8)) > 0))
})

test_that("the ratio's limits and invalid arguments", {
  expect_identical(bessel_i_ratio(c(0, Inf)), c(0, 1))
  expect_identical(bessel_i_ratio_d1(c(0, Inf)), c(0.5, 0))
  expect_equal(c(bessel_i_ratio_d2(0), bessel_i_ratio_d3(0), bessel_i_ratio_d4(0)),
               c(0, -3 / 8, 0))
  for (f in list(bessel_i_ratio, bessel_i_ratio_d1, bessel_i_ratio_d4,
                 bessel_i_ratio_inverse_d2)) {
    out <- f(c(NA, NaN, -1))
    expect_true(is.na(out[1]) && !is.nan(out[1]))
    expect_true(all(is.nan(out[2:3])))
  }
  expect_identical(bessel_i_ratio(numeric(0)), numeric(0))
})

test_that("the kernels give the same bits at any thread count", {
  k <- 10^seq(-6, 6, length.out = 2000)
  for (f in list(bessel_i_ratio, bessel_i_ratio_d3, bessel_i_ratio_inverse_d4)) {
    expect_identical(f(k, threads = 2L), f(k))
  }
  r <- seq(0.0005, 0.9995, length.out = 2000)
  expect_identical(bessel_i_ratio_inverse(r, threads = 2L),
                   bessel_i_ratio_inverse(r))
})

test_that("the inverse matches 150-digit preimages and rejects the boundary", {
  # columns: rho (a double), the exact kappa with A(kappa) = rho (mpmath).
  # rho is written in hexadecimal: a decimal literal can parse to a different
  # double on arm64 macOS, and near rho = 1 one unit in the last place moves
  # kappa by 1.2e-4 relative (the last row is 1 - 2^-40).
  inv <- rbind(
    c(0x1.19799812dea11p-40, 2.0e-12),
    c(0x1.3333333333333p-2, 6.292153761056903e-1),
    c(0x1.f0a3d70a3d70ap-1, 1.6928871205888453e+1),
    c(0x1.ff7ced916872bp-1, 5.0025037594098552e+2),
    c(0x1.fffffffffep-1, 5.4975581388825e+11))
  expect_lt(max(abs(bessel_i_ratio_inverse(inv[, 1]) / inv[, 2] - 1)), 4e-15)
  expect_true(all(is.na(bessel_i_ratio_inverse(c(0, 1, -0.5, 2, NA)))))
})

test_that("the inverse is located to the rounding of rho, down to tiny rho", {
  rho <- c(1e-15, 1e-12, 1e-8, seq(0.02, 0.98, by = 0.04))
  k <- bessel_i_ratio_inverse(rho)
  expect_true(all(abs(bessel_i_ratio(k) - rho) <= 4 * .Machine$double.eps * rho))
  expect_equal(k[1:2], 2 * rho[1:2], tolerance = 1e-12)
  one <- vapply(rho, bessel_i_ratio_inverse, numeric(1))
  expect_identical(k, one)
})

test_that("the inverse's derivatives match one numerical pass each", {
  skip_if_not_installed("numDeriv")
  d <- list(function(r) bessel_i_ratio_inverse(r),
            function(r) bessel_i_ratio_inverse_d1(bessel_i_ratio_inverse(r)),
            function(r) bessel_i_ratio_inverse_d2(bessel_i_ratio_inverse(r)),
            function(r) bessel_i_ratio_inverse_d3(bessel_i_ratio_inverse(r)),
            function(r) bessel_i_ratio_inverse_d4(bessel_i_ratio_inverse(r)))
  for (r in c(0.2, 0.6, 0.9, 0.99)) {
    for (m in 1:4) {
      expect_equal(d[[m + 1]](r), numDeriv::grad(d[[m]], r), tolerance = 1e-7,
                   info = sprintf("rho %g, order %d", r, m))
    }
  }
  # at a small rho numDeriv is the weak side (4e-7 at the fourth order);
  # there the four derivatives against mpmath's at 60 digits
  k <- bessel_i_ratio_inverse(0.001)
  got <- c(bessel_i_ratio_inverse_d1(k), bessel_i_ratio_inverse_d2(k),
           bessel_i_ratio_inverse_d3(k), bessel_i_ratio_inverse_d4(k))
  expect_equal(got, c(2.0000030000041667, 0.0060000166666999168,
                      6.0000500001662504, 0.10000066500240241),
               tolerance = 1e-13)
  # the first is the reciprocal of A'
  k <- c(1e-5, 0.7, 25, 1e5)
  expect_equal(bessel_i_ratio_inverse_d1(k), 1 / bessel_i_ratio_d1(k),
               tolerance = 1e-15)
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
