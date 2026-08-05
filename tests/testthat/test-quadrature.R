# The vectorized quadrature. The references are closed forms and
# stats::integrate, which is an independent implementation (QUADPACK through
# a different path, one theta at a time).

test_that("the Gauss-Kronrod pair is pinned by its defining property", {
  r <- gauss_kronrod15()
  expect_length(r$nodes, 15L)
  # a rule on [-1, 1] has weights summing to its length
  expect_equal(sum(r$wk), 2, tolerance = 1e-14)
  expect_equal(sum(r$wg), 2, tolerance = 1e-14)
  # symmetry
  expect_equal(r$nodes, -rev(r$nodes))
  expect_equal(r$wk, rev(r$wk))
  # the embedded Gauss rule uses seven of the fifteen nodes
  expect_identical(sum(r$wg > 0), 7L)

  # exactness: G7 to degree 13, K15 to degree 22. The odd degrees vanish by
  # symmetry, so the even ones carry the check.
  mono <- function(w, p) sum(w * r$nodes^p)
  for (p in seq(0, 12, by = 2)) {
    expect_equal(mono(r$wg, p), 2 / (p + 1), tolerance = 1e-13,
                 label = paste("G7 degree", p))
  }
  for (p in seq(0, 22, by = 2)) {
    expect_equal(mono(r$wk, p), 2 / (p + 1), tolerance = 1e-13,
                 label = paste("K15 degree", p))
  }
  # and one degree beyond, each rule must fail: exactness cannot be vacuous
  expect_false(isTRUE(all.equal(mono(r$wg, 14), 2 / 15, tolerance = 1e-10)))
  expect_false(isTRUE(all.equal(mono(r$wk, 24), 2 / 25, tolerance = 1e-10)))
})


test_that("finite intervals: many beta densities in one call", {
  a <- c(2, 0.8, 5, 1, 3)
  b <- c(5, 1.2, 5, 1, 0.7)
  f <- function(x, i) dbeta(x, a[i], b[i])
  got <- quad_vec(f, 0, rep(1, 5))
  expect_equal(got, rep(1, 5), tolerance = 1e-9)

  # the means, against the closed form a/(a+b)
  g <- function(x, i) x * dbeta(x, a[i], b[i])
  expect_equal(quad_vec(g, 0, rep(1, 5)), a / (a + b), tolerance = 1e-8)
})


test_that("semi-infinite and doubly infinite domains go through the maps", {
  shp <- c(0.5, 1, 2, 5, 50)
  f <- function(x, i) dgamma(x, shape = shp[i], rate = 1)
  expect_equal(quad_vec(f, 0, rep(Inf, 5)), rep(1, 5), tolerance = 1e-8)
  g <- function(x, i) x * dgamma(x, shape = shp[i], rate = 1)
  expect_equal(quad_vec(g, 0, rep(Inf, 5)), shp, tolerance = 1e-7)

  mu <- c(-3, 0, 10)
  sg <- c(0.5, 1, 4)
  h <- function(x, i) dnorm(x, mu[i], sg[i])
  expect_equal(quad_vec(h, rep(-Inf, 3), Inf), rep(1, 3), tolerance = 1e-9)
  m2 <- function(x, i) x^2 * dnorm(x, mu[i], sg[i])
  expect_equal(quad_vec(m2, rep(-Inf, 3), Inf), sg^2 + mu^2, tolerance = 1e-7)

  # a lower-infinite row: the left tail mass of a gaussian
  lt <- function(x, i) dnorm(x, 0, 1)
  expect_equal(quad_vec(lt, -Inf, 0), 0.5, tolerance = 1e-9)
})


test_that("rows of different kinds share one batch", {
  # finite, upper-infinite, lower-infinite and doubly infinite, together,
  # dispatched by i elementwise
  f <- function(x, i) {
    xv <- as.numeric(x)
    iv <- rep(i, times = ncol(x))
    out <- numeric(length(xv))
    out[iv == 1] <- dbeta(xv[iv == 1], 2, 3)
    out[iv == 2] <- dexp(xv[iv == 2])
    out[iv == 3] <- dnorm(xv[iv == 3], -1, 1)
    out[iv == 4] <- dcauchy(xv[iv == 4])
    matrix(out, nrow(x))
  }
  got <- quad_vec(f, lower = c(0, 0, -Inf, -Inf), upper = c(1, Inf, 0, Inf))
  want <- c(1, 1, pnorm(0, -1, 1), 1)
  expect_equal(got, want, tolerance = 1e-8)
})


test_that("quad_vec agrees with stats::integrate row by row", {
  set.seed(4)
  shp <- runif(20, 0.4, 30)
  rt <- runif(20, 0.2, 5)
  f <- function(x, i) sqrt(x) * dgamma(x, shape = shp[i], rate = rt[i])
  got <- quad_vec(f, 0, rep(Inf, 20))
  ref <- vapply(seq_len(20), function(j) {
    stats::integrate(function(x) sqrt(x) * dgamma(x, shp[j], rt[j]),
                     0, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
  expect_equal(got, ref, tolerance = 1e-7)
})


test_that("a row that cannot converge returns NA and is named", {
  # An endpoint singularity too harsh for the depth allowed: dbeta(x, 0.5, 1)
  # behaves like x^(-1/2) at zero, whose bisection error decays like
  # 2^(-depth/2), so at max_depth = 4 the error is SEEN but cannot be
  # brought under the budget -- and only ITS row fails. (A needle narrower
  # than every node spacing is not a usable probe here: no sampled rule can
  # estimate an error at points it never touches, so such a row converges
  # to the integral without the needle, exactly as stats::integrate does.)
  f <- function(x, i) {
    xv <- as.numeric(x)
    iv <- rep(i, times = ncol(x))
    out <- dunif(xv)
    sg <- iv == 2
    out[sg] <- dbeta(xv[sg], 0.5, 1)
    matrix(out, nrow(x))
  }
  expect_warning(
    got <- quad_vec(f, 0, c(1, 1), max_depth = 4L),
    "rows 2"
  )
  expect_equal(got[1], 1, tolerance = 1e-9)
  expect_true(is.na(got[2]))

  # and the same row converges once the depth is there: the refusal was
  # about the budget, not about the integrand
  g <- function(x, i) dbeta(x, 0.5, 1)
  expect_equal(quad_vec(g, 0, 1), 1, tolerance = 1e-8)
})


test_that("a corrupted rule is caught by what it integrates", {
  # the paired injection: a rule with one weight 5% wrong cannot integrate a
  # gamma density to one
  bad <- gauss_kronrod15()
  bad$wk[8] <- bad$wk[8] * 1.05
  f <- function(x, i) dgamma(x, shape = 3, rate = 1)
  got <- suppressWarnings(quad_vec(f, 0, Inf, rule = bad, max_depth = 6L))
  expect_false(isTRUE(all.equal(got, 1, tolerance = 1e-4)))
})


test_that("endpoints are validated", {
  expect_error(quad_vec(function(x, i) x, 1, 0), "smaller")
  expect_error(quad_vec(function(x, i) x, c(0, 2), c(1, 2)), "smaller")
})
