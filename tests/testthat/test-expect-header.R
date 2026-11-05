# Twins of inst/include/numericals7/n7_expect.h: the compiled routines
# return the same doubles as the R functions they transcribe.

test_that("the compiled mean() is R's", {
  set.seed(1)
  xs <- list(rnorm(1000), c(1, 2^-53, 2^-53), runif(7) * 1e300,
             c(.Machine$double.xmax, .Machine$double.xmax, -1e308),
             c(1, NaN, 2), c(Inf, 1), numeric(0) + 3)
  for (x in xs) expect_identical(n7_mean_probe(x), mean(x))
})

test_that("the compiled %% is R's", {
  set.seed(2)
  x <- c(rnorm(200, sd = 10), -pi, pi, 0, -0, 1e17, -1e-300, 2 * pi, -2 * pi)
  expect_identical(n7_fmod_probe(x, 2 * pi), x %% (2 * pi))
})

test_that("the probe distances are series_vec()'s", {
  s <- n7_series_steps_probe()
  block <- 64L
  nj <- ceiling(log(2^31 / block) / log1p(1 / 128))
  expect_identical(s, c(0, unique(floor(block * (1 + 1 / 128)^(0:nj)))))
})

test_that("a compiled row of quad_vec() is quad_vec()'s row", {
  cases <- list(
    list(f = function(x, i) dnorm(x), lo = c(-Inf, -1, 0), up = c(-1, 0, Inf)),
    list(f = function(x, i) dbeta(x, 0.8, 1.2), lo = 0, up = 1),
    list(f = function(x, i) 1 / sqrt(abs(x)), lo = c(-1, -1), up = c(1, 2)),
    list(f = function(x, i) 1 / abs(x), lo = -1, up = 1),
    list(f = function(x, i) sin(50 * x) * exp(-x^2), lo = -Inf, up = Inf),
    list(f = function(x, i) exp(-x) * x^3, lo = 0, up = Inf),
    list(f = function(x, i) ifelse(x == 0, NaN, 1), lo = -1, up = 1),
    list(f = function(x, i) 1 / (1 + x^2), lo = -Inf, up = 2)
  )
  for (cs in cases) {
    for (pr in list(c(1e-10, 1e-8, 4096), c(1e-12, 1e-10, 64))) {
      ref <- suppressWarnings(quad_vec(cs$f, cs$lo, cs$up, atol = pr[1], rtol = pr[2],
                                       max_panels = pr[3]))
      got <- n7_quad_rows_probe(cs$f, rep_len(cs$lo, length(ref)), rep_len(cs$up, length(ref)),
                                pr[1], pr[2], as.integer(pr[3]))
      expect_identical(got, ref)
    }
  }
})

test_that("a compiled row of series_vec() is series_vec()'s row", {
  lam <- c(0.5, 30, 1e4, 3)
  t1 <- function(k, i) k^2 * dpois(k, lam[i])
  expect_identical(n7_series_rows_probe(t1, 4L, 0L, 1e-12, 1e-10, 100000L),
                   series_vec(t1, 4L))
  t2 <- function(k, i) dnbinom(k, size = c(0.3, 50)[i], mu = c(200, 2)[i]) * log1p(k)
  expect_identical(n7_series_rows_probe(t2, 2L, 0L, 1e-12, 1e-10, 100000L),
                   series_vec(t2, 2L))
  t3 <- function(k, i) 0 * k
  expect_identical(n7_series_rows_probe(t3, 1L, 0L, 1e-12, 1e-10, 100000L),
                   series_vec(t3, 1L))
  t4 <- function(k, i) 1 / (k + 1)
  expect_identical(n7_series_rows_probe(t4, 1L, 0L, 1e-12, 1e-10, 2000L),
                   suppressWarnings(series_vec(t4, 1L, max_terms = 2000L)))
})

test_that("the compiled uniroot() finds uniroot()'s root", {
  fs <- list(function(x) x - 0.1, function(x) pnorm(x) - 0.3,
             function(x) x^3 - 2 * x - 5, function(x) log(x) + 3,
             function(x) exp(x) - 1e5)
  ivs <- list(c(-1, 1), c(-5, 5), c(0, 4), c(1e-6, 2), c(-3, 20))
  for (k in seq_along(fs)) for (tol in c(1e-14, .Machine$double.eps^0.5)) {
    ref <- uniroot(fs[[k]], ivs[[k]], tol = tol, maxiter = 1000)$root
    expect_identical(n7_uniroot_probe(fs[[k]], ivs[[k]][1], ivs[[k]][2], tol, 1000L), ref)
  }
})
