# The vectorized quadrature. The references are closed forms and
# stats::integrate, which is an independent implementation (QUADPACK through
# a different path, one theta at a time).

test_that("the Gauss-Kronrod pair is pinned by its defining property", {
  r <- gauss_kronrod15()
  expect_length(r$nodes, 15L)
  # a rule on [-1, 1] has weights summing to its length, and at full precision
  # to within the rounding of the sum: a few units in the last place, where the
  # fifteen-decimal table this replaced was 27 units short
  expect_lt(abs(sum(r$wk) - 2), 4 * .Machine$double.eps)
  expect_lt(abs(sum(r$wg) - 2), 4 * .Machine$double.eps)
  old_wkh <- c(0.022935322010529, 0.063092092629979, 0.104790010322250,
               0.140653259715525, 0.169004726639267, 0.190350578064785,
               0.204432940075298, 0.209482141084728)
  old_wk <- c(old_wkh[1:7], old_wkh[8], rev(old_wkh[1:7]))
  expect_gt(abs(sum(old_wk) - 2), 4 * .Machine$double.eps)
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
  # Asked for the accuracy it is held to. At the default rtol = 1e-8 the two
  # rows with an endpoint singularity stop 1.2e-9 and 2.6e-9 from 1, inside
  # that request, and this comparison passed only because expect_equal()
  # averages over the elements that are not exactly equal: with the
  # fifteen-decimal constants all five carried a bias and the average was
  # 7.5e-10, with the full-precision ones two rows are exactly 1 and the
  # average over the other three is 1.3e-9.
  got <- quad_vec(f, 0, rep(1, 5), rtol = 1e-10)
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


test_that("a relative budget below the rule's floor is refused before anything is evaluated", {
  calls <- 0L
  f <- function(x, i) {
    calls <<- calls + 1L
    matrix(1, nrow(x), ncol(x))
  }
  fl <- quad_floor(gauss_kronrod15())
  # the default rule's two weight sums agree up to their rounding, so its floor
  # is one or two units in the last place: the computed sums agree to the last
  # bit on x86_64 (2.2e-16) and differ by one unit in the last place of 2 on
  # the arm64 macOS runner (4.4e-16). A bound between that floor and the
  # 3.7e-15 of the fifteen-decimal table below asserts the rule, not the bit.
  bound <- 1e-15
  expect_lt(fl, bound)
  expect_lt(2 * .Machine$double.eps, bound)

  # with no absolute budget the call is refused, the floor is named, and the
  # integrand is never evaluated
  expect_error(quad_vec(f, 0, 1, atol = 0, rtol = fl / 2),
               sprintf("below the floor of the quadrature rule, %.3g", fl),
               fixed = TRUE)
  expect_identical(calls, 0L)
  # an absolute budget can end the refinement, so rtol = 0 stays valid beside one
  expect_equal(quad_vec(f, 0, 1, atol = 1e-10, rtol = 0), 1, tolerance = 1e-14)
  expect_gt(calls, 0L)

  # the floor is a property of the rule: the fifteen-decimal table this
  # package carried before 0.14.0 has its own, and a budget it cannot meet is
  # refused with it while the full-precision rule meets the same budget
  old_rule <- local({
    xh <- c(0.991455371120813, 0.949107912342759, 0.864864423359769,
            0.741531185599394, 0.586087235467691, 0.405845151377397,
            0.207784955007898, 0)
    wkh <- c(0.022935322010529, 0.063092092629979, 0.104790010322250,
             0.140653259715525, 0.169004726639267, 0.190350578064785,
             0.204432940075298, 0.209482141084728)
    wgh <- c(0, 0.129484966168870, 0, 0.279705391489277,
             0, 0.381830050505119, 0, 0.417959183673469)
    list(nodes = c(-xh[1:7], xh[8], rev(xh[1:7])),
         wk = c(wkh[1:7], wkh[8], rev(wkh[1:7])),
         wg = c(wgh[1:7], wgh[8], rev(wgh[1:7])))
  })
  fl_old <- quad_floor(old_rule)
  expect_gt(fl_old, bound)
  expect_error(quad_vec(f, 0, 1, atol = 0, rtol = 1e-15, rule = old_rule),
               sprintf("%.3g", fl_old), fixed = TRUE)
  expect_equal(quad_vec(f, 0, 1, atol = 0, rtol = 1e-15), 1, tolerance = 1e-15)
})


test_that("a row whose error estimate does not fall with bisection stops at max_panels", {
  # sin(1e9 x) on [0, 1]: until a panel is narrower than a period of 6.3e-9
  # every panel carries a share of the error, so the summed estimate does not
  # fall and the panel count doubles at every pass. Uncapped that is the
  # memory exhaustion the cap exists for, and nothing here runs uncapped.
  largest <- 0L
  f <- function(x, i) {
    largest <<- max(largest, nrow(x))
    sin(1e9 * x)
  }
  expect_warning(got <- quad_vec(f, c(0, 0), 1, max_panels = 256L),
                 "max_panels = 256")
  expect_true(all(is.na(got)))
  # a row splits at most the panels it holds, so f never sees more than two
  # caps' worth of fresh panels per row
  expect_lt(largest, 2L * 2L * 256L)

  # the cap is per row: a row that converges keeps the value it has alone
  g <- function(x, i) {
    xv <- as.numeric(x)
    iv <- rep(i, times = ncol(x))
    out <- cos(20 * xv)
    out[iv == 2] <- sin(1e9 * xv[iv == 2])
    matrix(out, nrow(x))
  }
  both <- suppressWarnings(quad_vec(g, c(0, 0), 1, max_panels = 256L))
  alone <- quad_vec(function(x, i) cos(20 * x), 0, 1, max_panels = 256L)
  # a tolerance and not an identity: the two calls multiply matrices of
  # different shape, and a BLAS may accumulate them in a different order
  expect_equal(both[1], alone, tolerance = 1e-13)
  expect_true(is.na(both[2]))

  # and a row that converges within the cap is untouched by it: cos(1000 x)
  # needs between 128 and 256 panels
  h <- function(x, i) cos(1000 * x)
  expect_identical(quad_vec(h, 0, 1, max_panels = 4096L),
                   quad_vec(h, 0, 1, max_panels = 1e6))
  expect_equal(quad_vec(h, 0, 1), sin(1000) / 1000, tolerance = 1e-6)

  expect_error(quad_vec(h, 0, 1, max_panels = 0), "max_panels")
  expect_error(quad_vec(h, 0, 1, max_panels = c(10, 20)), "max_panels")
})
