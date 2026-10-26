# log I, log K and their derivatives against two stored references
# (fixtures/log_bessel_ref.rds): mpmath at 50 digits on a grid that covers
# every branch of the input-range table (every series point of a 4300-point
# grid, 150 points of each other branch), and the values that numericals7
# 0.23.5 returned before the series was summed by its ratio and the u_k
# polynomials by Horner. Errors are in units of eps * max(1, |reference|).
# mpmath's besselk cancels two values of size exp(nu eta) at a large order
# and a large argument (nu = 925, x = 649 gives +61.9 for log K = -67.5), so
# the comparisons of log K leave out nu > 400 with x > 400.

ref <- readRDS(test_path("fixtures", "log_bessel_ref.rds"))
g <- ref$grid
eps <- .Machine$double.eps
err <- function(a, r) abs(a - r) / (eps * pmax(1, abs(r)))

test_that("log I and log K agree with mpmath on every branch", {
  li <- log_bessel_i(g$x, g$nu)
  e <- err(li, g$mp_lI)
  expect_lt(max(e[g$branch == 0L]), 64)
  expect_lt(max(e), 1024)
  k_ok <- is.finite(g$mp_lK) & !(g$nu > 400 & g$x > 400)
  expect_gt(sum(k_ok), 1800)
  expect_lt(max(err(log_bessel_k(g$x, g$nu), g$mp_lK)[k_ok]), 256)
})

test_that("log I is not less accurate than in 0.23.5", {
  e_new <- err(log_bessel_i(g$x, g$nu), g$mp_lI)
  e_old <- err(g$old_lI, g$mp_lI)
  # points where the new value is worse by more than a factor 2 and 4 eps
  expect_lte(sum(e_new > 2 * e_old & e_new > 4), 10)
  expect_lte(max(e_new[g$branch == 0L]), max(e_old[g$branch == 0L]) + 4)
})

test_that("log I and log K stay within a few eps of 0.23.5", {
  expect_lt(max(err(log_bessel_i(g$x, g$nu), g$old_lI)), 64)
  expect_lt(max(err(log_bessel_k(g$x, g$nu), g$old_lK)), 16)
  # The stored values were computed on Windows. At a large order the
  # uniform expansion carries an absolute error of about nu * eps (the
  # product nu * eta), so another platform's libm moves the result by that
  # much: on macOS 0.23.5 itself differs from the stored values by 1715 eps
  # on the scale above (x = 191.7, nu = 285.7, log I = 0.437). The random
  # set, which reaches nu = 1e4, is therefore compared on the scale
  # max(1, |reference|, nu).
  set.seed(7)
  xr <- exp(runif(5000, log(1e-4), log(1e6)))
  vr <- exp(runif(5000, log(1e-3), log(1e4)))
  err_nu <- function(a, r) abs(a - r) / (eps * pmax(1, abs(r), vr))
  expect_lt(max(err_nu(log_bessel_i(xr, vr), ref$random$lI)), 64)
  expect_lt(max(err_nu(log_bessel_k(xr, vr), ref$random$lK), na.rm = TRUE), 64)
})

test_that("the derivatives agree with mpmath", {
  di <- log_bessel_i_derivs(g$x, g$nu)
  dk <- log_bessel_k_derivs(g$x, g$nu)
  # d1 is exp of a difference of two logs, so its relative error scales
  # with the size of the log
  rel1 <- function(a, r, l) abs(a - r) / (abs(r) * eps * pmax(1, abs(l)))
  expect_lt(max(rel1(di$d1, g$mp_dI1, g$mp_lI)), 128)
  k_ok <- is.finite(g$mp_lK) & !(g$nu > 400 & g$x > 400)
  expect_lt(max(rel1(dk$d1[k_ok], g$mp_dK1[k_ok], g$mp_lK[k_ok])), 1024)
  # orders 2 to 4 come from the Bessel equation, which cancels terms of
  # order one; the comparison is kept to x <= 100, where the loss is small
  s <- g$x <= 100
  bound <- c(4e-9, 1e-7, 2e-5)
  for (q in 2:4) {
    r <- g[[paste0("mp_dI", q)]]
    ok <- s & is.finite(r) & r != 0
    expect_lt(max(abs(di[[q + 1]][ok] - r[ok]) / abs(r[ok])), bound[q - 1])
  }
})

test_that("the u_k table has the exponents that the Horner evaluation assumes", {
  for (k in seq_along(numericals7:::.bessel_uk)) {
    expect_identical(numericals7:::.bessel_uk[[k]]$e, as.numeric(seq(k, 3 * k, by = 2)))
  }
  bad <- numericals7:::.bessel_uk
  bad[[2]]$e <- bad[[2]]$e + 1
  expect_error(numericals7:::lb_set_uk_cpp(bad), "u_2")
  numericals7:::lb_set_uk_cpp(numericals7:::.bessel_uk)
})
