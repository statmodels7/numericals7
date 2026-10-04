# bessel_k() is R's bessel_k() compiled without its warnings and with its work
# array on the stack, so the check is identity with besselK() over a grid of
# orders and arguments that reaches both ends of the range, scaled and not.

test_that("the values are besselK()'s, bit for bit", {
  nu <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 7.3, 20, 63.5, 64, 100.2)
  x <- c(10^seq(-12, 3, by = 0.125), 704, 705.3, 706, 1e4)
  g <- expand.grid(x = x, nu = c(-rev(nu), nu))
  for (sc in c(FALSE, TRUE)) {
    expect_identical(suppressWarnings(bessel_k(g$x, g$nu, expon.scaled = sc)),
                     suppressWarnings(besselK(g$x, g$nu, expon.scaled = sc)))
  }
})

test_that("the edges and the thread count", {
  expect_true(is.nan(bessel_k(-1, 1)))
  expect_identical(bessel_k(c(NA, 1), c(1, NA)), c(NA_real_, NA_real_))
  expect_identical(bessel_k(0, 1), Inf)
  set.seed(2)
  x <- exp(runif(5000, -5, 5))
  nu <- runif(5000, 0, 6)
  expect_identical(bessel_k(x, nu, expon.scaled = TRUE, threads = 4L),
                   bessel_k(x, nu, expon.scaled = TRUE))
})
