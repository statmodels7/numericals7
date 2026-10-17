# fd_tensor() and fd_known(): one tensor-product stencil per component, from
# the highest declared analytic quantity.

counter <- function(f) {
  env <- new.env()
  env$n <- 0L
  list(f = function(x) { env$n <- env$n + 1L; f(x) }, calls = function() env$n)
}

test_that("the one-dimensional weights are Fornberg's, exactly symmetric", {
  for (a in c(2L, 4L, 6L)) for (m in 1:6) {
    pl <- numericals7:::fd_tensor_plan_cpp(matrix(m, 1, 1),
                                           list(matrix(0L, 1, 1)), a)
    off <- fd_offsets(m, a)$central
    w <- fd_weights(off, m)
    # the Vandermonde solve leaves the centre of an odd order at about 1e-12
    keep <- off != 0 | m %% 2L == 0L
    expect_equal(pl$term_w, w[keep], tolerance = 1e-10)
    expect_identical(as.integer(pl$node_off), as.integer(off[keep]))
    # an odd order drops the centre, which is an exact zero
    expect_identical(0L %in% pl$node_off, m %% 2L == 0L)
    expect_identical(pl$term_w, rev(pl$term_w) * (-1)^m)
  }
  # the seven-point fourth derivative, written out
  pl <- numericals7:::fd_tensor_plan_cpp(matrix(4L, 1, 1),
                                         list(matrix(0L, 1, 1)), 4L)
  expect_equal(pl$term_w, c(-1/6, 2, -13/2, 28/3, -13/2, 2, -1/6),
               tolerance = 1e-15)
})

test_that("a polynomial within the stencils' degree is reproduced exactly", {
  # degree at most 2R in each differenced direction: the stencils are exact,
  # and a wide step leaves no rounding to speak of
  P <- function(x) x$u^4 * x$v^3 + 2 * x$u^2 * x$v + x$v^4 - 3 * x$u^3
  x <- list(u = c(-0.7, 0.2, 1.3), v = c(0.4, -1.1, 0.9))
  u <- x$u; v <- x$v
  # at accuracy 4 every stencil up to order five is exact to degree four
  exact <- list(
    u_u_u = 24 * u * v^3 - 18, u_u_v = 36 * u^2 * v^2 + 4,
    u_v_v = 24 * u^3 * v, v_v_v = 6 * u^4 + 24 * v)
  got <- fd_tensor(P, x, order = 3, accuracy = 4, h_rel = 0.5)
  for (nm in names(exact)) expect_equal(got[[nm]], exact[[nm]], tolerance = 1e-10)
  exact5 <- list(u_u_u_u_u = 0 * u, u_u_u_u_v = 72 * v^2,
                 u_u_u_v_v = 144 * u * v, u_u_v_v_v = 72 * u^2,
                 u_v_v_v_v = 0 * u, v_v_v_v_v = 0 * u)
  got <- fd_tensor(P, x, order = 5, accuracy = 4, h_rel = 0.5)
  for (nm in names(exact5)) {
    expect_equal(got[[nm]], exact5[[nm]], tolerance = 1e-9)
  }
  # at accuracy 2 a first difference is exact to degree two
  Q <- function(x) x$u^2 * x$v^2 + 3 * x$u * x$v^2 - x$u^2 * x$v
  got <- fd_tensor(Q, x, order = 3, accuracy = 2, h_rel = 0.5)
  expect_equal(got$u_u_v, 4 * v - 2, tolerance = 1e-10)
  expect_equal(got$u_v_v, 4 * u + 6, tolerance = 1e-10)
  got <- fd_tensor(Q, x, order = 4, accuracy = 2, h_rel = 0.5)
  expect_equal(got$u_u_v_v, 4 + 0 * u, tolerance = 1e-10)
})

test_that("the components follow tuple_indices() and are named by coordinate", {
  f <- function(x) x$a * x$b * x$c
  x <- list(a = 1, b = 2, c = 3)
  for (k in 1:5) {
    got <- fd_tensor(f, x, order = k)
    want <- vapply(tuple_indices(3, k),
                   function(i) paste(c("a", "b", "c")[i], collapse = "_"),
                   character(1))
    expect_identical(names(got), want)
  }
  got <- fd_tensor(f, x, components = rbind(c(1, 1, 1), c(0, 2, 0)))
  expect_identical(names(got), c("a_b_c", "b_b"))
  expect_equal(got$a_b_c, 1, tolerance = 1e-8)
  expect_equal(got$b_b, 0, tolerance = 1e-8)
})

test_that("a smooth function is differentiated to the accuracy of the rule", {
  f <- function(x) exp(x$a) * sin(x$b) + x$a^2 * cos(x$b)
  x <- list(a = c(-0.5, 0.3, 1.1), b = c(0.2, -0.8, 1.4))
  a <- x$a; b <- x$b
  # partials of exp(a) sin(b) + a^2 cos(b)
  ex <- list(
    a_a_a_b = exp(a) * cos(b),
    a_a_b_b = -exp(a) * sin(b) - 2 * cos(b),
    a_a_a_a_b = exp(a) * cos(b),
    a_b_b_b_b = exp(a) * sin(b) + 2 * a * cos(b))
  err <- list()
  for (acc in c(2L, 4L)) {
    d4 <- fd_tensor(f, x, order = 4, accuracy = acc)
    d5 <- fd_tensor(f, x, order = 5, accuracy = acc)
    err[[acc]] <- max(abs(unlist(c(d4[names(ex)[1:2]], d5[names(ex)[3:4]])) -
                            unlist(ex)))
  }
  expect_lt(err[[2]], 1e-2)
  expect_lt(err[[4]], 1e-3)
  expect_lt(err[[4]], err[[2]])
})

test_that("each component starts from the declared quantity of least residual", {
  f <- function(x) x$a^3 * exp(x$b)
  g <- function(x) list(3 * x$a^2 * exp(x$b), x$a^3 * exp(x$b))
  H <- function(x) list(6 * x$a * exp(x$b), x$a^3 * exp(x$b),
                        3 * x$a^2 * exp(x$b))
  x <- list(a = c(0.5, 2), b = 0.3)
  got <- fd_tensor(f, x, order = 3,
                   known = list(fd_known(g, order = 1), fd_known(H, order = 2)))
  src <- attr(got, "fd_source")
  expect_identical(src$component, c("a_a_a", "a_a_b", "a_b_b", "b_b_b"))
  expect_identical(src$source, rep(2L, 4))
  expect_identical(src$residual, rep(1L, 4))
  # a_a_b differences the a_a entry along b; the most repeated coordinate
  # goes to the analytic part
  expect_identical(src$index, c(1L, 1L, 3L, 2L))
  expect_equal(got$a_a_b, 6 * x$a * exp(x$b), tolerance = 1e-8)

  # a quantity that leaves the full order serves only where nothing else does
  got <- fd_tensor(f, x, components = rbind(c(1, 0), c(3, 1)),
                   known = fd_known(H, order = 2))
  expect_identical(attr(got, "fd_source")$source, c(0L, 1L))
  expect_identical(attr(got, "fd_source")$residual, c(1L, 2L))
})

test_that("equal residuals go to the larger declaration, then the first one", {
  f <- function(x) x$y^2 * x$m^2
  x <- list(y = 0.7, m = 1.3)
  # both leave a residual of one; the second declares more derivatives
  dm <- fd_known(function(x) 2 * x$y^2 * x$m, indices = c(0, 1))
  dy <- fd_known(function(x) list(2 * x$y * x$m^2, 2 * x$m^2),
                 indices = rbind(c(1, 0), c(2, 0)))
  got <- fd_tensor(f, x, components = rbind(c(1, 1)), known = list(dm, dy))
  expect_identical(attr(got, "fd_source")$source, 2L)
  expect_equal(got$y_m, 4 * x$y * x$m, tolerance = 1e-8)
  one <- fd_known(function(x) 2 * x$y * x$m^2, indices = c(1, 0))
  # equal declarations: the first is taken
  one_b <- fd_known(function(x) 2 * x$y^2 * x$m, indices = c(0, 1))
  got <- fd_tensor(f, x, components = rbind(c(1, 1)), known = list(one, one_b))
  expect_identical(attr(got, "fd_source")$source, 1L)
  got <- fd_tensor(f, x, components = rbind(c(1, 1)), known = list(one_b, one))
  expect_identical(attr(got, "fd_source")$source, 1L)
  expect_equal(got$y_m, 4 * x$y * x$m, tolerance = 1e-8)
})

test_that("a node is evaluated once, whatever the components that share it", {
  cnt <- counter(function(x) stats::dnorm(1.2, x$mu, x$sigma, log = TRUE))
  invisible(fd_tensor(cnt$f, list(mu = 0.3, sigma = 1.2), order = 5,
                      accuracy = 2))
  # the distinct points of the six fifth-order tensor stencils
  expect_identical(cnt$calls(), 24L)
  cnt <- counter(function(x) x$a^2)
  g <- fd_known(function(x) list(2 * x$a, 0 * x$a), order = 1)
  ncall <- 0L
  g@fun <- function(x) { ncall <<- ncall + 1L; list(2 * x$a, 0 * x$a) }
  invisible(fd_tensor(cnt$f, list(a = 1, b = 2), order = 2, known = g,
                      accuracy = 2))
  expect_identical(cnt$calls(), 0L)
  expect_identical(ncall, 4L)
})

test_that("the step keeps every node inside the bounds", {
  f <- function(x) {
    if (any(x$s <= 0)) stop("left the domain")
    log(x$s)^2 + x$m * x$s
  }
  x <- list(m = 0.5, s = c(1e-3, 0.5, 3))
  s <- x$s
  exact <- (4 * log(s) - 6) / s^3
  for (acc in c(2L, 4L)) {
    got <- fd_tensor(f, x, order = 3, accuracy = acc,
                     bounds = list(s = c(0, Inf)))
    expect_true(all(is.finite(got$s_s_s)))
    expect_equal(got$s_s_s[2:3], exact[2:3], tolerance = 1e-4)
    # next to the bound the step is a quarter of the distance, and the
    # truncation is of that relative size
    expect_lt(abs(got$s_s_s[1] / exact[1] - 1), 0.25)
  }
  # without the bound the stencil leaves the domain
  expect_error(fd_tensor(f, x, order = 3), "left the domain")
})

test_that("the length of the values sets n when the coordinates are scalars", {
  y <- c(-1, 0, 2, 3)
  f <- function(x) stats::dnorm(y, x$mu, 1, log = TRUE)
  got <- fd_tensor(f, list(mu = 0.4), order = 1)
  expect_length(got$mu, 4L)
  expect_equal(got$mu, y - 0.4, tolerance = 1e-10)
  # a known quantity of length 1 is recycled to the common length
  got <- fd_tensor(f, list(mu = 0.4), components = rbind(1, 2),
                   known = fd_known(function(x) -1, order = 2))
  expect_identical(attr(got, "fd_source")$source, c(0L, 1L))
  expect_identical(got$mu_mu, rep(-1, 4))
})

test_that("the combination does not depend on the thread count", {
  set.seed(2)
  n <- 50000L
  vals <- lapply(1:6, function(i) rnorm(n))
  vals[[3]] <- 0.7
  ts <- c(0L, 3L, 6L)
  tv <- c(0L, 1L, 2L, 3L, 4L, 5L)
  tw <- c(0.5, -1.25, 2, 1 / 3, -0.1, 7)
  dens <- list(runif(n) + 0.5, 0.01)
  a <- numericals7:::fd_tensor_combine_cpp(vals, ts, tv, tw, dens, n, 1L)
  b <- numericals7:::fd_tensor_combine_cpp(vals, ts, tv, tw, dens, n, 4L)
  expect_identical(a, b)
  expect_equal(a[[1]], (0.5 * vals[[1]] - 1.25 * vals[[2]] + 2 * 0.7) / dens[[1]])
})

test_that("the plan is cached and returned unchanged", {
  f <- function(x) x$a * x$b^2
  x <- list(a = 1, b = 2)
  r1 <- fd_tensor(f, x, order = 3)
  r2 <- fd_tensor(f, x, order = 3)
  expect_identical(r1, r2)
})

test_that("malformed requests are rejected", {
  f <- function(x) x$a
  x <- list(a = 1)
  expect_error(fd_tensor(f, x), "exactly one of")
  expect_error(fd_tensor(f, x, order = 2, components = matrix(1L)), "exactly one of")
  expect_error(fd_tensor(f, x, order = 2, accuracy = 3), "even")
  expect_error(fd_tensor(f, x, order = 0), "positive whole")
  expect_error(fd_tensor(f, list(1), order = 1), "distinct names")
  expect_error(fd_tensor(f, list(a = 1:2, b = 1:3), order = 1), "length 1")
  expect_error(fd_tensor(f, x, components = matrix(1L, 1, 2)), "one column")
  expect_error(fd_tensor(f, x, order = 1, known = list(1)), "fd_known")
  expect_error(fd_tensor(f, x, order = 2,
                         known = fd_known(f, indices = c(1, 0))),
               "2 columns for 1")
  expect_error(fd_tensor(f, x, order = 2,
                         known = fd_known(function(x) list(1, 2), order = 1)),
               "returned 2 values for 1")
  expect_error(fd_tensor(f, x, order = 1, bounds = list(a = c(1, 0))),
               "lower < upper")
  expect_error(fd_known(f), "exactly one of")
  expect_error(fd_known(f, order = 1, indices = 1), "exactly one of")
  expect_error(fd_known(f, indices = c(-1, 0)), "non-negative")
})

test_that("a request whose key exceeds the length of a symbol is planned", {
  # 715 components over ten coordinates give a key of about 14 kB, more than
  # the 10000 bytes R allows a variable name.
  a <- seq(0.1, 1, length.out = 10)
  nm <- paste0("x", 1:10)
  x <- stats::setNames(as.list(seq(-0.5, 0.4, by = 0.1)), nm)
  f <- function(x) exp(sum(a * unlist(x)))
  comp <- numericals7:::counts_of_order(10L, 4L)
  k3 <- numericals7:::counts_of_order(10L, 3L)
  d3 <- fd_known(function(x) {
    v <- f(x)
    lapply(seq_len(nrow(k3)), function(r) prod(a^k3[r, ]) * v)
  }, indices = k3)
  got <- fd_tensor(f, x, components = comp, known = d3)
  want <- apply(comp, 1L, function(m) prod(a^m)) * f(x)
  expect_equal(unname(unlist(got)), want, tolerance = 1e-8)
  expect_identical(fd_tensor(f, x, components = comp, known = d3), got)
})
