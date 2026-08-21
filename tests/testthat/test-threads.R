# The thread policy: one object, constructed here at the root and passed as
# an argument from the fit entry points down to the kernels.

test_that("n_threads() validates its count", {
  expect_s3_class(n_threads(), "n_threads")
  expect_identical(thread_count(n_threads()), 1L)
  expect_identical(thread_count(n_threads(4)), 4L)
  expect_identical(thread_count(n_threads(2L)), 2L)
  expect_error(n_threads(0), "at least 1")
  expect_error(n_threads(-1), "at least 1")
  expect_error(n_threads(1.5), "whole number")
  expect_error(n_threads(c(1, 2)), "single")
  expect_error(n_threads("4"), "whole number")
  expect_error(n_threads(Inf), "whole number")
})

test_that("thread_count() rejects anything but the object", {
  expect_error(thread_count(4), "n_threads")
  expect_error(thread_count(NULL), "n_threads")
  expect_error(thread_count(list(threads = 2L)), "n_threads")
})

test_that("workers is carried, validated and read", {
  expect_identical(worker_count(n_threads()), 1L)
  expect_identical(worker_count(n_threads(2, workers = 4)), 4L)
  expect_identical(thread_count(n_threads(2, workers = 4)), 2L)
  expect_error(n_threads(1, workers = 0), "at least 1")
  expect_error(n_threads(1, workers = 1.5), "whole number")
  expect_error(worker_count(4), "n_threads")
  # an object stored before `workers` existed answers 1
  old <- structure(list(threads = 3L), class = "n_threads")
  expect_identical(worker_count(old), 1L)
})

test_that("print method names the sequential default", {
  expect_output(print(n_threads()), "sequential")
  expect_output(print(n_threads(3)), "n_threads\\(3\\)")
})

test_that("local_threads() sets for the frame and restores on exit", {
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  f <- function() {
    local_threads(n_threads(2))
    Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
  }
  expect_identical(f(), "2")
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS",
                              unset = NA_character_), old)

  # a previously SET value is restored, not unset
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "7")
  expect_identical(f(), "2")
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS"), "7")
  if (is.na(old)) Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS") else
    Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old)
})

test_that("local_threads() at 1 touches nothing", {
  Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
  f <- function() {
    local_threads(n_threads(1))
    Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "unset")
  }
  expect_identical(f(), "unset")
})

# The one elementwise kernel of this package that is threaded. Every branch
# of log I is the package's own arithmetic -- a series or a uniform
# asymptotic expansion, no call into Rmath -- so element i is computed and
# written by one thread and the answer is bit-identical at any count.
# Measured at n = 20000, the dearest branch (small argument) is 167 ms
# sequential against 52 at eight threads, and it is 80.8 per cent of a von
# Mises fit whose concentration is modelled.
test_that("the log I kernel does not depend on the thread count", {
  set.seed(41)
  n <- 5000
  for (rg in list(c(0.1, 5), c(0.1, 50), c(1, 1e4), c(1e4, 1e6))) {
    x <- runif(n, rg[1], rg[2])
    nu <- runif(n, 0, 3)
    ref <- log_bessel_i(x, nu, 1L)
    for (k in c(2L, 3L, 5L)) expect_identical(log_bessel_i(x, nu, k), ref)
  }
})

test_that("the threaded kernel still matches its R twin", {
  set.seed(42)
  x <- runif(300, 0.1, 40)
  nu <- runif(300, 0, 4)
  expect_equal(log_bessel_i(x, nu, 3L), .log_bessel_i_r(x, nu),
               tolerance = 1e-12)
})

# log K takes no count, and that is a refusal rather than an omission: its
# hybrid branch calls R's own scaled besselK, which can raise a warning, and
# a warning from a worker thread ends the session. An argument that would be
# read by nobody is worse than one that is absent.
test_that("log_bessel_k takes no thread count", {
  expect_false("threads" %in% names(formals(log_bessel_k)))
  expect_true("threads" %in% names(formals(log_bessel_i)))
})
