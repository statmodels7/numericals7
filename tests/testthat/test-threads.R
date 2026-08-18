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
