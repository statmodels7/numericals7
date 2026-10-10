# The thread policy of the toolkit, in one copy at the root. The object is
# PASSED as an argument from the fit entry points down to the compiled
# kernels: no package reads a setting that lives in another, which is the
# constraint the dependency graph imposes -- the kernels sit in packages
# statmodels7 cannot be imported by without a cycle.

#' How Many Threads and Worker Processes a Fit May Use
#'
#' @description
#' Constructs the parallelism policy an entry point such as
#' `statmodels7::statmod()` or `distributions7::fit_distrib()`
#' accepts through its `threads` argument. The default, one thread and
#' one process, takes exactly the sequential code path; a larger
#' `threads` lets the compiled per-observation kernels and the dense
#' assembly products run in parallel, and a larger `workers` distributes the
#' independent units of a fit over separate R processes: the folds of a
#' cross-validation and the combinations of a kinked path's product grid.
#'
#' @details
#' # Results independent of the thread count
#'
#' Every parallel region in the toolkit decomposes its work over the elements
#' of its output. Each accumulated value is therefore summed in full by one
#' thread, reductions are never split across threads, and a kernel returns the
#' same bits at any count. A fold's fit is a complete independent computation,
#' and the folds are collected in fold order however many processes ran them.
#'
#' Two properties of the drivers make that hold, and neither belongs to the
#' kernels themselves. A worker installs the calling thread's floating-point
#' environment before running its chunk, because without it some of the
#' platform's own math routines (among them R's `psigamma` at higher orders,
#' `bessel_k`, and `pgamma` and `pbeta` on the log scale) can return results
#' that differ in the last bits between threads. The sequential branch also
#' runs through the same compiled function as the parallel one, so a compiler
#' cannot optimize the two differently.
#'
#' One qualification. In the dense assembly products of \pkg{statmodels7},
#' raising the count first *engages* a threaded kernel where a BLAS expression
#' stood. That kernel is a second implementation of the same sum: bit-exact
#' against the reference BLAS R ships, and within the rounding of one dot
#' product against an optimized BLAS such as OpenBLAS or Accelerate, whose
#' accumulation order is its own.
#'
#' The guarantee is a design constraint on any kernel added later.
#'
#' # The sequential path
#'
#' At `threads = 1` no parallel region is entered and the process-level thread
#' setting is left unchanged. The code taken is the sequential one, not a
#' parallel path running on a single thread.
#'
#' Each kernel carries an internal threshold as well, measured where the cost of
#' opening a parallel region overtakes its gain. Below it the kernel stays
#' sequential regardless of the count.
#'
#' # Threads inside worker processes
#'
#' A fit running inside a worker process is sequential by construction, so
#' `workers = 4` opens up to four processes, each fitting on one thread, never
#' `4 * threads` of them.
#'
#' Worker processes are separate R sessions that load the installed packages,
#' which makes them safe for the S7 objects carried by a fold's fit; this is
#' the rule recorded by `optimizers7::multistart`. Starting a worker process
#' takes on the order of a second, so `workers` shortens a cross-validated path
#' only when each fold takes longer than that.
#'
#' # Printing
#'
#' `print()` shows the call that would rebuild the object, with `workers` shown
#' only when it is not 1, and appends `[sequential]` when both counts are 1. It
#' returns its argument invisibly.
#'
#' @param threads A single whole number, at least 1, and the default is `1`.
#'   The count is passed to the compiled kernels unchanged, and a kernel still
#'   runs sequentially on an input below its internal threshold. Any other
#'   value signals an error: a value below 1, a fractional value, `Inf`, `NA`,
#'   a character string and a vector of length other than one are all rejected
#'   with the same message. A value above the largest integer gives a coercion
#'   warning followed by R's error "missing value where TRUE/FALSE needed".
#' @param workers A single whole number, at least 1, and the default is `1`:
#'   how many R processes the independent units of a fit may use. `1` runs
#'   them in this process. Validated exactly as `threads` is.
#' @param x An `n_threads` object.
#' @param ... Unused, and present because `print()` requires it.
#'
#' @return `n_threads()` returns an object of class `"n_threads"`, a list of two
#'   components:
#'   \describe{
#'     \item{`threads`}{integer, the thread count, at least 1.}
#'     \item{`workers`}{integer, the worker-process count, at least 1.}
#'   }
#'   Read them with [thread_count()] and [worker_count()], which check the class
#'   first. `print()` returns its argument invisibly.
#'
#' @examples
#' # The default is sequential, and says so.
#' n_threads()
#'
#' # Threads within a fit, worker processes across folds, or both.
#' n_threads(4)
#' n_threads(1, workers = 4)
#' n_threads(4, workers = 2)
#'
#' # The counts come back out through the two readers.
#' p <- n_threads(4, workers = 2)
#' c(threads = thread_count(p), workers = worker_count(p))
#'
#' # A count must be a single whole number of at least one.
#' try(n_threads(0))
#' try(n_threads(2.5))
#'
#' @seealso [thread_count()] and [worker_count()] to read the counts,
#'   [local_threads()] to apply one for the duration of a fit.
#' @export
n_threads <- function(threads = 1, workers = 1) {
  whole <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 &&
      x == as.integer(x)
  }
  if (!whole(threads)) {
    stop("'threads' must be a single whole number, at least 1.",
         call. = FALSE)
  }
  if (!whole(workers)) {
    stop("'workers' must be a single whole number, at least 1.",
         call. = FALSE)
  }
  structure(list(threads = as.integer(threads),
                 workers = as.integer(workers)), class = "n_threads")
}

#' @rdname n_threads
#' @export
print.n_threads <- function(x, ...) {
  wrk <- if (x$workers != 1L) sprintf(", workers = %d", x$workers) else ""
  seq <- if (x$threads == 1L && x$workers == 1L) "  [sequential]" else ""
  cat(sprintf("n_threads(%d%s)%s\n", x$threads, wrk, seq))
  invisible(x)
}

#' Read the Thread Count Off a Policy
#'
#' @description
#' Returns the thread count an [n_threads()] object carries, as a single
#' integer. An entry point calls it once to validate the `threads` argument it
#' was given, then passes the plain count down to the kernels, which receive a
#' number and not the policy object.
#'
#' @param x An object returned by [n_threads()]. Any other value signals an
#'   error with a message naming the constructor, because a bare number passed
#'   to an entry point is the likely mistake and would otherwise be read as a
#'   policy.
#'
#' @return A single integer, at least 1.
#'
#' @examples
#' thread_count(n_threads())
#' thread_count(n_threads(8))
#'
#' # A bare count is not a policy.
#' try(thread_count(8))
#'
#' @seealso [n_threads()] for the object, [worker_count()] for the other count.
#' @export
thread_count <- function(x) {
  if (!inherits(x, "n_threads")) {
    stop("'threads' must be the object n_threads() returns, e.g. ",
         "threads = n_threads(4).", call. = FALSE)
  }
  x$threads
}

#' Read the Worker Count Off a Policy
#'
#' @description
#' Returns the worker-process count an [n_threads()] object carries, as a single
#' integer. A policy stored before `workers` existed has no such component, and
#' this function returns 1 for it, so that an object saved with an earlier fit
#' is read as sequential across folds.
#'
#' @param x An object returned by [n_threads()]. Any other value signals an
#'   error with a message naming the constructor.
#'
#' @return A single integer, at least 1. One for a policy carrying no `workers`
#'   component.
#'
#' @examples
#' worker_count(n_threads())
#' worker_count(n_threads(2, workers = 4))
#'
#' # A policy stored before workers existed still reads as sequential across
#' # folds. Only the readers tolerate one; print() on it throws.
#' worker_count(structure(list(threads = 3L), class = "n_threads"))
#'
#' @seealso [n_threads()] for the object, [thread_count()] for the other count.
#' @export
worker_count <- function(x) {
  if (!inherits(x, "n_threads")) {
    stop("'threads' must be the object n_threads() returns, e.g. ",
         "threads = n_threads(4).", call. = FALSE)
  }
  if (is.null(x$workers)) 1L else x$workers
}

#' Apply a Thread Policy for the Duration of One Fit
#'
#' @description
#' Sets \pkg{RcppParallel}'s thread count to the policy's for the calling
#' frame's lifetime and restores the previous process state when that frame
#' exits. Called once at the entry of a fit; the count itself still travels
#' to each kernel as an argument, and this call only sizes the worker pool from
#' which the parallel regions draw.
#'
#' @details
#' `RcppParallel::setThreadOptions()` writes a process-level variable, so a fit
#' that set it and returned would leave it moved for whatever runs next. The
#' restoration is therefore registered with `on.exit` in the caller's frame,
#' which fires however the fit leaves, an error included.
#'
#' At a count of 1 nothing is changed, so the sequential path does not touch
#' the parallel machinery.
#'
#' @param x An object returned by [n_threads()]. Any other value signals an
#'   error through [thread_count()].
#' @param frame The environment whose exit restores the previous setting,
#'   defaulting to the caller's. Pass an outer frame to hold the setting for
#'   longer than one call.
#'
#' @return The thread count, invisibly: a single integer, at least 1.
#'
#' @examples
#' f <- function() {
#'   local_threads(n_threads(2))
#'   Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
#' }
#' f()                                     # "2" inside the frame
#' Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "unset")   # restored
#'
#' @seealso [n_threads()], [thread_count()]
#' @export
local_threads <- function(x, frame = parent.frame()) {
  k <- thread_count(x)
  if (k > 1L) {
    old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
    RcppParallel::setThreadOptions(numThreads = k)
    cleanup <- if (is.na(old)) {
      function() Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
    } else {
      function() Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old)
    }
    do.call(base::on.exit,
            list(substitute(fn(), list(fn = cleanup)), add = TRUE),
            envir = frame)
  }
  invisible(k)
}
