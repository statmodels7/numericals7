# The thread policy of the toolkit, in one copy at the root. The object is
# PASSED as an argument from the fit entry points down to the compiled
# kernels: no package reads a setting that lives in another, which is the
# constraint the dependency graph imposes -- the kernels sit in packages
# statmodels7 cannot be imported by without a cycle.

#' How Many Threads and Worker Processes a Fit May Use
#'
#' @description
#' Constructs the parallelism policy an entry point such as
#' \code{statmodels7::statmod()} or \code{distributions7::fit_distrib()}
#' accepts through its \code{threads} argument. The default, one thread and
#' one process, takes exactly the sequential code path; a larger
#' \code{threads} lets the compiled per-observation kernels and the dense
#' assembly products run in parallel, and a larger \code{workers} fans the
#' independent fits of a cross-validation's folds out over separate R
#' processes.
#'
#' @details
#' The result does not depend on \code{threads} or on \code{workers}, bit
#' for bit, at any count. Every parallel region in the toolkit decomposes
#' its work over the elements of its output, so each accumulated value is
#' summed in full by one thread and no reduction is ever split across
#' threads; and a fold's fit is a complete, independent computation whose
#' results are collected in fold order whatever the number of processes.
#' The guarantee is also a design constraint on any kernel added later.
#'
#' At \code{threads = 1} nothing parallel is entered and no process-level
#' thread setting is touched: the code taken is the sequential path, not a
#' parallel one with a single thread. Below a per-kernel internal threshold,
#' measured where the cost of opening a parallel region overtakes its gain,
#' a kernel stays sequential whatever the count says.
#'
#' The two levels do not nest: a fit running inside a worker process is
#' sequential by construction, so \code{workers = 4} opens four processes
#' each fitting on one thread rather than \code{4 x threads} of them.
#' Worker processes are separate R sessions loading the installed packages,
#' which is what makes them safe for the S7 objects a fold's fit carries
#' (the same rule \code{optimizers7::multistart} records); starting one
#' costs on the order of a second, so \code{workers} pays on a
#' cross-validated path and is not worth asking for on a fit that takes
#' less than that.
#'
#' @param threads A single whole number, at least 1. \code{1}, the default,
#'   is sequential.
#' @param workers A single whole number, at least 1: how many R processes
#'   the independent fits of a cross-validation's folds may use. \code{1},
#'   the default, runs them in this process.
#'
#' @return An object of class \code{"n_threads"} carrying the counts.
#'
#' @examples
#' n_threads()
#' n_threads(4)
#' n_threads(1, workers = 4)
#' thread_count(n_threads(4))
#' worker_count(n_threads(1, workers = 4))
#'
#' @seealso \code{\link{thread_count}}, \code{\link{worker_count}},
#'   \code{\link{local_threads}}
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
#' @param x An \code{n_threads} object.
#' @param ... Unused.
#' @export
print.n_threads <- function(x, ...) {
  wrk <- if (x$workers != 1L) sprintf(", workers = %d", x$workers) else ""
  seq <- if (x$threads == 1L && x$workers == 1L) "  [sequential]" else ""
  cat(sprintf("n_threads(%d%s)%s\n", x$threads, wrk, seq))
  invisible(x)
}

#' Read the Count Off a Thread Policy
#'
#' @description
#' Returns the thread count an \code{\link{n_threads}} object carries, as a
#' single integer. Entry points call it once to validate their
#' \code{threads} argument and then pass the plain count down to the
#' kernels.
#'
#' @param x An object returned by \code{\link{n_threads}}.
#'
#' @return A single integer, at least 1.
#'
#' @examples
#' thread_count(n_threads())
#' thread_count(n_threads(8))
#'
#' @seealso \code{\link{n_threads}}
#' @export
thread_count <- function(x) {
  if (!inherits(x, "n_threads")) {
    stop("'threads' must be the object n_threads() returns, e.g. ",
         "threads = n_threads(4).", call. = FALSE)
  }
  x$threads
}

#' Read the Worker Count Off a Thread Policy
#'
#' @description
#' Returns the worker-process count an \code{\link{n_threads}} object
#' carries, as a single integer. An object built before \code{workers}
#' existed answers 1, so a stored policy keeps meaning what it meant.
#'
#' @param x An object returned by \code{\link{n_threads}}.
#'
#' @return A single integer, at least 1.
#'
#' @examples
#' worker_count(n_threads())
#' worker_count(n_threads(2, workers = 4))
#'
#' @seealso \code{\link{n_threads}}, \code{\link{thread_count}}
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
#' to each kernel as an argument, this call only sizes the worker pool the
#' parallel regions draw from.
#'
#' @details
#' \code{RcppParallel::setThreadOptions()} writes a process-level variable,
#' so a fit that set it and returned would leave it moved for whatever code
#' runs next; the restoration is registered with \code{on.exit} in the
#' caller's frame. At a count of 1 nothing is touched at all, which is what
#' keeps the sequential path exactly the sequential path.
#'
#' @param x An object returned by \code{\link{n_threads}}.
#' @param frame The environment whose exit restores the previous setting.
#'   Defaults to the caller's.
#'
#' @return The count, invisibly.
#'
#' @examples
#' f <- function() {
#'   local_threads(n_threads(2))
#'   Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
#' }
#' f()                                     # "2" inside the frame
#' Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "unset")   # restored
#'
#' @seealso \code{\link{n_threads}}, \code{\link{thread_count}}
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
