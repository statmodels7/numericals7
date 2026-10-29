# Apply a Thread Count for the Duration of One Fit

Validates a `threads` argument and sets RcppParallel's thread count to
it for the calling frame's lifetime, restoring the previous process
state when that frame exits. An entry point such as
`statmodels7::statmod()` or `distributions7::fit_distrib()` calls it
once at the start of a fit; the count itself still travels to each
kernel as an argument, and this call only sizes the worker pool from
which the parallel regions draw.

## Usage

``` r
local_threads(threads = 1, frame = parent.frame())
```

## Arguments

- threads:

  A single whole number, at least 1, and the default is `1`. The count
  is passed to the compiled kernels unchanged, and a kernel still runs
  sequentially on an input below its internal threshold. Any other value
  signals an error: a value below 1, a fractional value, `Inf`, `NA`, a
  character string and a vector of length other than one are all
  rejected with the same message.

- frame:

  The environment whose exit restores the previous setting, defaulting
  to the caller's. Pass an outer frame to hold the setting for longer
  than one call.

## Value

The thread count, invisibly: a single integer, at least 1.

## Threads and cores

`threads` is multithreading inside one R process: the compiled
per-observation kernels and the dense assembly products run on several
threads that share memory. It is distinct from `cores`, the argument of
the functions that distribute independent fits over separate R processes
(the folds of a cross-validation, the replicates of a bootstrap, the
starts of `optimizers7::multistart()`). The two counts are independent,
except that a fit running inside a separate R process is sequential by
construction, so `cores = 4` opens up to four processes and each of them
fits on one thread; it never opens `4 * threads` of them.

## Results independent of the thread count

Every parallel region in the toolkit decomposes its work over the
elements of its output. Each accumulated value is therefore summed in
full by one thread, reductions are never split across threads, and a
kernel returns the same bits at any count.

This holds because of two properties of the drivers, and not of the
kernels. A worker installs the calling thread's floating-point
environment before running its chunk, because without it some of the
platform's own math routines (among them R's `psigamma` at higher
orders, `bessel_k`, and `pgamma` and `pbeta` on the log scale) can
return results that differ in the last bits between threads. The
sequential branch also runs through the same compiled function as the
parallel one, so a compiler cannot optimize the two differently.

The dense assembly products of statmodels7 are an exception. There, a
count above 1 replaces a BLAS expression by a threaded kernel, which is
a second implementation of the same sum. Its result is bit-exact against
the reference BLAS shipped with R, and agrees within the rounding of one
dot product with an optimized BLAS such as OpenBLAS or Accelerate, which
accumulates in its own order.

The guarantee is a design constraint on any kernel added later.

## The sequential path

At `threads = 1` the code enters no parallel region and leaves the
process-level thread setting unchanged. The code taken is the sequential
one, not a parallel path running on a single thread.

Each kernel carries an internal threshold as well, measured where the
cost of opening a parallel region overtakes its gain. Below it the
kernel stays sequential regardless of the count.

## Restoration

[`RcppParallel::setThreadOptions()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)
writes a process-level variable, so a fit that set it and returned would
leave it moved for whatever runs next. The restoration is therefore
registered with `on.exit` in the caller's frame, which runs however the
fit exits, also on an error.

## Examples

``` r
f <- function() {
  local_threads(2)
  Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
}
f()                                     # "2" inside the frame
#> [1] "2"
Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "unset")   # restored
#> [1] "unset"

# A count must be a single whole number of at least one.
try(local_threads(0))
#> Error : 'threads' must be a single whole number, at least 1.
try(local_threads(2.5))
#> Error : 'threads' must be a single whole number, at least 1.
```
