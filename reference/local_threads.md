# Apply a Thread Policy for the Duration of One Fit

Sets RcppParallel's thread count to the policy's for the calling frame's
lifetime and restores the previous process state when that frame exits.
Called once at the entry of a fit; the count itself still travels to
each kernel as an argument, this call only sizes the worker pool the
parallel regions draw from.

## Usage

``` r
local_threads(x, frame = parent.frame())
```

## Arguments

- x:

  An object returned by
  [`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md).

- frame:

  The environment whose exit restores the previous setting. Defaults to
  the caller's.

## Value

The count, invisibly.

## Details

[`RcppParallel::setThreadOptions()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)
writes a process-level variable, so a fit that set it and returned would
leave it moved for whatever code runs next; the restoration is
registered with `on.exit` in the caller's frame. At a count of 1 nothing
is touched at all, which is what keeps the sequential path exactly the
sequential path.

## See also

[`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md),
[`thread_count`](https://statmodels7.github.io/numericals7/reference/thread_count.md)

## Examples

``` r
f <- function() {
  local_threads(n_threads(2))
  Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
}
f()                                     # "2" inside the frame
#> [1] "2"
Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "unset")   # restored
#> [1] "unset"
```
