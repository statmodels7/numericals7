# How Many Threads and Worker Processes a Fit May Use

Constructs the parallelism policy an entry point such as
`statmodels7::statmod()` or `distributions7::fit_distrib()` accepts
through its `threads` argument. The default, one thread and one process,
takes exactly the sequential code path; a larger `threads` lets the
compiled per-observation kernels and the dense assembly products run in
parallel, and a larger `workers` fans the independent fits of a
cross-validation's folds out over separate R processes.

## Usage

``` r
n_threads(threads = 1, workers = 1)

# S3 method for class 'n_threads'
print(x, ...)
```

## Arguments

- threads:

  A single whole number, at least 1. `1`, the default, is sequential.

- workers:

  A single whole number, at least 1: how many R processes the
  independent fits of a cross-validation's folds may use. `1`, the
  default, runs them in this process.

- x:

  An `n_threads` object.

- ...:

  Unused.

## Value

An object of class `"n_threads"` carrying the counts.

## Details

The count never changes what a kernel computes: every parallel region in
the toolkit decomposes its work over the elements of its output, so each
accumulated value is summed in full by one thread and no reduction is
ever split across threads – a kernel's result is bit for bit the same at
any count – and a fold's fit is a complete, independent computation
whose results are collected in fold order whatever the number of
processes. One qualification: where raising the count first ENGAGES a
threaded kernel in place of a BLAS expression (the dense assembly
products of statmodels7), the replacement is a second implementation of
the same sum, bit-exact against the reference BLAS R ships and within
the rounding of one dot product against an optimized one (OpenBLAS,
Accelerate), whose accumulation order is its own. The guarantee is also
a design constraint on any kernel added later.

At `threads = 1` nothing parallel is entered and no process-level thread
setting is touched: the code taken is the sequential path, not a
parallel one with a single thread. Below a per-kernel internal
threshold, measured where the cost of opening a parallel region
overtakes its gain, a kernel stays sequential whatever the count says.

The two levels do not nest: a fit running inside a worker process is
sequential by construction, so `workers = 4` opens four processes each
fitting on one thread rather than `4 x threads` of them. Worker
processes are separate R sessions loading the installed packages, which
is what makes them safe for the S7 objects a fold's fit carries (the
same rule `optimizers7::multistart` records); starting one costs on the
order of a second, so `workers` pays on a cross-validated path and is
not worth asking for on a fit that takes less than that.

## See also

[`thread_count`](https://statmodels7.github.io/numericals7/reference/thread_count.md),
[`worker_count`](https://statmodels7.github.io/numericals7/reference/worker_count.md),
[`local_threads`](https://statmodels7.github.io/numericals7/reference/local_threads.md)

## Examples

``` r
n_threads()
#> n_threads(1)  [sequential]
n_threads(4)
#> n_threads(4)
n_threads(1, workers = 4)
#> n_threads(1, workers = 4)
thread_count(n_threads(4))
#> [1] 4
worker_count(n_threads(1, workers = 4))
#> [1] 4
```
