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

  A single whole number, at least 1, and the default is `1`. The count
  is what the kernels are given and what they run on, not a ceiling they
  may come under. Anything else throws: a value below 1, a fractional
  value, `Inf`, a character string and a vector of length other than one
  are all rejected with the same message.

- workers:

  A single whole number, at least 1, and the default is `1`: how many R
  processes the independent fits of a cross-validation's folds may use.
  `1` runs them in this process. Validated exactly as `threads` is.

- x:

  An `n_threads` object.

- ...:

  Unused, and present because
  [`print()`](https://rdrr.io/r/base/print.html) requires it.

## Value

`n_threads()` returns an object of class `"n_threads"`, a list of two
components:

- `threads`:

  integer, the thread count, at least 1.

- `workers`:

  integer, the worker-process count, at least 1.

Read them with
[`thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.md)
and
[`worker_count()`](https://statmodels7.github.io/numericals7/reference/worker_count.md),
which check the class first.
[`print()`](https://rdrr.io/r/base/print.html) returns its argument
invisibly.

## The count does not change the answer

Every parallel region in the toolkit decomposes its work over the
elements of its output. Each accumulated value is therefore summed in
full by one thread, no reduction is ever split across threads, and a
kernel returns the same bits at any count. A fold's fit is a complete
independent computation, and the folds are collected in fold order
however many processes ran them.

Two properties of the drivers make that hold, and neither belongs to the
kernels themselves. A worker installs the calling thread's
floating-point environment before running its chunk: without it some of
the platform's own math routines return per-thread last bits, and R's
`psigamma` at higher orders, `bessel_k`, and `pgamma` and `pbeta` on the
log scale were each measured doing so. The sequential branch also runs
through the same compiled function the parallel one does, so a compiler
cannot optimize the two apart.

One qualification. In the dense assembly products of statmodels7,
raising the count first *engages* a threaded kernel where a BLAS
expression stood. That kernel is a second implementation of the same
sum: bit-exact against the reference BLAS R ships, and within the
rounding of one dot product against an optimized BLAS such as OpenBLAS
or Accelerate, whose accumulation order is its own.

The guarantee is a design constraint on any kernel added later.

## The sequential path is the sequential path

At `threads = 1` nothing parallel is entered and no process-level thread
setting is touched. The code taken is the sequential one, not a parallel
path running on a single thread.

Each kernel carries an internal threshold as well, measured where the
cost of opening a parallel region overtakes its gain. Below it the
kernel stays sequential whatever the count says.

## The two levels do not nest

A fit running inside a worker process is sequential by construction, so
`workers = 4` opens four processes each fitting on one thread, never
`4 * threads` of them.

Worker processes are separate R sessions that load the installed
packages, which makes them safe for the S7 objects a fold's fit carries;
this is the rule `optimizers7::multistart` records. Starting one costs
on the order of a second, so `workers` pays on a cross-validated path
and is not worth asking for on a fit shorter than that.

## Printing

[`print()`](https://rdrr.io/r/base/print.html) shows the call that would
rebuild the object, with `workers` shown only when it is not 1, and
appends `[sequential]` when both counts are 1. It returns its argument
invisibly.

## See also

[`thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.md)
and
[`worker_count()`](https://statmodels7.github.io/numericals7/reference/worker_count.md)
to read the counts,
[`local_threads()`](https://statmodels7.github.io/numericals7/reference/local_threads.md)
to apply one for the duration of a fit.

## Examples

``` r
# The default is sequential, and says so.
n_threads()
#> n_threads(1)  [sequential]

# Threads within a fit, worker processes across folds, or both.
n_threads(4)
#> n_threads(4)
n_threads(1, workers = 4)
#> n_threads(1, workers = 4)
n_threads(4, workers = 2)
#> n_threads(4, workers = 2)

# The counts come back out through the two readers.
p <- n_threads(4, workers = 2)
c(threads = thread_count(p), workers = worker_count(p))
#> threads workers 
#>       4       2 

# A count must be a single whole number of at least one.
try(n_threads(0))
#> Error : 'threads' must be a single whole number, at least 1.
try(n_threads(2.5))
#> Error : 'threads' must be a single whole number, at least 1.
```
