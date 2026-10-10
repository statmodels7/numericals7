# Read the Worker Count Off a Policy

Returns the worker-process count an
[`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
object carries, as a single integer. A policy stored before `workers`
existed has no such component, and this function returns 1 for it, so
that an object saved with an earlier fit is read as sequential across
folds.

## Usage

``` r
worker_count(x)
```

## Arguments

- x:

  An object returned by
  [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md).
  Any other value signals an error with a message naming the
  constructor.

## Value

A single integer, at least 1. One for a policy carrying no `workers`
component.

## See also

[`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
for the object,
[`thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.md)
for the other count.

## Examples

``` r
worker_count(n_threads())
#> [1] 1
worker_count(n_threads(2, workers = 4))
#> [1] 4

# A policy stored before workers existed still reads as sequential across
# folds. Only the readers tolerate one; print() on it throws.
worker_count(structure(list(threads = 3L), class = "n_threads"))
#> [1] 1
```
