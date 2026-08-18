# Read the Worker Count Off a Thread Policy

Returns the worker-process count an
[`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
object carries, as a single integer. An object built before `workers`
existed answers 1, so a stored policy keeps meaning what it meant.

## Usage

``` r
worker_count(x)
```

## Arguments

- x:

  An object returned by
  [`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md).

## Value

A single integer, at least 1.

## See also

[`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md),
[`thread_count`](https://statmodels7.github.io/numericals7/reference/thread_count.md)

## Examples

``` r
worker_count(n_threads())
#> [1] 1
worker_count(n_threads(2, workers = 4))
#> [1] 4
```
