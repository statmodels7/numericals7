# Read the Count Off a Thread Policy

Returns the thread count an
[`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
object carries, as a single integer. Entry points call it once to
validate their `threads` argument and then pass the plain count down to
the kernels.

## Usage

``` r
thread_count(x)
```

## Arguments

- x:

  An object returned by
  [`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md).

## Value

A single integer, at least 1.

## See also

[`n_threads`](https://statmodels7.github.io/numericals7/reference/n_threads.md)

## Examples

``` r
thread_count(n_threads())
#> [1] 1
thread_count(n_threads(8))
#> [1] 8
```
