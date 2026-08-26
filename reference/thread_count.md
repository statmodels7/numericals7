# Read the Thread Count Off a Policy

Returns the thread count an
[`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
object carries, as a single integer. An entry point calls it once to
validate the `threads` argument it was given, then passes the plain
count down to the kernels, which take a number and know nothing about
the policy object.

## Usage

``` r
thread_count(x)
```

## Arguments

- x:

  An object returned by
  [`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md).
  Anything else throws, with a message naming the constructor, because a
  bare number reaching an entry point is the likely mistake and would
  otherwise be read as a policy.

## Value

A single integer, at least 1.

## See also

[`n_threads()`](https://statmodels7.github.io/numericals7/reference/n_threads.md)
for the object,
[`worker_count()`](https://statmodels7.github.io/numericals7/reference/worker_count.md)
for the other count.

## Examples

``` r
thread_count(n_threads())
#> [1] 1
thread_count(n_threads(8))
#> [1] 8

# A bare count is not a policy.
try(thread_count(8))
#> Error : 'threads' must be the object n_threads() returns, e.g. threads = n_threads(4).
```
