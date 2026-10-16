# Helpers of the Tensor Stencil

These are internal functions of
[`fd_tensor()`](https://statmodels7.github.io/numericals7/reference/fd_tensor.md).
`tensor_plan()` returns the plan of a request from the session's cache,
computing it in C++ on the first call with that key. `is_whole()` tests
for a single finite whole number. `counts_of_order()` returns the
multi-indices of one total order as rows of counts, in the order of
[`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md).
`as_counts()` validates a matrix of multi-indices and stores it as
integers. `as_columns()` turns the value of a source function into a
list of double vectors, one per declared multi-index.

## Usage

``` r
tensor_plan(key, comp, kidx, accuracy)

is_whole(x)

counts_of_order(p, k)

as_counts(m, p, what)

as_columns(v, k, src)
```

## Arguments

- key:

  A string identifying the request: the accuracy, the number of
  coordinates, the components and the declarations.

- comp, m:

  An integer matrix of counts, one row per multi-index.

- kidx:

  A list of such matrices, of which the first is the function's own zero
  row.

- accuracy:

  The order of accuracy of the stencils.

- x:

  A value to test.

- p, k:

  The number of coordinates and a total order; in `as_columns()`, `k` is
  the number of declared multi-indices.

- what:

  The argument name used in an error message.

- v:

  The value returned by a source function.

- src:

  The position of the source: 1 for `f`, \\j + 1\\ for `known[[j]]`.

## Value

`tensor_plan()` a list of the plan's vectors and matrices; `is_whole()`
a single logical; `counts_of_order()` and `as_counts()` an integer
matrix; `as_columns()` a list of double vectors.
