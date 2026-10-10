# Multi-Indices for Derivatives of a Given Order

Enumerates the multi-indices that key a list of partial derivatives:
every non-decreasing tuple \\(i_1 \le \cdots \le i_k)\\ drawn from
`1:d`, where \\k\\ is the derivative order and \\d\\ the number of
variables. There is one tuple per distinct partial derivative, because a
mixed partial does not depend on the order in which the variables are
differentiated, so the count is \\\binom{d + k - 1}{k}\\.

## Usage

``` r
tuple_indices(d, order = 2L)
```

## Arguments

- d:

  The number of variables, a non-negative whole number. Zero gives an
  empty list. A negative value signals an error: at orders one and two
  it comes from [`seq_len()`](https://rdrr.io/r/base/seq.html) and
  concerns coercion to a non-negative integer, and at orders three and
  four it comes from
  [`utils::combn()`](https://rdrr.io/r/utils/combn.html).

- order:

  The derivative order \\k\\, one of `1`, `2`, `3` or `4`. Any other
  value signals an error, because the enumeration is written out only to
  order four.

## Value

A list of \\\binom{d + k - 1}{k}\\ integer vectors, each of length
`order`, each non-decreasing. An empty list when `d` is zero.

## The ordering

At order two the diagonal comes first, \\(1,1), (2,2), \dots, (d,d)\\,
and the off-diagonal pairs follow in lexicographic order, which is the
order used to index a Hessian. At orders three and four the enumeration
is plain lexicographic over non-decreasing tuples.

Every object in the toolkit holding derivatives over \\d\\ variables
keys its components by this enumeration, so two such objects can be
combined without reordering either. The order is fixed.

## Counts

\$\$\lvert T(d, k)\rvert = \binom{d + k - 1}{k},\$\$

the number of multisets of size \\k\\ drawn from \\d\\ symbols. At \\d =
3\\ the four orders give 3, 6, 10 and 15 tuples.

## See also

[`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md),
the other enumeration that a higher-order chain rule needs, and
[`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md)
for the ordered sums.

## Examples

``` r
# Second-order tuples of three variables: the diagonal first, then the
# off-diagonal pairs. A Hessian consumer indexes in this order.
tuple_indices(3, 2)
#> [[1]]
#> [1] 1 1
#> 
#> [[2]]
#> [1] 2 2
#> 
#> [[3]]
#> [1] 3 3
#> 
#> [[4]]
#> [1] 1 2
#> 
#> [[5]]
#> [1] 1 3
#> 
#> [[6]]
#> [1] 2 3
#> 

# One tuple per distinct partial derivative, so the count is the number of
# multisets of size k drawn from d symbols.
vapply(1:4, function(k) length(tuple_indices(3, k)), integer(1))
#> [1]  3  6 10 15
choose(3 + 1:4 - 1, 1:4)
#> [1]  3  6 10 15

# Every tuple is non-decreasing, so no derivative is enumerated twice.
all(vapply(tuple_indices(4, 3), function(i) !is.unsorted(i), logical(1)))
#> [1] TRUE

# With no variables there are no derivatives.
tuple_indices(0, 2)
#> list()
```
