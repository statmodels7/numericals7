# The Index Tuples of a Given Width

Every multi-index of a derivative of the given order over `d` variables:
the tuples \\(i_1 \le \dots \le i_k)\\ that key a list of partial
derivatives.

## Usage

``` r
tuple_indices(d, order = 2L)
```

## Arguments

- d:

  The number of variables.

- order:

  The derivative order, 1 to 4.

## Value

A list of integer vectors of length `order`.

## Details

The ordering is part of the contract. At order two the diagonal comes
first and the off-diagonal pairs follow, because that is the order a
Hessian consumer indexes by; at orders three and four the enumeration is
the lexicographic one over non-decreasing tuples. Anything that holds
derivatives over \\d\\ variables – a jet, a constrained parameter – keys
its components by this enumeration, so nothing has to be reordered when
they meet.

## See also

[`jet_layout`](https://statmodels7.github.io/numericals7/reference/jet_layout.md),
[`set_partitions`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)

## Examples

``` r
tuple_indices(2, 2)
#> [[1]]
#> [1] 1 1
#> 
#> [[2]]
#> [1] 2 2
#> 
#> [[3]]
#> [1] 1 2
#> 
lengths(tuple_indices(3, 4))[1:3]
#> [1] 4 4 4
# there are choose(d + k - 1, k) tuples of width k
length(tuple_indices(3, 4)) == choose(3 + 4 - 1, 4)
#> [1] TRUE
```
