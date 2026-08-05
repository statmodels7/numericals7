# The Set Partitions of the First n Integers

Every way of splitting `1:n` into disjoint non-empty blocks, which is
what a chain rule of order \\n\\ sums over.

## Usage

``` r
set_partitions(n)
```

## Arguments

- n:

  A positive integer, at most four here.

## Value

A list of partitions, each a list of integer vectors.

## Details

Built by the standard recursion: the partitions of `1:n` are obtained
from those of `1:(n-1)` by placing `n` into each existing block in turn
and then into a block of its own. There are 1, 2, 5 and 15 of them for
\\n = 1, \dots, 4\\, the Bell numbers.

The blocks index **positions** rather than variables, which is what
makes a repeated variable count with the right multiplicity without any
bookkeeping of its own – the same device
[`jet_mul`](https://statmodels7.github.io/numericals7/reference/jet_mul.md)
uses with subsets of positions.

## See also

[`jet_compose`](https://statmodels7.github.io/numericals7/reference/jet_compose.md),
[`tuple_indices`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)

## Examples

``` r
set_partitions(3)
#> [[1]]
#> [[1]][[1]]
#> [1] 1 2 3
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 1 2
#> 
#> [[2]][[2]]
#> [1] 3
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] 1 3
#> 
#> [[3]][[2]]
#> [1] 2
#> 
#> 
#> [[4]]
#> [[4]][[1]]
#> [1] 1
#> 
#> [[4]][[2]]
#> [1] 2 3
#> 
#> 
#> [[5]]
#> [[5]][[1]]
#> [1] 1
#> 
#> [[5]][[2]]
#> [1] 2
#> 
#> [[5]][[3]]
#> [1] 3
#> 
#> 
lengths(lapply(1:4, set_partitions))
#> [1]  1  2  5 15
```
