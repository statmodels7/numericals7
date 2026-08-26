# Set Partitions of the First n Integers

Enumerates every way of splitting `1:n` into disjoint non-empty blocks.
A chain rule of order \\n\\ contributes one term per partition, so this
is the index set of Faà di Bruno's formula and of the Bell-polynomial
identities built on it. The number of partitions is the Bell number
\\B_n\\: 1, 2, 5, 15, 52 and 203 for \\n\\ from one to six.

## Usage

``` r
set_partitions(n)
```

## Arguments

- n:

  A positive whole number. The recursion has no base case below one, so
  zero and negative values recurse until the stack overflows instead of
  throwing.

## Value

A list of \\B_n\\ partitions. Each partition is a list of integer
vectors, its blocks, which between them contain each of `1:n` exactly
once. The storage is integer whether `n` is given as `4` or as `4L`, as
it is for
[`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
and
[`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md),
so the three enumerations agree and a block may be compared against
`1:n` with [`identical()`](https://rdrr.io/r/base/identical.html).

## How they are built

By the standard recursion. The partitions of `1:n` come from those of
`1:(n-1)` by putting \\n\\ into each existing block in turn, and then
into a block of its own.

The cost is therefore \\B_n\\, which grows faster than any exponential:
\\B_8\\ is 4140 and \\B\_{10}\\ is 115975. Four is as high as the
toolkit's derivatives go, where the sum has fifteen terms.

## The blocks index positions

A block holds positions within a multi-index, not variables. A variable
repeated in that multi-index therefore counts with its correct
multiplicity and needs no bookkeeping of its own: differentiating three
times in one variable and once in another sums over the same fifteen
partitions as four distinct variables do, and what differs is which
derivative each block names.
[`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
supplies the multi-index the positions point into.

## References

Constantine, G. M. and Savits, T. H. (1996). A multivariate Faà di Bruno
formula with applications. *Transactions of the American Mathematical
Society* **348**, 503-520.

## See also

[`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md)
for the multi-indices the blocks index into, and
[`compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.md)
for the ordered sums.

## Examples

``` r
# The five partitions of {1, 2, 3}.
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

# The counts are the Bell numbers.
lengths(lapply(1:6, set_partitions))
#> [1]   1   2   5  15  52 203

# Every partition covers 1:n exactly once, the blocks being disjoint. The
# comparison is `identical()` because the blocks are integer however `n` was
# given, which is the convention all three enumerations follow.
all(vapply(set_partitions(4), function(p) identical(sort(unlist(p)), 1:4),
           logical(1)))
#> [1] TRUE

# A fourth-order chain rule sums fifteen terms, one per partition, and the
# number of blocks in each is the number of factors in that term.
table(lengths(set_partitions(4)))
#> 
#> 1 2 3 4 
#> 1 7 6 1 
```
