# Every Way to Write an Integer as an Ordered Sum

The weak compositions of `n` into `k` parts: every vector of `k`
non-negative integers summing to `n`, one per row.

## Usage

``` r
compositions(n, k)
```

## Arguments

- n:

  The total, a non-negative integer.

- k:

  The number of parts, a positive integer.

## Value

An integer matrix with `k` columns.

## Details

The set enumerated is

\$\$\mathcal{C}(n, k) = \Bigl\\(c_1, \dots, c_k) \in \mathbb{N}\_0^{k} :
\textstyle\sum\_{j=1}^{k} c_j = n\Bigr\\, \qquad \lvert\mathcal{C}(n,
k)\rvert = \binom{n + k - 1}{k - 1},\$\$

the count being the number of ways to place \\k - 1\\ dividers among
\\n\\ units.

Built by recursion on the number of parts, which is what keeps the
result in a fixed order and avoids generating and filtering a full grid.
There are `choose(n + k - 1, k - 1)` of them, so the enumeration is only
practical for a moderate size: at `n = 20` and `k = 5` it is 10626 rows,
and at `k = 10` it is 10015005.

A discrete distribution on a fixed total – a multinomial – has exactly
this set as its support, which is what makes its expectations exact
sums.

## See also

[`set_partitions`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)

## Examples

``` r
compositions(3, 2)
#>      [,1] [,2]
#> [1,]    0    3
#> [2,]    1    2
#> [3,]    2    1
#> [4,]    3    0
nrow(compositions(5, 3)) == choose(5 + 3 - 1, 3 - 1)
#> [1] TRUE
```
