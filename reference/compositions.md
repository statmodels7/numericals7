# Weak Compositions of an Integer

Enumerates the weak compositions of \\n\\ into \\k\\ parts: every vector
of \\k\\ non-negative whole numbers summing to \\n\\, returned one per
row of an integer matrix. *Weak* means a part may be zero, so the set is
the whole lattice simplex of counts and includes its faces. There are
\\\binom{n + k - 1}{k - 1}\\ rows.

## Usage

``` r
compositions(n, k)
```

## Arguments

- n:

  The total, a non-negative whole number. Neither argument is validated.
  A negative `n` returns a matrix that answers no question, since the
  recursion walks `0:n` and that sequence runs downwards.

- k:

  The number of parts, a positive whole number. Zero recurses until the
  stack overflows.

## Value

An integer matrix with `k` columns and \\\binom{n + k - 1}{k - 1}\\
rows, every row summing to `n`.

## The set

\$\$\mathcal{C}(n, k) = \Bigl\\(c_1, \dots, c_k) \in \mathbb{N}\_0^{k} :
\textstyle\sum\_{j=1}^{k} c_j = n\Bigr\\, \qquad \lvert\mathcal{C}(n,
k)\rvert = \binom{n + k - 1}{k - 1}.\$\$

The count is the number of ways to place \\k - 1\\ dividers among \\n\\
units.

## Order and cost

Rows come out in lexicographic order, the first part ascending from zero
to \\n\\. The recursion runs on the number of parts, so the rows are
produced in that order directly and no full grid is built and filtered.

The whole matrix is materialized and the count grows quickly, which
bounds the practical size: \\n = 20\\ with \\k = 5\\ is 10626 rows, and
the same \\n\\ with \\k = 10\\ is 10015005.

## Where the set turns up

It is exactly the support of a multinomial with total \\n\\ over \\k\\
categories, and of any other distribution on a fixed total. An
expectation under such a law is therefore a finite sum over these rows,
evaluated exactly, where a continuous family would need a quadrature.

## See also

[`set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
and
[`tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md),
the other two enumerations collected here.

## Examples

``` r
# The four ways to split 3 into two ordered non-negative parts.
compositions(3, 2)
#>      [,1] [,2]
#> [1,]    0    3
#> [2,]    1    2
#> [3,]    2    1
#> [4,]    3    0

# Every row sums to the total, and the row count is the divider formula.
cs <- compositions(5, 3)
all(rowSums(cs) == 5)
#> [1] TRUE
nrow(cs) == choose(5 + 3 - 1, 3 - 1)
#> [1] TRUE

# Zero parts are allowed, so the faces of the simplex are included. Of the
# fifteen rows, twelve touch a face and three are strictly positive.
cs4 <- compositions(4, 3)
c(rows = nrow(cs4),
  with_a_zero = sum(apply(cs4, 1, function(r) any(r == 0))),
  positive = choose(4 - 1, 3 - 1))
#>        rows with_a_zero    positive 
#>          15          12           3 

# One part, or a total of zero, leaves a single composition.
compositions(3, 1)
#>      [,1]
#> [1,]    3
compositions(0, 3)
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
```
