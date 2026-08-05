# The Bookkeeping a Jet Needs

The index tuples up to fourth order over \\d\\ variables, together with
the lookup from a tuple to its position, computed once and shared by
every jet in a calculation.

## Usage

``` r
jet_layout(d)
```

## Arguments

- d:

  The number of variables.

## Value

A list with `d`, the tuple lists `tuples`, and the environment `pos`
mapping a sorted tuple to its order and position.

## Details

The tuples come from
[`tuple_indices`](https://statmodels7.github.io/numericals7/reference/tuple_indices.md),
the one enumeration the toolkit shares, so a jet's components are keyed
exactly as any other holder of derivatives over the same variables and
nothing has to be reordered on the way out.

## See also

[`jet_var`](https://statmodels7.github.io/numericals7/reference/jet_var.md),
[`Ops.jet`](https://statmodels7.github.io/numericals7/reference/Ops.jet.md)

## Examples

``` r
lay <- jet_layout(2)
lay$d
#> [1] 2
length(lay$tuples[[2]])
#> [1] 3
```
