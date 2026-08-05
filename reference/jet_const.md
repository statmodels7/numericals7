# A Constant Jet

A number with every derivative zero.

## Usage

``` r
jet_const(v, lay)
```

## Arguments

- v:

  The value.

- lay:

  A layout from
  [`jet_layout`](https://statmodels7.github.io/numericals7/reference/jet_layout.md);
  taken from the jet itself unless given.

## Value

A jet.

## See also

[`jet_var`](https://statmodels7.github.io/numericals7/reference/jet_var.md),
[`jet_layout`](https://statmodels7.github.io/numericals7/reference/jet_layout.md)

## Examples

``` r
lay <- jet_layout(2)
x <- jet_var(1, list(2, 1, 0, 0, 0), lay)
acc <- jet_const(0, lay)   # a running sum starts at a constant
(acc + x * 3)$v
#> [1] 6
```
