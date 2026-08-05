# A Jet for One Variable

The jet of a scalar map of a single free value: its value and its four
derivatives sit in the pure components of that variable, everything else
being zero.

## Usage

``` r
jet_var(k, dv, lay)
```

## Arguments

- k:

  The index of the free value.

- dv:

  A list of five numbers: the value and four derivatives.

- lay:

  A layout from
  [`jet_layout`](https://statmodels7.github.io/numericals7/reference/jet_layout.md);
  taken from the jet itself unless given.

## Value

A jet.

## See also

[`jet_layout`](https://statmodels7.github.io/numericals7/reference/jet_layout.md),
[`Ops.jet`](https://statmodels7.github.io/numericals7/reference/Ops.jet.md)

## Examples

``` r
# seed two variables, then write a map in ordinary R
lay <- jet_layout(2)
x <- jet_var(1, list(1.3, 1, 0, 0, 0), lay)
y <- jet_var(2, list(0.8, 1, 0, 0, 0), lay)
z <- x / gamma(1 + 1 / y)
z$v
#> [1] 1.147393
z$d[[1]]
#> [1] 0.8826101 1.0264623
```
