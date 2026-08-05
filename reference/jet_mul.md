# Product of Jets

Multiplies two jets, propagating every derivative exactly.

## Usage

``` r
jet_mul(a, b, lay = a$lay)
```

## Arguments

- a, b:

  Jets over the same layout.

- lay:

  A layout from
  [`jet_layout`](https://statmodels7.github.io/numericals7/reference/jet_layout.md);
  taken from the jet itself unless given.

## Value

A jet.

## Details

The Leibniz rule again, in its scalar form: the component of the product
at a tuple \\T\\ is \\\sum\_{S \subseteq T} a_S b\_{T \setminus S}\\,
the sum running over subsets of *positions* so that a repeated variable
is counted with the right multiplicity without a bookkeeping of its own.
It is the same enumeration parameters7's Leibniz machinery uses for a
matrix product.
