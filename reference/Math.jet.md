# Mathematical Functions of a Jet

The mathematical functions a jet understands, so that a map written in
ordinary R differentiates itself.

## Usage

``` r
# S3 method for class 'jet'
Math(x, ...)
```

## Arguments

- x:

  A jet.

- ...:

  Further arguments; refused, since a second argument would change what
  is being differentiated.

## Value

A jet.

## Details

Each is one call to
[`jet_compose`](https://statmodels7.github.io/numericals7/reference/jet_compose.md)
with the function's own five derivatives, so extending the vocabulary
means writing five derivatives and not a chain rule.

A function that is not smooth is refused rather than approximated:
`abs`, `floor`, `ceiling`, `round` and `sign` have no fourth derivative
to carry, and returning one would be a wrong answer in the shape of a
right one.

## See also

[`Ops.jet`](https://statmodels7.github.io/numericals7/reference/Ops.jet.md),
[`jet_compose`](https://statmodels7.github.io/numericals7/reference/jet_compose.md)
