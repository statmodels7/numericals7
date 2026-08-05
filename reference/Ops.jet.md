# Arithmetic on Jets

The arithmetic operators applied to jets, so that a map is written in
ordinary R and its derivatives come out exact.

## Usage

``` r
# S3 method for class 'jet'
Ops(e1, e2)
```

## Arguments

- e1, e2:

  Jets, or a jet and a number.

## Value

A jet.

## Details

This is what makes a jet worth having. A quantity written as
`mu / gamma(1 + 1 / sigma)` carries every partial derivative to fourth
order with it, and whoever wrote the expression need not know that: the
operators dispatch, and nothing about the derivative machinery appears
in the expression.

A number mixed with a jet is treated as a constant. Multiplication and
division by one take the cheap route rather than building a constant jet
and running the full Leibniz rule. A jet raised to a numeric power goes
through
[`jet_pow`](https://statmodels7.github.io/numericals7/reference/jet_pow.md);
a jet raised to a jet through \\\exp(b\log a)\\.

Comparison and logical operators are refused. A jet is a value with its
derivatives attached, and a branch taken on it would silently discard
them: whichever side the comparison chose, the result would carry the
derivatives of that side alone and claim to be the derivative of the
whole expression.

## See also

[`Math.jet`](https://statmodels7.github.io/numericals7/reference/Math.jet.md),
[`jet_compose`](https://statmodels7.github.io/numericals7/reference/jet_compose.md)
