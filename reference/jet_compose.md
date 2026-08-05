# A Smooth Function of a Jet

Applies a scalar function to a jet, given that function's own
derivatives at the jet's value, and propagates every partial derivative
exactly.

## Usage

``` r
jet_compose(a, fd, lay = a$lay)
```

## Arguments

- a:

  A jet.

- fd:

  A numeric vector of five values: \\f\\ and its first four derivatives,
  all evaluated at `a$v`.

- lay:

  A layout from
  [`jet_layout`](https://statmodels7.github.io/numericals7/reference/jet_layout.md);
  taken from the jet itself unless given.

## Value

A jet.

## Details

This is Faa di Bruno in the form the multivariate case takes when the
inner function is scalar valued: for a tuple \\T\\ of positions, \$\$(f
\circ a)\_T = \sum\_{\pi} f^{(\|\pi\|)}(a)\prod\_{B \in \pi}
a\_{T\[B\]}\$\$ the sum running over the set partitions of the positions
of \\T\\.

Every transcendental a jet needs is one call to this with the right five
numbers, so the exponential, the logarithm, a power and the gamma
function need no derivative machinery of their own. Extending the
vocabulary is a matter of writing five derivatives, not of writing a
chain rule.

## See also

[`jet_mul`](https://statmodels7.github.io/numericals7/reference/jet_mul.md),
[`set_partitions`](https://statmodels7.github.io/numericals7/reference/set_partitions.md)
