# The u_k Polynomials of the Uniform Asymptotic Expansion

Holds the coefficients of the polynomials \\u_k(t)\\ for \\k = 1, \dots,
13\\, which enter the large-order uniform asymptotic expansions of the
modified Bessel functions (DLMF 10.41). They are what
[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
and
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md)
sum in their large-order branches, and the truncation depth \\K\\ that
names each of those branches is how many of these polynomials it uses.

## Format

A list of 13 elements, one per \\k\\ in order. Each is a list of two
equal-length vectors:

- `e`:

  the exponents of \\t\\ that appear, ascending.

- `c`:

  their coefficients, alternating in sign.

## Value

The list described under Format. It is a stored constant, never called,
and is documented so that the expansion can be checked against its
source.

## Details

Each \\u_k\\ is a polynomial in \\t\\ with only every other power
present, so it is stored as the exponents that appear and their
coefficients. They satisfy the recurrence

\$\$u\_{k+1}(t) = \tfrac{1}{2} t^2 (1 - t^2)\\ u_k'(t) + \tfrac{1}{8}
\int_0^{t} (1 - 5 s^2)\\ u_k(s)\\\mathrm{d}s, \qquad u_0(t) = 1,\$\$

DLMF 10.41.9, and were generated from it in exact rational arithmetic
and converted to double only at the end. The tests re-run the recurrence
numerically and compare, so a mistyped digit fails rather than degrading
an expansion quietly.

The polynomials grow: \\u_1\\ has two terms and \\u\_{13}\\ has
fourteen, with coefficients reaching \\10^{12}\\.

## References

Olver, F. W. J., et al. (2024). *NIST Digital Library of Mathematical
Functions*, section 10.41. <https://dlmf.nist.gov/10.41>.

## See also

[`log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.md)
and
[`log_bessel_k()`](https://statmodels7.github.io/numericals7/reference/log_bessel_k.md),
whose large-order branches sum these.
