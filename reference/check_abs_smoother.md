# Check a Smoother of the Absolute Value Numerically

Compares a smoother against the properties every smoother of \\\lvert
u\rvert\\ must have, and each derivative order against one numerical
differentiation of the analytic order below it. Returns one row per
check with the worst error and a pass or fail. Write a smoother of your
own and this says whether its five derivatives are right.

## Usage

``` r
check_abs_smoother(smoother, width = NULL, tol = 1e-06, verbose = TRUE)
```

## Arguments

- smoother:

  An
  [`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md).
  Anything else is rejected.

- width:

  The width to check at, a single positive number, or `NULL` (the
  default) for the smoother's own. A smoother carrying none is checked
  at `0.5`, a scale chosen to exercise the formulas and saying nothing
  about any data.

- tol:

  The relative error above which a row is reported as `FAILED`. A single
  positive number, `1e-6` by default. The worst error the three shipped
  smoothers reach is `1.3e-10`, four orders under it, so the default
  separates a correct smoother from one whose derivative is wrong in its
  fourth digit.

- verbose:

  `TRUE`, the default, prints a header naming the smoother and the
  width, then the table without row names, and returns invisibly.
  `FALSE` returns the table.

## Value

A data frame with one row per check and three columns: `check`
(character), `max_error` (numeric) and `status` (character, `"OK"` or
`"FAILED"`). Ten rows, or twelve when the smoother declares a
`tau_correction`. Returned invisibly when `verbose` is `TRUE`.

## The rows

Ten checks, and two more for a smoother declaring a scale correction:

|  |  |
|----|----|
| row | what it asserts |
| `s is even` | \\s(-u) = s(u)\\, relative to the size of \\s\\ |
| `s' is odd` | \\s'(-u) = -s'(u)\\ |
| `\|s'\| bounded by one` | the smooth sign never leaves \\\[-1, 1\]\\ |
| `s convex` | \\s'' \ge 0\\ on the grid |
| `matches \|u\| in the tails` | the excess ten intrinsic widths out is at most a fifth of the excess at the kink |
| `order k vs numDeriv on order k-1` | five rows, \\k\\ from 1 to 5 |
| `tau_correction bounded by tau` | the correction never increases a scale |
| `tau_correction is the identity at width zero` | it vanishes as the smoothing does |

Order \\k\\ goes against numDeriv applied **once** to the analytic order
\\k-1\\, so a wrong derivative is caught against an independent route
while the reference never degenerates into a difference of differences.
Over the three shipped smoothers at four widths each, every row passes
and the worst error is `1.3e-10`.

The tail row is deliberately loose, at a fifth of the excess at the
kink, because
[`smooth_hyperbolic()`](https://statmodels7.github.io/numericals7/reference/smooth_hyperbolic.md)'s
polynomial tails would fail a tight one and are a legitimate choice: its
excess is \\c/(2\lvert u\rvert)\\, a twentieth of the kink's value ten
widths out.

## Where the grid goes

The grid spans the transition and the tails, placed by the smoother's
own intrinsic scale \\s(0)\\, which is of order the transition width for
any smoother of the absolute value, so no reading of what the width
parameter means is needed. Every point sits off zero and off the seam: a
piecewise smoother has measure-zero points where a one-sided derivative
is read, and a difference straddling one compares nothing.

## What it catches

This is not decoration. Writing
[`smooth_quintic()`](https://statmodels7.github.io/numericals7/reference/smooth_quintic.md),
a factor-2 error in its own \\s'''\\ was caught here before anything
shipped: the order-3 row read a relative `0.5` and the order-4 row
`1.0`, the wrong third derivative becoming the reference for the fourth.
The examples below reproduce it.

## Errors

numDeriv is in `Suggests` and every derivative row needs it, so the
function stops when it is not installed.

## See also

[`abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.md)
for the contract being checked,
[`smoother_deriv()`](https://statmodels7.github.io/numericals7/reference/smoother_deriv.md)
for the derivatives it reads, `penalties7::check_penalty()` for the
sibling on a penalty.

## Examples

``` r
res <- check_abs_smoother(smooth_probit(h = 0.3))
#> check_abs_smoother: probit at h = 0.3
#>                                         check    max_error status
#>                                     s is even 6.023629e-17     OK
#>                                     s' is odd 1.110223e-16     OK
#>                           |s'| bounded by one 0.000000e+00     OK
#>                                      s convex 0.000000e+00     OK
#>                      matches |u| in the tails 0.000000e+00     OK
#>                order 1 vs numDeriv on order 0 5.906771e-11     OK
#>                order 2 vs numDeriv on order 1 6.359873e-11     OK
#>                order 3 vs numDeriv on order 2 1.267688e-10     OK
#>                order 4 vs numDeriv on order 3 1.657041e-11     OK
#>                order 5 vs numDeriv on order 4 4.669206e-11     OK
#>                 tau_correction bounded by tau 0.000000e+00     OK
#>  tau_correction is the identity at width zero 0.000000e+00     OK
all(res$status == "OK")
#> [1] TRUE

# Ten rows for a smoother with no scale correction, twelve with one.
nrow(check_abs_smoother(smooth_quintic(h = 0.3), verbose = FALSE))
#> [1] 10
nrow(check_abs_smoother(smooth_probit(h = 0.3), verbose = FALSE))
#> [1] 12

# The error this caught while the quintic was being written: its own
# third derivative with a 4 in the denominator where the algebra gives 2.
broken <- smooth_quintic(h = 0.3)
fns <- broken@s
fns[[4]] <- function(u, width) {
  t <- u / width
  ifelse(abs(u) < width, -15 * t * (1 - t^2) / (4 * width^2), 0)
}
broken@s <- fns
bad <- check_abs_smoother(broken, verbose = FALSE)
bad[bad$status != "OK", c("check", "max_error")]
#>                            check max_error
#> 8 order 3 vs numDeriv on order 2       0.5
#> 9 order 4 vs numDeriv on order 3       1.0
```
