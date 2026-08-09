## Test environments

* local Windows 11, R 4.6.0
* GitHub Actions:
  * macOS latest, R release
  * Windows latest, R release
  * Ubuntu latest, R devel
  * Ubuntu latest, R release
  * Ubuntu latest, R oldrel-1

## R CMD check results

0 errors | 0 warnings | 1 note

The note is

```
Maintainer: 'Giovanni Tinervia <giovannitinervia9@gmail.com>'
New submission
```

which is expected for a first submission.

## Notes for the reviewer

The package compiles one C++ translation unit through Rcpp. Its kernels are
paired with pure R implementations of the same recursions, kept in the test
suite and compared against the compiled route at machine precision, so a
platform on which the two disagree fails the tests rather than returning a
different answer quietly.

`log_bessel_i()` and `log_bessel_k()` implement the uniform asymptotic
expansions of Plesner, Sorensen and Hauberg (arXiv:2409.08729), cited on
their help pages.

## Downstream dependencies

None on CRAN. This is the first package of a toolkit whose other members
depend on it and are not yet submitted.
