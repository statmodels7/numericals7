.onLoad <- function(libname, pkgname) {
  # the compiled log-Bessel kernels read the u_k polynomial table once; the
  # R twins read the same list, so the two routes share the table and
  # nothing else
  lb_set_uk_cpp(.bessel_uk)
}
