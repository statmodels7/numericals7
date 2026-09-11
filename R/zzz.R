.onLoad <- function(libname, pkgname) {
  # the compiled log-Bessel kernels read the u_k polynomial table once; the
  # R twins read the same list, so the two routes share the table and
  # nothing else
  lb_set_uk_cpp(.bessel_uk)
  # abs_smoother's print method sits on a base generic, and S7 registers such
  # a method for an installed package only here; pkgload registers it anyway,
  # which is how its absence would go unnoticed
  S7::methods_register()
}
