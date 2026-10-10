#' The Modified Bessel Function of the Third Kind
#'
#' @description
#' Computes \eqn{K_\nu(x)}, or \eqn{e^x K_\nu(x)} when `expon.scaled = TRUE`,
#' for real orders \eqn{\nu}. The values are those of [base::besselK()], from a
#' compiled function that never reaches the R API, so that it may be called
#' from a worker thread. Its C entry point is `n7_bessel_k`, resolved with
#' `R_GetCCallable("numericals7", "n7_bessel_k")` and taking the arguments
#' `(x, alpha, expo)` that R's C function `bessel_k()` takes, `expo` being 1
#' for the function and 2 for the scaled function.
#'
#' @details
#' The computation is R's own: the routine RKBESL of W. J. Cody and L. Stoltz,
#' after J. B. Campbell's implementation of Temme's algorithm, as R carries it,
#' copied with the warnings R raises from it removed and with its work array
#' on the stack (allocated with `malloc()` beyond 64 orders) instead of in R's
#' memory pool. A warning calls into the R API,
#' and R's allocator is not safe to call from a worker thread.
#'
#' @param x A numeric vector of non-negative arguments.
#' @param nu A numeric vector of orders, recycled against `x`; the order
#'   enters as \eqn{|\nu|}.
#' @param expon.scaled If `TRUE`, \eqn{e^x K_\nu(x)}; otherwise, the default,
#'   \eqn{K_\nu(x)}.
#' @param threads A single positive integer, how many threads the
#'   computation may use. Defaults to `1L`.
#'
#' @return A numeric vector of the recycled length of `x` and `nu`. A negative
#'   `x` gives `NaN`.
#'
#' @references
#' Campbell, J. B. (1980). On Temme's algorithm for the modified Bessel
#' function of the third kind. *ACM Transactions on Mathematical Software* 6,
#' 581-586.
#'
#' @seealso [base::besselK()], whose values these are, and [log_bessel_k()]
#'   for the logarithm where the function itself overflows or underflows.
#'
#' @examples
#' x <- c(0.01, 0.5, 3, 40)
#' identical(bessel_k(x, nu = 1.5), besselK(x, nu = 1.5))
#'
#' # the scaled function at a large argument, where K itself underflows
#' bessel_k(1000, nu = 2, expon.scaled = TRUE)
#'
#' @export
bessel_k <- function(x, nu, expon.scaled = FALSE, threads = 1L) {
  bessel_k_cpp(as.double(x), as.double(nu), isTRUE(expon.scaled),
               as.integer(threads))
}
