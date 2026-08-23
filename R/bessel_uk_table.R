# u_k(t) polynomials of the uniform asymptotic expansion, k = 1..13,
# generated in exact rational arithmetic from the recurrence
# u_{k+1} = t^2(1-t^2)/2 u_k' + (1/8) int_0^t (1-5s^2) u_k ds
# (DLMF 10.41.4); coefficients converted to double at the end.
#' The u_k Polynomials of the Uniform Asymptotic Expansion
#'
#' @description
#' Holds the coefficients of the polynomials \eqn{u_k(t)} for
#' \eqn{k = 1, \dots, 13}, which enter the large-order uniform asymptotic
#' expansions of the modified Bessel functions (DLMF 10.41). They are what
#' [log_bessel_i()] and [log_bessel_k()] sum in their large-order branches, and
#' the truncation depth \eqn{K} that names each of those branches is how many of
#' these polynomials it uses.
#'
#' @details
#' Each \eqn{u_k} is a polynomial in \eqn{t} with only every other power
#' present, so it is stored as the exponents that appear and their
#' coefficients. They satisfy the recurrence
#'
#' \deqn{u_{k+1}(t) = \tfrac{1}{2} t^2 (1 - t^2)\, u_k'(t)
#'       + \tfrac{1}{8} \int_0^{t} (1 - 5 s^2)\, u_k(s)\,\mathrm{d}s,
#'       \qquad u_0(t) = 1,}
#'
#' DLMF 10.41.9, and were generated from it in exact rational arithmetic and
#' converted to double only at the end. The tests re-run the recurrence
#' numerically and compare, so a mistyped digit fails rather than degrading an
#' expansion quietly.
#'
#' The polynomials grow: \eqn{u_1} has two terms and \eqn{u_{13}} has fourteen,
#' with coefficients reaching \eqn{10^{12}}.
#'
#' @format A list of 13 elements, one per \eqn{k} in order. Each is a list of
#'   two equal-length vectors:
#'   \describe{
#'     \item{`e`}{the exponents of \eqn{t} that appear, ascending.}
#'     \item{`c`}{their coefficients, alternating in sign.}
#'   }
#' @return The list described under Format. It is a stored constant, never
#'   called, and is documented so that the expansion can be checked against its
#'   source.
#'
#' @references
#' Olver, F. W. J., et al. (2024). *NIST Digital Library of Mathematical
#' Functions*, section 10.41. <https://dlmf.nist.gov/10.41>.
#'
#' @seealso [log_bessel_i()] and [log_bessel_k()], whose large-order branches
#'   sum these.
#'
#' @keywords internal
#' @name dot-bessel_uk
#' @aliases .bessel_uk
NULL

.bessel_uk <- list(
  list(e = c(1, 3),
       c = c(0.125, -0.20833333333333334)),
  list(e = c(2, 4, 6),
       c = c(0.0703125, -0.4010416666666667, 0.3342013888888889)),
  list(e = c(3, 5, 7, 9),
       c = c(0.0732421875, -0.8912109375, 1.8464626736111112, -1.0258125964506173)),
  list(e = c(4, 6, 8, 10, 12),
       c = c(0.112152099609375, -2.3640869140625, 8.78912353515625, -11.207002616222994, 4.669584423426247)),
  list(e = c(5, 7, 9, 11, 13, 15),
       c = c(0.22710800170898438, -7.368794359479632, 42.53499874538846, -91.81824154324002, 84.63621767460073, -28.212072558200244)),
  list(e = c(6, 8, 10, 12, 14, 16, 18),
       c = c(0.5725014209747314, -26.491430486951554, 218.1905117442116, -699.5796273761325, 1059.9904525279999, -765.2524681411817, 212.57013003921713)),
  list(e = c(7, 9, 11, 13, 15, 17, 19, 21),
       c = c(1.7277275025844574, -108.09091978839466, 1200.9029132163525, -5305.646978613403, 11655.393336864534, -13586.550006434138, 8061.722181737309, -1919.457662318407)),
  list(e = c(8, 10, 12, 14, 16, 18, 20, 22, 24),
       c = c(6.074042001273483, -493.915304773088, 7109.514302489364, -41192.65496889755, 122200.46498301746, -203400.17728041555, 192547.00123253153, -96980.59838863752, 20204.29133096615)),
  list(e = c(9, 11, 13, 15, 17, 19, 21, 23, 25, 27),
       c = c(24.380529699556064, -2499.8304818112097, 45218.76898136273, -331645.1724845636, 1268365.2733216248, -2813563.226586534, 3763271.297656404, -2998015.9185381066, 1311763.6146629772, -242919.18790055133)),
  list(e = c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
       c = c(110.01714026924674, -13886.08975371704, 308186.4046126624, -2785618.1280864547, 13288767.166421818, -37567176.66076335, 66344512.27472903, -74105148.21153265, 50952602.49266464, -19706819.118432228, 3284469.853072038)),
  list(e = c(11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33),
       c = c(551.3358961220206, -84005.43360302408, 2243768.1779224495, -24474062.72573873, 142062907.7975331, -495889784.2750303, 1106842816.8230145, -1621080552.1083372, 1553596899.57058, -939462359.6815784, 325573074.18576574, -49329253.66450996)),
  list(e = c(12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36),
       c = c(3038.090510922384, -549842.3275722887, 17395107.553978164, -225105661.88941526, 1559279864.8792574, -6563293792.619285, 17954213731.1556, -33026599749.800724, 41280185579.753975, -34632043388.158775, 18688207509.295826, -5866481492.051847, 814789096.1183121)),
  list(e = c(13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39),
       c = c(18257.755474293175, -3871833.442572613, 143157876.71888897, -2167164983.223795, 17634730606.83497, -87867072178.02327, 287900649906.1506, -645364869245.3765, 1008158106865.3821, -1098375156081.2233, 819218669548.5773, -399096175224.4665, 114498237732.0258, -14679261247.695616))
)
