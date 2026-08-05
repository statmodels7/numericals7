#' @include enumerate.R
NULL

# A jet is a number carried together with all its partial derivatives up to
# fourth order in the free values. Adding and multiplying jets propagates
# those derivatives exactly, so a quantity built from a recursion of sums and
# products -- which is what the Levinson-Durbin recursion behind an AR(p)
# parameter is -- comes out analytic at every order with no formula
# transcribed. The alternative, expanding the recursion by hand, is pages of
# algebra whose only check would be the very finite differences the toolkit
# refuses to trust at fourth order.


#' The Bookkeeping a Jet Needs
#'
#' @description
#' The index tuples up to fourth order over \eqn{d} variables, together with
#' the lookup from a tuple to its position, computed once and shared by every
#' jet in a calculation.
#'
#' @details
#' The tuples come from \code{\link{tuple_indices}}, the one enumeration the
#' toolkit shares, so a jet's components are keyed exactly as any other holder
#' of derivatives over the same variables and nothing has to be reordered on
#' the way out.
#'
#' @param d The number of variables.
#'
#' @return A list with \code{d}, the tuple lists \code{tuples}, and the
#'   environment \code{pos} mapping a sorted tuple to its order and position.
#'
#' @seealso \code{\link{jet_var}}, \code{\link{Ops.jet}}
#'
#' @examples
#' lay <- jet_layout(2)
#' lay$d
#' length(lay$tuples[[2]])
#'
#' @export
jet_layout <- function(d) {
  tuples <- lapply(1:4, function(o) tuple_indices(d, o))
  # The empty tuple gets no key: R has no zero-length variable name, and the
  # lookups short-circuit on it before reaching the environment anyway.
  pos <- new.env(parent = emptyenv())
  for (o in 1:4) {
    for (i in seq_along(tuples[[o]])) {
      assign(paste(sort(tuples[[o]][[i]]), collapse = ","), c(o, i), envir = pos)
    }
  }
  list(d = d, tuples = tuples, pos = pos)
}


#' A Constant Jet
#'
#' @description
#' A number with every derivative zero.
#'
#' @param v The value.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#'
#' @return A jet.
#'
#' @seealso \code{\link{jet_var}}, \code{\link{jet_layout}}
#'
#' @examples
#' lay <- jet_layout(2)
#' x <- jet_var(1, list(2, 1, 0, 0, 0), lay)
#' acc <- jet_const(0, lay)   # a running sum starts at a constant
#' (acc + x * 3)$v
#'
#' @export
jet_const <- function(v, lay) {
  structure(
    list(v = v,
         d = lapply(1:4, function(o) numeric(length(lay$tuples[[o]]))),
         lay = lay),
    class = "jet"
  )
}


#' A Jet for One Variable
#'
#' @description
#' The jet of a scalar map of a single free value: its value and its four
#' derivatives sit in the pure components of that variable, everything else
#' being zero.
#'
#' @param k The index of the free value.
#' @param dv A list of five numbers: the value and four derivatives.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#'
#' @return A jet.
#'
#' @seealso \code{\link{jet_layout}}, \code{\link{Ops.jet}}
#'
#' @examples
#' # seed two variables, then write a map in ordinary R
#' lay <- jet_layout(2)
#' x <- jet_var(1, list(1.3, 1, 0, 0, 0), lay)
#' y <- jet_var(2, list(0.8, 1, 0, 0, 0), lay)
#' z <- x / gamma(1 + 1 / y)
#' z$v
#' z$d[[1]]
#'
#' @export
jet_var <- function(k, dv, lay) {
  out <- jet_const(dv[[1L]], lay)
  for (o in 1:4) {
    slot <- get(paste(rep(k, o), collapse = ","), envir = lay$pos)
    out$d[[o]][slot[2L]] <- dv[[o + 1L]]
  }
  out
}


#' Sum of Jets
#'
#' @description
#' Adds two jets, or a jet and a number, componentwise.
#'
#' @param a,b Jets, or a jet and a single number.
#'
#' @return A jet.
#'
#' @keywords internal
jet_add <- function(a, b) {
  structure(
    list(v = a$v + b$v,
         d = lapply(1:4, function(o) a$d[[o]] + b$d[[o]]),
         lay = a$lay),
    class = "jet"
  )
}


#' Product of Jets
#'
#' @description
#' Multiplies two jets, propagating every derivative exactly.
#'
#' @details
#' The Leibniz rule again, in its scalar form: the component of the product at
#' a tuple \eqn{T} is \eqn{\sum_{S \subseteq T} a_S b_{T \setminus S}}, the sum
#' running over subsets of \emph{positions} so that a repeated variable is
#' counted with the right multiplicity without a bookkeeping of its own. It is
#' the same enumeration parameters7's Leibniz machinery uses for a matrix
#' product.
#'
#' @param a,b Jets over the same layout.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#'
#' @return A jet.
#'
#' @keywords internal
jet_mul <- function(a, b, lay = a$lay) {
  grab <- function(x, tup) {
    if (!length(tup)) return(x$v)
    slot <- get(paste(sort(tup), collapse = ","), envir = lay$pos)
    x$d[[slot[1L]]][slot[2L]]
  }
  out <- jet_const(a$v * b$v, lay)
  for (o in 1:4) {
    tt <- lay$tuples[[o]]
    vals <- numeric(length(tt))
    positions <- seq_len(o)
    for (i in seq_along(tt)) {
      t <- tt[[i]]
      acc <- 0
      for (m in 0:(2^o - 1L)) {
        take <- as.logical(bitwAnd(m, bitwShiftL(1L, positions - 1L)) > 0L)
        acc <- acc + grab(a, t[take]) * grab(b, t[!take])
      }
      vals[i] <- acc
    }
    out$d[[o]] <- vals
  }
  out
}


#' A Jet Times a Constant
#'
#' @description
#' Scales a jet by a plain number, which is the common case and avoids
#' building a constant jet for it.
#'
#' @param a A jet.
#' @param c A single number.
#'
#' @return A jet.
#'
#' @keywords internal
jet_cmul <- function(a, c) {
  structure(list(v = a$v * c, d = lapply(a$d, function(x) x * c), lay = a$lay),
            class = "jet")
}


#' A Smooth Function of a Jet
#'
#' @description
#' Applies a scalar function to a jet, given that function's own derivatives at
#' the jet's value, and propagates every partial derivative exactly.
#'
#' @details
#' This is Faa di Bruno in the form the multivariate case takes when the inner
#' function is scalar valued: for a tuple \eqn{T} of positions,
#' \deqn{(f \circ a)_T = \sum_{\pi} f^{(|\pi|)}(a)\prod_{B \in \pi} a_{T[B]}}
#' the sum running over the set partitions of the positions of \eqn{T}.
#'
#' Every transcendental a jet needs is one call to this with the right five
#' numbers, so the exponential, the logarithm, a power and the gamma function
#' need no derivative machinery of their own. Extending the vocabulary is a
#' matter of writing five derivatives, not of writing a chain rule.
#'
#' @param a A jet.
#' @param fd A numeric vector of five values: \eqn{f} and its first four
#'   derivatives, all evaluated at \code{a$v}.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#'
#' @return A jet.
#'
#' @seealso \code{\link{jet_mul}}, \code{\link{set_partitions}}
#'
#' @keywords internal
jet_compose <- function(a, fd, lay = a$lay) {
  grab <- function(tup) {
    slot <- get(paste(sort(tup), collapse = ","), envir = lay$pos)
    a$d[[slot[1L]]][slot[2L]]
  }
  out <- jet_const(fd[[1L]], lay)
  for (o in 1:4) {
    parts <- set_partitions(o)
    tt <- lay$tuples[[o]]
    vals <- numeric(length(tt))
    for (i in seq_along(tt)) {
      t <- tt[[i]]
      acc <- 0
      for (p in parts) {
        term <- fd[[length(p) + 1L]]
        if (term == 0) next
        for (b in p) term <- term * grab(t[b])
        acc <- acc + term
      }
      vals[i] <- acc
    }
    out$d[[o]] <- vals
  }
  out
}


#' The Exponential of a Jet
#'
#' @description \eqn{e^a}, whose four derivatives are all \eqn{e^a}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_exp <- function(a, lay = a$lay) {
  e <- exp(a$v)
  jet_compose(a, rep(e, 5L), lay)
}


#' The Logarithm of a Jet
#'
#' @description \eqn{\log a}, with derivatives \eqn{(-1)^{k-1}(k-1)!/a^k}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_log <- function(a, lay = a$lay) {
  v <- a$v
  jet_compose(a, c(log(v), 1 / v, -1 / v^2, 2 / v^3, -6 / v^4), lay)
}


#' A Power of a Jet
#'
#' @description
#' \eqn{a^p} for any real exponent, which covers the square root, the
#' reciprocal and the square without a routine for each.
#' @param a A jet.
#' @param p The exponent.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_pow <- function(a, p, lay = a$lay) {
  v <- a$v
  fd <- c(v^p,
          p * v^(p - 1),
          p * (p - 1) * v^(p - 2),
          p * (p - 1) * (p - 2) * v^(p - 3),
          p * (p - 1) * (p - 2) * (p - 3) * v^(p - 4))
  jet_compose(a, fd, lay)
}


#' The Reciprocal of a Jet
#'
#' @description \eqn{1/a}, the common case of \code{\link{jet_pow}}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_pow}}
#' @keywords internal
jet_inv <- function(a, lay = a$lay) jet_pow(a, -1, lay)


#' The Square Root of a Jet
#'
#' @description \eqn{\sqrt{a}}, the common case of \code{\link{jet_pow}}.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_pow}}
#' @keywords internal
jet_sqrt <- function(a, lay = a$lay) jet_pow(a, 0.5, lay)


#' The Quotient of Two Jets
#'
#' @description \eqn{a/b}, as a product with the reciprocal.
#' @param a,b Jets over the same layout.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_mul}}
#' @keywords internal
jet_div <- function(a, b, lay = a$lay) jet_mul(a, jet_inv(b, lay), lay)


#' The Log-Gamma of a Jet
#'
#' @description
#' \eqn{\log\Gamma(a)}, whose derivatives are the polygamma functions.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_lgamma <- function(a, lay = a$lay) {
  v <- a$v
  jet_compose(a, c(lgamma(v), digamma(v), trigamma(v),
                   psigamma(v, 2L), psigamma(v, 3L)), lay)
}


#' The Gamma Function of a Jet
#'
#' @description
#' \eqn{\Gamma(a)}, formed as \eqn{\exp\log\Gamma(a)} so that the derivatives
#' come from the polygamma functions rather than from a second table.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_lgamma}}
#' @keywords internal
jet_gamma <- function(a, lay = a$lay) jet_exp(jet_lgamma(a, lay), lay)


#' The Digamma of a Jet
#'
#' @description \eqn{\psi(a)}, with the polygamma functions above it.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_digamma <- function(a, lay = a$lay) {
  v <- a$v
  jet_compose(a, c(digamma(v), trigamma(v), psigamma(v, 2L),
                   psigamma(v, 3L), psigamma(v, 4L)), lay)
}


#' The Trigamma of a Jet
#'
#' @description \eqn{\psi'(a)}, with the polygamma functions above it.
#' @param a A jet.
#' @param lay A layout from \code{\link{jet_layout}}; taken from the jet
#'   itself unless given.
#' @return A jet.
#' @seealso \code{\link{jet_compose}}
#' @keywords internal
jet_trigamma <- function(a, lay = a$lay) {
  v <- a$v
  jet_compose(a, c(trigamma(v), psigamma(v, 2L), psigamma(v, 3L),
                   psigamma(v, 4L), psigamma(v, 5L)), lay)
}


#' Arithmetic on Jets
#'
#' @description
#' The arithmetic operators applied to jets, so that a map is written in
#' ordinary R and its derivatives come out exact.
#'
#' @details
#' This is what makes a jet worth having. A quantity written as
#' \code{mu / gamma(1 + 1 / sigma)} carries every partial derivative to fourth
#' order with it, and whoever wrote the expression need not know that: the
#' operators dispatch, and nothing about the derivative machinery appears in
#' the expression.
#'
#' A number mixed with a jet is treated as a constant. Multiplication and
#' division by one take the cheap route rather than building a constant jet and
#' running the full Leibniz rule. A jet raised to a numeric power goes through
#' \code{\link{jet_pow}}; a jet raised to a jet through
#' \eqn{\exp(b\log a)}.
#'
#' Comparison and logical operators are refused. A jet is a value with its
#' derivatives attached, and a branch taken on it would silently discard them:
#' whichever side the comparison chose, the result would carry the derivatives
#' of that side alone and claim to be the derivative of the whole expression.
#'
#' @param e1,e2 Jets, or a jet and a number.
#'
#' @return A jet.
#'
#' @seealso \code{\link{Math.jet}}, \code{\link{jet_compose}}
#'
#' @keywords internal
#' @export
Ops.jet <- function(e1, e2) {
  op <- .Generic
  if (missing(e2)) {
    if (op == "+") return(e1)
    if (op == "-") return(jet_cmul(e1, -1))
    stop(sprintf("unary '%s' is not defined for a jet.", op), call. = FALSE)
  }
  if (!op %in% c("+", "-", "*", "/", "^")) {
    stop(sprintf(paste0(
      "'%s' is not defined for a jet. A jet is a value with its derivatives\n",
      "  attached, so a branch taken on it would keep one side's derivatives\n",
      "  and report them as the whole expression's."
    ), op), call. = FALSE)
  }
  j1 <- inherits(e1, "jet")
  j2 <- inherits(e2, "jet")
  lay <- if (j1) e1$lay else e2$lay
  as_jet <- function(x) if (inherits(x, "jet")) x else jet_const(x, lay)

  switch(op,
    "+" = jet_add(as_jet(e1), as_jet(e2)),
    "-" = jet_add(as_jet(e1), jet_cmul(as_jet(e2), -1)),
    "*" = if (!j1) jet_cmul(e2, e1) else if (!j2) jet_cmul(e1, e2) else
      jet_mul(e1, e2, lay),
    "/" = if (!j2) jet_cmul(e1, 1 / e2) else jet_div(as_jet(e1), e2, lay),
    "^" = if (j2) jet_exp(jet_mul(e2, jet_log(as_jet(e1), lay), lay), lay) else
      jet_pow(as_jet(e1), e2, lay)
  )
}


#' Mathematical Functions of a Jet
#'
#' @description
#' The mathematical functions a jet understands, so that a map written in
#' ordinary R differentiates itself.
#'
#' @details
#' Each is one call to \code{\link{jet_compose}} with the function's own five
#' derivatives, so extending the vocabulary means writing five derivatives and
#' not a chain rule.
#'
#' A function that is not smooth is refused rather than approximated:
#' \code{abs}, \code{floor}, \code{ceiling}, \code{round} and \code{sign} have
#' no fourth derivative to carry, and returning one would be a wrong answer in
#' the shape of a right one.
#'
#' @param x A jet.
#' @param ... Further arguments; refused, since a second argument would change
#'   what is being differentiated.
#'
#' @return A jet.
#'
#' @seealso \code{\link{Ops.jet}}, \code{\link{jet_compose}}
#'
#' @keywords internal
#' @export
Math.jet <- function(x, ...) {
  if (length(list(...))) {
    stop(sprintf("'%s' takes no further arguments on a jet.", .Generic),
         call. = FALSE)
  }
  v <- x$v
  fd <- switch(.Generic,
    exp = rep(exp(v), 5L),
    log = c(log(v), 1 / v, -1 / v^2, 2 / v^3, -6 / v^4),
    sqrt = NULL,
    gamma = NULL,
    lgamma = c(lgamma(v), digamma(v), trigamma(v),
               psigamma(v, 2L), psigamma(v, 3L)),
    digamma = c(digamma(v), trigamma(v), psigamma(v, 2L),
                psigamma(v, 3L), psigamma(v, 4L)),
    trigamma = c(trigamma(v), psigamma(v, 2L), psigamma(v, 3L),
                 psigamma(v, 4L), psigamma(v, 5L)),
    sin = c(sin(v), cos(v), -sin(v), -cos(v), sin(v)),
    cos = c(cos(v), -sin(v), -cos(v), sin(v), cos(v)),
    expm1 = c(expm1(v), exp(v), exp(v), exp(v), exp(v)),
    log1p = c(log1p(v), 1 / (1 + v), -1 / (1 + v)^2,
              2 / (1 + v)^3, -6 / (1 + v)^4),
    stop(sprintf(paste0(
      "'%s' is not defined for a jet. A function without four derivatives\n",
      "  has nothing to carry; one that has them is added by writing those\n",
      "  five numbers in Math.jet()."
    ), .Generic), call. = FALSE)
  )
  if (.Generic == "sqrt") return(jet_pow(x, 0.5, x$lay))
  if (.Generic == "gamma") return(jet_gamma(x, x$lay))
  jet_compose(x, fd, x$lay)
}


#' Print a Jet
#'
#' @description
#' Shows the value and how many variables the derivatives are taken in, rather
#' than the whole list of components, which is thirty-five numbers even in two
#' variables.
#'
#' @param x A jet.
#' @param ... Unused.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{Ops.jet}}
#'
#' @keywords internal
#' @export
print.jet <- function(x, ...) {
  cat(sprintf("<jet> value %s, derivatives to order 4 in %d variable%s\n",
              format(x$v, digits = 7), x$lay$d,
              if (x$lay$d == 1L) "" else "s"))
  cat("  gradient: ", paste(format(x$d[[1]], digits = 5), collapse = "  "),
      "\n", sep = "")
  invisible(x)
}
