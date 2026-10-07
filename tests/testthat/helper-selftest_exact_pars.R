# A joint self test at parameters known exactly: a precision of 1e12 on every parameter makes each replicate's
# draw the fit's own vector to within 1e-6, so only the processes and the data are drawn fresh.

library(SPoRC)
library(testthat)

#' sdreport stand in for a joint self test whose parameters are not drawn
#'
#' @param obj Model built at the truth, its parameter vector the one each replicate runs on.
#' @param random Names of the parameters the refits integrate out.
#'
#' @return \code{par.fixed}, \code{par.random} and \code{jointPrecision}, in the order of
#'   \code{obj$env$last.par.best}.
exact_pars_sd_rep <- function(obj, random = NULL) {
  full <- obj$env$last.par.best
  is_random <- names(full) %in% random
  list(par.fixed = full[!is_random], par.random = full[is_random],
       jointPrecision = Matrix::.sparseDiagonal(length(full), 1e12, shape = "s"))
}
