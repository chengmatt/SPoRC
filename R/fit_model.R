#' Fit a SPoRC RTMB model
#'
#' Constructs an RTMB automatic differentiation function via
#' \code{RTMB::MakeADFun}, optimises it with \code{stats::nlminb}, and
#' optionally refines the solution with Newton steps using the analytic
#' Hessian. The best parameter vector (\code{obj$env$last.par.best}),
#' optimiser output, and model report are attached to the returned object.
#'
#' @param data Named list of model data as constructed by the
#'   \code{Setup_Mod_*} family of functions.
#' @param parameters Named list of parameter starting values.
#' @param mapping Named list of factor maps controlling parameter sharing
#'   and fixing.
#' @param random Character vector of parameter names to integrate out as
#'   random effects. \code{NULL} (default) fits a fixed-effects-only model.
#' @param newton_loops Integer. Number of Newton refinement steps applied
#'   after \code{nlminb} convergence to reduce gradient magnitudes. Each
#'   step solves \eqn{\Delta\theta = -H^{-1} g} and updates the objective.
#'   Default \code{3}. Errors and warnings are caught silently via
#'   \code{tryCatch}.
#' @param silent Logical. If \code{TRUE}, suppresses RTMB and optimiser
#'   console output. Default \code{FALSE}.
#' @param do_optim Logical. If \code{TRUE} (default), runs \code{nlminb}
#'   and Newton refinement. If \code{FALSE}, returns the un-optimised
#'   \code{MakeADFun} object only.
#' @param nlminb_control Named list of control parameters passed to
#'   \code{stats::nlminb}. Default
#'   \code{list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)}.
#' @param ... Additional arguments forwarded to \code{RTMB::MakeADFun}.
#' @param do_internal_comp_osa Logical. If \code{TRUE}, allows OSA residuals for composition datasets.
#' Default \code{FALSE}.
#' @param do_internal_conv_tag_osa Logical. If \code{TRUE}, allows OSA residuals for tagging datasets.
#' Default \code{FALSE}.
#'
#' @return The RTMB \code{ADFun} object with three additional fields:
#'   \code{$optim} (the \code{nlminb} output list), and \code{$rep} (the
#'   model report evaluated at \code{obj$env$last.par.best}).
#'
#' @importFrom stats nlminb optimHess
#' @export fit_model
#' @family Utility
#'
#' @examples
#' \dontrun{
#' obj <- fit_model(data, parameters, mapping, random = NULL, newton_loops = 3)
#' obj$rep$SSB
#' }
fit_model <- function(data,
                      parameters,
                      mapping,
                      random = NULL,
                      newton_loops = 3,
                      silent = FALSE,
                      do_optim = TRUE,
                      do_internal_comp_osa = FALSE,
                      do_internal_conv_tag_osa = FALSE,
                      nlminb_control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15),
                      ...
                      ) {

  # Flag for OSAs
  data$do_internal_comp_osa <- do_internal_comp_osa
  data$do_internal_conv_tag_osa <- do_internal_conv_tag_osa

  # make AD model function
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), parameters = parameters,
                         map = mapping, random = random, silent = silent, ...)

  if(do_optim == TRUE) {
    # Now, optimize the function
    optim <- stats::nlminb(obj$par, obj$fn, obj$gr,
                           control = nlminb_control)
    # newton steps
    try_improve <- tryCatch(expr =
                              for(i in 1:newton_loops) {
                                g = as.numeric(obj$gr(optim$par))
                                h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                                optim$par = optim$par - solve(h,g)
                                optim$objective = obj$fn(optim$par)
                              }
                            , error = function(e){e}, warning = function(w){w})

    # save optim
    obj$optim <- optim
  }

  # save report
  obj$rep <- obj$report(obj$env$last.par.best)

  return(obj)
}
