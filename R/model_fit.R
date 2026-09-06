# Stage 2 of 3: objective function
#
# Wraps RTMB::MakeADFun around SPoRC_rtmb, optimizes with nlminb plus Newton refinement, and
# returns the fitted object every post fit routine reads. cmb curries the data list in.

#' Fit a SPoRC RTMB model
#'
#' Constructs an RTMB automatic differentiation function via
#' \code{RTMB::MakeADFun}, optimizes it with \code{stats::nlminb}, and
#' optionally refines the solution with Newton steps using the analytic
#' Hessian. The best parameter vector (\code{obj$env$last.par.best}),
#' optimizer output, and model report are attached to the returned object.
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
#'   \code{tryCatch}, so a step that fails leaves the \code{nlminb} solution
#'   in place without a message.
#'   \eqn{H} comes from the AD tape (\code{obj$he}) for fixed-effects models,
#'   which is exact and costs a single call. Random-effects models fall back
#'   to \code{\link[stats]{optimHess}} differencing the gradient, since RTMB
#'   does not implement a tape Hessian when random effects are present.
#'   Refinement stops early if \eqn{H} comes back non-finite, which happens on
#'   models that have not converged, where second derivatives can be undefined
#'   at parameter values the objective and gradient still evaluate at. The
#'   \code{nlminb} solution is kept in that case.
#' @param silent Logical. If \code{TRUE}, suppresses RTMB and optimizer
#'   console output. Default \code{FALSE}.
#' @param do_optim Logical. If \code{TRUE} (default), runs \code{nlminb}
#'   and Newton refinement. If \code{FALSE}, returns the un-optimized
#'   \code{MakeADFun} object only.
#' @param nlminb_control Named list of control parameters passed to
#'   \code{stats::nlminb}. Default
#'   \code{list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15)}.
#' @param lower Numeric vector of lower bounds for \code{obj$par} (the
#'   estimated parameter vector, i.e. after mapping and random-effects
#'   marginalization), passed to \code{stats::nlminb} and used to clamp each
#'   Newton refinement step. \code{NULL} (default) is unbounded
#'   (\code{-Inf} for every element).
#' @param upper Numeric vector of upper bounds for \code{obj$par}, passed to
#'   \code{stats::nlminb} and used to clamp each Newton refinement step.
#'   \code{NULL} (default) is unbounded (\code{Inf} for every element).
#' @param model Function with signature \code{function(pars, data)} passed to
#'   \code{RTMB::MakeADFun} via \code{\link{cmb}}. Default \code{\link{SPoRC_rtmb}}.
#'   Allows non-SPoRC RTMB models to be fit with the same optimization and
#'   Newton-refinement routines.
#' @param ... Additional arguments forwarded to \code{RTMB::MakeADFun}.
#'
#' @return The RTMB \code{ADFun} object with additional fields: \code{$optim}
#'   (the \code{nlminb} output list, with \code{$lower}/\code{$upper} recording
#'   the bounds used), \code{$rep} (the model report evaluated at
#'   \code{obj$env$last.par.best}), and \code{$data}, \code{$parameters},
#'   \code{$mapping}, \code{$random}.
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
fit_model <- function(
  data,
  parameters,
  mapping,
  random = NULL,
  newton_loops = 3,
  silent = FALSE,
  do_optim = TRUE,
  nlminb_control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15),
  lower = NULL,
  upper = NULL,
  model = SPoRC_rtmb,
  ...
) {

  # check par n map len
  check_par_map_lengths(parameters, mapping)

  # the deviation penalties key on the map mirrors in the data list, so refresh
  # them from the map actually being handed to MakeADFun
  data <- sync_dev_map_data(data, mapping)

  # make AD model function
  obj <- RTMB::MakeADFun(
    cmb(model, data),
    parameters = parameters,
    map = mapping,
    random = random,
    silent = silent,
    ...
  )

  if(do_optim == TRUE) {

    # set bounds on optimizaiton
    if(is.null(lower)) lower <- rep(-Inf, length(obj$par))
    if(is.null(upper)) upper <- rep(Inf, length(obj$par))

    # Now, optimize the function
    optim <- stats::nlminb(
      obj$par,
      obj$fn,
      obj$gr,
      control = nlminb_control,
      lower = lower,
      upper = upper
    )

    # newton steps
    try_improve <- tryCatch(expr =
                              for(i in 1:newton_loops) {
                                g = as.numeric(obj$gr(optim$par))

                                # the tape gives the Hessian exactly in one call; optimHess
                                # differences the gradient once per parameter

                                if(is.null(random)) {
                                  h = as.matrix(obj$he(optim$par)) # analytical hessian from obj
                                } else {
                                  h = stats::optimHess(optim$par, fn = obj$fn, gr = obj$gr)
                                } # end if else for hessian source

                                # some second derivatives are undefined where the objective and
                                # gradient still are, showing up as a non-finite Hessian
                                if(!all(is.finite(h))) break

                                new_par = optim$par - solve(h,g)
                                optim$par = pmax(lower, pmin(upper, new_par)) # keep Newton step within bounds
                                optim$objective = obj$fn(optim$par)
                              }
                            , error = function(e){e}, warning = function(w){w})

    # record bounds used alongside optim output
    optim$lower <- lower
    optim$upper <- upper

    # save optim
    obj$optim <- optim
  }

  # make sure to rename if naming got stripped by nlminb
  if (is.null(names(obj$env$last.par.best)) && !is.null(names(obj$par)) &&
      length(obj$env$last.par.best) == length(obj$par)) {
    names(obj$env$last.par.best) <- names(obj$par)
  }

  # save report
  obj$rep <- obj$report(obj$env$last.par.best)

  # attach inputs used to construct the model for reference, reuse, and serialization
  obj$data <- data
  obj$parameters <- parameters
  obj$mapping <- mapping
  obj$random <- random

  return(obj)
}

#' Combine a parameter function and a data list for RTMB
#'
#' Returns a closure that calls \code{f(p, d)}, allowing the data list to be
#' fixed at construction time so that \code{RTMB::MakeADFun} receives a
#' single-argument objective function of the form \code{function(p)}.
#'
#' @param f Function with signature \code{function(pars, data)}, typically
#'   \code{\link{SPoRC_rtmb}}.
#' @param d Named list of model data passed as the second argument to
#'   \code{f} on every call.
#'
#' @return A single-argument function \code{function(p)} equivalent to
#'   \code{f(p, d)}.
#'
#' @keywords internal
cmb <- function(f, d) {
  function(p) f(p, d)
}
