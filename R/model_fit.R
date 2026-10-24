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
#'   Default \code{3}. A step that would raise the largest absolute gradient
#'   is not taken and refinement stops there, as happens when the Hessian is
#'   not positive definite. Errors and warnings are caught silently via
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
#' @param newton_hessian Character. Where each refinement step's Hessian comes
#'   from. \code{"exact"} (default) is the behaviour above. \code{"diff"}
#'   differences the exact gradient (\code{\link[stats]{optimHess}}, steps of
#'   1e-5), for models whose tape Hessian is too large to build, and in place of
#'   the Newton step reruns \code{nlminb} with each parameter in its own
#'   standard errors (the parameters scaled by the Cholesky factor of the
#'   inverse Hessian). A differenced Hessian is not accurate enough for a raw
#'   Newton step when the curvature spans many orders of magnitude, but it
#'   makes the surface round for \code{nlminb}, which then converges in a few
#'   iterations. The last Hessian is returned as \code{$hessian}, for
#'   \code{RTMB::sdreport(obj, hessian.fixed = obj$hessian)}.
#' @param silent Logical. If \code{TRUE}, suppresses RTMB and optimizer
#'   console output. Default \code{FALSE}.
#' @param do_optim Logical. If \code{TRUE} (default), runs \code{nlminb}
#'   and Newton refinement. If \code{FALSE}, returns the un-optimized
#'   \code{MakeADFun} object only.
#' @param nlminb_control Named list of control parameters passed to
#'   \code{stats::nlminb}. Default
#'   \code{list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-12)}.
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
#' @param do_sdrep Logical. If \code{TRUE}, runs \code{RTMB::sdreport(obj)} at
#'   the fitted parameters, returns it as \code{$sd_rep}, and runs
#'   \code{\link{post_optim_sanity_checks}} on it. Default \code{FALSE}, which
#'   skips both. \code{sdreport} differences the gradient for its Hessian, so on
#'   an ill-conditioned fit it can report a Hessian that is not positive definite
#'   where the tape Hessian is; rerun it as
#'   \code{RTMB::sdreport(obj, hessian.fixed = obj$he())} before reading that
#'   check as a real failure.
#' @param post_optim_sanity_checks_control Named list of the tolerances passed to
#'   \code{\link{post_optim_sanity_checks}} under \code{do_sdrep = TRUE}:
#'   \code{gradient_tol} (largest absolute gradient), \code{se_tol} (largest
#'   standard error) and \code{corr_tol} (largest absolute correlation between
#'   two parameters). Default
#'   \code{list(gradient_tol = 1e-3, se_tol = 100, corr_tol = 0.99)}. A list
#'   naming only some of them, such as \code{list(gradient_tol = 1e-4)}, takes
#'   the defaults for the rest. Any other element, named or not, stops the fit
#'   with an error before it starts. The checks print a message for each one
#'   that fails, or a single message when all pass, whatever \code{silent} is
#'   set to.
#' @param ... Additional arguments forwarded to \code{RTMB::MakeADFun}.
#'
#' @return The RTMB \code{ADFun} object with additional fields: \code{$optim}
#'   (the \code{nlminb} output list, with \code{$lower}/\code{$upper} recording
#'   the bounds used), \code{$rep} (the model report evaluated at
#'   \code{obj$env$last.par.best}), \code{$data}, \code{$parameters},
#'   \code{$mapping}, \code{$random}, \code{$provenance} (the package
#'   version and commit that built the fit, see \code{\link{fit_provenance}}),
#'   \code{$hessian} (the last refinement step's Hessian; \code{NULL} if
#'   none was taken), and under \code{do_sdrep = TRUE} \code{$sd_rep} (the
#'   \code{sdreport}) and \code{$post_optim_sanity_checks} (\code{TRUE} if
#'   every check passed, \code{FALSE} otherwise). Both are \code{NULL} under
#'   \code{do_sdrep = FALSE}.
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
  newton_hessian = "exact",
  silent = FALSE,
  do_optim = TRUE,
  nlminb_control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-12),
  lower = NULL,
  upper = NULL,
  model = SPoRC_rtmb,
  do_sdrep = FALSE,
  post_optim_sanity_checks_control = list(gradient_tol = 1e-3, se_tol = 100, corr_tol = 0.99),
  ...
) {

  # check par n map len
  check_par_map_lengths(parameters, mapping)

  # modifyList below would drop a misspelled tolerance without a word, so stop on it before fitting
  tol_names <- names(post_optim_sanity_checks_control)
  if(is.null(tol_names)) tol_names <- rep("", length(post_optim_sanity_checks_control))
  bad_tol <- setdiff(tol_names, c("gradient_tol", "se_tol", "corr_tol"))
  if(length(bad_tol) > 0) stop("post_optim_sanity_checks_control takes only gradient_tol, se_tol and corr_tol, not: ",
                               paste0("'", bad_tol, "'", collapse = ", "))

  # movement, tagging and the dsem are optional, so figure out what each needs when not set up
  filled <- fill_optional_modules(data, parameters, mapping)
  data <- filled$data
  parameters <- filled$par
  mapping <- filled$map

  # refresh data mapping stuff if mapping declared after setup
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
    if(!newton_hessian %in% c("exact", "diff")) stop("newton_hessian must be exact or diff")
    tryCatch(expr =
               for(i in seq_len(newton_loops)) { # 1:0 would run two steps when asked for none
                 g = as.numeric(obj$gr(optim$par))

                 # differenced Hessian: nlminb again with each parameter in its own standard errors, par + C y, where the
                 # surface is round; the remaining gain is tiny, so a tight tolerance
                 if(newton_hessian == "diff") {
                   h = stats::optimHess(optim$par, fn = obj$fn, gr = obj$gr, control = list(ndeps = rep(1e-5, length(g))))
                   if(!all(is.finite(h))) break
                   obj$hessian = h
                   C = t(chol(solve(h)))
                   x0 = optim$par
                   scaled = stats::nlminb(0 * x0, function(y) obj$fn(x0 + C %*% y), function(y) as.numeric(obj$gr(x0 + C %*% y) %*% C),
                                          control = utils::modifyList(nlminb_control, list(rel.tol = 1e-14)))
                   optim$par[] = pmax(lower, pmin(upper, as.numeric(x0 + C %*% scaled$par)))
                   optim$objective = obj$fn(optim$par)
                   next
                 }

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

                 obj$hessian = h
                 new_par = pmax(lower, pmin(upper, optim$par - solve(h,g))) # keep Newton step within bounds
                 if(max(abs(obj$gr(new_par))) >= max(abs(g))) break # a step that raises the gradient (an indefinite Hessian) is not taken
                 optim$par = new_par
                 optim$objective = obj$fn(optim$par)
               }
             , error = function(e) {e}, warning = function(w) {w})

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

  # if doing sdrep
  if(do_sdrep) obj$sd_rep <- RTMB::sdreport(obj)

  # attach inputs used to construct the model for reference, reuse, and serialization
  obj$data <- data
  obj$parameters <- parameters
  obj$mapping <- mapping
  obj$random <- random

  # input the fit with the package build
  obj$provenance <- fit_provenance()

  # check convergence, tolerances left out of the control list take defaults
  if(do_sdrep) {
    check_control <- utils::modifyList(list(gradient_tol = 1e-3, se_tol = 100, corr_tol = 0.99), post_optim_sanity_checks_control)
    obj$post_optim_sanity_checks <- post_optim_sanity_checks(sd_rep = obj$sd_rep, rep = obj$rep,
                                                             gradient_tol = check_control$gradient_tol,
                                                             se_tol = check_control$se_tol,
                                                             corr_tol = check_control$corr_tol
    )
  }

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

#' Record which SPoRC build produced a fit
#'
#' Reads the package version and, when the install came from GitHub through
#' remotes or pak, the commit it was built from. A package loaded with
#' \code{devtools::load_all} from a git checkout reports that checkout's HEAD
#' instead, and a plain local install has neither, so the commit is \code{NA}.
#' A saved fit then says which code made it, so a later package version can
#' tell which report names to expect.
#'
#' @return Named list: \code{package_version}, \code{commit_sha} (40 character
#'   hash or \code{NA}), \code{commit_source} (\code{"github"},
#'   \code{"local_git"}, or \code{"unknown"}), \code{r_version},
#'   \code{rtmb_version}, and \code{fit_time}.
#'
#' @keywords internal
fit_provenance <- function() {
  sha <- utils::packageDescription("SPoRC")$RemoteSha
  source <- if(is.null(sha)) "unknown" else "github"

  # load_all carries no remote fields and points system.file at inst/, so look for the
  # checkout there and one level up, and ask git for its HEAD
  if(is.null(sha)) {
    root <- system.file(package = "SPoRC")
    repo <- c(root, dirname(root))
    repo <- repo[file.exists(file.path(repo, ".git"))]
    if(length(repo) > 0 && nzchar(Sys.which("git"))) {
      head <- tryCatch(
        suppressWarnings(system2("git", c("-C", shQuote(repo[1]), "rev-parse", "HEAD"),
                                 stdout = TRUE, stderr = FALSE)),
        error = function(e) character(0)
      )
      if(length(head) == 1 && grepl("^[0-9a-f]{40}$", head)) {
        sha <- head
        source <- "local_git"
      }
    }
  }

  list(
    package_version = as.character(utils::packageVersion("SPoRC")),
    commit_sha = if(is.null(sha)) NA_character_ else sha,
    commit_source = source,
    r_version = R.version.string,
    rtmb_version = as.character(utils::packageVersion("RTMB")),
    fit_time = Sys.time()
  )
}
