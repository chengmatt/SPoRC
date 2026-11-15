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
#' With a projection set up by \code{\link{Setup_Mod_Projection}}, the model is
#' fitted with the projection off, since it adds nothing to the likelihood, then
#' rebuilt at the estimates with the projection on. RTMB records a model once, at
#' its starting values, and a rule or cap that compares two numbers keeps whichever
#' side held then, so the rebuild is what puts the estimated stock on the right
#' side. \code{do_sdrep = TRUE} then covers the projected quantities.
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
#'   from. \code{"exact"} (default) is the behavior above. \code{"diff"}
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
#'   Newton-refinement routines. A \code{data} list without the dims
#'   \code{\link{Setup_Mod_Dim}} writes (\code{years}, \code{n_regions},
#'   \code{ages}, \code{lens}, \code{n_sexes}, \code{n_fish_fleets} and
#'   \code{n_srv_fleets}) is taken as such a model, and the SPoRC setup steps
#'   are skipped: filling in movement and tagging when they were not set up,
#'   and refreshing the copies of the map kept in the data list. A SPoRC data
#'   list with a model of your own, such as an edited copy of
#'   \code{SPoRC_rtmb}, still takes them.
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
#' @param tape_config Named list of how RTMB records the model, passed to
#'   \code{RTMB::TapeConfig()} before it is built. \code{NULL} (default)
#'   leaves the session's settings as they are. Elements, each a single string:
#'   \describe{
#'     \item{\code{vectorize}}{\code{"enable"} records arithmetic, math
#'       functions, sums and sparse matrix products as one operation per
#'       vector rather than one per element. \code{"disable"} is RTMB's default.}
#'     \item{\code{atomic}}{\code{"enable"} records each matrix product as one
#'       operation. \code{"disable"} records matrix products and multivariate
#'       normal densities element by element; RTMB's default records the
#'       densities as one operation and the products element by element, so
#'       \code{"disable"} also changes the densities.}
#'     \item{\code{comparison}}{What happens when the model compares two
#'       quantities it is differentiating: \code{"forbid"} (RTMB's default)
#'       stops with an error, \code{"tape"} records the comparison as an
#'       indicator, and \code{"allow"}, for testing only, decides it once at
#'       the starting values.}
#'   }
#'   \code{"NA"} leaves an element as set. Each recording gives the same
#'   objective and gradient, to rounding. A vectorized or atomic operation is
#'   taken as making every output depend on every input, so with
#'   \code{random} set, \code{vectorize = "enable"} makes the Hessian of the
#'   random effects dense, and \code{atomic = "enable"} does too when one
#'   product spans many random effects; each Laplace gradient is then several
#'   times slower. Time a fit both ways before keeping either.
#' @param tmb_config Named list of TMB settings for the models RTMB builds,
#'   passed to \code{TMB::config(..., DLL = "RTMB")} before the model is
#'   built. \code{NULL} (default) leaves them as they are. Any name
#'   \code{TMB::config(DLL = "RTMB")} returns is allowed, each a single 0 or 1
#'   (or a count for \code{nthreads}). The ones that change a fit:
#'   \describe{
#'     \item{\code{tmbad.sparse_hessian_compress}}{\code{1} compresses the
#'       sparse Hessian of the random effects, which takes less memory and can
#'       be slower. Read only when \code{random} is set. Default \code{0}.}
#'     \item{\code{tmbad.atomic_sparse_log_determinant}}{\code{0} computes the
#'       Laplace log determinant from the plain factorization rather than TMB's
#'       single operation, which has no second derivative.
#'       \code{\link{Setup_Mod_DSEM}} sets \code{0} itself when a moderated
#'       arrow needs it, and setting \code{1} then gives NaN. Default \code{1}.}
#'     \item{\code{optimize.instantly}}{\code{0} skips removing repeated and
#'       unused operations from the recorded model, which is slower but keeps
#'       the operations in the order they were recorded, so each can be traced
#'       to the line that recorded it. Default \code{1}.}
#'   }
#'   Both lists are in force for the builds, Newton steps and \code{sdreport}
#'   inside \code{fit_model}, and the earlier settings are restored when it
#'   returns or stops. The diagnostics that refit the model
#'   (\code{\link{do_retrospective}}, \code{\link{do_jitter}},
#'   \code{\link{do_likelihood_profile}}) do not record / save these settings
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
  tape_config = NULL,
  tmb_config = NULL,
  ...
) {

  # check par n map len
  check_par_map_lengths(parameters, mapping)

  # check tape settings before fitting, so a misspelled name or value stops here rather than inside RTMB or TMB
  check_fit_tape_settings(tape_config = tape_config, tmb_config = tmb_config)

  # modifyList below would drop a misspelled tolerance without a word, so stop on it before fitting
  tol_names <- names(post_optim_sanity_checks_control)
  if(is.null(tol_names)) tol_names <- rep("", length(post_optim_sanity_checks_control))
  bad_tol <- setdiff(tol_names, c("gradient_tol", "se_tol", "corr_tol"))
  if(length(bad_tol) > 0) stop("post_optim_sanity_checks_control takes only gradient_tol, se_tol and corr_tol, not: ",
                               paste0("'", bad_tol, "'", collapse = ", "))

  # a data list without the dims Setup_Mod_Dim has always written is a model of the user's own, which has none of the setup fields below
  sporc_data <- all(c("years", "n_regions", "ages", "lens", "n_sexes", "n_fish_fleets", "n_srv_fleets") %in% names(data))

  if(sporc_data) {

    # movement, tagging and the dsem are optional, so figure out what each needs when not set up
    filled <- fill_optional_modules(data, parameters, mapping)
    data <- filled$data
    parameters <- filled$par
    mapping <- filled$map

    # refresh data mapping stuff if mapping declared after setup
    data <- sync_dev_map_data(data, mapping)

  } # end if sporc_data

  # the projection adds nothing to the likelihood, so the model is fitted without it
  proj_on <- !is.null(data$proj_n_yrs) && data$proj_n_yrs > 0
  data_fit <- data # new data object for fitting only
  if(proj_on) data_fit$proj_n_yrs <- 0 # set projection years to zero here if actually projection since the devs aren't actually fit

  # tape settings are global to the session, so set them for the builds, newton steps and sdreport below and put the old ones back on exit
  if(!is.null(tape_config)) {
    old_tape_config <- RTMB::TapeConfig() # no arguments reads the current settings
    on.exit(RTMB::TapeConfig(old_tape_config), add = TRUE) # a list as the only argument restores all seven
    do.call(RTMB::TapeConfig, tape_config)
  } # end if tape_config

  if(!is.null(tmb_config)) {
    old_tmb_config <- TMB::config(DLL = "RTMB") # no settings reads the current ones
    on.exit(do.call(TMB::config, c(old_tmb_config, list(DLL = "RTMB"))), add = TRUE)
    do.call(TMB::config, c(tmb_config, list(DLL = "RTMB")))
  } # end if tmb_config

  # make AD model function
  obj <- RTMB::MakeADFun(
    cmb(model, data_fit),
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
                   optim$objective = as.numeric(obj$fn(optim$par))
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
                 optim$objective = as.numeric(obj$fn(optim$par))
               }
             , error = function(e) {e}, warning = function(w) {w})

    # return obj estimate here
    optim$objective = as.numeric(obj$fn(optim$par))

    # record bounds used alongside optim output
    optim$lower <- lower
    optim$upper <- upper

    # save optim
    obj$optim <- optim
  }

  # make sure to rename if naming got stripped by nlminb
  if (is.null(names(obj$env$last.par.best)) && !is.null(names(obj$par)) && length(obj$env$last.par.best) == length(obj$par)) {
    names(obj$env$last.par.best) <- names(obj$par)
  }

  # RTMB records arithmetic only and an if in the an HCR or a cap is decided by R on the plain numbers (ad_value) while
  # recording, so the record keeps the estiamtes bascially held at the starting values (B/B40 above one, catch over the cap).
  # rebuild at the estimates
  if(proj_on) {
    obj_proj <- RTMB::MakeADFun(
      cmb(model, data), # use actual data (fitting + projection options to do the projections)
      parameters = obj$env$parList(),
      map = mapping,
      random = random,
      silent = silent,
      ...
    )
    obj_proj$optim <- obj$optim
    obj_proj$hessian <- obj$hessian # the likelihood's, which the projection leaves as it was
    obj <- obj_proj # replace obj with obj_proj here
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

#' Check the tape settings passed to fit_model
#'
#' \code{RTMB::TapeConfig} stops on a bad value only once the fit has begun,
#' and \code{TMB::config} keeps a misspelled name without a word, so both lists
#' are checked before anything is built and the error names the argument.
#'
#' @param tape_config \code{NULL} or a named list of \code{comparison},
#'   \code{atomic} and \code{vectorize}, as in \code{\link{fit_model}}.
#' @param tmb_config \code{NULL} or a named list of TMB settings, as in
#'   \code{\link{fit_model}}.
#'
#' @return \code{NULL}, invisibly. Stops on any name or value not allowed.
#'
#' @keywords internal
check_fit_tape_settings <- function(tape_config,
                                    tmb_config) {

  # rtmb tape settings, each with the values RTMB::TapeConfig accepts
  if(!is.null(tape_config)) {
    tape_values <- list(comparison = c("NA", "forbid", "tape", "allow"), # how a comparison of two model quantities is recorded
                        atomic = c("NA", "enable", "disable"), # matrix products and multivariate normal densities
                        vectorize = c("NA", "disable", "enable")) # arithmetic, math functions, sums and sparse products

    if(!is.list(tape_config)) stop("tape_config must be NULL or a named list, such as list(vectorize = 'enable').")
    tape_names <- names(tape_config)
    if(is.null(tape_names)) tape_names <- rep("", length(tape_config))
    bad_tape_names <- setdiff(tape_names, names(tape_values))
    if(length(bad_tape_names) > 0) stop("tape_config takes only comparison, atomic and vectorize, not: ",
                                        paste0("'", bad_tape_names, "'", collapse = ", "))

    for(i in seq_along(tape_config)) {
      setting_name <- tape_names[i] # comparison, atomic or vectorize
      setting_value <- tape_config[[i]] # value asked for
      if(!(is.character(setting_value) && length(setting_value) == 1 && setting_value %in% tape_values[[setting_name]])) {
        stop("tape_config$", setting_name, " must be one of ",
             paste0("'", tape_values[[setting_name]], "'", collapse = ", "), ".")
      }
    } # end i loop
  } # end if tape_config

  # tmb settings, named as TMB::config reads them back for the RTMB library
  if(!is.null(tmb_config)) {
    tmb_names_allowed <- names(TMB::config(DLL = "RTMB")) # no settings reads the current ones

    if(!is.list(tmb_config)) stop("tmb_config must be NULL or a named list, such as list(tmbad.sparse_hessian_compress = 1).")
    tmb_names <- names(tmb_config)
    if(is.null(tmb_names)) tmb_names <- rep("", length(tmb_config))
    bad_tmb_names <- setdiff(tmb_names, tmb_names_allowed)
    if(length(bad_tmb_names) > 0) stop("tmb_config takes only ", paste0("'", sort(tmb_names_allowed), "'", collapse = ", "),
                                       ", not: ", paste0("'", bad_tmb_names, "'", collapse = ", "))

    for(i in seq_along(tmb_config)) {
      setting_value <- tmb_config[[i]] # value asked for, 0 or 1 for a switch, a count for nthreads
      if(!((is.numeric(setting_value) || is.logical(setting_value)) && length(setting_value) == 1 && !is.na(setting_value))) {
        stop("tmb_config$", tmb_names[i], " must be a single number or TRUE/FALSE.")
      }
    } # end i loop
  } # end if tmb_config

  return(invisible(NULL))

} # end function

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
