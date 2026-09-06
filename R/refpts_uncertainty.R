# Stage 3 of 3: post fit
#
# Uncertainty for reference points. Get_Reference_Points takes estimated quantities in as fixed data, so its
# output has none; a reference point is a deterministic function of the fit. Albertsen and Trijoulet (2020).

#' Flatten a reference point array into a named vector
#'
#' @param x Numeric vector or array of reference point values.
#' @param prefix Character. Name stem, e.g. \code{"f_ref_pt"}.
#'
#' @return Named numeric vector. Length one input keeps the bare prefix.
#'
#' @keywords internal
flatten_refpt <- function(x, prefix) {
  x <- as.array(x)
  if(length(x) == 1) return(stats::setNames(as.numeric(x), prefix))
  idx <- expand.grid(lapply(dim(x), seq_len)) # index grid, one column per dimension
  stats::setNames(as.numeric(x), paste0(prefix, "_", apply(idx, 1, paste, collapse = "_")))
}

#' Joint covariance of the fitted parameters and states
#'
#' Selectivity is read from the terminal year, so a fit with random effects needs
#' the joint precision rather than the fixed effect block. The precision is left
#' sparse and never inverted.
#'
#' @param obj Fitted RTMB model object from \code{\link{fit_model}}.
#' @param sd_rep Optional \code{sdreport} to reuse. Needs \code{jointPrecision}
#'   when the fit has random effects.
#'
#' @return List with \code{kind} (\code{"cov"} or \code{"prec"}), \code{mat}, and
#'   \code{n}.
#'
#' @keywords internal
refpt_par_cov <- function(obj, sd_rep = NULL) {

  m <- length(obj$env$last.par.best)

  if(is.null(obj$random)) {

    sr <- if(!is.null(sd_rep)) sd_rep else if(!is.null(obj$sdrep)) obj$sdrep else RTMB::sdreport(obj)

    if(is.null(sr$cov.fixed) || any(!is.finite(sr$cov.fixed)))
      stop("cov.fixed is missing or non-finite, so the Hessian is not positive definite. ",
           "Try sdreport(obj, hessian.fixed = obj$he(obj$env$last.par.best)).", call. = FALSE)

    if(nrow(sr$cov.fixed) != m)
      stop("cov.fixed is ", nrow(sr$cov.fixed), " x ", ncol(sr$cov.fixed), " but the fit has ",
           m, " parameters. The sdreport does not belong to this model.", call. = FALSE)

    return(list(kind = "cov", mat = sr$cov.fixed, n = m))

  } # end fixed effects only

  sr <- if(!is.null(sd_rep)) sd_rep else RTMB::sdreport(obj, getJointPrecision = TRUE)

  if(is.null(sr$jointPrecision))
    stop("This fit has random effects. Pass sd_rep = sdreport(obj, getJointPrecision = TRUE).",
         call. = FALSE)

  if(nrow(sr$jointPrecision) != m)
    stop("jointPrecision is ", nrow(sr$jointPrecision), " x ", ncol(sr$jointPrecision),
         " but the fit has ", m, " parameters and states. The sdreport does not belong to this model.",
         call. = FALSE)

  list(kind = "prec", mat = sr$jointPrecision, n = m)

}

#' Delta method covariance of derived quantities
#'
#' Evaluates \eqn{d \Sigma d'}. A precision matrix is applied by solving, so the
#' dense inverse is never built.
#'
#' @param cov_obj Output of \code{\link{refpt_par_cov}}.
#' @param d Numeric matrix \code{[n_quantity, n_par]} of sensitivities.
#'
#' @return Numeric matrix \code{[n_quantity, n_quantity]}.
#'
#' @keywords internal
refpt_quad_form <- function(cov_obj, d) {
  if(cov_obj$kind == "cov") return(d %*% cov_obj$mat %*% t(d)) # if covariance
  as.matrix(d %*% Matrix::solve(cov_obj$mat, t(d))) # if precision
}

#' Draw parameter deviations from the joint covariance
#'
#' @param cov_obj Output of \code{\link{refpt_par_cov}}.
#' @param n_draw Integer. Number of draws.
#'
#' @return Numeric matrix \code{[n_draw, n_par]} of mean zero deviations.
#'
#' @keywords internal
refpt_draw_devs <- function(cov_obj, n_draw) {

  m <- cov_obj$n
  z <- matrix(stats::rnorm(n_draw * m), nrow = m, ncol = n_draw)

  # sparse Cholesky of the precision, so the dense inverse is not constructed
  if(cov_obj$kind == "prec") {
    ch <- Matrix::Cholesky(cov_obj$mat, super = TRUE)
    z <- Matrix::solve(ch, z, system = "Lt")
    return(t(as.matrix(Matrix::solve(ch, z, system = "Pt"))))
  }

  # a nearly singular covariance fails the Cholesky on rounding alone
  R <- tryCatch(chol(cov_obj$mat), error = function(e) {
    ev <- eigen(cov_obj$mat, symmetric = TRUE)
    t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors))
  })

  t(z) %*% R

}

#' Evaluate reference point quantities at a trial parameter vector
#'
#' @param obj Fitted RTMB model object.
#' @param p Numeric vector ordered as \code{obj$env$last.par.best}.
#' @param refpt_args List forwarded to \code{\link{Get_Reference_Points}}.
#' @param extra_quantities Optional \code{function(rep, refpts)} returning a named
#'   vector of extra positive quantities.
#' @param keep Optional character vector of quantity names to retain.
#'
#' @return Named numeric vector on the log scale.
#'
#' @keywords internal
eval_refpt_log_quantities <- function(obj, p, refpt_args, extra_quantities = NULL, keep = NULL) {

  rep_p <- obj$report(p) # rebuild derived quantities at p
  rp <- suppressMessages(do.call(Get_Reference_Points, c(list(data = obj$data, rep = rep_p), refpt_args)))

  out <- c(flatten_refpt(rp$f_ref_pt, "f_ref_pt"),
           flatten_refpt(rp$b_ref_pt, "b_ref_pt"),
           flatten_refpt(rp$virgin_b_ref_pt, "virgin_b_ref_pt"),
           flatten_refpt(rp$pop_b_ref_pt, "pop_b_ref_pt"),
           flatten_refpt(rp$virgin_pop_b_ref_pt, "virgin_pop_b_ref_pt"))

  if(!is.null(extra_quantities)) {
    ex <- extra_quantities(rep_p, rp)
    if(is.null(names(ex)) || any(names(ex) == ""))
      stop("extra_quantities must return a fully named numeric vector.", call. = FALSE)
    out <- c(out, ex)
  }

  if(!is.null(keep)) out <- out[keep]

  log(out)

}

#' Confidence intervals for reference points
#'
#' Attaches uncertainty to \code{\link{Get_Reference_Points}}, which returns point
#' estimates only.
#'
#' @details
#' Write \eqn{\mathbf{p}} for everything the fit estimates and
#' \eqn{x^\ast = \log F^\ast} for the reference point. The reference point has no
#' closed form, so its sensitivity \eqn{d_j = \partial x^\ast / \partial p_j} is
#' taken by central differences of the solved value, using a step of
#' \code{rel_step} times each parameter's own standard error. The delta method then
#' gives \eqn{\mathrm{Var}(x^\ast) = d \Sigma d'} with \eqn{\Sigma} the joint
#' covariance of the fit. Intervals are built on the log scale and exponentiated,
#' so they stay positive and are asymmetric.
#'
#' \code{method = "mvn"} skips the linearization and draws parameter vectors from
#' \eqn{\Sigma} instead, re-solving the reference point at each. Use it to check
#' the delta method when a reference point sits near \eqn{F_{crash}} or the stock
#' recruit curve is depensatory.
#'
#' Only estimated quantities contribute. Anything the model holds fixed, including
#' maturity, fixed \eqn{M}, fixed steepness, and anything turned off through
#' \code{mapping}, gives exactly zero, so an interval can look tight simply because
#' what drives it was fixed.
#'
#' @param obj Fitted RTMB model object from \code{\link{fit_model}}.
#' @param SPR_x Numeric. Target spawning potential ratio fraction.
#' @param t_spawn Numeric. Spawning season fraction elapsed before spawning. Default = 0.
#' @param sex_ratio_f Numeric array \code{[n_pop, n_regions]}. Defaults to 0.5.
#' @param calc_rec_st_yr Integer. First year of the mean recruitment window. Default = 1.
#' @param rec_age Integer. Recruitment lag in years. Default = 1.
#' @param type Character. \code{"single_region"} or \code{"multi_region"}.
#' @param what Character. Reference point method, as in \code{\link{Get_Reference_Points}}.
#' @param n_avg_yrs Integer. Terminal years averaged for demographic rates. Default = 1.
#' @param local_bh_msy_newton_steps Integer. Newton steps for \code{"local_MSY"}. Default = 6.
#' @param is_discard_fleet Integer vector \code{[n_fish_fleets]}. Defaults to all zeros.
#' @param sd_rep Optional \code{sdreport} to reuse. Default \code{NULL}.
#' @param method Character. \code{"delta"}, \code{"mvn"}, or \code{"both"}. Default \code{"delta"}.
#' @param rel_step Numeric. Difference step as a multiple of each parameter's
#'   standard error. Default 1e-3. Vary it and confirm \code{d} is unchanged.
#' @param min_step Numeric. Floor on the absolute step. Default 1e-7.
#' @param n_draw Integer. Draws for \code{"mvn"} and \code{"both"}. Default 300.
#' @param seed Integer or \code{NULL}. Seed for the draws. Default \code{NULL}.
#' @param level Numeric. Confidence level. Default 0.95.
#' @param extra_quantities Optional \code{function(rep, refpts)} returning a named
#'   vector of extra positive quantities to propagate through, e.g. stock status.
#' @param par_subset Optional character vector of parameter names to perturb. This is the lever for
#'   a large model if taking too long. Everything left out is asserted to have exactly zero effect,
#'   which is a claim about the model rather than a shortcut. Make sure to check any subset
#'   against \code{method = "mvn"}, whose draws perturb everything.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{refpts}}{Delta method estimates, log scale SE, bounds, and CV.}
#'     \item{\code{mvn}}{Draw based estimates and quantiles, when requested.}
#'     \item{\code{d}}{Sensitivities \code{[n_quantity, n_par]}.}
#'     \item{\code{log_cov}}{Covariance of the log quantities \code{[n_quantity, n_quantity]},
#'       needed for anything combining two of them.}
#'     \item{\code{draws}}{Log scale draws \code{[n_draw, n_quantity]}, when requested.}
#'     \item{\code{n_par}}{Parameters perturbed.}
#'     \item{\code{n_zero_par}}{How many of those had no effect on any quantity.}
#'     \item{\code{step}}{Absolute step used per parameter.}
#'   }
#'
#' @references
#' Albertsen, C.M. and Trijoulet, V. (2020). Model-based estimates of reference
#' points in an age-based state-space stock assessment model. Fisheries Research
#' 230, 105618.
#'
#' @export Get_Reference_Point_Uncertainty
#' @family Reference Points and Projections
#'
#' @examples
#' \dontrun{
#' data("dusky_rtmb_model")
#'
#' rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model, SPR_x = 0.4,
#'                                       type = "single_region", what = "SPR")
#' rp$refpts
#'
#' # have stock status through, so the correlation with the reference point is kept
#' status <- function(rep, refpts) {
#'   ssb <- rep$SSB[1, 1, dim(rep$SSB)[3]]
#'   c(SSB_terminal = ssb, status = ssb / as.numeric(refpts$b_ref_pt))
#' }
#'
#' rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model, SPR_x = 0.4,
#'                                       type = "single_region", what = "SPR",
#'                                       method = "both", extra_quantities = status)
#' }
Get_Reference_Point_Uncertainty <- function(obj,
                                            SPR_x = NULL,
                                            t_spawn = 0,
                                            sex_ratio_f = NULL,
                                            calc_rec_st_yr = 1,
                                            rec_age = 1,
                                            type,
                                            what,
                                            n_avg_yrs = 1,
                                            local_bh_msy_newton_steps = 6,
                                            is_discard_fleet = NULL,
                                            sd_rep = NULL,
                                            method = "delta",
                                            rel_step = 1e-3,
                                            min_step = 1e-7,
                                            n_draw = 300,
                                            seed = NULL,
                                            level = 0.95,
                                            extra_quantities = NULL,
                                            par_subset = NULL) {

  method <- match.arg(method, c("delta", "mvn", "both"))

  if(is.null(obj$data) || is.null(obj$report) || is.null(obj$env$last.par.best))
    stop("obj must be a fitted model object from fit_model, not a report list.", call. = FALSE)

  data <- obj$data
  if(is.null(sex_ratio_f)) sex_ratio_f <- array(0.5, dim = c(data$n_pop, data$n_regions))
  if(is.null(is_discard_fleet)) is_discard_fleet <- array(0, dim = data$n_fish_fleets)

  refpt_args <- list(
    SPR_x = SPR_x,
    t_spawn = t_spawn,
    sex_ratio_f = sex_ratio_f,
    calc_rec_st_yr = calc_rec_st_yr,
    rec_age = rec_age,
    type = type,
    what = what,
    n_avg_yrs = n_avg_yrs,
    local_bh_msy_newton_steps = local_bh_msy_newton_steps,
    is_discard_fleet = is_discard_fleet
  )

  p <- obj$env$last.par.best
  m <- length(p)
  cov_obj <- refpt_par_cov(obj, sd_rep)

  # steps scale with each parameter's own standard error, so every parameter is
  # probed over the range its uncertainty spans
  se_p <- if(cov_obj$kind == "cov") sqrt(pmax(diag(cov_obj$mat), 0)) else sqrt(pmax(Matrix::diag(Matrix::solve(cov_obj$mat)), 0))
  step <- pmax(rel_step * se_p, min_step)

  # base case w/ no perturbation
  base_raw <- eval_refpt_log_quantities(obj, p, refpt_args, extra_quantities)
  keep <- names(base_raw)[is.finite(base_raw)] # get values in the base case

  if(length(keep) == 0)
    stop("No reference point quantity was finite and positive at the estimates.", call. = FALSE)

  base <- base_raw[keep]

  idx <- seq_len(m)
  if(!is.null(par_subset)) {
    idx <- which(names(p) %in% par_subset)
    if(length(idx) == 0) stop("par_subset matched none of the fitted parameter names.", call. = FALSE)
  }

  # central differences of the solved reference point. The criterion curvature is
  # never formed, since the measured movement already has it divided out.
  d <- matrix(0, nrow = length(keep), ncol = m, dimnames = list(keep, names(p)))

  for(k in seq_along(idx)) {

    j <- idx[k]

    # subtract step
    pp <- p
    pp[j] <- pp[j] + step[j]

    # add step
    pm <- p
    pm[j] <- pm[j] - step[j]

    # what reference points give after perturbation
    up <- eval_refpt_log_quantities(obj, pp, refpt_args, extra_quantities, keep)
    dn <- eval_refpt_log_quantities(obj, pm, refpt_args, extra_quantities, keep)

    if(any(!is.finite(up)) || any(!is.finite(dn)))
      stop("Reference point solve failed while perturbing ", names(p)[j], ". Try a smaller rel_step.", call. = FALSE)

    # get sensitivity to perturbaiton (the jacobian)
    d[, j] <- (up - dn) / (2 * step[j])

  } # end k loop

  log_cov <- refpt_quad_form(cov_obj, d) # get covariance of reference point (weighted variance sum basically)
  dimnames(log_cov) <- list(keep, keep)
  log_se <- sqrt(pmax(diag(log_cov), 0))
  z <- stats::qnorm(1 - (1 - level) / 2) # get cis

  out <- list(
    refpts = NULL,
    mvn = NULL,
    d = d,
    log_cov = log_cov,
    draws = NULL,
    n_par = length(idx),
    n_zero_par = sum(apply(d, 2, function(col) all(col == 0))),
    step = step
  )

  if(method %in% c("delta", "both")) {
    out$refpts <- data.frame(quantity = keep,
                             est = as.numeric(exp(base)),
                             log_se = as.numeric(log_se),
                             lwr = as.numeric(exp(base - z * log_se)),
                             upr = as.numeric(exp(base + z * log_se)),
                             cv = as.numeric(sqrt(exp(log_se^2) - 1)),
                             row.names = NULL)
  }

  if(method %in% c("mvn", "both")) {

    if(!is.null(seed)) set.seed(seed)
    devs <- refpt_draw_devs(cov_obj, n_draw)

    draws <- matrix(NA_real_, nrow = n_draw, ncol = length(keep), dimnames = list(NULL, keep))
    for(b in seq_len(n_draw)) {
      draws[b, ] <- tryCatch(eval_refpt_log_quantities(obj, p + devs[b, ], refpt_args, extra_quantities, keep),
                             error = function(e) rep(NA_real_, length(keep)))
    } # end b loop

    # a draw far enough into the tail can fail to solve, so keep only the usable ones
    ok <- stats::complete.cases(draws)
    draws <- draws[ok, , drop = FALSE]
    out$draws <- draws
    out$mvn <- data.frame(
      quantity = keep,
      est = as.numeric(exp(base)),
      log_sd = apply(draws, 2, stats::sd),
      lwr = as.numeric(exp(apply(draws, 2, stats::quantile, probs = (1 - level) / 2))),
      upr = as.numeric(exp(apply(draws, 2, stats::quantile, probs = 1 - (1 - level) / 2))),
      n_draw = nrow(draws),
      row.names = NULL
    )

  } # end mvn draws

  out

}
