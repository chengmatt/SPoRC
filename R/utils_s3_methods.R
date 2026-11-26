# Shared helpers
#
# Every S3 method in the package: print, summary, logLik, coef, vcov, confint, as.data.frame and plot for a
# fit (sporc_fit), and plot for retrospective results (sporc_retro) and OSA residuals (sporc_osa).

# Fitted model (sporc_fit) ------------------------------------------------

#' Print, summarize, compare and plot a fitted SPoRC model
#'
#' \code{\link{fit_model}} returns a list of class \code{sporc_fit}. Typing its
#' name prints a short description of the fit in place of the whole list,
#' \code{summary()} tabulates the estimates by parameter, and \code{logLik()}
#' returns the maximized log likelihood, so \code{AIC(fit_1, fit_2)} compares
#' fits. \code{coef()}, \code{vcov()} and \code{confint()} return the fixed
#' effect estimates, their covariance matrix and their confidence intervals, and
#' \code{as.data.frame()} returns one reported quantity as a long data frame.
#' \code{plot()} draws one of the package's plots for this fit alone; to
#' compare fits, call the plotting function itself with a list of them. The fit
#' is still a list: \code{fit$rep}, \code{fit$sd_rep} and
#' \code{RTMB::sdreport(fit)} work as before, and \code{unclass(fit)} prints
#' every element.
#'
#' \code{print} says whether \code{nlminb} converged and gives the largest
#' absolute gradient at the estimates, both stored with every optimized fit.
#' The Hessian check and the standard errors are read from \code{$sd_rep}, so
#' they are shown only for a fit run with \code{do_sdrep = TRUE}.
#'
#' @param x A fit returned by \code{\link{fit_model}}.
#' @param object A fit returned by \code{\link{fit_model}}.
#' @param digits Number of significant digits printed for the objective and
#'   the largest absolute gradient. Default \code{6}.
#' @param random For \code{coef}, \code{TRUE} returns the random effects with
#'   the fixed effects, in the order RTMB stores them. Default \code{FALSE}.
#' @param parm For \code{confint}, which estimates to return: parameter names
#'   (every estimated value of each, e.g. \code{"ln_RecDevs"}) or positions in
#'   \code{coef(object)}. Default is every fixed effect.
#' @param level Confidence level for \code{confint}. Default \code{0.95}.
#' @param row.names Not used.
#' @param optional Not used.
#' @param what For \code{as.data.frame}, name of the quantity in \code{x$rep}
#'   to return, e.g. \code{"SSB"} or \code{"NAA"}.
#' @param labels For \code{as.data.frame}, optional named list of labels for
#'   populations, regions, seasons, sexes and fleets, as in
#'   \code{\link{label_dims}}.
#' @param type Which plot to draw, each from the function named:
#'   \code{"ts"} (the default, \code{\link{get_ts_plot}}, which needs
#'   \code{do_sdrep = TRUE}), \code{"selex"} (\code{\link{get_selex_plot}}),
#'   \code{"biological"} (\code{\link{get_biological_plot}}),
#'   \code{"data_fitted"} (\code{\link{get_data_fitted_plot}}),
#'   \code{"nLL"} (\code{\link{get_nLL_plot}}), \code{"idx_fits"}
#'   (\code{\link{get_idx_fits_plot}}), \code{"catch_fits"}
#'   (\code{\link{get_catch_fits_plot}}) or \code{"at_age_fits"}
#'   (\code{\link{get_at_age_fits_plot}}).
#' @param model_name Name the fit is given in the plot legend. Default
#'   \code{"Model 1"}.
#' @param ... For \code{plot}, further arguments of the plotting function, such
#'   as \code{do_ci}, \code{Selex_Type}, \code{year_indx} or
#'   \code{data_source}. Not used by the other methods.
#'
#' @return \code{plot} returns what the plotting function returns, a
#'   \code{ggplot} or a list of them. \code{print} returns the fit invisibly. \code{summary} returns a
#'   data frame with one row per parameter name: \code{Par}, \code{Type}
#'   (\code{"fixed"} or \code{"random"}), \code{N} (number of estimated
#'   values), \code{Est_Min} and \code{Est_Max} (range of the estimates, on
#'   the scale the parameter is estimated on), \code{SE_Max} (largest standard
#'   error) and \code{Abs_Grad_Max} (largest absolute gradient, \code{NA} for
#'   random effects). \code{\link{get_par_est_info}} gives one row per
#'   estimated value. \code{logLik} returns minus the objective the fit
#'   minimized, with \code{df} the number of estimated fixed effects. That
#'   objective includes every penalty and prior in the model, and with
#'   \code{random} set it is the marginal. No number of observations is set,
#'   so \code{BIC()} is not available. \code{coef} returns the fixed effect
#'   estimates on the scale each is estimated on (mostly log), named by
#'   parameter, so a name repeats once per estimated value. \code{vcov} returns
#'   their covariance matrix from the sdreport, and \code{confint} a matrix of
#'   lower and upper limits on that same scale, the estimate plus and minus a
#'   normal quantile times its standard error. Both need
#'   \code{do_sdrep = TRUE} and leave out random effects, whose standard errors
#'   are in \code{fit$sd_rep}. \code{as.data.frame} returns one row per value
#'   of the quantity, with a column for each dim (\code{Region}, \code{Year},
#'   \code{Age}, ...) and the value in \code{value}. A fit run with
#'   \code{do_optim = FALSE} has no estimates, and \code{summary},
#'   \code{logLik}, \code{coef} and \code{confint} stop on it.
#'
#' @name sporc_fit
#' @family Utility
#'
#' @examples
#' \dontrun{
#' fit <- fit_model(data, parameters, mapping, do_sdrep = TRUE)
#' fit
#' summary(fit)
#' AIC(fit)
#' coef(fit)
#' confint(fit, parm = "ln_global_R0")
#' ssb_df <- as.data.frame(fit, what = "SSB") # columns Pop, Region, Year, value
#' plot(fit, type = "selex")
#' }
#'
#' @rdname sporc_fit
#' @export
print.sporc_fit <- function(x,
                            digits = 6,
                            ...
                            ) {

  # package version and time of the fit
  cat("SPoRC fit (version ", x$provenance$package_version, ", fitted ", format(x$provenance$fit_time, "%Y-%m-%d %H:%M"), ")\n", sep = "")

  # model dims, which a model of the user's own does not have
  if(!is.null(x$data$years)) {

    cat("  years:         ", min(x$data$years), "to", max(x$data$years), paste0("(", length(x$data$years), ")"), "\n")

    # c() drops a dim an older data list does not have
    dim_values <- c(n_pop = x$data$n_pop,
                    n_regions = x$data$n_regions,
                    n_seas = x$data$n_seas,
                    n_ages = length(x$data$ages),
                    n_lens = length(x$data$lens),
                    n_sexes = x$data$n_sexes,
                    n_fish_fleets = x$data$n_fish_fleets,
                    n_srv_fleets = x$data$n_srv_fleets)
    cat("  dimensions:    ", paste0(names(dim_values), " = ", dim_values, collapse = ", "), "\n")

  } # end if SPoRC data list

  # number of fixed effects, and the random effects by name
  n_random <- length(x$env$random)
  random_text <- "no random effects"
  if(n_random > 0) random_text <- paste0(n_random, " random effects (", paste(x$random, collapse = ", "), ")")
  cat("  parameters:    ", length(x$par), "fixed effects,", random_text, "\n")

  if(is.null(x$optim)) {

    # nothing was estimated, so the report is at the starting values
    cat("  optimizer:      not run (do_optim = FALSE)\n")
    if(!is.null(x$rep$jnLL)) cat("  jnLL:          ", format(x$rep$jnLL, digits = digits), "at the starting values\n")

  } else {

    # nlminb returns code 0 when it stopped at a minimum, anything else when it gave up
    converged_text <- "not converged"
    if(x$optim$convergence == 0) converged_text <- "converged"
    cat("  optimizer:     ", paste0(converged_text, " (nlminb code ", x$optim$convergence, ": ", x$optim$message, "), ", x$optim$iterations, " iterations"), "\n")

    # with random effects the fit minimized the marginal, not jnLL
    if(is.null(x$random)) cat("  jnLL:          ", format(x$optim$objective, digits = digits), "\n")
    else cat("  marginal nLL:  ", format(x$optim$objective, digits = digits), "\n")

    # largest absolute gradient at the estimates and the parameter it belongs to
    abs_grad <- abs(x$optim$gradient)
    max_grad_idx <- which.max(abs_grad)
    cat("  max gradient:  ", format(abs_grad[max_grad_idx], digits = digits), paste0("(", names(x$par)[max_grad_idx], ")"), "\n") # a newton step drops the names on optim$par

  } # end if else optimized

  if(is.null(x$sd_rep)) {

    cat("  Hessian:        not checked (do_sdrep = FALSE)\n")

  } else {

    if(x$sd_rep$pdHess) cat("  Hessian:        positive definite\n")
    else cat("  Hessian:        not positive definite\n")

    # result of post_optim_sanity_checks, which fit_model runs with the sdreport
    if(isTRUE(x$post_optim_sanity_checks)) cat("  sanity checks:  passed\n")
    if(isFALSE(x$post_optim_sanity_checks)) cat("  sanity checks:  failed\n")

  } # end if else sdreport

  return(invisible(x))
} # end function

#' @rdname sporc_fit
#' @export
summary.sporc_fit <- function(object,
                              ...
                              ) {

  if(is.null(object$optim)) stop("This fit was run with do_optim = FALSE, so there are no estimates to summarize.")

  # every estimated value at the optimum, fixed and random effects in the order RTMB stores them
  est_val <- object$env$last.par.best
  random_idx <- object$env$random # where the random effects sit in est_val, NULL without any
  fixed_idx <- setdiff(seq_along(est_val), random_idx)

  par_type <- rep("fixed", length(est_val))
  par_type[random_idx] <- "random"

  # gradient at the estimates, which random effects do not have
  abs_grad_val <- rep(NA, length(est_val))
  abs_grad_val[fixed_idx] <- abs(object$optim$gradient)

  # standard errors need the sdreport
  se_val <- rep(NA, length(est_val))
  if(!is.null(object$sd_rep)) {
    se_val[fixed_idx] <- sqrt(diag(object$sd_rep$cov.fixed))
    if(length(random_idx) > 0) se_val[random_idx] <- sqrt(object$sd_rep$diag.cov.random) # diag.cov.random is a variance
  } # end if sdreport

  # one row per parameter name
  par_names <- unique(names(est_val))
  par_table <- data.frame(Par = par_names,
                          Type = NA,
                          N = NA,
                          Est_Min = NA,
                          Est_Max = NA,
                          SE_Max = NA,
                          Abs_Grad_Max = NA)

  for(i in seq_along(par_names)) {

    par_idx <- which(names(est_val) == par_names[i]) # every estimated value of this parameter

    par_table$Type[i] <- par_type[par_idx[1]]
    par_table$N[i] <- length(par_idx)
    par_table$Est_Min[i] <- min(est_val[par_idx])
    par_table$Est_Max[i] <- max(est_val[par_idx])

    # left as NA where there is no standard error or gradient to report
    if(!all(is.na(se_val[par_idx]))) par_table$SE_Max[i] <- max(se_val[par_idx], na.rm = TRUE)
    if(!all(is.na(abs_grad_val[par_idx]))) par_table$Abs_Grad_Max[i] <- max(abs_grad_val[par_idx], na.rm = TRUE)

  } # end i loop

  return(par_table)
} # end function

#' @rdname sporc_fit
#' @importFrom stats logLik
#' @export
logLik.sporc_fit <- function(object,
                             ...
                             ) {

  if(is.null(object$optim)) stop("This fit was run with do_optim = FALSE, so it has no maximized likelihood.")

  # minus the minimized objective, the marginal when random effects were integrated out
  loglik_val <- -object$optim$objective
  attr(loglik_val, "df") <- length(object$optim$par) # estimated fixed effects
  class(loglik_val) <- "logLik"

  return(loglik_val)
} # end function

#' @rdname sporc_fit
#' @importFrom stats coef
#' @export
coef.sporc_fit <- function(object,
                           random = FALSE,
                           ...
                           ) {

  if(is.null(object$optim)) stop("This fit was run with do_optim = FALSE, so there are no estimates.")

  # fixed effects at the estimates
  coef_val <- object$optim$par
  names(coef_val) <- names(object$par) # a newton step drops the names on optim$par

  # fixed and random effects together, in the order RTMB stores them
  if(random) coef_val <- object$env$last.par.best

  return(coef_val)
} # end function

#' @rdname sporc_fit
#' @importFrom stats vcov
#' @export
vcov.sporc_fit <- function(object,
                           ...
                           ) {

  if(is.null(object$sd_rep)) stop("vcov needs the sdreport (object$sd_rep). Fit with do_sdrep = TRUE.")

  return(object$sd_rep$cov.fixed) # fixed effects only, rows and columns named by parameter
} # end function

#' @rdname sporc_fit
#' @importFrom stats confint
#' @export
confint.sporc_fit <- function(object,
                              parm,
                              level = 0.95,
                              ...
                              ) {

  est_val <- coef(object) # fixed effects, on the scale they are estimated on
  se_val <- sqrt(diag(vcov(object)))

  # every estimated value, unless parm picks some
  par_idx <- seq_along(est_val)

  if(!missing(parm)) {

    if(is.character(parm)) {

      missing_names <- setdiff(parm, names(est_val))
      if(length(missing_names) > 0) stop("parm names parameters that were not estimated as fixed effects: ", paste(missing_names, collapse = ", "), ". See unique(names(coef(object))).")
      par_idx <- which(names(est_val) %in% parm) # every value of a named parameter

    } else {

      if(!all(parm %in% seq_along(est_val))) stop("parm must be parameter names, or positions in coef(object) from 1 to ", length(est_val), ".")
      par_idx <- parm

    } # end if else names or positions

  } # end if parm given

  # normal quantiles for the lower and upper limit
  tail_prob <- (1 - level) / 2
  lower_val <- est_val[par_idx] + stats::qnorm(tail_prob) * se_val[par_idx]
  upper_val <- est_val[par_idx] + stats::qnorm(1 - tail_prob) * se_val[par_idx]

  # one row per estimated value, columns labeled the way stats::confint labels them
  ci_table <- cbind(lower_val, upper_val)
  ci_labels <- paste(format(100 * c(tail_prob, 1 - tail_prob), trim = TRUE, scientific = FALSE, digits = 3), "%")
  dimnames(ci_table) <- list(names(est_val)[par_idx], ci_labels)

  return(ci_table)
} # end function

#' @rdname sporc_fit
#' @export
as.data.frame.sporc_fit <- function(x,
                                    row.names = NULL,
                                    optional = FALSE,
                                    what,
                                    labels = list(),
                                    ...
                                    ) {

  if(missing(what) || length(what) != 1 || !what %in% names(x$rep)) stop("what must name one quantity in the report, e.g. what = 'SSB'. See names(x$rep).")

  rep_val <- x$rep[what] # kept as a list of one, which label_dims reads

  # a vector has no dim, so give it one for melt to make a column of (e.g. Year)
  if(is.null(dim(rep_val[[what]])) && length(rep_val[[what]]) > 1) rep_val[[what]] <- array(data = rep_val[[what]], dim = length(rep_val[[what]]))

  # name the dims, which needs the years and ages of a SPoRC data list
  if(!is.null(x$data$years) && !is.null(x$data$ages)) rep_val <- label_dims(x = rep_val, data = x$data, labels = labels)

  # one row per value, one column per dim
  rep_df <- reshape2::melt(rep_val[[what]])

  return(rep_df)
} # end function

#' @rdname sporc_fit
#' @export
plot.sporc_fit <- function(x,
                           type = "ts",
                           model_name = "Model 1",
                           ...
                           ) {

  if(!type %in% c("ts", "selex", "biological", "data_fitted", "nLL", "idx_fits", "catch_fits", "at_age_fits"))
    stop("type must be one of: ts, selex, biological, data_fitted, nLL, idx_fits, catch_fits, at_age_fits")

  # the time series plot reads its standard errors from the sdreport
  if(type == "ts" && is.null(x$sd_rep)) stop("x needs the sdreport (x$sd_rep) for type = 'ts'. Fit with do_sdrep = TRUE, or choose another type.")

  # the plotting functions compare models, so each takes a list with one element per model
  data <- list(x$data)
  rep <- list(x$rep)
  sd_rep <- list(x$sd_rep)

  if(type == "ts") out <- get_ts_plot(rep = rep, sd_rep = sd_rep, model_names = model_name, ...)
  if(type == "selex") out <- get_selex_plot(rep = rep, model_names = model_name, ...)
  if(type == "biological") out <- get_biological_plot(data = data, rep = rep, model_names = model_name)
  if(type == "data_fitted") out <- get_data_fitted_plot(data = data, model_names = model_name)
  if(type == "nLL") out <- get_nLL_plot(data = data, rep = rep, model_names = model_name)
  if(type == "idx_fits") out <- get_idx_fits_plot(data = data, rep = rep, model_names = model_name)
  if(type == "catch_fits") out <- get_catch_fits_plot(data = data, rep = rep, model_names = model_name)
  if(type == "at_age_fits") out <- get_at_age_fits_plot(data = data, rep = rep, model_names = model_name, ...)

  return(out)
} # end function

# Retrospective results (sporc_retro) -------------------------------------

#' @rdname get_retrospective_plot
#' @param x Retrospective results from \code{\link{do_retrospective}}, the
#'   data frame itself or \code{$retro_df} under \code{return_models = TRUE}.
#'   \code{plot(x, Rec_Age)} draws the same three plots.
#' @param ... Not used.
#' @export
plot.sporc_retro <- function(x,
                             Rec_Age,
                             ...
                             ) {

  return(get_retrospective_plot(retro_output = x, Rec_Age = Rec_Age))
} # end function

# OSA residuals (sporc_osa) -----------------------------------------------

#' @rdname plot_resids
#' @param x Residuals from \code{\link{get_osa}}. \code{plot(x)} draws the same
#'   two plots.
#' @param ... Not used.
#' @export
plot.sporc_osa <- function(x,
                           ...
                           ) {

  return(plot_resids(osa_results = x))
} # end function
