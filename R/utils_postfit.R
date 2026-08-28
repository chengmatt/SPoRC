# Shared helpers
#
# Helpers for working with a fitted object: convergence checks, parameter
# estimate tables, the parameter list for a refit, marginal AIC, and pulling a
# report out of MCMC draws.

#' Run post-optimization convergence checks on a fitted SPoRC model
#'
#' Evaluates four convergence criteria on the output of a fitted RTMB model:
#' (1) finiteness of the joint negative log-likelihood, (2) maximum absolute
#' gradient of fixed-effect parameters, (3) positive-definiteness of the
#' Hessian, and (4) finiteness and magnitude of parameter standard errors and
#' pairwise correlations. A diagnostic message is printed for each failed
#' check identifying the offending parameter or parameter pair.
#'
#' @param sd_rep \code{sdreport} object returned by \code{RTMB::sdreport}.
#'   Must contain \code{$gradient.fixed}, \code{$par.fixed},
#'   \code{$pdHess}, and \code{$cov.fixed}.
#' @param rep Report list returned by \code{obj$rep} after fitting a SPoRC
#'   model. Must contain \code{$jnLL}.
#' @param gradient_tol Numeric. Maximum tolerated absolute gradient for any
#'   fixed-effect parameter. Default \code{1e-3}; values above this threshold
#'   suggest the optimizer did not reach a local minimum.
#' @param se_tol Numeric. Maximum tolerated parameter standard error. Default
#'   \code{100}; very large SEs indicate poorly identified parameters.
#' @param corr_tol Numeric. Maximum tolerated absolute pairwise parameter
#'   correlation. Default \code{0.99}; values approaching 1 indicate
#'   near-redundant parameters.
#'
#' @return Logical. \code{TRUE} if all four checks pass; \code{FALSE} if any
#'   check fails. In either case, informative messages are printed via
#'   \code{message}.
#'
#' @export post_optim_sanity_checks
#' @family Utility
post_optim_sanity_checks <- function(sd_rep,
                                     rep,
                                     gradient_tol = 1e-3,
                                     se_tol = 100,
                                     corr_tol = 0.99
                                     ) {

  passed_post_sanity_checks <- TRUE

  # check likelihoods are all finite and not NA
  if(!all(is.finite(rep$jnLL))) {
    message("Found Inf in joint log-likelihood, model is not converged!")
    passed_post_sanity_checks <- F
  }

  # check maximum absolute gradients
  max_abs_grad_ndx <- which.max(abs(sd_rep$gradient.fixed))
  max_abs_grad <- abs(sd_rep$gradient.fixed)[max_abs_grad_ndx]
  if(gradient_tol < max_abs_grad) {
    message("Parameter: ", names(sd_rep$par.fixed)[max_abs_grad_ndx], " had absolute gradient = ", max_abs_grad,
            " which was greater than tolerance ", gradient_tol,". This indicates potential non-convergence according to the tolerance.\n")
    passed_post_sanity_checks <- F
  }

  # check hessian
  if(!sd_rep$pdHess) {
    message("Hessian is not positive definite, model is not converged!")
    passed_post_sanity_checks <- F
  }

  # check if standard errors are finite (if finite, then check other stuff)
  if(!all(is.finite(sqrt(diag(sd_rep$cov.fixed))))) {
    message("Found non finite elements in standard errors of parameters, model is not converged!")
    passed_post_sanity_checks <- F
  } else {
    # check if standard errors are big
    if(max(sqrt(diag(sd_rep$cov.fixed))) > se_tol) {
      message("Parameter: ", names(diag(sd_rep$cov.fixed))[which.max(sqrt(diag(sd_rep$cov.fixed)))], " has a standard error = ",
              max(sqrt(diag(sd_rep$cov.fixed))), " which was greated than tolerance ", se_tol, ". This indicates potential non-convergence according to the tolerance. \n")
      passed_post_sanity_checks <- F
    }

    # check if correlations are big
    corr_mat <- cov2cor(sd_rep$cov.fixed)
    diag(corr_mat) <- "Same" # set diagonal to "Same" to remove from max calculations

    # reshape to dataframe
    corr_df <- reshape2::melt(corr_mat) %>%
      dplyr::filter(value != 'Same') %>%
      dplyr::mutate(value = as.numeric(value))

    if(max(abs(corr_df$value)) > corr_tol) {
      message("Parameter pairs: ", corr_df$Var1[which.max(abs(corr_df$value))], " and ", corr_df$Var2[which.max(abs(corr_df$value))], " have a correlation of ", max(abs(corr_df$value)), ". This indicates potential non-convergence according to the tolerance.")
      passed_post_sanity_checks <- F
    }
  }

  if(passed_post_sanity_checks) {
    message("Successfully passed post-optim-sanity checks\n")
  }

  return(passed_post_sanity_checks)

}

#' Extract and tabulate parameter estimates and metadata from a fitted SPoRC model
#'
#' Joins the parameter list, factor map, and \code{sdreport} to produce two
#' tidy data frames: one for estimated parameters (with initial values,
#' posterior estimates, standard errors, and absolute gradients) and one for
#' fixed (non-estimated) parameters. Parameter indices are re-sequenced to
#' match the sequential numbering used internally by RTMB's \code{sdreport},
#' and map entries of \code{NA} (i.e., parameters fixed via \code{mapping})
#' are labeled \code{"NE"} (not estimated).
#'
#' @param parameters Named list of parameter starting values passed to
#'   \code{RTMB::MakeADFun}.
#' @param mapping Named list of factor maps passed to \code{RTMB::MakeADFun}.
#'   Parameters absent from \code{mapping} are treated as freely estimated.
#'   Elements with \code{NA} factor levels are treated as fixed.
#' @param sd_rep \code{sdreport} object returned by \code{RTMB::sdreport}.
#'   Must contain \code{$par.fixed}, \code{$par.random}, \code{$cov.fixed},
#'   \code{$diag.cov.random}, and \code{$gradient.fixed}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{est_pars}}{Data frame of estimated parameters, with columns
#'       \code{Par} (parameter name), \code{Est_Val} (estimated value on the
#'       native scale), \code{SE_Val} (standard error), \code{Abs_Grad_Val}
#'       (absolute gradient; \code{NA} for random effects), \code{Init_Val}
#'       (starting value), and dimension columns \code{Dim1}, \code{Dim2},
#'       \ldots indicating the array indices of each element.}
#'     \item{\code{non_est_pars}}{Data frame of non-estimated (fixed) parameters
#'       with the same dimension and \code{Init_Val} columns, plus \code{map =
#'       "NE"} indicating they were excluded from estimation.}
#'   }
#'
#' @export get_par_est_info
#' @importFrom stringr str_replace
#' @family Utility
get_par_est_info <- function(parameters, mapping, sd_rep) {

  # get parameter names
  par_names <- reshape2::melt(parameters) %>%
    dplyr::rename_with(~stringr::str_replace(., "^Var(\\d+)$", "Dim\\1")) %>%
    dplyr::rename(Init_Val = value, Par = L1) %>%
    dplyr::group_by(Par) %>%
    # unique parameters based on dimensions of parameter list
    dplyr::mutate(Par_Num = paste(Par, dplyr::row_number(), sep = "_"))

  # get mapping names
  map_names <- reshape2::melt(mapping) %>%
    dplyr::rename(map = value, Par = L1) %>%
    dplyr::group_by(Par) %>%
    dplyr::mutate(Par_Num = paste(Par, dplyr::row_number(), sep = "_"),
                  map = as.character(map), # make character
                  map = ifelse(is.na(map), 'NE', map) # denote NAs as NE (for not estimated instead)
    )

  # join parameter names and mapping names
  par_map_names <- par_names %>%
    dplyr::left_join(map_names, by = c("Par", "Par_Num")) %>%
    dplyr::group_by(Par) %>%
    # make sure to turn NAs (they are estimated, but just were not in the mapping list)
    dplyr::mutate(map = ifelse(is.na(map), as.character(dplyr::row_number()), map))

  # Make sure mapping numbers are sequential to match up with the sdreport
  par_map_names <- par_map_names %>%
    # filter(str_detect(Par, "ln_F_mean")) %>%
    dplyr::group_by(Par) %>%
    dplyr::mutate(map = ifelse(map != 'NE', as.character(cumsum(map != 'NE')), 'NE'),
                  Par_Num_map = paste(Par, map, sep = "_")) # now, make a variable that is consistent with numbering in sdreport

  # now, get estimater parameter names and values
  est_names <- data.frame(Par = c(names(sd_rep$par.fixed), names(sd_rep$par.random)),
                          Est_Val = c(sd_rep$par.fixed, sd_rep$par.random),
                          SE_Val = c(sqrt(diag(sd_rep$cov.fixed)), sd_rep$diag.cov.random),
                          Abs_Grad_Val = c(abs(as.vector(sd_rep$gradient.fixed)), rep(NA, length(sd_rep$diag.cov.random)))) %>%
    dplyr::group_by(Par) %>%
    dplyr::mutate(Par_Num_map = paste(Par, dplyr::row_number(), sep = '_'))

  # join to get estimated parameters, along with initial starting values
  estimated_pars <- est_names %>%
    dplyr::left_join(par_map_names, by = c("Par", "Par_Num_map")) %>%
    dplyr::select(-c(Par_Num))

  # also get non-estimated parameters
  non_estimated_pars <- par_map_names %>%
    dplyr::filter(map == 'NE')

  return(
    list(est_pars = estimated_pars,
         non_est_pars = non_estimated_pars)
  )
}

#' Populate a parameter list with optimized values from sdreport
#'
#' Replaces starting values in \code{parameters} with the corresponding
#' optimized estimates from \code{sd_rep}, respecting the factor-map sharing
#' structure in \code{mapping}. For each parameter: if a map exists, factor
#' level integers are used to index into the estimated value vector so that
#' shared elements receive the same optimized value and \code{NA}-mapped
#' (fixed) elements are left unchanged. Parameters absent from \code{mapping}
#' are treated as fully estimated and filled in sequentially. Random effects
#' are sourced from \code{sd_rep$par.random}; all other estimated parameters
#' from \code{sd_rep$par.fixed}.
#'
#' @param parameters Named list of parameter starting values passed to
#'   \code{RTMB::MakeADFun}.
#' @param mapping Named list of factor maps passed to \code{RTMB::MakeADFun}.
#' @param sd_rep \code{sdreport} object returned by \code{RTMB::sdreport}.
#'   Must contain \code{$par.fixed}, \code{$par.random}, and
#'   \code{$cov.fixed}.
#' @param random Character vector of parameter names declared as random
#'   effects in \code{RTMB::MakeADFun}. Used to route extraction to
#'   \code{sd_rep$par.random} rather than \code{sd_rep$par.fixed}.
#'
#' @return The \code{parameters} list with all estimated elements replaced by
#'   their optimized values. Fixed (\code{NA}-mapped) elements retain their
#'   original starting values.
#'
#' @keywords internal
get_optim_param_list <- function(parameters, mapping, sd_rep, random) {

  est_param_names <- names(c(sd_rep$par.fixed, sd_rep$par.random)) # get estimated parameter names

  for (param_name in names(parameters)) {
    map_name <- param_name
    if (map_name %in% names(mapping) && param_name %in% est_param_names) {
      # checking to see if vector
      if(!is.vector(parameters[[param_name]])) {
        param_map <- array(mapping[[map_name]], dim = dim(parameters[[param_name]])) # not a vector
      } else {
        param_map <- mapping[[map_name]] # vector
      }

      est_values <- if(param_name %in% random) sd_rep$par.random[names(sd_rep$par.random) == param_name] else sd_rep$par.fixed[names(sd_rep$par.fixed) == param_name] # Get estimated values for this parameter
      param_map_int <- as.integer(mapping[[map_name]])  # codes: 1=level1, 2=level2, etc.
      if (!is.vector(parameters[[param_name]])) param_map_int <- array(param_map_int, dim = dim(parameters[[param_name]]))
      non_na <- !is.na(param_map_int)
      parameters[[param_name]][non_na] <- est_values[param_map_int[non_na]]

      # No mapping - estimated by default
    } else if (param_name %in% est_param_names) {
      if(!is.vector(parameters[[param_name]])) parameters[[param_name]] <- array(sd_rep$par.fixed[est_param_names == param_name], dim = dim(parameters[[param_name]])) # not a vector
         else parameters[[param_name]] <- sd_rep$par.fixed[est_param_names == param_name] # vector
         }

  }

  return(parameters)
}

#' Extract model report quantities from MCMC posterior samples
#'
#' Discards warmup iterations, collapses all chains into a single matrix of
#' posterior draws, evaluates the RTMB report function at each draw in
#' parallel, and returns the requested report components as tidy
#' \code{data.table}s with a \code{posterior_sample} index column.
#'
#' @param rtmb_obj An RTMB \code{ADFun} object with a \code{$report()}
#'   method, as returned by \code{RTMB::MakeADFun}.
#' @param mcmc_obj An \code{adnuts} or \code{SparseNUTS} posterior object containing
#'   \code{$samples} (array \code{[n_iter × n_chain × n_param]}) and
#'   \code{$warmup} (number of warmup iterations to discard).
#' @param what Character vector. Names of components in the model report
#'   (i.e., quantities passed to \code{RTMB::REPORT} inside
#'   \code{\link{SPoRC_rtmb}}) to extract from each posterior draw.
#' @param n_cores Integer. Number of parallel workers to use via
#'   \code{future::multisession}.
#'
#' @return Named list of \code{data.table}s, one per element of \code{what}.
#'   Each table is the row-bound result of \code{reshape2::melt} applied to
#'   the report component across all post-warmup draws, with an additional
#'   integer column \code{posterior_sample} identifying the draw index
#'   (1 to \code{n_iter × n_chain}).
#'
#' @export get_model_rep_from_mcmc
#' @family Model Diagnostics
#'
#' @examples
#' \dontrun{
#' model_reports <- get_model_rep_from_mcmc(
#'   rtmb_obj, mcmc_obj,
#'   what = c("SSB", "Rec"),
#'   n_cores = 4
#' )
#' }
get_model_rep_from_mcmc <- function(rtmb_obj, mcmc_obj, what, n_cores) {

  # discard warmup samples
  mcmc_obj$samples <- mcmc_obj$samples[-c(1:mcmc_obj$warmup),,]

  # define dimensions
  n_iter <- dim(mcmc_obj$samples)[1]
  n_chain <- dim(mcmc_obj$samples)[2]
  n_param <- dim(mcmc_obj$samples)[3]

  # collapse chains and iterations into a single posterior draw
  samples_collapsed <- matrix(mcmc_obj$samples, nrow = n_iter * n_chain, ncol = n_param)
  colnames(samples_collapsed) <- dimnames(mcmc_obj$samples)[[3]] # rename columns
  what_list <- vector("list", length(what)) # create list to store
  names(what_list) <- what # name list

  future::plan(future::multisession, workers = n_cores)
  all_results <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(samples_collapsed)) # progress
    future.apply::future_lapply(1:nrow(samples_collapsed), function(idx) {
      tmp_rep <- rtmb_obj$report(par = samples_collapsed[idx, ])
      what_results <- vector("list", length(what)) # empty list
      names(what_results) <- what
      for (w in seq_along(what)) {
        tmp_what <- reshape2::melt(tmp_rep[[what[w]]]) # reshape2 in df
        tmp_what$posterior_sample <- idx # get posterior idx
        what_results[[w]] <- tmp_what # in[ut]
      } # end w loop
      p()  # upd prog
      what_results
    }, future.seed = TRUE)
  })

  # combine results
  what_list <- lapply(what, function(w) data.table::rbindlist(lapply(all_results, `[[`, w)))
  names(what_list) <- what # rename

  return(what_list)
}

#' Compute the corrected marginal AIC (AICc)
#'
#' Calculates the corrected Akaike Information Criterion
#' \eqn{\mathrm{AICc} = p k + 2 \ell + 2k(k+1)/(n-k-1)}, where \eqn{k} is
#' the number of estimated parameters, \eqn{\ell} is the negative
#' log-likelihood at the optimum, \eqn{p} is the penalty multiplier, and
#' \eqn{n} is the sample size. When \eqn{n = \infty} the small-sample
#' correction term vanishes and the result reduces to standard marginal AIC.
#' Compatible with output from both \code{optim} (\code{$value}) and
#' \code{nlminb} (\code{$objective}).
#'
#' @param opt Named list of optimizer output. Must contain \code{$par} and
#'   either \code{$objective} (e.g., \code{nlminb}) or \code{$value}
#'   (e.g., \code{optim}).
#' @param p Numeric. Penalty multiplier per parameter. Default \code{2}
#'   (standard AIC).
#' @param n Numeric. Sample size used for the small-sample correction.
#'   Default \code{Inf} (no correction, equivalent to AIC).
#'
#' @return Numeric. The AICc value.
#'
#' @export marg_AIC
marg_AIC <- function(opt, p = 2, n = Inf){
  k <- length(opt[["par"]])
  if(all(c("par","objective") %in% names(opt))) negloglike <- opt[["objective"]]
  if(all(c("par","value") %in% names(opt))) negloglike <- opt[["value"]]
  Return <- p * k + 2 * negloglike + 2 * k * (k + 1) / (n - k - 1)
  return(Return)
}
