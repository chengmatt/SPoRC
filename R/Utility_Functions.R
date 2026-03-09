#' ggplot2 theme for SPoRC plots
#'
#' A clean \code{theme_bw}-based ggplot2 theme with enlarged text elements
#' suitable for publication-quality SPoRC diagnostic and results figures.
#'
#' @return A \code{ggplot2} theme object.
#'
#' @import ggplot2
#' @export theme_sablefish
#' @family Plotting
theme_sablefish <- function() {
   theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 17),
          title = element_text(size = 21, color = 'black'),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 17, color = 'black'),
          legend.text = element_text(size = 15, color = "black"),
          legend.title = element_text(size = 17, color = 'black'))
}


#' Construct an AR(1) correlation matrix
#'
#' Builds an \eqn{n \times n} correlation matrix whose \eqn{(i,j)} element
#' equals \eqn{\rho^{|i-j|}}, corresponding to a stationary AR(1) process
#' with autocorrelation parameter \eqn{\rho}.
#'
#' @param n Integer. Matrix dimension (number of bins, ages, or lengths).
#' @param rho Numeric. AR(1) autocorrelation parameter in \eqn{(-1, 1)}.
#'   Values close to 1 produce strong positive correlation between adjacent
#'   bins; values close to 0 approach the identity matrix.
#'
#' @return Numeric \eqn{n \times n} correlation matrix.
#'
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' get_AR1_CorrMat(10, 0.5)
#' }
get_AR1_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the correlation based on the lag distance
      corrMatrix[i, j] <- rho^(abs(i - j))
    } # end i
  } # end j
  return(corrMatrix)
}

#' Construct a constant (exchangeable) correlation matrix
#'
#' Builds an \eqn{n \times n} correlation matrix with 1 on the diagonal and
#' \eqn{\rho} on all off-diagonal elements, corresponding to a compound
#' symmetry (exchangeable) covariance structure. Used in SPoRC to model
#' constant correlation across sexes in the 2D logistic-normal composition
#' likelihood (\code{comp_like = 4}).
#'
#' @param n Integer. Matrix dimension (typically \code{n_sexes}).
#' @param rho Numeric. Off-diagonal correlation in \eqn{(-1, 1)}. A value
#'   of 0 produces the identity matrix; a value approaching 1 produces
#'   near-perfect correlation across sexes.
#'
#' @return Numeric \eqn{n \times n} correlation matrix.
#'
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' get_Constant_CorrMat(2, 0.5)
#' }
get_Constant_CorrMat <- function(n, rho) {
  corrMatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i != j) corrMatrix[i, j] <- rho
      else corrMatrix[i, j] <- 1
    } # end i
  } # end j
  return(corrMatrix)
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

#' Append a message to the global messages list
#'
#' Concatenates its arguments into a single string and appends the result to
#' \code{messages_list} in the calling environment via \code{<<-}. Used
#' internally to accumulate validation and setup messages for deferred display.
#'
#' @param ... Character strings passed to \code{paste(..., sep = "")}.
#'
#' @return \code{NULL} invisibly. Side effect: \code{messages_list} is
#'   updated in the parent environment.
#'
#' @keywords internal
collect_message <- function(...) {
  messages_list <<- c(messages_list, paste(..., sep = ""))
}

#' Check whether an R package is installed
#'
#' Tests for package availability without loading it, using
#' \code{system.file} to locate the package directory.
#'
#' @param pkg Character string. Name of the package to check.
#'
#' @return Logical. \code{TRUE} if the package is installed and findable;
#'   \code{FALSE} otherwise.
#'
#' @keywords internal
is_package_available <- function(pkg) {
  nzchar(system.file(package = pkg))
}

#' Convert a target catch to fishing mortality for a single fleet via bisection
#'
#' Uses interval bisection to find the scalar fishing mortality \eqn{F} that
#' produces a predicted catch (Baranov catch equation, summed over ages and
#' sexes in biomass) equal to \code{catch}. Intended for closed-loop MSE
#' harvest control rules where a TAC must be translated into an \eqn{F} for
#' the operating model.
#'
#' @param f_guess Numeric. Initial \eqn{F} guess (not used directly by the
#'   bisection algorithm but retained for API consistency).
#' @param catch Numeric. Target catch in biomass units.
#' @param NAA Numeric matrix \code{[n_ages × n_sexes]}. Numbers-at-age at
#'   the start of the time step.
#' @param WAA Numeric matrix \code{[n_ages × n_sexes]}. Weight-at-age.
#' @param natmort Numeric matrix \code{[n_ages × n_sexes]}. Instantaneous
#'   natural mortality rate.
#' @param fish_sel Numeric matrix \code{[n_ages × n_sexes]}. Fishery
#'   selectivity (scaled to a maximum of 1).
#' @param n.iter Integer. Number of bisection iterations. Default \code{20};
#'   approximately \eqn{\log_2((ub - lb) / \epsilon)} iterations are required
#'   for tolerance \eqn{\epsilon}.
#' @param lb Numeric. Lower bound of the \eqn{F} search interval.
#'   Default \code{0}.
#' @param ub Numeric. Upper bound of the \eqn{F} search interval.
#'   Default \code{2}.
#'
#' @return Scalar numeric. The \eqn{F} value at the final bisection midpoint
#'   that most closely produces \code{catch}.
#'
#'
#' @export catch_to_F_singlefleet
#' @family Closed Loop Simulations
catch_to_F_singlefleet <- function(f_guess,
                                    catch,
                                    NAA,
                                    WAA,
                                    natmort,
                                    fish_sel,
                                    n.iter = 20,
                                    lb = 0,
                                    ub = 2) {

  range <- vector(length=2) # F range
  range[1] <- lb # Lower bound
  range[2] <- ub # Upper bound

  for(i in 1:n.iter) {

    # Get midpoint of range
    midpoint <- mean(range)

    # Caclulate baranov's
    FAA <- (midpoint * fish_sel)
    ZAA <- FAA + natmort
    pred_catch <- sum((FAA / ZAA * NAA * (1 - exp(-ZAA))) * WAA)

    if(pred_catch < catch) {
      range[1] <- midpoint
      range[2] <- range[2]
    }else {
      range[1] <- range[1]
      range[2] <- midpoint
    }

  } # end i loop

  return(midpoint)
}

#' Convert target catches to fishing mortality rates for multiple fleets
#'
#' Solves for the vector of fleet-specific fishing mortality rates
#' \eqn{\mathbf{F} = (F_1, \ldots, F_k)} that simultaneously satisfy the
#' Baranov catch equations for all fleets, given a vector of target catches.
#' Total mortality at age accounts for contributions from all fleets:
#' \eqn{Z_a = M_a + \sum_f F_f s_{a,f}}. The system of equations is solved
#' via \code{nleqslv::nleqslv}. Intended for closed-loop MSE harvest control
#' rules with multiple interacting fishery fleets.
#'
#' @param target_catch Numeric vector \code{[n_fleets]}. Target catch in
#'   biomass units for each fleet.
#' @param NAA Numeric matrix \code{[n_ages × n_sexes]}. Numbers-at-age.
#' @param WAA Numeric matrix \code{[n_ages × n_sexes]}. Weight-at-age.
#' @param natmort Numeric matrix \code{[n_ages × n_sexes]}. Instantaneous
#'   natural mortality rate.
#' @param fish_sel Numeric array \code{[n_ages × n_sexes × n_fleets]}.
#'   Fishery selectivity for each fleet.
#' @param f_init Numeric scalar or vector \code{[n_fleets]}. Starting values
#'   for the \eqn{F} solver. If a scalar is supplied it is recycled across
#'   all fleets. Default \code{0.05}.
#' @param control Named list of control parameters passed to
#'   \code{nleqslv::nleqslv}. Default \code{list(btol = 1e-6)}.
#'
#' @return Numeric vector \code{[n_fleets]} of solved fishing mortality rates,
#'   one per fleet.
#'
#'
#' @export catch_to_F_multifleet
#' @family Closed Loop Simulations
catch_to_F_multifleet <- function(target_catch, NAA, WAA, natmort, fish_sel,
                               f_init = 0.05, control = list(btol = 1e-6)) {

  n_fleets <- length(target_catch)

  # Expand f_init if scalar
  if(length(f_init) == 1) f_init <- rep(f_init, n_fleets)

  # Function to minimize: difference between predicted and target catch for all fleets
  catch_diff <- function(f_vec) {
    pred_catches <- numeric(n_fleets)

    for(f in 1:n_fleets) {
      # F-at-age for this fleet
      FAA <- f_vec[f] * fish_sel[, , f]

      # Total Z includes F from ALL fleets
      ZAA_total <- natmort
      for(ff in 1:n_fleets) {
        ZAA_total <- ZAA_total + f_vec[ff] * fish_sel[, , ff]
      }

      # Predicted catch for this fleet (Baranov catch equation)
      pred_catches[f] <- sum((FAA / ZAA_total * NAA * (1 - exp(-ZAA_total))) * WAA)
    }

    return(pred_catches - target_catch)  # Difference from target
  }

  # Solve for F vector
  result <- nleqslv::nleqslv(f_init, catch_diff, control = control)

  return(result$x)
}

#' Run post-optimisation convergence checks on a fitted SPoRC model
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
#'   suggest the optimiser did not reach a local minimum.
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

#' Safely extract a named element from a TMB report object
#'
#' Returns the named element if it exists and is non-\code{NULL}; returns
#' \code{0} otherwise. Used to guard against missing or \code{NULL} report
#' fields when accumulating likelihood components.
#'
#' @param obj Named list, typically \code{obj$rep} from a fitted RTMB model.
#' @param name Character string. Name of the element to extract.
#'
#' @return The value of \code{obj[[name]]} if present and non-\code{NULL};
#'   \code{0} otherwise.
#'
#' @keywords internal
safe_extract <- function(obj, name) {
  if (name %in% names(obj) && !is.null(obj[[name]])) {
    return(obj[[name]])
  } else {
    return(0)
  }
}

#' Extract and tabulate parameter estimates and metadata from a fitted SPoRC model
#'
#' Joins the parameter list, factor map, and \code{sdreport} to produce two
#' tidy data frames: one for estimated parameters (with initial values,
#' posterior estimates, standard errors, and absolute gradients) and one for
#' fixed (non-estimated) parameters. Parameter indices are re-sequenced to
#' match the sequential numbering used internally by RTMB's \code{sdreport},
#' and map entries of \code{NA} (i.e., parameters fixed via \code{mapping})
#' are labelled \code{"NE"} (not estimated).
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
#' @family Utility
get_par_est_info <- function(parameters, mapping, sd_rep) {

  # get parameter names
  par_names <- reshape2::melt(parameters) %>%
    dplyr::rename_with(~str_replace(., "^Var(\\d+)$", "Dim\\1")) %>%
    dplyr::rename(Init_Val = value, Par = L1) %>%
    dplyr::group_by(Par) %>%
    # unique parameters based on dimensions of parameter list
    dplyr::mutate(Par_Num = paste(Par, row_number(), sep = "_"))

  # get mapping names
  map_names <- reshape2::melt(mapping) %>%
    dplyr::rename(map = value, Par = L1) %>%
    dplyr::group_by(Par) %>%
    dplyr::mutate(Par_Num = paste(Par, row_number(), sep = "_"),
                  map = as.character(map), # make character
                  map = ifelse(is.na(map), 'NE', map) # denote NAs as NE (for not estimated instead)
    )

  # join parameter names and mapping names
  par_map_names <- par_names %>%
    dplyr::left_join(map_names, by = c("Par", "Par_Num")) %>%
    dplyr::group_by(Par) %>%
    # make sure to turn NAs (they are estimated, but just were not in the mapping list)
    dplyr::mutate(map = ifelse(is.na(map), as.character(row_number()), map))

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
    dplyr::mutate(Par_Num_map = paste(Par, row_number(), sep = '_'))

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

#' Populate a parameter list with optimised values from sdreport
#'
#' Replaces starting values in \code{parameters} with the corresponding
#' optimised estimates from \code{sd_rep}, respecting the factor-map sharing
#' structure in \code{mapping}. For each parameter: if a map exists, factor
#' level integers are used to index into the estimated value vector so that
#' shared elements receive the same optimised value and \code{NA}-mapped
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
#'   their optimised values. Fixed (\code{NA}-mapped) elements retain their
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

#' Extend an array along its year dimension
#'
#' Appends additional year slices to an array along a specified dimension,
#' using one of several fill strategies. Used in SPoRC to extend biological,
#' selectivity, and sample-size arrays from the conditioning period into
#' projection years before running closed-loop MSE simulations.
#'
#' @param arr Array of any dimensionality to extend.
#' @param n_years Integer. Total number of years in the extended array
#'   (i.e., the new size of dimension \code{yr_dim}). Must be greater than
#'   \code{dim(arr)[yr_dim]}.
#' @param yr_dim Integer. Index of the dimension in \code{arr} corresponding
#'   to years.
#' @param fill Character string or numeric. Fill strategy for the appended
#'   year slices:
#'   \describe{
#'     \item{\code{"zeros"}}{Fill with zeros.}
#'     \item{\code{"last"}}{Repeat the last year slice that contains at least
#'       one non-\code{NA}, non-\code{NaN} value. If no valid slice exists,
#'       fills with \code{NA}.}
#'     \item{\code{"mean"}}{Fill with the per-element mean across years,
#'       excluding zeros, \code{NA}, and \code{NaN} values. Elements with no
#'       valid values are set to zero.}
#'     \item{\code{"F_pattern"}}{Fill with zeros; signals to downstream
#'       functions (e.g., \code{\link{predict_sim_fish_iss_fmort}}) that
#'       sample sizes should be dynamically updated based on projected fishing
#'       mortality during the closed-loop simulation.}
#'     \item{Numeric scalar or array}{Fill all appended slices with the
#'       supplied constant value, recycled via \code{array()}.}
#'   }
#'
#' @return Array with the same dimensions as \code{arr} except that
#'   \code{dim(result)[yr_dim] == dim(arr)[yr_dim] + n_years}, formed by
#'   binding \code{arr} and the fill array along \code{yr_dim} via
#'   \code{abind::abind}.
#'
#' @keywords internal
extend_years <- function(arr, n_years, yr_dim, fill = "zeros") {
  new_dims <- dim(arr); new_dims[yr_dim] <- n_years
  if(fill %in% c("zeros", "F_pattern")) {
    fill_array <- array(0, dim = new_dims)
  } else if(fill == "last") {
    # Get last non-NaN year slice along yr_dim
    # First, find the last year index that contains at least some non-NaN values
    last_valid_idx <- NULL
    for(i in dim(arr)[yr_dim]:1) {
      indices <- rep(list(quote(expr=)), length(dim(arr)))
      indices[[yr_dim]] <- i
      year_slice <- do.call(`[`, c(list(arr), indices, drop = FALSE))
      # check if this slice has any non-NaN values
      if(any(!is.na(year_slice) & !is.nan(year_slice))) {
        last_valid_idx <- i
        break
      }
    }
    # use NA if no valid year found
    if(is.null(last_valid_idx)) {
      fill_array <- array(NA, dim = new_dims)
    } else {
      # Get the last valid year slice
      indices <- rep(list(quote(expr=)), length(dim(arr)))
      indices[[yr_dim]] <- last_valid_idx
      last_year_slice <- do.call(`[`, c(list(arr), indices, drop = FALSE))

      # repeat slice n_years times
      fill_array <- array(0, dim = new_dims)
      for(i in 1:n_years) {
        fill_indices <- rep(list(quote(expr=)), length(dim(arr)))
        fill_indices[[yr_dim]] <- i
        fill_array <- do.call(`[<-`, c(list(fill_array), fill_indices, list(last_year_slice)))
      }
    }
  } else if (fill == "mean") {
    # get mean along the year dimension, excluding zeros and NaNs
    margins <- setdiff(seq_along(dim(arr)), yr_dim)
    # get mean excluding zeros and NaNs
    mean_slice <- apply(arr, margins, function(x) {
      valid_values <- x[!is.na(x) & !is.nan(x) & x != 0]
      if(length(valid_values) == 0) return(0) else mean(valid_values)
    })
    # extend mean_slice along year dimension
    fill_array <- array(mean_slice, dim = new_dims)
  } else if (is.numeric(fill)) {
    fill_array <- array(fill, dim = new_dims)
  }
  return(abind::abind(arr, fill_array, along = yr_dim))
}

#' Set data indicators to unused for specified years
#'
#' Sets \code{Use*} binary indicator arrays to \code{0} for the specified
#' years across one or more data types, and optionally removes conventional
#' tag cohorts released in those years. Used in MSE closed-loop simulations
#' and retrospective analyses to exclude future or withheld data from the
#' estimation model without modifying the underlying observation arrays.
#' Only years present in \code{data$years} are affected; out-of-range values
#' in \code{unused_years} are silently ignored.
#'
#' @param data Named list of model data as constructed by the
#'   \code{Setup_Mod_*} family of functions.
#' @param unused_years Integer vector. Year indices (relative to
#'   \code{data$years}) to mark as unused. Values not present in
#'   \code{1:length(data$years)} are dropped.
#' @param what Character vector. Data types to modify. Any combination of:
#'   \describe{
#'     \item{\code{"Catch"}}{Sets \code{UseCatch[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishIdx"}}{Sets \code{UseFishIdx[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishAgeComps"}}{Sets \code{UseFishAgeComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"FishLenComps"}}{Sets \code{UseFishLenComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"SrvIdx"}}{Sets \code{UseSrvIdx[, unused_years, , ] <- 0}.}
#'     \item{\code{"SrvAgeComps"}}{Sets \code{UseSrvAgeComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"SrvLenComps"}}{Sets \code{UseSrvLenComps[, unused_years, , ] <- 0}.}
#'     \item{\code{"conv_tagging"}}{Removes tag cohorts whose release year
#'       falls in \code{unused_years} from \code{conv_tagged_fish},
#'       \code{obs_conv_tag_fish_recap}, \code{conv_tag_release_indicator},
#'       and updates \code{n_conv_tag_cohorts}. Only applied when
#'       \code{any(data$use_conv_fish_tagging == 1)}.}
#'   }
#'
#' @return The modified \code{data} list with \code{Use*} indicators set to
#'   \code{0} for \code{unused_years} and, if applicable, tagging cohorts
#'   from those years removed.
#'
#' @export set_data_indicator_unused
#' @family Utility
set_data_indicator_unused <- function(data,
                                      unused_years,
                                      what = c('Catch', "FishIdx",
                                               "FishAgeComps", "FishLenComps",
                                               "SrvIdx", "SrvAgeComps", "SrvLenComps",
                                               "conv_tagging")) {

  # figure out year dimensions
  data_years <- 1:length(data$years)
  unused_years <- unused_years[which(unused_years %in% data_years)]

  if(length(unused_years) > 0) {
    # set to not use
    if("Catch" %in% what) data$UseCatch[,unused_years,,] <- 0
    if("FishIdx" %in% what) data$UseFishIdx[,unused_years,,] <- 0
    if("FishAgeComps" %in% what) data$UseFishAgeComps[,unused_years,,] <- 0
    if("FishLenComps" %in% what) data$UseFishLenComps[,unused_years,,] <- 0
    if("SrvIdx" %in% what) data$UseSrvIdx[,unused_years,,] <- 0
    if("SrvAgeComps" %in% what) data$UseSrvAgeComps[,unused_years,,] <- 0
    if("SrvLenComps" %in% what) data$UseSrvLenComps[,unused_years,,] <- 0
  }

  # modify tagging stuff
  if(any(data$use_conv_fish_tagging == 1) && "conv_tagging" %in% what) {
    tags_to_remove <- which(data$conv_tag_release_indicator[,2] %in% unused_years)
    if(length(tags_to_remove) > 0) {
      data$conv_tagged_fish <- data$Tagged_Fish[-tags_to_remove,,,,drop=FALSE]
      data$obs_conv_tag_fish_recap <- data$Obs_Tag_Recap[,,-tags_to_remove,,,,,,drop=FALSE]
      data$conv_tag_release_indicator <- data$conv_tag_release_indicator[-tags_to_remove,,drop=FALSE]
      data$n_conv_tag_cohorts <- nrow(data$tag_release_indicator)
    }
  }

  return(data)
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
#' @param adnuts_obj An \code{adnuts} posterior object containing
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
#'   rtmb_obj, adnuts_obj,
#'   what = c("SSB", "Rec"),
#'   n_cores = 4
#' )
#' }
get_model_rep_from_mcmc <- function(rtmb_obj, adnuts_obj, what, n_cores) {

  # discard warmup samples
  adnuts_obj$samples <- adnuts_obj$samples[-c(1:adnuts_obj$warmup),,]

  # define dimensions
  n_iter <- dim(adnuts_obj$samples)[1]
  n_chain <- dim(adnuts_obj$samples)[2]
  n_param <- dim(adnuts_obj$samples)[3]

  # collapse chains and iterations into a single posterior draw
  samples_collapsed <- matrix(adnuts_obj$samples, nrow = n_iter * n_chain, ncol = n_param)
  colnames(samples_collapsed) <- dimnames(adnuts_obj$samples)[[3]] # rename columns
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

#' Transform a real-valued parameter to the interval (-1, 1)
#'
#' Applies the scaled logistic transformation
#' \eqn{2 / (1 + e^{-2x}) - 1} to map an unconstrained real value to
#' \eqn{(-1, 1)}, suitable for parameterising correlation coefficients.
#' Used in SPoRC to back-transform raw correlation parameters before
#' constructing AR(1) and constant covariance matrices for logistic-normal
#' composition likelihoods.
#'
#' @param x Numeric. Unconstrained real-valued parameter.
#'
#' @return Numeric. Transformed value in \eqn{(-1, 1)}.
#'
#'
#' @export rho_trans
#' @family Utility
rho_trans <- function(x){
  2/(1+ exp(-2 * x)) - 1 # constraint between -1 and 1
}

#' Construct a logistic-normal covariance matrix
#'
#' Builds the covariance matrix \eqn{\Sigma} used in logistic-normal
#' composition likelihoods for a given correlation structure. Three structures
#' are supported, matching the \code{comp_like} codes used throughout SPoRC:
#' iid (\code{2}), AR(1) across bins (\code{3}), and AR(1) across bins with
#' constant correlation across sexes via a Kronecker product (\code{4}).
#'
#' @param comp_like Integer. Covariance structure: \code{2} = iid (diagonal
#'   \eqn{\theta^2 I}), \code{3} = AR(1) across bins
#'   (\eqn{\theta^2 C_{\text{AR1}}}), \code{4} = Kronecker product of
#'   constant sex correlation and AR(1) bin correlation
#'   (\eqn{\theta^2 (C_{\text{sex}} \otimes C_{\text{AR1}})}).
#' @param n_bins Integer. Number of composition categories (ages or lengths).
#'   The resulting matrix has dimension \code{n_bins} for \code{comp_like}
#'   \code{2} and \code{3}, or \code{n_bins × n_sexes} for \code{comp_like
#'   = 4}.
#' @param n_sexes Integer. Number of sexes. Required for \code{comp_like = 4};
#'   ignored otherwise.
#' @param theta Numeric. Marginal standard deviation \eqn{\theta > 0}
#'   controlling the overall scale of \eqn{\Sigma}.
#' @param corr_b Numeric. AR(1) correlation across bins in \eqn{(-1, 1)}.
#'   Required for \code{comp_like} \code{3} and \code{4}; ignored for
#'   \code{comp_like = 2}.
#' @param corr_s Numeric. Constant (exchangeable) correlation across sexes in
#'   \eqn{(-1, 1)}. Required for \code{comp_like = 4}; ignored otherwise.
#'
#' @return Numeric covariance matrix \eqn{\Sigma} of dimension
#'   \code{n_bins × n_bins} (\code{comp_like} \code{2}, \code{3}) or
#'   \code{(n_bins × n_sexes) × (n_bins × n_sexes)} (\code{comp_like = 4}).
#'
#'
#' @export get_logistN_Sigma
#' @family Utility
#'
#' @examples
#' \dontrun{
#' # iid
#' get_logistN_Sigma(comp_like = 2, n_bins = 5, n_sexes = NULL, theta = 0.5)
#'
#' # AR(1) across bins
#' get_logistN_Sigma(comp_like = 3, n_bins = 5, n_sexes = NULL,
#'                   theta = 0.5, corr_b = 0.3)
#'
#' # AR(1) across bins x constant across sexes
#' get_logistN_Sigma(comp_like = 4, n_bins = 5, n_sexes = 2,
#'                   theta = 0.5, corr_b = 0.3, corr_s = 0.2)
#' }
get_logistN_Sigma <- function(comp_like,
                              n_bins,
                              n_sexes,
                              theta,
                              corr_b = NULL,
                              corr_s = NULL
                              ) {

  # iid
  if(comp_like == 2) Sigma <- diag(rep(theta^2, n_bins))

  # 1dar1 across
  if(comp_like == 3) {
    # Construct Sigma matrix
    LN_corr_b <- corr_b # correlation by age / length
    Sigma <- get_AR1_CorrMat(n_bins, LN_corr_b)
    Sigma <- Sigma * theta^2
  }

  # 2dar1 across
  if(comp_like == 4) {
    # Construct Sigma matrix
    LN_corr_b <- corr_b
    LN_corr_s <- corr_s
    Sigma <- kronecker(get_Constant_CorrMat(n_sexes, LN_corr_s), get_AR1_CorrMat(n_bins, LN_corr_b))
    Sigma <- Sigma * theta^2
  }

  return(Sigma)
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
#' @param opt Named list of optimiser output. Must contain \code{$par} and
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

#' Convert character or numeric input to numeric codes
#'
#' Maps a character vector to integer codes via a named lookup list, or
#' passes numeric input through unchanged. Arrays and matrices are flattened,
#' converted element-wise, and restored to their original dimensions.
#' Unrecognised character values raise an informative error listing both the
#' invalid inputs and the valid options.
#'
#' @param x Character vector, numeric vector, or array to convert.
#' @param lookup Named list mapping valid character strings to numeric codes
#'   (e.g., \code{list("none" = 999, "multinomial" = 0)}).
#'
#' @return Numeric vector or array of the same shape as \code{x}.
#'
#' @keywords internal
convert_to_numeric <- function(x, lookup) {

  # Return numberic if already numeric
  if (is.numeric(x)) {
    return(x)
  }

  # if character, return numeric and convert
  if (is.character(x)) {
    result <- lookup[x]
    if (any(is.na(result))) {
      invalid <- x[is.na(lookup[x])]
      stop("Invalid character input: ", paste(invalid, collapse = ", "),
           "\nValid options: ", paste(names(lookup), collapse = ", "))
    }
    return(unlist(result))
  }

  # Handle arrays/matrices
  if (is.array(x)) {
    dims <- dim(x)
    result <- convert_to_numeric(as.vector(x), lookup)
    return(array(result, dim = dims))
  }

  stop("Input must be numeric or character")
}
