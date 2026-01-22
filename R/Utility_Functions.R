#' ggplot theme for sablefish
#'
#' @return ggplot theme
#' @export theme_sablefish
#' @import ggplot2
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


#' Function to fill in an n x n correlation AR(1) matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return correlation matrix for an ar1 process
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_AR1_CorrMat(10, 0.5)
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

#' Constant correlation matrix
#'
#' @param n Number of bins
#' @param rho correaltion parameter
#'
#' @return constant correlation matrix
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' mat <- get_Constant_CorrMat(10, 0.5)
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

#' For combining a parameter and data list in RTMB so a data object can be explicitly defined
#'
#' @param f Parameter list
#' @param d Data list
#' @keywords internal
cmb <- function(f, d) {
  function(p) f(p, d)
}

#' Helper function to collect messages
#'
#' @param ... character vector of messages
#' @keywords internal
collect_message <- function(...) {
  messages_list <<- c(messages_list, paste(..., sep = ""))
}

#' Helper function to check package availbility
#'
#' @param pkg package name character
#' @keywords internal
is_package_available <- function(pkg) {
  nzchar(system.file(package = pkg))
}

#' Go from TAC to Fishing Mortality using bisection for when a single fishery fleet exists
#'
#' @param f_guess Initial guess of F
#' @param catch Provided catch values
#' @param NAA Numbers, dimensioned by ages, and sexes
#' @param WAA Weight, dimensioned by ages and sexes
#' @param natmort Natural mortality dimensioned by ages and sex
#' @param fish_sel Fishery selectivity, dimesnioned by ages and sex
#' @param n.iter Number of iterations for bisection
#' @param lb Lower bound of F
#' @param ub Upper bound of F
#'
#' @returns Fishing mortality values for a single fleet
#' @export bisection_F
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

#' Solve for fishing mortality rates that achieve target catches for multiple fleets
#'
#' @param target_catch Numeric vector of target catch values for each fleet
#' @param NAA Matrix of numbers-at-age (ages x sexes)
#' @param WAA Matrix of weight-at-age (ages x sexes)
#' @param natmort Matrix of natural mortality (ages x sexes)
#' @param fish_sel 3D array of fishery selectivity (ages x sexes x fleets)
#' @param f_init Initial guess for F values (scalar or vector)
#' @param control List of control parameters for nleqslv
#' @return Numeric vector of F values for each fleet
#' @export solve_multifleet_F
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

#' Post Optimization Model Convergence Checks
#'
#' @param sd_rep sd report list from a `SPoRC` model
#' @param rep report list from a `SPoRC` model
#' @param gradient_tol Value for maximum gradient tolerance to use
#' @param se_tol Value for maximum standard error tolerance to use
#' @param corr_tol Value for maximum correlation tolerance to use
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

#' Helper function for extracting elements from TMB report
#'
#' @param obj TMB report object
#' @param name Name of object to be extracted
#' @keywords internal
safe_extract <- function(obj, name) {
  if (name %in% names(obj) && !is.null(obj[[name]])) {
    return(obj[[name]])
  } else {
    return(0)
  }
}

#' Helper function for extracting parameter information and names from TMB
#'
#' @param parameters Parameter list from setting up TMB object
#' @param mapping Mapping list from setting up TMB object
#' @param sd_rep SD Report from TMB obj
#'
#' @returns A list of dataframes for estimated and non-estimated parameter values.
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

#' Title Populates an unoptimized parameter list with optimized values
#'
#' @param parameters Parameter list
#' @param mapping Mapping list
#' @param sd_rep Sd report list from RTMB
#' @param random Character vector of whether random effects were estimated
#'
#' @returns Updated parameter list with associated optimized parameters corresponding to the same mapping indices
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
      unique_map_values <- unique(param_map[!is.na(param_map)]) # Get unique non-NA mapping values

      # Create estimated parameter array/vector by mapping unique values to estimates
      if(!is.vector(parameters[[param_name]])) {
        param_est <- array(NA, dim = dim(parameters[[param_name]]))
      } else {
        param_est <- rep(NA, length(parameters[[param_name]]))
      }

      # mapg unique values back to the parameter structure
      for (i in seq_along(unique_map_values)) {
        param_est[param_map == unique_map_values[i]] <- est_values[i]
      }

      # if mapping is NA use original, else use estimated (shared values)
      parameters[[param_name]][!is.na(param_map)] <- param_est[!is.na(param_map)]

      # No mapping - estimated by default
    } else if (param_name %in% est_param_names) {
      if(!is.vector(parameters[[param_name]])) parameters[[param_name]] <- array(sd_rep$par.fixed[est_param_names == param_name], dim = dim(parameters[[param_name]])) # not a vector
         else parameters[[param_name]] <- sd_rep$par.fixed[est_param_names == param_name] # vector
         }

  }

  return(parameters)
}

#' Extend an array along a year dimension
#'
#' This function takes an array and extends it along the specified year dimension.
#' The extension can either be filled with zeros or by repeating the last year slice.
#'
#' @param arr Array to extend. Can have any number of dimensions.
#' @param n_years Integer. The total number of years to extend the array to.
#' @param yr_dim Integer. The dimension of `arr` that corresponds to years.
#' @param fill Character or Numeric (scalar or array). How to fill the extended years:
#'   - `"zeros"`: fill with zeros
#'   - `"last"`: repeat the last year slice that is not a NaN or NA value. If all values are NA, then the array gets populated with NAs.
#'   - `"mean"`: take mean of time series
#'   - `"F_pattern"`: Used for fishery input sample sizes, dynamically fills in sample sizes based on fishing mortality values specified in a closed loop simulation
#'   -  Numeric: Constant scalar or array
#' @return An array extended along the `yr_dim` dimension.
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

#' Set Data Indicators to Unused for Specified Years
#'
#' @param data Data list for RTMB model
#' @param unused_years Integer vector specifying which years to mark as unused. Only years present in \code{data$years} are considered.
#' @param what Character vector specifying which data types to modify. Possible values include:
#'   \describe{
#'     \item{"Catch"}{Catch data indicators.}
#'     \item{"FishIdx"}{Fishery index data indicators.}
#'     \item{"FishAgeComps"}{Fishery age composition data indicators.}
#'     \item{"FishLenComps"}{Fishery length composition data indicators.}
#'     \item{"SrvIdx"}{Survey index data indicators.}
#'     \item{"SrvAgeComps"}{Survey age composition data indicators.}
#'     \item{"SrvLenComps"}{Survey length composition data indicators.}
#'     \item{"Tagging"}{Tagging data and associated cohorts.}
#'   }
#'
#' @returns The modified \code{data} object, with indicators set to 0 for the specified years and tagging cohorts removed if relevant.
#' @export set_data_indicator_unused
#' @family Utility
set_data_indicator_unused <- function(data,
                                      unused_years,
                                      what = c('Catch', "FishIdx",
                                               "FishAgeComps", "FishLenComps",
                                               "SrvIdx", "SrvAgeComps", "SrvLenComps",
                                               "Tagging")) {

  # figure out year dimensions
  data_years <- 1:length(data$years)
  unused_years <- unused_years[which(unused_years %in% data_years)]

  if(length(unused_years) > 0) {
    # set to not use
    if("Catch" %in% what) data$UseCatch[,unused_years,] <- 0
    if("FishIdx" %in% what) data$UseFishIdx[,unused_years,] <- 0
    if("FishAgeComps" %in% what) data$UseFishAgeComps[,unused_years,] <- 0
    if("FishLenComps" %in% what) data$UseFishLenComps[,unused_years,] <- 0
    if("SrvIdx" %in% what) data$UseSrvIdx[,unused_years,] <- 0
    if("SrvAgeComps" %in% what) data$UseSrvAgeComps[,unused_years,] <- 0
    if("SrvLenComps" %in% what) data$UseSrvLenComps[,unused_years,] <- 0
  }

  # modify tagging stuff
  if(data$UseTagging == 1 && "Tagging" %in% what) {
    tags_to_remove <- which(data$tag_release_indicator[,2] %in% unused_years)
    if(length(tags_to_remove) > 0) {
      data$Tagged_Fish <- data$Tagged_Fish[-tags_to_remove,,,drop=FALSE]
      data$Obs_Tag_Recap <- data$Obs_Tag_Recap[,-tags_to_remove,,,,drop=FALSE]
      data$tag_release_indicator <- data$tag_release_indicator[-tags_to_remove,,drop=FALSE]
      data$n_tag_cohorts <- nrow(data$tag_release_indicator)
    }
  }

  return(data)
}

#' Extract model report from MCMC posterior samples
#'
#' This function collapses MCMC chains from an RTMB/ADNUTS object,
#' generates model reports for each posterior draw, and extracts
#' specified components of the report.
#'
#' @param rtmb_obj An RTMB object created via `ADFun`.
#' @param adnuts_obj An `adnuts` object containing MCMC samples.
#' @param what Character vector specifying the names of components
#'   in the model report to extract.
#' @param n_cores Number of cores to use
#'
#' @return A named list of `data.table`s, one for each element in
#'   `what`. Each table contains the melted report component across
#'   all posterior samples, with an additional column
#'   `posterior_sample` indicating the MCMC draw index.
#' @family Model Diagnostics
#' @examples
#' \dontrun{
#' model_reports <- get_model_rep_from_mcmc(rtmb_obj, adnuts_obj,
#'                                          what = c("SSB", "Rec"))
#' }
#'
#' @export get_model_rep_from_mcmc
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

#' Title Constrains value between -1 and 1
#'
#' @param x Numeric value to constrain
#'
#' @returns Constrained value between -1 and 1
#' @export rho_trans
#' @family Utility
rho_trans <- function(x){
  2/(1+ exp(-2 * x)) - 1 # constraint between -1 and 1
}

#' Construct logistic-normal covariance matrix
#'
#' Helper function to generate the covariance matrix (\eqn{\Sigma}) used in
#' logistic-normal composition models. The structure depends on the
#' specification of \code{comp_like}:
#' \itemize{
#'   \item \code{comp_like = 2}: independent and identically distributed (iid)
#'     across categories (\eqn{n_\mathrm{categories}}).
#'   \item \code{comp_like = 3}: first-order autoregressive (AR1) correlation
#'     across categories (\eqn{n_\mathrm{categories}}).
#'   \item \code{comp_like = 4}: two-dimensional AR1 correlation across
#'     categories and sexes (\eqn{n_\mathrm{categories} \times n_\mathrm{sexes}}).
#' }
#'
#' @param comp_like Integer specifying the logistic-normal correlation structure:
#'   \itemize{
#'     \item 2 = iid across categories
#'     \item 3 = AR1 across categories
#'     \item 4 = AR1 across categories and sexes
#'   }
#' @param n_bins Number of composition categories (e.g., age or length bins).
#'   For \code{comp_like = 2, 3}, the covariance matrix is dimensioned
#'   \code{n_bins}. For \code{comp_like = 4}, it is dimensioned
#'   \code{n_bins * n_sexes}.
#' @param n_sexes Number of sexes. Required when \code{comp_like = 4}.
#' @param theta Standard deviation parameter controlling the overall scale
#'   of the covariance.
#' @param corr_b Correlation parameter across categories, in the interval
#'   \eqn{(-1, 1)}. Used when \code{comp_like = 3} or \code{4}.
#' @param corr_s Correlation parameter across sexes, in the interval
#'   \eqn{(-1, 1)}. Used when \code{comp_like = 4}.
#'
#' @return A covariance matrix \eqn{\Sigma} with dimension:
#'   \itemize{
#'     \item \code{n_bins} (\code{comp_like = 2, 3})
#'     \item \code{n_bins * n_sexes} (\code{comp_like = 4})
#'   }
#'
#' @export get_logistN_Sigma
#' @family Utility
#' @examples \dontrun{
#' n_cat <- 5
#' n_sexes <- 2
#'
#' # iid example (categories only)
#' get_logistN_Sigma(comp_like = 2, n_bins = n_cat, n_sexes = NULL, theta = 0.5)
#'
#' # AR1 across categories
#' get_logistN_Sigma(comp_like = 3, n_bins = n_cat, n_sexes = NULL, theta = 0.5,
#'                   corr_b = 0.3)
#'
#' # AR1 across categories and sexes
#' get_logistN_Sigma(comp_like = 4, n_bins = n_cat, n_sexes = n_sexes, theta = 0.5,
#'                   corr_b = 0.3, corr_s = 0.2)
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

#' Calculate the Corrected marginal AIC (AICc) from Optimization Results
#'
#' Computes the corrected marginal Akaike Information Criterion (AICc)
#' for model selection using optimization results. It supports objects returned
#' from different optimizers, such as `optim` or `nlminb`.
#'
#' @param opt A list containing optimization results. Must include either:
#'   \itemize{
#'     \item `"par"` and `"objective"` (e.g., from `optim`), or
#'     \item `"par"` and `"value"` (e.g., from `nlminb`)
#'   }
#' @param p Numeric. Penalty multiplier for the number of parameters. Default is 2.
#' @param n Numeric. Sample size. Default is `Inf`.
#'
#' @return Numeric. The corrected AIC (AICc) value.
#' @export marg_AIC
marg_AIC <- function(opt, p = 2, n = Inf){
  k <- length(opt[["par"]])
  if(all(c("par","objective") %in% names(opt))) negloglike <- opt[["objective"]]
  if(all(c("par","value") %in% names(opt))) negloglike <- opt[["value"]]
  Return <- p * k + 2 * negloglike + 2 * k * (k + 1) / (n - k - 1)
  return(Return)
}

