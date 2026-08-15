# Stage 1 of 3: model setup
#
# Catch and fishing mortality inputs, plus discard mortality rate.
# Setup_Mod_Catch_and_F is the entry point and calls every do_*_mapping helper in
# this file to build the parameter maps for F, its deviations, the catch and
# discard observation error terms, and the discard mortality rate.
# Setup_Sim_Fishing, the operating model counterpart, lives in setup_sim_fleets.R.

#' Map sigmaF (fishing mortality process error SD) parameters
#'
#' Constructs the \code{ln_sigmaF} factor map used by the TMB/RTMB objective
#' function to share or fix the log-scale standard deviation of fishing mortality
#' process error across regions, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigmaF_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigmaF}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × season × fleet combination.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per region × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per region × season.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons; unique per fleet.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets; unique per season.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets; unique per region.}
#'     \item{\code{"est_shared_r_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigmaF} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigmaF} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigmaF))}. Each element is an integer
#'   estimation index for shared or estimated configurations, or \code{NA} when
#'   \code{sigmaF_spec = "fix"}.
#'
#' @keywords internal
do_sigmaF_mapping <- function(input_list, sigmaF_spec) {

  # Sigma F -----------------------------------------------------------------
  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaF <- build_shared_spec_map(
    dims = dims, spec = sigmaF_spec,
    dim_abbrev = c(r = "region", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigmaF is specified as: ", sigmaF_spec)

  return(input_list)
}

#' Map AR1 correlation parameter for fishing mortality deviations
#'
#' Constructs the \code{Fdev_rho} factor map. \code{Fdev_rho} is only
#' meaningful when \code{Fdev_model = "ar1"} (see
#' \code{\link{Setup_Mod_Catch_and_F}}); for any other \code{Fdev_model}, all
#' \code{Fdev_rho} parameters are mapped to \code{NA} regardless of
#' \code{Fdev_rho_spec}, since they are unused by \code{\link{Get_Fdev_PE_loglik}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param Fdev_rho_spec Character string controlling the sharing and
#'   estimation structure for \code{Fdev_rho}, following the same convention
#'   as \code{\link{do_sigmaF_mapping}}'s \code{sigmaF_spec}: one of
#'   \code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_seas"},
#'   \code{"est_shared_f"}, \code{"est_shared_r_seas"}, \code{"est_shared_r_f"},
#'   \code{"est_shared_seas_f"}, \code{"est_shared_r_seas_f"}, or \code{"fix"}.
#'
#' @return The input \code{input_list} with \code{$map$Fdev_rho} set to a
#'   factor vector of length \code{prod(dim(par$Fdev_rho))}.
#'
#' @keywords internal
do_Fdev_rho_mapping <- function(input_list, Fdev_rho_spec) {

  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  if(input_list$data$Fdev_model != 3) { # only AR1 uses Fdev_rho
    input_list$map$Fdev_rho <- factor(rep(NA, prod(dims)))
  } else {
    input_list$map$Fdev_rho <- build_shared_spec_map(
      dims = dims, spec = Fdev_rho_spec,
      dim_abbrev = c(r = "region", seas = "season", f = "fleet")
    )
  }

  collect_message("Fdev_rho is specified as: ", Fdev_rho_spec)

  return(input_list)
}

#' Map sigma_C (catch observation error SD) parameters
#'
#' Constructs the \code{ln_sigmaC} factor map used by the TMB/RTMB objective
#' function to share or fix the log-scale standard deviation of catch observation
#' error across regions, years, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigmaC_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigmaC}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions.}
#'     \item{\code{"est_shared_y"}}{Shared across years.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years, and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years, and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons, and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigmaC} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigmaC} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigmaC))}. Each element is an integer
#'   estimation index, or \code{NA} when \code{sigmaC_spec = "fix"}.
#'
#' @keywords internal
do_sigmaC_mapping <- function(input_list, sigmaC_spec) {

  # Sigma C -----------------------------------------------------------------
  dims <- c(region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaC <- build_shared_spec_map(
    dims = dims, spec = sigmaC_spec,
    dim_abbrev = c(r = "region", y = "year", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigmaC is specified as: ", sigmaC_spec)

  return(input_list)
}

#' Map fishing mortality parameters
#'
#' Constructs the \code{ln_F_devs} and \code{ln_F_mean} factor maps, assigning
#' unique estimation indices to cells where catch data are used
#' (\code{UseCatch == 1}) and mapping cells without catch data to \code{NA}.
#' This ensures that fishing mortality parameters are only estimated for
#' dimensions with observed catch. \code{ln_F_devs} is resolved per
#' region-year-season-fleet cell, while \code{ln_F_mean} is resolved per
#' region-season-fleet cell and is estimated whenever that cell is fished in
#' at least one year.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$UseCatch} and \code{$data$UseCatch_pop} to be populated by
#'   \code{\link{Setup_Mod_Catch_and_F}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_F_devs} and
#'   \code{$map$ln_F_mean} set to factor vectors. Cells with catch are assigned
#'   sequential integer indices; cells without catch are \code{NA}.
#'
#' @keywords internal
do_Fmort_mapping <- function(input_list) {

  dims <- c(region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  # Estimate F devs if aggregated catch OR any pop-specific catch is used, or
  # if the aggregate catch observation is missing (NA) rather than a true
  # recorded zero; fishing is assumed to have continued through a missing
  # observation, whereas a recorded zero (or no catch data used at all) with
  # no missing observation indicates a true closure
  has_catch <- input_list$data$UseCatch == 1 |
    apply(input_list$data$UseCatch_pop == 1, c(2,3,4,5), any) |
    is.na(input_list$data$ObsCatch)

  F_dev_map <- build_pe_map(dims, share_over = character(0))
  F_dev_map[!has_catch] <- NA

  input_list$map$ln_F_devs <- factor(as.vector(F_dev_map))

  # ln_F_mean only enters the objective through cells that are fished in at
  # least one year, the same condition that keeps Fmort free rather than
  # pinned to zero. A region-season-fleet cell that is closed in every year
  # (e.g. a fleet that only operates in its own region/season) contributes
  # nothing to the likelihood, so estimating it leaves a flat gradient
  F_mean_dims <- dims[c("region", "season", "fleet")]
  F_mean_active <- apply(has_catch, c(1,3,4), any)

  F_mean_map <- build_pe_map(F_mean_dims, share_over = character(0))
  F_mean_map[!F_mean_active] <- NA

  input_list$map$ln_F_mean <- factor(as.vector(F_mean_map))

  # mirror the deviation map into the data list so the process error penalty can
  # key on the cells that are actually estimated. A deviation mapped off by hand
  # after setup is then neither estimated nor penalized
  input_list$data$map_ln_F_devs <- array(as.numeric(input_list$map$ln_F_devs), dim = dims)

  return(input_list)
}

#' Map population-specific catch observation error SD parameters
#'
#' Constructs the \code{ln_sigmaC_pop} factor map used by the TMB/RTMB
#' objective function to share or fix the log-scale standard deviation of
#' population-specific catch observation error across populations, regions,
#' years, seasons, and fleets. All cells within a shared group are assigned
#' the same estimation index.
#'
#' The sharing specification encodes which dimensions are collapsed into a
#' single parameter via underscore-separated tokens. For example,
#' \code{"est_shared_pop_r"} shares across populations and regions (one
#' parameter per year × season × fleet combination), while
#' \code{"est_shared_r_f"} shares across regions and fleets (one parameter
#' per population × year × season combination).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions. Requires
#'   \code{$data$n_pop}, \code{$data$n_regions}, \code{$data$years},
#'   \code{$data$n_seas}, \code{$data$n_fish_fleets}, and
#'   \code{$par$ln_sigmaC_pop} to be populated before calling.
#' @param sigmaC_pop_spec Character string controlling the sharing and
#'   estimation structure for \code{ln_sigmaC_pop}. One of:
#'   \describe{
#'     \item{\code{"fix"}}{All parameters fixed at starting values
#'       (mapped to \code{NA}).}
#'     \item{\code{"est_all"}}{Unique parameter per population × region ×
#'       year × season × fleet cell.}
#'     \item{\code{"est_shared_pop"}}{Shared across populations; unique per
#'       region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per
#'       population × year × season × fleet.}
#'     \item{\code{"est_shared_y"}}{Shared across years; unique per
#'       population × region × season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per
#'       population × region × year × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per
#'       population × region × year × season.}
#'     \item{\code{"est_shared_pop_r"}}{Shared across populations and regions.}
#'     \item{\code{"est_shared_pop_y"}}{Shared across populations and years.}
#'     \item{\code{"est_shared_pop_seas"}}{Shared across populations and seasons.}
#'     \item{\code{"est_shared_pop_f"}}{Shared across populations and fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_pop_r_y"}}{Shared across populations, regions,
#'       and years.}
#'     \item{\code{"est_shared_pop_r_seas"}}{Shared across populations, regions,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_r_f"}}{Shared across populations, regions,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_y_seas"}}{Shared across populations, years,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_y_f"}}{Shared across populations, years,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_seas_f"}}{Shared across populations, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years,
#'       and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years,
#'       and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas"}}{Shared across populations,
#'       regions, years, and seasons.}
#'     \item{\code{"est_shared_pop_r_y_f"}}{Shared across populations, regions,
#'       years, and fleets.}
#'     \item{\code{"est_shared_pop_r_seas_f"}}{Shared across populations,
#'       regions, seasons, and fleets.}
#'     \item{\code{"est_shared_pop_y_seas_f"}}{Shared across populations, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Shared across regions, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas_f"}}{Single parameter shared across
#'       all dimensions.}
#'   }
#'
#'
#' @keywords internal
do_sigmaC_pop_mapping <- function(input_list, sigmaC_pop_spec) {

  dims <- c(pop    = input_list$data$n_pop,
            region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaC_pop <- build_shared_spec_map(
    dims = dims, spec = sigmaC_pop_spec,
    dim_abbrev = c(pop = "pop", r = "region", y = "year", seas = "season", f = "fleet")
  )

  collect_message("sigmaC_pop is specified as: ", sigmaC_pop_spec)
  return(input_list)
}

#' Map sigma_dmr (discard mortality process error SD) parameters
#'
#' Constructs the \code{ln_sigma_dmr} factor map used by the TMB/RTMB objective
#' function to share or fix the logit-scale standard deviation of discard mortality
#' process error across regions, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigma_dmr_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigma_dmr}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × season × fleet combination.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per region × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per region × season.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons; unique per fleet.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets; unique per season.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets; unique per region.}
#'     \item{\code{"est_shared_r_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigma_dmr} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigma_dmr} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigma_dmr))}. Each element is an integer
#'   estimation index for shared or estimated configurations, or \code{NA} when
#'   \code{sigma_dmr_spec = "fix"}.
#'
#' @keywords internal
do_sigma_dmr_mapping <- function(input_list, sigma_dmr_spec) {

  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigma_dmr <- build_shared_spec_map(
    dims = dims, spec = sigma_dmr_spec,
    dim_abbrev = c(r = "region", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigma_dmr is specified as: ", sigma_dmr_spec)

  return(input_list)
}

#' Map sigmaD (discard mortality rate observation error SD) parameters
#'
#' Constructs the \code{ln_sigmaD} factor map used by the TMB/RTMB objective
#' function to share or fix the log-scale standard deviation of discard mrotality observation
#' error across regions, years, seasons, and fleets. All cells within a shared
#' group are assigned the same estimation index.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions.
#' @param sigmaD_spec Character string controlling the sharing and estimation
#'   structure for \code{ln_sigmaD}. One of:
#'   \describe{
#'     \item{\code{"est_all"}}{Unique parameter per region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions.}
#'     \item{\code{"est_shared_y"}}{Shared across years.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years, and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years, and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons, and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Single parameter shared across all dimensions.}
#'     \item{\code{"fix"}}{All \code{ln_sigmaD} parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_sigmaD} set to a factor
#'   vector of length \code{prod(dim(par$ln_sigmaD))}. Each element is an integer
#'   estimation index, or \code{NA} when \code{sigmaD_spec = "fix"}.
#'
#' @keywords internal
do_sigmaD_mapping <- function(input_list, sigmaD_spec) {

  dims <- c(region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaD <- build_shared_spec_map(
    dims = dims, spec = sigmaD_spec,
    dim_abbrev = c(r = "region", y = "year", seas = "season", f = "fleet")
  )

  # Print Message
  collect_message("sigmaD is specified as: ", sigmaD_spec)

  return(input_list)
}

#' Map discard mortality deviation parameters
#'
#' Constructs the \code{logit_dmr_devs} factor map, assigning unique estimation
#' indices to region-year-season-fleet cells that are fished and mapping true
#' closures to \code{NA}. A cell is fished when aggregated or any
#' population-specific catch is used, or when the aggregate catch observation is
#' missing (\code{NA}) rather than a recorded zero, which is the same condition
#' under which the objective computes a non-zero \code{dmr}.
#'
#' Discard observations are not required for a deviation to be estimable.
#' \code{dmr} enters the likelihood through total mortality (\code{ZAA}), so it
#' is informed by retained catch, indices, and compositions in any fished cell
#' where retention is less than one. \code{\link{get_dmr_penalty}} keys on the
#' same condition, so the estimated and penalized sets coincide.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$UseCatch}, \code{$data$UseCatch_pop}, and
#'   \code{$data$ObsCatch} to be populated by
#'   \code{\link{Setup_Mod_Catch_and_F}}.
#' @param dmr_dev_spec Character string specifying whether to estimate or fix
#'   deviations. Currently supports \code{"est_all"} and \code{"fix"}.
#'
#' @return The input \code{input_list} with \code{$map$logit_dmr_devs} set to a
#'   factor vector. Fished cells are assigned unique integer indices; true
#'   closures are \code{NA}.
#'
#' @keywords internal
do_dmr_dev_mapping <- function(input_list, dmr_dev_spec) {

  valid_specs <- c("fix", "est_all")
  if(!dmr_dev_spec %in% valid_specs)
    stop("dmr_dev_spec '", dmr_dev_spec, "' not recognized. Valid options: ",
         paste(valid_specs, collapse = ", "))

  dims <- c(region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  if(dmr_dev_spec == "fix") {
    input_list$map$logit_dmr_devs <- factor(rep(NA, length(input_list$par$logit_dmr_devs)))
    input_list$data$map_logit_dmr_devs <- array(NA_real_, dim = dims)
    collect_message("dmr_devs is specified as: fix")
    return(input_list)
  }

  # Estimate a deviation wherever dmr enters the dynamics, which is every cell
  # the objective does not treat as a true closure: aggregated OR any
  # population-specific catch is used, or the aggregate catch observation is
  # missing (NA) rather than a recorded zero. dmr is identified through total
  # mortality (ZAA), so a cell is informative whenever it is fished and
  # retention is less than one, with or without discard observations.
  # get_dmr_penalty() keys on the resulting map, so every estimated deviation is
  # penalized and no deviation fixed at zero contributes to the penalty
  has_catch <- input_list$data$UseCatch == 1 |
    apply(input_list$data$UseCatch_pop == 1, c(2,3,4,5), any) |
    is.na(input_list$data$ObsCatch)

  dmr_dev_map <- build_pe_map(dims, share_over = character(0))
  dmr_dev_map[!has_catch] <- NA

  input_list$map$logit_dmr_devs <- factor(as.vector(dmr_dev_map))

  # mirror the map into the data list so a deviation mapped off by hand after
  # setup is neither estimated nor penalized
  input_list$data$map_logit_dmr_devs <- array(as.numeric(input_list$map$logit_dmr_devs), dim = dims)

  collect_message("dmr_devs is specified as: est_all")
  return(input_list)
}

#' Map population-specific discard observation error SD parameters
#'
#' Constructs the \code{ln_sigmaD_pop} factor map used by the TMB/RTMB
#' objective function to share or fix the log-scale standard deviation of
#' population-specific discard observation error across populations, regions,
#' years, seasons, and fleets. All cells within a shared group are assigned
#' the same estimation index.
#'
#' The sharing specification encodes which dimensions are collapsed into a
#' single parameter via underscore-separated tokens. For example,
#' \code{"est_shared_pop_r"} shares across populations and regions (one
#' parameter per year × season × fleet combination), while
#' \code{"est_shared_r_f"} shares across regions and fleets (one parameter
#' per population × year × season combination).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists, as constructed by upstream setup functions. Requires
#'   \code{$data$n_pop}, \code{$data$n_regions}, \code{$data$years},
#'   \code{$data$n_seas}, \code{$data$n_fish_fleets}, and
#'   \code{$par$ln_sigmaD_pop} to be populated before calling.
#' @param sigmaD_pop_spec Character string controlling the sharing and
#'   estimation structure for \code{ln_sigmaD_pop}. One of:
#'   \describe{
#'     \item{\code{"fix"}}{All parameters fixed at starting values
#'       (mapped to \code{NA}).}
#'     \item{\code{"est_all"}}{Unique parameter per population × region ×
#'       year × season × fleet cell.}
#'     \item{\code{"est_shared_pop"}}{Shared across populations; unique per
#'       region × year × season × fleet.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per
#'       population × year × season × fleet.}
#'     \item{\code{"est_shared_y"}}{Shared across years; unique per
#'       population × region × season × fleet.}
#'     \item{\code{"est_shared_seas"}}{Shared across seasons; unique per
#'       population × region × year × fleet.}
#'     \item{\code{"est_shared_f"}}{Shared across fleets; unique per
#'       population × region × year × season.}
#'     \item{\code{"est_shared_pop_r"}}{Shared across populations and regions.}
#'     \item{\code{"est_shared_pop_y"}}{Shared across populations and years.}
#'     \item{\code{"est_shared_pop_seas"}}{Shared across populations and seasons.}
#'     \item{\code{"est_shared_pop_f"}}{Shared across populations and fleets.}
#'     \item{\code{"est_shared_r_y"}}{Shared across regions and years.}
#'     \item{\code{"est_shared_r_seas"}}{Shared across regions and seasons.}
#'     \item{\code{"est_shared_r_f"}}{Shared across regions and fleets.}
#'     \item{\code{"est_shared_y_seas"}}{Shared across years and seasons.}
#'     \item{\code{"est_shared_y_f"}}{Shared across years and fleets.}
#'     \item{\code{"est_shared_seas_f"}}{Shared across seasons and fleets.}
#'     \item{\code{"est_shared_pop_r_y"}}{Shared across populations, regions,
#'       and years.}
#'     \item{\code{"est_shared_pop_r_seas"}}{Shared across populations, regions,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_r_f"}}{Shared across populations, regions,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_y_seas"}}{Shared across populations, years,
#'       and seasons.}
#'     \item{\code{"est_shared_pop_y_f"}}{Shared across populations, years,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_seas_f"}}{Shared across populations, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_r_y_seas"}}{Shared across regions, years,
#'       and seasons.}
#'     \item{\code{"est_shared_r_y_f"}}{Shared across regions, years,
#'       and fleets.}
#'     \item{\code{"est_shared_r_seas_f"}}{Shared across regions, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_y_seas_f"}}{Shared across years, seasons,
#'       and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas"}}{Shared across populations,
#'       regions, years, and seasons.}
#'     \item{\code{"est_shared_pop_r_y_f"}}{Shared across populations, regions,
#'       years, and fleets.}
#'     \item{\code{"est_shared_pop_r_seas_f"}}{Shared across populations,
#'       regions, seasons, and fleets.}
#'     \item{\code{"est_shared_pop_y_seas_f"}}{Shared across populations, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_r_y_seas_f"}}{Shared across regions, years,
#'       seasons, and fleets.}
#'     \item{\code{"est_shared_pop_r_y_seas_f"}}{Single parameter shared across
#'       all dimensions.}
#'   }
#'
#'
#' @keywords internal
do_sigmaD_pop_mapping <- function(input_list, sigmaD_pop_spec) {

  dims <- c(pop    = input_list$data$n_pop,
            region = input_list$data$n_regions,
            year   = length(input_list$data$years),
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$ln_sigmaD_pop <- build_shared_spec_map(
    dims = dims, spec = sigmaD_pop_spec,
    dim_abbrev = c(pop = "pop", r = "region", y = "year", seas = "season", f = "fleet")
  )

  collect_message("sigmaD_pop is specified as: ", sigmaD_pop_spec)
  return(input_list)
}

#' Map discard mortality process error SD parameters
#'
#' Constructs the \code{ln_sigma_dmr} factor map used in the TMB/RTMB
#' objective function to share or fix log-scale standard deviations of
#' discard mortality process error across regions, seasons, and fleets.
#'
#' The mapping assigns integer indices to parameter groups defined by
#' \code{dmr_mean_spec}. Cells within the same group share a single
#' estimated parameter.
#'
#' @param input_list Named list containing \code{$data}, \code{$par}, and
#'   \code{$map}.
#'
#' @param dmr_mean_spec Character string specifying sharing structure.
#'   Options include \code{"est_all"}, \code{"est_shared_r"},
#'   \code{"est_shared_seas"}, \code{"est_shared_f"}, and combinations thereof,
#'   or \code{"fix"}.
#'
#' @return Updated \code{input_list} with \code{$map$ln_sigma_dmr}.
#'
#' @keywords internal
do_dmr_mean_mapping <- function(input_list, dmr_mean_spec) {

  dims <- c(region = input_list$data$n_regions,
            season = input_list$data$n_seas,
            fleet  = input_list$data$n_fish_fleets)

  input_list$map$logit_dmr_mean <- build_shared_spec_map(
    dims = dims, spec = dmr_mean_spec,
    dim_abbrev = c(r = "region", seas = "season", f = "fleet")
  )

  collect_message("dmr_mean is specified as: ", dmr_mean_spec)
  return(input_list)
}

#' Set up fishing mortality, discard mortality, and catch observation inputs
#'
#' Populates \code{input_list} with observed catch, catch usage indicators,
#' fishing mortality parameters (\code{ln_F_mean}, \code{ln_F_devs}), and
#' observation/process error structures (\code{ln_sigmaC}, \code{ln_sigmaC_pop},
#' \code{ln_sigmaF}). Also populates discard observations, discard mortality
#' rate parameters (\code{logit_dmr_mean}, \code{logit_dmr_devs}), and
#' discard observation/process error structures (\code{ln_sigmaD},
#' \code{ln_sigmaD_pop}, \code{ln_sigma_dmr}).
#' Must be called after \code{\link{Setup_Mod_Biologicals}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsCatch Observed aggregated catch array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{catch_units}. For a
#'   cell with \code{UseCatch == 0} (and no population-specific catch
#'   used), an \code{NA} entry here is treated as a genuinely missing
#'   observation; fishing is assumed to have continued and \code{Fmort}/
#'   \code{ln_F_devs} are estimated normally for that year, whereas a
#'   true recorded value (typically \code{0}) is treated as a real
#'   closure: \code{Fmort} is forced to zero and no deviation is estimated.
#'   See \code{\link{Get_Fdev_PE_loglik}}.
#' @param ObsCatch_pop Observed population-specific catch array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{catch_units}.
#' @param UseCatch Binary indicator array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]} controlling which
#'   aggregated catch observations enter the likelihood and whether
#'   \code{ln_F_devs} are estimated for each cell. \code{1} = use;
#'   \code{0} = exclude, unless \code{ObsCatch} is \code{NA} at that cell
#'   (see \code{ObsCatch} above), in which case \code{ln_F_devs} is still
#'   estimated as an ordinary active year despite not being fit against an
#'   observation.
#' @param UseCatch_pop Binary indicator array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]} controlling
#'   which population-specific catch observations enter the likelihood.
#'   \code{1} = use; \code{0} = exclude.
#' @param catch_units Character array \code{[n_fish_fleets]} specifying catch
#'   units per fleet. \code{"biom"} = biomass (default); \code{"abd"} =
#'   abundance. Converted internally to \code{0}/\code{1} integer codes.
#' @param Use_F_pen Integer flag for applying a fishing mortality penalty to
#'   penalise large deviations in \code{ln_F_devs}. \code{1} = apply
#'   (default); \code{0} = do not apply.
#' @param sigmaC_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaC} (aggregated catch observation error SD). Default
#'   \code{"fix"} holds \code{ln_sigmaC} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the convention \code{"est_shared_<dims>"} where \code{<dims>} is
#'   an underscore-separated list of dimensions to collapse: \code{"r"}
#'   (regions), \code{"y"} (years), \code{"seas"} (seasons), \code{"f"}
#'   (fleets), or any combination (e.g., \code{"est_shared_r_y"},
#'   \code{"est_shared_r_y_seas_f"}). Use \code{"est_all"} for a fully
#'   independent parameter per cell. A warning is issued if \code{"fix"} is
#'   selected without providing a starting value in \code{...}.
#' @param sigmaC_pop_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaC_pop} (population-specific catch observation error SD).
#'   Default \code{"fix"} holds \code{ln_sigmaC_pop} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaC_spec} but with an additional
#'   population dimension: e.g., \code{"est_shared_pop"} shares across
#'   populations, \code{"est_shared_pop_r"} shares across populations and
#'   regions, and \code{"est_shared_pop_r_y_seas_f"} collapses all dimensions
#'   into a single parameter. A warning is issued if \code{"fix"} is selected
#'   without providing a starting value in \code{...}.
#' @param sigmaF_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaF} (fishing mortality process error SD). Default
#'   \code{"fix"} holds \code{ln_sigmaF} at its starting value (\code{log(1)},
#'   i.e., \eqn{\sigma_F = 1}, unless overridden via \code{...}). A warning is
#'   issued if \code{"fix"} is selected without providing a starting value
#'   in \code{...}.
#' @param Fdev_pen_center Where the fishing mortality deviation penalty is
#'   centred. \code{"fixed"} (default) centres on zero, constraining both the
#'   level and the spread of the deviations. \code{"own_mean"} centres on the
#'   mean of the estimated deviations, penalizing only their spread and leaving
#'   the level free, which is what a sum of squares about the series' own mean
#'   amounts to. Under a mean-plus-deviations parameterization the level is
#'   already carried by \code{ln_F_mean}, so \code{"own_mean"} avoids
#'   penalizing it twice; note that it also leaves \code{ln_F_mean} and the
#'   deviations' level mutually unidentified unless one of them is fixed,
#'   which \code{ln_F_mean_spec = "fix"} does.
#' @param ln_F_mean_spec Character string, matched by exact name only because it
#'   sits after \code{...}. \code{"est"} (default, the previous
#'   and only behaviour) or \code{"fix"}. \code{"fix"} maps \code{ln_F_mean}
#'   off at its starting value, which defaults to \code{0} under this spec
#'   unless supplied through \code{...}, so the deviations carry all of log
#'   fishing mortality: \code{F = exp(ln_F_devs)}, where it follows a free annual log-F
#'   parameterization. It must be paired with \code{Fdev_pen_center = "own_mean"}
#'   (penalize only the spread about the deviations' own mean),
#'   \code{Fdev_model = "rw"}, or \code{Use_F_pen = 0}: an \code{"iid"} or
#'   \code{"ar1"} penalty centred on a fixed zero mean would shrink the
#'   deviations toward \code{F = 1}, so that combination is rejected at setup.
#'   \code{"est"} keeps the mean-plus-deviations form, where the \code{"iid"}
#'   penalty shrinks each year toward the estimated average F.
#' @param Fdev_model Character string specifying the process error structure
#'   for \code{ln_F_devs}. One of \code{"iid"} (default; independent
#'   deviations), \code{"rw"} (random walk; the first catch-active year per
#'   region/season/fleet is initialized with a diffuse \eqn{N(0,5)} prior),
#'   or \code{"ar1"} (first-order autoregressive; the first catch-active year
#'   is drawn from its stationary marginal distribution, and \code{Fdev_rho_spec}
#'   controls the AR1 correlation parameter). Catch-active years do not need
#'   to be contiguous for \code{"rw"} or \code{"ar1"}: the transition between
#'   two active years spanning a gap of \eqn{d} closed years is taken over
#'   the elapsed gap directly (the same marginal transition as estimating
#'   deviations for the closed years and integrating them out, without
#'   actually estimating them), see \code{\link{Get_Fdev_PE_loglik}}.
#'   A warning is issued if \code{"rw"} or \code{"ar1"} is selected but
#'   \code{Use_F_pen = 0} (the penalty is never evaluated, so the process
#'   structure has no effect), \code{sigmaF_spec = "fix"} (the process error
#'   SD is not estimated), or (for \code{"ar1"}) \code{Fdev_rho_spec =
#'   "fix"} (the correlation is not estimated), any of these may be
#'   intentional, but are common oversights when switching away from
#'   \code{"iid"}.
#' @param Fdev_rho_spec Character string specifying the sharing structure for
#'   the AR1 correlation parameter \code{Fdev_rho}, following the same
#'   convention as \code{sigmaF_spec}. Only used when \code{Fdev_model =
#'   "ar1"}; ignored (and mapped entirely to \code{NA}) otherwise.
#' @param ObsDiscard Observed aggregated discard array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{discard_units}.
#'   Default: \code{NULL} (no discard observations).
#' @param UseDiscard Binary indicator array
#'   \code{[n_regions x n_years x n_seas x n_fish_fleets]} controlling which
#'   aggregated discard observations enter the likelihood. \code{1} = use;
#'   \code{0} = exclude. Default: all zeros.
#' @param discard_units Character array \code{[n_fish_fleets]} specifying
#'   discard units per fleet. \code{"abd"} = abundance (\code{0}),
#'   \code{"biom"} = biomass (\code{1}), \code{"abd_frac"} = abundance
#'   fraction (\code{2}), \code{"biom_frac"} = biomass fraction (\code{3},
#'   default). Converted internally to integer codes.
#' @param UseDiscard_pop Binary indicator array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]} controlling
#'   which population-specific discard observations enter the likelihood.
#'   \code{1} = use; \code{0} = exclude. Default: all zeros.
#' @param ObsDiscard_pop Observed population-specific discard array
#'   \code{[n_pop x n_regions x n_years x n_seas x n_fish_fleets]}.
#'   Values should be in the units specified by \code{discard_units}.
#'   Default: \code{NULL} (no population-specific discard observations).
#' @param Use_dmr_pen Integer flag for applying a discard mortality rate
#'   penalty to penalise large deviations in \code{logit_dmr_devs}.
#'   \code{1} = apply; \code{0} = do not apply (default). Must be \code{1}
#'   when \code{dmr_dev_spec = "est_all"} and \code{0} when
#'   \code{dmr_dev_spec = "fix"}.
#' @param sigmaD_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaD} (aggregated discard observation error SD). Default
#'   \code{"fix"} holds \code{ln_sigmaD} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaC_spec}. A warning is issued
#'   if \code{"fix"} is selected without providing a starting value in
#'   \code{...}.
#' @param sigmaD_pop_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaD_pop} (population-specific discard observation error SD).
#'   Default \code{"fix"} holds \code{ln_sigmaD_pop} at its starting value
#'   (\code{log(0.01)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaC_pop_spec}. A warning is
#'   issued if \code{"fix"} is selected without providing a starting value
#'   in \code{...}.
#' @param sigma_dmr_spec Character string specifying the sharing structure for
#'   \code{ln_sigma_dmr} (discard mortality rate process error SD). Default
#'   \code{"fix"} holds \code{ln_sigma_dmr} at its starting value
#'   (\code{log(1)} unless overridden via \code{...}). Sharing options
#'   follow the same convention as \code{sigmaF_spec}. A warning is issued
#'   if \code{"fix"} is selected without providing a starting value in
#'   \code{...}.
#' @param dmr_mean_spec Character string specifying the sharing/estimation
#'   structure for \code{logit_dmr_mean} (logit-scale mean discard mortality
#'   rate). Default \code{"fix"} holds at its starting value (\code{0},
#'   i.e., DMR = 0.5 on the natural scale, unless overridden via \code{...}).
#'   See \code{\link{do_dmr_mean_mapping}} for sharing options.
#' @param dmr_dev_spec Character string specifying the sharing/estimation
#'   structure for \code{logit_dmr_devs} (logit-scale annual discard mortality
#'   rate deviations). Default \code{"fix"} holds deviations at zero
#'   (unless overridden via \code{...}). Use \code{"est_all"} to estimate a
#'   deviation in every fished cell; requires \code{Use_dmr_pen = 1}. See
#'   \code{\link{do_dmr_dev_mapping}} for sharing options.
#' @param ... Optional starting value overrides for catch and discard related parameters.
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated. Key additions:
#'   \describe{
#'     \item{\code{$data}}{
#'       \code{ObsCatch}, \code{ObsCatch_pop}, \code{UseCatch},
#'       \code{UseCatch_pop}, \code{Use_F_pen}, \code{catch_units},
#'       \code{Fdev_model},
#'       \code{ObsDiscard}, \code{ObsDiscard_pop}, \code{UseDiscard},
#'       \code{UseDiscard_pop}, \code{Use_dmr_pen}, \code{discard_units}.}
#'     \item{\code{$par}}{
#'       \code{ln_sigmaC}, \code{ln_sigmaC_pop}, \code{ln_sigmaF},
#'       \code{Fdev_rho},
#'       \code{ln_F_mean}, \code{ln_F_devs},
#'       \code{ln_sigmaD}, \code{ln_sigmaD_pop}, \code{ln_sigma_dmr},
#'       \code{logit_dmr_mean}, \code{logit_dmr_devs}.}
#'     \item{\code{$map}}{
#'       \code{ln_sigmaC}, \code{ln_sigmaC_pop}, \code{ln_sigmaF},
#'       \code{Fdev_rho},
#'       \code{ln_F_mean}, \code{ln_F_devs},
#'       \code{ln_sigmaD}, \code{ln_sigmaD_pop}, \code{ln_sigma_dmr},
#'       \code{logit_dmr_mean}, \code{logit_dmr_devs}.}
#'   }
#'
#' @export Setup_Mod_Catch_and_F
#' @family Model Setup
Setup_Mod_Catch_and_F <- function(input_list,

                                  # Retained Catch Stuff
                                  ObsCatch,
                                  UseCatch,
                                  catch_units = array("biom", dim = c(input_list$data$n_fish_fleets)),
                                  UseCatch_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                                  length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets )),
                                  ObsCatch_pop = NULL,
                                  Use_F_pen = 1,
                                  sigmaC_spec = "fix",
                                  sigmaC_pop_spec = 'fix',
                                  sigmaF_spec = "fix",
                                  Fdev_model = "iid",
                                  Fdev_pen_center = "fixed",
                                  Fdev_rho_spec = "fix",

                                  # Discarded Catch Stuff
                                  ObsDiscard = NULL,
                                  UseDiscard = array(0, dim = c(input_list$data$n_regions,
                                                                length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets )),
                                  discard_units = array("biom_frac", dim = c(input_list$data$n_fish_fleets)),
                                  UseDiscard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                                    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets )),
                                  ObsDiscard_pop = NULL,
                                  Use_dmr_pen = 0,
                                  sigmaD_spec = "fix",
                                  sigmaD_pop_spec = 'fix',
                                  sigma_dmr_spec = "fix",
                                  dmr_mean_spec = 'fix',
                                  dmr_dev_spec = 'fix',
                                  ...,
                                  ln_F_mean_spec = "est") {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Catch_and_F <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Catch objects
  check_data_dimensions(ObsCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsCatch')
  check_data_dimensions(UseCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseCatch')

  if(any(UseCatch_pop == 1)) {
    check_data_dimensions(ObsCatch_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsCatch_pop')
    check_data_dimensions(UseCatch_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseCatch_pop')
  }

  if(any(UseDiscard == 1)) {
    check_data_dimensions(ObsDiscard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsDiscard')
    check_data_dimensions(UseDiscard, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseDiscard')
  }

  if(any(UseDiscard_pop == 1)) {
    check_data_dimensions(ObsDiscard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsDiscard_pop')
    check_data_dimensions(UseDiscard_pop, n_pop = input_list$data$n_pop, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseDiscard_pop')
  }

  # Fishing Mortality checking
  if(!Use_F_pen %in% c(0,1)) stop("Use_F_pen incorrectly specified. Either set at 0 (don't use F penalty) or 1 (use F penalty)")
  else collect_message("Fishing mortality penalty is: ", ifelse(Use_F_pen == 0, 'Not Used', "Used"))
  if(sigmaC_spec == "fix" && !("ln_sigmaC" %in% names(starting_values))) warning("sigmaC is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaC_pop_spec == "fix" && !("ln_sigmaC_pop" %in% names(starting_values))) warning("sigmaC_pop is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaF_spec == "fix" && !("ln_sigmaF" %in% names(starting_values))) warning("sigmaF_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")

  # Fdev_model checking
  if(!Fdev_model %in% c("iid", "rw", "ar1")) stop("Fdev_model incorrectly specified. Must be one of 'iid', 'rw', or 'ar1'")
  else collect_message("Fdev_model is specified as: ", Fdev_model)

  if(Fdev_model %in% c("rw", "ar1") && Use_F_pen == 0)
    warning("Fdev_model = '", Fdev_model, "' but Use_F_pen = 0; the fishing mortality deviation penalty is never evaluated, so the ", Fdev_model, " process error structure has no effect on the model. Set Use_F_pen = 1 to actually apply it.")

  if(Fdev_model %in% c("rw", "ar1") && sigmaF_spec == "fix")
    warning("Fdev_model = '", Fdev_model, "' but sigmaF_spec = 'fix'; the process error standard deviation (ln_sigmaF) driving the ", Fdev_model, " process is not being estimated. This may be intentional (e.g. fixing sigma at a known value), but if not, consider estimating ln_sigmaF via sigmaF_spec.")

  if(Fdev_model == "ar1" && Fdev_rho_spec == "fix")
    warning("Fdev_model = 'ar1' but Fdev_rho_spec = 'fix'; the AR1 correlation parameter (Fdev_rho) is not being estimated. This may be intentional (e.g. fixing rho at a known value), but if not, consider estimating Fdev_rho via Fdev_rho_spec.")

  # Discard Mortality checking
  if(!Use_dmr_pen %in% c(0,1)) stop("Use_dmr_pen incorrectly specified. Either set at 0 (don't use D penalty) or 1 (use D penalty)")
  else collect_message("Discard mortality penalty is: ", ifelse(Use_dmr_pen == 0, 'Not Used', "Used"))
  if(sigmaD_spec == "fix" && !("ln_sigmaD" %in% names(starting_values))) warning("sigmaD is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaD_pop_spec == "fix" && !("ln_sigmaD_pop" %in% names(starting_values))) warning("sigmaD_pop is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigma_dmr_spec == "fix" && !("ln_sigma_dmr" %in% names(starting_values))) warning("sigma_dmr_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")

  # Validation checks for dmr_dev_spec and Use_dmr_pen consistency
  if(dmr_dev_spec == "est_all" && Use_dmr_pen == 0)
    warning("dmr_dev_spec is 'est_all' but Use_dmr_pen is 0. Estimating dmr deviations without a penalty will likely cause convergence issues.")

  if(dmr_dev_spec == "fix" && Use_dmr_pen == 1)
    stop("dmr_dev_spec is 'fix' but Use_dmr_pen is 1. Cannot apply a penalty on deviations that are not estimated.")

  # Catch units
  catch_units[catch_units == 'abd'] <- 0
  catch_units[catch_units == 'biom'] <- 1
  catch_units <- array(as.numeric(catch_units), dim = c(input_list$data$n_fish_fleets)) # convert to numeric array

  # Discard units
  discard_units[discard_units == 'abd'] <- 0
  discard_units[discard_units == 'biom'] <- 1
  discard_units[discard_units == 'abd_frac'] <- 2
  discard_units[discard_units == 'biom_frac'] <- 3
  discard_units <- array(as.numeric(discard_units), dim = c(input_list$data$n_fish_fleets)) # convert to numeric array

  # Populate Data List ------------------------------------------------------

  # Retained Catch Stuff
  input_list$data$ObsCatch <- ObsCatch
  input_list$data$UseCatch <- UseCatch
  input_list$data$ObsCatch_pop <- ObsCatch_pop
  input_list$data$UseCatch_pop <- UseCatch_pop
  input_list$data$Use_F_pen <- Use_F_pen
  input_list$data$catch_units <- catch_units
  input_list$data$Fdev_model <- match(Fdev_model, c("iid", "rw", "ar1")) # 1 = iid, 2 = rw, 3 = ar1
  if(!Fdev_pen_center %in% c("fixed", "own_mean")) stop("Fdev_pen_center must be fixed or own_mean")
  input_list$data$Fdev_pen_center <- convert_to_numeric(Fdev_pen_center, list(fixed = 0, own_mean = 1))
  if(!ln_F_mean_spec %in% c("est", "fix")) stop("ln_F_mean_spec must be est or fix")
  collect_message("ln_F_mean is specified as: ", ln_F_mean_spec)

  # A fixed zero mean with a penalty centred on that mean shrinks the deviations
  # toward F = 1, which is a prior nobody intends, so the combination is rejected
  # rather than warned about.
  if(ln_F_mean_spec == "fix" && Use_F_pen == 1 && Fdev_pen_center == "fixed" && Fdev_model %in% c("iid", "ar1"))
    stop("ln_F_mean_spec = 'fix' with a zero-centered '", Fdev_model, "' penalty shrinks the deviations toward F = 1. Pair it with Fdev_pen_center = 'own_mean', Fdev_model = 'rw', or Use_F_pen = 0.")

  # Under own-mean centering the penalty no longer pins the level of log F, so
  # an estimated ln_F_mean and the deviations' mean trade off along an exactly
  # flat ridge unless something else still reads the mean: "prop" initialization
  # F does, but absolute-rate initialization ("abs") and a free initial age
  # structure (init_age_strc = 4) do not.
  if(Fdev_pen_center == "own_mean" && ln_F_mean_spec == "est" &&
     (isTRUE(input_list$data$init_F_form == 1) || isTRUE(input_list$data$init_age_strc == 4)))
    warning("Fdev_pen_center = 'own_mean' with an estimated ln_F_mean leaves the mean and the deviations' level mutually unidentified in this configuration, since nothing else reads ln_F_mean. Consider ln_F_mean_spec = 'fix'.")

  # Discarded Catch Stuff
  input_list$data$ObsDiscard <- ObsDiscard
  input_list$data$UseDiscard <- UseDiscard
  input_list$data$ObsDiscard_pop <- ObsDiscard_pop
  input_list$data$UseDiscard_pop <- UseDiscard_pop
  input_list$data$Use_dmr_pen <- Use_dmr_pen
  input_list$data$discard_units <- discard_units

  # Populate Parameter List -------------------------------------------------

  # Catch observation error
  if("ln_sigmaC" %in% names(starting_values)) input_list$par$ln_sigmaC <- starting_values$ln_sigmaC
  else input_list$par$ln_sigmaC <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  if("ln_sigmaC_pop" %in% names(starting_values)) input_list$par$ln_sigmaC_pop <- starting_values$ln_sigmaC_pop
  else input_list$par$ln_sigmaC_pop <- array(log(0.01), dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Process error fishing deviations
  if("ln_sigmaF" %in% names(starting_values)) input_list$par$ln_sigmaF <- starting_values$ln_sigmaF
  else input_list$par$ln_sigmaF <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # AR1 correlation for fishing mortality deviations (only used when Fdev_model = "ar1")
  if("Fdev_rho" %in% names(starting_values)) input_list$par$Fdev_rho <- starting_values$Fdev_rho
  else input_list$par$Fdev_rho <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Log mean fishing mortality. A fixed mean defaults to zero so the deviations
  # are log F outright.
  if("ln_F_mean" %in% names(starting_values)) input_list$par$ln_F_mean <- starting_values$ln_F_mean
  else input_list$par$ln_F_mean <- array(if(ln_F_mean_spec == "fix") 0 else log(0.1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Log fishing deviations
  if("ln_F_devs" %in% names(starting_values)) input_list$par$ln_F_devs <- starting_values$ln_F_devs
  else input_list$par$ln_F_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Discard observation error
  if("ln_sigmaD" %in% names(starting_values)) input_list$par$ln_sigmaD <- starting_values$ln_sigmaD
  else input_list$par$ln_sigmaD <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  if("ln_sigmaD_pop" %in% names(starting_values)) input_list$par$ln_sigmaD_pop <- starting_values$ln_sigmaD_pop
  else input_list$par$ln_sigmaD_pop <- array(log(0.01), dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Process error discard deviations
  if("ln_sigma_dmr" %in% names(starting_values)) input_list$par$ln_sigma_dmr <- starting_values$ln_sigma_dmr
  else input_list$par$ln_sigma_dmr <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Logit mean discard mortality
  if("logit_dmr_mean" %in% names(starting_values)) input_list$par$logit_dmr_mean <- starting_values$logit_dmr_mean
  else input_list$par$logit_dmr_mean <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Logit discard mortality deviations
  if("logit_dmr_devs" %in% names(starting_values)) input_list$par$logit_dmr_devs <- starting_values$logit_dmr_devs
  else input_list$par$logit_dmr_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------

  # Retained Catch Stuff
  input_list <- do_sigmaC_mapping(input_list, sigmaC_spec)
  input_list <- do_sigmaC_pop_mapping(input_list, sigmaC_pop_spec)
  input_list <- do_sigmaF_mapping(input_list, sigmaF_spec)
  input_list <- do_Fdev_rho_mapping(input_list, Fdev_rho_spec)
  input_list <- do_Fmort_mapping(input_list)
  # Free log-F parameterization: the mean is fixed at its starting value and the
  # deviations carry all of log F.
  if(ln_F_mean_spec == "fix") input_list$map$ln_F_mean <- factor(rep(NA, length(input_list$par$ln_F_mean)))

  # Discard Catch Stuff
  input_list <- do_sigmaD_mapping(input_list, sigmaD_spec)
  input_list <- do_sigmaD_pop_mapping(input_list, sigmaD_pop_spec)
  input_list <- do_sigma_dmr_mapping(input_list, sigma_dmr_spec)
  input_list <- do_dmr_mean_mapping(input_list, dmr_mean_spec)
  input_list <- do_dmr_dev_mapping(input_list, dmr_dev_spec)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

