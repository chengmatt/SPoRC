# Stage 1 of 3: model setup
#
# Catch and fishing mortality inputs, plus discard mortality rate. Setup_Mod_Catch_and_F builds the
# maps for F, its deviations, the catch and discard observation error, and the discard mortality rate.

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
    dims = dims,
    spec = sigmaF_spec,
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
      dims = dims,
      spec = Fdev_rho_spec,
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
    dims = dims,
    spec = sigmaC_spec,
    dim_abbrev = c(r = "region", y = "year", seas = "season", f = "fleet"),
    use = input_list$data$UseCatch,
    what = "sigmaC"
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

  # estimate F devs if any catch is used, or if the aggregate catch is missing rather than a
  # recorded zero: fishing is assumed to continue through a missing observation
  has_catch <- input_list$data$UseCatch == 1 |
    apply(input_list$data$UseCatch_pop == 1, c(2,3,4,5), any) |
    is.na(input_list$data$ObsCatch)

  # a fleet fitting catch at age keeps its observations in a different array and is fished wherever
  # any age or sex is fit. without this the fishing mortality level is mapped off silently
  if(!is.null(input_list$data$UseCatchAA) && any(input_list$data$UseCatchAA == 1)) {
    has_catch <- has_catch | apply(input_list$data$UseCatchAA == 1, c(1,2,3,6), any)
  }
  if(!is.null(input_list$data$UseCatchAA_pop) && any(input_list$data$UseCatchAA_pop == 1)) {
    has_catch <- has_catch | apply(input_list$data$UseCatchAA_pop == 1, c(2,3,4,7), any)
  }

  # a fleet reporting one annual total still fishes every season of a year it reports in, so
  # spread that year's observation across the seasons rather than leaving the rest unfished
  # an input list built by hand, or by an older version, has no seasonal reporting fields
  or_seasonal <- function(x) if(is.null(x)) rep(0, dims[["fleet"]]) else x
  seas_agg <- rbind(or_seasonal(input_list$data$Catch_seas_Type),
                    or_seasonal(input_list$data$Catch_pop_seas_Type),
                    or_seasonal(input_list$data$CatchAA_seas_Type),
                    or_seasonal(input_list$data$CatchAA_pop_seas_Type))

  for(f in which(apply(seas_agg == 1, 2, any))) {
    fished_year <- apply(has_catch[,,,f, drop = FALSE], c(1,2), any) # region by year
    for(seas in seq_len(dims[["season"]])) has_catch[,,seas,f] <- fished_year
  } # end f loop

  F_dev_map <- build_pe_map(dims, share_over = character(0))
  F_dev_map[!has_catch] <- NA

  input_list$map$ln_F_devs <- factor(as.vector(F_dev_map))

  # ln_F_mean only enters through cells fished in at least one year. a region-season-fleet cell
  # closed in every year contributes nothing, so estimating it leaves a flat gradient
  F_mean_dims <- dims[c("region", "season", "fleet")]
  F_mean_active <- apply(has_catch, c(1,3,4), any)

  F_mean_map <- build_pe_map(F_mean_dims, share_over = character(0))
  F_mean_map[!F_mean_active] <- NA

  input_list$map$ln_F_mean <- factor(as.vector(F_mean_map))

  # mirror the deviation map into the data list so the process error penalty keys on the cells that
  # are estimated. a deviation mapped off by hand after setup is then neither estimated nor penalized
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
    dims = dims,
    spec = sigmaC_pop_spec,
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
    dims = dims,
    spec = sigma_dmr_spec,
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
    dims = dims,
    spec = sigmaD_spec,
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

  # estimate a deviation wherever dmr enters the dynamics: any catch used, or the aggregate catch
  # missing rather than a recorded zero. get_dmr_penalty() keys on the resulting map
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
    dims = dims,
    spec = sigmaD_pop_spec,
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
    dims = dims,
    spec = dmr_mean_spec,
    dim_abbrev = c(r = "region", seas = "season", f = "fleet")
  )

  collect_message("dmr_mean is specified as: ", dmr_mean_spec)
  return(input_list)
}

#' Set up fishing mortality, discard mortality, and catch observation inputs
#'
#' Sets the observed catch and discards with their use flags, the fishing
#' mortality parameters (\code{ln_F_mean}, \code{ln_F_devs}) and their observation
#' and process error, the catch and discard at age data sources, and the discard
#' mortality rate parameters. Call after \code{\link{Setup_Mod_Biologicals}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map} and
#'   \code{$verbose}.
#' @param ObsCatch Observed aggregated catch array \code{[n_regions x n_years x
#'   n_seas x n_fish_fleets]} in the units \code{catch_units} names. Where
#'   \code{UseCatch == 0} and no population-specific catch is used, an \code{NA}
#'   here is a missing observation: fishing is assumed to have continued and
#'   \code{Fmort} and \code{ln_F_devs} are estimated as usual. A recorded value,
#'   typically \code{0}, is a real closure: \code{Fmort} is forced to zero and no
#'   deviation is estimated. See \code{\link{Get_Fdev_PE_loglik}}.
#' @param ObsCatch_pop Observed population-specific catch array \code{[n_pop x
#'   n_regions x n_years x n_seas x n_fish_fleets]}, in \code{catch_units}.
#' @param UseCatch Binary array dimensioned like \code{ObsCatch} controlling which
#'   aggregated catch observations enter the likelihood and whether
#'   \code{ln_F_devs} is estimated in each cell. \code{0} excludes the
#'   observation, unless \code{ObsCatch} is \code{NA} there, in which case the
#'   deviation is still estimated.
#' @param UseCatch_pop Binary array dimensioned like \code{ObsCatch_pop}.
#' @param catch_units Character array \code{[n_fish_fleets]}: \code{"biom"}
#'   (default) or \code{"abd"}, stored as \code{0}/\code{1}.
#' @param Use_F_pen Integer flag for the fishing mortality penalty on
#'   \code{ln_F_devs}. \code{1} (default) applies it.
#' @param sigmaC_spec Sharing structure for \code{ln_sigmaC}, the aggregated catch
#'   observation error sd. \code{"fix"} (default) holds it at its starting value,
#'   \code{log(0.01)} unless supplied through \code{...}, and warns when no
#'   starting value was given. Estimated options are
#'   \code{"est_shared_<dims>"} over any of \code{"r"} (regions), \code{"y"}
#'   (years), \code{"seas"} and \code{"f"} (fleets), e.g.
#'   \code{"est_shared_r_y_seas_f"}, or \code{"est_all"} for one parameter per
#'   cell.
#' @param sigmaC_pop_spec Sharing structure for \code{ln_sigmaC_pop}, as
#'   \code{sigmaC_spec} with an added population dim, e.g.
#'   \code{"est_shared_pop_r"} or \code{"est_shared_pop_r_y_seas_f"}.
#' @param sigmaF_spec Sharing structure for \code{ln_sigmaF}, the fishing mortality
#'   process error sd, following \code{sigmaC_spec}. \code{"fix"} (default) holds
#'   it at \code{log(1)} unless supplied through \code{...}, and warns.
#' @param Fdev_pen_center Where the fishing mortality deviation penalty is
#'   centered. \code{"fixed"} (default) centers on zero, constraining the level
#'   and the spread. \code{"own_mean"} centers on the deviations' own mean,
#'   penalizing only their spread; the level is then already set by
#'   \code{ln_F_mean}, so it is not penalized twice, but the two are mutually
#'   unidentified unless one is fixed, which \code{ln_F_mean_spec = "fix"} does.
#' @param ln_F_mean_spec \code{"est"} (default) or \code{"fix"}, matched by exact
#'   name only because it sits after \code{...}. \code{"fix"} maps
#'   \code{ln_F_mean} off at its starting value, \code{0} unless supplied through
#'   \code{...}, so the deviations hold all of log fishing mortality,
#'   \code{F = exp(ln_F_devs)}. It must be paired with
#'   \code{Fdev_pen_center = "own_mean"}, \code{Fdev_model = "rw"} or
#'   \code{Use_F_pen = 0}: an \code{"iid"} or \code{"ar1"} penalty centered on a
#'   fixed zero would shrink the deviations toward \code{F = 1}, so that
#'   combination is rejected at setup.
#' @param Fdev_model Process error on \code{ln_F_devs}: \code{"iid"} (default),
#'   \code{"rw"} (the first catch-active year per region, season and fleet takes a
#'   diffuse \eqn{N(0,5)}), or \code{"ar1"} (that year is drawn from the
#'   stationary marginal, with \code{Fdev_rho_spec} setting the correlation).
#'   Catch-active years need not be contiguous under \code{"rw"} or \code{"ar1"}:
#'   the transition across a gap of \eqn{d} closed years is taken over the elapsed
#'   gap, the same marginal as estimating the closed years and integrating them
#'   out. See \code{\link{Get_Fdev_PE_loglik}}. Warns under \code{Use_F_pen = 0}
#'   (the penalty is never evaluated), \code{sigmaF_spec = "fix"}, or, for
#'   \code{"ar1"}, \code{Fdev_rho_spec = "fix"}.
#' @param Fdev_rho_spec Sharing structure for \code{Fdev_rho}, following
#'   \code{sigmaF_spec}. Only read under \code{Fdev_model = "ar1"} and mapped
#'   entirely to \code{NA} otherwise.
#' @param ObsDiscard Observed aggregated discard array \code{[n_regions x n_years x
#'   n_seas x n_fish_fleets]} in \code{discard_units}. Default \code{NULL}.
#' @param UseDiscard Binary array dimensioned like \code{ObsDiscard}. Default all
#'   zeros.
#' @param discard_units Character array \code{[n_fish_fleets]}: \code{"abd"}
#'   (\code{0}), \code{"biom"} (\code{1}), \code{"abd_frac"} (\code{2}) or
#'   \code{"biom_frac"} (\code{3}, default).
#' @param UseDiscard_pop Binary array \code{[n_pop x n_regions x n_years x n_seas x
#'   n_fish_fleets]}. Default all zeros.
#' @param ObsDiscard_pop Observed population-specific discard array, same dims, in
#'   \code{discard_units}. Default \code{NULL}.
#' @param Use_dmr_pen Integer flag for the penalty on \code{logit_dmr_devs}.
#'   Default \code{0}. Must be \code{1} under \code{dmr_dev_spec = "est_all"} and
#'   \code{0} under \code{"fix"}.
#' @param sigmaD_spec,sigmaD_pop_spec Sharing structures for \code{ln_sigmaD} and
#'   \code{ln_sigmaD_pop}, the discard observation error sds, following
#'   \code{sigmaC_spec} and \code{sigmaC_pop_spec}. \code{"fix"} (default) holds
#'   them at \code{log(0.01)} and warns when no starting value was given.
#' @param sigma_dmr_spec Sharing structure for \code{ln_sigma_dmr}, the discard
#'   mortality rate process error sd, following \code{sigmaF_spec}. \code{"fix"}
#'   (default) holds it at \code{log(1)} and warns.
#' @param dmr_mean_spec Sharing structure for \code{logit_dmr_mean}. \code{"fix"}
#'   (default) holds it at \code{0}, a rate of 0.5 on the natural scale. See
#'   \code{\link{do_dmr_mean_mapping}}.
#' @param dmr_dev_spec Sharing structure for \code{logit_dmr_devs}. \code{"fix"}
#'   (default) holds the deviations at zero; \code{"est_all"} estimates one in
#'   every fished cell and requires \code{Use_dmr_pen = 1}. See
#'   \code{\link{do_dmr_dev_mapping}}.
#' @param ObsCatchAA Observed catch at age \code{[n_regions, n_years, n_seas,
#'   n_obs_ages, n_sexes, n_fish_fleets]}, the ages being the columns of the
#'   fleet's ageing error matrix, through which the predicted catch at each model
#'   age is read before it is compared. The sex dim is required whatever the fleet
#'   reports: a data source summed over sexes has its observation in sex slot one.
#'   Supplying this fits the catch at age directly, every age its own lognormal
#'   observation, in place of an aggregated catch with compositions, which is the
#'   native form for ICES age-structured assessments. The exact factorization of an
#'   at-age observation into a total and a composition holds for Poisson and
#'   multinomial but not lognormal, so a fleet must use one or the other and
#'   supplying both is an error. \code{NULL} (default) keeps the fleet on
#'   aggregated catch.
#' @param UseCatchAA Integer array shaped like \code{ObsCatchAA}, \code{1} where an
#'   observation is fit. A cell that is not fit is also not fished, so this governs
#'   closures the way \code{UseCatch} does.
#' @param sigmaCAA_key Integer array \code{[n_obs_ages, n_sexes, n_fish_fleets]}
#'   coupling the catch at age observation error, the key matrix ICES assessments
#'   use. Equal entries share a parameter and \code{NA} excludes one. The sex dim
#'   is required; a key coupling the sexes repeats its entries across them. Along
#'   ages, \code{1 2 3 4 5} gives one sd per age, \code{1 1 2 2 2} gives sds by age
#'   group, and \code{1 1 1 1 1} gives one for the fleet. Defaults to one parameter
#'   per fleet. A parameter informed by fewer than two observations is refused,
#'   since an sd with a single observation drives the likelihood to negative
#'   infinity rather than failing outright.
#' @param sigmaCAA_spec \code{"est"} (default) or \code{"fix"}. Starting values go
#'   through \code{...} as \code{ln_sigmaCAA}.
#' @param ObsDiscardAA,UseDiscardAA Observed discard at age and its use flags,
#'   shaped like \code{ObsCatchAA} and read through the same fishery ageing error.
#' @param ObsDiscardAA_pop,UseDiscardAA_pop,ObsCatchAA_pop,UseCatchAA_pop
#'   Population-specific counterparts, with a leading population dim.
#' @param sigmaCAA_pop_key,sigmaDAA_key,sigmaDAA_pop_key Integer arrays coupling
#'   the observation error for the population-specific catch, the discards and the
#'   population-specific discards, following \code{sigmaCAA_key}.
#'   \code{sigmaDAA_key} is \code{[n_obs_ages, n_sexes, n_fish_fleets]}; the two
#'   population-specific keys take a leading population dim.
#' @param sigmaCAA_pop_spec,sigmaDAA_spec,sigmaDAA_pop_spec \code{"est"} or
#'   \code{"fix"}.
#' @param AgeObsCorr_catch,AgeObsCorr_discard,AgeObsCorr_catch_pop,AgeObsCorr_discard_pop
#'   Correlation across ages within a cell, one setting for every fleet or one per
#'   fleet. \code{"iid"} (default) treats ages as independent, \code{"1dar1"}
#'   correlates them as an AR(1) in age distance, \code{"us"} estimates an
#'   unstructured correlation, and \code{"2dar1"} correlates over ages and years
#'   jointly through a separable AR(1), which needs the fleet's observed ages and
#'   years to form a complete grid. A cell with one observed age falls back to
#'   independent. The population-specific data sources have their own settings. The
#'   index data sources are set in \code{\link{Setup_Mod_FishIdx_and_Comps}} and
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}.
#' @param rho_catch_spec,rho_discard_spec,rho_catch_pop_spec,rho_discard_pop_spec
#'   How each data source's correlation parameters are shared, in the spec strings
#'   \code{sigmaF_spec} uses. The correlations sit over region, sex and fleet, with
#'   a leading population dim for the population-specific data sources, so
#'   \code{"est_shared_r_s"} gives one per fleet, \code{"est_shared_s"} one per
#'   region and fleet, \code{"est_shared_r_s_f"} a single value, \code{"est_all"}
#'   one per cell, and \code{"fix"} holds them. \code{NULL} (default) takes
#'   \code{"est_shared_r_s"}, or \code{"est_shared_p_r_s"} for the population data
#'   sources. The spec governs the across-age correlation, the across-year
#'   correlation and the unstructured matrix together, so fleets sharing under
#'   \code{"us"} share a whole matrix. A region, sex or population a fleet never
#'   observes has no parameter.
#' @param ObsCatchAA_SE,ObsDiscardAA_SE,ObsCatchAA_pop_SE,ObsDiscardAA_pop_SE
#'   Reported standard errors shaped like their observation array, read only when
#'   that data source's \code{sigma_form} asks for them.
#' @param Catch_seas_Type,Catch_pop_seas_Type,Discard_seas_Type,Discard_pop_seas_Type,CatchAA_seas_Type,CatchAA_pop_seas_Type,DiscardAA_seas_Type,DiscardAA_pop_seas_Type
#'   Whether a seasonal model reports this data source once a season or once a
#'   year, one value for every fleet or one per fleet. \code{"spltSeas"} (default)
#'   fits the observation against the prediction for the season it sits in;
#'   \code{"aggSeas"} sums the prediction over the year's seasons and fits one
#'   observation, which is how a fleet that lands catch all year but reports one
#'   annual total is usually recorded. Under \code{"aggSeas"} the observation stays
#'   in the season it was placed in and exactly one season per region and year may
#'   be on in the matching \code{Use} array, since more than one would be fit
#'   against the same year total. Fishing mortality is still estimated season by
#'   season, so a fleet with one annual observation and free seasonal deviations
#'   leaves the split between seasons unidentified: share the deviations or fix the
#'   seasonal pattern.
#' @param CatchAA_Type,DiscardAA_Type,CatchAA_pop_Type,DiscardAA_pop_Type Which
#'   dims the fleet reports separately, in the composition vocabulary, as one
#'   setting for every fleet, one per fleet, or year and fleet specifications such
#'   as \code{"spltRaggS_Year_1-20_Fleet_1"}. \code{"agg"} sums over regions and
#'   sexes, \code{"spltRaggS"} (default) splits regions and sums over sexes,
#'   \code{"aggRspltS"} does the reverse, and \code{"spltRspltS"} splits both. An
#'   observation summed over a dim belongs in slot one of it.
#' @param CatchAA_LikeType,DiscardAA_LikeType,CatchAA_pop_LikeType,DiscardAA_pop_LikeType
#'   \code{"lognormal"} (default) or \code{"normal"}, one setting for every fleet
#'   or one per fleet.
#' @param CatchAA_sigma_form,DiscardAA_sigma_form,CatchAA_pop_sigma_form,DiscardAA_pop_sigma_form
#'   Where the observation error comes from. \code{"none"} (default) uses the
#'   estimated parameter alone, \code{"data"} the reported standard errors alone,
#'   and \code{"est_additive"} or \code{"est_quadrature"} both. Naming
#'   \code{"data"} holds the parameter fixed, since nothing reads it.
#' @param ... Optional starting values for the catch and discard parameters.
#'
#' @return \code{input_list} with \code{$data}, \code{$par} and \code{$map}
#'   updated. \code{$data} gains \code{ObsCatch}, \code{ObsCatch_pop},
#'   \code{UseCatch}, \code{UseCatch_pop}, \code{Use_F_pen}, \code{catch_units},
#'   \code{Fdev_model}, \code{ObsDiscard}, \code{ObsDiscard_pop},
#'   \code{UseDiscard}, \code{UseDiscard_pop}, \code{Use_dmr_pen} and
#'   \code{discard_units}. \code{$par} and \code{$map} both gain
#'   \code{ln_sigmaC}, \code{ln_sigmaC_pop}, \code{ln_sigmaF}, \code{Fdev_rho},
#'   \code{ln_F_mean}, \code{ln_F_devs}, \code{ln_sigmaD}, \code{ln_sigmaD_pop},
#'   \code{ln_sigma_dmr}, \code{logit_dmr_mean} and \code{logit_dmr_devs}.
#'
#' @export Setup_Mod_Catch_and_F
#' @family Model Setup
Setup_Mod_Catch_and_F <- function(input_list,

                                  # Retained Catch Stuff
                                  ObsCatch,
                                  ObsCatchAA = NULL,
                                  UseCatchAA = NULL,
                                  ObsCatchAA_SE = NULL,
                                  sigmaCAA_key = NULL,
                                  sigmaCAA_spec = "est",
                                  ObsDiscardAA = NULL,
                                  UseDiscardAA = NULL,
                                  ObsDiscardAA_SE = NULL,
                                  ObsDiscardAA_pop = NULL,
                                  UseDiscardAA_pop = NULL,
                                  ObsDiscardAA_pop_SE = NULL,
                                  ObsCatchAA_pop = NULL,
                                  UseCatchAA_pop = NULL,
                                  ObsCatchAA_pop_SE = NULL,
                                  sigmaCAA_pop_key = NULL,
                                  sigmaCAA_pop_spec = "est",
                                  sigmaDAA_key = NULL,
                                  sigmaDAA_spec = "est",
                                  sigmaDAA_pop_key = NULL,
                                  sigmaDAA_pop_spec = "est",
                                  CatchAA_Type = "spltRaggS",
                                  CatchAA_pop_Type = "spltRaggS",
                                  DiscardAA_Type = "spltRaggS",
                                  DiscardAA_pop_Type = "spltRaggS",
                                  Catch_seas_Type = NULL,
                                  Catch_pop_seas_Type = NULL,
                                  Discard_seas_Type = NULL,
                                  Discard_pop_seas_Type = NULL,
                                  CatchAA_seas_Type = NULL,
                                  CatchAA_pop_seas_Type = NULL,
                                  DiscardAA_seas_Type = NULL,
                                  DiscardAA_pop_seas_Type = NULL,
                                  CatchAA_LikeType = "lognormal",
                                  CatchAA_pop_LikeType = "lognormal",
                                  DiscardAA_LikeType = "lognormal",
                                  DiscardAA_pop_LikeType = "lognormal",
                                  CatchAA_sigma_form = "none",
                                  CatchAA_pop_sigma_form = "none",
                                  DiscardAA_sigma_form = "none",
                                  DiscardAA_pop_sigma_form = "none",
                                  AgeObsCorr_catch = "iid",
                                  AgeObsCorr_catch_pop = "iid",
                                  AgeObsCorr_discard = "iid",
                                  AgeObsCorr_discard_pop = "iid",
                                  rho_catch_spec = NULL,
                                  rho_catch_pop_spec = NULL,
                                  rho_discard_spec = NULL,
                                  rho_discard_pop_spec = NULL,
                                  UseCatch,
                                  catch_units = array("biom", dim = c(input_list$data$n_fish_fleets)),
                                  UseCatch_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                                  length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
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
                                                                length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                  discard_units = array("biom_frac", dim = c(input_list$data$n_fish_fleets)),
                                  UseDiscard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
                                                                    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets)),
                                  ObsDiscard_pop = NULL,
                                  Use_dmr_pen = 0,
                                  sigmaD_spec = "fix",
                                  sigmaD_pop_spec = 'fix',
                                  sigma_dmr_spec = "fix",
                                  dmr_mean_spec = 'fix',
                                  dmr_dev_spec = 'fix',
                                  ...,
                                  ln_F_mean_spec = "est") {

  messages_list <<- character(0) # string to attach to for printing messages # nolint: object_usage_linter.
  starting_values <- list(...)
  if(input_list$store_config) input_list$config$Setup_Mod_Catch_and_F <- mget(names(formals()))[-1]

  # Input Validation --------------------------------------------------------

  # Catch objects
  check_data_dimensions(
    ObsCatch,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_fish_fleets = input_list$data$n_fish_fleets,
    what = 'ObsCatch'
  )
  check_data_dimensions(
    UseCatch,
    n_regions = input_list$data$n_regions,
    n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas,
    n_fish_fleets = input_list$data$n_fish_fleets,
    what = 'UseCatch'
  )

  if(any(UseCatch_pop == 1)) {
    check_data_dimensions(
      ObsCatch_pop,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_fish_fleets = input_list$data$n_fish_fleets,
      what = 'ObsCatch_pop'
    )
    check_data_dimensions(
      UseCatch_pop,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_fish_fleets = input_list$data$n_fish_fleets,
      what = 'UseCatch_pop'
    )
  }

  if(any(UseDiscard == 1)) {
    check_data_dimensions(
      ObsDiscard,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_fish_fleets = input_list$data$n_fish_fleets,
      what = 'ObsDiscard'
    )
    check_data_dimensions(
      UseDiscard,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_fish_fleets = input_list$data$n_fish_fleets,
      what = 'UseDiscard'
    )
  }

  if(any(UseDiscard_pop == 1)) {
    check_data_dimensions(
      ObsDiscard_pop,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_fish_fleets = input_list$data$n_fish_fleets,
      what = 'ObsDiscard_pop'
    )
    check_data_dimensions(
      UseDiscard_pop,
      n_pop = input_list$data$n_pop,
      n_regions = input_list$data$n_regions,
      n_years = length(input_list$data$years),
      n_seas = input_list$data$n_seas,
      n_fish_fleets = input_list$data$n_fish_fleets,
      what = 'UseDiscard_pop'
    )
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

  # catch at age. a fleet fits either the aggregated catch with compositions or the catch at age,
  # never both: they are the same information stated twice
  input_list <- do_at_age_data_setup(input_list, ObsCatchAA, UseCatchAA, ObsCatchAA_SE,
                                     "CatchAA", "n_fish_fleets")
  input_list <- do_at_age_data_setup(input_list, ObsCatchAA_pop, UseCatchAA_pop, ObsCatchAA_pop_SE,
                                     "CatchAA", "n_fish_fleets", pop = TRUE)
  input_list <- do_at_age_data_setup(input_list, ObsDiscardAA, UseDiscardAA, ObsDiscardAA_SE,
                                     "DiscardAA", "n_fish_fleets")
  input_list <- do_at_age_data_setup(input_list, ObsDiscardAA_pop, UseDiscardAA_pop, ObsDiscardAA_pop_SE,
                                     "DiscardAA", "n_fish_fleets", pop = TRUE)

  use_catch_aa <- rep(0, input_list$data$n_fish_fleets)
  use_discard_aa <- rep(0, input_list$data$n_fish_fleets)
  for(f in 1:input_list$data$n_fish_fleets) {

    if(any(input_list$data$UseCatchAA[,,,,,f] == 1) ||
       any(input_list$data$UseCatchAA_pop[,,,,,,f] == 1)) {
      use_catch_aa[f] <- 1
      if(any(UseCatch[,,,f] == 1)) {
        stop("Fishery fleet ", f, " has both aggregated catch and catch at age in use. ",
             "A fleet fits one or the other: the two are the same information stated twice, ",
             "and the factorization of an at-age observation into a total and a composition ",
             "is exact for multinomial counts but not for the lognormal used here.")
      }
    }

    if(any(input_list$data$UseDiscardAA[,,,,,f] == 1) ||
       any(input_list$data$UseDiscardAA_pop[,,,,,,f] == 1)) {
      use_discard_aa[f] <- 1
      # a fraction is a property of the whole catch, not of one age, so the
      # at-age data source has discards in numbers or weight
      if(!discard_units[f] %in% c(0, 1)) {
        stop("Fishery fleet ", f, " fits discards at age with discard_units '",
             c("abd_frac", "biom_frac")[discard_units[f] - 1], "'. A discard fraction is a ",
             "property of the catch as a whole rather than of a single age. Use 'abd' or ",
             "'biom' for a fleet reporting discards at age.")
      }
    }
  } # end f loop

  input_list$data$use_catch_aa <- use_catch_aa
  input_list$data$use_discard_aa <- use_discard_aa

  # whether each data source is fit once a season or once a year as a season total
  n_fish <- input_list$data$n_fish_fleets
  input_list$data$Catch_seas_Type <- parse_seas_agg_spec(Catch_seas_Type, "Catch_seas_Type", n_fish)
  input_list$data$Catch_pop_seas_Type <- parse_seas_agg_spec(Catch_pop_seas_Type, "Catch_pop_seas_Type", n_fish)
  input_list$data$Discard_seas_Type <- parse_seas_agg_spec(Discard_seas_Type, "Discard_seas_Type", n_fish)
  input_list$data$Discard_pop_seas_Type <- parse_seas_agg_spec(Discard_pop_seas_Type, "Discard_pop_seas_Type", n_fish)
  input_list$data$CatchAA_seas_Type <- parse_seas_agg_spec(CatchAA_seas_Type, "CatchAA_seas_Type", n_fish)
  input_list$data$CatchAA_pop_seas_Type <- parse_seas_agg_spec(CatchAA_pop_seas_Type, "CatchAA_pop_seas_Type", n_fish)
  input_list$data$DiscardAA_seas_Type <- parse_seas_agg_spec(DiscardAA_seas_Type, "DiscardAA_seas_Type", n_fish)
  input_list$data$DiscardAA_pop_seas_Type <- parse_seas_agg_spec(DiscardAA_pop_seas_Type, "DiscardAA_pop_seas_Type", n_fish)

  # an aggregated observation is compared against the whole year, so only one season may be fit
  check_seas_agg_use(UseCatch, input_list$data$Catch_seas_Type, "UseCatch")
  check_seas_agg_use(UseCatch_pop, input_list$data$Catch_pop_seas_Type, "UseCatch_pop")
  check_seas_agg_use(UseDiscard, input_list$data$Discard_seas_Type, "UseDiscard")
  check_seas_agg_use(UseDiscard_pop, input_list$data$Discard_pop_seas_Type, "UseDiscard_pop")
  check_seas_agg_use(input_list$data$UseCatchAA, input_list$data$CatchAA_seas_Type, "UseCatchAA")
  check_seas_agg_use(input_list$data$UseCatchAA_pop, input_list$data$CatchAA_pop_seas_Type, "UseCatchAA_pop")
  check_seas_agg_use(input_list$data$UseDiscardAA, input_list$data$DiscardAA_seas_Type, "UseDiscardAA")
  check_seas_agg_use(input_list$data$UseDiscardAA_pop, input_list$data$DiscardAA_pop_seas_Type, "UseDiscardAA_pop")

  for(f in 1:n_fish) {
    if(input_list$data$Catch_seas_Type[f] == 1) collect_message("Aggregate catch for fishery fleet ", f, " is fit as a season total")
    if(input_list$data$CatchAA_seas_Type[f] == 1) collect_message("Catch at age for fishery fleet ", f, " is fit as a season total")
    if(input_list$data$Discard_seas_Type[f] == 1) collect_message("Discards for fishery fleet ", f, " is fit as a season total")
    if(input_list$data$DiscardAA_seas_Type[f] == 1) collect_message("Discards at age for fishery fleet ", f, " is fit as a season total")
  } # end f loop

  # how each data source is reported, what density it uses, and where its error comes from
  input_list <- do_at_age_type_setup(input_list, CatchAA_Type, "CatchAA", "n_fish_fleets", "UseCatchAA")
  input_list <- do_at_age_type_setup(input_list, CatchAA_pop_Type, "CatchAA", "n_fish_fleets", "UseCatchAA_pop", pop = TRUE)
  input_list <- do_at_age_type_setup(input_list, DiscardAA_Type, "DiscardAA", "n_fish_fleets", "UseDiscardAA")
  input_list <- do_at_age_type_setup(input_list, DiscardAA_pop_Type, "DiscardAA", "n_fish_fleets", "UseDiscardAA_pop", pop = TRUE)

  input_list <- do_at_age_like_setup(input_list, CatchAA_LikeType, CatchAA_sigma_form, "CatchAA", "n_fish_fleets")
  input_list <- do_at_age_like_setup(input_list, CatchAA_pop_LikeType, CatchAA_pop_sigma_form, "CatchAA", "n_fish_fleets", pop = TRUE)
  input_list <- do_at_age_like_setup(input_list, DiscardAA_LikeType, DiscardAA_sigma_form, "DiscardAA", "n_fish_fleets")
  input_list <- do_at_age_like_setup(input_list, DiscardAA_pop_LikeType, DiscardAA_pop_sigma_form, "DiscardAA", "n_fish_fleets", pop = TRUE)

  input_list <- do_age_corr_setup(input_list, AgeObsCorr_catch, "catch", "n_fish_fleets",
                                  "UseCatchAA", starting_values, rho_catch_spec)
  input_list <- do_age_corr_setup(input_list, AgeObsCorr_catch_pop, "catch", "n_fish_fleets",
                                  "UseCatchAA_pop", starting_values, rho_catch_pop_spec, pop = TRUE)
  input_list <- do_age_corr_setup(input_list, AgeObsCorr_discard, "discard", "n_fish_fleets",
                                  "UseDiscardAA", starting_values, rho_discard_spec)
  input_list <- do_age_corr_setup(input_list, AgeObsCorr_discard_pop, "discard", "n_fish_fleets",
                                  "UseDiscardAA_pop", starting_values, rho_discard_pop_spec, pop = TRUE)

  input_list <- do_key_mapping(input_list, sigmaCAA_key,
                               at_age_sigma_spec(sigmaCAA_spec, CatchAA_sigma_form, any(use_catch_aa == 1)),
                               "ln_sigmaCAA", "n_fish_fleets", "UseCatchAA", starting_values)
  input_list <- do_key_mapping(input_list, sigmaCAA_pop_key,
                               at_age_sigma_spec(sigmaCAA_pop_spec, CatchAA_pop_sigma_form, any(input_list$data$UseCatchAA_pop == 1)),
                               "ln_sigmaCAA_pop", "n_fish_fleets", "UseCatchAA_pop", starting_values, pop = TRUE)
  input_list <- do_key_mapping(input_list, sigmaDAA_key,
                               at_age_sigma_spec(sigmaDAA_spec, DiscardAA_sigma_form, any(use_discard_aa == 1)),
                               "ln_sigmaDAA", "n_fish_fleets", "UseDiscardAA", starting_values)
  input_list <- do_key_mapping(input_list, sigmaDAA_pop_key,
                               at_age_sigma_spec(sigmaDAA_pop_spec, DiscardAA_pop_sigma_form, any(input_list$data$UseDiscardAA_pop == 1)),
                               "ln_sigmaDAA_pop", "n_fish_fleets", "UseDiscardAA_pop", starting_values, pop = TRUE)

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

  # a fixed zero mean with a penalty centered on that mean shrinks the deviations toward F = 1,
  # a prior nobody intends, so the combination is rejected rather than warned about
  if(ln_F_mean_spec == "fix" && Use_F_pen == 1 && Fdev_pen_center == "fixed" && Fdev_model %in% c("iid", "ar1"))
    stop("ln_F_mean_spec = 'fix' with a zero-centered '", Fdev_model, "' penalty shrinks the deviations toward F = 1. Pair it with Fdev_pen_center = 'own_mean', Fdev_model = 'rw', or Use_F_pen = 0.")

  # under own-mean centering the penalty no longer sets the level of log F, so ln_F_mean and the
  # deviations' mean trade along a flat ridge unless "prop" initialization F reads the mean
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
  input_list$par$ln_sigmaC <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_sigmaC <- use_starting_value(input_list$par$ln_sigmaC, starting_values, "ln_sigmaC")

  input_list$par$ln_sigmaC_pop <- array(log(0.01), dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_sigmaC_pop <- use_starting_value(input_list$par$ln_sigmaC_pop, starting_values, "ln_sigmaC_pop")

  # Process error fishing deviations
  input_list$par$ln_sigmaF <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_sigmaF <- use_starting_value(input_list$par$ln_sigmaF, starting_values, "ln_sigmaF")

  # AR1 correlation for fishing mortality deviations (only used when Fdev_model = "ar1")
  input_list$par$Fdev_rho <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$Fdev_rho <- use_starting_value(input_list$par$Fdev_rho, starting_values, "Fdev_rho")

  # Log mean fishing mortality. A fixed mean defaults to zero so the deviations
  # are log F outright.
  input_list$par$ln_F_mean <- array(if(ln_F_mean_spec == "fix") 0 else log(0.1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_F_mean <- use_starting_value(input_list$par$ln_F_mean, starting_values, "ln_F_mean")

  # Log fishing deviations
  input_list$par$ln_F_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_F_devs <- use_starting_value(input_list$par$ln_F_devs, starting_values, "ln_F_devs")

  # Discard observation error
  input_list$par$ln_sigmaD <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_sigmaD <- use_starting_value(input_list$par$ln_sigmaD, starting_values, "ln_sigmaD")

  input_list$par$ln_sigmaD_pop <- array(log(0.01), dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_sigmaD_pop <- use_starting_value(input_list$par$ln_sigmaD_pop, starting_values, "ln_sigmaD_pop")

  # Process error discard deviations
  input_list$par$ln_sigma_dmr <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$ln_sigma_dmr <- use_starting_value(input_list$par$ln_sigma_dmr, starting_values, "ln_sigma_dmr")

  # Logit mean discard mortality
  input_list$par$logit_dmr_mean <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$logit_dmr_mean <- use_starting_value(input_list$par$logit_dmr_mean, starting_values, "logit_dmr_mean")

  # Logit discard mortality deviations
  input_list$par$logit_dmr_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))
  input_list$par$logit_dmr_devs <- use_starting_value(input_list$par$logit_dmr_devs, starting_values, "logit_dmr_devs")

  # Mapping Options ---------------------------------------------------------

  # Retained Catch Stuff
  input_list <- do_sigmaC_mapping(input_list, sigmaC_spec)
  input_list <- do_sigmaC_pop_mapping(input_list, sigmaC_pop_spec)
  input_list <- do_sigmaF_mapping(input_list, sigmaF_spec)
  input_list <- do_Fdev_rho_mapping(input_list, Fdev_rho_spec)
  input_list <- do_Fmort_mapping(input_list)
  # Free log-F parameterization: the mean is fixed at its starting value and the
  # deviations have all of log F.
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
