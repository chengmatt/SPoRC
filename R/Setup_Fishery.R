#' Set up fishing process inputs for the operating model simulation
#'
#' Populates \code{sim_list} with all fishery-related inputs needed by the
#' operating model: fishing mortality schedules, selectivity, catchability,
#' catch observation error, fishery indices, and age and length composition
#' simulation settings (likelihoods, sample sizes, overdispersion, and
#' correlation parameters). Must be called after \code{\link{Setup_Sim_Dim}}.
#'
#' Most array arguments accept either numeric codes or their character string
#' equivalents (e.g., \code{"biom"} instead of \code{1}); these are converted
#' internally via \code{convert_to_numeric()}.
#'
#' @param sim_list Simulation list returned by \code{\link{Setup_Sim_Dim}}.
#'
#' @section Fishing mortality:
#' \describe{
#'   \item{\code{Fmort_input}}{Instantaneous fishing mortality array
#'     \code{[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]}.
#'     Default: \code{0.1} for all cells.}
#'   \item{\code{init_F_val}}{Numeric vector of length \code{n_seas} giving the
#'     initial fishing mortality used during equilibrium initialisation for the
#'     dominant fleet. Default: \code{rep(0, n_seas)} (unfished initialisation).}
#'   \item{\code{fish_sel_input}}{Fishery selectivity-at-age array
#'     \code{[n_regions × n_yrs × n_ages × n_sexes × n_fish_fleets × n_sims]}.
#'     No default; must be provided.}
#'   \item{\code{fish_q_input}}{Fishery catchability array
#'     \code{[n_regions × n_yrs × n_fish_fleets × n_sims]}.
#'     Default: \code{1} for all cells.}
#'   \item{\code{catch_units}}{Integer array \code{[n_fish_fleets]} specifying
#'     the units in which catch is recorded per fleet.
#'     \code{0}/\code{"abd"} = abundance; \code{1}/\code{"biom"} = biomass
#'     (default).}
#' }
#'
#' @section Catch observation error:
#' \describe{
#'   \item{\code{ln_sigmaC}}{Log-scale standard deviation of lognormal catch
#'     observation error, array
#'     \code{[n_regions × n_yrs × n_seas × n_fish_fleets]}.
#'     Default: \code{log(0.02)}.}
#' }
#'
#' @section Fishery index:
#' \describe{
#'   \item{\code{ObsFishIdx_SE}}{Standard error of the observed fishery index
#'     (lognormal scale), array
#'     \code{[n_regions × n_yrs × n_seas × n_fish_fleets]}.
#'     Default: \code{0.2}.}
#'   \item{\code{fish_idx_type}}{Integer array \code{[n_regions × n_fish_fleets]}
#'     specifying whether each fleet's index is an abundance or biomass index.
#'     \code{0}/\code{"abd"} = abundance; \code{1}/\code{"biom"} = biomass
#'     (default).}
#' }
#'
#' @section Fishery age composition:
#' \describe{
#'   \item{\code{FishAgeComps_Type}}{Integer array \code{[n_yrs × n_fish_fleets]}
#'     controlling how age compositions are structured before simulation.
#'     \code{0}/\code{"agg"} = aggregated across regions and sexes;
#'     \code{1}/\code{"spltRspltS"} = split by region and sex;
#'     \code{2}/\code{"spltRjntS"} = split by region, joint across sexes (default);
#'     \code{999}/\code{"none"} = not simulated.}
#'   \item{\code{comp_fishage_like}}{Integer vector \code{[n_fish_fleets]}
#'     specifying the composition likelihood used to simulate age observations.
#'     \code{0}/\code{"Multinomial"} (default);
#'     \code{1}/\code{"Dirichlet-Multinomial"};
#'     \code{2}/\code{"iid-Logistic-Normal"};
#'     \code{3}/\code{"1d-Logistic-Normal"} (AR1 by age);
#'     \code{4}/\code{"2d-Logistic-Normal"} (AR1 by age, constant by sex).}
#'   \item{\code{ISS_FishAgeComps}}{Input sample sizes for age composition
#'     simulation, array
#'     \code{[n_regions × n_yrs × n_seas × n_sexes × n_fish_fleets × n_sims]}.
#'     Interpreted as the Dirichlet-Multinomial or multinomial sample size
#'     depending on \code{comp_fishage_like}. Default: \code{100}.}
#'   \item{\code{ln_FishAge_theta}}{Log-scale overdispersion parameters for
#'     fleet- region- and sex-specific compositions, array
#'     \code{[n_regions × n_sexes × n_fish_fleets]}.
#'     Only used for Dirichlet-Multinomial and logistic-normal likelihoods.
#'     Default: \code{log(1)}.}
#'   \item{\code{ln_FishAge_theta_agg}}{Log-scale overdispersion parameters
#'     for aggregated compositions, vector \code{[n_fish_fleets]}.
#'     Default: \code{log(1)}.}
#'   \item{\code{FishAge_corr_pars}}{Correlation parameters for 1D/2D
#'     logistic-normal likelihoods, array
#'     \code{[n_regions × n_sexes × n_fish_fleets × 2]}.
#'     The two trailing elements correspond to the age AR1 coefficient and
#'     the sex correlation. Default: \code{0.01}.}
#'   \item{\code{FishAge_corr_pars_agg}}{Correlation parameters for aggregated
#'     compositions under options 3–4, vector \code{[n_fish_fleets]}.
#'     Default: \code{0.01}.}
#' }
#'
#' @section Fishery length composition:
#' \describe{
#'   \item{\code{FishLenComps_Type}}{Integer array \code{[n_yrs × n_fish_fleets]}
#'     with the same coding as \code{FishAgeComps_Type} but applied to length
#'     compositions. Default: \code{2} (split by region, joint across sexes).}
#'   \item{\code{comp_fishlen_like}}{Integer vector \code{[n_fish_fleets]}
#'     with the same likelihood codes as \code{comp_fishage_like} applied to
#'     length compositions. Default: \code{0} (Multinomial).}
#'   \item{\code{ISS_FishLenComps}}{Input sample sizes for length composition
#'     simulation, array
#'     \code{[n_regions × n_yrs × n_seas × n_sexes × n_fish_fleets × n_sims]}.
#'     Default: \code{100}.}
#'   \item{\code{ln_FishLen_theta}}{Log-scale overdispersion parameters for
#'     length compositions, array
#'     \code{[n_regions × n_sexes × n_fish_fleets]}.
#'     Default: \code{log(1)}.}
#'   \item{\code{ln_FishLen_theta_agg}}{Log-scale overdispersion for aggregated
#'     length compositions, vector \code{[n_fish_fleets]}.
#'     Default: \code{log(1)}.}
#'   \item{\code{FishLen_corr_pars}}{Correlation parameters for 1D/2D
#'     logistic-normal length likelihoods, array
#'     \code{[n_regions × n_sexes × n_fish_fleets × 2]}.
#'     Default: \code{0.01}.}
#'   \item{\code{FishLen_corr_pars_agg}}{Correlation parameters for aggregated
#'     length compositions, vector \code{[n_fish_fleets]}.
#'     Default: \code{0.01}.}
#' }
#'
#' @return The input \code{sim_list} with all fishery fields appended under their
#'   respective names (e.g., \code{$Fmort}, \code{$fish_sel}, \code{$ln_sigmaC},
#'   \code{$comp_fishage_like}, \code{$ISS_FishAgeComps}, etc.). Character-coded
#'   inputs are converted to their integer equivalents before storage.
#'
#'
#' @export Setup_Sim_Fishing
#' @family Simulation Setup
Setup_Sim_Fishing <- function(sim_list,
                              ln_sigmaC = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              catch_units = array(1, dim = c(sim_list$n_fish_fleets)),
                              init_F_val = rep(0, sim_list$n_seas),
                              Fmort_input = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets, sim_list$n_sims)),
                              fish_sel_input,
                              fish_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ObsFishIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
                              fish_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
                              comp_fishage_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishAge_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
                              comp_fishlen_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishLen_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets))
                              ) {

  # Convert character inputs to numeric codes
  catch_units <- convert_to_numeric(catch_units,  list(abd = 0, biom = 1))
  fish_idx_type <- convert_to_numeric(fish_idx_type, list(abd = 0, biom = 1))
  comp_fishage_like <- convert_to_numeric(comp_fishage_like, list(Multinomial = 0,  `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  comp_fishlen_like <- convert_to_numeric(comp_fishlen_like, list(Multinomial = 0, `Dirichlet-Multinomial` = 1, `iid-Logistic-Normal` = 2, `1d-Logistic-Normal` = 3, `2d-Logistic-Normal` = 4))
  FishAgeComps_Type <- convert_to_numeric(FishAgeComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))
  FishLenComps_Type <- convert_to_numeric(FishLenComps_Type,  list(agg = 0, spltRspltS = 1, spltRjntS = 2, none = 999))

  # Validate dimensions of all input parameters
  check_sim_dimensions(ln_sigmaC, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_sigmaC")
  check_sim_dimensions(Fmort_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "Fmort_input")
  check_sim_dimensions(fish_sel_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "fish_sel_input")
  check_sim_dimensions(fish_q_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "fish_q_input")
  check_sim_dimensions(ObsFishIdx_SE, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ObsFishIdx_SE")

  # Validate fishery age composition parameters
  check_sim_dimensions(comp_fishage_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_like")
  check_sim_dimensions(ISS_FishAgeComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_FishAgeComps")
  check_sim_dimensions(ln_FishAge_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_theta")
  check_sim_dimensions(ln_FishAge_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishAge_theta_agg")
  check_sim_dimensions(FishAge_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_corr_pars_agg")
  check_sim_dimensions(FishAge_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "FishAge_corr_pars")
  check_sim_dimensions(FishAgeComps_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets,
                       what = "FishAgeComps_Type")

  # Validate fishery length composition parameters
  check_sim_dimensions(comp_fishlen_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishlen_like")
  check_sim_dimensions(ISS_FishLenComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs, n_seas = sim_list$n_seas,
                       n_sexes = sim_list$n_sexes, n_fish_fleets = sim_list$n_fish_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_FishLenComps")
  check_sim_dimensions(ln_FishLen_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_theta")
  check_sim_dimensions(ln_FishLen_theta_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "ln_FishLen_theta_agg")
  check_sim_dimensions(FishLen_corr_pars_agg, n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_corr_pars_agg")
  check_sim_dimensions(FishLen_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "FishLen_corr_pars")
  check_sim_dimensions(FishLenComps_Type, n_years = sim_list$n_yrs, n_fish_fleets = sim_list$n_fish_fleets,
                       what = "FishLenComps_Type")

  # output variables into list
  sim_list$Fmort <- Fmort_input # input fishing mortality pattern
  sim_list$catch_units <- catch_units # catch units
  sim_list$ln_sigmaC <- ln_sigmaC # Observation sd for catch
  sim_list$init_F <- init_F_val # initial F value
  sim_list$fish_sel <- fish_sel_input # fishery selectivity
  sim_list$fish_q <- fish_q_input # fishery catchability
  sim_list$ObsFishIdx_SE <- ObsFishIdx_SE # fishery index SE
  sim_list$fish_idx_type <- fish_idx_type # fishery index type

  # Fishery age compositions
  sim_list$comp_fishage_like <- comp_fishage_like
  sim_list$ISS_FishAgeComps <- ISS_FishAgeComps
  sim_list$ln_FishAge_theta <- ln_FishAge_theta
  sim_list$ln_FishAge_theta_agg <- ln_FishAge_theta_agg
  sim_list$FishAge_corr_pars_agg <- FishAge_corr_pars_agg
  sim_list$FishAge_corr_pars <- FishAge_corr_pars
  sim_list$FishAgeComps_Type <- FishAgeComps_Type

  # Fishery length compositions
  sim_list$comp_fishlen_like <- comp_fishlen_like
  sim_list$ISS_FishLenComps <- ISS_FishLenComps
  sim_list$ln_FishLen_theta <- ln_FishLen_theta
  sim_list$ln_FishLen_theta_agg <- ln_FishLen_theta_agg
  sim_list$FishLen_corr_pars_agg <- FishLen_corr_pars_agg
  sim_list$FishLen_corr_pars <- FishLen_corr_pars
  sim_list$FishLenComps_Type <- FishLenComps_Type

  return(sim_list)
}



#' Map sigma_F (fishing mortality process error SD) parameters
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
  map_sigmaF <- input_list$par$ln_sigmaF
  n_regions <- input_list$data$n_regions
  n_seas <- input_list$data$n_seas
  n_fish_fleets <- input_list$data$n_fish_fleets

  shared_specs <- c("est_shared_r", "est_shared_seas", "est_shared_f",
                    "est_shared_r_seas", "est_shared_r_f", "est_shared_seas_f",
                    "est_shared_r_seas_f")

  # In do_sigmaF_mapping, after defining shared_specs:
  valid_specs <- c(shared_specs, "fix", "est_all")
  if(!sigmaF_spec %in% valid_specs) stop("sigmaF_spec '", sigmaF_spec, "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  if(sigmaF_spec %in% shared_specs) {

    counter <- 1
    for(r in 1:n_regions) {
      for(seas in 1:n_seas) {
        for(f in 1:n_fish_fleets) {

          # --- Single dimension sharing ---
          if(sigmaF_spec == "est_shared_r" && r == 1) {
            map_sigmaF[, seas, f] <- counter; counter <- counter + 1
          }
          if(sigmaF_spec == "est_shared_seas" && seas == 1) {
            map_sigmaF[r, , f] <- counter; counter <- counter + 1
          }
          if(sigmaF_spec == "est_shared_f" && f == 1) {
            map_sigmaF[r, seas, ] <- counter; counter <- counter + 1
          }

          # --- Two dimension sharing ---
          if(sigmaF_spec == "est_shared_r_seas" && r == 1 && seas == 1) {
            map_sigmaF[, , f] <- counter; counter <- counter + 1
          }
          if(sigmaF_spec == "est_shared_r_f" && r == 1 && f == 1) {
            map_sigmaF[, seas, ] <- counter; counter <- counter + 1
          }
          if(sigmaF_spec == "est_shared_seas_f" && seas == 1 && f == 1) {
            map_sigmaF[r, , ] <- counter; counter <- counter + 1
          }

          # --- Three dimension sharing (single parameter) ---
          if(sigmaF_spec == "est_shared_r_seas_f" && r == 1 && seas == 1 && f == 1) {
            map_sigmaF[, , ] <- counter; counter <- counter + 1
          }

        } # end f
      } # end seas
    } # end r

    input_list$map$ln_sigmaF <- factor(map_sigmaF)
  }

  # Fixing sigmaF
  if(sigmaF_spec == "fix") input_list$map$ln_sigmaF <- factor(rep(NA, length(input_list$par$ln_sigmaF)))
  # Estimating all sigmaF
  if(sigmaF_spec == "est_all") input_list$map$ln_sigmaF <- factor(1:length(input_list$par$ln_sigmaF))

  # Print Message
  collect_message("sigmaF is specified as: ", sigmaF_spec)


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
  map_sigmaC <- input_list$par$ln_sigmaC
  n_regions <- input_list$data$n_regions
  n_years <- length(input_list$data$years)
  n_seas <- input_list$data$n_seas
  n_fish_fleets <- input_list$data$n_fish_fleets

  shared_specs <- c("est_shared_r", "est_shared_y", "est_shared_seas", "est_shared_f",
                    "est_shared_r_y", "est_shared_r_seas", "est_shared_r_f",
                    "est_shared_y_seas", "est_shared_y_f", "est_shared_seas_f",
                    "est_shared_r_y_seas", "est_shared_r_y_f", "est_shared_r_seas_f", "est_shared_y_seas_f",
                    "est_shared_r_y_seas_f")

  # In do_sigmaF_mapping, after defining shared_specs:
  valid_specs <- c(shared_specs, "fix", "est_all")
  if(!sigmaC_spec %in% valid_specs) stop("sigmaF_spec '", sigmaC_spec, "' not recognized. Valid options: ", paste(valid_specs, collapse = ", "))

  if(sigmaC_spec %in% shared_specs) {

    counter <- 1
    for(r in 1:n_regions) {
      for(y in 1:n_years) {
        for(seas in 1:n_seas) {
          for(f in 1:n_fish_fleets) {

            # --- Single dimension sharing ---
            if(sigmaC_spec == "est_shared_r" && r == 1) {
              map_sigmaC[, y, seas, f] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_y" && y == 1) {
              map_sigmaC[r, , seas, f] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_seas" && seas == 1) {
              map_sigmaC[r, y, , f] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_f" && f == 1) {
              map_sigmaC[r, y, seas, ] <- counter; counter <- counter + 1
            }

            # --- Two dimension sharing ---
            if(sigmaC_spec == "est_shared_r_y" && r == 1 && y == 1) {
              map_sigmaC[, , seas, f] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_r_seas" && r == 1 && seas == 1) {
              map_sigmaC[, y, , f] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_r_f" && r == 1 && f == 1) {
              map_sigmaC[, y, seas, ] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_y_seas" && y == 1 && seas == 1) {
              map_sigmaC[r, , , f] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_y_f" && y == 1 && f == 1) {
              map_sigmaC[r, , seas, ] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_seas_f" && seas == 1 && f == 1) {
              map_sigmaC[r, y, , ] <- counter; counter <- counter + 1
            }

            # --- Three dimension sharing ---
            if(sigmaC_spec == "est_shared_r_y_seas" && r == 1 && y == 1 && seas == 1) {
              map_sigmaC[, , , f] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_r_y_f" && r == 1 && y == 1 && f == 1) {
              map_sigmaC[, , seas, ] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_r_seas_f" && r == 1 && seas == 1 && f == 1) {
              map_sigmaC[, y, , ] <- counter; counter <- counter + 1
            }
            if(sigmaC_spec == "est_shared_y_seas_f" && y == 1 && seas == 1 && f == 1) {
              map_sigmaC[r, , , ] <- counter; counter <- counter + 1
            }

            # --- Four dimension sharing (single parameter) ---
            if(sigmaC_spec == "est_shared_r_y_seas_f" && r == 1 && y == 1 && seas == 1 && f == 1) {
              map_sigmaC[, , , ] <- counter; counter <- counter + 1
            }

          } # end f
        } # end seas
      } # end y
    } # end r

    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }

  # Fixing sigmaC
  if(sigmaC_spec == "fix") input_list$map$ln_sigmaC <- factor(rep(NA, length(input_list$par$ln_sigmaC)))
  # Estimating all sigmaC
  if(sigmaC_spec == "est_all") input_list$map$ln_sigmaC <- factor(1:length(input_list$par$ln_sigmaC))

  # Print Message
  collect_message("sigmaC is specified as: ", sigmaC_spec)


  return(input_list)
}

#' Map fishing mortality deviation parameters
#'
#' Constructs the \code{ln_F_devs} factor map, assigning unique estimation
#' indices to region–year–season–fleet cells where catch data are used
#' (\code{UseCatch == 1}) and mapping cells without catch data to \code{NA}.
#' This ensures that \code{ln_F_devs} parameters are only estimated for
#' dimensions with observed catch.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{$data$UseCatch} to be populated by
#'   \code{\link{Setup_Mod_Catch_and_F}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_F_devs} set to a
#'   factor vector. Cells with catch are assigned sequential integer indices;
#'   cells without catch are \code{NA}.
#'
#' @keywords internal
do_Fmort_mapping <- function(input_list) {

  # Mapping for fishing mortality deviations
  F_dev_map <- input_list$par$ln_F_devs # initialize for mapping
  F_dev_map[] <- NA
  F_dev_counter <- 1

  for(r in 1:input_list$data$n_regions) {
    for(y in 1:length(input_list$data$years)) {
      for(seas in 1:input_list$data$n_seas) {
        for(f in 1:input_list$data$n_fish_fleets) {

          # if no catch, don't estimate devs
          if(input_list$data$UseCatch[r,y,seas,f] == 0) {
            F_dev_map[r,y,seas,f] <- NA
          }

          # if have catch, estimate f devs
          if(input_list$data$UseCatch[r,y,seas,f] == 1) {
            F_dev_map[r,y,seas,f] <- F_dev_counter
            F_dev_counter <- F_dev_counter + 1
          } # end if

        } # end f
      } # end seas loop
    } # end y
  } # end r

  input_list$map$ln_F_devs <- factor(F_dev_map)

  return(input_list)
}

##' Set up fishing mortality and catch observation inputs
#'
#' Populates \code{input_list} with observed catch, catch usage indicators,
#' fishing mortality parameters (\code{ln_F_mean}, \code{ln_F_devs}), and
#' observation/process error structures (\code{ln_sigmaC}, \code{ln_sigmaF}).
#' Must be called after \code{\link{Setup_Mod_Biologicals}}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsCatch Observed catch array \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   Values should be in the units specified by \code{catch_units}.
#' @param UseCatch Binary indicator array \code{[n_regions × n_years × n_seas × n_fish_fleets]}
#'   controlling which catch observations enter the likelihood and whether
#'   \code{ln_F_devs} are estimated for each cell. \code{1} = use; \code{0} = exclude
#'   (corresponding \code{ln_F_devs} will be mapped to \code{NA}).
#' @param catch_units Character array \code{[n_fish_fleets]} specifying catch
#'   units per fleet. \code{"biom"} = biomass (default); \code{"abd"} = abundance.
#'   Converted internally to \code{0}/\code{1} integer codes.
#' @param Use_F_pen Integer flag for applying a fishing mortality penalty to
#'   penalise large deviations in \code{ln_F_devs}. \code{1} = apply (default);
#'   \code{0} = do not apply.
#' @param sigmaC_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaC} (catch observation error SD). Default \code{"fix"} holds
#'   \code{ln_sigmaC} at its starting value (\code{log(0.01)} unless overridden
#'   via \code{...}). See \code{\link{do_sigmaC_mapping}} for all sharing options
#'   (\code{"est_shared_r"}, \code{"est_shared_y"}, \code{"est_all"}, etc.).
#'   A warning is issued if \code{"fix"} is selected without providing a starting
#'   value in \code{...}.
#' @param sigmaF_spec Character string specifying the sharing structure for
#'   \code{ln_sigmaF} (fishing mortality process error SD). Default \code{"fix"}
#'   holds \code{ln_sigmaF} at its starting value (\code{log(1)}, i.e.,
#'   \eqn{\sigma_F = 1}, unless overridden via \code{...}). See
#'   \code{\link{do_sigmaF_mapping}} for all sharing options. A warning is issued
#'   if \code{"fix"} is selected without providing a starting value in \code{...}.
#' @param ... Optional starting value overrides, passed by name. Recognised
#'   arguments:
#'   \describe{
#'     \item{\code{ln_sigmaC}}{Array of log-scale starting values for catch
#'       observation error, dimensioned
#'       \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'       Default: \code{log(0.01)}.}
#'     \item{\code{ln_sigmaF}}{Array of log-scale starting values for fishing
#'       mortality process error, dimensioned
#'       \code{[n_regions × n_seas × n_fish_fleets]}.
#'       Default: \code{log(1)}.}
#'     \item{\code{ln_F_mean}}{Array of log-scale mean fishing mortality starting
#'       values, dimensioned \code{[n_regions × n_seas × n_fish_fleets]}.
#'       Default: \code{log(0.1)}.}
#'     \item{\code{ln_F_devs}}{Array of log-scale annual fishing mortality
#'       deviation starting values, dimensioned
#'       \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'       Default: \code{0}.}
#'   }
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated. Key additions: \code{$data$ObsCatch},
#'   \code{$data$UseCatch}, \code{$data$Use_F_pen}, \code{$data$catch_units},
#'   \code{$par$ln_sigmaC}, \code{$par$ln_sigmaF}, \code{$par$ln_F_mean},
#'   \code{$par$ln_F_devs}, \code{$map$ln_sigmaC}, \code{$map$ln_sigmaF},
#'   \code{$map$ln_F_devs}.
#'
#'
#' @export Setup_Mod_Catch_and_F
#' @family Model Setup
Setup_Mod_Catch_and_F <- function(input_list,
                                  ObsCatch,
                                  catch_units = array("biom", dim = c(input_list$data$n_fish_fleets)),
                                  UseCatch,
                                  Use_F_pen = 1,
                                  sigmaC_spec = "fix",
                                  sigmaF_spec = "fix",
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Catch objects
  check_data_dimensions(ObsCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsCatch')
  check_data_dimensions(UseCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseCatch')

  # Indicators for whether catch is aggregated across regions
  if(any(UseCatch == 0)) collect_message("User specified catch for some years, seasons, and fleets to not be fit to, and ln_F_devs will not be estimated for those dimensions")

  # Fishing Mortality checking
  if(!Use_F_pen %in% c(0,1)) stop("Use_F_pen incorrectly specified. Either set at 0 (don't use F penalty) or 1 (use F penalty)")
  else collect_message("Fishing mortality penalty is: ", ifelse(Use_F_pen == 0, 'Not Used', "Used"))

  if(sigmaC_spec == "fix" && !("ln_sigmaC" %in% names(starting_values))) warning("sigmaC is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaF_spec == "fix" && !("ln_sigmaF" %in% names(starting_values))) warning("sigmaF_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")

  # Catch units
  catch_units[catch_units == 'abd'] <- 0
  catch_units[catch_units == 'biom'] <- 1
  catch_units <- array(as.numeric(catch_units), dim = c(input_list$data$n_fish_fleets)) # convert to numeric array

  # Populate Data List ------------------------------------------------------

  input_list$data$ObsCatch <- ObsCatch
  input_list$data$UseCatch <- UseCatch
  input_list$data$Use_F_pen <- Use_F_pen
  input_list$data$catch_units <- catch_units

  # Populate Parameter List -------------------------------------------------

  # Catch observation error
  if("ln_sigmaC" %in% names(starting_values)) input_list$par$ln_sigmaC <- starting_values$ln_sigmaC
  else input_list$par$ln_sigmaC <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Process error fishing deviations for regional catch
  if("ln_sigmaF" %in% names(starting_values)) input_list$par$ln_sigmaF <- starting_values$ln_sigmaF
  else input_list$par$ln_sigmaF <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Log mean fishing mortality
  if("ln_F_mean" %in% names(starting_values)) input_list$par$ln_F_mean <- starting_values$ln_F_mean
  else input_list$par$ln_F_mean <- array(log(0.1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Log fishing deviations
  if("ln_F_devs" %in% names(starting_values)) input_list$par$ln_F_devs <- starting_values$ln_F_devs
  else input_list$par$ln_F_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------
  input_list <- do_sigmaC_mapping(input_list, sigmaC_spec)
  input_list <- do_sigmaF_mapping(input_list, sigmaF_spec)
  input_list <- do_Fmort_mapping(input_list)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Map fishery age composition overdispersion parameters
#'
#' Constructs factor maps for \code{ln_FishAge_theta} (fleet- region- and
#' sex-specific overdispersion) and \code{ln_FishAge_theta_agg} (aggregated
#' overdispersion) based on the composition type and likelihood specified in
#' \code{$data$FishAgeComps_Type} and \code{$data$FishAgeComps_LikeType}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' (\code{LikeType == 0}) or with no observed age compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishAgeComps_Type}, \code{FishAgeComps_LikeType},
#'   and \code{UseFishAgeComps} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishAge_theta} and
#'   \code{$map$ln_FishAge_theta_agg} set to factor vectors. Active parameters
#'   receive sequential integer indices; inactive parameters are \code{NA}.
#'
#' @keywords internal
do_FishAge_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishage_agg <- 1
  counter_fishage <- 1

  # initialize array to set up mapping
  map_FishAge_theta <- input_list$par$ln_FishAge_theta
  map_FishAge_theta_agg <- input_list$par$ln_FishAge_theta_agg
  map_FishAge_theta[] <- NA
  map_FishAge_theta_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # get unique fishery comp types
    fishage_comp_type <- unique(input_list$data$FishAgeComps_Type[,f])

    # If aggregated (ages)
    if(any(fishage_comp_type == 0) && input_list$data$FishAgeComps_LikeType[f] != 0) {
      map_FishAge_theta_agg[f] <- counter_fishage_agg
      counter_fishage_agg <- counter_fishage_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(fishage_comp_type == 1) && input_list$data$FishAgeComps_LikeType[f] != 0) {
          map_FishAge_theta[r,s,f] <- counter_fishage
          counter_fishage <- counter_fishage + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(fishage_comp_type == 2) && input_list$data$FishAgeComps_LikeType[f] != 0 && s == 1) {
          map_FishAge_theta[r,1,f] <- counter_fishage
          counter_fishage <- counter_fishage + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any age comps for a given fleet
    if(input_list$data$FishAgeComps_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps[,,,f]) == 0) {
      map_FishAge_theta[,,f] <- NA
      map_FishAge_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishAge_theta <- factor(map_FishAge_theta)
  input_list$map$ln_FishAge_theta_agg <- factor(map_FishAge_theta_agg)

  return(input_list)
}

#' Map fishery length composition overdispersion parameters
#'
#' Analogous to \code{\link{do_FishAge_theta_mapping}} but for length
#' compositions. Constructs factor maps for \code{ln_FishLen_theta} and
#' \code{ln_FishLen_theta_agg} based on \code{FishLenComps_Type},
#' \code{FishLenComps_LikeType}, and \code{UseFishLenComps}.
#' Parameters are mapped to \code{NA} for fleets using multinomial likelihoods
#' or with no observed length compositions.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists. Requires \code{FishLenComps_Type}, \code{FishLenComps_LikeType},
#'   and \code{UseFishLenComps} to be set by
#'   \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' @return The input \code{input_list} with \code{$map$ln_FishLen_theta} and
#'   \code{$map$ln_FishLen_theta_agg} set to factor vectors.
#'
#' @keywords internal
do_FishLen_theta_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_agg <- 1
  counter_fishlen <- 1

  # initialize array to set up mapping
  map_FishLen_theta <- input_list$par$ln_FishLen_theta
  map_FishLen_theta_agg <- input_list$par$ln_FishLen_theta_agg
  map_FishLen_theta[] <- NA
  map_FishLen_theta_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # get unique fishery comp types
    fishlen_comp_type <- unique(input_list$data$FishLenComps_Type[,f])

    # If aggregated (ages)
    if(any(fishlen_comp_type == 0) && input_list$data$FishLenComps_LikeType[f] != 0) {
      map_FishLen_theta_agg[f] <- counter_fishlen_agg
      counter_fishlen_agg <- counter_fishlen_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(fishlen_comp_type == 1) && input_list$data$FishLenComps_LikeType[f] != 0) {
          map_FishLen_theta[r,s,f] <- counter_fishlen
          counter_fishlen <- counter_fishlen + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(fishlen_comp_type == 2) && input_list$data$FishLenComps_LikeType[f] != 0 && s == 1) {
          map_FishLen_theta[r,1,f] <- counter_fishlen
          counter_fishlen <- counter_fishlen + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any lenght comps for a given fleet
    if(input_list$data$FishLenComps_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps[,,,f]) == 0) {
      map_FishLen_theta[,,f] <- NA
      map_FishLen_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishLen_theta <- factor(map_FishLen_theta)
  input_list$map$ln_FishLen_theta_agg <- factor(map_FishLen_theta_agg)

  return(input_list)
}

#' Map fishery age composition correlation parameters
#'
#' Constructs factor maps for \code{FishAge_corr_pars} (region- and sex-specific
#' AR1 and sex correlation parameters) and \code{FishAge_corr_pars_agg}
#' (aggregated correlation parameters) for 1D and 2D logistic-normal age
#' composition likelihoods. Parameters are activated only when
#' \code{FishAgeComps_LikeType %in% c(3, 4)} (1D or 2D logistic-normal);
#' all other likelihoods map correlation parameters to \code{NA}.
#'
#' For the 2D logistic-normal (\code{LikeType == 4}), both trailing elements of
#' the \code{[..., 2]} dimension are activated — element 1 for the age AR1
#' coefficient and element 2 for the sex correlation (skipped when
#' \code{n_sexes == 1}).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#'
#' @return The input \code{input_list} with \code{$map$FishAge_corr_pars} and
#'   \code{$map$FishAge_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishAge_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishage_corr <- 1
  counter_fishage_corr_agg <- 1

  # initialize array to set up mapping
  map_FishAge_corr_pars <- input_list$par$FishAge_corr_pars
  map_FishAge_corr_pars_agg <- input_list$par$FishAge_corr_pars_agg
  map_FishAge_corr_pars[] <- NA
  map_FishAge_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$FishAgeComps_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps[,,,f]) == 0) {
      map_FishAge_corr_pars[,,f,] <- NA
      map_FishAge_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique fishery comp types
    fishage_comp_type <- unique(input_list$data$FishAgeComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(fishage_comp_type == 0) && input_list$data$FishAgeComps_LikeType[f] != 0) {
      if(input_list$data$FishAgeComps_LikeType[f] == 3) {
        map_FishAge_corr_pars_agg[f] <- counter_fishage_corr_agg
        counter_fishage_corr_agg <- counter_fishage_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(fishage_comp_type == 1) && input_list$data$FishAgeComps_LikeType[f] != 0) {
          if(input_list$data$FishAgeComps_LikeType[f] == 3) {
            map_FishAge_corr_pars[r,s,f,1] <- counter_fishage_corr
            counter_fishage_corr <- counter_fishage_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(fishage_comp_type == 2) && input_list$data$FishAgeComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$FishAgeComps_LikeType[f] == 3) {
            map_FishAge_corr_pars[r,1,f,1] <- counter_fishage_corr
            counter_fishage_corr <- counter_fishage_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$FishAgeComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishAge_corr_pars[r,1,f,i] <- counter_fishage_corr
              counter_fishage_corr <- counter_fishage_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$FishAge_corr_pars_agg <- factor(map_FishAge_corr_pars_agg)
  input_list$map$FishAge_corr_pars <- factor(map_FishAge_corr_pars)

  return(input_list)
}

#' Map fishery length composition correlation parameters
#'
#' Analogous to \code{\link{do_FishAge_corr_pars_mapping}} but for length
#' compositions. Constructs factor maps for \code{FishLen_corr_pars} and
#' \code{FishLen_corr_pars_agg} for 1D and 2D logistic-normal length composition
#' likelihoods (\code{FishLenComps_LikeType %in% c(3, 4)}).
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#'
#' @return The input \code{input_list} with \code{$map$FishLen_corr_pars} and
#'   \code{$map$FishLen_corr_pars_agg} set to factor vectors.
#'
#' @keywords internal
do_FishLen_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_fishlen_corr <- 1
  counter_fishlen_corr_agg <- 1

  # initialize array to set up mapping
  map_FishLen_corr_pars <- input_list$par$FishLen_corr_pars
  map_FishLen_corr_pars_agg <- input_list$par$FishLen_corr_pars_agg
  map_FishLen_corr_pars[] <- NA
  map_FishLen_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$FishLenComps_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps[,,,f]) == 0) {
      map_FishLen_corr_pars[,,f,] <- NA
      map_FishLen_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique fishery comp types
    fishlen_comp_type <- unique(input_list$data$FishLenComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(fishlen_comp_type == 0) && input_list$data$FishLenComps_LikeType[f] != 0) {
      if(input_list$data$FishLenComps_LikeType[f] == 3) {
        map_FishLen_corr_pars_agg[f] <- counter_fishlen_corr_agg
        counter_fishlen_corr_agg <- counter_fishlen_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(fishlen_comp_type == 1) && input_list$data$FishLenComps_LikeType[f] != 0) {
          if(input_list$data$FishLenComps_LikeType[f] == 3) {
            map_FishLen_corr_pars[r,s,f,1] <- counter_fishlen_corr
            counter_fishlen_corr <- counter_fishlen_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(fishlen_comp_type == 2) && input_list$data$FishLenComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$FishLenComps_LikeType[f] == 3) {
            map_FishLen_corr_pars[r,1,f,1] <- counter_fishlen_corr
            counter_fishlen_corr <- counter_fishlen_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$FishLenComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_FishLen_corr_pars[r,1,f,i] <- counter_fishlen_corr
              counter_fishlen_corr <- counter_fishlen_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$FishLen_corr_pars_agg <- factor(map_FishLen_corr_pars_agg)
  input_list$map$FishLen_corr_pars <- factor(map_FishLen_corr_pars)

  return(input_list)
}

#' Set up fishery index, age composition, and length composition inputs
#'
#' Populates \code{input_list} with observed fishery indices, age compositions,
#' and length compositions along with their usage indicators, likelihood types,
#' composition structure types, input sample sizes, and overdispersion and
#' correlation parameter starting values and mappings. Must be called after
#' \code{\link{Setup_Mod_Catch_and_F}}.
#'
#' When \code{ISS_FishAgeComps} or \code{ISS_FishLenComps} are \code{NULL},
#' input sample sizes are derived automatically by summing the observed
#' composition arrays within each year–fleet–season–region cell, consistent
#' with the specified \code{FishAgeComps_Type} or \code{FishLenComps_Type}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists, as returned by upstream setup functions.
#' @param ObsFishIdx Observed fishery CPUE or biomass index array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#' @param ObsFishIdx_SE Standard errors of \code{ObsFishIdx} on the log scale,
#'   same dimensions as \code{ObsFishIdx}.
#' @param UseFishIdx Binary indicator array \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = include index in the likelihood; \code{0} = exclude.
#' @param fish_idx_type Character vector of length \code{n_fish_fleets} specifying
#'   the index type for each fleet. \code{"biom"} = biomass; \code{"abd"} =
#'   abundance; \code{"none"} = no index for this fleet.
#' @param ObsFishAgeComps Observed fishery age composition array
#'   \code{[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]}.
#'   Values may be raw counts or proportions; if proportions, supply
#'   \code{ISS_FishAgeComps} explicitly.
#' @param UseFishAgeComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit age compositions; \code{0} = exclude.
#' @param ISS_FishAgeComps Input sample size array
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), sample sizes are computed automatically by
#'   summing \code{ObsFishAgeComps} within each year–fleet–season–region cell
#'   according to \code{FishAgeComps_Type}.
#' @param ObsFishLenComps Observed fishery length composition array
#'   \code{[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]}.
#'   Only required when \code{input_list$data$fit_lengths == 1}.
#' @param UseFishLenComps Binary indicator array
#'   \code{[n_regions × n_years × n_seas × n_fish_fleets]}.
#'   \code{1} = fit length compositions; \code{0} = exclude.
#' @param ISS_FishLenComps Input sample size array for length compositions,
#'   \code{[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]}.
#'   If \code{NULL} (default), derived automatically from \code{ObsFishLenComps}.
#' @param FishAgeComps_LikeType Character vector of length \code{n_fish_fleets}
#'   specifying the likelihood for fishery age compositions. Options:
#'   \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#'   \code{"iid-Logistic-Normal"}, \code{"1d-Logistic-Normal"},
#'   \code{"2d-Logistic-Normal"}, \code{"none"}.
#' @param FishLenComps_LikeType Same as \code{FishAgeComps_LikeType} but for
#'   length compositions.
#' @param FishAgeComps_Type Character vector defining the age composition
#'   structure (aggregation level) for each fleet and time period. Each element
#'   must follow the format \code{"<type>_Year_<start>-<end>_Fleet_<f>"} or
#'   \code{"<type>_Year_<start>-terminal_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"agg"}}{Aggregated across regions and sexes
#'       (incompatible with \code{"2d-Logistic-Normal"}).}
#'     \item{\code{"spltRspltS"}}{Split by region and sex.}
#'     \item{\code{"spltRjntS"}}{Split by region, summed jointly across sexes.}
#'     \item{\code{"none"}}{No composition data for this fleet and period.}
#'   }
#'   Example: \code{c("spltRjntS_Year_1-10_Fleet_1", "agg_Year_11-terminal_Fleet_1")}.
#' @param FishLenComps_Type Same format and options as \code{FishAgeComps_Type}
#'   but applied to length compositions.
#' @param ... Optional starting value overrides for overdispersion and
#'   correlation parameters. Recognised arguments:
#'   \describe{
#'     \item{\code{ln_FishAge_theta}}{Log-scale overdispersion starting values
#'       \code{[n_regions × n_sexes × n_fish_fleets]}. Default: \code{0}.}
#'     \item{\code{ln_FishAge_theta_agg}}{Aggregated overdispersion starting
#'       values \code{[n_fish_fleets]}. Default: \code{0}.}
#'     \item{\code{FishAge_corr_pars}}{Correlation parameter starting values
#'       \code{[n_regions × n_sexes × n_fish_fleets × 2]}. Default: \code{0.01}.}
#'     \item{\code{FishAge_corr_pars_agg}}{Aggregated correlation starting values
#'       \code{[n_fish_fleets]}. Default: \code{0.01}.}
#'     \item{\code{ln_FishLen_theta}}{Length composition overdispersion
#'       \code{[n_regions × n_sexes × n_fish_fleets]}. Default: \code{0}.}
#'     \item{\code{ln_FishLen_theta_agg}}{Aggregated length overdispersion
#'       \code{[n_fish_fleets]}. Default: \code{0}.}
#'     \item{\code{FishLen_corr_pars}}{Length correlation parameters
#'       \code{[n_regions × n_sexes × n_fish_fleets × 2]}. Default: \code{0.01}.}
#'     \item{\code{FishLen_corr_pars_agg}}{Aggregated length correlation
#'       \code{[n_fish_fleets]}. Default: \code{0.01}.}
#'   }
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated with all fishery index and composition fields, including
#'   computed or supplied ISS arrays, integer-coded likelihood and composition
#'   type matrices, overdispersion parameters, and their factor maps.
#'
#'
#' @export Setup_Mod_FishIdx_and_Comps
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_FishIdx_and_Comps <- function(input_list,
                                        ObsFishIdx,
                                        ObsFishIdx_SE,
                                        fish_idx_type,
                                        UseFishIdx,
                                        ObsFishAgeComps,
                                        UseFishAgeComps,
                                        ISS_FishAgeComps = NULL,
                                        ObsFishLenComps,
                                        UseFishLenComps,
                                        ISS_FishLenComps = NULL,
                                        FishAgeComps_LikeType,
                                        FishLenComps_LikeType,
                                        FishAgeComps_Type,
                                        FishLenComps_Type,
                                        ...
                                        ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation ---------------------------------------------------------

  # Fishery indices
  check_data_dimensions(ObsFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx')
  check_data_dimensions(ObsFishIdx_SE, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_SE')
  check_data_dimensions(UseFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishIdx')
  if(!all(fish_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for fish_idx_type. Should be either abd, biom, or none")

  # Fishery compositions
  check_data_dimensions(ObsFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps')
  check_data_dimensions(UseFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps')
  check_data_dimensions(UseFishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_seas = input_list$data$n_seas, n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsFishLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps')
  if(!is.null(ISS_FishAgeComps)) check_data_dimensions(ISS_FishAgeComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps')
  if(!is.null(ISS_FishLenComps)) check_data_dimensions(ISS_FishLenComps, n_regions = input_list$data$n_regions, n_seas = input_list$data$n_seas, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps')
  check_data_dimensions(FishAgeComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_LikeType')
  check_data_dimensions(FishLenComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_LikeType')
  if(!all(FishAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")


  # Fishery Index Options ---------------------------------------------------

  fish_idx_type_vals <- array(NA, dim = c(input_list$data$n_fish_fleets))
  for(f in 1:input_list$data$n_fish_fleets) {
    if(fish_idx_type[f] == 'biom') fish_idx_type_vals[f] <- 1 # biomass
    if(fish_idx_type[f] == 'abd') fish_idx_type_vals[f] <- 0 # abundance
    if(fish_idx_type[f] == 'none') fish_idx_type_vals[f] <- 999 # none
    collect_message(paste("Fishery Index", "for fishery fleet", f, "specified as:" , fish_idx_type[f]))
  } # end f loop


  # Fishery Age Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishage_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishAgeComps_LikeType[f] == 'none') comp_fishage_like_vals <- c(comp_fishage_like_vals, 999)
    if(FishAgeComps_LikeType[f] == "Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 0)
    if(FishAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishage_like_vals <- c(comp_fishage_like_vals, 1)
    if(FishAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 2)
    if(FishAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 3)
    if(FishAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_fishage_like_vals <- c(comp_fishage_like_vals, 4)
    collect_message(paste("Fishery Age Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishAgeComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishAgeComps_Type)) {

    # Extract out components from list
    tmp <- FishAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishAgeComps_Type. This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # Composition type
    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishage_like_vals[fleet] == 4) stop("Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishAgeComps_Type_Mat))) stop("FishAgeComps_Type is returning an NA. Did you update the year range of FishAgeComps_Type?")

  # Fishery Length Composition Options -----------------------------------------

  # Specifying composition likelihood
  comp_fishlen_like_vals <- vector()
  for(f in 1:input_list$data$n_fish_fleets) {
    if(FishLenComps_LikeType[f] == 'none') comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 999)
    if(FishLenComps_LikeType[f] == "Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 0)
    if(FishLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 1)
    if(FishLenComps_LikeType[f] == "iid-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 2)
    if(FishLenComps_LikeType[f] == "1d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 3)
    if(FishLenComps_LikeType[f] == "2d-Logistic-Normal") comp_fishlen_like_vals <- c(comp_fishlen_like_vals, 4)
    collect_message(paste("Fishery Length Composition Likelihoods", "for fishery fleet", f, "specified as:" , FishLenComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  FishLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(FishLenComps_Type)) {

    # Extract out components from list
    tmp <- FishLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # define composition types
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("FishLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for FishLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

    # get year ranges
    if(!str_detect(tmp, "terminal")) { # if not terminal year
      year_range <- as.numeric(unlist(strsplit(tmp_vec[3], "-")))
      years <- year_range[1]:year_range[2] # get sequence of years
    } else { # if terminal year
      year_range <- unlist(strsplit(tmp_vec[3], '-'))[1] # get year range
      years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
    }

    # define composition types
    if(comps_type_tmp == "agg") {
      if(comp_fishlen_like_vals[fleet] == 4) stop("Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    FishLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(FishLenComps_Type_Mat))) stop("FishLenComps_Type_Mat is returning an NA. Did you update the year range of FishLenComps_Type_Mat?")


  # ISS Munging -------------------------------------------------------------

  # Fishery Ages
  if(is.null(ISS_FishAgeComps)) {
    collect_message("No ISS is specified for FishAgeComps. ISS weighting is calculated by summing up values from ObsFishAgeComps each year")
    ISS_FishAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0) or joint across sexes
          if(FishAgeComps_Type_Mat[y,f] == 0) ISS_FishAgeComps[1,y,seas,1,f] <- sum(ObsFishAgeComps[,y,seas,,,f])
          # if split by region and sex
          if(FishAgeComps_Type_Mat[y,f] == 1) ISS_FishAgeComps[,y,seas,,f] <- apply(ObsFishAgeComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishAgeComps_Type_Mat[y,f] == 2) ISS_FishAgeComps[,y,seas,1,f] <- apply(ObsFishAgeComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps)) {
    collect_message("No ISS is specified for FishLenComps. ISS weighting is calculated by summing up values from ObsFishLenComps each year")
    ISS_FishLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        for(seas in 1:input_list$data$n_seas) {
          # if aggregated across sexes and regions (0)
          if(FishLenComps_Type_Mat[y,f] == 0) ISS_FishLenComps[1,y,seas,1,f] <- sum(ObsFishLenComps[,y,seas,,,f])
          # if split by region and sex
          if(FishLenComps_Type_Mat[y,f] == 1) ISS_FishLenComps[,y,seas,,f] <- apply(ObsFishLenComps[,y,seas,,,f, drop = FALSE], c(1,4), sum)
          # if split by region, joint by sex
          if(FishLenComps_Type_Mat[y,f] == 2) ISS_FishLenComps[,y,seas,1,f] <- apply(ObsFishLenComps[,y,seas,,,f, drop = FALSE], 1, sum)
        } # end seas loop
      } # end f loop
    } # end y loop
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_FishAgeComps <- ISS_FishAgeComps
  input_list$data$ISS_FishLenComps <- ISS_FishLenComps
  input_list$data$ObsFishIdx <- ObsFishIdx
  input_list$data$ObsFishIdx_SE <- ObsFishIdx_SE
  input_list$data$UseFishIdx <- UseFishIdx
  input_list$data$fish_idx_type <- fish_idx_type_vals
  input_list$data$ObsFishAgeComps <- ObsFishAgeComps
  input_list$data$UseFishAgeComps <- UseFishAgeComps
  input_list$data$ObsFishLenComps <- ObsFishLenComps
  input_list$data$UseFishLenComps <- UseFishLenComps
  input_list$data$FishAgeComps_LikeType <- comp_fishage_like_vals
  input_list$data$FishLenComps_LikeType <- comp_fishlen_like_vals
  input_list$data$FishAgeComps_Type <- FishAgeComps_Type_Mat
  input_list$data$FishLenComps_Type <- FishLenComps_Type_Mat

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the fishery age comps
  if("ln_FishAge_theta" %in% names(starting_values)) input_list$par$ln_FishAge_theta <- starting_values$ln_FishAge_theta
  else input_list$par$ln_FishAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery age comps
  if("FishAge_corr_pars" %in% names(starting_values)) input_list$par$FishAge_corr_pars <- starting_values$FishAge_corr_pars
  else input_list$par$FishAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishAge_theta_agg" %in% names(starting_values)) input_list$par$ln_FishAge_theta_agg <- starting_values$ln_FishAge_theta_agg
  else input_list$par$ln_FishAge_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  # aggregated correlation parameters
  if("FishAge_corr_pars_agg" %in% names(starting_values)) input_list$par$FishAge_corr_pars_agg <- starting_values$FishAge_corr_pars_agg
  else input_list$par$FishAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Dispersion parameters for fishery length comps
  if("ln_FishLen_theta" %in% names(starting_values)) input_list$par$ln_FishLen_theta <- starting_values$ln_FishLen_theta
  else input_list$par$ln_FishLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # logistic normal correlation parameters for fishery length comps
  if("FishLen_corr_pars" %in% names(starting_values)) input_list$par$FishLen_corr_pars <- starting_values$FishLen_corr_pars
  else input_list$par$FishLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_fish_fleets, 2))

  # aggregated
  if("ln_FishLen_theta_agg" %in% names(starting_values)) input_list$par$ln_FishLen_theta_agg <- starting_values$ln_FishLen_theta_agg
  else input_list$par$ln_FishLen_theta_agg <- array(0, dim = c(input_list$data$n_fish_fleets))

  if("FishLen_corr_pars_agg" %in% names(starting_values)) input_list$par$FishLen_corr_pars_agg <- starting_values$FishLen_corr_pars_agg
  else input_list$par$FishLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_fish_fleets))

  # Mapping Options ---------------------------------------------------------

  input_list <- do_FishAge_theta_mapping(input_list)
  input_list <- do_FishLen_theta_mapping(input_list)
  input_list <- do_FishAge_corr_pars_mapping(input_list)
  input_list <- do_FishLen_corr_pars_mapping(input_list)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Map fishery selectivity fixed-effect parameters
#'
#' Constructs the factor map for \code{ln_fish_fixed_sel_pars} (e.g., \eqn{a_{50}},
#' \eqn{k}, \eqn{a_{max}}), controlling whether selectivity shape parameters are
#' estimated independently or shared across regions, sexes, or fleets. Cells with
#' no catch data (\code{UseCatch == 0}) are automatically mapped to \code{NA}.
#'
#' Fleet sharing (\code{"est_shared_f_x"}) is handled in a second pass after all
#' base fleet mappings are established, copying the reference fleet's index
#' assignments into the sharing fleet.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param fish_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}.
#'   Each element specifies the estimation structure for one fleet. Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate parameters per region × sex × block.}
#'     \item{\code{"est_shared_r"}}{Parameters shared across regions; unique per sex × block.}
#'     \item{\code{"est_shared_s"}}{Parameters shared across sexes; unique per region × block.}
#'     \item{\code{"est_shared_r_s"}}{Parameters shared across regions and sexes; unique per block.}
#'     \item{\code{"est_shared_f_x"}}{Copy parameters from fleet \code{x} (e.g.,
#'       \code{"est_shared_f_2"} shares with fleet 2). Fleet \code{x} must not
#'       itself use \code{"est_shared_f_y"}.}
#'     \item{\code{"fix"}}{All parameters fixed at starting values (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_fish_fixed_sel_pars}
#'   set to a factor vector.
#'
#' @keywords internal
do_fish_fixed_sel_pars_mapping <- function(input_list, fish_fixed_sel_pars_spec) {

  # Initialize counter and mapping array for fixed effects fishery selectivity
  fish_fixed_sel_pars_counter <- 1
  map_fish_fixed_sel_pars <- input_list$par$ln_fish_fixed_sel_pars
  map_fish_fixed_sel_pars[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate Options
    if(!fish_fixed_sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix", "est_shared_s") &&
       !stringr::str_detect(fish_fixed_sel_pars_spec[f], "est_shared_f_\\d+"))
      stop("fish_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")

    # Skip fleet sharing specs in first pass
    if(stringr::str_detect(fish_fixed_sel_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # Only add a counter if caatches are avaliable in some years for a given region and fleet combination
      if(sum(input_list$data$UseCatch[r,,,f]) > 0) {

        # Extract number of fishery selectivity blocks
        fishsel_blocks_tmp <- unique(as.vector(input_list$data$fish_sel_blocks[r,,f]))

        for(s in 1:input_list$data$n_sexes) {
          for(b in 1:length(fishsel_blocks_tmp)) {

            block_years <- which(input_list$data$fish_sel_blocks[r,,f] == fishsel_blocks_tmp[b]) # figure out block years
            sel_model_this_block <- unique(input_list$data$fish_sel_model[r, block_years, f]) # get selectivity form for a given block
            if(length(sel_model_this_block) > 1) stop("Block ", fishsel_blocks_tmp[b], " for fleet ", f, " region ", r, " has multiple selectivity models assigned to it")

            # determine maximum selectivity parameters
            if(sel_model_this_block == 2) max_sel_pars <- 1 # exponential
            if(sel_model_this_block %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(sel_model_this_block == 4) max_sel_pars <- 6 # double normal

            for(i in 1:max_sel_pars) {

              # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
              if(fish_fixed_sel_pars_spec[f] == "est_all") {
                map_fish_fixed_sel_pars[r,i,b,s,f] <- fish_fixed_sel_pars_counter
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fish_fixed_sel_pars_spec[f] == 'est_shared_r' && r == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  if(fishsel_blocks_tmp[b] %in% input_list$data$fish_sel_blocks[rr,,f]) {
                    map_fish_fixed_sel_pars[rr, i, b, s, f] <- fish_fixed_sel_pars_counter
                  } # end if
                } # end rr loop
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fish_fixed_sel_pars_spec[f] == 'est_shared_s' && s == 1) {
                for(ss in 1:input_list$data$n_sexes) {
                  map_fish_fixed_sel_pars[r, i, b, ss, f] <- fish_fixed_sel_pars_counter
                } # end ss loop
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fish_fixed_sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  for(ss in 1:input_list$data$n_sexes) {
                    if(fishsel_blocks_tmp[b] %in% input_list$data$fish_sel_blocks[rr,,f]) {
                      map_fish_fixed_sel_pars[rr, i, b, ss, f] <- fish_fixed_sel_pars_counter
                    } # end if
                  } # end ss loop
                } #end rr loop
                fish_fixed_sel_pars_counter <- fish_fixed_sel_pars_counter + 1
              } # end if

            } # end i loop
          } # end b loop
        } # end s loop
      } # end if statement
    } # end r loop

    # fix all parameters
    if(fish_fixed_sel_pars_spec[f] == "fix") map_fish_fixed_sel_pars[,,,,f] <- NA
    collect_message("fish_fixed_sel_pars_spec is specified as: ", fish_fixed_sel_pars_spec[f], " for fishery fleet ", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(fish_fixed_sel_pars_spec[f], "est_shared_f")) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(fish_fixed_sel_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(fish_fixed_sel_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_fish_fixed_sel_pars[,,,,f] <- map_fish_fixed_sel_pars[,,,,flt_shared]
      collect_message("fish_fixed_sel_pars_spec is specified as: ", fish_fixed_sel_pars_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_fish_fixed_sel_pars <- factor(map_fish_fixed_sel_pars)
  return(input_list)
}


#' Map fishery catchability parameters
#'
#' Constructs the factor map for \code{ln_fish_q}, controlling whether
#' catchability parameters are estimated independently per region and time block
#' or shared across regions. Cells with no fishery index observations
#' (\code{UseFishIdx == 0}) are automatically mapped to \code{NA}.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param fish_q_spec Character vector of length \code{n_fish_fleets}. Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate catchability per region × block × fleet.}
#'     \item{\code{"est_shared_r"}}{Single catchability shared across regions,
#'       unique per block × fleet.}
#'     \item{\code{"fix"}}{All catchability parameters fixed (mapped to \code{NA}).}
#'   }
#'
#' @return The input \code{input_list} with \code{$map$ln_fish_q} set to a
#'   factor vector.
#'
#' @keywords internal
do_fish_q_mapping <- function(input_list, fish_q_spec) {

  # Initialize counter and mapping array for fishery catchability
  fish_q_counter <- 1
  map_fish_q <- input_list$par$ln_fish_q
  map_fish_q[] <- NA

  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate options
    if(!is.null(fish_q_spec)) {
      if(!fish_q_spec[f] %in% c("est_all", "est_shared_r", "fix"))
        stop("fish_q_spec not correctly specfied. Should be one of these: est_all, est_shared_r, fix")
    }

    for(r in 1:input_list$data$n_regions) {

      if(sum(input_list$data$UseFishIdx[r,,,f]) == 0) {
        map_fish_q[r,,f] <- NA # fix parameters if we are not using fishery indices for these fleets and regions
      } else {

        # Extract number of fishery catchability blocks
        fishq_blocks_tmp <- unique(as.vector(input_list$data$fish_q_blocks[r,,f]))

        for(b in 1:length(fishq_blocks_tmp)) {

          # Estimate for all regions
          if(fish_q_spec[f] == 'est_all') {
            map_fish_q[r,b,f] <- fish_q_counter
            fish_q_counter <- fish_q_counter + 1
          } # end if

          # Estimate but share q across regions
          if(fish_q_spec[f] == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              if(fishq_blocks_tmp[b] %in% input_list$data$fish_q_blocks[rr,,f]) {
                map_fish_q[rr, b, f] <- fish_q_counter
              } # end if
            } # end rr loop
            fish_q_counter <- fish_q_counter + 1
          } # end if

        } # end b loop
      } # end else loop
    } # end r loop

    # fix all parameters
    if(fish_q_spec[f] == 'fix') map_fish_q[,,f] <- NA
    collect_message("fish_q_spec is specified as: ", fish_q_spec[f], " for fishery fleet ", f)
  } # end f loop

  # input into mapping list
  input_list$map$ln_fish_q <- factor(map_fish_q)

  return(input_list)
}

#' Map fishery selectivity process error hyperparameters
#'
#' Constructs the factor map for \code{fishsel_pe_pars}, which contains the
#' variance and correlation hyperparameters governing continuous time-varying
#' selectivity. The set of active parameters depends on the time-variation type
#' (\code{cont_tv_fish_sel}): iid/random-walk forms use up to 2 parameters
#' (log-sigma); 3D GMRF forms use up to 4 (partial correlations for age, year,
#' cohort dimensions plus log-sigma); the 2D AR1 form uses 3 (bin AR1, year AR1,
#' log-sigma). Correlation components can be selectively suppressed via
#' \code{corr_opt_semipar}.
#'
#' Fleet sharing (\code{"est_shared_f_x"}) is handled in a second pass.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param fishsel_pe_pars_spec Character vector of length \code{n_fish_fleets}.
#'   Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate hyperparameters per region × sex.}
#'     \item{\code{"est_shared_r"}}{Shared across regions; unique per sex.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes; unique per region.}
#'     \item{\code{"est_shared_r_s"}}{Shared across regions and sexes.}
#'     \item{\code{"est_shared_f_x"}}{Copy hyperparameters from fleet \code{x}.}
#'     \item{\code{"fix"} or \code{"none"}}{All parameters fixed (mapped to \code{NA}).}
#'   }
#' @param corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   specifying which correlation components to suppress for semi-parametric
#'   models. Valid values per fleet: \code{NA} (no suppression),
#'   \code{"corr_zero_y"}, \code{"corr_zero_b"}, \code{"corr_zero_y_b"},
#'   \code{"corr_zero_c"}, \code{"corr_zero_y_c"}, \code{"corr_zero_b_c"},
#'   \code{"corr_zero_y_b_c"}. Cohort options (\code{"corr_zero_c"}, etc.) are
#'   only valid for 3D GMRF forms and will error if applied to the 2D AR1
#'   (\code{cont_tv_fish_sel == 5}).
#'
#' @return The input \code{input_list} with \code{$map$fishsel_pe_pars} set to a
#'   factor vector. Index numbering is reset after any correlation suppression to
#'   maintain contiguous integer indices.
#'
#' @keywords internal
do_fishsel_pe_pars_mapping <- function(input_list, fishsel_pe_pars_spec, corr_opt_semipar) {

  # Initialize counter and mapping array for fishery process errors
  fishsel_pe_pars_counter <- 1 # initalize counter
  map_fishsel_pe_pars <- input_list$par$fishsel_pe_pars # initalize array
  map_fishsel_pe_pars[] <- NA

  # Fishery process error parameters
  for(f in 1:input_list$data$n_fish_fleets) {

    # Validate options
    if(!is.null(fishsel_pe_pars_spec)) {
      if(!fishsel_pe_pars_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s") &&
         !stringr::str_detect(fishsel_pe_pars_spec[f], "est_shared_f_\\d+"))
        stop("fishsel_pe_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    }

    # Skip fleet sharing specs in first pass
    if(!is.null(fishsel_pe_pars_spec)) if(stringr::str_detect(fishsel_pe_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data$cont_tv_fish_sel[r,f] == 0 || sum(input_list$data$UseCatch[r,,,f]) == 0) {
        map_fishsel_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation

        # Figure out max number of selectivity parameters for a given region and fleet
        if(unique(input_list$data$fish_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
        if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
        if(unique(input_list$data$fish_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal

        for(s in 1:input_list$data$n_sexes) {

          # If iid time-variation or random walk for this fleet
          if(input_list$data$cont_tv_fish_sel[r,f] %in% c(1,2)) {

            for(i in 1:max_sel_pars) {

              # either fixing parameters or not used for a given fleet
              if(fishsel_pe_pars_spec[f] %in% c("none", "fix")) map_fishsel_pe_pars[r,i,s,f] <- NA

              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == "est_all") {
                map_fishsel_pe_pars[r,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_fishsel_pe_pars[,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_fishsel_pe_pars[r,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_fishsel_pe_pars[,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

            } # end i loop
          } # end iid or random walk variation

          # If 3d gmrf or 2dar1
          if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_bin, 2 = pcorr_year, 4 = log_sigma)
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4) && input_list$data$Selex_Type == 1) stop("Cohort-based selectivity deviations are specified, but selectivity is specified as length-based. Please choose another deviation form!")

            for(i in idx) {

              # either fixing parameters or not used for a given fleet
              if(fishsel_pe_pars_spec[f] %in% c("none", "fix")) map_fishsel_pe_pars[r,i,s,f] <- NA

              # Estimating all process error parameters
              if(fishsel_pe_pars_spec[f] == "est_all") {
                map_fishsel_pe_pars[r,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_fishsel_pe_pars[,i,s,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_fishsel_pe_pars[r,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(fishsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_fishsel_pe_pars[,i,,f] <- fishsel_pe_pars_counter
                fishsel_pe_pars_counter <- fishsel_pe_pars_counter + 1
              }

            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(corr_opt_semipar)) {

              opt <- input_list$data$cont_tv_fish_sel[r,f] # get random effects options

              # Validate options
              if(!corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_b", "corr_zero_y_b", "corr_zero_c", "corr_zero_y_c", "corr_zero_b_c", "corr_zero_y_b_c"))
                stop("corr_opt_semipar not correctly specfied. Should be one of these: corr_zero_y, corr_zero_b, corr_zero_y_b, corr_zero_c, corr_zero_y_c, corr_zero_b_c, corr_zero_y_b_c, NA")
              if(opt == 5 && corr_opt_semipar[f] %in% c("corr_zero_c","corr_zero_y_c","corr_zero_b_c","corr_zero_y_b_c"))
                stop("Invalid corr_opt_semipar for 2dar1 (opt=5): cohort correlations are not allowed.")

              if (opt %in% c(3,4,5)) {
                # 2d and 3d options
                if (corr_opt_semipar[f] == "corr_zero_y")    map_fishsel_pe_pars[,2,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_b")    map_fishsel_pe_pars[,1,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b")  map_fishsel_pe_pars[,1:2,,f]   <- NA
              }

              if(opt %in% c(3,4)) {
                # 3d gmrf options only (adds the cohort dimension)
                if (corr_opt_semipar[f] == "corr_zero_c")      map_fishsel_pe_pars[,3,,f]   <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_c")    map_fishsel_pe_pars[,2:3,,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_b_c")    map_fishsel_pe_pars[,c(1,3),,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b_c")  map_fishsel_pe_pars[,1:3,,f] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_fishsel_pe_pars))
              map_fishsel_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("corr_opt_semipar is specified as: ", corr_opt_semipar[f], "for fishery fleet", f)

            }
          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(fishsel_pe_pars_spec[f] == "fix") map_fishsel_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(fishsel_pe_pars_spec)) collect_message("fishsel_pe_pars_spec is specified as: ", fishsel_pe_pars_spec[f], "for fishery fleet", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(fishsel_pe_pars_spec[f], "est_shared_f") && !is.null(fishsel_pe_pars_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(fishsel_pe_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(fishsel_pe_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_fishsel_pe_pars[,,,f] <- map_fishsel_pe_pars[,,,flt_shared]
      collect_message("fishsel_pe_pars_spec is specified as: ", fishsel_pe_pars_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$fishsel_pe_pars <- factor(map_fishsel_pe_pars)

  return(input_list)
}

#' Map fishery selectivity deviation parameters
#'
#' Constructs the factor map for \code{ln_fishsel_devs}, the annual deviations in
#' continuous time-varying fishery selectivity. For iid and random-walk forms,
#' the deviation dimension corresponds to selectivity parameters (up to 6 for
#' double-normal); for semi-parametric forms (3D GMRF, 2D AR1), it corresponds
#' to age or length bins. Cells with no time-variation
#' (\code{cont_tv_fish_sel == 0}) or no catch data are mapped to \code{NA}.
#'
#' Age-sharing (\code{"est_shared_a"} and related options) groups bins into
#' blocks defined by \code{fishsel_devs_shared_ages}, reducing the number of
#' estimated deviation series. Fleet sharing (\code{"est_shared_f_x"}) is handled
#' in a second pass. The resulting integer map is also stored as
#' \code{$data$map_ln_fishsel_devs} for use in the C++ objective function.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, and \code{$map}
#'   sublists.
#' @param fish_sel_devs_spec Character vector of length \code{n_fish_fleets}.
#'   Options:
#'   \describe{
#'     \item{\code{"est_all"}}{Separate deviation series per region × sex × bin.}
#'     \item{\code{"est_shared_r"}}{Shared across regions.}
#'     \item{\code{"est_shared_s"}}{Shared across sexes.}
#'     \item{\code{"est_shared_r_s"}}{Shared across regions and sexes.}
#'     \item{\code{"est_shared_a"}}{Shared across age/bin groups defined by
#'       \code{fishsel_devs_shared_ages}.}
#'     \item{\code{"est_shared_r_a"}}{Shared across regions and age groups.}
#'     \item{\code{"est_shared_a_s"}}{Shared across age groups and sexes.}
#'     \item{\code{"est_shared_r_a_s"}}{Shared across regions, age groups, and sexes.}
#'     \item{\code{"est_shared_f_x"}}{Copy deviation map from fleet \code{x}.}
#'     \item{\code{"fix"} or \code{"none"}}{All deviations fixed at zero (mapped to \code{NA}).}
#'   }
#'   Age-sharing options (\code{"est_shared_a"}, etc.) are only valid with
#'   semi-parametric time-varying forms (\code{cont_tv_fish_sel %in% c(3,4,5)}).
#' @param fishsel_devs_shared_ages List of integer vectors defining age/bin
#'   groupings for age-sharing options. Each element groups bins that share a
#'   single deviation series, e.g., \code{list(1:5, 6:10, 11:30)}.
#'   Only used when \code{fish_sel_devs_spec} includes \code{"est_shared_a"}.
#'
#' @return The input \code{input_list} with \code{$map$ln_fishsel_devs} set to a
#'   factor vector and \code{$data$map_ln_fishsel_devs} set to the corresponding
#'   integer array (for use in the C++ template).
#'
#' @keywords internal
do_fishsel_devs_mapping <- function(input_list, fish_sel_devs_spec, fishsel_devs_shared_ages) {

  # Initialize counter and mapping array for fishery selectivity deviations
  fishsel_devs_counter <- 1
  map_fishsel_devs <- input_list$par$ln_fishsel_devs
  map_fishsel_devs[] <- NA

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:input_list$data$n_fish_fleets) {

      # Validate options
      if(!is.null(fish_sel_devs_spec)) {
        if(!fish_sel_devs_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "est_shared_a", "est_shared_r_a", "est_shared_r_a_s", "est_shared_a_s") &&
           !stringr::str_detect(fish_sel_devs_spec[f], "est_shared_f_\\d+"))
          stop("fish_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, est_shared_a, est_shared_r_a, est_shared_r_a_s, est_shared_r_s, fix, or est_shared_f_# (where # is fleet number)")
        if(fish_sel_devs_spec[f] %in% c("est_shared_a", "est_shared_r_a", "est_shared_r_a_s", "est_shared_a_s") &&
           !input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) stop("Sharing age deviations with iid or random walk parametric forms is not supported!")
       }

      # Skip fleet sharing specs in first pass
      if(!is.null(fish_sel_devs_spec)) if(stringr::str_detect(fish_sel_devs_spec[f], "est_shared_f")) next

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:(length(input_list$data$years) + input_list$data$n_proj_yrs_devs)) {

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data$cont_tv_fish_sel[r,f] == 0 || sum(input_list$data$UseCatch[r,,,f]) == 0) {
            map_fishsel_devs[r,y,,s,f] <- NA
          } else {

            # Figure out max number of selectivity parameters for a given region and fleet
            if(unique(input_list$data$fish_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
            if(unique(input_list$data$fish_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(unique(input_list$data$fish_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal

            # If iid or random walk time-variation for this fleet
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(1,2)) {

              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(fish_sel_devs_spec[f] == 'est_all') {
                  map_fishsel_devs[r,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(fish_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_fishsel_devs[,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(fish_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_fishsel_devs[r,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(fish_sel_devs_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                  map_fishsel_devs[,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data$cont_tv_fish_sel[r,f] %in% c(3,4,5)) {

              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and bins
                if(fish_sel_devs_spec[f] == 'est_all') {
                  map_fishsel_devs[r,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across regions
                if(fish_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_fishsel_devs[,y,i,s,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes
                if(fish_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_fishsel_devs[r,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes and regions
                if(fish_sel_devs_spec[f] == 'est_shared_r_s' && s == 1 && r == 1) {
                  map_fishsel_devs[,y,i,,f] <- fishsel_devs_counter
                  fishsel_devs_counter <- fishsel_devs_counter + 1
                }

                if(fish_sel_devs_spec[f] == 'est_shared_a') {
                  for(k in 1:length(fishsel_devs_shared_ages)) {
                    map_fishsel_devs[r,y,fishsel_devs_shared_ages[[k]],s,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

                if(fish_sel_devs_spec[f] == 'est_shared_r_a' && r == 1) {
                  for(k in 1:length(fishsel_devs_shared_ages)) {
                    map_fishsel_devs[,y,fishsel_devs_shared_ages[[k]],s,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

                if(fish_sel_devs_spec[f] == 'est_shared_a_s' && s == 1) {
                  for(k in 1:length(fishsel_devs_shared_ages)) {
                    map_fishsel_devs[r,y,fishsel_devs_shared_ages[[k]],,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

                if(fish_sel_devs_spec[f] == 'est_shared_r_a_s' && s == 1 && r == 1) {
                  for(k in 1:length(fishsel_devs_shared_ages)) {
                    map_fishsel_devs[,y,fishsel_devs_shared_ages[[k]],,f] <- fishsel_devs_counter
                    fishsel_devs_counter <- fishsel_devs_counter + 1
                  } # end k loop
                }

              } # end i loop
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

      if(!is.null(fish_sel_devs_spec)) collect_message("fish_sel_devs_spec is specified as: ", fish_sel_devs_spec[f], "for fishery fleet", f, "and region ", r)

    } # end f loop
  } # end r loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_fish_fleets) {
    if(stringr::str_detect(fish_sel_devs_spec[f], "est_shared_f") && !is.null(fish_sel_devs_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(fish_sel_devs_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_fish_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(fish_sel_devs_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_fishsel_devs[,,,,f] <- map_fishsel_devs[,,,,flt_shared]
      collect_message("fish_sel_devs_spec is specified as: ", fish_sel_devs_spec[f], " for fishery fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_fishsel_devs <- factor(map_fishsel_devs)
  input_list$data$map_ln_fishsel_devs <- array(as.numeric(input_list$map$ln_fishsel_devs), dim = dim(input_list$par$ln_fishsel_devs))

  return(input_list)
}


#' Set up fishery selectivity and catchability specifications
#'
#' Configures all aspects of fishery selectivity and catchability for the
#' estimation model: functional forms, time blocks, continuous time-varying
#' structures, process error hyperparameters, annual deviations, and
#' catchability blocks and estimation structure. Must be called after
#' \code{\link{Setup_Mod_FishIdx_and_Comps}}.
#'
#' Selectivity time-variation and blocked selectivity are mutually exclusive
#' within a fleet — specifying both for the same fleet will raise an error.
#'
#' @param input_list Named list with \code{$data}, \code{$par}, \code{$map},
#'   and \code{$verbose} sublists.
#' @param fish_sel_model Character vector specifying the selectivity functional
#'   form for each fleet (and optionally each time block). Each element must
#'   follow one of:
#'   \itemize{
#'     \item \code{"<model>_Fleet_<f>"} — single form for all years of fleet \code{f}.
#'     \item \code{"<model>_Fleet_<f>_Block_<b>"} — form specific to block \code{b}
#'       of fleet \code{f}, as defined in \code{fish_sel_blocks}.
#'   }
#'   Available models:
#'   \describe{
#'     \item{\code{"logist1"}}{Logistic with \eqn{a_{50}} and slope \eqn{k} (2 parameters).}
#'     \item{\code{"logist2"}}{Logistic with \eqn{a_{50}} and \eqn{a_{95}} (2 parameters).}
#'     \item{\code{"gamma"}}{Dome-shaped gamma with \eqn{a_{max}} and \eqn{\delta} (2 parameters).}
#'     \item{\code{"exponential"}}{Exponential with a single power parameter (1 parameter).}
#'     \item{\code{"dbnrml"}}{Double-normal with 6 parameters.}
#'   }
#'   See the model equations vignette for mathematical definitions.
#' @param cont_tv_fish_sel Character vector of length \code{n_fish_fleets}
#'   specifying continuous time-varying selectivity per fleet. Each element
#'   must be \code{"<type>_Fleet_<f>"}. Valid types:
#'   \describe{
#'     \item{\code{"none"}}{No continuous time-variation (default).}
#'     \item{\code{"iid"}}{IID annual deviations on selectivity parameters.}
#'     \item{\code{"rw"}}{Random walk in selectivity parameters over time.}
#'     \item{\code{"3dmarg"}}{3D GMRF with marginal variance parameterisation.}
#'     \item{\code{"3dcond"}}{3D GMRF with conditional variance parameterisation.}
#'     \item{\code{"2dar1"}}{2D separable AR1 in bin and year dimensions.}
#'   }
#'   If any fleet has \code{cont_tv_fish_sel != "none"}, both
#'   \code{fishsel_pe_pars_spec} and \code{fish_sel_devs_spec} must also be
#'   provided.
#' @param fish_sel_blocks Character vector defining discrete selectivity time
#'   blocks per fleet. Each element follows \code{"Block_<b>_Year_<s>-<e>_Fleet_<f>"}
#'   or \code{"Block_<b>_Year_<s>-terminal_Fleet_<f>"}. Use
#'   \code{"none_Fleet_<f>"} (default) for a single constant block. Blocks must
#'   be non-overlapping and together span all model years for the specified fleet.
#'   Mutually exclusive with \code{cont_tv_fish_sel != "none"} for the same fleet.
#' @param fish_q_blocks Character vector defining catchability time blocks per
#'   fleet, using the same format as \code{fish_sel_blocks}. Default
#'   \code{"none_Fleet_<f>"} gives a single constant block.
#' @param fish_fixed_sel_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying how fixed-effect selectivity parameters are estimated. See
#'   \code{\link{do_fish_fixed_sel_pars_mapping}} for all options
#'   (\code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"},
#'   \code{"est_shared_r_s"}, \code{"est_shared_f_x"}, \code{"fix"}).
#' @param fish_q_spec Character vector of length \code{n_fish_fleets} specifying
#'   catchability estimation structure. See \code{\link{do_fish_q_mapping}} for
#'   options (\code{"est_all"}, \code{"est_shared_r"}, \code{"fix"}).
#' @param fishsel_pe_pars_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for selectivity process error
#'   hyperparameters. Required when any fleet has continuous time-variation.
#'   See \code{\link{do_fishsel_pe_pars_mapping}} for all options.
#' @param fish_sel_devs_spec Character vector of length \code{n_fish_fleets}
#'   specifying the estimation structure for annual selectivity deviations.
#'   Required when any fleet has continuous time-variation. See
#'   \code{\link{do_fishsel_devs_mapping}} for all options including age-sharing
#'   options for semi-parametric forms.
#' @param fishsel_devs_shared_ages List of integer vectors grouping age or length
#'   bins that share a single deviation series. Only used when
#'   \code{fish_sel_devs_spec} contains one of the \code{"est_shared_a"} variants.
#'   Example: \code{list(1:5, 6:10, 11:30)}.
#' @param corr_opt_semipar Character vector of length \code{n_fish_fleets}
#'   controlling which correlation components to suppress in semi-parametric
#'   (3D GMRF or 2D AR1) time-varying selectivity. Set to \code{NA} (default)
#'   for no suppression. See \code{\link{do_fishsel_pe_pars_mapping}} for valid
#'   suppression codes. Cohort-correlation options are invalid for \code{"2dar1"}.
#' @param cont_tv_fish_sel_penalty Logical. If \code{TRUE} (default), applies a
#'   penalty on the continuous time-varying selectivity deviations to regularise
#'   the process.
#' @param Use_fish_q_prior Integer flag. \code{1} = apply lognormal priors to
#'   catchability; \code{0} = no priors (default). Requires \code{fish_q_prior}.
#' @param fish_q_prior Data frame of catchability prior hyperparameters. Required
#'   columns: \code{region}, \code{fleet}, \code{block} (block index), \code{mu}
#'   (prior mean on natural scale), \code{sd} (prior SD on log scale). Each row
#'   specifies a \eqn{\text{Normal}(\log(\mu), \sigma)} prior for one catchability
#'   parameter. Only used when \code{Use_fish_q_prior = 1}.
#' @param Use_fish_selex_prior Integer flag. \code{1} = apply lognormal priors to
#'   selectivity parameters; \code{0} = no priors (default). Requires
#'   \code{fish_selex_prior}.
#' @param fish_selex_prior Data frame of selectivity prior hyperparameters.
#'   Required columns: \code{region}, \code{fleet}, \code{block}, \code{sex},
#'   \code{par} (parameter index within the functional form), \code{mu}, \code{sd}.
#'   Only used when \code{Use_fish_selex_prior = 1}.
#' @param ... Optional starting value overrides, passed by name. Recognised
#'   arguments:
#'   \describe{
#'     \item{\code{ln_fish_fixed_sel_pars}}{Array dimensioned
#'       \code{[n_regions × max_pars × max_blocks × n_sexes × n_fish_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{ln_fish_q}}{Array dimensioned
#'       \code{[n_regions × max_q_blocks × n_fish_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{fishsel_pe_pars}}{Array dimensioned
#'       \code{[n_regions × 4 × n_sexes × n_fish_fleets]}.
#'       Default: \code{0}.}
#'     \item{\code{ln_fishsel_devs}}{Array dimensioned
#'       \code{[n_regions × (n_years + n_proj_yrs_devs) × n_bins × n_sexes × n_fish_fleets]}.
#'       Default: \code{0}.}
#'   }
#'
#' @return The input \code{input_list} with \code{$data}, \code{$par}, and
#'   \code{$map} updated. Key additions include the parsed integer arrays for
#'   \code{cont_tv_fish_sel}, \code{fish_sel_blocks}, \code{fish_sel_model}, and
#'   \code{fish_q_blocks}; starting value arrays for all four parameter groups;
#'   and their corresponding factor maps.
#'
#'
#' @export Setup_Mod_Fishsel_and_Q
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Fishsel_and_Q <- function(input_list,
                                    cont_tv_fish_sel = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fish_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fish_sel_model,
                                    Use_fish_q_prior = 0,
                                    fish_q_prior = NA,
                                    fish_q_blocks = paste("none_Fleet_", 1:input_list$data$n_fish_fleets, sep = ''),
                                    fishsel_pe_pars_spec = NULL,
                                    fish_fixed_sel_pars_spec = NULL,
                                    fish_q_spec = NULL,
                                    fish_sel_devs_spec = NULL,
                                    corr_opt_semipar = NULL,
                                    Use_fish_selex_prior = 0,
                                    fish_selex_prior = NULL,
                                    cont_tv_fish_sel_penalty = TRUE,
                                    fishsel_devs_shared_ages = NULL,
                                    ...
                                    ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Selectivity Type
  if(is.null(input_list$data$Selex_Type)) stop("Selectivity type (age or length-based) has not been specified yet! Make sure to first specify biological inputs with Setup_Mod_Biologicals.")

  # Continuous Selectivity Deviations
  if(!is.null(fishsel_pe_pars_spec)) if(length(fishsel_pe_pars_spec) != input_list$data$n_fish_fleets) stop("fishsel_pe_pars_spec is not length n_fish_fleets")
  if(!is.null(fish_sel_devs_spec)) if(length(fish_sel_devs_spec) != input_list$data$n_fish_fleets) stop("fish_sel_devs_spec is not length n_fish_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_fish_fleets) stop("corr_opt_semipar is not length n_fish_fleets")

  # Catchability Priors
  if(!Use_fish_q_prior %in% c(0,1)) stop("Values for Use_fish_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking catchability priors
  if(Use_fish_q_prior == 1) {
    required_cols <- c("region", "fleet", "block", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(fish_q_prior))
    if(length(missing_cols) > 0) {
      stop("fish_q_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Fishery Catchability priors are: ", ifelse(Use_fish_q_prior == 0, "Not Used", "Used"))

  # Selectivity Priors
  if(!Use_fish_selex_prior %in% c(0,1)) stop("Values for Use_fish_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_fish_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(fish_selex_prior))
    if(length(missing_cols) > 0) {
      stop("fish_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Fishery Selectivity priors are: ", ifelse(Use_fish_selex_prior == 0, "Not Used", "Used"))

  # Continuous Time-Varying Selectivity Options -----------------------------
  cont_tv_fish_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_fish_sel)) {
    # Extract out components from list
    tmp <- cont_tv_fish_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for cont_tv_fish_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_fish_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_fish_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous fishery time-varying selectivity specified as: ", cont_tv_type, " for fishery fleet ", fleet)
  }

  if(any(cont_tv_fish_sel_mat > 0) && is.null(fishsel_pe_pars_spec) && is.null(fish_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but fishsel_pe_pars_spec and/or fish_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  fish_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_sel_blocks)) {

    # Extract out components from list
    tmp <- fish_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_sel_blocks_arr[,,fleet] <- 1 # input only 1 fishery time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # extract fleet index

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fish_sel_blocks_arr[,years,fleet] <- block_val
    }

  }

  if(any(is.na(fish_sel_blocks_arr))) stop("Fishery Selectivtiy Blocks are returning an NA. Did you forget to specify the year range of fish_sel_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Selectivity Time Blocks for fishery", f, "is specified at:", length(unique(fish_sel_blocks_arr[,,f]))))

  # Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml"), num = c(0,1,2,3,4)) # set up values we can map to
  fish_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))

  for(i in 1:length(fish_sel_model)) {

    # Extract out fishery selectivity components from vector
    tmp_sel_form <- fish_sel_model[i]
    tmp_sel_form_vec <- unlist(strsplit(tmp_sel_form, "_")) # split string
    sel_form <- tmp_sel_form_vec[1] # get selectivity type

    # get fleet index
    tmp_fleet <- if(length(tmp_sel_form_vec) == 3) as.numeric(tmp_sel_form_vec[3]) else as.numeric(tmp_sel_form_vec[5]) # fleet index changes if block is included in character vector
    # get block index
    tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[3]) else NULL

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("fish_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_fish_fleets)) stop("Invalid fleet specified for fish_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) fish_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
    else fish_sel_model_arr[,which(fish_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)]
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Fishery selectivity functional form specified as:", sel_form, " for fishery fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_fish_fleets) {
    has_blocks <- length(unique(fish_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_fish_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Blocked Catchability Options --------------------------------------------
  fish_q_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
  for(i in 1:length(fish_q_blocks)) {
    # Extract out components from list
    tmp <- fish_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Validate options
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Fishery Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      fish_q_blocks_arr[,,fleet] <- 1 # input only 1 fishery catchability time block
    }

    if(tmp_vec[1] == "Block") {

      block_val <- as.numeric(tmp_vec[2]) # get block value
      fleet <- as.numeric(tmp_vec[6]) # get fleet number

      # get year ranges
      if(!str_detect(tmp, "terminal")) { # if not terminal year
        year_range <- as.numeric(unlist(strsplit(tmp_vec[4], "-")))
        years <- year_range[1]:year_range[2] # get sequence of years
      } else { # if terminal year
        year_range <- unlist(strsplit(tmp_vec[4], '-'))[1] # get year range
        years <- as.numeric(year_range):length(input_list$data$years) # get sequence of years
      }

      fish_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }

  if(any(is.na(fish_q_blocks))) stop("Fishery Catchability Blocks are returning an NA. Did you forget to specify the year range of fish_q_blocks?")
  for(f in 1:input_list$data$n_fish_fleets) collect_message(paste("Fishery Catchability Time Blocks for fishery", f, "is specified at:", length(unique(fish_q_blocks_arr[,,f]))))

  # Populate Data List ------------------------------------------------------

  input_list$data$cont_tv_fish_sel <- cont_tv_fish_sel_mat
  input_list$data$cont_tv_fish_sel_penalty <- cont_tv_fish_sel_penalty
  input_list$data$fish_sel_blocks <- fish_sel_blocks_arr
  input_list$data$fish_sel_model <- fish_sel_model_arr
  input_list$data$fish_q_blocks <- fish_q_blocks_arr
  input_list$data$fish_q_prior <- fish_q_prior
  input_list$data$Use_fish_q_prior <- Use_fish_q_prior
  input_list$data$Use_fish_selex_prior <- Use_fish_selex_prior
  input_list$data$fish_selex_prior <- fish_selex_prior

  # Populate Parameter List -------------------------------------------------

  # Figure out number of selectivity parameters for a given functional form
  unique_fishsel_vals <- unique(as.vector(input_list$data$fish_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_fishsel_vals)) {
    if(unique_fishsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_fishsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_fishsel_vals[i] %in% c(4)) sel_pars_vec[i] <- 6 # double normal
  } # end i loop

  # figure out maximum number of fishery selectivity blocks for a given reigon and fleet
  max_fishsel_blks <- max(apply(input_list$data$fish_sel_blocks, c(1,3), FUN = function(x) length(unique(x))))
  # maximum number of selectivity parameters across all forms
  max_fishsel_pars <- max(sel_pars_vec)
  if("ln_fish_fixed_sel_pars" %in% names(starting_values)) input_list$par$ln_fish_fixed_sel_pars <- starting_values$ln_fish_fixed_sel_pars
  else input_list$par$ln_fish_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_fishsel_pars, max_fishsel_blks, input_list$data$n_sexes, input_list$data$n_fish_fleets))

  # Fishery catchability
  max_fishq_blks <- max(apply(input_list$data$fish_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of fishery catchability blocks for a given reigon and fleet
  if("ln_fish_q" %in% names(starting_values)) input_list$par$ln_fish_q <- starting_values$ln_fish_q
  else input_list$par$ln_fish_q <- array(0, dim = c(input_list$data$n_regions, max_fishq_blks, input_list$data$n_fish_fleets))

  # Fishery selectivity process error parameters
  if("fishsel_pe_pars" %in% names(starting_values)) input_list$par$fishsel_pe_pars <- starting_values$fishsel_pe_pars
  else input_list$par$fishsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_fishsel_pars, 4), input_list$data$n_sexes, input_list$data$n_fish_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Fishery selectivity deviations
  if(input_list$data$Selex_Type == 0) bins <- length(input_list$data$ages) # age based deviations
  if(input_list$data$Selex_Type == 1) bins <- length(input_list$data$lens) # length based deviations
  if("ln_fishsel_devs" %in% names(starting_values)) input_list$par$ln_fishsel_devs <- starting_values$ln_fishsel_devs
  else input_list$par$ln_fishsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_fish_fleets))


  # Mapping Options ---------------------------------------------------------
  input_list <- do_fish_fixed_sel_pars_mapping(input_list, fish_fixed_sel_pars_spec)
  input_list <- do_fish_q_mapping(input_list, fish_q_spec)
  input_list <- do_fishsel_pe_pars_mapping(input_list, fishsel_pe_pars_spec, corr_opt_semipar)
  input_list <- do_fishsel_devs_mapping(input_list, fish_sel_devs_spec, fishsel_devs_shared_ages)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}


