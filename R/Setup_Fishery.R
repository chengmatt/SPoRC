#' Setup values and dimensions of fishing processes
#'
#' @param sim_list Simulation list object from `Setup_Sim_Dim()`
#' @param ln_sigmaC Observation error for catch
#'   [n_regions × n_yrs × n_fish_fleets]
#'   (default: `log(0.02)`)
#' @param init_F_val Initial fishing mortality value
#'   (default: `0`)
#' @param Fmort_input Fishing mortality input array
#'   [n_regions × n_yrs × n_fish_fleets × n_sims]
#'   (default: `0.1`)
#' @param fish_sel_input Fishery selectivity array
#'   [n_regions × n_yrs × n_ages × n_sexes × n_fish_fleets × n_sims]
#'   (no default, must be provided)
#' @param fish_q_input Fishery catchability array
#'   [n_regions × n_yrs × n_fish_fleets × n_sims]
#'   (default: `1`)
#'
#' @param ObsFishIdx_SE Observation error of fishery index
#'   [n_regions × n_yrs × n_fish_fleets]
#'   (default: `0.2`)
#' @param fish_idx_type Array of index types [n_regions x n_fish_fleets]
#'   (default: all `1` = biomass index)
#'   \itemize{
#'     \item \code{0}: Abundance index
#'     \item \code{1}: Biomass index
#'   }
#'
#' @param comp_fishage_like Vector [n_fish_fleets] specifying likelihood for simulating age comps
#'   (default: all `0` = multinomial)
#'   \itemize{
#'     \item \code{0}: Multinomial
#'     \item \code{1}: Dirichlet-Multinomial
#'     \item \code{2}: Logistic Normal iid
#'     \item \code{3}: Logistic Normal 1dar1
#'     \item \code{4}: Logistic Normal 2d correlation (constant by sex, 1dar1 by age)
#'   }
#' @param ISS_FishAgeComps Input sample sizes
#'   [n_regions × n_yrs × n_sexes × n_fish_fleets × n_sims]
#'   (default: `100`)
#' @param ln_FishAge_theta Overdispersion parameters
#'   [n_regions × n_sexes × n_fish_fleets]
#'   (default: `log(1)`)
#' @param ln_FishAge_theta_agg Overdispersion parameters for aggregated comps
#'   [n_fish_fleets]
#'   (default: `log(1)`)
#' @param FishAge_corr_pars_agg Correlation parameters (agg.) for options 3–4
#'   [n_fish_fleets]
#'   (default: `0.01`)
#' @param FishAge_corr_pars Correlation parameters
#'   [n_regions × n_sexes × n_fish_fleets x 2]
#'   (default: `0.01`)
#' @param FishAgeComps_Type Array [n_yrs × n_fish_fleets]
#'   (default: `2` = joint by sex, split by region)
#'   \itemize{
#'     \item \code{0}: Aggregated
#'     \item \code{1}: Split by sex and region
#'     \item \code{2}: Joint by sex, split by region
#'     \item \code{999}: Not simulated
#'   }
#'
#' @param comp_fishlen_like Vector [n_fish_fleets] specifying likelihood for simulating length comps
#'   (default: all `0` = multinomial)
#'   \itemize{
#'     \item \code{0}: Multinomial
#'     \item \code{1}: Dirichlet-Multinomial
#'     \item \code{2}: Logistic Normal iid
#'     \item \code{3}: Logistic Normal 1dar1
#'     \item \code{4}: Logistic Normal 2d correlation (constant by sex, 1dar1 by length)
#'   }
#' @param ISS_FishLenComps Input sample sizes
#'   [n_regions × n_yrs × n_sexes × n_fish_fleets × n_sims]
#'   (default: `100`)
#' @param ln_FishLen_theta Overdispersion parameters
#'   [n_regions × n_sexes × n_fish_fleets x 2]
#'   (default: `log(1)`)
#' @param ln_FishLen_theta_agg Overdispersion parameters for aggregated comps
#'   [n_fish_fleets]
#'   (default: `log(1)`)
#' @param FishLen_corr_pars_agg Correlation parameters (agg.) for options 3–4
#'   [n_fish_fleets]
#'   (default: `0.01`)
#' @param FishLen_corr_pars Correlation parameters
#'   [n_regions × n_sexes × n_fish_fleets]
#'   (default: `0.01`)
#' @param FishLenComps_Type Array [n_yrs × n_fish_fleets]
#'   (default: `2` = joint by sex, split by region)
#'   \itemize{
#'     \item \code{0}: Aggregated
#'     \item \code{1}: Split by sex and region
#'     \item \code{2}: Joint by sex, split by region
#'     \item \code{999}: Not simulated
#'   }
#' @param catch_units Units of catch - Array [n_regions × n_fish_fleets]
#'   \itemize{
#'     \item \code{0}: Abundance
#'     \item \code{1}: Biomass (default)
#'   }
#' @export Setup_Sim_Fishing
#' @family Simulation Setup
Setup_Sim_Fishing <- function(sim_list,
                              ln_sigmaC = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets)),
                              catch_units = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
                              init_F_val = 0,
                              Fmort_input = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              fish_sel_input,
                              fish_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ObsFishIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_fish_fleets)),
                              fish_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
                              comp_fishage_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishAge_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishAge_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
                              comp_fishlen_like = rep(0, sim_list$n_fish_fleets),
                              ISS_FishLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
                              ln_FishLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets)),
                              ln_FishLen_theta_agg = rep(log(1), sim_list$n_fish_fleets),
                              FishLen_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
                              FishLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_fish_fleets, 2)),
                              FishLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets))
                              ) {

  # Validate dimensions of all input parameters
  check_sim_dimensions(ln_sigmaC, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ln_sigmaC")
  check_sim_dimensions(catch_units, n_regions = sim_list$n_regions, n_fish_fleets = sim_list$n_fish_fleets, what = "catch_units")
  check_sim_dimensions(Fmort_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "Fmort_input")
  check_sim_dimensions(fish_sel_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "fish_sel_input")
  check_sim_dimensions(fish_q_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_fish_fleets = sim_list$n_fish_fleets, n_sims = sim_list$n_sims, what = "fish_q_input")
  check_sim_dimensions(ObsFishIdx_SE, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_fish_fleets = sim_list$n_fish_fleets, what = "ObsFishIdx_SE")
  check_sim_dimensions(fish_idx_type, n_regions = sim_list$n_regions, n_fish_fleets = sim_list$n_fish_fleets, what = "fish_idx_type")

  # Validate fishery age composition parameters
  check_sim_dimensions(comp_fishage_like, n_fish_fleets = sim_list$n_fish_fleets, what = "comp_fishage_like")
  check_sim_dimensions(ISS_FishAgeComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
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
  check_sim_dimensions(ISS_FishLenComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
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



#' Helper function to setup sigma F and sigma F agg mapping
#'
#' @param input_list Input list
#' @param sigmaF_spec Character vector for sigmaF specification
#' @param sigmaF_agg_spec Character vector for sigmaF_agg specification
#' @keywords internal
do_sigmaF_mapping <- function(input_list, sigmaF_spec, sigmaF_agg_spec) {

  # Sigma F -----------------------------------------------------------------
  map_sigmaF <- input_list$par$ln_sigmaF # initialize

  # Same sigmaF across fleets, but unique across regions
  if(sigmaF_spec == "est_shared_f") {
    map_sigmaF[1:input_list$data$n_regions,] <- 1:input_list$data$n_regions
    input_list$map$ln_sigmaF <- factor(map_sigmaF)
  }

  # Same sigmaF across regions, but unique across fleets
  if(sigmaF_spec == "est_shared_r") {
    map_sigmaF[,1:input_list$data$n_fish_fleets] <- 1:input_list$data$n_fish_fleets
    input_list$map$ln_sigmaF <- factor(map_sigmaF)
  }

  # Same sigmaF across regions and fleets
  if(sigmaF_spec == "est_shared_r_f") input_list$map$ln_sigmaF <- factor(rep(1, length(input_list$par$ln_sigmaF)))

  # Fixing sigmaF
  if(sigmaF_spec == "fix") input_list$map$ln_sigmaF <- factor(rep(NA, length(input_list$par$ln_sigmaF)))

  # Estimating all sigmaF
  if(sigmaF_spec == "est_all") input_list$map$ln_sigmaF <- factor(1:length(input_list$par$ln_sigmaF))

  # Print Message
  collect_message("sigmaF is specified as: ", sigmaF_spec)


  # Sigma F Aggregated ------------------------------------------------------

  # Validate Options
  if(input_list$data$est_all_regional_F == 1 && sigmaF_agg_spec != 'fix') stop("Fishing mortality is specified to be estimated for all regions, but sigmaF_agg_spec is not specified at `fix`!")

  # Process error for aggregated fishing mortality / catch options
  if(input_list$data$est_all_regional_F == 0) {

    # Same sigmaF across fleets
    if(sigmaF_agg_spec == "est_shared_f") input_list$map$ln_sigmaF <- factor(rep(1, length(input_list$par$ln_sigmaF_agg)))

    # Fixing ln_sigmaF_agg
    if(sigmaF_agg_spec == "fix") input_list$map$ln_sigmaF_agg <- factor(rep(NA, length(input_list$par$ln_sigmaF_agg)))

    # Estimating all ln_sigmaF_agg
    if(sigmaF_agg_spec == "est_all") input_list$map$ln_sigmaF_agg <- factor(1:length(input_list$par$ln_sigmaF_agg))

  } # end if some fishing mortality deviations are aggregated

  # If all fishing mortality deviations are regional, then don't estimate this
  if(input_list$data$est_all_regional_F == 1) input_list$map$ln_sigmaF_agg <- factor(rep(NA, length(input_list$par$ln_sigmaF_agg)))

  # Print message
  collect_message("sigmaF_agg is specified as: ", sigmaF_agg_spec)

  return(input_list)
}

#' Helper function to setup sigma C and sigma C agg mapping
#'
#' @param input_list Input list
#' @param sigmaC_spec Character vector for sigmaC specification
#' @param sigmaC_agg_spec Character vector for sigmaC_agg specification
#' @keywords internal
do_sigmaC_mapping <- function(input_list, sigmaC_spec, sigmaC_agg_spec) {

  # Sigma C -----------------------------------------------------------------
  map_sigmaC <- input_list$par$ln_sigmaC # initialize

  # Same sigmaC across fleets, but unique across regions and years
  if(sigmaC_spec == "est_shared_f") {
    map_sigmaC[1:input_list$data$n_regions,1:length(input_list$data$years),] <- 1:(input_list$data$n_regions * length(input_list$data$years))
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across regions, but unique across fleets and years
  if(sigmaC_spec == "est_shared_r") {
    map_sigmaC[,1:length(input_list$data$years),1:input_list$data$n_fish_fleets] <- 1:(input_list$data$n_fish_fleets * length(input_list$data$years))
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across years, but unique across regions and fleets
  if(sigmaC_spec == "est_shared_y") {
    counter <- 1
    for(r in 1:input_list$data$n_regions) {
      for(f in 1:input_list$data$n_fish_fleets) {
        map_sigmaC[r, , f] <- counter
        counter <- counter + 1
      }
    }
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across regions and fleets, but unique across years
  if(sigmaC_spec == "est_shared_r_f") {
    map_sigmaC[,1:length(input_list$data$years),] <- 1:length(input_list$data$years)
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across fleets and years, but unique across regions
  if(sigmaC_spec == "est_shared_f_y") {
    for(r in 1:input_list$data$n_regions) map_sigmaC[r, , ] <- r
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across regions and years, but unique across fleets
  if(sigmaC_spec == "est_shared_r_y") {
    for(f in 1:input_list$data$n_fish_fleets) map_sigmaC[, , f] <- f
    input_list$map$ln_sigmaC <- factor(map_sigmaC)
  }
  # Same sigmaC across regions, years, and fleets
  if(sigmaC_spec == "est_shared_r_y_f") {
    input_list$map$ln_sigmaC <- factor(rep(1, length(input_list$par$ln_sigmaC)))
  }
  # Fixing sigmaC
  if(sigmaC_spec == "fix") {
    input_list$map$ln_sigmaC <- factor(rep(NA, length(input_list$par$ln_sigmaC)))
  }
  # Estimating all sigmaC
  if(sigmaC_spec == "est_all") {
    input_list$map$ln_sigmaC <- factor(1:length(input_list$par$ln_sigmaC))
  }
  # Print Message
  collect_message("sigmaC is specified as: ", sigmaC_spec)


  # Sigma C Aggregated ------------------------------------------------------
  map_sigmaC_agg <- input_list$par$ln_sigmaC_agg # initialize

  # Validate Options
  if(input_list$data$est_all_regional_F == 1 && sigmaC_agg_spec != 'fix') stop("Catch is avaliable in all regions and periods, but sigmaC_agg_spec is not specified at `fix`!")

  # Process error for aggregated fishing mortality / catch options
  if(input_list$data$est_all_regional_F == 0) {
    # Same sigmaC across fleets, but unique across regions and years
    if(sigmaC_agg_spec == "est_shared_f") {
      map_sigmaC_agg[1:length(input_list$data$years),] <- 1:(length(input_list$data$years))
      input_list$map$ln_sigmaC_agg <- factor(map_sigmaC_agg)
    }
    # Same sigmaC_agg across years, but unique across fleets
    if(sigmaC_agg_spec == "est_shared_y") {
      for(f in 1:input_list$data$n_fish_fleets) {
        map_sigmaC_agg[, f] <- f
      }
      input_list$map$ln_sigmaC_agg <- factor(map_sigmaC_agg)
    }
    # Same sigmaC_agg across years and fleets
    if(sigmaC_agg_spec == "est_shared_y_f") {
      input_list$map$ln_sigmaC_agg <- factor(rep(1, length(input_list$par$ln_sigmaC_agg)))
    }
    # Fixing sigmaC_agg
    if(sigmaC_agg_spec == "fix") {
      input_list$map$ln_sigmaC_agg <- factor(rep(NA, length(input_list$par$ln_sigmaC_agg)))
    }
    # Estimating all sigmaC_agg
    if(sigmaC_agg_spec == "est_all") {
      input_list$map$ln_sigmaC_agg <- factor(1:length(input_list$par$ln_sigmaC_agg))
    }
  } # end if some fishing mortality deviations are aggregated

  # If all fishing mortality deviations are regional, then don't estimate this
  if(input_list$data$est_all_regional_F == 1) input_list$map$ln_sigmaC_agg <- factor(rep(NA, length(input_list$par$ln_sigmaC_agg)))

  # Print message
  collect_message("sigmaC_agg is specified as: ", sigmaC_agg_spec)

  return(input_list)
}

#' Helper function to setup fmort mapping
#'
#' @param input_list Input list
#' @keywords internal
do_Fmort_mapping <- function(input_list) {

  # Mapping for fishing mortality deviations
  F_dev_map <- input_list$par$ln_F_devs # initialize for mapping
  F_dev_map[] <- NA
  F_dev_counter <- 1

  for(r in 1:input_list$data$n_regions) {
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {

        # if we are not using catch or its specified as aggregated catch, then turn off estimation of these f devs
        if(input_list$data$UseCatch[r,y,f] == 0) {
          F_dev_map[r,y,f] <- NA
        }

        # if we are using regional catch or its specified as regional catch, then estimate the f devs
        if(input_list$data$UseCatch[r,y,f] == 1) {
          F_dev_map[r,y,f] <- F_dev_counter
          F_dev_counter <- F_dev_counter + 1
        } # end if

      } # end f
    } # end y
  } # end r

  # If we have regional aggregated catch but some F devs are not regional
  if(any(input_list$dat$Catch_Type == 0) && input_list$data$est_all_regional_F == 0) {
    Catch_Type_map <- input_list$data$Catch_Type # initialize for mapping
    # Loop through to find years and fleets with aggregated catch
    for(f in 1:input_list$data$n_fish_fleets) {
      agg_tmp <- which(Catch_Type_map[,f] == 0) # years with aggregated catch for a given fleet
      F_dev_map[,agg_tmp,f] <- NA # specify as NA for 0s
      F_dev_map[,-agg_tmp,f] <- 1:length(F_dev_map[,-agg_tmp,f]) # specify as unique values for others
    } # end f loop
  }

  input_list$map$ln_F_devs <- factor(F_dev_map)

  # If we are estimating some aggregated F devs
  if(input_list$data$est_all_regional_F == 0) {
    input_list$map$ln_F_devs_AggCatch <- factor(1:length(input_list$par$ln_F_devs_AggCatch))
    input_list$map$ln_F_mean_AggCatch <- factor(1:input_list$data$n_fish_fleets)
  }

  # If we are estimating regional F devs for all
  if(input_list$data$est_all_regional_F == 1) {
    input_list$map$ln_F_devs_AggCatch <- factor(array(NA, dim = dim(input_list$data$Catch_Type[rowSums(input_list$data$Catch_Type) == 0, , drop = FALSE])))
    input_list$map$ln_F_mean_AggCatch <- factor(rep(NA, input_list$data$n_fish_fleets))
  }

  return(input_list)
}

#' Setup fishing mortality and catch observations
#'
#' @param input_list A list containing data, parameters, and map lists used by the model.
#'
#' @param ObsCatch Numeric array of observed catches, dimensioned \code{[n_regions, n_years, n_fish_fleets]}.
#'
#' @param Catch_Type Integer matrix with dimensions \code{[n_years, n_fish_fleets]}, specifying catch data types:
#' \itemize{
#'   \item \code{0}: Use aggregated catch data for the year.
#'   \item \code{1}: Use region-specific catch data for the year.
#' }
#'
#' @param UseCatch Indicator array \code{[n_regions, n_years, n_fish_fleets]} specifying whether to include catch data in the fit:
#' \itemize{
#'   \item \code{0}: Do not use catch data.
#'   \item \code{1}: Use catch data and fit.
#' }
#'
#' @param Use_F_pen Integer flag indicating whether to apply a fishing mortality penalty:
#' \itemize{
#'   \item \code{0}: Do not apply penalty.
#'   \item \code{1}: Apply penalty.
#' }
#'
#' @param est_all_regional_F Integer flag indicating whether all regional fishing mortality deviations are estimated:
#' \itemize{
#'   \item \code{0}: Some fishing mortality deviations are aggregated across regions.
#'   \item \code{1}: All fishing mortality deviations are regional.
#' }
#'
#'
#' @param sigmaC_spec Character string specifying observation error structure for catch data. Default behavior fixes \code{sigmaC} at a starting value of \code{1e-3} (log-scale \code{ln_sigmaC = log(1e-3)}) for all regions, years, and fleets. Other options include:
#' \itemize{
#'   \item \code{"est_shared_f"}: Estimate \code{sigmaC} shared across fishery fleets, unique by region and year.
#'   \item \code{"est_shared_r"}: Estimate \code{sigmaC} shared across regions, unique by fleet and year.
#'   \item \code{"est_shared_y"}: Estimate \code{sigmaC} shared across years, unique by region and fleet.
#'   \item \code{"est_shared_r_f"}: Estimate \code{sigmaC} shared across regions and fleets, unique by year.
#'   \item \code{"est_shared_f_y"}: Estimate \code{sigmaC} shared across fleets and years, unique by region.
#'   \item \code{"est_shared_r_y"}: Estimate \code{sigmaC} shared across regions and years, unique by fleet.
#'   \item \code{"est_shared_r_y_f"}: Estimate single \code{sigmaC} shared across regions, years, and fleets.
#'   \item \code{"fix"}: Fix \code{sigmaC} at the starting value.
#'   \item \code{"est_all"}: Estimate separate \code{sigmaC} for each region, year, and fleet combination.
#' }
#'
#' @param sigmaF_spec Character string specifying process error structure for fishing mortality. Default fixes \code{sigmaF} at \code{1} on the log scale (i.e., \code{ln_sigmaF = 0}). Other options include:
#' \itemize{
#'   \item \code{"est_shared_f"}: Estimate \code{sigmaF} shared across fishery fleets.
#'   \item \code{"est_shared_r"}: Estimate \code{sigmaF} shared across regions but unique by fleet.
#'   \item \code{"est_shared_r_f"}: Estimate \code{sigmaF} shared across regions and fleets.
#'   \item \code{"fix"}: Fix \code{sigmaF} at the starting value.
#'   \item \code{"est_all"}: Estimate separate \code{sigmaF} for each region and fleet.
#' }
#'
#' @param sigmaF_agg_spec Character string specifying process error structure for aggregated fishing mortality. Default fixes \code{sigmaF_agg} at the starting value (log-scale \code{ln_sigmaF_agg}). Other options include:
#' \itemize{
#'   \item \code{"est_shared_f"}: Estimate \code{sigmaF_agg} shared across fishery fleets.
#'   \item \code{"fix"}: Fix at the starting value.
#'   \item \code{"est_all"}: Estimate separate parameters for each fishery fleet.
#' }
#'
#' @param ... Additional arguments specifying starting values for \code{ln_sigmaC}, \code{ln_sigmaF}, and \code{ln_sigmaF_agg}.
#' @param sigmaC_agg_spec Character string specifying process error structure for aggregated catch observation error. Default fixes \code{sigmaC_agg} at the starting value (log-scale \code{sigmaC_agg}). Other options include:
#' \itemize{
#'   \item \code{"est_shared_f"}: Estimate \code{sigmaC_agg} shared across fishery fleets, unique by year.
#'   \item \code{"est_shared_y"}: Estimate \code{sigmaC_agg} shared across years, unique by fleet.
#'   \item \code{"est_shared_y_f"}: Estimate single \code{sigmaC_agg} shared across years and fleets.
#'   \item \code{"fix"}: Fix at the starting value.
#'   \item \code{"est_all"}: Estimate separate parameters for each year and fleet combination.
#' }
#' @param catch_units Catch units - Array dimensioned by n_regions x n_fish_fleets
#' \itemize{
#'   \item \code{"abd"}: Catch units in abundance
#'   \item \code{"biom"}: Catch units in biomass (default)
#' }
#'
#' @export Setup_Mod_Catch_and_F
#' @family Model Setup
Setup_Mod_Catch_and_F <- function(input_list,
                                  ObsCatch,
                                  catch_units = array("biom", dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets)),
                                  Catch_Type = array(1, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets)),
                                  UseCatch,
                                  Use_F_pen = 1,
                                  est_all_regional_F = 1,
                                  sigmaC_spec = "fix",
                                  sigmaC_agg_spec = "fix",
                                  sigmaF_spec = "fix",
                                  sigmaF_agg_spec = "fix",
                                  ...
                                  ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Catch objects
  check_data_dimensions(ObsCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsCatch')
  check_data_dimensions(Catch_Type, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'Catch_Type')
  check_data_dimensions(UseCatch, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseCatch')
  check_data_dimensions(catch_units, n_regions = input_list$data$n_region, n_fish_fleets = input_list$data$n_fish_fleets, what = 'catch_units')

  # Indicators for whether catch is aggregated across regions
  if(est_all_regional_F == 0 && any(unique(Catch_Type) == 0)) collect_message("Catch is aggregated by region in some years, with a separate aggregated ln_F_Mean and ln_F_devs estimated in those years")
  if(est_all_regional_F == 1 && any(unique(Catch_Type) == 0)) collect_message("Catch is aggregated by region in some years, with a region specific ln_F_Mean and ln_F_devs estiamted, where these fishing mortalities are estimated using information from data (age and indices) in subsequent years")
  if(est_all_regional_F == 1 && any(unique(Catch_Type) == 1)) collect_message("Catch is region specific, with region specific ln_F_Mean and ln_F_devs")
  if(any(UseCatch == 0)) collect_message("User specified catch for some years and fleets to not be fit to, and ln_F_devs will not be estimated for those dimensions")

  # Fishing Mortality checking
  if(!est_all_regional_F %in% c(0,1)) stop("est_all_regional_F incorrectly specified. Either set at 0 (not all regional Fs are estiamted) or 1 (all regional Fs are estimated)")
  else collect_message("Fishing mortality is estimated: ", ifelse(est_all_regional_F == 0, 'Not For All Regions', "For All Regions"))
  if(!Use_F_pen %in% c(0,1)) stop("Use_F_pen incorrectly specified. Either set at 0 (don't use F penalty) or 1 (use F penalty)")
  else collect_message("Fishing mortality penalty is: ", ifelse(Use_F_pen == 0, 'Not Used', "Used"))

  if(sigmaC_spec == "fix" && !("ln_sigmaC" %in% names(starting_values))) warning("sigmaC is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaC_agg_spec == "fix" && !("ln_sigmaC_agg" %in% names(starting_values))) warning("sigmaF_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaF_spec == "fix" && !("ln_sigmaF" %in% names(starting_values))) warning("sigmaF_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")
  if(sigmaF_agg_spec == "fix" && !("ln_sigmaF_agg" %in% names(starting_values))) warning("sigmaF_agg_spec is specified as fix, but no starting values / fixed values are provided. Either do this post-hoc, or use the ... argument if you do not want to use default values")

  # Catch units
  catch_units[catch_units == 'abd'] <- 0
  catch_units[catch_units == 'biom'] <- 1
  catch_units <- array(as.numeric(catch_units), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets)) # convert to numeric array

  # Populate Data List ------------------------------------------------------

  input_list$data$ObsCatch <- ObsCatch
  input_list$data$Catch_Type <- Catch_Type
  input_list$data$UseCatch <- UseCatch
  input_list$data$est_all_regional_F <- est_all_regional_F
  input_list$data$Use_F_pen <- Use_F_pen
  input_list$data$catch_units <- catch_units

  # Populate Parameter List -------------------------------------------------

  # Catch observation error
  if("ln_sigmaC" %in% names(starting_values)) input_list$par$ln_sigmaC <- starting_values$ln_sigmaC
  else input_list$par$ln_sigmaC <- array(log(0.01), dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))

  # Catch observation error for aggregated catch
  if("ln_sigmaC_agg" %in% names(starting_values)) input_list$par$ln_sigmaC_agg <- starting_values$ln_sigmaC_agg
  else input_list$par$ln_sigmaC_agg <- array(log(0.01), dim = c(length(input_list$data$years), input_list$data$n_fish_fleets))

  # Process error fishing deviations for regional catch
  if("ln_sigmaF" %in% names(starting_values)) input_list$par$ln_sigmaF <- starting_values$ln_sigmaF
  else input_list$par$ln_sigmaF <- array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))

  # Process error fishing deviations for aggregated catch
  if("ln_sigmaF_agg" %in% names(starting_values)) input_list$par$ln_sigmaF <- starting_values$ln_sigmaF
  else input_list$par$ln_sigmaF_agg <- rep(log(1), input_list$data$n_fish_fleets)

  # Log mean fishing mortality
  if("ln_F_mean" %in% names(starting_values)) input_list$par$ln_F_mean <- starting_values$ln_F_mean
  else input_list$par$ln_F_mean <- array(log(0.1), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))

  # Log fishing deviations
  if("ln_F_devs" %in% names(starting_values)) input_list$par$ln_F_devs <- starting_values$ln_F_devs
  else input_list$par$ln_F_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))

  # Log mean fishing mortality for aggregated catch
  if("ln_F_mean_AggCatch" %in% names(starting_values)) input_list$par$ln_F_mean_AggCatch <- starting_values$ln_F_mean_AggCatch
  else input_list$par$ln_F_mean_AggCatch <- rep(log(0.1), dim(input_list$data$Catch_Type[rowSums(input_list$data$Catch_Type) == 0, , drop = FALSE])[2])

  # Log fishing deviations for aggregated catch
  if("ln_F_devs_AggCatch" %in% names(starting_values)) input_list$par$ln_F_devs_AggCatch <- starting_values$ln_F_devs_AggCatch
  else input_list$par$ln_F_devs_AggCatch <- array(0, dim = dim(input_list$data$Catch_Type[rowSums(input_list$data$Catch_Type) == 0, , drop = FALSE]))


  # Mapping Options ---------------------------------------------------------
  input_list <- do_sigmaC_mapping(input_list, sigmaC_spec, sigmaC_agg_spec)
  input_list <- do_sigmaF_mapping(input_list, sigmaF_spec, sigmaF_agg_spec)
  input_list <- do_Fmort_mapping(input_list)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}

#' Helper function to set up fishery age overdispersion parameters
#'
#' @param input_list Input list
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
    if(input_list$data$FishAgeComps_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps[,,f]) == 0) {
      map_FishAge_theta[,,f] <- NA
      map_FishAge_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishAge_theta <- factor(map_FishAge_theta)
  input_list$map$ln_FishAge_theta_agg <- factor(map_FishAge_theta_agg)

  return(input_list)
}

#' Helper function to set up fishery length overdispersion parameters
#'
#' @param input_list Input list
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
    if(input_list$data$FishLenComps_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps[,,f]) == 0) {
      map_FishLen_theta[,,f] <- NA
      map_FishLen_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_FishLen_theta <- factor(map_FishLen_theta)
  input_list$map$ln_FishLen_theta_agg <- factor(map_FishLen_theta_agg)

  return(input_list)
}

#' Helper function to set up fishery age overdispersion correlation parameters
#'
#' @param input_list Input list
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
    if(input_list$data$FishAgeComps_LikeType[f] == 0 || sum(input_list$data$UseFishAgeComps[,,f]) == 0) {
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

#' Helper function to set up fishery length overdispersion correlation parameters
#'
#' @param input_list Input list
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
    if(input_list$data$FishLenComps_LikeType[f] == 0 || sum(input_list$data$UseFishLenComps[,,f]) == 0) {
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

#' Setup observed fishery indices and composition data (age and length comps)
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param ObsFishIdx Observed fishery index data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_fish_fleets]}.
#'
#' @param ObsFishIdx_SE Standard errors associated with \code{ObsFishIdx},
#' also dimensioned \code{[n_regions, n_years, n_fish_fleets]}.
#'
#' @param UseFishIdx Logical or binary indicator array (\code{[n_regions, n_years, n_fish_fleets]})
#' specifying whether to include a fishery index in the likelihood (\code{1}) or ignore it (\code{0}).
#'
#' @param ObsFishAgeComps Observed fishery age composition data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_ages, n_sexes, n_fish_fleets]}. Values should reflect counts or proportions
#' (not required to sum to 1, but should be on a comparable scale).
#'
#' @param UseFishAgeComps Indicator array (\code{[n_regions, n_years, n_fish_fleets]}) specifying whether
#' to fit fishery age composition data (\code{1}) or ignore it (\code{0}).
#'
#' @param ObsFishLenComps Observed fishery length composition data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_lens, n_sexes, n_fish_fleets]}. Values should reflect counts or proportions.
#'
#' @param UseFishLenComps Indicator array (\code{[n_regions, n_years, n_fish_fleets]}) specifying whether
#' to fit fishery length composition data (\code{1}) or ignore it (\code{0}).
#'
#' @param FishAgeComps_LikeType Character vector of length \code{n_fish_fleets} specifying the likelihood
#' type used for fishery age composition data. Options include \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#' and \code{"iid-Logistic-Normal"}. Use \code{"none"} to omit the likelihood.
#'
#' @param FishLenComps_LikeType Same as \code{FishAgeComps_LikeType}, but for fishery length composition data.
#'
#' @param FishAgeComps_Type Character vector specifying how age compositions are structured by fleet and year range.
#' Options include:
#' \itemize{
#'   \item \code{"agg"}: Aggregated across regions and sexes.
#'   \item \code{"spltRspltS"}: Split by region and by sex (compositions sum to 1 within region-sex group).
#'   \item \code{"spltRjntS"}: Split by region but summed jointly across sexes.
#'   \item \code{"none"}: No composition data used.
#' }
#' Format each element as \code{"<type>_Year_<start>-<end>_Fleet_<fleet number>"}
#' (e.g., \code{"agg_Year_1-10_Fleet_1"}).
#'
#' @param FishLenComps_Type Same as \code{FishAgeComps_Type}, but for length compositions.
#' @param fish_idx_type Character vector of length \code{n_fish_fleets} specifying the type of index data.
#' Options are \code{"abd"} for abundance, \code{"biom"} for biomass, and \code{"none"} if no index is available.
#'
#' @param ISS_FishAgeComps Input sample size for age compositions, array dimensioned
#' \code{[n_regions, n_years, n_sexes, n_fish_fleets]}. Required if observed age comps are normalized
#' (i.e., sum to 1), to correctly scale the contribution to the likelihood.
#'
#' @param ISS_FishLenComps Same as \code{ISS_FishAgeComps}, but for length compositions.
#'
#' @param ... Additional arguments specifying starting values for overdispersion parameters
#' (e.g., \code{ln_FishAge_theta}, \code{ln_FishLen_theta}, \code{ln_FishAge_theta_agg}, \code{ln_FishLen_theta_agg}).
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
                                        ISS_FishAgeComps,
                                        ObsFishLenComps,
                                        UseFishLenComps,
                                        ISS_FishLenComps,
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
  check_data_dimensions(ObsFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx')
  check_data_dimensions(ObsFishIdx_SE, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishIdx_SE')
  check_data_dimensions(UseFishIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishIdx')
  check_data_dimensions(fish_idx_type, n_fish_fleets = input_list$data$n_fish_fleets, what = 'fish_idx_type')
  if(!all(fish_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for fish_idx_type. Should be either abd, biom, or none")

  # Fishery compositions
  check_data_dimensions(ObsFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishAgeComps')
  check_data_dimensions(UseFishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishAgeComps')
  check_data_dimensions(UseFishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_fish_fleets = input_list$data$n_fish_fleets, what = 'UseFishLenComps')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsFishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ObsFishLenComps')
  if(!is.null(ISS_FishAgeComps)) check_data_dimensions(ISS_FishAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishAgeComps')
  if(!is.null(ISS_FishLenComps)) check_data_dimensions(ISS_FishLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_fish_fleets = input_list$data$n_fish_fleets, what = 'ISS_FishLenComps')
  check_data_dimensions(FishAgeComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishAgeComps_LikeType')
  check_data_dimensions(FishLenComps_LikeType, n_fish_fleets = input_list$data$n_fish_fleets, what = 'FishLenComps_LikeType')
  if(!all(FishAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(FishLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for FishLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")


  # Fishery Index Options ---------------------------------------------------

  fish_idx_type_vals <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets))
  for(f in 1:ncol(fish_idx_type_vals)) {
    if(fish_idx_type[f] == 'biom') fish_idx_type_vals[,f] <- 1 # biomass
    if(fish_idx_type[f] == 'abd') fish_idx_type_vals[,f] <- 0 # abundance
    if(fish_idx_type[f] == 'none') fish_idx_type_vals[,f] <- 999 # none
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
    ISS_FishAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        # if aggregated across sexes and regions (0) or joint across sexes
        if(FishAgeComps_Type_Mat[y,f] == 0) ISS_FishAgeComps[1,y,1,f] <- sum(ObsFishAgeComps[,y,,,f])
        # if split by region and sex
        if(FishAgeComps_Type_Mat[y,f] == 1) ISS_FishAgeComps[,y,,f] <- apply(ObsFishAgeComps[,y,,,f, drop = FALSE], c(1,4), sum)
        # if split by region, joint by sex
        if(FishAgeComps_Type_Mat[y,f] == 2) ISS_FishAgeComps[,y,1,f] <- apply(ObsFishAgeComps[,y,,,f, drop = FALSE], 1, sum)
      } # end f loop
    } # end y loop
  }

  # Fishery Lengths
  if(is.null(ISS_FishLenComps)) {
    collect_message("No ISS is specified for FishLenComps. ISS weighting is calculated by summing up values from ObsFishLenComps each year")
    ISS_FishLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_fish_fleets) {
        # if aggregated across sexes and regions (0)
        if(FishLenComps_Type_Mat[y,f] == 0) ISS_FishLenComps[1,y,1,f] <- sum(ObsFishLenComps[,y,,,f])
        # if split by region and sex
        if(FishLenComps_Type_Mat[y,f] == 1) ISS_FishLenComps[,y,,f] <- apply(ObsFishLenComps[,y,,,f, drop = FALSE], c(1,4), sum)
        # if split by region, joint by sex
        if(FishLenComps_Type_Mat[y,f] == 2) ISS_FishLenComps[,y,1,f] <- apply(ObsFishLenComps[,y,,,f, drop = FALSE], 1, sum)
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

#' Helper function to set up fishery selctivity fixed effects mapping
#'
#' @param input_list Input list
#' @param fish_fixed_sel_pars_spec Character vector specifying fishery selectivity fixed effects parameterization
#' @keywords internal
#' Helper function to set up fishery selctivity fixed effects mapping
#'
#' @param input_list Input list
#' @param fish_fixed_sel_pars_spec Character vector specifying fishery selectivity fixed effects parameterization
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
      if(sum(input_list$data$UseCatch[r,,f]) > 0) {

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


#' Helper function to set up fishery catchability mapping
#'
#' @param input_list Input list
#' @param fish_q_spec Character vector specifying fishery catchability parameterization
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

      if(sum(input_list$data$UseFishIdx[r,,f]) == 0) {
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

#' Helper function to set up fishery process error mapping
#'
#' @param input_list Input list
#' @param corr_opt_semipar Character vector specifying correlation parameter options
#' @param fishsel_pe_pars_spec Character vector specifying fishery process error parameterization
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
      if(input_list$data$cont_tv_fish_sel[r,f] == 0 || sum(input_list$data$UseCatch[r,,f]) == 0) {
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

#' Helper function to set up fishery selectivity deviations mapping
#'
#' @param input_list Input list
#' @param fishsel_devs_shared_ages List object for specifying which ages are shared when selectivity deviations are semi-parametric (e.g., list(1:5, 6:10, 11:30) specifies that ages 1-5, 6-10, and 11-30 have the same deviations.)
#' @param fish_sel_devs_spec Character vector specifying fishery selectivity deviations parameterization
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
          if(input_list$data$cont_tv_fish_sel[r,f] == 0 || sum(input_list$data$UseCatch[r,,f]) == 0) {
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

#' Setup fishery selectivity and catchability specifications
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param cont_tv_fish_sel Character vector specifying the form of continuous time-varying selectivity for each fishery fleet.
#' The vector must be length \code{n_fish_fleets}, and each element must follow the structure:
#' \code{"<time variation type>_Fleet_<fleet number>"}.
#'
#' Valid time variation types include:
#' \itemize{
#'   \item \code{"none"}: No continuous time variation (default)
#'   \item \code{"iid"}: Independent and identically distributed deviations across years.
#'   \item \code{"rw"}: Random walk in time.
#'   \item \code{"3dmarg"}: 3D marginal time-varying selectivity.
#'   \item \code{"3dcond"}: 3D conditional time-varying selectivity.
#'   \item \code{"2dar1"}: Two-dimensional AR1 process.
#' }
#'
#' For example:
#' \itemize{
#'   \item \code{"iid_Fleet_1"} applies an iid time-varying structure to Fleet 1.
#'   \item \code{"none_Fleet_2"} means no time variation is used for Fleet 2.
#' }
#'
#' @param fish_sel_blocks Character vector specifying the fishery selectivity blocks for each region and fleet.
#'
#' Each element must follow one of the following structures:
#' \itemize{
#'   \item `"Block_<block number>_Year_<start>-<end>_Fleet_<fleet number>"`
#'   \item `"Block_<block number>_Year_<start>-terminal_Fleet_<fleet number>"`
#'   \item `"none_Fleet_<fleet number>"`
#' }
#'
#' This argument defines how fishery selectivity varies over time for each fleet:
#' \itemize{
#'   \item \code{"Block_..."} entries specify discrete time blocks during which selectivity parameters are assumed constant.
#'   \item \code{"none_..."} entries indicate that selectivity is constant across all years for the specified fleet.
#' }
#'
#' If time-block-based selectivity is specified for a fleet (via \code{fish_sel_blocks}), its corresponding continuous selectivity option (in \code{cont_tv_fish_sel}) must be set to \code{"none_Fleet_<fleet number>"}. The two approaches—blocked and continuous time-varying selectivity—are mutually exclusive.
#' The default for each fleet is \code{"none_Fleet_x"} (i.e., no selectivity blocks).
#'
#' @param fish_sel_model Character vector specifying the fishery selectivity functional form for each fleet, and optionally by time block.
#'
#' Each element must follow one of the following structures:
#' \itemize{
#'   \item \code{"<selectivity model>_Fleet_<fleet number>"}
#'   \item \code{"<selectivity model>_Fleet_<fleet number>_Block_<block number>"}
#' }
#'
#' The first form applies a single selectivity model across all years for the specified fleet.
#' The second form allows the user to assign a distinct selectivity model to a specific time block, as defined in \code{fish_sel_blocks}.
#'
#' Available selectivity model types include:
#' \itemize{
#'   \item \code{"logist1"} — Logistic function with parameters \code{a50} and \code{k}.
#'   \item \code{"logist2"} — Logistic function with parameters \code{a50} and \code{a95}.
#'   \item \code{"gamma"} — Dome-shaped gamma function with parameters \code{amax} and \code{delta}.
#'   \item \code{"exponential"} — Exponential function with a power parameter.
#'   \item \code{"dbnrml"} — Double-normal function with six parameters.
#' }
#' If multiple selectivity time blocks are specified for a fleet (using \code{fish_sel_blocks}), then the corresponding selectivity model for each block must be explicitly defined using the \code{"<model>_Block_<block>_Fleet_<fleet>"} format.
#' If blocks are not defined for a fleet, use the \code{"<model>_Fleet_<fleet number>"} format only.
#' For mathematical definitions and implementation details of each selectivity form, refer to the model equations vignette.
#'
#' @param fish_q_blocks Character vector specifying fishery catchability (q) blocks for each fleet.
#' Each element must follow the structure: \code{"Block_<block number>_Year_<start>-<end>_Fleet_<fleet number>"}
#' or \code{"none_Fleet_<fleet number>"}. Default is "none_Fleet_x".
#'
#' This allows users to define time-varying catchability blocks independently of selectivity blocks.
#' The blocks must be non-overlapping and sequential in time within each fleet.
#'
#' For example:
#' \itemize{
#'   \item \code{"Block_1_Year_1-35_Fleet_1"} assigns block 1 to Fleet 1 for years 1–35.
#'   \item \code{"Block_2_Year_36-56_Fleet_1"} continues with block 2 for years 36–56.
#'   \item \code{"Block_3_Year_57-terminal_Fleet_1"} assigns block 3 from year 57 to the terminal year for Fleet 1.
#'   \item \code{"none_Fleet_2"} indicates no catchability blocks are used for Fleet 2.
#' }
#'
#' Internally, these specifications are converted to a \code{[n_regions, n_years, n_fish_fleets]} array,
#' where each block is mapped to the appropriate years and fleets.
#' @param fishsel_pe_pars_spec Character string specifying how process error parameters for fishery selectivity
#' are estimated across regions and sexes. This is only relevant if \code{cont_tv_fish_sel} is not set to \code{"none"};
#' otherwise, all process error parameters are treated as fixed.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates separate process error parameters for each region and sex.
#'   \item \code{"est_shared_r"}: Shares process error parameters across regions (sex-specific parameters are still estimated).
#'   \item \code{"est_shared_s"}: Shares process error parameters across sexes (region-specific parameters are still estimated).
#'   \item \code{"est_shared_r_s"}: Shares process error parameters across both regions and sexes, estimating a single set of parameters.
#'   \item \code{"est_shared_f_x"}: Shares process error parameters with another fleet, where \code{x} is the fleet number to share with.
#'     This option forces multiple fleets to have identical process error variance and correlation structures for their
#'     time-varying selectivity. For example, \code{"est_shared_f_2"} means the current fleet will use the same
#'     process error parameters as fleet 2. The reference fleet (fleet x) must use one of the other sharing options
#'     and cannot itself be sharing with another fleet.
#'   \item \code{"fix"} or \code{"none"}: Does not estimate process error parameters; all are treated as fixed.
#' }
#' @param fish_fixed_sel_pars_spec Character string specifying the structure for estimating
#' fixed-effect parameters of the fishery selectivity model (e.g., a50, k, amax).
#' This controls whether selectivity parameters are estimated separately or shared across regions and sexes.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates separate fixed-effect selectivity parameters for each region and sex.
#'   \item \code{"est_shared_r"}: Shares parameters across regions (sex-specific parameters are still estimated).
#'   \item \code{"est_shared_s"}: Shares parameters across sexes (region-specific parameters are still estimated).
#'   \item \code{"est_shared_r_s"}: Shares parameters across both regions and sexes, estimating a single set of fixed-effect parameters.
#'   \item \code{"est_shared_f_x"}: Shares fixed-effect selectivity parameters with another fleet, where \code{x} is the fleet number to share with.
#'     This option forces multiple fleets to have identical selectivity curves by using the same underlying parameters
#'     (e.g., same a50, k, amax values). For example, \code{"est_shared_f_2"} means the current fleet will use the same
#'     fixed-effect selectivity parameters as fleet 2. The reference fleet (fleet x) must use one of the other sharing options
#'     and cannot itself be sharing with another fleet.
#'   \item \code{"fix"}: Fixes all selectivity parameters to their initial values (no estimation).
#'   \item \code{"none"}: No selectivity parameters are estimated (equivalent to \code{"fix"}).
#' }
#' @param fish_q_spec Character string specifying the structure of fishery catchability (\code{q}) estimation
#' across regions. This controls whether separate or shared parameters are used.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates separate catchability parameters for each region.
#'   \item \code{"est_shared_r"}: Estimates a single catchability parameter shared across all regions.
#' }
#' @param fish_sel_devs_spec Character string specifying the structure of process error deviations
#' in time-varying fishery selectivity dimensioned by the number of fishery fleets. This determines how deviations are estimated across regions and sexes.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates a separate deviation time series for each region and sex.
#'   \item \code{"est_shared_r"}: Shares deviations across regions (sex-specific deviations are still estimated).
#'   \item \code{"est_shared_s"}: Shares deviations across sexes (region-specific deviations are still estimated).
#'   \item \code{"est_shared_r_s"}: Shares deviations across both regions and sexes, estimating a single deviation time series.
#'   \item \code{"est_shared_a"}: Shares deviations across age blocks.
#'   \item \code{"est_shared_r_a"}: Shares deviations across regions and age shared blocks.
#'   \item \code{"est_shared_a_s"}: Shares deviations across age shared blocks and sexes.
#'   \item \code{"est_shared_r_a_s"}: Shares deviations across regions, age shared blocks, and sexes.
#'   \item \code{"est_shared_f_x"}: Shares deviations with another fleet, where \code{x} is the fleet number to share with.
#'     This option allows multiple fleets to use identical deviation parameters, reducing the number of parameters
#'     to estimate. For example, \code{"est_shared_f_2"} means the current fleet will use the same deviation
#'     parameters as fleet 2. The reference fleet (fleet x) must use one of the other sharing options
#'     (\code{"est_all"}, \code{"est_shared_r"}, \code{"est_shared_s"}, or \code{"est_shared_r_s"})
#'     and cannot itself be sharing with another fleet.
#'   \item \code{"fix"}: Fixes all deviation parameters to zero (no time-variation).
#'   \item \code{"none"}: No deviation parameters are estimated (equivalent to \code{"fix"}).
#' }
#'
#' This argument is only used when a continuous time-varying selectivity form is specified (e.g., via \code{cont_tv_fish_sel}).
#' @param corr_opt_semipar Character string specifying which correlation structures to suppress
#'   when using semi-parametric time-varying selectivity models. Only used if \code{cont_tv_sel}
#'   is set to one of \code{"3dmarg"}, \code{"3dcond"}, or \code{"2dar1"}.
#'
#'   This option allows users to turn off estimation of specific correlation components in the
#'   time-varying selectivity model. This can improve stability or enforce assumptions about
#'   independence in the temporal or age structure.
#'
#'   Available options:
#'   \itemize{
#'     \item \code{"corr_zero_y"}: Sets year (temporal) correlations to 0.
#'     \item \code{"corr_zero_b"}: Sets age correlations to 0.
#'     \item \code{"corr_zero_y_b"}: Sets both year and bin correlations to 0.
#'     \item \code{"corr_zero_c"}: Sets cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_y_c"}: Sets year and cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_b_c"}: Sets bin (age) and cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_y_b_c"}: Sets all correlations (year, bin (age), and cohort) to 0.
#'       Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}; equivalent to an iid structure.
#'   }
#'
#' These correlation-suppression flags are ignored when \code{cont_tv_sel} is set to any other value.
#' @param Use_fish_q_prior Integer (0 or 1). Flag to enable/disable fishery catchability priors.
#'   When set to 1, applies log-normal priors to fishery selectivity parameters as specified
#'   in \code{fish_q_prior}. When set to 0, no priors are applied.
#' @param fish_q_prior Data frame containing prior specifications for fishery catchability parameters.
#'   Must include columns: \code{region} (region index), \code{fleet} (fleet index),
#'   \code{block} (time block index), \code{mu} (prior mean on natural scale), and \code{sd} (prior standard deviation on log scale).
#'   Each row specifies a log-normal prior N(log(mu), sd) for a given catchability parameter.
#'   Only parameters with rows in this data frame will have priors applied.
#' @param Use_fish_selex_prior Integer (0 or 1). Flag to enable/disable fishery selectivity priors.
#'   When set to 1, applies log-normal priors to fishery selectivity parameters as specified
#'   in \code{fish_selex_prior}. When set to 0, no priors are applied.
#' @param fish_selex_prior Data frame containing prior specifications for fishery selectivity parameters.
#'   Must include columns: \code{region} (region index), \code{fleet} (fleet index),
#'   \code{block} (time block index), \code{sex} (sex index), \code{par} (parameter index),
#'   \code{mu} (prior mean on natural scale), and \code{sd} (prior standard deviation on log scale).
#'   Each row specifies a log-normal prior N(log(mu), sd) for one selectivity parameter.
#'   Only parameters with rows in this data frame will have priors applied.
#' @param cont_tv_fish_sel_penalty Whether or not continuous fishery time varying selectivity penalties are applied (if cont_tv_fish_sel > 0)
#' @param fishsel_devs_shared_ages List object for specifying which ages are shared when selectivity deviations are semi-parametric (e.g., list(1:5, 6:10, 11:30) specifies that ages 1-5, 6-10, and 11-30 have the same deviations.)
#'
#' @export Setup_Mod_Fishsel_and_Q
#'
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
    tmp_fleet <- as.numeric(tmp_sel_form_vec[3])
    # get block index
    tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[5])

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


