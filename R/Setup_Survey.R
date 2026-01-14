#' Setup values for survey parameterization
#'
#' @param sim_list Simulation list object from `Setup_Sim_Dim()`
#' @param ObsSrvIdx_SE Survey index observation error
#'   [n_regions × n_yrs × n_srv_fleets]
#'   (default: `0.2`)
#' @param srv_sel_input Survey selectivity array
#'   [n_regions × n_yrs × n_ages × n_sexes × n_srv_fleets × n_sims]
#'   (no default, must be provided)
#' @param srv_q_input Survey catchability array
#'   [n_regions × n_yrs × n_srv_fleets × n_sims]
#'   (default: `1`)
#' @param t_srv Survey timing fraction
#'   [n_regions × n_srv_fleets]
#'   (default: `0`)
#' @param srv_idx_type Array of index types [n_regions x n_srv_fleets]
#'   (default: all `1` = biomass index)
#'   \itemize{
#'     \item \code{0}: Abundance index
#'     \item \code{1}: Biomass index
#'   }
#' @param comp_srvage_like Vector [n_srv_fleets] specifying likelihood for simulating age comps
#'   (default: all `0` = multinomial)
#'   \itemize{
#'     \item \code{0}: Multinomial
#'     \item \code{1}: Dirichlet-Multinomial
#'     \item \code{2}: Logistic Normal iid
#'     \item \code{3}: Logistic Normal 1dar1
#'     \item \code{4}: Logistic Normal 2d correlation (constant by sex, 1dar1 by age)
#'   }
#' @param ISS_SrvAgeComps Input sample sizes
#'   [n_regions × n_yrs × n_sexes × n_srv_fleets × n_sims]
#'   (default: `100`)
#' @param ln_SrvAge_theta Overdispersion parameters
#'   [n_regions × n_sexes × n_srv_fleets]
#'   (default: `log(1)`)
#' @param ln_SrvAge_theta_agg Overdispersion parameters for aggregated comps
#'   [n_srv_fleets]
#'   (default: `log(1)`)
#' @param SrvAge_corr_pars_agg Correlation parameters (agg.) for options 3–4
#'   [n_srv_fleets]
#'   (default: `0.01`)
#' @param SrvAge_corr_pars Correlation parameters
#'   [n_regions × n_sexes × n_srv_fleets x 2]
#'   (default: `0.01`)
#' @param SrvAgeComps_Type Array [n_yrs × n_srv_fleets]
#'   (default: `2` = joint by sex, split by region)
#'   \itemize{
#'     \item \code{0}: Aggregated
#'     \item \code{1}: Split by sex and region
#'     \item \code{2}: Joint by sex, split by region
#'     \item \code{999}: Not simulated
#'   }
#' @param comp_srvlen_like Vector [n_srv_fleets] specifying likelihood for simulating length comps
#'   (default: all `0` = multinomial)
#'   \itemize{
#'     \item \code{0}: Multinomial
#'     \item \code{1}: Dirichlet-Multinomial
#'     \item \code{2}: Logistic Normal iid
#'     \item \code{3}: Logistic Normal 1dar1
#'     \item \code{4}: Logistic Normal 2d correlation (constant by sex, 1dar1 by length)
#'   }
#' @param ISS_SrvLenComps Input sample sizes
#'   [n_regions × n_yrs × n_sexes × n_srv_fleets × n_sims]
#'   (default: `100`)
#' @param ln_SrvLen_theta Overdispersion parameters
#'   [n_regions × n_sexes × n_srv_fleets]
#'   (default: `log(1)`)
#' @param ln_SrvLen_theta_agg Overdispersion parameters for aggregated comps
#'   [n_srv_fleets]
#'   (default: `log(1)`)
#' @param SrvLen_corr_pars_agg Correlation parameters (agg.) for options 3–4
#'   [n_srv_fleets]
#'   (default: `0.01`)
#' @param SrvLen_corr_pars Correlation parameters
#'   [n_regions × n_sexes × n_srv_fleets x 2]
#'   (default: `0.01`)
#' @param SrvLenComps_Type Array [n_yrs × n_srv_fleets]
#'   (default: `2` = joint by sex, split by region)
#'   \itemize{
#'     \item \code{0}: Aggregated
#'     \item \code{1}: Split by sex and region
#'     \item \code{2}: Joint by sex, split by region
#'     \item \code{999}: Not simulated
#'   }
#'
#' @export Setup_Sim_Survey
#' @family Simulation Setup
Setup_Sim_Survey <- function(ObsSrvIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_srv_fleets)),
                             sim_list,
                             srv_sel_input,
                             srv_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_srv_fleets, sim_list$n_sims)),
                             t_srv = array(0, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)),
                             srv_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)),
                             comp_srvage_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvAge_theta_agg = rep(log(1), sim_list$n_srv_fleets),
                             SrvAge_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
                             SrvAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
                             comp_srvlen_like = rep(0, sim_list$n_srv_fleets),
                             ISS_SrvLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
                             ln_SrvLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets)),
                             ln_SrvLen_theta_agg = rep(log(1), sim_list$n_srv_fleets),
                             SrvLen_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
                             SrvLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes, sim_list$n_srv_fleets, 2)),
                             SrvLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets))
                             ) {

  # Validate dimensions of all input parameters
  check_sim_dimensions(srv_sel_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_ages = sim_list$n_ages, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "srv_sel_input")
  check_sim_dimensions(srv_q_input, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_srv_fleets = sim_list$n_srv_fleets, n_sims = sim_list$n_sims, what = "srv_q_input")
  check_sim_dimensions(ObsSrvIdx_SE, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "ObsSrvIdx_SE")
  check_sim_dimensions(srv_idx_type, n_regions = sim_list$n_regions, n_srv_fleets = sim_list$n_srv_fleets, what = "srv_idx_type")
  check_sim_dimensions(t_srv, n_regions = sim_list$n_regions, n_srv_fleets = sim_list$n_srv_fleets, what = "t_srv")

  # Validate survey age composition parameters
  check_sim_dimensions(comp_srvage_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvage_like")
  check_sim_dimensions(ISS_SrvAgeComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_SrvAgeComps")
  check_sim_dimensions(ln_SrvAge_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_theta")
  check_sim_dimensions(ln_SrvAge_theta_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvAge_theta_agg")
  check_sim_dimensions(SrvAge_corr_pars_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_corr_pars_agg")
  check_sim_dimensions(SrvAge_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "SrvAge_corr_pars")
  check_sim_dimensions(SrvAgeComps_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets,
                       what = "SrvAgeComps_Type")

  # Validate suvey length composition parameters
  check_sim_dimensions(comp_srvlen_like, n_srv_fleets = sim_list$n_srv_fleets, what = "comp_srvlen_like")
  check_sim_dimensions(ISS_SrvLenComps, n_regions = sim_list$n_regions, n_years = sim_list$n_yrs,
                       n_sexes = sim_list$n_sexes, n_srv_fleets = sim_list$n_srv_fleets,
                       n_sims = sim_list$n_sims, what = "ISS_SrvLenComps")
  check_sim_dimensions(ln_SrvLen_theta, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_theta")
  check_sim_dimensions(ln_SrvLen_theta_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "ln_SrvLen_theta_agg")
  check_sim_dimensions(SrvLen_corr_pars_agg, n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_corr_pars_agg")
  check_sim_dimensions(SrvLen_corr_pars, n_regions = sim_list$n_regions, n_sexes = sim_list$n_sexes,
                       n_srv_fleets = sim_list$n_srv_fleets, what = "SrvLen_corr_pars")
  check_sim_dimensions(SrvLenComps_Type, n_years = sim_list$n_yrs, n_srv_fleets = sim_list$n_srv_fleets,
                       what = "SrvLenComps_Type")

  # output into list
  sim_list$srv_sel <- srv_sel_input
  sim_list$srv_q <- srv_q_input
  sim_list$ObsSrvIdx_SE <- ObsSrvIdx_SE
  sim_list$t_srv <- t_srv
  sim_list$srv_idx_type <- srv_idx_type

  # Survey age compositions
  sim_list$comp_srvage_like <- comp_srvage_like
  sim_list$ISS_SrvAgeComps <- ISS_SrvAgeComps
  sim_list$ln_SrvAge_theta <- ln_SrvAge_theta
  sim_list$ln_SrvAge_theta_agg <- ln_SrvAge_theta_agg
  sim_list$SrvAge_corr_pars_agg <- SrvAge_corr_pars_agg
  sim_list$SrvAge_corr_pars <- SrvAge_corr_pars
  sim_list$SrvAgeComps_Type <- SrvAgeComps_Type

  # Survey length compositions
  sim_list$comp_srvlen_like <- comp_srvlen_like
  sim_list$ISS_SrvLenComps <- ISS_SrvLenComps
  sim_list$ln_SrvLen_theta <- ln_SrvLen_theta
  sim_list$ln_SrvLen_theta_agg <- ln_SrvLen_theta_agg
  sim_list$SrvLen_corr_pars_agg <- SrvLen_corr_pars_agg
  sim_list$SrvLen_corr_pars <- SrvLen_corr_pars
  sim_list$SrvLenComps_Type <- SrvLenComps_Type

  return(sim_list)

} # end function

#' Helper function to set up survey age overdispersion parameters
#'
#' @param input_list Input list
#' @keywords internal
do_SrvAge_theta_mapping <- function(input_list) {

  # setup counters
  counter_srvage_agg <- 1
  counter_srvage <- 1

  # initialize array to set up mapping
  map_SrvAge_theta <- input_list$par$ln_SrvAge_theta
  map_SrvAge_theta_agg <- input_list$par$ln_SrvAge_theta_agg
  map_SrvAge_theta[] <- NA
  map_SrvAge_theta_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # get unique survey comp types
    srvage_comp_type <- unique(input_list$data$SrvAgeComps_Type[,f])

    # If aggregated (ages)
    if(any(srvage_comp_type == 0) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
      map_SrvAge_theta_agg[f] <- counter_srvage_agg
      counter_srvage_agg <- counter_srvage_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(srvage_comp_type == 1) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
          map_SrvAge_theta[r,s,f] <- counter_srvage
          counter_srvage <- counter_srvage + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(srvage_comp_type == 2) && input_list$data$SrvAgeComps_LikeType[f] != 0 && s == 1) {
          map_SrvAge_theta[r,1,f] <- counter_srvage
          counter_srvage <- counter_srvage + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any age comps for a given fleet
    if(input_list$data$SrvAgeComps_LikeType[f] == 0 || sum(input_list$data$UseSrvAgeComps[,,f]) == 0) {
      map_SrvAge_theta[,,f] <- NA
      map_SrvAge_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_SrvAge_theta <- factor(map_SrvAge_theta)
  input_list$map$ln_SrvAge_theta_agg <- factor(map_SrvAge_theta_agg)

  return(input_list)
}

#' Helper function to set up survey length overdispersion parameters
#'
#' @param input_list Input list
#' @keywords internal
do_SrvLen_theta_mapping <- function(input_list) {

  # setup counters
  counter_srvlen_agg <- 1
  counter_srvlen <- 1

  # initialize array to set up mapping
  map_SrvLen_theta <- input_list$par$ln_SrvLen_theta
  map_SrvLen_theta_agg <- input_list$par$ln_SrvLen_theta_agg
  map_SrvLen_theta[] <- NA
  map_SrvLen_theta_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # get unique survey comp types
    srvlen_comp_type <- unique(input_list$data$SrvLenComps_Type[,f])

    # If aggregated (ages)
    if(any(srvlen_comp_type == 0) && input_list$data$SrvLenComps_LikeType[f] != 0) {
      map_SrvLen_theta_agg[f] <- counter_srvlen_agg
      counter_srvlen_agg <- counter_srvlen_agg + 1 # aggregated
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # if split by sex and region
        if(any(srvlen_comp_type == 1) && input_list$data$SrvLenComps_LikeType[f] != 0) {
          map_SrvLen_theta[r,s,f] <- counter_srvlen
          counter_srvlen <- counter_srvlen + 1 # split by sex and region
        }

        # joint by sex, split by region
        if(any(srvlen_comp_type == 2) && input_list$data$SrvLenComps_LikeType[f] != 0 && s == 1) {
          map_SrvLen_theta[r,1,f] <- counter_srvlen
          counter_srvlen <- counter_srvlen + 1 # joint by sex, split by region
        }

      } # end s loop
    } # end r loop

    # If we are using a multinomial or there aren't any lenght comps for a given fleet
    if(input_list$data$SrvLenComps_LikeType[f] == 0 || sum(input_list$data$UseSrvLenComps[,,f]) == 0) {
      map_SrvLen_theta[,,f] <- NA
      map_SrvLen_theta_agg[f] <- NA
    }

  } # end f loop

  # Input into mapping list
  input_list$map$ln_SrvLen_theta <- factor(map_SrvLen_theta)
  input_list$map$ln_SrvLen_theta_agg <- factor(map_SrvLen_theta_agg)

  return(input_list)
}

#' Helper function to set up survey age overdispersion correlation parameters
#'
#' @param input_list Input list
#' @keywords internal
do_SrvAge_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_srvage_corr <- 1
  counter_srvage_corr_agg <- 1

  # initialize array to set up mapping
  map_SrvAge_corr_pars <- input_list$par$SrvAge_corr_pars
  map_SrvAge_corr_pars_agg <- input_list$par$SrvAge_corr_pars_agg
  map_SrvAge_corr_pars[] <- NA
  map_SrvAge_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$SrvAgeComps_LikeType[f] == 0 || sum(input_list$data$UseSrvAgeComps[,,f]) == 0) {
      map_SrvAge_corr_pars[,,f,] <- NA
      map_SrvAge_corr_pars_agg[f] <- NA
      next # skip if none
    }

    # get unique survey comp types
    srvage_comp_type <- unique(input_list$data$SrvAgeComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(srvage_comp_type == 0) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
      if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
        map_SrvAge_corr_pars_agg[f] <- counter_srvage_corr_agg
        counter_srvage_corr_agg <- counter_srvage_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(srvage_comp_type == 1) && input_list$data$SrvAgeComps_LikeType[f] != 0) {
          if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
            map_SrvAge_corr_pars[r,s,f,1] <- counter_srvage_corr
            counter_srvage_corr <- counter_srvage_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(srvage_comp_type == 2) && input_list$data$SrvAgeComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$SrvAgeComps_LikeType[f] == 3) {
            map_SrvAge_corr_pars[r,1,f,1] <- counter_srvage_corr
            counter_srvage_corr <- counter_srvage_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$SrvAgeComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_SrvAge_corr_pars[r,1,f,i] <- counter_srvage_corr
              counter_srvage_corr <- counter_srvage_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$SrvAge_corr_pars_agg <- factor(map_SrvAge_corr_pars_agg)
  input_list$map$SrvAge_corr_pars <- factor(map_SrvAge_corr_pars)

  return(input_list)
}

#' Helper function to set up survey length overdispersion correlation parameters
#'
#' @param input_list Input list
#' @keywords internal
do_SrvLen_corr_pars_mapping <- function(input_list) {

  # setup counters
  counter_srvlen_corr <- 1
  counter_srvlen_corr_agg <- 1

  # initialize array to set up mapping
  map_SrvLen_corr_pars <- input_list$par$SrvLen_corr_pars
  map_SrvLen_corr_pars_agg <- input_list$par$SrvLen_corr_pars_agg
  map_SrvLen_corr_pars[] <- NA
  map_SrvLen_corr_pars_agg[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # No overdispersion parameters estimated
    if(input_list$data$SrvLenComps_LikeType[f] == 0 || sum(input_list$data$UseSrvLenComps[,,f]) == 0) {
      map_SrvLen_corr_pars[,,f,] <- NA
      map_SrvLen_corr_pars_agg[f] <- NA
      next # skip if none should be estimated
    }

    # get unique survey comp types
    srvlen_comp_type <- unique(input_list$data$SrvLenComps_Type[,f])

    # Aggregated Correlation Parameters
    if(any(srvlen_comp_type == 0) && input_list$data$SrvLenComps_LikeType[f] != 0) {
      if(input_list$data$SrvLenComps_LikeType[f] == 3) {
        map_SrvLen_corr_pars_agg[f] <- counter_srvlen_corr_agg
        counter_srvlen_corr_agg <- counter_srvlen_corr_agg + 1 # aggregated
      }
    }

    # Loop through to make sure mapping stuff off correctly
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {

        # Split by region and sex
        if(any(srvlen_comp_type == 1) && input_list$data$SrvLenComps_LikeType[f] != 0) {
          if(input_list$data$SrvLenComps_LikeType[f] == 3) {
            map_SrvLen_corr_pars[r,s,f,1] <- counter_srvlen_corr
            counter_srvlen_corr <- counter_srvlen_corr + 1
          }
        }

        # Joint by sex, split by region
        if(any(srvlen_comp_type == 2) && input_list$data$SrvLenComps_LikeType[f] != 0 && s == 1) {

          # 1dar1 correlation
          if(input_list$data$SrvLenComps_LikeType[f] == 3) {
            map_SrvLen_corr_pars[r,1,f,1] <- counter_srvlen_corr
            counter_srvlen_corr <- counter_srvlen_corr + 1
          }

          # 2dar1 correlation
          if(input_list$data$SrvLenComps_LikeType[f] == 4) {
            for(i in 1:2) {
              if(i == 2 && input_list$data$n_sexes == 1) next # skip if we only have 1 sex
              map_SrvLen_corr_pars[r,1,f,i] <- counter_srvlen_corr
              counter_srvlen_corr <- counter_srvlen_corr + 1
            } # end i
          } # end if

        }
      } # end s loop
    } # end r loop

  } # end f loop

  # Input into mapping list
  input_list$map$SrvLen_corr_pars_agg <- factor(map_SrvLen_corr_pars_agg)
  input_list$map$SrvLen_corr_pars <- factor(map_SrvLen_corr_pars)

  return(input_list)
}

#' Setup observed survey indices and composition data (age and length comps)
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param ObsSrvIdx Observed survey index data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_srv_fleets]}.
#'
#' @param ObsSrvIdx_SE Standard errors associated with \code{ObsSrvIdx},
#' also dimensioned \code{[n_regions, n_years, n_srv_fleets]}.
#'
#' @param UseSrvIdx Logical or binary indicator array (\code{[n_regions, n_years, n_srv_fleets]})
#' specifying whether to include a survey index in the likelihood (\code{1}) or ignore it (\code{0}).
#'
#' @param ObsSrvAgeComps Observed survey age composition data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_ages, n_sexes, n_srv_fleets]}. Values should reflect counts or proportions
#' (not required to sum to 1, but should be on a comparable scale).
#'
#' @param UseSrvAgeComps Indicator array (\code{[n_regions, n_years, n_srv_fleets]}) specifying whether
#' to fit survey age composition data (\code{1}) or ignore it (\code{0}).
#'
#' @param ObsSrvLenComps Observed survey length composition data as a numeric array with dimensions
#' \code{[n_regions, n_years, n_lens, n_sexes, n_srv_fleets]}. Values should reflect counts or proportions.
#'
#' @param UseSrvLenComps Indicator array (\code{[n_regions, n_years, n_srv_fleets]}) specifying whether
#' to fit survey length composition data (\code{1}) or ignore it (\code{0}).
#'
#' @param SrvAgeComps_LikeType Character vector of length \code{n_srv_fleets} specifying the likelihood
#' type used for survey age composition data. Options include \code{"Multinomial"}, \code{"Dirichlet-Multinomial"},
#' and \code{"iid-Logistic-Normal"}. Use \code{"none"} to omit the likelihood.
#'
#' @param SrvLenComps_LikeType Same as \code{SrvAgeComps_LikeType}, but for survey length composition data.
#'
#' @param SrvAgeComps_Type Character vector specifying how age compositions are structured by fleet and year range.
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
#' @param SrvLenComps_Type Same as \code{SrvAgeComps_Type}, but for length compositions.
#' @param srv_idx_type Character vector of length \code{n_srv_fleets} specifying the type of index data.
#' Options are \code{"abd"} for abundance, \code{"biom"} for biomass, and \code{"none"} if no index is available.
#'
#' @param ISS_SrvAgeComps Input sample size for age compositions, array dimensioned
#' \code{[n_regions, n_years, n_sexes, n_srv_fleets]}. Required if observed age comps are normalized
#' (i.e., sum to 1), to correctly scale the contribution to the likelihood.
#'
#' @param ISS_SrvLenComps Same as \code{ISS_SrvAgeComps}, but for length compositions.
#'
#' @param ... Additional arguments specifying starting values for overdispersion parameters
#' (e.g., \code{ln_SrvAge_theta}, \code{ln_SrvLen_theta}, \code{ln_SrvAge_theta_agg}, \code{ln_SrvLen_theta_agg}).
#'
#' @export Setup_Mod_SrvIdx_and_Comps
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_SrvIdx_and_Comps <- function(input_list,
                                       ObsSrvIdx,
                                       ObsSrvIdx_SE,
                                       UseSrvIdx,
                                       srv_idx_type,
                                       ObsSrvAgeComps,
                                       UseSrvAgeComps,
                                       ObsSrvLenComps,
                                       UseSrvLenComps,
                                       ISS_SrvAgeComps,
                                       ISS_SrvLenComps,
                                       SrvAgeComps_LikeType,
                                       SrvLenComps_LikeType,
                                       SrvAgeComps_Type,
                                       SrvLenComps_Type,
                                       ...
                                       ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)


  # Input Validation --------------------------------------------------------

  # Survey Indices
  check_data_dimensions(ObsSrvIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx')
  check_data_dimensions(ObsSrvIdx_SE, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvIdx_SE')
  check_data_dimensions(UseSrvIdx, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvIdx')
  check_data_dimensions(srv_idx_type, n_srv_fleets = input_list$data$n_srv_fleets, what = 'srv_idx_type')
  if(!all(srv_idx_type %in% c("biom", "abd", "none"))) stop("Invalid specification for srv_idx_type. Should be either abd, biom, or none")

  # Survey compositions
  check_data_dimensions(ObsSrvAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvAgeComps')
  check_data_dimensions(UseSrvAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvAgeComps')
  check_data_dimensions(UseSrvLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_srv_fleets = input_list$data$n_srv_fleets, what = 'UseSrvLenComps')
  if(input_list$data$fit_lengths == 1) check_data_dimensions(ObsSrvLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_lens = length(input_list$data$lens), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ObsSrvLenComps')
  if(!is.null(ISS_SrvAgeComps)) check_data_dimensions(ISS_SrvAgeComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvAgeComps')
  if(!is.null(ISS_SrvLenComps)) check_data_dimensions(ISS_SrvLenComps, n_regions = input_list$data$n_regions, n_years = length(input_list$data$years), n_sexes = input_list$data$n_sexes, n_srv_fleets = input_list$data$n_srv_fleets, what = 'ISS_SrvLenComps')
  check_data_dimensions(SrvAgeComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvAgeComps_LikeType')
  check_data_dimensions(SrvLenComps_LikeType, n_srv_fleets = input_list$data$n_srv_fleets, what = 'SrvLenComps_LikeType')
  if(!all(SrvAgeComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvAgeComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")
  if(!all(SrvLenComps_LikeType %in% c("none", "Multinomial", "Dirichlet-Multinomial", "iid-Logistic-Normal", "1d-Logistic-Normal", "2d-Logistic-Normal")))
    stop("Invalid specification for SrvLenComps_LikeType Should be either none, Multinomial, Dirichlet-Multinomial, iid-Logistic-Normal, 1d-Logistic-Normal, 2d-Logistic-Normal")


  # Survey Index Options ----------------------------------------------------

  srv_idx_type_vals <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets))
  for(f in 1:ncol(srv_idx_type_vals)) {
    if(srv_idx_type[f] == 'biom') srv_idx_type_vals[,f] <- 1 # biomass
    if(srv_idx_type[f] == 'abd') srv_idx_type_vals[,f] <- 0 # abundance
    if(srv_idx_type[f] == 'none') srv_idx_type_vals[,f] <- 999 # none
    collect_message(paste("Survey Index", "for survey fleet", f, "specified as:" , srv_idx_type[f]))
  } # end f loop


  # Survey Age Composition Options ------------------------------------------

  comp_srvage_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvAgeComps_LikeType[f] == 'none') comp_srvage_like_vals <- c(comp_srvage_like_vals, 999)
    if(SrvAgeComps_LikeType[f] == "Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 0)
    if(SrvAgeComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvage_like_vals <- c(comp_srvage_like_vals, 1)
    if(SrvAgeComps_LikeType[f] == "iid-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 2)
    if(SrvAgeComps_LikeType[f] == "1d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 3)
    if(SrvAgeComps_LikeType[f] == "2d-Logistic-Normal") comp_srvage_like_vals <- c(comp_srvage_like_vals, 4)
    collect_message(paste("Survey Age Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvAgeComps_LikeType[f]))
  } # end f loop

  # Specifying composition type
  SrvAgeComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvAgeComps_Type)) {

    # Extract out components from list
    tmp <- SrvAgeComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvAgeComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvAgeComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(comp_srvage_like_vals[fleet] == 4) stop("Age composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvAgeComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvAgeComps_Type_Mat))) stop("SrvAgeComps_Type_Mat is returning an NA. Did you update the year range of SrvAgeComps_Type_Mat?")

  # Survey Length Composition Options ---------------------------------------

  comp_srvlen_like_vals <- vector()
  for(f in 1:input_list$data$n_srv_fleets) {
    if(SrvLenComps_LikeType[f] == 'none') comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 999)
    if(SrvLenComps_LikeType[f] == "Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 0)
    if(SrvLenComps_LikeType[f] == "Dirichlet-Multinomial") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 1)
    if(SrvLenComps_LikeType[f] == "iid-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 2)
    if(SrvLenComps_LikeType[f] == "1d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 3)
    if(SrvLenComps_LikeType[f] == "2d-Logistic-Normal") comp_srvlen_like_vals <- c(comp_srvlen_like_vals, 4)
    collect_message(paste("Survey Length Composition Likelihoods", "for survey fleet", f, "specified as:" , SrvLenComps_LikeType[f]))
  } # end f loop

  SrvLenComps_Type_Mat <- array(NA, dim = c(length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(SrvLenComps_Type)) {

    # Extract out components from list
    tmp <- SrvLenComps_Type[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    comps_type_tmp <- tmp_vec[1] # get composition type
    fleet <- as.numeric(tmp_vec[5]) # extract fleet index

    # Checking character string
    if(!comps_type_tmp %in% c("agg", "spltRspltS", "spltRjntS", 'none')) stop("SrvLenComps_Type not specified correctly. Must be one of: agg, spltRspltS, spltRjntS, none")
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for SrvLenComps_Type This needs to be specified as CompType_Year_x-y_Fleet_x")

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
      if(comp_srvlen_like_vals[fleet] == 4) stop("Length composition likelihood specified as 2d-Logistic-Normal, but composition type is aggregated. This is not valid.")
      comps_type_val <- 0
    }
    if(comps_type_tmp == "spltRspltS") comps_type_val <- 1
    if(comps_type_tmp == "spltRjntS") comps_type_val <- 2
    if(comps_type_tmp == "none") comps_type_val <- 999

    # input into matrix
    SrvLenComps_Type_Mat[years,fleet] <- comps_type_val
  } # end i

  if(any(is.na(SrvLenComps_Type_Mat))) stop("SrvLenComps_Type_Mat is returning an NA. Did you update the year range of SrvLenComps_Type_Mat?")

  # ISS Munging -------------------------------------------------------------

  # Survey Ages
  if(is.null(ISS_SrvAgeComps)) {
    collect_message("No ISS is specified for SrvAgeComps. ISS weighting is calculated by summing up values from ObsSrvAgeComps each year")
    ISS_SrvAgeComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        # if aggregated across sexes and regions (0)
        if(SrvAgeComps_Type_Mat[y,f] == 0) ISS_SrvAgeComps[1,y,1,f] <- sum(srvObsSrvAgeComps[,y,,,f])
        # if split by region and sex
        if(SrvAgeComps_Type_Mat[y,f] == 1) ISS_SrvAgeComps[,y,,f] <- apply(srvObsSrvAgeComps[,y,,,f, drop = FALSE], c(1,4), sum)
        # if split by region, joint by sex
        if(SrvAgeComps_Type_Mat[y,f] == 2) ISS_SrvAgeComps[,y,1,f] <- apply(srvObsSrvAgeComps[,y,,,f, drop = FALSE], 1, sum)
      } # end f loop
    } # end y loop
  }

  # Survey Lengths
  if(is.null(ISS_SrvLenComps)) {
    collect_message("No ISS is specified for SrvLenComps. ISS weighting is calculated by summing up values from ObsSrvLenComps each year")
    ISS_SrvLenComps <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    for(y in 1:length(input_list$data$years)) {
      for(f in 1:input_list$data$n_srv_fleets) {
        # if aggregated across sexes and regions (0)
        if(SrvLenComps_Type_Mat[y,f] == 0) ISS_SrvLenComps[1,y,1,f] <- sum(ObsSrvLenComps[,y,,,f])
        # if split by region and sex
        if(SrvLenComps_Type_Mat[y,f] == 1) ISS_SrvLenComps[,y,,f] <- apply(ObsSrvLenComps[,y,,,f, drop = FALSE], c(1,4), sum)
        # if split by region, joint by sex
        if(SrvLenComps_Type_Mat[y,f] == 2) ISS_SrvLenComps[,y,1,f] <- apply(ObsSrvLenComps[,y,,,f, drop = FALSE], 1, sum)
      } # end f loop
    } # end y loop
  }

  # Populate Data List ------------------------------------------------------

  input_list$data$ISS_SrvAgeComps <- ISS_SrvAgeComps
  input_list$data$ISS_SrvLenComps <- ISS_SrvLenComps
  input_list$data$ObsSrvIdx <- ObsSrvIdx
  input_list$data$ObsSrvIdx_SE <- ObsSrvIdx_SE
  input_list$data$UseSrvIdx <- UseSrvIdx
  input_list$data$ObsSrvAgeComps <- ObsSrvAgeComps
  input_list$data$UseSrvAgeComps <- UseSrvAgeComps
  input_list$data$ObsSrvLenComps <- ObsSrvLenComps
  input_list$data$UseSrvLenComps <- UseSrvLenComps
  input_list$data$SrvAgeComps_LikeType <- comp_srvage_like_vals
  input_list$data$SrvLenComps_LikeType <- comp_srvlen_like_vals
  input_list$data$SrvAgeComps_Type <- SrvAgeComps_Type_Mat
  input_list$data$SrvLenComps_Type <- SrvLenComps_Type_Mat
  input_list$data$srv_idx_type <- srv_idx_type_vals

  # Populate Parameter List -------------------------------------------------

  # Dispersion parameters for the survey age comps
  if("ln_SrvAge_theta" %in% names(starting_values)) input_list$par$ln_SrvAge_theta <- starting_values$ln_SrvAge_theta
  else input_list$par$ln_SrvAge_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey age comps
  if("SrvAge_corr_pars" %in% names(starting_values)) input_list$par$SrvAge_corr_pars <- starting_values$SrvAge_corr_pars
  else input_list$par$SrvAge_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated
  if("ln_SrvAge_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvAge_theta_agg <- starting_values$ln_SrvAge_theta_agg
  else input_list$par$ln_SrvAge_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  # aggregated correlation parameters
  if("SrvAge_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvAge_corr_pars_agg <- starting_values$SrvAge_corr_pars_agg
  else input_list$par$SrvAge_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Dispersion parameters for survey length comps
  if("ln_SrvLen_theta" %in% names(starting_values)) input_list$par$ln_SrvLen_theta <- starting_values$ln_SrvLen_theta
  else input_list$par$ln_SrvLen_theta <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # logistic normal correlation parameters for survey length comps
  if("SrvLen_corr_pars" %in% names(starting_values)) input_list$par$SrvLen_corr_pars <- starting_values$SrvLen_corr_pars
  else input_list$par$SrvLen_corr_pars <- array(0.01, dim = c(input_list$data$n_regions, input_list$data$n_sexes, input_list$data$n_srv_fleets, 2))

  # aggregated
  if("ln_SrvLen_theta_agg" %in% names(starting_values)) input_list$par$ln_SrvLen_theta_agg <- starting_values$ln_SrvLen_theta_agg
  else input_list$par$ln_SrvLen_theta_agg <- array(0, dim = c(input_list$data$n_srv_fleets))

  if("SrvLen_corr_pars_agg" %in% names(starting_values)) input_list$par$SrvLen_corr_pars_agg <- starting_values$SrvLen_corr_pars_agg
  else input_list$par$SrvLen_corr_pars_agg <- array(0.01, dim = c(input_list$data$n_srv_fleets))

  # Mapping Options ---------------------------------------------------------

  input_list <- do_SrvAge_theta_mapping(input_list)
  input_list <- do_SrvLen_theta_mapping(input_list)
  input_list <- do_SrvAge_corr_pars_mapping(input_list)
  input_list <- do_SrvLen_corr_pars_mapping(input_list)

  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}


#' Helper function to set up survey selctivity fixed effects mapping
#'
#' @param input_list Input list
#' @param srv_fixed_sel_pars_spec Character vector specifying survey selectivity fixed effects parameterization
#' @keywords internal
#' Helper function to set up survey selctivity fixed effects mapping
#'
#' @param input_list Input list
#' @param srv_fixed_sel_pars_spec Character vector specifying survey selectivity fixed effects parameterization
#' @keywords internal
do_srv_fixed_sel_pars_mapping <- function(input_list, srv_fixed_sel_pars_spec) {

  # Initialize counter and mapping array for fixed effects survey selectivity
  srv_fixed_sel_pars_counter <- 1
  map_srv_fixed_sel_pars <- input_list$par$ln_srv_fixed_sel_pars
  map_srv_fixed_sel_pars[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # Validate Options
    if(!srv_fixed_sel_pars_spec[f] %in% c("est_all", "est_shared_r", "est_shared_r_s", "fix", "est_shared_s") &&
       !stringr::str_detect(srv_fixed_sel_pars_spec[f], "est_shared_f_\\d+"))
      stop("srv_fixed_sel_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")

    # Skip fleet sharing specs in first pass
    if(stringr::str_detect(srv_fixed_sel_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # Only add a counter if caatches are avaliable in some years for a given region and fleet combination
      if(sum(input_list$data$UseSrvIdx[r,,f]) > 0) {

        # Extract number of srvery selectivity blocks
        srvsel_blocks_tmp <- unique(as.vector(input_list$data$srv_sel_blocks[r,,f]))

        for(s in 1:input_list$data$n_sexes) {
          for(b in 1:length(srvsel_blocks_tmp)) {

            block_years <- which(input_list$data$srv_sel_blocks[r,,f] == srvsel_blocks_tmp[b]) # figure out block years
            sel_model_this_block <- unique(input_list$data$srv_sel_model[r, block_years, f]) # get selectivity form for a given block
            if(length(sel_model_this_block) > 1) stop("Block ", srvsel_blocks_tmp[b], " for fleet ", f, " region ", r, " has multiple selectivity models assigned to it")

            # determine maximum selectivity parameters
            if(sel_model_this_block == 2) max_sel_pars <- 1 # exponential
            if(sel_model_this_block %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(sel_model_this_block == 4) max_sel_pars <- 6 # double normal

            for(i in 1:max_sel_pars) {

              # Estimate all selectivity fixed effects parameters within the constraints of the defined blocks
              if(srv_fixed_sel_pars_spec[f] == "est_all") {
                map_srv_fixed_sel_pars[r,i,b,s,f] <- srv_fixed_sel_pars_counter
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srv_fixed_sel_pars_spec[f] == 'est_shared_r' && r == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  if(srvsel_blocks_tmp[b] %in% input_list$data$srv_sel_blocks[rr,,f]) {
                    map_srv_fixed_sel_pars[rr, i, b, s, f] <- srv_fixed_sel_pars_counter
                  } # end if
                } # end rr loop
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srv_fixed_sel_pars_spec[f] == 'est_shared_s' && s == 1) {
                for(ss in 1:input_list$data$n_sexes) {
                  map_srv_fixed_sel_pars[r, i, b, ss, f] <- srv_fixed_sel_pars_counter
                } # end ss loop
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srv_fixed_sel_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                for(rr in 1:input_list$data$n_regions) {
                  for(ss in 1:input_list$data$n_sexes) {
                    if(srvsel_blocks_tmp[b] %in% input_list$data$srv_sel_blocks[rr,,f]) {
                      map_srv_fixed_sel_pars[rr, i, b, ss, f] <- srv_fixed_sel_pars_counter
                    } # end if
                  } # end ss loop
                } #end rr loop
                srv_fixed_sel_pars_counter <- srv_fixed_sel_pars_counter + 1
              } # end if

            } # end i loop
          } # end b loop
        } # end s loop
      } # end if statement
    } # end r loop

    # fix all parameters
    if(srv_fixed_sel_pars_spec[f] == "fix") map_srv_fixed_sel_pars[,,,,f] <- NA
    collect_message("srv_fixed_sel_pars_spec is specified as: ", srv_fixed_sel_pars_spec[f], " for survey fleet ", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_srv_fleets) {
    if(stringr::str_detect(srv_fixed_sel_pars_spec[f], "est_shared_f")) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(srv_fixed_sel_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_srv_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(srv_fixed_sel_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_srv_fixed_sel_pars[,,,,f] <- map_srv_fixed_sel_pars[,,,,flt_shared]
      collect_message("srv_fixed_sel_pars_spec is specified as: ", srv_fixed_sel_pars_spec[f], " for survey fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_srv_fixed_sel_pars <- factor(map_srv_fixed_sel_pars)

  return(input_list)
}


#' Helper function to set up survey catchability mapping
#'
#' @param input_list Input list
#' @param srv_q_spec Character vector specifying survey catchability parameterization
#' @keywords internal
do_srv_q_mapping <- function(input_list, srv_q_spec) {

  # Initialize counter and mapping array for survey catchability
  srv_q_counter <- 1
  map_srv_q <- input_list$par$ln_srv_q
  map_srv_q[] <- NA

  for(f in 1:input_list$data$n_srv_fleets) {

    # Validate options
    if(!is.null(srv_q_spec)) {
      if(!srv_q_spec[f] %in% c("est_all", "est_shared_r", "fix"))
        stop("srv_q_spec not correctly specfied. Should be one of these: est_all, est_shared_r, fix")
    }

    for(r in 1:input_list$data$n_regions) {

      if(sum(input_list$data$UseSrvIdx[r,,f]) == 0) {
        map_srv_q[r,,f] <- NA # fix parameters if we are not using survey indices for these fleets and regions
      } else {

        # Extract number of survey catchability blocks
        srvq_blocks_tmp <- unique(as.vector(input_list$data$srv_q_blocks[r,,f]))

        for(b in 1:length(srvq_blocks_tmp)) {

          # Estimate for all regions
          if(srv_q_spec[f] == 'est_all') {
            map_srv_q[r,b,f] <- srv_q_counter
            srv_q_counter <- srv_q_counter + 1
          } # end if

          # Estimate but share q across regions
          if(srv_q_spec[f] == 'est_shared_r' && r == 1) {
            for(rr in 1:input_list$data$n_regions) {
              if(srvq_blocks_tmp[b] %in% input_list$data$srv_q_blocks[rr,,f]) {
                map_srv_q[rr, b, f] <- srv_q_counter
              } # end if
            } # end rr loop
            srv_q_counter <- srv_q_counter + 1
          } # end if

        } # end b loop
      } # end else loop
    } # end r loop

    # fix all parameters
    if(srv_q_spec[f] == 'fix') map_srv_q[,,f] <- NA
    collect_message("srv_q_spec is specified as: ", srv_q_spec[f], " for survey fleet ", f)
  } # end f loop

  # input into mapping list
  input_list$map$ln_srv_q <- factor(map_srv_q)

  return(input_list)
}

#' Helper function to set up survey process error mapping
#'
#' @param input_list Input list
#' @param corr_opt_semipar Character vector specifying correlation parameter options
#' @param srvsel_pe_pars_spec Character vector specifying survey process error parameterization
#'
#' @keywords internal
do_srvsel_pe_pars_mapping <- function(input_list, srvsel_pe_pars_spec, corr_opt_semipar) {

  # Initialize counter and mapping array for survey process errors
  srvsel_pe_pars_counter <- 1 # initalize counter
  map_srvsel_pe_pars <- input_list$par$srvsel_pe_pars # initalize array
  map_srvsel_pe_pars[] <- NA

  # Survey process error parameters
  for(f in 1:input_list$data$n_srv_fleets) {

    # Validate options
    if(!is.null(srvsel_pe_pars_spec)) {
      if(!srvsel_pe_pars_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s") &&
         !stringr::str_detect(srvsel_pe_pars_spec[f], "est_shared_f_\\d+"))
        stop("srvsel_pe_pars_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, fix, or est_shared_f_# (where # is fleet number)")
    }

    # Skip fleet sharing specs in first pass
    if(!is.null(srvsel_pe_pars_spec)) if(stringr::str_detect(srvsel_pe_pars_spec[f], "est_shared_f")) next

    for(r in 1:input_list$data$n_regions) {

      # if no time-variation, then fix all parameters for this fleet
      if(input_list$data$cont_tv_srv_sel[r,f] == 0 || sum(input_list$data$UseSrvIdx[r,,f]) == 0) {
        map_srvsel_pe_pars[r,,,f] <- NA
      } else { # if we have time-variation

        # Figure out max number of selectivity parameters for a given region and fleet
        if(unique(input_list$data$srv_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
        if(unique(input_list$data$srv_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
        if(unique(input_list$data$srv_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal

        for(s in 1:input_list$data$n_sexes) {

          # If iid time-variation or random walk for this fleet
          if(input_list$data$cont_tv_srv_sel[r,f] %in% c(1,2)) {

            for(i in 1:max_sel_pars) {

              # either fixing parameters or not used for a given fleet
              if(srvsel_pe_pars_spec[f] %in% c("none", "fix")) map_srvsel_pe_pars[r,i,s,f] <- NA

              # Estimating all parameters separately (unique for each region, sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == "est_all") {
                map_srvsel_pe_pars[r,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_srvsel_pe_pars[,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_srvsel_pe_pars[r,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_srvsel_pe_pars[,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

            } # end i loop
          } # end iid or random walk variation

          # If 3d gmrf or 2dar1
          if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) {

            # Set up indexing to loop through
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4)) idx = 1:4 # 3dgmrf (1 = pcorr_age, 2 = pcorr_year, 3= pcorr_cohort, 4 = log_sigma)
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(5)) idx = c(1,2,4) # 2dar1 (1 = pcorr_bin, 2 = pcorr_year, 4 = log_sigma)
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4) && input_list$data$Selex_Type == 1) stop("Cohort-based selectivity deviations are specified, but selectivity is specified as length-based. Please choose another deviation form!")

            for(i in idx) {

              # either fixing parameters or not used for a given fleet
              if(srvsel_pe_pars_spec[f] %in% c("none", "fix")) map_srvsel_pe_pars[r,i,s,f] <- NA

              # Estimating all process error parameters
              if(srvsel_pe_pars_spec[f] == "est_all") {
                map_srvsel_pe_pars[r,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              } # end est_all

              # Estimating process error parameters shared across regions (but unique for each sex, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r' && r == 1) {
                map_srvsel_pe_pars[,i,s,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across sexes (but unique for each region, fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_s' && s == 1) {
                map_srvsel_pe_pars[r,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

              # Estimating process error parameters shared across regions and sexes (but unique for each fleet, parameter)
              if(srvsel_pe_pars_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                map_srvsel_pe_pars[,i,,f] <- srvsel_pe_pars_counter
                srvsel_pe_pars_counter <- srvsel_pe_pars_counter + 1
              }

            } # end i loop

            # Options to set correaltions to 0 for 3dgmrf
            if(!is.null(corr_opt_semipar)) {

              opt <- input_list$data$cont_tv_srv_sel[r,f] # get random effects options

              # Validate options
              if(!corr_opt_semipar[f] %in% c(NA, "corr_zero_y", "corr_zero_b", "corr_zero_y_b", "corr_zero_c", "corr_zero_y_c", "corr_zero_b_c", "corr_zero_y_b_c"))
                stop("corr_opt_semipar not correctly specfied. Should be one of these: corr_zero_y, corr_zero_b, corr_zero_y_b, corr_zero_c, corr_zero_y_c, corr_zero_b_c, corr_zero_y_b_c, NA")
              if(opt == 5 && corr_opt_semipar[f] %in% c("corr_zero_c","corr_zero_y_c","corr_zero_b_c","corr_zero_y_b_c"))
                stop("Invalid corr_opt_semipar for 2dar1 (opt=5): cohort correlations are not allowed.")

              if (opt %in% c(3,4,5)) {
                # 2d and 3d options
                if (corr_opt_semipar[f] == "corr_zero_y")    map_srvsel_pe_pars[,2,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_b")    map_srvsel_pe_pars[,1,,f]     <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b")  map_srvsel_pe_pars[,1:2,,f]   <- NA
              }

              if(opt %in% c(3,4)) {
                # 3d gmrf options only (adds the cohort dimension)
                if (corr_opt_semipar[f] == "corr_zero_c")      map_srvsel_pe_pars[,3,,f]   <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_c")    map_srvsel_pe_pars[,2:3,,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_b_c")    map_srvsel_pe_pars[,c(1,3),,f] <- NA
                if (corr_opt_semipar[f] == "corr_zero_y_b_c")  map_srvsel_pe_pars[,1:3,,f] <- NA
              }

              # Reset numbering for mapping off correlation parameters for clarity
              non_na_positions <- which(!is.na(map_srvsel_pe_pars))
              map_srvsel_pe_pars[non_na_positions] <- seq_along(non_na_positions)
              collect_message("corr_opt_semipar is specified as: ", corr_opt_semipar[f], "for survey fleet", f)

            }
          } # end if 3d gmrf marginal or conditional variance

          # fix all parameters
          if(srvsel_pe_pars_spec[f] == "fix") map_srvsel_pe_pars[r,,s,f] <- NA

        } # end s loop
      } # end else
    } # end r loop

    if(!is.null(srvsel_pe_pars_spec)) collect_message("srvsel_pe_pars_spec is specified as: ", srvsel_pe_pars_spec[f], "for survey fleet", f)

  } # end f loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_srv_fleets) {
    if(stringr::str_detect(srvsel_pe_pars_spec[f], "est_shared_f") && !is.null(srvsel_pe_pars_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(srvsel_pe_pars_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_srv_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(srvsel_pe_pars_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_srvsel_pe_pars[,,,f] <- map_srvsel_pe_pars[,,,flt_shared]
      collect_message("srvsel_pe_pars_spec is specified as: ", srvsel_pe_pars_spec[f], " for survey fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$srvsel_pe_pars <- factor(map_srvsel_pe_pars)

  return(input_list)
}

#' Helper function to set up survey selectivity deviations mapping
#'
#' @param input_list Input list
#' @param srvsel_devs_shared_ages List object for specifying which ages are shared when selectivity deviations are semi-parametric (e.g., list(1:5, 6:10, 11:30) specifies that ages 1-5, 6-10, and 11-30 have the same deviations.)
#' @param srv_sel_devs_spec Character vector specifying survey selectivity deviations parameterization
#'
#' @keywords internal
do_srvsel_devs_mapping <- function(input_list, srv_sel_devs_spec, srvsel_devs_shared_ages) {

  # Initialize counter and mapping array for survey selectivity deviations
  srvsel_devs_counter <- 1
  map_srvsel_devs <- input_list$par$ln_srvsel_devs
  map_srvsel_devs[] <- NA

  for(r in 1:input_list$data$n_regions) {
    for(f in 1:input_list$data$n_srv_fleets) {

      # Validate options
      if(!is.null(srv_sel_devs_spec)) {
        if(!srv_sel_devs_spec[f] %in% c("fix", "none", "est_all", "est_shared_r", "est_shared_s", "est_shared_r_s", "est_shared_a", "est_shared_r_a", "est_shared_r_a_s", "est_shared_a_s") &&
           !stringr::str_detect(srv_sel_devs_spec[f], "est_shared_f_\\d+"))
          stop("srv_sel_devs_spec not correctly specfied. Should be one of these: est_all, est_shared_r, est_shared_r_s, est_shared_s, est_shared_a, est_shared_r_a, est_shared_r_a_s, est_shared_r_s, fix, or est_shared_f_# (where # is fleet number)")
        if(srv_sel_devs_spec[f] %in% c("est_shared_a", "est_shared_r_a", "est_shared_r_a_s", "est_shared_a_s") &&
           !input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) stop("Sharing age deviations with iid or random walk parametric forms is not supported!")
      }

      # Skip fleet sharing specs in first pass
      if(!is.null(srv_sel_devs_spec)) if(stringr::str_detect(srv_sel_devs_spec[f], "est_shared_f")) next

      for(s in 1:input_list$data$n_sexes) {
        for(y in 1:(length(input_list$data$years) + input_list$data$n_proj_yrs_devs)) {

          # if no time-variation, then fix all parameters for this fleet
          if(input_list$data$cont_tv_srv_sel[r,f] == 0 || sum(input_list$data$UseSrvIdx[r,,f]) == 0) {
            map_srvsel_devs[r,y,,s,f] <- NA
          } else {

            # Figure out max number of selectivity parameters for a given region and fleet
            if(unique(input_list$data$srv_sel_model[r,,f]) %in% 2) max_sel_pars <- 1 # exponential
            if(unique(input_list$data$srv_sel_model[r,,f]) %in% c(0,1,3)) max_sel_pars <- 2 # logistic or gamma
            if(unique(input_list$data$srv_sel_model[r,,f]) == 4) max_sel_pars <- 6 # double normal

            # If iid or random walk time-variation for this fleet
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(1,2)) {

              for(i in 1:max_sel_pars) {
                # Estimating all selectivity deviations across regions, sexes, fleets, and parameter
                if(srv_sel_devs_spec[f] == 'est_all') {
                  map_srvsel_devs[r,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating selectivity deviations across sexes, fleets, and parameters, but shared across regions
                if(srv_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_srvsel_devs[,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating selectivity deviations across regions, fleets, and parameters, but shared across sexes
                if(srv_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_srvsel_devs[r,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating selectivity deviations across fleets, and parameters, but shared across sexes and regions
                if(srv_sel_devs_spec[f] == 'est_shared_r_s' && r == 1 && s == 1) {
                  map_srvsel_devs[,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

              } # end i loop
            } # end iid or random walk variation

            # If 3d gmrf for this fleet
            if(input_list$data$cont_tv_srv_sel[r,f] %in% c(3,4,5)) {

              for(i in 1:length(input_list$data$ages)) {
                # Estimating all selectivity deviations across regions, years and bins
                if(srv_sel_devs_spec[f] == 'est_all') {
                  map_srvsel_devs[r,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across regions
                if(srv_sel_devs_spec[f] == 'est_shared_r' && r == 1) {
                  map_srvsel_devs[,y,i,s,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes
                if(srv_sel_devs_spec[f] == 'est_shared_s' && s == 1) {
                  map_srvsel_devs[r,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                # Estimating all selectivity deviations across years and bins, but shared across sexes and regions
                if(srv_sel_devs_spec[f] == 'est_shared_r_s' && s == 1 && r == 1) {
                  map_srvsel_devs[,y,i,,f] <- srvsel_devs_counter
                  srvsel_devs_counter <- srvsel_devs_counter + 1
                }

                if(srv_sel_devs_spec[f] == 'est_shared_a') {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[r,y,srvsel_devs_shared_ages[[k]],s,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

                if(srv_sel_devs_spec[f] == 'est_shared_r_a' && r == 1) {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[,y,srvsel_devs_shared_ages[[k]],s,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

                if(srv_sel_devs_spec[f] == 'est_shared_a_s' && s == 1) {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[r,y,srvsel_devs_shared_ages[[k]],,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

                if(srv_sel_devs_spec[f] == 'est_shared_r_a_s' && s == 1 && r == 1) {
                  for(k in 1:length(srvsel_devs_shared_ages)) {
                    map_srvsel_devs[,y,srvsel_devs_shared_ages[[k]],,f] <- srvsel_devs_counter
                    srvsel_devs_counter <- srvsel_devs_counter + 1
                  } # end k loop
                }

              } # end i loop
            } # end 3d gmrf

          } # end else
        } # end y loop
      } # end s loop

      if(!is.null(srv_sel_devs_spec)) collect_message("srv_sel_devs_spec is specified as: ", srv_sel_devs_spec[f], "for survey fleet", f, "and region ", r)

    } # end f loop
  } # end r loop

  # Handle fleet sharing after all base mappings are established
  for(f in 1:input_list$data$n_srv_fleets) {
    if(stringr::str_detect(srv_sel_devs_spec[f], "est_shared_f") && !is.null(srv_sel_devs_spec)) {
      # extract fleet sharing index
      flt_shared <- as.numeric(unlist(stringr::str_extract_all(srv_sel_devs_spec[f], "\\d+")))

      # Validate options here
      if(flt_shared > input_list$data$n_srv_fleets || flt_shared < 1) stop("Fleet sharing specification 'est_shared_f", flt_shared, "' for fleet ", f, " references invalid fleet number.")
      if(stringr::str_detect(srv_sel_devs_spec[flt_shared], "est_shared_f")) stop("Fleet ", f, " cannot share with fleet ", flt_shared, " because fleet ", flt_shared, " is self-sharing parameters, which does not make sense.")

      # Copy mapping from reference fleet
      map_srvsel_devs[,,,,f] <- map_srvsel_devs[,,,,flt_shared]
      collect_message("srv_sel_devs_spec is specified as: ", srv_sel_devs_spec[f], " for survey fleet ", f, " (sharing with fleet ", flt_shared, ")")
    } # end if statement
  } # end f loop

  # input into mapping list
  input_list$map$ln_srvsel_devs <- factor(map_srvsel_devs)
  input_list$data$map_ln_srvsel_devs <- array(as.numeric(input_list$map$ln_srvsel_devs), dim = dim(input_list$par$ln_srvsel_devs))

  return(input_list)
}


#' Setup survey selectivity and catchability specifications
#'
#' @param input_list List containing a data list, parameter list, and map list
#' @param cont_tv_srv_sel Character vector specifying the form of continuous time-varying selectivity for each survey fleet.
#' The vector must be length \code{n_srv_fleets}, and each element must follow the structure:
#' \code{"<time variation type>_Fleet_<fleet number>"}.
#'
#' Valid time variation types include:
#' \itemize{
#'   \item \code{"none"}: No continuous time variation. (default)
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
#' @param srv_sel_blocks Character vector specifying the survey selectivity blocks for each region and fleet.
#'
#' Each element must follow one of the following structures:
#' \itemize{
#'   \item `"Block_<block number>_Year_<start>-<end>_Fleet_<fleet number>"`
#'   \item `"Block_<block number>_Year_<start>-terminal_Fleet_<fleet number>"`
#'   \item `"none_Fleet_<fleet number>"`
#' }
#'
#' This argument defines how survey selectivity varies over time for each fleet:
#' \itemize{
#'   \item \code{"Block_..."} entries specify discrete time blocks during which selectivity parameters are assumed constant.
#'   \item \code{"none_..."} entries indicate that selectivity is constant across all years for the specified fleet.
#' }
#'
#' If time-block-based selectivity is specified for a fleet (via \code{srv_sel_blocks}), its corresponding continuous selectivity option (in \code{cont_tv_srv_sel}) must be set to \code{"none_Fleet_<fleet number>"}. The two approaches—blocked and continuous time-varying selectivity—are mutually exclusive.
#' The default for each fleet is \code{"none_Fleet_x"} (i.e., no selectivity blocks).
#'
#' @param srv_sel_model Character vector specifying the survey selectivity functional form for each fleet, and optionally by time block.
#'
#' Each element must follow one of the following structures:
#' \itemize{
#'   \item \code{"<selectivity model>_Fleet_<fleet number>"}
#'   \item \code{"<selectivity model>_Fleet_<fleet number>_Block_<block number>"}
#' }
#'
#' The first form applies a single selectivity model across all years for the specified fleet.
#' The second form allows the user to assign a distinct selectivity model to a specific time block, as defined in \code{srv_sel_blocks}.
#'
#' Available selectivity model types include:
#' \itemize{
#'   \item \code{"logist1"} — Logistic function with parameters \code{a50} and \code{k}.
#'   \item \code{"logist2"} — Logistic function with parameters \code{a50} and \code{a95}.
#'   \item \code{"gamma"} — Dome-shaped gamma function with parameters \code{amax} and \code{delta}.
#'   \item \code{"exponential"} — Exponential function with a power parameter.
#'   \item \code{"dbnrml"} — Double-normal function with six parameters.
#' }
#' If multiple selectivity time blocks are specified for a fleet (using \code{srv_sel_blocks}), then the corresponding selectivity model for each block must be explicitly defined using the \code{"<model>_Block_<block>_Fleet_<fleet>"} format.
#' If blocks are not defined for a fleet, use the \code{"<model>_Fleet_<fleet number>"} format only.
#' For mathematical definitions and implementation details of each selectivity form, refer to the model equations vignette.
#'
#' @param srv_q_blocks Character vector specifying survey catchability (q) blocks for each fleet.
#' Each element must follow the structure: \code{"Block_<block number>_Year_<start>-<end>_Fleet_<fleet number>"}
#' or \code{"none_Fleet_<fleet number>"}.
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
#' Internally, these specifications are converted to a \code{[n_regions, n_years, n_srv_fleets]} array,
#' where each block is mapped to the appropriate years and fleets.
#' @param srvsel_pe_pars_spec Character string specifying how process error parameters for survey selectivity
#' are estimated across regions and sexes. This is only relevant if \code{cont_tv_srv_sel} is not set to \code{"none"};
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
#' @param srv_fixed_sel_pars_spec Character string specifying the structure for estimating
#' fixed-effect parameters of the survey selectivity model (e.g., a50, k, amax).
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
#' @param srv_q_spec Character string specifying the structure of survey catchability (\code{q}) estimation
#' across regions. This controls whether separate or shared parameters are used.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates separate catchability parameters for each region.
#'   \item \code{"est_shared_r"}: Estimates a single catchability parameter shared across all regions.
#' }
#' @param srv_sel_devs_spec Character string specifying the structure of process error deviations
#' in time-varying survey selectivity dimensioned by the number of survey fleets. This determines how deviations are estimated across regions and sexes.
#'
#' Available options include:
#' \itemize{
#'   \item \code{"est_all"}: Estimates a separate deviation time series for each region and sex.
#'   \item \code{"est_shared_r"}: Shares deviations across regions (sex-specific deviations are still estimated).
#'   \item \code{"est_shared_s"}: Shares deviations across sexes (region-specific deviations are still estimated).
#'   \item \code{"est_shared_r_s"}: Shares deviations across both regions and sexes, estimating a single deviation time series.
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
#' This argument is only used when a continuous time-varying selectivity form is specified (e.g., via \code{cont_tv_srv_sel}).
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
#'     \item \code{"corr_zero_b"}: Sets bin correlations to 0.
#'     \item \code{"corr_zero_y_b"}: Sets both year and bin correlations to 0.
#'     \item \code{"corr_zero_c"}: Sets cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_y_c"}: Sets year and cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_b_c"}: Sets bin (age) and cohort correlations to 0. Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}.
#'     \item \code{"corr_zero_y_b_c"}: Sets all correlations (year, bin (age), and cohort) to 0.
#'       Only valid for \code{cont_tv_sel} = \code{"3dmarg"} or \code{"3dcond"}; equivalent to an iid structure.
#'   }
#'
#' These correlation-suppression flags are ignored when \code{cont_tv_sel} is set to any other value.
#' @param Use_srv_q_prior Integer (0 or 1). Flag to enable/disable survey catchability priors.
#'   When set to 1, applies log-normal priors to survey selectivity parameters as specified
#'   in \code{srv_q_prior}. When set to 0, no priors are applied.
#' @param srv_q_prior Data frame containing prior specifications for survey catchability parameters.
#'   Must include columns: \code{region} (region index), \code{fleet} (fleet index),
#'   \code{block} (time block index), \code{mu} (prior mean on natural scale), and \code{sd} (prior standard deviation on log scale).
#'   Each row specifies a log-normal prior N(log(mu), sd) for a given catchability parameter.
#'   Only parameters with rows in this data frame will have priors applied.
#' @param ... Additional arguments specifying starting values for survey selectivity and catchability parameters (srvsel_pe_pars, ln_srvsel_devs, ln_srv_fixed_sel_pars, ln_srv_q, srv_q_coeff)
#' @param srv_q_formula A named list of formulas specifying environmental covariate relationships
#'   for each region and survey fleet. Each element should be named using the convention
#'   `"Region_<region>_Fleet_<fleet>"` and contain a formula object using covariate names present in
#'   `srv_q_cov_dat`. The formula determines how environmental covariates influence survey catchability.
#'   If `NULL`, no environmental covariate effects are included.
#' @param srv_q_cov_dat A named list containing time series vectors (typically by year) of environmental
#'   covariates used in the `srv_q_formula`. Each entry should be a numeric vector of length equal to the
#'   number of years, and names must match the variable names used in the formulas. If `NULL`, survey
#'   catchability is assumed to be time-invariant (i.e., not influenced by environmental variables).
#' @param Use_srv_selex_prior Integer (0 or 1). Flag to enable/disable survey selectivity priors.
#'   When set to 1, applies log-normal priors to survey selectivity parameters as specified
#'   in \code{srv_selex_prior}. When set to 0, no priors are applied.
#' @param srv_selex_prior Data frame containing prior specifications for survey selectivity parameters.
#'   Must include columns: \code{region} (region index), \code{fleet} (fleet index),
#'   \code{block} (time block index), \code{sex} (sex index), \code{par} (parameter index),
#'   \code{mu} (prior mean on natural scale), and \code{sd} (prior standard deviation on log scale).
#'   Each row specifies a log-normal prior N(log(mu), sd) for one selectivity parameter.
#'   Only parameters with rows in this data frame will have priors applied.
#'
#' If both `srv_q_formula` and `srv_q_cov_dat` are non-`NULL`, the model constructs time-varying design matrices
#' for each region and fleet based on the provided formulas and environmental covariates. A coefficient array
#' (`srv_q_coeff`) and a mapping array (`map_srv_q_coeff`) are created to estimate and track the associated
#' regression coefficients. The design matrix is stored in `srv_q_env`, a 4D array indexed by
#' \code{[region, year, fleet, covariate]}.
#'
#' If either argument is `NULL`, environmental covariate effects are excluded and survey catchability is treated
#' as constant over time.
#'
#' \strong{Important:} All covariate time series in `srv_q_cov_dat` must:
#' \itemize{
#'   \item Be numeric vectors with a length equal to the number of years in the model.
#'   \item Align to the same years across all covariates.
#'   \item Contain no missing values; users must impute or interpolate missing covariate values prior to use. For years in which the index is not used, values can be set at 0.
#' }
#'
#' Covariates that are defined but not used in any formula can be filled with zeros (e.g., \code{rep(0, n_yrs)}).
#' This avoids issues with list structure but does not affect the design matrix or model results.
#'
#' \strong{Example formulas:}
#' \itemize{
#'   \item \code{"Region_1_Fleet_1"} = \code{~ 0 + poly(env1_r1_f1, 3) + env2_r1_f1} uses a 3rd-degree
#'         polynomial for \code{env1_r1_f1} and a linear term for \code{env2_r1_f1}.
#'   \item \code{"Region_2_Fleet_1"} = \code{~ 0 + env1_r2_f1 + env2_r2_f1} includes additive effects of
#'         two covariates.
#'   \item \code{"Region_3_Fleet_2"} = \code{~ NULL} disables environmental covariates for that fleet-region.
#' }
#' @param t_srv Survey timing in fractions (n_regions * n_srv_fleets; default is 0.5)
#' @param cont_tv_srv_sel_penalty Whether or not to apply continuous time-varying selectivity penalties (if cont_tv_srv_sel > 0)
#' @param srvsel_devs_shared_ages List object for specifying which ages are shared when selectivity deviations are semi-parametric (e.g., list(1:5, 6:10, 11:30) specifies that ages 1-5, 6-10, and 11-30 have the same deviations.)
#'
#' @export Setup_Mod_Srvsel_and_Q
#' @importFrom stringr str_detect
#' @family Model Setup
Setup_Mod_Srvsel_and_Q <- function(input_list,
                                   cont_tv_srv_sel = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                   srv_sel_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                   srv_sel_model,
                                   Use_srv_q_prior = 0,
                                   srv_q_prior = NA,
                                   srv_q_blocks = paste("none_Fleet_", 1:input_list$data$n_srv_fleets, sep = ''),
                                   srvsel_pe_pars_spec = NULL,
                                   srv_fixed_sel_pars_spec,
                                   srv_q_spec = NULL,
                                   srv_sel_devs_spec = NULL,
                                   corr_opt_semipar = NULL,
                                   srv_q_formula = NULL,
                                   srv_q_cov_dat = NULL,
                                   Use_srv_selex_prior = 0,
                                   srv_selex_prior = NULL,
                                   t_srv = array(0.5, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets)),
                                   cont_tv_srv_sel_penalty = TRUE,
                                   srvsel_devs_shared_ages = NULL,
                                   ...
                                   ) {

  messages_list <<- character(0) # string to attach to for printing messages
  starting_values <- list(...)

  # Input Validation --------------------------------------------------------

  # Selectivity
  # Selectivity Type
  if(is.null(input_list$data$Selex_Type)) stop("Selectivity type (age or length-based) has not been specified yet! Make sure to first specify biological inputs with Setup_Mod_Biologicals.")

  # Continuous Selectivity Deviations
  if(!is.null(srvsel_pe_pars_spec)) if(length(srvsel_pe_pars_spec) != input_list$data$n_srv_fleets) stop("srvsel_pe_pars_spec is not length n_srv_fleets")
  if(!is.null(srv_sel_devs_spec)) if(length(srv_sel_devs_spec) != input_list$data$n_srv_fleets) stop("srv_sel_devs_spec is not length n_srv_fleets")
  if(!is.null(corr_opt_semipar)) if(length(corr_opt_semipar) != input_list$data$n_srv_fleets) stop("corr_opt_semipar is not length n_srv_fleets")

  # Catchability Priors
  if(!Use_srv_q_prior %in% c(0,1)) stop("Values for Use_srv_q_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking catchability priors
  if(Use_srv_q_prior == 1) {
    required_cols <- c("region", "fleet", "block", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(srv_q_prior))
    if(length(missing_cols) > 0) {
      stop("srv_q_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Survey Catchability priors are: ", ifelse(Use_srv_q_prior == 0, "Not Used", "Used"))

  # Selectivity Priors
  if(!Use_srv_selex_prior %in% c(0,1)) stop("Values for Use_srv_selex_prior are not valid. They are == 0 (don't use prior), or == 1 (use prior)")
  # Checking selectivity priors
  if(Use_srv_selex_prior == 1) {
    required_cols <- c("region", "fleet", "block", "sex", "par", "mu", "sd")
    missing_cols <- setdiff(required_cols, names(srv_selex_prior))
    if(length(missing_cols) > 0) {
      stop("srv_selex_prior is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  collect_message("Survey Selectivity priors are: ", ifelse(Use_srv_selex_prior == 0, "Not Used", "Used"))


  # Continuous Time-Varying Selectivity Options -----------------------------
  # define for continuous time-varying selectivity
  cont_tv_srv_sel_mat <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets))
  cont_tv_map <- data.frame(type = c("none", "iid", "rw", "3dmarg", "3dcond", "2dar1"), num = c(0,1,2,3,4,5)) # set up values we map to

  for(i in 1:length(cont_tv_srv_sel)) {
    # Extract out components from list
    tmp <- cont_tv_srv_sel[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))
    cont_tv_type <- tmp_vec[1] # get continuous selex type
    fleet <- as.numeric(tmp_vec[3]) # extract fleet index

    # Validate options
    if(!fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for cont_tv_srv_sel This needs to be specified as timevarytype_Fleet_x")
    if(!cont_tv_type %in% c(cont_tv_map$type)) stop("cont_tv_srv_sel is not correctly specified. This needs to be one of these: none, iid, rw, 3dmarg, 3dcond, 2dar1 (the timevarytypes) and specified as timevarytype_Fleet_x")

    # Input options
    cont_tv_srv_sel_mat[,fleet] <- cont_tv_map$num[which(cont_tv_map$type == cont_tv_type)]
    collect_message("Continuous survey time-varying selectivity specified as: ", cont_tv_type, " for survey fleet ", fleet)
  }

  if(any(cont_tv_srv_sel_mat > 0) && is.null(srvsel_pe_pars_spec) && is.null(srv_sel_devs_spec)) stop("Continuous time-varying selectivity specified, but srvsel_pe_pars_spec and/or srv_sel_devs_spec is NULL (i.e., not specified)!")

  # Blocked Time-Varying Selectivity Options --------------------------------
  srv_sel_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_sel_blocks)) {

    # Extract out components from list
    tmp <- srv_sel_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Selectivity Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_sel_blocks_arr[,,fleet] <- 1 # input only 1 survey time block
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

      srv_sel_blocks_arr[,years,fleet] <- block_val
    }
  }

  if(any(is.na(srv_sel_blocks_arr))) stop("Survey Selectivity Blocks are returning an NA. Did you update the year range of srv_sel_blocks?")
  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Selectivity Time Blocks for survey", f, "is specified at:", length(unique(srv_sel_blocks_arr[,,f]))))

  # Selectivity Functional Forms --------------------------------------------
  sel_map <- data.frame(sel = c('logist1', "gamma", "exponential", "logist2", "dbnrml"), num = c(0,1,2,3,4)) # set up values we can map to
  srv_sel_model_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_sel_model)) {

    # Extract out survey selectivity components from vector
    tmp_sel_form <- srv_sel_model[i]
    tmp_sel_form_vec <- unlist(strsplit(tmp_sel_form, "_")) # split string
    sel_form <- tmp_sel_form_vec[1] # get selectivity type

    # get fleet index
    tmp_fleet <- as.numeric(tmp_sel_form_vec[3])
    # get block index
    tmp_block <- if(length(tmp_sel_form_vec) == 5) as.numeric(tmp_sel_form_vec[5])

    # validate options
    if(!sel_form %in% c(sel_map$sel)) stop("srv_sel_model is not correctly specified. This needs to be one of these: logist1, gamma, exponential, logist2, dbnrml (the seltypes) and specified as seltype_Fleet_x")
    if(!tmp_fleet %in% c(1:input_list$data$n_srv_fleets)) stop("Invalid fleet specified for srv_sel_model This needs to be specified as seltype_Fleet_x or seltype_Fleet_x_Block_x (if blocks are specified to change for a fleet)")

    # Input options
    if(is.null(tmp_block)) srv_sel_model_arr[,,tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)] # same selectivity form across blocks
    else srv_sel_model_arr[,which(srv_sel_blocks_arr[,,tmp_fleet] == tmp_block),tmp_fleet] <- sel_map$num[which(sel_map$sel == sel_form)]
    rm(tmp_block) # remove tmp block to start next loop
    collect_message("Survey selectivity functional form specified as:", sel_form, " for survey fleet ", tmp_fleet)
  }

  # Validate that blocks and continuous time-variation aren't both specified for same fleet
  for(f in 1:input_list$data$n_srv_fleets) {
    has_blocks <- length(unique(srv_sel_blocks_arr[1,,f])) > 1
    has_cont_tv <- cont_tv_srv_sel_mat[1,f] != 0  # 0 = "none"
    if(has_blocks && has_cont_tv) {
      stop("Fleet ", f, " has both selectivity blocks and continuous time-varying selectivity specified. ",
           "These are mutually exclusive - choose one approach to time-variation.")
    }
  }

  # Blocked Catchability Options --------------------------------------------
  srv_q_blocks_arr <- array(NA, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
  for(i in 1:length(srv_q_blocks)) {

    # Extract out components from list
    tmp <- srv_q_blocks[i]
    tmp_vec <- unlist(strsplit(tmp, "_"))

    # Vakudate option
    if(!tmp_vec[1] %in% c("none", "Block")) stop("Survey Catchability Blocks not correctly specified. This should be either none_Fleet_x or Block_x_Year_x-y_Fleet_x")

    # extract out fleets if constant
    if(tmp_vec[1] == "none") {
      fleet <- as.numeric(tmp_vec[3]) # get fleet number
      srv_q_blocks_arr[,,fleet] <- 1 # input only 1 survey catchability time block
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
      srv_q_blocks_arr[,years,fleet] <- block_val # input catchability time block
    }
  }
  if(any(is.na(srv_q_blocks_arr))) stop("Survey Catchability Blocks are returning an NA. Did you update the year range of srv_q_blocks?")
  for(f in 1:input_list$data$n_srv_fleets) collect_message(paste("Survey Catchability Time Blocks for survey", f, "is specified at:", length(unique(srv_q_blocks_arr[,,f]))))

  # Covariate Catchability Options ------------------------------------------
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) collect_message("Using covariates to predict survey catchability")

  # Figure out the total number of regression coefficients that could be estimated
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) {
    n_srv_q_cov <- max(sapply(names(srv_q_formula), function(key) { # sapply to extract names from formula
      tmp_formula <- srv_q_formula[[key]] # get formula
      var_names <- all.vars(tmp_formula) # get var names
      tmp_dat <- data.frame(srv_q_cov_dat[var_names]) # make dataframe
      ncol(model.matrix(tmp_formula, data = tmp_dat)) # figure out number of columns for formula (max number of coefficients to estimate)
    }))
  } else {
    do_srv_q_cov <- 0 # Indicator for whether covariates are included into survey catchability
    n_srv_q_cov <- 1 # dummy to initialize the array
  }

  # Catchability covariate containers
  srv_q_cov <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets, n_srv_q_cov)) # environmental time series
  srv_q_coeff <- array(0, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets, n_srv_q_cov)) # coefficients to be estimated
  map_srv_q_coeff <- array(NA, dim = c(input_list$data$n_regions, input_list$data$n_srv_fleets, n_srv_q_cov)) # coefficients to be mapped off

  # Loop through to map stuff off and populate containers
  if(!is.null(srv_q_cov_dat) && !is.null(srv_q_formula)) {

    # Validate covariate length
    cov_lengths <- lengths(srv_q_cov_dat)

    # Validate options
    # Check all covariates are the same length
    if (length(unique(cov_lengths)) != 1) stop("All covariates in 'srv_q_cov_dat' must have the same length. If some years are missing data, either impute some value, or set at 0 (if it is not used in the calculation).")
    # Check that length matches the model year structure
    if (unique(cov_lengths) != length(input_list$data$years)) stop(paste0("Covariate length mismatch: expected ",  length(input_list$data$years),  " years but got ", unique(cov_lengths),  "."))

    do_srv_q_cov <- 1 # Indicator for whether covariates are included into survey catchability
    coeff_counter <- 0 # setup counter for mapping

    for(r in 1:n_regions) {
      for(f in 1:n_srv_fleets) {

        # Get key to index
        key <- paste(paste("Region", r, sep = "_"), "_Fleet_", f, sep = "")
        # get temporary formula
        tmp_formula <- srv_q_formula[[key]]
        # extract variable names
        var_names <- all.vars(tmp_formula)
        if(length(var_names) == 0) next # skip if no variables
        # get environmental covariates from environmental data list, based on model formula
        tmp_dat <- data.frame(srv_q_cov_dat[var_names])
        # Generate design matrix
        tmp_design_mat <- model.matrix(tmp_formula, data = tmp_dat)
        # store covariate effects into container
        srv_q_cov[r,,f,1:ncol(tmp_design_mat)] <- tmp_design_mat

        # setup mapping - assign unique counter values for each coefficient
        for(i in 1:ncol(tmp_design_mat)) {
          coeff_counter <- coeff_counter + 1
          map_srv_q_coeff[r,f,i] <- coeff_counter
        } # end i loop

      } # end sf loop
    } # end r loop
  } # if using covariates

  # Populate Data List ------------------------------------------------------
  input_list$data$cont_tv_srv_sel <- cont_tv_srv_sel_mat
  input_list$data$cont_tv_srv_sel_penalty <- cont_tv_srv_sel_penalty
  input_list$data$srv_sel_blocks <- srv_sel_blocks_arr
  input_list$data$srv_sel_model <- srv_sel_model_arr
  input_list$data$srv_q_blocks <- srv_q_blocks_arr
  input_list$data$srv_q_prior <- srv_q_prior
  input_list$data$Use_srv_q_prior <- Use_srv_q_prior
  input_list$data$do_srv_q_cov <- do_srv_q_cov
  input_list$data$srv_q_cov <- srv_q_cov
  input_list$data$Use_srv_selex_prior <- Use_srv_selex_prior
  input_list$data$srv_selex_prior <- srv_selex_prior
  input_list$data$t_srv <- t_srv


  # Populate Parameter List -------------------------------------------------
  # Figure out number of selectivity parameters for a given functional form
  unique_srvsel_vals <- unique(as.vector(input_list$data$srv_sel_model))
  sel_pars_vec <- vector() # create empty vector to populate

  for(i in 1:length(unique_srvsel_vals)) {
    if(unique_srvsel_vals[i] %in% c(2)) sel_pars_vec[i] <- 1 # exponential
    if(unique_srvsel_vals[i] %in% c(0,1,3)) sel_pars_vec[i] <- 2 # logistic or gamma
    if(unique_srvsel_vals[i] == 4) sel_pars_vec[i] <- 6 # double normal
  } # end i loop

  max_srvsel_blks <- max(apply(input_list$data$srv_sel_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey selectivity blocks for a given reigon and fleet
  max_srvsel_pars <- max(sel_pars_vec) # maximum number of selectivity parameters across all forms
  if("ln_srv_fixed_sel_pars" %in% names(starting_values)) input_list$par$ln_srv_fixed_sel_pars <- starting_values$ln_srv_fixed_sel_pars
  else input_list$par$ln_srv_fixed_sel_pars <- array(0, dim = c(input_list$data$n_regions, max_srvsel_pars, max_srvsel_blks, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability
  max_srvq_blks <- max(apply(input_list$data$srv_q_blocks, c(1,3), FUN = function(x) length(unique(x)))) # figure out maximum number of survey catchability blocks for a given reigon and fleet
  if("ln_srv_q" %in% names(starting_values)) input_list$par$ln_srv_q <- starting_values$ln_srv_q
  else input_list$par$ln_srv_q <- array(0, dim = c(input_list$data$n_regions, max_srvq_blks, input_list$data$n_srv_fleets))

  # Survey selectivity process error parameters
  if("srvsel_pe_pars" %in% names(starting_values)) input_list$par$srvsel_pe_pars <- starting_values$srvsel_pe_pars
  else input_list$par$srvsel_pe_pars <- array(0, dim = c(input_list$data$n_regions, max(max_srvsel_pars, 4), input_list$data$n_sexes, input_list$data$n_srv_fleets)) # dimensioned 4 as the max number of pars for process errors (e.g., sigmas), and then just map off if not using

  # Survey selectivity deviations
  if(input_list$data$Selex_Type == 0) bins <- length(input_list$data$ages) # age based deviations
  if(input_list$data$Selex_Type == 1) bins <- length(input_list$data$lens) # length based deviations
  if("ln_srvsel_devs" %in% names(starting_values)) input_list$par$ln_srvsel_devs <- starting_values$ln_srvsel_devs
  else input_list$par$ln_srvsel_devs <- array(0, dim = c(input_list$data$n_regions, length(input_list$data$years) + input_list$data$n_proj_yrs_devs, bins, input_list$data$n_sexes, input_list$data$n_srv_fleets))

  # Survey catchability covariate effects
  if("srv_q_coeff" %in% names(starting_values)) input_list$par$srv_q_coeff <- starting_values$srv_q_coeff
  else input_list$par$srv_q_coeff <- srv_q_coeff # input parameter array

  # Mapping Options ---------------------------------------------------------
  input_list$map$srv_q_coeff <- factor(map_srv_q_coeff) # set up mapping for catchability covariate
  input_list <- do_srv_fixed_sel_pars_mapping(input_list, srv_fixed_sel_pars_spec)
  input_list <- do_srv_q_mapping(input_list, srv_q_spec)
  input_list <- do_srvsel_pe_pars_mapping(input_list, srvsel_pe_pars_spec, corr_opt_semipar)
  input_list <- do_srvsel_devs_mapping(input_list, srv_sel_devs_spec, srvsel_devs_shared_ages)


  # Print Messages ----------------------------------------------------------
  if(input_list$verbose) for(msg in messages_list) message(msg)

  return(input_list)
}
