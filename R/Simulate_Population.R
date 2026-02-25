#' Constructs simulation objects in a new simulation environment for use in simulation functions
#'
#' @param sim_list Simulation list objects
#'
#' @returns A new simulation environment with objects from sim_list
#' @export Setup_sim_env
#'
#' @examples
#' \dontrun{
#' sim_env <- Setup_sim_env(sim_list)
#' }
#' @family Simulation Setup
Setup_sim_env <- function(sim_list) {

  sim_env <- new.env(parent = parent.frame()) # define new environment for simulation

  # Get SPoRC functions in simulation environment
  sim_env$generate_initial_age_structure <- generate_initial_age_structure
  sim_env$generate_recruitment <- generate_recruitment
  sim_env$apply_pop_dy <- apply_pop_dy
  sim_env$generate_fishery_catch_comp_idx <- generate_fishery_catch_comp_idx
  sim_env$generate_survey_comp_idx <- generate_survey_comp_idx
  sim_env$release_conv_tags <- release_conv_tags
  sim_env$generate_fishery_conv_tags_recap <- generate_fishery_conv_tags_recap
  sim_env$Get_Det_Recruitment <- Get_Det_Recruitment
  sim_env$Get_Init_NAA <- Get_Init_NAA
  sim_env$Get_Tagging_Mortality <- Get_Tagging_Mortality
  sim_env$predict_sim_fish_iss_fmort <- predict_sim_fish_iss_fmort
  sim_env$rho_trans <- rho_trans
  sim_env$simulate_comps <- simulate_comps
  sim_env$simulate_conv_tag_fish_recaptures <- simulate_conv_tag_fish_recaptures

  # output into simulation environment
  list2env(sim_list, envir = sim_env)

  return(sim_env)
}

#' Simulate Age or Length Compositions
#'
#' Generates Observed fish compositions by age or length for a given region, year, fleet, and simulation iteration.
#' Supports multinomial, Dirichlet-multinomial, and logistic-normal likelihoods, with optional ageing error applied.
#'
#' @param r Integer. Region index.
#' @param y Integer. Year index.
#' @param f Integer. Fleet index.
#' @param sim Integer. Simulation iteration index.
#' @param Exp Array. Expected compositions (age or length) with dimensions [region, year, category, sex, fleet, sim].
#' @param ISS Array. Sample size (integer) for the Observed compositions with dimensions [region, year, sex, fleet, sim].
#' @param AgeingError Array. Ageing error matrix for each year, dimensions [year, category, category, sim].
#' @param comp_like Integer vector. Composition likelihood type per fleet: 0 = multinomial, 1 = Dirichlet-multinomial, 2-4 = logistic-normal.
#' @param ln_theta Array. Log-variance parameter for compositions per region, sex, and fleet, dimensions [region, sex, fleet].
#' @param corr_pars Array. Correlation parameters for logistic-normal likelihood, dimensions [region, sex, fleet, ?].
#' @param ln_theta_agg Numeric vector. Log-variance parameter for aggregated compositions per fleet.
#' @param corr_pars_agg Numeric vector. Correlation parameters for aggregated logistic-normal likelihood per fleet.
#' @param comp_type Integer array. Composition type: 0 = aggregated across regions, 1 = split by sex, 2 = joint across sexes, dimensions [year, fleet].
#' @param n_sexes Integer. Number of sexes.
#' @param n_regions Integer. Number of regions.
#' @param n_cat Integer. Number of categories (ages or lengths).
#' @param Obs Array. Observed compositions array to fill, same dimensions as `Exp`.
#' @param age_or_len Integer. Flag to indicate if ageing error should be applied: 0 = apply ageing error (for ages), 1 = do not apply (for lengths).
#' @param seas Intege. Seasonal index
#' @param pop_specific Boolean on whether or not composition data are population-specific
#' @param n_pop Integer. Number of populations.
#'
#' @return Array of Observed compositions with the same dimensions as `Obs`, updated with simulated Observations.
#'
#' @details
#' The function handles three cases based on `comp_type`:
#' 1. Split by sex (comp_type = 1): compositions are simulated separately for each sex.
#' 2. Joint compositions across sexes (comp_type = 2): compositions simulated jointly and multiplied by a kronecker matrix for logistic-normal or Dirichlet-multinomial likelihoods.
#' 3. Aggregated across regions (comp_type = 0): only applied in the last region and averages across regions and sexes.
#'
#' The function normalizes expected compositions, applies the selected likelihood (`comp_like`), and multiplies by `AgeingError` when applicable.
#'
#' @keywords internal
simulate_comps <- function(r,
                           y,
                           f,
                           seas,
                           sim,
                           Exp,
                           ISS,
                           AgeingError,
                           comp_like,
                           ln_theta,
                           corr_pars,
                           ln_theta_agg,
                           corr_pars_agg,
                           comp_type,
                           n_sexes,
                           n_pop,
                           n_regions,
                           n_cat,
                           Obs,
                           age_or_len = 0,
                           pop_specific = FALSE) {

  if(comp_type[y,f] == 999 || comp_like[f] == 999) return(Obs)

  # helper functions
  get_expected <- function(prob_vec) prob_vec / sum(prob_vec)
  apply_error <- function(mat, age_or_len, AgeingError) {
    if(age_or_len == 0) return(mat %*% AgeingError)
    if(age_or_len == 1) return(mat)
  }

  if(age_or_len == 0) {
    if(comp_type[y,f] %in% c(0,1)) age_error_mat <- AgeingError[y,,,sim] # aggregated or split
    if(comp_type[y,f] == 2) age_error_mat <- kronecker(diag(n_sexes), AgeingError[y,,,sim]) # joint
  } else if(age_or_len == 1) age_error_mat <- NULL # length compositions

  if(pop_specific == FALSE) {

    # Split by sex
    if(comp_type[y,f] == 1) {
      for(s in 1:n_sexes) {

        tmp_prob <- apply(Exp[,r,y,seas,,s,f,sim, drop = FALSE], 5, sum) # extract compositions

        # multinomial
        if(comp_like[f] == 0) {
          Obs[r,y,seas,,s,f,sim] <- array(
            apply_error(as.vector(
              stats::rmultinom(n = 1, ISS[r,y,seas,s,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
            dim = dim(Obs[r,y,seas,,s,f,sim, drop = FALSE])
          )

          # dirichlet-multinomial
        } else if(comp_like[f] == 1) {
          Obs[r,y,seas,,s,f,sim] <- array(
            apply_error(as.vector(
              rdirM(
                n = 1,
                N = ISS[r,y,seas,s,f,sim],
                alpha = (exp(ln_theta[r,s,f]) * ISS[r,y,seas,s,f,sim]) * get_expected(tmp_prob)
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[r,y,seas,,s,f,sim, drop = FALSE])
          )

          # logistic normal
        } else if(comp_like[f] %in% 2:4) {
          Obs[r,y,seas,,s,f,sim] <- array(
            apply_error(as.vector(
              rlogistnormal(
                exp = get_expected(tmp_prob),
                pars = c(exp(ln_theta[r,s,f]), rho_trans(corr_pars[r,s,f,])),
                comp_like = comp_like[f],
                n_sexes = n_sexes
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[r,y,seas,,s,f,sim, drop = FALSE])
          )
        }

      } # end s loop
    } # end split by sex

    # Joint compositions
    if(comp_type[y,f] == 2) {

      tmp_prob <- apply(Exp[,r,y,seas,,,f,sim, drop = FALSE], c(5,6), sum) # extract compositions

      # multinomial
      if(comp_like[f] == 0) {
        Obs[r,y,seas,,,f,sim] <- array(
          apply_error(as.vector(stats::rmultinom(1, ISS[r,y,seas,1,f,sim], get_expected(tmp_prob))),
                      age_or_len, age_error_mat),
          dim = dim(Obs[r,y,seas,,,f,sim, drop = FALSE])
        )

        # dirichlet-multinomial
      } else if(comp_like[f] == 1) {
        Obs[r,y,seas,,,f,sim] <- array(
          apply_error(as.vector(
            rdirM(
              n = 1,
              N = ISS[r,y,seas,1,f,sim],
              alpha = (exp(ln_theta[r,1,f]) * ISS[r,y,seas,1,f,sim]) * get_expected(tmp_prob)
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[r,y,seas,,,f,sim, drop = FALSE])
        )

        # logistic normal
      } else if(comp_like[f] %in% 2:4) {
        Obs[r,y,seas,,,f,sim] <- array(
          apply_error(as.vector(
            rlogistnormal(
              exp = get_expected(tmp_prob),
              pars = c(exp(ln_theta[r,1,f]), rho_trans(corr_pars[r,1,f,])),
              comp_like = comp_like[f],
              n_sexes = n_sexes
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[r,y,seas,,,f,sim, drop = FALSE])
        )
      }

    } # end joint compositions

    # Aggregated comps across regions
    if(r == n_regions && comp_type[y,f] == 0) {

      # extract compositions
      tmp_prob <- apply(Exp[,,y,seas,,,f,sim, drop = FALSE], 5, sum)
      tmp_prob <- tmp_prob / sum(tmp_prob)

      # multinomial
      if(comp_like[f] == 0) {
        Obs[1,y,seas,,1,f,sim] <- array(
          apply_error(as.vector(stats::rmultinom(1, ISS[1,y,seas,1,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
          dim = dim(Obs[1,y,seas,,1,f,sim, drop = FALSE])
        )

        # dirichlet-multinomial
      } else if(comp_like[f] == 1) {
        Obs[1,y,seas,,1,f,sim] <- array(
          apply_error(as.vector(
            rdirM(
              n = 1,
              N = ISS[1,y,seas,1,f,sim],
              alpha = (exp(ln_theta_agg[f]) * ISS[1,y,seas,1,f,sim]) * get_expected(tmp_prob)
            )
          ), age_or_len, age_error_mat),
          dim = dim(Obs[1,y,seas,,1,f,sim, drop = FALSE])
        )

        # logistic normal
      } else if(comp_like[f] %in% 2:4) {
        Obs[1,y,seas,,1,f,sim] <- array(
          apply_error(as.vector(
            rlogistnormal(
              exp = get_expected(tmp_prob),
              pars = c(exp(ln_theta_agg[f]), rho_trans(corr_pars_agg[f])),
              comp_like = comp_like[f],
              n_sexes = n_sexes
            )
          )),
          dim = dim(Obs[1,y,seas,,1,f,sim, drop = FALSE])
        )
      }
    }

  } # end if not population specific

  return(Obs)
}

#' Simulate conventional tag recaptures for fishery fleets
#'
#' Simulates observed tag recaptures from fishery fleets given predicted recapture
#' arrays, supporting multiple likelihood structures (Poisson, Negative Binomial,
#' Multinomial, Dirichlet-Multinomial) and flexible attribute reporting levels
#' (population, age, sex). The function marginalizes over unreported dimensions
#' based on the specified \code{tag_recaptures_attr} string.
#'
#' @param conv_fish_tag_like Integer specifying the likelihood for tag recaptures.
#'   \itemize{
#'     \item \code{0} Poisson
#'     \item \code{1} Negative Binomial
#'     \item \code{2} Multinomial, release conditioned
#'     \item \code{3} Multinomial, recovery conditioned
#'     \item \code{4} Dirichlet-Multinomial, release conditioned
#'     \item \code{5} Dirichlet-Multinomial, recovery conditioned
#'   }
#' @param tag_recaptures_attr Character string specifying which biological attributes
#'   are recorded at recapture. Constructed from any combination of \code{"p"} (population),
#'   \code{"a"} (age), and \code{"s"} (sex), joined by underscores. Region and fleet are
#'   always retained. Supported values: \code{"p_a_s"}, \code{"a_s"}, \code{"p_a"},
#'   \code{"p_s"}, \code{"a"}, \code{"s"}, \code{"p"}, \code{"none"}. Dimensions not
#'   present in the string are marginalized out and assigned to index 1 in the output array.
#' @param conv_tagged_fish Array of tagged fish at release with dimensions
#'   \code{[tc, pop, region, sex, sim]}. Used to determine the number of tags at liberty
#'   for release-conditioned likelihoods.
#' @param pred_conv_tag_fish_recap Array of predicted tag recaptures with dimensions
#'   \code{[year, season, cohort, pop, region, age, sex, fleet, sim]}.
#' @param obs_conv_tag_fish_recap Array of observed tag recaptures with the same
#'   dimensions as \code{pred_conv_tag_fish_recap}. Simulated values are written into
#'   this array and the updated array is returned.
#' @param ln_conv_fish_tag_theta Numeric. Log of the overdispersion parameter used
#'   in Negative Binomial (size = \code{exp(ln_conv_fish_tag_theta)}) and
#'   Dirichlet-Multinomial (\eqn{\theta = \exp(\text{ln\_conv\_fish\_tag\_theta})}) likelihoods.
#'   Ignored for Poisson and Multinomial.
#' @param ry Integer. Recovery year index into the recapture arrays.
#' @param rseas Integer. Recovery season index into the recapture arrays.
#' @param tc Integer. Tag cohort index, used to index into \code{pred_conv_tag_fish_recap}
#'   and \code{conv_tagged_fish}.
#' @param sim Integer. Simulation replicate index.
#' @param n_pop Integer. Number of populations.
#' @param n_regions Integer. Number of regions.
#' @param n_ages Integer. Number of age classes.
#' @param n_sexes Integer. Number of sexes.
#' @param n_fish_fleets Integer. Number of fishery fleets.
#'
#' @returns The \code{obs_conv_tag_fish_recap} array with simulated recaptures filled in
#'   at indices \code{[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim]},
#'   where marginalized dimensions are fixed at index 1.
#'
#' @details
#' For release-conditioned likelihoods (\code{2}, \code{4}), the total number of tags
#' at liberty is taken from \code{conv_tagged_fish} and predicted recaptures are expressed
#' as proportions. A "not recaptured" bin is appended to complete the probability vector
#' before drawing, then removed before assignment.
#'
#' For recovery-conditioned likelihoods (\code{3}, \code{5}), the total number of
#' recaptures is taken directly from the sum of \code{pred_conv_tag_fish_recap} and
#' no "not recaptured" bin is needed.
#'
#' Marginalization is performed by reshaping the flat predicted recapture vector into
#' a 5-dimensional array of shape \code{c(n_pop, n_regions, n_ages, n_sexes, n_fish_fleets)}
#' and summing over dimensions absent from \code{tag_recaptures_attr}.
#'
#' @keywords internal
simulate_conv_tag_fish_recaptures <- function(conv_fish_tag_like,
                                              tag_recaptures_attr,
                                              conv_tagged_fish,
                                              pred_conv_tag_fish_recap,
                                              obs_conv_tag_fish_recap,
                                              ln_conv_fish_tag_theta,
                                              ry,
                                              rseas,
                                              tc,
                                              sim,
                                              n_pop,
                                              n_regions,
                                              n_ages,
                                              n_sexes,
                                              n_fish_fleets
                                              ) {

  # get full dimensions of tag recaptures we simulate
  full_dims  <- c(n_pop, n_regions, n_ages, n_sexes, n_fish_fleets)
  attr_parts <- strsplit(tag_recaptures_attr, "_")[[1]]

  # Which of the 5 free dims (pop=1, region=2, age=3, sex=4, fleet=5) to retain
  # Region and fleet are always kept; pop/age/sex kept only if present in attr string
  keep_dims <- c(
    if("p" %in% attr_parts) 1,
    2,                            # region always kept
    if("a" %in% attr_parts) 3,
    if("s" %in% attr_parts) 4,
    5                             # fleet always kept
  )

  # get obs slice indices: full range if dim is kept, fixed at 1 if marginalized out
  pop_idx <- if("p" %in% attr_parts) seq_len(n_pop)   else 1
  age_idx <- if("a" %in% attr_parts) seq_len(n_ages)  else 1
  sex_idx <- if("s" %in% attr_parts) seq_len(n_sexes) else 1
  reg_idx <- seq_len(n_regions)   # always full
  flt_idx <- seq_len(n_fish_fleets) # always full

  # Function to marginalize tag recaptures. If dims == 5, then keep dimensions, otherwise,
  # marginalize and retain the kept dimensions
  marginalize <- function(vals) {
    tmp <- array(vals, dim = full_dims)
    if(length(keep_dims) < 5) apply(tmp, keep_dims, sum) else tmp
  }

  # Poisson or Neg Bin
  if(conv_fish_tag_like %in% c(0, 1)) {
    lambda <- marginalize(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim]) # get lambda / mu parameter
    # input and simulate
    obs_conv_tag_fish_recap[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim] <-
      if(conv_fish_tag_like == 0) {
        rpois(n = length(lambda), lambda = lambda)
      } else {
        rnbinom(n = length(lambda), mu = lambda, size = exp(ln_conv_fish_tag_theta))
      }
  }

  # Multinomial or Dirichlet Multinomial (Release conditioned)
  if(conv_fish_tag_like %in% c(2, 4)) {
    tmp_n_tags_rel <- round(sum(conv_tagged_fish[tc, , , , sim])) # get sample size
    tmp_recap      <- marginalize(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim] / tmp_n_tags_rel) # marginalize
    tmp_probs      <- c(tmp_recap, 1 - sum(tmp_recap)) # get non-recaptured state

    # simualte
    tmp_sim_recap <-
      if(conv_fish_tag_like == 2) {
        stats::rmultinom(1, tmp_n_tags_rel, tmp_probs)
      } else {
        rdirM(n = 1, N = tmp_n_tags_rel, exp(ln_conv_fish_tag_theta) * tmp_n_tags_rel * tmp_probs)
      }

    # Drop the "not recaptured" bin and restore array shape
    tmp_sim_recap <- array(tmp_sim_recap[-length(tmp_sim_recap)], dim(tmp_recap))
    obs_conv_tag_fish_recap[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim] <- tmp_sim_recap
  }

  # Multinomial or Dirichlet Multinomial (Recovery conditioned)
  if(conv_fish_tag_like %in% c(3, 5)) {
    tmp_n_tags_recap <- round(sum(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim])) # get sample size
    tmp_probs        <- marginalize(pred_conv_tag_fish_recap[ry, rseas, tc, , , , , , sim] / tmp_n_tags_recap)

    # simulate
    tmp_sim_recap <-
      if(conv_fish_tag_like == 3) {
        stats::rmultinom(1, tmp_n_tags_recap, c(tmp_probs))
      } else {
        rdirM(n = 1, N = tmp_n_tags_recap, exp(ln_conv_fish_tag_theta) * tmp_n_tags_recap * c(tmp_probs))
      }

    # input
    obs_conv_tag_fish_recap[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim] <- tmp_sim_recap
  }

  return(obs_conv_tag_fish_recap)

}

#' Title Initialize Age Structure for Simulation
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation Environment
#' @keywords internal
generate_initial_age_structure <- function(y,
                                           sim,
                                           sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {
    tmp_ln_init_devs <- NULL
    for(p in 1:n_pop) {

      # reset deviations for each population (draws for each popn)
      if(n_pop > 1) tmp_ln_init_devs <- NULL

      for(r in 1:n_regions) {

        # if local DD and n_pop = 1, reset deviations for each region (draws for each region, but if n_pop > 1,
        # shares deviations across regions withn a given population)
        if(n_pop == 1 && init_dd == 0) tmp_ln_init_devs <- NULL

        if(exists("ln_InitDevs_input")) { # if exists in environment, then use input
          tmp_ln_init_devs <- ln_InitDevs_input[p,r,,sim]
        } else { # simulate new initial age devs otherwise

          # index ln_sigmaR, if n_pop == 1 and local DD, use region specific rates,
          # otherwise use 1 (i.e., n_pop > 1 or init_dd == 1)
          sigma_idx <- ifelse(n_pop == 1 && init_dd == 0, r, 1)

          # simulate initial age deviations
          if(is.null(tmp_ln_init_devs)) tmp_ln_init_devs <- stats::rnorm(n_ages-1, 0, exp(ln_sigmaR[1,p,sigma_idx]))
        }

        # input age deviations
        if(R0[p,r,y,sim] != 0) {
          sim_env$ln_InitDevs[p,r,,sim] <- tmp_ln_init_devs
        } else sim_env$ln_InitDevs[p,r,,sim] <- 0

      } # end r loop
    } # end p loop


    # Get initial fished NAA
    Init_Fished_NAA = Get_Init_NAA(
      init_age_strc = init_age_strc, # initial age structure
      init_iter = n_ages * 5, # if init_age_strc == 0, number of iterations to run
      n_pop = n_pop, # number of populations
      n_regions = n_regions, # regions
      n_sexes = n_sexes, # sexes
      n_ages = n_ages, # ages
      n_seas = n_seas, # seasons
      seasdur = seasdur,  # fracion of time in season
      natmort = array(natmort[,,1,,,sim], dim = c(n_pop, n_regions, n_ages, n_sexes)), # natural mortality in first year
      init_F = init_F, # initial F applied (0 for unfished)
      fish_sel = array(fish_sel[,1,,,,sim], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)), # fishery selectivity in first year
      R0_r = array(R0[,,1,sim], dim = c(n_pop, n_regions)), # regional mean or virgin recruitment
      sexratio = array(sexratio[,,1,,sim], dim = c(n_pop, n_regions, n_sexes)), # sex ratio in first year
      Movement = array(Movement[,,,1,,,,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
      do_recruits_move = do_recruits_move, # whether recruits move
      ln_InitDevs = array(ln_InitDevs[,,,sim], dim = c(n_pop, n_regions, n_ages - 1)) # initial deviations
    )

    # Get initial unfished NAA
    Init_Unfished_NAA = Get_Init_NAA(
      init_age_strc = init_age_strc, # initial age structure
      init_iter = n_ages * 5, # if init_age_strc == 0, number of iterations to run
      n_pop = n_pop, # number of populations
      n_regions = n_regions, # regions
      n_sexes = n_sexes, # sexes
      n_ages = n_ages, # ages
      natmort = array(natmort[,,1,,,sim], dim = c(n_pop, n_regions, n_ages, n_sexes)), # natural mortality in first year
      init_F = rep(0, n_seas), # initial F applied (0 for unfished)
      n_seas = n_seas, # seasons
      seasdur = seasdur,  # fracion of time in season
      fish_sel = array(fish_sel[,1,,,,sim], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)), # fishery selectivity in first year
      R0_r = array(R0[,,1,sim], dim = c(n_pop, n_regions)), # regional mean or virgin recruitment
      sexratio = array(sexratio[,,1,,sim], dim = c(n_pop, n_regions, n_sexes)), # sex ratio in first year
      Movement = array(Movement[,,,1,,,,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
      do_recruits_move = do_recruits_move, # whether recruits move
      ln_InitDevs = array(ln_InitDevs[,,,sim], dim = c(n_pop, n_regions, n_ages - 1)) # initial deviations
    )

    # Input into model arrays and assign back to simulation environment (first year and first season)
    sim_env$NAA[,,1,1,,,sim] = Init_Fished_NAA
    sim_env$NAA0[,,1,1,,,sim] = Init_Unfished_NAA

  })

}

#' Title Generate Recruitment for Simulation
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation Environment
#' @keywords internal
generate_recruitment <- function(y,
                                 sim,
                                 sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {

    # Get deterministic recruitment
    tmp_det_rec <- Get_Det_Recruitment(recruitment_model = recruitment_opt,
                                       rec_dd = rec_dd,
                                       y = y,
                                       rec_lag = rec_lag,
                                       R0 = apply(R0[,,y,sim, drop = FALSE], 1, sum), # sum to get global R0
                                       Rec_Prop =       t(apply(R0[,,y,sim, drop = FALSE], c(1), function(x) x / sum(x))), # get R0 proportion
                                       h = array(h[,,y,sim], dim = c(n_pop, n_regions)),
                                       n_pop = n_pop,
                                       n_regions = n_regions,
                                       n_ages = n_ages,

                                       # Note: Using first year and female quantities to compute unfished SSB0
                                       sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,1,1,sim], dim = c(n_pop, n_regions)),
                                       WAA = array(WAA[,,1,,,1,sim], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                       MatAA = array(MatAA[,,1,,,1,sim], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                       natmort = array(natmort[,,1,,1,sim], dim = c(n_pop, n_regions, n_ages)),
                                       Movement = array(Movement[,,,1,,,1,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                       sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,1,,1,sim], dim = c(n_pop, n_regions, n_regions, n_ages)),
                                       SSB_vals = array(SSB[,,,sim], dim = c(n_pop, n_regions, n_yrs)),
                                       t_spawn = t_spawn,
                                       n_seas = n_seas,
                                       spawn_seas = spawn_seas,
                                       seasdur = seasdur,
                                       do_recruits_move = do_recruits_move,
                                       init_F = init_F, # initial F applied
                                       fish_sel = array(fish_sel[,1,,1,1,sim], dim = c(n_regions, n_ages)) # fishery selectivity in first year
    )


    # if Rec_input exists and year index is within bounds
    use_rec_input <- exists("Rec_input") && (y <= dim(Rec_input)[3])
    tmp_ln_rec_devs <- NULL

    for(p in 1:n_pop) {

      # reset deviations for each population (draws for each popn)
      if(n_pop > 1) tmp_ln_rec_devs <- NULL

      for(r in 1:n_regions) {

        # if local DD and n_pop = 1, reset deviations for each region (draws for each region, but if n_pop > 1,
        # shares deviations across regions withn a given population)
        if(n_pop == 1 && rec_dd == 0) tmp_ln_rec_devs <- NULL

        if(use_rec_input) {
          # recruitment input
          for(s in 1:n_sexes) sim_env$NAA[p,r,y,1,1,s,sim] <- Rec_input[p,r,y,sim] * sexratio[p,r,y,s,sim]
        } else {

          # index ln_sigmaR, if n_pop == 1 and local DD, use region specific rates,
          # otherwise use p (i.e., n_pop > 1 or rec_dd == 1)
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, p)

          # simulate rec deviations
          if(is.null(tmp_ln_rec_devs)) tmp_ln_rec_devs <- stats::rnorm(1, 0, exp(ln_sigmaR[2,p,sigma_idx]))

          if(R0[p,r,y,sim] != 0) {
            sim_env$ln_RecDevs[p,r,y,sim] <- tmp_ln_rec_devs
          } else sim_env$ln_RecDevs[p,r,y,sim] <- 0

          # compute rec
          for(s in 1:n_sexes) sim_env$NAA[p,r,y,1,1,s,sim] <-
            tmp_det_rec[p,r] * exp(sim_env$ln_RecDevs[p,r,y,sim] -
            exp(ln_sigmaR[2,p,sigma_idx])^2/2) * sexratio[p,r,y,s,sim]
        }

        sim_env$Rec[p,r,y,sim] <- sum(NAA[p,r,y,1,1,,sim]) # Save recruitment estimates
        sim_env$NAA0[p,r,y,1,1,,sim] = NAA[p,r,y,1,1,,sim] # populate unfished NAA

      } # end r loop
    } # end p loop
  })
}

#' Title Apply Population Dynamics (Movement, Mortality, and SSB) in Simulation
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation Environment
#' @keywords internal
apply_pop_dy <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {

    for(seas in 1:n_seas) {

      # Mortality and Ageing
      tmp_Fmort <- array(Fmort[,y,seas,,sim], dim = c(n_regions, n_fish_fleets))
      tmp_fish_sel  <- array(fish_sel[,y,,,,sim], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets))
      tmp_natmort <- array((natmort[,,y,,,sim] * seasdur[seas]), dim = c(n_pop, n_regions, n_ages, n_sexes))
      tmp_FAA   <- apply(sweep(tmp_fish_sel, c(1,4), tmp_Fmort, "*"), c(1,2,3), sum)
      sim_env$ZAA[,,y,seas,,,sim] <- sweep(tmp_natmort, c(2,3,4), tmp_FAA, "+")

      # Movement
      for(p in 1:n_pop) {
        if(n_regions > 1) {
          if(do_recruits_move == 0) { # Recruits don't move
            for(a in 2:n_ages) { # apply movement after ageing processes - start movement at age 2
              for(s in 1:n_sexes) {
                sim_env$NAA[p,,y,seas,a,s,sim] <- t(NAA[p,,y,seas,a,s,sim]) %*% Movement[p,,,y,seas,a,s,sim] # Fished
                sim_env$NAA0[p,,y,seas,a,s,sim] <- t(NAA0[p,,y,seas,a,s,sim]) %*% Movement[p,,,y,seas,a,s,sim] # Unfished
              } # end s loop
            } # end a loop
          } # end if recruits don't move

          if(do_recruits_move == 1) { # Recruits move here
            for(a in 1:n_ages) {
              for(s in 1:n_sexes) {
                sim_env$NAA[p,,y,seas,a,s,sim] <- t(NAA[p,,y,seas,a,s,sim]) %*% Movement[p,,,y,seas,a,s,sim] # Fished
                sim_env$NAA0[p,,y,seas,a,s,sim] <- t(NAA0[p,,y,seas,a,s,sim]) %*% Movement[p,,,y,seas,a,s,sim] # Unfished
              } # end s loop
            } # end a loop
          } # end if
        } # only compute if spatial
      } # end p loop

      # print(NAA[1,,y,1,1,1,sim])

      if(seas < n_seas) { # Within year seasonal mortality
        sim_env$NAA[,,y,seas + 1,1:n_ages,,sim] = NAA[,,y,seas,1:n_ages,,sim] * exp(-ZAA[,,y,seas,1:n_ages,,sim]) # fished
        sim_env$NAA0[,,y,seas + 1,1:n_ages,,sim] <- NAA0[,,y,seas,1:n_ages,,sim] * exp(-(tmp_natmort[,,1:n_ages,])) # unfished
      } else {
        # Advance into the next year, season 1
        sim_env$NAA[,,y+1,1,2:n_ages,,sim] = NAA[,,y,seas,1:(n_ages-1),,sim] * exp(-ZAA[,,y,seas,1:(n_ages-1),,sim]) # fished
        sim_env$NAA[,,y+1,1,n_ages,,sim] = NAA[,,y+1,1,n_ages,,sim] + NAA[,,y,seas,n_ages,,sim] * exp(-ZAA[,,y,seas,n_ages,,sim]) # Acuumulate plus group (fished)
        sim_env$NAA0[,,y+1,1,2:n_ages,,sim] = NAA0[,,y,seas,1:(n_ages-1),,sim] * exp(-(tmp_natmort[,,1:(n_ages - 1),])) # fished
        sim_env$NAA0[,,y+1,1,n_ages,,sim] = NAA0[,,y+1,1,n_ages,,sim] + NAA0[,,y,seas,n_ages,,sim] * exp(-(tmp_natmort[,,n_ages,])) # Acuumulate plus group (fished)
      }

      # Compute Biomass Quantities
      if(seas == spawn_seas) {

        # Get NAA for spawning
        tmp_NAA_spawn <- NAA[,,y,spawn_seas,,,sim, drop = FALSE]
        tmp_NAA0_spawn <- NAA0[,,y,spawn_seas,,,sim, drop = FALSE]

        # If we we are natal homing with 1 season
        if(n_seas == 1 && n_pop > 1) {
          # Get NAA during spawning
          for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
            tmp_NAA_spawn[p,,1,1,a,s,1] <- tmp_NAA_spawn[p,,1,1,a,s,1] %*% sgl_seas_spawning_movement[p,,,y,a,s,sim]
            tmp_NAA0_spawn[p,,1,1,a,s,1] <- tmp_NAA0_spawn[p,,1,1,a,s,1] %*% sgl_seas_spawning_movement[p,,,y,a,s,sim]
          } # end s loop
        }


        # Total Biomass
        sim_env$Total_Biom[,, y, sim] <- apply(tmp_NAA_spawn *
                                                 WAA[,, y, spawn_seas, , , sim,drop = FALSE] *
                                                 exp(-ZAA[,,y,spawn_seas,,,sim,drop = FALSE] * t_spawn), c(1,2), sum)

        # Spawning Stock Biomass
        sim_env$SSB[,, y, sim] <- apply(tmp_NAA_spawn *
                                          WAA[,, y, spawn_seas, , 1, sim,drop = FALSE] *
                                          MatAA[,, y, spawn_seas, , 1, sim,drop = FALSE] *
                                          exp(-ZAA[,, y, spawn_seas, , 1, sim,drop = FALSE] * t_spawn), c(1,2), sum)

        # Get dynamic B0
        SSB0_array <- tmp_NAA0_spawn *  WAA[,,  y, spawn_seas, , 1, sim, drop = FALSE] * MatAA[,,y, spawn_seas, , 1, sim, drop = FALSE]
        mort_spawn <- exp(-natmort[,, y, , 1, sim, drop = FALSE] * t_spawn * seasdur[spawn_seas])
        mort_spawn <- array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
        sim_env$Dynamic_SSB0[,,y,sim] <- apply(SSB0_array * mort_spawn, c(1,2), sum) # Dynamic B0

        if(n_sexes == 1) { # If single sex model, multiply SSB calculations by 0.5
          sim_env$SSB[,,y,sim] <- SSB[,,y,sim] * 0.5
          sim_env$Dynamic_SSB0[,,y,sim] <- Dynamic_SSB0[,,y,sim] * 0.5
        }

      } # if season = spawning season
    } # end seas loop
  })
}

#' Title Generate Fishery Catches, Comps, and Indices in Simulation
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation Environment
#' @keywords internal
generate_fishery_catch_comp_idx <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {
    for(seas in 1:n_seas) {
      for(r in 1:n_regions) {
        for(f in 1:n_fish_fleets) {

          for(p in 1:n_pop) {
            # Baranov's catch equation
            sim_env$CAA[p,r,y,seas,,,f,sim] <- (Fmort[r,y,seas,f,sim] * fish_sel[r,y,,,f,sim]) / ZAA[p,r,y,seas,,,sim] *
                                                NAA[p,r,y,seas,,,sim] * (1 - exp(-ZAA[p,r,y,seas,,,sim]))
            if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) for(s in 1:n_sexes) sim_env$CAL[p,r,y,seas,,s,f,sim] <- SizeAgeTrans[p,r,y,seas,,,s,sim] %*% CAA[p,r,y,seas,,s,f,sim] # Catch at length
          } # end p loop

          # Catch
          if(catch_units[f] == 0) sim_env$TrueCatch[r,y,seas,f,sim] <- sum(CAA[,r,y,seas,,,f,sim]) # abundance
          if(catch_units[f] == 1) sim_env$TrueCatch[r,y,seas,f,sim] <- sum(CAA[,r,y,seas,,,f,sim] * WAA_fish[,r,y,seas,,,f,sim]) # biomass
          sim_env$ObsCatch[r,y,seas,f,sim] <- TrueCatch[r,y,seas,f,sim] * exp(stats::rnorm(1, 0, exp(ln_sigmaC[r,y,seas,f]))) # Observed Catch w/ lognormal deviations

          # Fishery Index
          tmp_expl_abd <- sweep(NAA[,r,y,seas,,,sim, drop = F], c(5,6), fish_sel[r,y,,,f,sim, drop = F], "*") # get exploitable abundance
          tmp_expl_biom <- sweep(tmp_expl_abd, c(1,5,6), WAA_fish[,r,y,seas,,,f,sim, drop = F], "*") # get exploitable abundance
          if(fish_idx_type[f] == 0) sim_env$TrueFishIdx[r,y,seas,f,sim] <- fish_q[r,y,f,sim] * sum(tmp_expl_abd) # True Fishery Index (abundance)
          if(fish_idx_type[f] == 1) sim_env$TrueFishIdx[r,y,seas,f,sim] <- fish_q[r,y,f,sim] * sum(tmp_expl_biom) # True Fishery Index (biomass)
          sim_env$ObsFishIdx[r,y,seas,f,sim] <- TrueFishIdx[r,y,seas,f,sim] * exp(stats::rnorm(1, 0, ObsFishIdx_SE[r,y,seas,f])) # Observed Fishery index w/ lognormal deviations

          # Fishery Compositions
          if(Fmort[r,y,seas,f,sim] > 0) { # only simulate if Fishing Mortality > 0

            # Age Compositions (Dynamic ISS based on feedback fishing mortality)
            if(exists("ISS_FishAgeComps_fill") && isTRUE(ISS_FishAgeComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              sim_env$ISS_FishAgeComps[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishAgeComps, Fmort = Fmort, y = y, sim = sim, seas = seas)
            }

            # Length Compositions (Dynamic ISS based on feedback fishing mortality)
            if(exists("ISS_FishLenComps_fill") && isTRUE(ISS_FishLenComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              sim_env$ISS_FishLenComps[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishLenComps, Fmort = Fmort, y = y, sim = sim, seas = seas)
            }

            # Sample fishery ages
            sim_env$ObsFishAgeComps <- simulate_comps(r = r,
                                                      y = y,
                                                      seas = seas,
                                                      f = f,
                                                      sim = sim,
                                                      Exp = CAA,
                                                      ISS = ISS_FishAgeComps,
                                                      AgeingError = AgeingError,
                                                      comp_like = comp_fishage_like,
                                                      ln_theta = ln_FishAge_theta,
                                                      ln_theta_agg = ln_FishAge_theta_agg,
                                                      corr_pars = FishAge_corr_pars,
                                                      corr_pars_agg = FishAge_corr_pars_agg,
                                                      comp_type = FishAgeComps_Type,
                                                      n_sexes = n_sexes,
                                                      n_regions = n_regions,
                                                      n_cat = n_ages,
                                                      Obs = ObsFishAgeComps,
                                                      age_or_len = 0,
                                                      pop_specific = FALSE)

            # Sample fishery lengths
            if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) {
              sim_env$ObsFishLenComps <- simulate_comps(r = r,
                                                        y = y,
                                                        seas = seas,
                                                        f = f,
                                                        sim = sim,
                                                        Exp = CAL,
                                                        ISS = ISS_FishLenComps,
                                                        AgeingError = NULL,
                                                        comp_like = comp_fishlen_like,
                                                        ln_theta = ln_FishLen_theta,
                                                        ln_theta_agg = ln_FishLen_theta_agg,
                                                        corr_pars = FishLen_corr_pars,
                                                        corr_pars_agg = FishLen_corr_pars_agg,
                                                        comp_type = FishLenComps_Type,
                                                        n_sexes = n_sexes,
                                                        n_regions = n_regions,
                                                        n_cat = n_lens,
                                                        Obs = ObsFishLenComps,
                                                        age_or_len = 1,
                                                        pop_specific = FALSE)

            } # end if size age transition if availiable
          } # end if Fmort > 0

        } # end f loop
      } # end r loop
    } # end seas loop
  })

}


#' Title Generate Survey Comps and Indices in Simulation
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation Environment
#' @keywords internal
generate_survey_comp_idx <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {

    for(seas in 1:n_seas) {
      for(r in 1:n_regions) {
        for(sf in 1:n_srv_fleets) {

          for(p in 1:n_pop) {
            # Survey Ages Indexed (midpoint year)
            sim_env$SrvIAA[p,r,y,seas,,,sf,sim] <- NAA[p,r,y,seas,,,sim] * srv_sel[r,y,,,sf,sim] * exp(-t_srv[r,seas,sf] * ZAA[p,r,y,seas,,,sim])
            if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) for(s in 1:n_sexes) sim_env$SrvIAL[p,r,y,seas,,s,sf,sim] <- SizeAgeTrans[p,r,y,seas,,,s,sim] %*% SrvIAA[p,r,y,seas,,s,sf,sim] # Survey index at length
          } # end p loop

          # Survey Index
          if(srv_idx_type[sf] == 0) sim_env$TrueSrvIdx[r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[,r,y,seas,,,sf,sim]) # True Survey Index (abundance)
          if(srv_idx_type[sf] == 1) sim_env$TrueSrvIdx[r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[,r,y,seas,,,sf,sim] * WAA_srv[,r,y,seas,,,sf,sim]) # True Survey Index (biomass)
          sim_env$ObsSrvIdx[r,y,seas,sf,sim] <- TrueSrvIdx[r,y,seas,sf,sim] * exp(stats::rnorm(1, 0, ObsSrvIdx_SE[r,y,seas,sf])) # Observed survey index w/ lognormal deviations

          # Survey Compositions
          # Sample survey ages
          sim_env$ObsSrvAgeComps <- simulate_comps(r = r,
                                                   y = y,
                                                   f = sf,
                                                   seas = seas,
                                                   sim = sim,
                                                   Exp = SrvIAA,
                                                   ISS = ISS_SrvAgeComps,
                                                   AgeingError = AgeingError,
                                                   comp_like = comp_srvage_like,
                                                   ln_theta = ln_SrvAge_theta,
                                                   ln_theta_agg = ln_SrvAge_theta_agg,
                                                   corr_pars = SrvAge_corr_pars,
                                                   corr_pars_agg = SrvAge_corr_pars_agg,
                                                   comp_type = SrvAgeComps_Type,
                                                   n_sexes = n_sexes,
                                                   n_regions = n_regions,
                                                   n_cat = n_ages,
                                                   Obs = ObsSrvAgeComps,
                                                   age_or_len = 0,
                                                   pop_specific = FALSE)

          # Sample survey lengths
          if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) {
            sim_env$ObsSrvLenComps <- simulate_comps(r = r,
                                                     y = y,
                                                     f = sf,
                                                     seas = seas,
                                                     sim = sim,
                                                     Exp = SrvIAL,
                                                     ISS = ISS_SrvLenComps,
                                                     AgeingError = NULL,
                                                     comp_like = comp_srvlen_like,
                                                     ln_theta = ln_SrvLen_theta,
                                                     ln_theta_agg = ln_SrvLen_theta_agg,
                                                     corr_pars = SrvLen_corr_pars,
                                                     corr_pars_agg = SrvLen_corr_pars_agg,
                                                     comp_type = SrvLenComps_Type,
                                                     n_sexes = n_sexes,
                                                     n_regions = n_regions,
                                                     n_cat = n_lens,
                                                     Obs = ObsSrvLenComps,
                                                     age_or_len = 1,
                                                     pop_specific = FALSE)

          } # end if size age transition if availiable

        } # end sf loop
      } # end r loop
    } # end seas loop

  })
}

#' Title Generate Conventional Tag Releases
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation Environment
#' @keywords internal
release_conv_tags <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {
    for(seas in 1:n_seas) {
      for(r in 1:n_regions) {

        # Get indices for tag cohorts in the current year and region
        tag_rel <- which(tag_release_indicator[,1] == r & tag_release_indicator[,2] == y & tag_release_indicator[,3] == seas) # Get tag cohort (release event)

        # Release Tags if any events
        if(length(tag_rel) != 0) {

          # Tag Indexing
          tr <- tag_release_indicator[tag_rel,1] # tag release region
          ty <- tag_release_indicator[tag_rel,2] # tag release year
          tseas <- tag_release_indicator[tag_rel,3] # tag release season
          tplat <- tag_release_platform[tag_rel, 1] # get tagging platform information
          tplat_f <- as.numeric(tag_release_platform[tag_rel, 2]) # get tagging fleet

          # distribute tags by survey
          if(tplat[1] == 'survey') {
            if(!exists("n_tags_rel_input")) {
              n_tags_rel <- round(sum(SrvIAA[,tr,ty,tseas,,,tplat_f,sim]) / sum(SrvIAA[,,ty,tseas,,,tplat_f,sim]) * n_tags) # get region specific tags
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel] # use input tags by cohort if availiable
            }
            tmp_props <- SrvIAA[, tr, ty, tseas, , , tplat_f, sim] / sum(SrvIAA[, tr, ty, tseas, , , tplat_f, sim]) # get proportions by population, age, and sex
          }

          # distribute tags by fishery
          if(tplat[1] == 'fishery') {
            if(!exists("n_tags_rel_input")) {
              n_tags_rel <- round(sum(CAA[,tr,ty,tseas,,,tplat_f,sim]) / sum(CAA[,,ty,tseas,,,tplat_f,sim]) * n_tags) # get region specific tags
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel] # use input tags by cohort if availiable
            }
            tmp_props <- CAA[, tr, ty, tseas, , , tplat_f, sim] / sum(CAA[, tr, ty, tseas, , , tplat_f, sim]) # get proportions by population, age, and sex
          }

          # distribute tags by population
          if(tplat[1] == 'population') {
            if(!exists("n_tags_rel_input")) {
              n_tags_rel <- round(sum(NAA[,tr,ty,tseas,,,sim]) / sum(NAA[,,ty,tseas,,,sim]) * n_tags) # get region specific tags
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel] # use input tags by cohort if availiable
            }
            tmp_props <- NAA[, tr, ty, tseas, , , sim] / sum(NAA[, tr, ty, tseas, , , sim]) # get proportions by population, age, and sex
          }

          # multiply by tags in each region and distributa cross ages and sexes
          sim_env$conv_tagged_fish[tag_rel, , , , sim] <- array(round(tmp_props * n_tags_rel), dim = c(n_pop, n_ages, n_sexes))

        } # end if no tag releases

      } # end r loop
    } # end seas loop
  })
}

#' Title Generate Conventional Tag Recaptures from Fisheries
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation Environment
#' @keywords internal
generate_fishery_conv_tags_recap <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env,{

      for(rseas in 1:n_seas) {
        for(tc in 1:n_tag_rel_events) {

          # get indexing
          tr <- tag_release_indicator[tc,1] # tag release region
          ty <- tag_release_indicator[tc,2] # tag release year
          tseas <- tag_release_indicator[tc,3] # tag release seasons

          # Skipping stuff if hasn't occurred yet, or if max liberty
          if(y < ty || (y == ty && rseas < tseas)) next
          ry <- y - ty + 1 # get tag liberty
          if(ry > max_liberty) next # skip if max liberty

          # get fishing mortality
          tmp_F <- array(Fmort[, y, rseas, , sim] , dim = c(n_regions, 1, n_ages, n_sexes, n_fish_fleets, 1))
          tmp_FAA <- tmp_F * fish_sel[, y, , , , sim, drop = FALSE]

          # get total mortality
          tmp_natmort <- array(natmort[,,y,,,sim], dim = c(n_pop, n_regions, 1, n_ages, n_sexes, 1))
          tmp_ZAA <- sweep((tmp_natmort * seasdur[rseas]), c(2,3,4,5), apply(tmp_FAA, 1:4, sum), "+") + (exp(ln_conv_tag_shed) * seasdur[rseas])

          # Discount with tagging time (t_tagging) if it doesn't happen at the start of the season / year
          if(ry == 1 && rseas == tseas) {
            if(t_tagging != 1) tmp_ZAA <- tmp_ZAA * t_tagging
            # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
            sim_env$conv_tag_fish_avail[1, rseas, tc, , tr, , , sim] <- array(conv_tagged_fish[tc, , , , sim] * exp(-exp(ln_init_conv_tag_mort)), dim = c(n_pop, n_ages, n_sexes))
          }

          # get temporary survival value
          tmp_SAA <- exp(-tmp_ZAA)

          # Move tagged fish around (skip only in first release year + tagging season when tagging occurs mid-season)
          if(t_tagging == 1 || ry != 1 || rseas != tseas) {
            for(p in 1:n_pop) {
              # Movement of tag cohorts
              if(do_recruits_move == 0) {
                for(a in 2:n_ages) for(s in 1:n_sexes) {
                  sim_env$conv_tag_fish_avail[ry, rseas, tc, p, , a, s, sim] <-
                    t(conv_tag_fish_avail[ry, rseas, tc, p, , a, s, sim]) %*%
                    Movement[p, , , y, rseas, a, s, sim]
                }
              } else { # if recruits move
                for(a in 1:n_ages) for(s in 1:n_sexes) {
                  sim_env$conv_tag_fish_avail[ry, rseas, tc, p, , a, s, sim] <-
                    t(conv_tag_fish_avail[ry, rseas, tc, p, , a, s, sim]) %*%
                    Movement[p, , , y, rseas, a, s, sim]
                } # end s loop
              } # end else
            } # end p loop
          } # end if

          # Apply mortality and ageing to tagged fish
          if(rseas < n_seas) {

            # Season mortality within a given year, advance to next season same year/age
            sim_env$conv_tag_fish_avail[ry, rseas + 1, tc, , , , , sim] <-
              conv_tag_fish_avail[ry, rseas, tc, , , , , sim] *
              tmp_SAA[,,1,,,1]

          } else {

            # End of year mortality and age advancement (end of season)
            sim_env$conv_tag_fish_avail[ry + 1, 1, tc, , , 2:n_ages, , sim] <-
              conv_tag_fish_avail[ry, n_seas, tc, , , 1:(n_ages-1), , sim] *
              tmp_SAA[,,1,1:(n_ages - 1),,1]

            # Accumulate plus group
            sim_env$conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, , sim] <-
              conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, , sim] +
              conv_tag_fish_avail[ry, n_seas, tc, , , n_ages, , sim] *
              tmp_SAA[,,1,n_ages,,1]
          }

          # # Apply Baranov's to get predicted recaptures
          for(f in 1:n_fish_fleets) {
            for(p in 1:n_pop) {
              sim_env$pred_conv_tag_fish_recap[ry,rseas,tc,p,,,,f,sim] <- conv_tag_fish_reporting[,y,f,sim] *
                (tmp_FAA[,1,,,f,1] / tmp_ZAA[p,,1,,,1]) *
                conv_tag_fish_avail[ry,rseas,tc,p,,,,sim] *
                (1 - tmp_SAA[p,,1,,,1])
            } # end p loop
          } # end f loop

          # Simulate Tag Recoveries
          sim_env$obs_conv_tag_fish_recap <- simulate_conv_tag_fish_recaptures(
            conv_fish_tag_like = conv_fish_tag_like,
            tag_recaptures_attr = conv_fish_tag_attr,
            conv_tagged_fish = conv_tagged_fish,
            pred_conv_tag_fish_recap = pred_conv_tag_fish_recap,
            obs_conv_tag_fish_recap = obs_conv_tag_fish_recap,
            ln_conv_fish_tag_theta = ln_conv_fish_tag_theta,
            ry = ry,
            rseas = rseas,
            tc = tc,
            sim = sim,
            n_pop = n_pop,
            n_regions = n_regions,
            n_ages = n_ages,
            n_sexes = n_sexes,
            n_fish_fleets = n_fish_fleets
          )

        } # end tc loop
      } # end rseas loop
  })
}


#' Run Annual Cycle in Simulation Environment
#'
#' @param y Year index
#' @param sim Simulation index
#' @param sim_env Simulation environment will all the necessary elements to run the annual cycle
#' @export run_annual_cycle
#' @importFrom stats rnorm rmultinom
#' @family Simulation Setup
run_annual_cycle <- function(y,
                             sim,
                             sim_env) {

  if(y == 1) {
    generate_initial_age_structure(y = 1, sim, sim_env) # Initialize age structure
    generate_recruitment(y = 1, sim, sim_env) # Get recruitment in the first year
  }

  apply_pop_dy(y, sim, sim_env) # Apply population dynamics (movement, mortality, and biomass calculations)
  generate_fishery_catch_comp_idx(y, sim, sim_env) # Get Fishery Catches, Compositions, and Indices
  generate_survey_comp_idx(y, sim, sim_env) # Get Fishery Catches, Compositions, and Indices
  release_conv_tags(y, sim, sim_env) # Release conventional tags

  if(sim_env$use_conv_fish_tagging == 1) generate_fishery_conv_tags_recap(y, sim, sim_env) # Generate fishery conventional tag recaptures
  if(y < sim_env$n_yrs) generate_recruitment(y = y + 1, sim, sim_env) # Get recruitment in the following year


  return(invisible(NULL))

}

#' Simulates a static spatial, sex, and age-structured population (no feedback loop)
#'
#' @param output_path path to output simulation objects
#' @param sim_list Simulation list objects
#'
#' @returns a list object with a bunch of simulated values and outputs
#' @export Simulate_Pop_Static
#' @family Simulation Setup
Simulate_Pop_Static <- function(sim_list,
                                output_path = NULL) {

  # Setup simulation environment
  sim_env <- Setup_sim_env(sim_list)

  # Start Simulation
  for (sim in 1:sim_env$n_sims) {
    for (y in 1:sim_env$n_yrs) {
      # Run annual cycle here
      run_annual_cycle(y = y, sim = sim, sim_env = sim_env)
    } # end y loop
  } # end sim loop

  # Output simulation outputs as a list
  sim_out <- list(init_F = sim_env$init_F,
                  Fmort = sim_env$Fmort,
                  ln_sigmaC = sim_env$ln_sigmaC,
                  fish_sel = sim_env$fish_sel,
                  fish_q = sim_env$fish_q,
                  ln_RecDevs = sim_env$ln_RecDevs,
                  ln_InitDevs = sim_env$ln_InitDevs,
                  natmort = sim_env$natmort,
                  ZAA = sim_env$ZAA,
                  sexratio = sim_env$sexratio,
                  R0 = sim_env$R0,
                  Rec = sim_env$Rec,
                  WAA = sim_env$WAA,
                  WAA_fish = sim_env$WAA_fish,
                  WAA_srv = sim_env$WAA_srv,
                  MatAA = sim_env$MatAA,
                  h = sim_env$h,
                  do_recruits_move = sim_env$do_recruits_move,
                  ln_sigmaR = sim_env$ln_sigmaR,
                  Movement = sim_env$Movement,
                  NAA = sim_env$NAA,
                  NAA0 = sim_env$NAA0,
                  Dynamic_SSB0 = sim_env$Dynamic_SSB0,
                  SSB = sim_env$SSB,
                  t_spawn = sim_env$t_spawn,
                  Total_Biom = sim_env$Total_Biom,
                  TrueCatch = sim_env$TrueCatch,
                  ObsCatch = sim_env$ObsCatch,
                  ObsFishIdx = sim_env$ObsFishIdx,
                  TrueFishIdx = sim_env$TrueFishIdx,
                  ObsFishIdx_SE = sim_env$ObsFishIdx_SE,
                  CAA = sim_env$CAA,
                  CAL = sim_env$CAL,
                  ObsFishAgeComps = sim_env$ObsFishAgeComps,
                  ObsFishLenComps = sim_env$ObsFishLenComps,
                  ObsSrvIdx = sim_env$ObsSrvIdx,
                  TrueSrvIdx = sim_env$TrueSrvIdx,
                  ObsSrvIdx_SE = sim_env$ObsSrvIdx_SE,
                  SrvIAA = sim_env$SrvIAA,
                  SrvIAL = sim_env$SrvIAL,
                  srv_sel = sim_env$srv_sel,
                  srv_q = sim_env$srv_q,
                  ObsSrvAgeComps = sim_env$ObsSrvAgeComps,
                  ObsSrvLenComps = sim_env$ObsSrvLenComps,
                  tag_release_indicator = as.matrix(sim_env$tag_release_indicator),
                  conv_tag_fish_reporting = sim_env$conv_tag_fish_reporting,
                  conv_tagged_fish = sim_env$conv_tagged_fish,
                  ln_init_conv_tag_mort = sim_env$ln_init_conv_tag_mort,
                  ln_conv_tag_shed = sim_env$ln_conv_tag_shed,
                  conv_tag_fish_avail = sim_env$conv_tag_fish_avail,
                  use_conv_fish_tagging = sim_env$use_conv_fish_tagging,
                  pred_conv_tag_fish_recap = sim_env$pred_conv_tag_fish_recap,
                  obs_conv_tag_fish_recap = sim_env$obs_conv_tag_fish_recap,
                  SizeAgeTrans = if(!is.null(sim_env$SizeAgeTrans)) sim_env$SizeAgeTrans else NULL,
                  AgeingError = sim_env$AgeingError,
                  ISS_FishAgeComps = sim_env$ISS_FishAgeComps,
                  ISS_FishLenComps = sim_env$ISS_FishLenComps,
                  ISS_SrvAgeComps = sim_env$ISS_SrvAgeComps,
                  ISS_SrvLenComps = sim_env$ISS_SrvLenComps,
                  n_sims = sim_env$n_sims,
                  n_regions = sim_env$n_regions,
                  n_pop = sim_env$n_pop,
                  n_years = sim_env$n_yrs,  # duplicated to ensure backwards compatbility
                  n_yrs = sim_env$n_yrs, # duplicated to ensure backwards compatbility
                  n_ages = sim_env$n_ages,
                  n_seas = sim_env$n_seas,
                  seasdur = sim_env$seasdur,
                  spawn_seas = sim_env$spawn_seas,
                  n_lens = if(!is.null(sim_env$n_lens)) sim_env$n_lens else NULL,
                  n_sexes = sim_env$n_sexes,
                  n_fish_fleets = sim_env$n_fish_fleets,
                  n_srv_fleets = sim_env$n_srv_fleets
                  )

  # save RDS file
  if(!is.null(output_path)) saveRDS(sim_out, file = output_path)

  return(sim_out)

} # end function


#' Conduct a Simulation Self Test
#'
#' This function runs a self test of the fitted RTMB model by simulating new
#' datasets under the fitted parameters, refitting the model, and comparing
#' estimated outputs to the true values used for simulation. It can be run
#' sequentially or in parallel.
#'
#' @param data A list containing model data from an RTMB object.
#' @param parameters A list of fitted parameter values from an RTMB object.
#' @param mapping A list specifying parameter mappings from an RTMB object.
#' @param random Character vector specifying random effects.
#' @param rep A list of report values from an RTMB object (`$rep`).
#' @param sd_rep An `sdreport` object from RTMB summarizing parameter uncertainty.
#' @param n_sims Integer. Number of simulation replicates to run.
#' @param newton_loops Integer. Number of Newton loops used in model fitting (default: `3`).
#' @param do_sdrep Logical. If `TRUE`, compute `sdreport` for each fitted replicate (default: `FALSE`).
#' @param do_par Logical. If `TRUE`, run simulations in parallel (default: `FALSE`).
#' @param n_cores Integer. Number of cores to use for parallelization (default: `NULL` = detect automatically).
#' @param output_path Optional file path. If provided, the simulated datasets are written to this location.
#' @param what Character vector. Names of report elements in `rep` to extract and store for each replicate.
#' @family Simulation Setup
#' @return
#' A list with elements corresponding to the requested `what` values, each containing
#' an array of simulation results across replicates. If `do_sdrep = TRUE`, an additional
#' element `"sd_rep"` is included with the list of `sdreport` objects (or `NA` if a replicate fails).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run a simple self test with 10 simulations, extracting SSB
#' res <- simulation_self_test(
#'   data = model$data,
#'   parameters = model$parameters,
#'   mapping = model$mapping,
#'   random = model$random,
#'   rep = model$rep,
#'   sd_rep = model$sd_rep,
#'   n_sims = 10,
#'   what = "SSB"
#' )
#'
#' str(res$SSB) # look at simulated SSB arrays
#' }
simulation_self_test <- function(data,
                                 parameters,
                                 mapping,
                                 random,
                                 rep,
                                 sd_rep,
                                 n_sims,
                                 newton_loops = 3,
                                 do_sdrep = FALSE,
                                 do_par = FALSE,
                                 n_cores = NULL,
                                 output_path = NULL,
                                 what = c('SSB', 'Rec')
                                 ) {

  missing_names <- setdiff(what, names(rep))
  if(length(missing_names) > 0)  stop(paste("The following elements in 'what' are not found in rep:",  paste(missing_names, collapse = ", ")))
  optim_parameters_list <- get_optim_param_list(parameters, mapping, sd_rep, random) # get optimized parameters in original list format

  # Modify any data weights that are NA to 0
  if(any(is.na(data$Wt_Catch))) data$Wt_Catch[is.na(data$Wt_Catch)] <- 0
  if(any(is.na(data$Wt_FishAgeComps))) data$Wt_FishAgeComps[is.na(data$Wt_FishAgeComps)] <- 0
  if(any(is.na(data$Wt_FishLenComps))) data$Wt_FishLenComps[is.na(data$Wt_FishLenComps)] <- 0
  if(any(is.na(data$Wt_FishIdx))) data$Wt_FishIdx[is.na(data$Wt_FishIdx)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps))) data$Wt_SrvAgeComps[is.na(data$Wt_SrvAgeComps)] <- 0
  if(any(is.na(data$Wt_SrvLenComps))) data$Wt_SrvLenComps[is.na(data$Wt_SrvLenComps)] <- 0
  if(any(is.na(data$Wt_SrvIdx))) data$Wt_SrvIdx[is.na(data$Wt_SrvIdx)] <- 0
  if(any(is.na(data$Wt_Tagging))) data$Wt_Tagging[is.na(data$Wt_Tagging)] <- 0

  # Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims, # number of simulations
                            n_yrs = length(data$years), # number of years
                            n_regions = data$n_regions,  # number of regions
                            n_ages = length(data$ages), # number of ages
                            n_obs_ages = dim(data$ObsFishAgeComps)[4], # number of observed ages
                            n_lens = length(data$lens), # number of lengths
                            n_sexes = data$n_sexes, # number of sexes
                            n_fish_fleets = data$n_fish_fleets, # number of fishery fleets
                            n_srv_fleets = data$n_srv_fleets, # number of survey fleets
                            # Seasonal stuff
                            n_seas = data$n_seas,
                            seasdur = data$seasdur
  )

  # Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list)

  # Setup Fishing Processes -------------------------------------------------
  ln_sigmaC <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch)) ln_sigmaC[r,,,f] <- log(exp(optim_parameters_list$ln_sigmaC[r,,,f]) / sqrt(data$Wt_Catch[r,,,f]))
    else ln_sigmaC[r,,,f] <- log(exp(optim_parameters_list$ln_sigmaC[r,,,f]) / sqrt(data$Wt_Catch))
  }

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(sim_list = sim_list, # update simulate list
                                ln_sigmaC = ln_sigmaC, # sigmaC
                                catch_units = data$catch_units, # catch units
                                Fmort_input = replicate(n = sim_list$n_sims, rep$Fmort[,1:length(data$years),,,drop = FALSE]), # use fishing mortality from report
                                fish_sel_input = replicate(n = sim_list$n_sims, rep$fish_sel[,1:length(data$years),,,,drop = FALSE]), # use fishery selectivity from report
                                fish_q_input = replicate(n = sim_list$n_sims, rep$fish_q[,1:length(data$years),,drop = FALSE]), # use fishery catchability from report
                                ObsFishIdx_SE = data$ObsFishIdx_SE / sqrt(data$Wt_FishIdx), # fishery index uncertainty
                                fish_idx_type = data$fish_idx_type, # fishery index type
                                init_F_val = rep$init_F,

                                # fishery age composition specifications
                                comp_fishage_like = data$FishAgeComps_LikeType, # age comps likelihood
                                FishAgeComps_Type = data$FishAgeComps_Type, # age comps structure
                                ISS_FishAgeComps = replicate(sim_list$n_sims, data$ISS_FishAgeComps[,,,,,drop = F] * data$Wt_FishAgeComps), # input sample size
                                ln_FishAge_theta = optim_parameters_list$ln_FishAge_theta[,,,drop = F] , # overdispersion
                                ln_FishAge_theta_agg = optim_parameters_list$ln_FishAge_theta_agg, # aggregated overdispersion for fishery age comps
                                FishAge_corr_pars_agg = optim_parameters_list$FishAge_corr_pars_agg, # correaltion parameters for aggregated fishery age comps
                                FishAge_corr_pars = optim_parameters_list$FishAge_corr_pars[,,,,drop = F], # correlation parameters for fishery age comps

                                # fishery length composition specifications
                                comp_fishlen_like = data$FishLenComps_LikeType, # length comp likelihood
                                FishLenComps_Type = data$FishLenComps_Type, # length comps structure
                                ISS_FishLenComps = replicate(sim_list$n_sims, data$ISS_FishLenComps[,,,,,drop = F] * data$Wt_FishLenComps), # input sample size
                                ln_FishLen_theta = optim_parameters_list$ln_FishLen_theta[,,,drop = F], # overdispersion
                                ln_FishLen_theta_agg = optim_parameters_list$ln_FishLen_theta_agg, # aggregated overdispersion
                                FishLen_corr_pars_agg = optim_parameters_list$FishLen_corr_pars_agg, # correaltion parameters for aggregated fish len comps
                                FishLen_corr_pars = optim_parameters_list$FishLen_corr_pars[,,,,drop = F] # correlation parameters for fish len comps
  )


  # Setup Survey Processes --------------------------------------------------
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(n = sim_list$n_sims, rep$srv_sel[,1:length(data$years),,,,drop = FALSE]),
    srv_q_input = replicate(n = sim_list$n_sims, rep$srv_q[,1:length(data$years),,drop = FALSE]),
    ObsSrvIdx_SE = data$ObsSrvIdx_SE / sqrt(data$Wt_SrvIdx), # survey observation error
    srv_idx_type = data$srv_idx_type,
    t_srv = data$t_srv,

    # survey age composition specifications
    comp_srvage_like = data$SrvAgeComps_LikeType, # age comps likelihood
    SrvAgeComps_Type = data$SrvAgeComps_Type, # age comps structure
    ISS_SrvAgeComps = replicate(sim_list$n_sims, data$ISS_SrvAgeComps[,,,,,drop = F] * data$Wt_SrvAgeComps), # input sample size
    ln_SrvAge_theta = optim_parameters_list$ln_SrvAge_theta[,,,drop = F] , # overdispersion
    ln_SrvAge_theta_agg = optim_parameters_list$ln_SrvAge_theta_agg, # aggregated overdispersion for survey age comps
    SrvAge_corr_pars_agg = optim_parameters_list$SrvAge_corr_pars_agg, # correaltion parameters for aggregated survey age comps
    SrvAge_corr_pars = optim_parameters_list$SrvAge_corr_pars[,,,,drop = F], # correlation parameters for survey age comps

    # survey length composition specifications
    comp_srvlen_like = data$SrvLenComps_LikeType, # length comp likelihood
    SrvLenComps_Type = data$SrvLenComps_Type, # length comps structure
    ISS_SrvLenComps = replicate(sim_list$n_sims, data$ISS_SrvLenComps[,,,,,drop = F] * data$Wt_SrvLenComps), # input sample size
    ln_SrvLen_theta = optim_parameters_list$ln_SrvLen_theta[,,,drop = F], # overdispersion
    ln_SrvLen_theta_agg = optim_parameters_list$ln_SrvLen_theta_agg, # aggregated overdispersion
    SrvLen_corr_pars_agg = optim_parameters_list$SrvLen_corr_pars_agg, # correaltion parameters for aggregated srv len comps
    SrvLen_corr_pars = optim_parameters_list$SrvLen_corr_pars[,,,,drop = F] # correlation parameters for srv len comps
  )

  # Setup Biological Dynamics -----------------------------------------------
  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list, # simualtion list
    natmort_input = replicate(n = sim_list$n_sims, rep$natmort[,1:length(data$years),,,drop = FALSE]), # natuyral mortality
    WAA_input = replicate(n = sim_list$n_sims, data$WAA[,1:length(data$years),,,,drop = FALSE]), # weight at age
    WAA_fish_input = replicate(n = sim_list$n_sims, data$WAA_fish[,1:length(data$years),,,,,drop = FALSE]), # fishery weight at age
    WAA_srv_input = replicate(n = sim_list$n_sims, data$WAA_srv[,1:length(data$years),,,,,drop = FALSE]), # survey weight at age
    MatAA_input = replicate(n = sim_list$n_sims, data$MatAA[,1:length(data$years),,,,drop = FALSE]), # maturity at age
    AgeingError_input = replicate(n = sim_list$n_sims, data$AgeingError[1:length(data$years),,,drop = FALSE]), # ageing error
    SizeAgeTrans_input = replicate(n = sim_list$n_sims, data$SizeAgeTrans[,1:length(data$years),,,,,drop = FALSE]) # size age transition matrix
  )

  # Movement
  sim_list$Movement <- replicate(n = sim_list$n_sims, rep$Movement[,,1:length(data$years),,,,drop = FALSE])

  # Setup Recruitment Processes ---------------------------------------------
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    spawn_seas = data$spawn_seas, # spawning season
    do_recruits_move = data$do_recruits_move, # whether recruits move
    t_spawn = data$t_spawn, # spawn timing
    init_age_strc = data$init_age_strc, # initilaizing age structure
    h_input = replicate(n = sim_list$n_sims, array(rep$h_trans, dim = c(sim_list$n_regions, sim_list$n_yrs))), # steepness
    R0_input = replicate(n = sim_list$n_sims, expr = array(rep$R0 * rep$Rec_trans_prop, dim = c(sim_list$n_regions, sim_list$n_yrs))), # R0
    sexratio_input = replicate(n = sim_list$n_sims, expr = rep$sexratio[,1:length(data$years),,drop = FALSE]), # sex ratio
    ln_sigmaR = optim_parameters_list$ln_sigmaR / sqrt(data$Wt_Rec), # ln_sigmaR
    Rec_input = replicate(n = sim_list$n_sims, expr = rep$Rec[,1:length(data$years),drop = FALSE]), # recruitment time series
    ln_InitDevs_input = replicate(sim_list$n_sims, optim_parameters_list$ln_InitDevs) # init devs
  )

  # Setup Tagging -----------------------------------------------------------
  if(!is.na(sum(data$conv_tagged_fish))) n_tags_rel_input <- apply(data$conv_tagged_fish, 1, sum) else n_tags_rel_input <- NA
  if(exists("tag_release_indicator", data)) tag_release_indicator <- data$tag_release_indicator  else tag_release_indicator <- NA
  conv_tag_fish_reporting_input <- if(!is.null(rep$conv_tag_fish_reporting)) replicate(n = sim_list$n_sims, rep$conv_tag_fish_reporting) else NULL

  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list, # simulation list
    max_liberty = data$max_tag_liberty, # maximum tag liberty
    t_tagging = data$t_tagging, # time of tagging
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,  # number of tags to release per event
    tag_release_indicator = tag_release_indicator,  # tag release indicator
    ln_init_conv_tag_mort = optim_parameters_list$ln_init_conv_tag_mort,  # inital tagging mortality
    ln_conv_tag_shed = optim_parameters_list$ln_conv_tag_shed, # chronic tag shedding
    conv_tag_fish_reporting_input = conv_tag_fish_reporting_input, # tag reporting rates
    use_conv_fish_tagging = data$use_conv_fish_tagging, # whether or not tagging is used / simulated
    tag_selex = data$tag_selex, # tag selectivity type
    tag_natmort = data$tag_natmort, # tag natural mortality type
    conv_fish_tag_like = data$conv_fish_tag_likeType, # tag likelihood
    ln_conv_fish_tag_theta = parameters$ln_conv_fish_tag_theta # tag overdispersion
  )


  # Run Simulation ----------------------------------------------------------

  # storage
  store_res_list <- vector("list", length(what) + 1) # get list
  names(store_res_list) <- c(what, "sd_rep") # name list
  for(j in 1:length(what)) store_res_list[[j]] <- vector("list", n_sims) # stick in n_sims lists into storage

  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = output_path) # get simulated datasets

  if(do_par == FALSE) {

    for(i in 1:n_sims) {

      tryCatch({

        # set up data stuff
        tmp_data <- data # set up temporary data list
        tmp_pars <- parameters # set up temporary parameter list
        tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx)) # new fish index
        tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx)) # new srv index
        tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch)) # new catch
        tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps)) # new fishery ages
        tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps)) # new srv ages
        tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps)) # new fishery lens
        tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps)) # new survey lens

        # setup tagging data stuff if tagging is done
        if(tmp_data$use_conv_fish_tagging == 1) {
          tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,i], dim = dim(tmp_data$conv_tagged_fish)) # new tagged fish
          tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap)) # new tag recaps
          tmp_data$tag_release_indicator <- sim_obj$tag_release_indicator # release indicator
        }

        # reset weights
        tmp_data$Wt_Rec <- 1
        tmp_data$Wt_Tagging <- 1
        tmp_data$Wt_Catch[] <- 1
        tmp_data$Wt_FishAgeComps[] <- 1
        tmp_data$Wt_FishIdx <- 1
        tmp_data$Wt_FishLenComps[] <- 1
        tmp_data$Wt_SrvAgeComps[] <- 1
        tmp_data$Wt_SrvIdx <- 1
        tmp_data$Wt_SrvLenComps[] <- 1

        # input simulated uncertainty
        tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
        tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
        tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
        tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,,i]
        tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,,i]
        tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,,i]
        tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,,i]

        # Fit model
        obj <- fit_model(
          data = tmp_data,
          parameters = tmp_pars,
          mapping = mapping,
          random = random,
          newton_loops = newton_loops,
          silent = TRUE
        )

        # Populate results into store list
        for(j in 1:length(what)) store_res_list[[j]][[i]] <- obj$rep[[what[j]]]

        if(do_sdrep == TRUE) {
          tryCatch({
            obj$sd_rep <- RTMB::sdreport(obj)
            store_res_list[[length(what) + 1]][[i]] <- obj$sd_rep # input sd report
          }, error = function(e) {
            store_res_list[[length(what) + 1]][[i]] <- NA
          })
        }

      }, error = function(e) {
        # Skip failed simulations
        for(j in 1:length(what)) store_res_list[[j]][[i]] <- NA
        if(do_sdrep == TRUE) store_res_list[[length(what) + 1]][[i]] <- NA
      })

    } # end i loop

    # Convert result lists to array
    for(j in 1:length(what)) store_res_list[[j]] <- simplify2array(store_res_list[[j]])

  } # not doing parallelization

  if(do_par == TRUE) {

    # initialize cores
    if(is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    options(future.globals.maxSize = 5e3 * 1024^2)  # increase parrlalel size
    future::plan(future::multisession, workers = n_cores) # set up cores
    progressr::with_progress({
      p <- progressr::progressor(along = 1:n_sims) # progress bar

      sim_results <- future.apply::future_lapply(1:n_sims, function(i) {

        tryCatch({

          # set up data stuff
          tmp_data <- data # set up temporary data list
          tmp_pars <- parameters # set up temporary parameter list
          tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx)) # new fish index
          tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx)) # new srv index
          tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch)) # new catch
          tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps)) # new fishery ages
          tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps)) # new srv ages
          tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps)) # new fishery lens
          tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps)) # new survey lens

          # setup tagging data stuff if tagging is done
          if(tmp_data$use_conv_fish_tagging == 1) {
            tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,i], dim = dim(tmp_data$conv_tagged_fish)) # new tagged fish
            tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap)) # new tag recaps
            tmp_data$tag_release_indicator <- sim_obj$tag_release_indicator # release indicator
          }

          # reset weights
          tmp_data$Wt_Rec <- 1
          tmp_data$Wt_Tagging <- 1
          tmp_data$Wt_Catch[] <- 1
          tmp_data$Wt_FishAgeComps[] <- 1
          tmp_data$Wt_FishIdx <- 1
          tmp_data$Wt_FishLenComps[] <- 1
          tmp_data$Wt_SrvAgeComps[] <- 1
          tmp_data$Wt_SrvIdx <- 1
          tmp_data$Wt_SrvLenComps[] <- 1

          # input simulated uncertainty
          tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
          tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
          tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
          tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,i]
          tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,i]
          tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,i]
          tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,i]

          # Fit model
          obj <- fit_model(
            data = tmp_data,
            parameters = tmp_pars,
            mapping = mapping,
            random = random,
            newton_loops = newton_loops,
            silent = TRUE
          )

          # Extract what we need and return
          result <- list()
          for(j in 1:length(what)) result[[what[j]]] <- obj$rep[[what[j]]]

          if(do_sdrep == TRUE) {
            tryCatch({
              obj$sd_rep <- RTMB::sdreport(obj) # get sdreport
              result[[length(what) + 1]] <- obj$sd_rep # input sd report
            }, error = function(e) {
              result[[length(what) + 1]] <- NA
            })
          }

          p() # update progress
          return(result)

        }, error = function(e) {
          # Skip failed simulations
          result <- list()
          for(j in 1:length(what)) result[[what[j]]] <- NA
          if(do_sdrep == TRUE) result[[length(what) + 1]] <- NA

          p() # update progress
          return(result)
        })

      }, future.seed = TRUE)

      future::plan(future::sequential)  # Reset
    })

    # Populate results from parallel run
    for(i in 1:n_sims) for(j in 1:length(what)) store_res_list[[j]][[i]] <- sim_results[[i]][[what[j]]]
    if(do_sdrep == TRUE) for(i in 1:n_sims) store_res_list[[length(what) + 1]][[i]] <- sim_results[[i]][[length(what) + 1]]
    for(j in 1:length(what)) store_res_list[[j]] <- simplify2array(store_res_list[[j]])  # Convert lists to array
  }

  return(store_res_list)

}

#' Extract simulation data into SPoRC format
#'
#' This function subsets and reshapes biological, tagging, fishery, and survey
#' data from a simulation environment for use in SPoRC analyses.
#'
#' @param sim_env A simulation environment / object (list or environment) containing
#'   arrays of biological quantities, tagging information, fishery data,
#'   and survey data.
#' @param y Integer. Number of years to retain (subset from `1:y`).
#' @param sim Integer. Simulation replicate index to extract.
#' @family Simulation Setup
#' @return A named list with the following elements:
#' \describe{
#'   \item{WAA}{Weight-at-age array [region × year × seas x age × sex].}
#'   \item{MatAA}{Maturity-at-age array [region × year × seas x age × sex].}
#'   \item{SizeAgeTrans}{Size–age transition array [region × year × seas x length × age × sex].}
#'   \item{AgeingError}{Ageing error matrix [year × age × error].}
#'   \item{tag_release_indicator}{Tag release indicators (or `NULL` if tagging not used).}
#'   \item{obs_conv_tag_fish_recap}{Observed tag recapture array (or `NULL`).}
#'   \item{conv_tagged_fish}{Tagged fish counts (or `NULL`).}
#'   \item{n_tag_cohorts}{Number of tag release cohorts (or `NULL`).}
#'   \item{ObsCatch}{Observed fishery catch array [region × year × seas x fleet].}
#'   \item{ln_sigmaC}{Log Fishery Catch SD [region × year × seas x fleet].}
#'   \item{UseCatch}{Binary indicator array for catch data availability.}
#'   \item{ObsFishIdx}{Observed fishery index array [region × year × seas x fleet].}
#'   \item{ObsFishIdx_SE}{Standard error for fishery index array.}
#'   \item{UseFishIdx}{Binary indicator array for fishery indices.}
#'   \item{ObsFishAgeComps}{Observed fishery age composition array.}
#'   \item{ObsFishLenComps}{Observed fishery length composition array.}
#'   \item{ISS_FishAgeComps}{Implied sample sizes for fishery age compositions.}
#'   \item{ISS_FishLenComps}{Implied sample sizes for fishery length compositions.}
#'   \item{UseFishAgeComps}{Binary indicator array for fishery age comps.}
#'   \item{UseFishLenComps}{Binary indicator array for fishery length comps.}
#'   \item{ObsSrvIdx}{Observed survey index array [region × year × seas x fleet].}
#'   \item{ObsSrvIdx_SE}{Standard error for survey index array.}
#'   \item{UseSrvIdx}{Binary indicator array for survey indices.}
#'   \item{ObsSrvAgeComps}{Observed survey age composition array.}
#'   \item{ObsSrvLenComps}{Observed survey length composition array.}
#'   \item{ISS_SrvAgeComps}{Implied sample sizes for survey age compositions.}
#'   \item{ISS_SrvLenComps}{Implied sample sizes for survey length compositions.}
#'   \item{UseSrvAgeComps}{Binary indicator array for survey age comps.}
#'   \item{UseSrvLenComps}{Binary indicator array for survey length comps.}
#' }
#'
#' @export simulation_data_to_SPoRC
simulation_data_to_SPoRC <- function(sim_env,
                                     y,
                                     sim) {

    # Biologicals
    WAA <- array(sim_env$WAA[,,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
    WAA_fish <- array(sim_env$WAA_fish[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_fish_fleets))
    WAA_srv <- array(sim_env$WAA_srv[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes, sim_env$n_srv_fleets))
    MatAA <- array(sim_env$MatAA[,,1:y,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_ages, sim_env$n_sexes))
    SizeAgeTrans <- if(!is.null(sim_env$SizeAgeTrans)) {
      array(sim_env$SizeAgeTrans[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_ages, sim_env$n_sexes))
      } else NULL
    AgeingError <- sim_env$AgeingError[1:y,,,sim]

    # Tagging
    if(sim_env$use_conv_fish_tagging == 1) {
      keep_tag_cohorts <- which(sim_env$tag_release_indicator[,2] %in% 1:y)
      tag_release_indicator <- sim_env$tag_release_indicator[keep_tag_cohorts,,drop = FALSE]
      obs_conv_tag_fish_recap <- array(sim_env$obs_conv_tag_fish_recap[,,keep_tag_cohorts,,,,,,sim], dim = dim(sim_env$obs_conv_tag_fish_recap)[-length(dim(sim_env$obs_conv_tag_fish_recap))])
      conv_tagged_fish <- array(sim_env$conv_tagged_fish[keep_tag_cohorts,,,,sim], dim = c(dim(sim_env$conv_tagged_fish)[-length(dim(sim_env$conv_tagged_fish))]))
      n_tag_cohorts <- nrow(tag_release_indicator)
    } else {
      tag_release_indicator = obs_conv_tag_fish_recap = conv_tagged_fish = n_tag_cohorts = NULL
    }

    # Fishery Catches
    ObsCatch <- array(sim_env$ObsCatch[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
    ln_sigmaC <- array(sim_env$ln_sigmaC[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
    UseCatch <- array(0, dim = dim(ObsCatch))
    UseCatch[!is.na(ObsCatch) & ObsCatch > 0] <- 1

    # Fishery Indices
    ObsFishIdx <- array(sim_env$ObsFishIdx[, 1:y,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
    ObsFishIdx_SE <- array(sim_env$ObsFishIdx_SE[, 1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
    UseFishIdx <- array(0, dim = dim(ObsFishIdx))
    UseFishIdx[!is.na(ObsFishIdx) & ObsFishIdx > 0] <- 1

    # Fishery Compositions
    ObsFishAgeComps <- array(sim_env$ObsFishAgeComps[, 1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
    ObsFishLenComps <- if(!is.null(sim_env$n_lens)) {
      array(sim_env$ObsFishLenComps[, 1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
    } else NULL
    ISS_FishAgeComps <- array(sim_env$ISS_FishAgeComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
    ISS_FishLenComps <- array(sim_env$ISS_FishLenComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
    UseFishAgeComps <- apply(ObsFishAgeComps, c(1,2,3,6), sum) # sum across, regions, years, seasons, and fleets
    UseFishAgeComps[!is.na(UseFishAgeComps) & UseFishAgeComps > 0] <- 1
    if(!is.null(sim_env$n_lens)) {
      UseFishLenComps <- apply(ObsFishLenComps, c(1,2,3,6), sum)
      UseFishLenComps[!is.na(UseFishLenComps) & UseFishLenComps > 0] <- 1
    } else UseFishLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

    # Survey Indices
    ObsSrvIdx <- array(sim_env$ObsSrvIdx[, 1:y,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
    ObsSrvIdx_SE <- array(sim_env$ObsSrvIdx_SE[, 1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
    UseSrvIdx <- array(0, dim = dim(ObsSrvIdx))
    UseSrvIdx[!is.na(ObsSrvIdx) & ObsSrvIdx > 0] <- 1

    # Survey Compositions
    ObsSrvAgeComps <- array(sim_env$ObsSrvAgeComps[, 1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
    ObsSrvLenComps <- if(!is.null(sim_env$n_lens)) {
      array(sim_env$ObsSrvLenComps[, 1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
    } else NULL
    ISS_SrvAgeComps <- array(sim_env$ISS_SrvAgeComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
    ISS_SrvLenComps <- array(sim_env$ISS_SrvLenComps[, 1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
    UseSrvAgeComps <- apply(ObsSrvAgeComps, c(1,2,3,6), sum) # sum across, regions, years, seasons, and fleets
    UseSrvAgeComps[!is.na(UseSrvAgeComps) & UseSrvAgeComps > 0] <- 1
    if(!is.null(sim_env$n_lens)) {
      UseSrvLenComps <- apply(ObsSrvLenComps, c(1,2,3,6), sum)
      UseSrvLenComps[!is.na(UseSrvLenComps) & UseSrvLenComps > 0] <- 1
    } else UseSrvLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

    return(list(WAA = WAA,
                WAA_fish = WAA_fish,
                WAA_srv = WAA_srv,
                MatAA = MatAA,
                SizeAgeTrans = SizeAgeTrans,
                AgeingError = AgeingError,
                use_conv_fish_tagging = sim_env$use_conv_fish_tagging,
                tag_release_indicator = tag_release_indicator,
                obs_conv_tag_fish_recap = obs_conv_tag_fish_recap,
                conv_tagged_fish = conv_tagged_fish,
                n_tag_cohorts = n_tag_cohorts,
                ObsCatch = ObsCatch,
                ln_sigmaC = ln_sigmaC,
                UseCatch = UseCatch,
                ObsFishIdx = ObsFishIdx,
                ObsFishIdx_SE = ObsFishIdx_SE,
                UseFishIdx = UseFishIdx,
                ObsFishAgeComps = ObsFishAgeComps,
                ObsFishLenComps = ObsFishLenComps,
                ISS_FishAgeComps = ISS_FishAgeComps,
                ISS_FishLenComps = ISS_FishLenComps,
                UseFishAgeComps = UseFishAgeComps,
                UseFishLenComps = UseFishLenComps,
                ObsSrvIdx = ObsSrvIdx,
                ObsSrvIdx_SE = ObsSrvIdx_SE,
                UseSrvIdx = UseSrvIdx,
                ObsSrvAgeComps = ObsSrvAgeComps,
                ObsSrvLenComps = ObsSrvLenComps,
                ISS_SrvAgeComps = ISS_SrvAgeComps,
                ISS_SrvLenComps = ISS_SrvLenComps,
                UseSrvAgeComps = UseSrvAgeComps,
                UseSrvLenComps = UseSrvLenComps
    ))

}

#' Predict ISS fishery compositions under fishing mortality
#'
#' Uses historical ISS fishery compositions and fishing mortality rates
#' to estimate ISS compositions in the projection year. Compositions are
#' scaled relative to the historical maximum fishing mortality with
#' linear interpolation between the minimum and maximum observed ISS values.
#' If historical values are not available, defaults to the mean or zero.
#'
#' @param ISS_FishComps Array of ISS fishery compositions with dimensions
#'   `[region, year, seas, sex, fleet, sim]`.
#' @param Fmort Array of fishing mortality rates with dimensions
#'   `[region, year, seas, fleet, sim]`.
#' @param y Integer, projection year index for prediction.
#' @param sim Integer, simulation index.
#' @param seas Integer, season
#'
#' @returns Array with predicted ISS values for year `y` and season `seas`.
#' @keywords internal
#'
predict_sim_fish_iss_fmort <- function(ISS_FishComps,
                                       Fmort,
                                       y,
                                       seas,
                                       sim
                                       ) {

  # dimensions
  dims <- dim(ISS_FishComps)
  n_regions <- dims[1]
  n_seas <- dims[3]
  n_sexes <- dims[4]
  n_fish_fleets <- dims[4]

  # extract temp vars
  tmp_iss <- ISS_FishComps[, 1:(y-1), seas , , , sim, drop = FALSE]
  tmp_fmort <- Fmort[, 1:(y-1), seas ,, sim, drop = FALSE]

  # container
  iss_container <- array(0, dim = c(n_regions, length(1:y), 1, n_sexes, n_fish_fleets))
  iss_container[, 1:(y-1), seas, , ] <- ISS_FishComps[, 1:(y-1), seas , , , sim] # fill in values back

  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      for(f in 1:n_fish_fleets) {
        # get ISS and Fmort
        iss_vec <- tmp_iss[r, , 1, s, f, ]
        fmort_vec <- tmp_fmort[r, , 1, f, ]
        # remove zeros/NAs
        valid_idx <- which(iss_vec > 0 & !is.na(iss_vec) & !is.na(fmort_vec))
        if(length(valid_idx) > 0) {
          iss_valid <- iss_vec[valid_idx]
          fmort_valid <- fmort_vec[valid_idx]
          # min and max ISS from conditioning period
          min_iss <- min(iss_valid)
          max_iss <- max(iss_valid)
          # max Fmort from conditioning period
          max_fmort_hist <- max(fmort_valid)
          # new Fmort
          fmort_new <- Fmort[r, y, seas, f, sim]
          # scale ISS proportionally to Fmort relative to historical max
          if(max_fmort_hist > 0 && fmort_new >= 0) {
            # linear scaling between min and max ISS
            scaling_factor <- min(fmort_new / max_fmort_hist, 1)  # cap scaling at 1
            new_iss <- min_iss + scaling_factor * (max_iss - min_iss) # linear scaling
          } else {
            # use mean ISS if conditions not met ...
            new_iss <- mean(iss_valid)
          }
          iss_container[r, y, 1, s, f] <- new_iss
        } else {
          iss_container[r, y, 1, s, f] <- 0
        }
      } # end f loop
    } # end s loop
  } # end r loop

  return(iss_container)
}
