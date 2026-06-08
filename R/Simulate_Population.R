#' Construct and populate a simulation execution environment
#'
#' Creates a new R environment populated with all objects from \code{sim_list}
#' and binds the SPoRC simulation functions required by
#' \code{\link{run_annual_cycle}}. Isolating the simulation state in a
#' dedicated environment prevents name collisions with the calling frame and
#' allows \code{with()} / \code{<<-} assignment patterns used internally by
#' the annual-cycle helpers to modify shared state without polluting the
#' global workspace.
#'
#' @param sim_list Named list returned by \code{\link{Setup_Sim_Rec}} (or the
#'   last upstream setup function called). All elements are copied into the
#'   new environment via \code{list2env}.
#'
#' @return A new environment (parent = calling frame) containing every element
#'   of \code{sim_list} as a named object, plus bound references to the
#'   following SPoRC simulation functions: \code{generate_initial_age_structure},
#'   \code{generate_recruitment}, \code{apply_pop_dy},
#'   \code{generate_fishery_catch_comp_idx}, \code{generate_survey_comp_idx},
#'   \code{release_conv_tags}, \code{generate_fishery_conv_tags_recap},
#'   \code{Get_Det_Recruitment}, \code{Get_Init_NAA},
#'   \code{predict_sim_fish_iss_fmort}, \code{rho_trans},
#'   \code{simulate_comps}, \code{simulate_conv_tag_fish_recaptures}.
#'
#'
#' @export Setup_sim_env
#' @family Simulation Setup
#'
#' @examples
#' \dontrun{
#' sim_env <- Setup_sim_env(sim_list)
#' }
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
  sim_env$predict_sim_fish_iss_fmort <- predict_sim_fish_iss_fmort
  sim_env$rho_trans <- rho_trans
  sim_env$simulate_comps <- simulate_comps
  sim_env$simulate_conv_tag_fish_recaptures <- simulate_conv_tag_fish_recaptures

  # output into simulation environment
  list2env(sim_list, envir = sim_env)

  return(sim_env)
}

#' Simulate age or length compositions
#'
#' Draws observed composition samples (by age or length) for a single
#' region–year–fleet–season–simulation cell, supporting multinomial,
#' Dirichlet-multinomial, and logistic-normal likelihoods. Ageing error is
#' optionally applied post-draw. Three composition aggregation structures are
#' handled: sex-split (\code{comp_type = 1}), joint across sexes
#' (\code{comp_type = 2}), and spatially aggregated across all regions
#' (\code{comp_type = 0}). The sentinel value \code{comp_type = 999} or
#' \code{comp_like = 999} causes the function to return \code{Obs} unchanged.
#'
#' When \code{pop_specific = TRUE}, compositions are simulated separately for
#' each population, extending all relevant inputs (e.g., sample size,
#' dispersion, and correlation parameters) to include a population dimension.
#' In this case, aggregation across regions (\code{comp_type = 0}) is performed
#' within population, and results are written to \code{Obs[p, ...]}.
#'
#' For joint compositions (\code{comp_type = 2}), the Kronecker product
#' \code{diag(n_sexes) ⊗ AgeingError} is used to apply ageing error across
#' the combined age–sex vector. For aggregated compositions
#' (\code{comp_type = 0}), the draw is only executed when \code{r == n_regions}
#' (i.e., on the final region pass), and uses region- and sex-marginalised
#' expected proportions.
#'
#' @param r Integer. Region index.
#' @param y Integer. Year index.
#' @param f Integer. Fleet index (fishery or survey).
#' @param seas Integer. Season index.
#' @param sim Integer. Simulation replicate index.
#' @param Exp Array. Expected compositions
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_cat × n_sexes × n_fleets × n_sims]}.
#' @param ISS Array. Integer sample sizes
#'   \code{[n_regions × n_yrs × n_seas × n_sexes × n_fleets × n_sims]}.
#'   Used when \code{pop_specific = FALSE}.
#' @param AgeingError Array. Ageing error transition matrices
#'   \code{[n_yrs × n_obs_ages × n_ages × n_sims]}. Ignored when
#'   \code{age_or_len = 1}.
#' @param comp_like Integer vector \code{[n_fleets]}. Likelihood type per
#'   fleet: \code{0} = multinomial, \code{1} = Dirichlet-multinomial,
#'   \code{2}–\code{4} = logistic-normal variants.
#' @param ln_theta Array. Log overdispersion or log-variance parameters
#'   \code{[n_regions × n_sexes × n_fleets]}. Used when
#'   \code{pop_specific = FALSE}.
#' @param corr_pars Array. Correlation parameters for logistic-normal
#'   likelihoods \code{[n_regions × n_sexes × n_fleets × n_corr_pars]}.
#' @param ln_theta_agg Numeric vector \code{[n_fleets]}. Log overdispersion
#'   for spatially aggregated compositions (\code{comp_type = 0}).
#' @param corr_pars_agg Numeric vector \code{[n_fleets]}. Correlation
#'   parameter(s) for aggregated logistic-normal compositions.
#' @param comp_type Integer matrix \code{[n_yrs × n_fleets]}. Aggregation
#'   structure: \code{0} = aggregated across regions, \code{1} = split by
#'   sex, \code{2} = joint across sexes, \code{999} = no data (skip).
#' @param n_sexes Integer. Number of sexes.
#' @param n_pop Integer. Number of populations.
#' @param n_regions Integer. Number of regions.
#' @param n_cat Integer. Number of composition categories (ages or lengths).
#' @param Obs Array. Observed compositions container with the same dimensions
#'   as \code{Exp}. Simulated values are written in-place.
#' @param pop_specific Logical. If \code{TRUE}, simulate compositions
#'   separately for each population using population-specific inputs.
#' @param age_or_len Integer. Indicator for composition type:
#'   \code{0} = age compositions (apply ageing error),
#'   \code{1} = length compositions (no ageing error).
#'
#' @param ISS_pop Array. Population-specific sample sizes
#'   \code{[n_pop × n_regions × n_yrs × n_seas × n_sexes × n_fleets × n_sims]}.
#'   Used when \code{pop_specific = TRUE}.
#' @param pop_comp_like Integer vector \code{[n_fleets]}. Likelihood type per
#'   fleet for population-specific compositions.
#' @param pop_comp_type Integer matrix \code{[n_yrs × n_fleets]}. Aggregation
#'   structure for population-specific compositions.
#' @param ln_pop_theta Array. Log overdispersion parameters
#'   \code{[n_pop × n_regions × n_sexes × n_fleets]}.
#' @param pop_corr_pars Array. Correlation parameters for logistic-normal
#'   likelihoods
#'   \code{[n_pop × n_regions × n_sexes × n_fleets × n_corr_pars]}.
#' @param ln_pop_theta_agg Numeric array \code{[n_pop × n_fleets]}. Log
#'   overdispersion for population-specific aggregated compositions.
#' @param pop_corr_pars_agg Numeric array \code{[n_pop × n_fleets]}.
#'   Correlation parameter(s) for population-specific aggregated
#'   logistic-normal compositions.
#'
#' @return The \code{Obs} array with simulated composition draws filled in at
#'   the appropriate slice. When \code{pop_specific = FALSE}, values are written
#'   to \code{[r, y, seas, , , f, sim]}; when \code{pop_specific = TRUE}, values
#'   are written to \code{[p, r, y, seas, , , f, sim]}. All other slices are
#'   unchanged.
#'
#' @keywords internal
simulate_comps <- function(r,
                           y,
                           f,
                           seas,
                           sim,
                           Exp,
                           ISS  = NULL,
                           AgeingError,
                           comp_like  = NULL,
                           ln_theta  = NULL,
                           corr_pars  = NULL,
                           ln_theta_agg  = NULL,
                           corr_pars_agg  = NULL,
                           comp_type = NULL,
                           n_sexes,
                           n_pop = NULL,
                           n_regions,
                           n_cat,
                           Obs,
                           pop_specific = FALSE,
                           ISS_pop = NULL,
                           pop_comp_like = NULL,
                           pop_comp_type = NULL,
                           ln_pop_theta = NULL,
                           pop_corr_pars = NULL,
                           ln_pop_theta_agg = NULL,
                           pop_corr_pars_agg = NULL,
                           age_or_len = 0) {

  if(!pop_specific && (comp_type[y,f] == 999 || comp_like[f] == 999)) return(Obs)
  if(pop_specific && (pop_comp_type[y,f] == 999 || pop_comp_like[f] == 999)) return(Obs)

  # helper functions
  get_expected <- function(prob_vec) prob_vec / sum(prob_vec)
  apply_error <- function(mat, age_or_len, AgeingError) {
    if(age_or_len == 0) return(mat %*% AgeingError)
    if(age_or_len == 1) return(mat)
  }

  if(!pop_specific) {
    if(age_or_len == 0) {
      if(comp_type[y,f] %in% c(0,1)) age_error_mat <- AgeingError[y,,,sim]
      if(comp_type[y,f] == 2) age_error_mat <- kronecker(diag(n_sexes), AgeingError[y,,,sim])
    }
  } else {
    if(age_or_len == 0) {
      if(pop_comp_type[y,f] %in% c(0,1)) age_error_mat <- AgeingError[y,,,sim]
      if(pop_comp_type[y,f] == 2) age_error_mat <- kronecker(diag(n_sexes), AgeingError[y,,,sim])
    }
  }

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
            ), age_or_len, age_error_mat),
            dim = dim(Obs[1,y,seas,,1,f,sim, drop = FALSE])
          )
        }
      }
    } # end if not pop-specific

  if(pop_specific == TRUE) {
    for(p in 1:n_pop) {
      # Split by sex
      if(pop_comp_type[y,f] == 1) {
        for(s in 1:n_sexes) {

          tmp_prob <- Exp[p,r,y,seas,,s,f,sim, drop = FALSE]# extract compositions

          # multinomial
          if(pop_comp_like[f] == 0) {
            Obs[p,r,y,seas,,s,f,sim] <- array(
              apply_error(as.vector(
                stats::rmultinom(n = 1, ISS_pop[p,r,y,seas,s,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
              dim = dim(Obs[p,r,y,seas,,s,f,sim, drop = FALSE])
            )

            # dirichlet-multinomial
          } else if(pop_comp_like[f] == 1) {
            Obs[p,r,y,seas,,s,f,sim] <- array(
              apply_error(as.vector(
                rdirM(
                  n = 1,
                  N = ISS_pop[p,r,y,seas,s,f,sim],
                  alpha = (exp(ln_pop_theta[p,r,s,f]) * ISS_pop[p,r,y,seas,s,f,sim]) * get_expected(tmp_prob)
                )
              ), age_or_len, age_error_mat),
              dim = dim(Obs[p,r,y,seas,,s,f,sim, drop = FALSE])
            )

            # logistic normal
          } else if(pop_comp_like[f] %in% 2:4) {
            Obs[p,r,y,seas,,s,f,sim] <- array(
              apply_error(as.vector(
                rlogistnormal(
                  exp = get_expected(tmp_prob),
                  pars = c(exp(ln_pop_theta[p,r,s,f]), rho_trans(pop_corr_pars[p,r,s,f,])),
                  comp_like = pop_comp_like[f],
                  n_sexes = n_sexes
                )
              ), age_or_len, age_error_mat),
              dim = dim(Obs[p,r,y,seas,,s,f,sim, drop = FALSE])
            )
          }

        } # end s loop
      } # end split by sex

      # Joint compositions
      if(pop_comp_type[y,f] == 2) {

        tmp_prob <- Exp[p,r,y,seas,,,f,sim, drop = FALSE] # extract compositions

        # multinomial
        if(pop_comp_like[f] == 0) {
          Obs[p,r,y,seas,,,f,sim] <- array(
            apply_error(as.vector(stats::rmultinom(1, ISS_pop[p,r,y,seas,1,f,sim], get_expected(tmp_prob))),
                        age_or_len, age_error_mat),
            dim = dim(Obs[p,r,y,seas,,,f,sim, drop = FALSE])
          )

          # dirichlet-multinomial
        } else if(pop_comp_like[f] == 1) {
          Obs[p,r,y,seas,,,f,sim] <- array(
            apply_error(as.vector(
              rdirM(
                n = 1,
                N = ISS_pop[p,r,y,seas,1,f,sim],
                alpha = (exp(ln_pop_theta[p,r,1,f]) * ISS_pop[p,r,y,seas,1,f,sim]) * get_expected(tmp_prob)
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,r,y,seas,,,f,sim, drop = FALSE])
          )

          # logistic normal
        } else if(pop_comp_like[f] %in% 2:4) {
          Obs[p,r,y,seas,,,f,sim] <- array(
            apply_error(as.vector(
              rlogistnormal(
                exp = get_expected(tmp_prob),
                pars = c(exp(ln_pop_theta[p,r,1,f]), rho_trans(pop_corr_pars[p,r,1,f,])),
                comp_like = pop_comp_like[f],
                n_sexes = n_sexes
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,r,y,seas,,,f,sim, drop = FALSE])
          )
        }

      } # end joint compositions

      # Aggregated comps across regions
      if(r == n_regions && pop_comp_type[y,f] == 0) {

        # extract compositions
        tmp_prob <- apply(Exp[p,,y,seas,,,f,sim, drop = FALSE], 5, sum)
        tmp_prob <- tmp_prob / sum(tmp_prob)

        # multinomial
        if(pop_comp_like[f] == 0) {
          Obs[p,1,y,seas,,1,f,sim] <- array(
            apply_error(as.vector(stats::rmultinom(1, ISS_pop[p,1,y,seas,1,f,sim], get_expected(tmp_prob))), age_or_len, age_error_mat),
            dim = dim(Obs[p,1,y,seas,,1,f,sim, drop = FALSE])
          )

          # dirichlet-multinomial
        } else if(pop_comp_like[f] == 1) {
          Obs[p,1,y,seas,,1,f,sim] <- array(
            apply_error(as.vector(
              rdirM(
                n = 1,
                N = ISS_pop[p,1,y,seas,1,f,sim],
                alpha = (exp(ln_pop_theta_agg[p,f]) * ISS_pop[p,1,y,seas,1,f,sim]) * get_expected(tmp_prob)
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,1,y,seas,,1,f,sim, drop = FALSE])
          )

          # logistic normal
        } else if(pop_comp_like[f] %in% 2:4) {
          Obs[p,1,y,seas,,1,f,sim] <- array(
            apply_error(as.vector(
              rlogistnormal(
                exp = get_expected(tmp_prob),
                pars = c(exp(ln_pop_theta_agg[p,f]), rho_trans(pop_corr_pars_agg[p,f])),
                comp_like = pop_comp_like[f],
                n_sexes = n_sexes
              )
            ), age_or_len, age_error_mat),
            dim = dim(Obs[p,1,y,seas,,1,f,sim, drop = FALSE])
          )
        }
      }
    } # end p loop
  } # end if pop-specific

  return(Obs)
}

#' Simulate conventional tag recaptures for fishery fleets
#'
#' Draws observed tag recapture counts for a single liberty–season–cohort cell
#' from predicted recapture arrays, supporting six likelihood structures:
#' Poisson, negative binomial, and release- or recovery-conditioned
#' multinomial and Dirichlet-multinomial. Dimensions absent from
#' \code{tag_recaptures_attr} are marginalised by summing over them, and all
#' recaptures are placed into index 1 of the corresponding dimension in the
#' output array.
#'
#' For release-conditioned likelihoods (\code{2}, \code{4}), predicted
#' recaptures are expressed as proportions of total tags released. A
#' "not-recaptured" bin is appended to complete the probability vector before
#' drawing and removed before assignment. For recovery-conditioned likelihoods
#' (\code{3}, \code{5}), the draw is conditioned on the total predicted
#' recapture count with no not-recaptured bin needed. The overdispersion
#' parameter \code{ln_conv_fish_tag_theta} governs the negative-binomial size
#' parameter and the Dirichlet-multinomial concentration scaling, and is
#' ignored for Poisson and multinomial likelihoods.
#'
#' @param conv_fish_tag_like Integer. Likelihood for tag recaptures:
#'   \code{0} = Poisson, \code{1} = negative binomial,
#'   \code{2} = multinomial (release-conditioned),
#'   \code{3} = multinomial (recovery-conditioned),
#'   \code{4} = Dirichlet-multinomial (release-conditioned),
#'   \code{5} = Dirichlet-multinomial (recovery-conditioned).
#' @param tag_recaptures_attr Character string specifying which biological
#'   dimensions are attended in the recapture likelihood. Built from any
#'   combination of \code{"p"} (population), \code{"a"} (age), and \code{"s"}
#'   (sex), joined by underscores. Region and fleet are always retained.
#'   Unattended dimensions are marginalised and output into index 1.
#' @param conv_tagged_fish Array of released tagged fish
#'   \code{[n_conv_tag_cohorts × n_pop × n_ages × n_sexes × n_sims]}. Used
#'   as the release sample size for release-conditioned likelihoods.
#' @param pred_conv_tag_fish_recap Array of predicted recaptures
#'   \code{[conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop ×
#'   n_regions × n_ages × n_sexes × n_fish_fleets × n_sims]}.
#' @param obs_conv_tag_fish_recap Array of observed recaptures with the same
#'   dimensions as \code{pred_conv_tag_fish_recap}. Simulated values are
#'   written in-place at the \code{[ry, rseas, tc, ...]} slice.
#' @param ln_conv_fish_tag_theta Numeric. Log overdispersion: negative
#'   binomial size = \code{exp(ln_conv_fish_tag_theta)}; Dirichlet-multinomial
#'   concentration = \code{exp(ln_conv_fish_tag_theta) × N × p}.
#' @param ry Integer. Years-at-liberty index (first dimension of recapture
#'   arrays).
#' @param rseas Integer. Recovery season index.
#' @param tc Integer. Tag cohort index.
#' @param sim Integer. Simulation replicate index.
#' @param n_pop Integer. Number of populations.
#' @param n_regions Integer. Number of regions.
#' @param n_ages Integer. Number of age classes.
#' @param n_sexes Integer. Number of sexes.
#' @param n_fish_fleets Integer. Number of fishery fleets.
#'
#' @return The \code{obs_conv_tag_fish_recap} array with simulated recaptures
#'   filled in at \code{[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx,
#'   flt_idx, sim]}. Marginalised dimensions are fixed at index 1.
#'
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

#' Marginalise conventional fishery tag arrays across unattended dimensions
#'
#' Collapses population, age, and/or sex dimensions of a tag count array
#' \code{[n_pop × n_ages × n_sexes]} by summing over dimensions absent from
#' \code{tag_recaptures_attr}, placing the result into index 1 of the
#' corresponding dimension and zeroing all other indices. Region and fleet are
#' not handled here (they are managed at the calling level). The function is
#' used to align the release array \code{conv_tagged_fish} with the attended
#' resolution of the recapture likelihood.
#'
#' @param vals Numeric vector or array of tag counts, interpreted as a
#'   \code{[n_pop × n_ages × n_sexes]} array.
#' @param tag_recaptures_attr Character string specifying attended dimensions.
#'   Same format as \code{conv_fish_tag_attr} in
#'   \code{\link{Setup_Sim_Tagging}}: any combination of \code{"p"},
#'   \code{"a"}, \code{"s"} joined by underscores.
#' @param n_pop Integer. Number of populations.
#' @param n_ages Integer. Number of age classes.
#' @param n_sexes Integer. Number of sexes.
#'
#' @return Array \code{[n_pop × n_ages × n_sexes]} with unattended dimensions
#'   summed into index 1 and all other indices set to zero.
#'
#'
#' @keywords internal
marginalize_conv_fish_tags <- function(vals,
                                       tag_recaptures_attr,
                                       n_pop,
                                       n_ages,
                                       n_sexes) {

  full_dims  <- c(n_pop, n_ages, n_sexes)
  attr_parts <- strsplit(tag_recaptures_attr, "_")[[1]]
  tmp <- array(vals, dim = full_dims)  # temporary array

  if (!("p" %in% attr_parts) && n_pop > 1) {
    summed        <- apply(tmp, c(2, 3), sum) # collapse pop
    tmp[]         <- 0 # zero out stuff
    tmp[1, , ]  <- summed # input into 1st pop
  }
  if (!("a" %in% attr_parts) && n_ages > 1) {
    summed        <- apply(tmp, c(1, 3), sum) # collapse ages
    tmp[]         <- 0 # zero out
    tmp[, 1, ]  <- summed # input into 1st age
  }
  if (!("s" %in% attr_parts) && n_sexes > 1) {
    summed        <- apply(tmp, c(1, 2), sum)  # collapse sexes
    tmp[]         <- 0 # zero out
    tmp[, , 1]  <- summed # input into 1st sex
  }

  return(tmp)
}

#' Initialise age structure for a simulation replicate
#'
#' Simulates or reads in initial age deviations and calls
#' \code{\link{Get_Init_NAA}} to compute both the fished and unfished
#' equilibrium numbers-at-age for year 1 and season 1. Results are written
#' directly into the simulation environment arrays \code{NAA} and
#' \code{NAA0}. This function is called once per simulation replicate at
#' \code{y = 1} by \code{\link{run_annual_cycle}}.
#'
#' Initial deviation sharing follows the same logic as the estimation model:
#' deviations are drawn once per population when \code{n_pop > 1}, or once
#' per region when \code{n_pop = 1} and \code{init_dd = 0} (local
#' density-dependence). If \code{ln_InitDevs_input} exists in the simulation
#' environment, those values are used directly rather than simulating new
#' draws. Populations with \code{R0 = 0} receive zero deviations. The
#' equilibrium solver uses \code{init_iter = n_ages × 5} iterations.
#'
#' @param y Integer. Year index (must be \code{1}).
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by
#'   \code{\link{Setup_sim_env}}. Modified in place: \code{$ln_InitDevs},
#'   \code{$NAA[,,1,1,,,sim]}, and \code{$NAA0[,,1,1,,,sim]} are updated.
#'
#' @return \code{invisible(NULL)}. All modifications are made by reference
#'   within \code{sim_env}.
#'
#'
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

          # get init devs devs
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])
          # simulate initial age deviations
          if(is.null(tmp_ln_init_devs)) tmp_ln_init_devs <- stats::rnorm(n_ages-1, -exp(ln_sigmaR[1,p,sigma_idx])^2/2, exp(ln_sigmaR[1,p,sigma_idx]))
        }

        # input age deviations
        if(R0[p,r,1,sim] != 0) {
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
      n_fish_fleets = n_fish_fleets, # number of fishery fleets
      seasdur = seasdur,  # fracion of time in season
      natmort = array(natmort[,,1,,,sim], dim = c(n_pop, n_regions, n_ages, n_sexes)), # natural mortality in first year
      init_F = init_F, # initial F applied (0 for unfished)\
      dmr = array(dmr[,1,,,sim], dim = c(n_regions, n_seas, n_fish_fleets)), # discard mortality rate
      fish_sel = array(fish_sel[,,1,,,,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # total fishery selectivity in first year
      ret_sel = array(ret_sel[,,1,,,,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # retained selectivity in first year
      R0_r = if(use_rinit == 0) array(R0[,,1,sim], dim = c(n_pop, n_regions)) else array(rinit[,,sim], dim = c(n_pop, n_regions)), # regional mean or virgin recruitment
      rec_seas_prop = array(rec_seas_prop[,,sim], dim = c(n_pop, n_seas)), # recruitment seasonal apportionment
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
      init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)), # initial F applied (0 for unfished)
      dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)), # dmr applied (0 for unfished)
      n_seas = n_seas, # seasons
      n_fish_fleets = n_fish_fleets, # number of fishery fleets
      seasdur = seasdur,  # fracion of time in season
      rec_seas_prop = array(rec_seas_prop[,,sim], dim = c(n_pop, n_seas)), # recruitment seasonal apportionment
      fish_sel = array(fish_sel[,,1,,,,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # total fishery selectivity in first year
      ret_sel = array(ret_sel[,,1,,,,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # retained selectivity in first year
      R0_r = if(use_rinit == 0) array(R0[,,1,sim], dim = c(n_pop, n_regions)) else array(rinit[,,sim], dim = c(n_pop, n_regions)), # regional mean or virgin recruitment
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

#' Generate recruitment for a simulation year
#'
#' Computes total recruitment for year \code{y} by first obtaining
#' deterministic expected recruitment from \code{\link{Get_Det_Recruitment}}
#' and then multiplying by lognormal deviations (bias-corrected). Recruitment
#' is apportioned across sexes and seasons and written into
#' \code{sim_env$NAA[p, r, y, 1, 1, s, sim]} (first season, age-1 slot).
#' Unfished NAA (\code{NAA0}) is synchronised to match fished NAA at
#' recruitment. If \code{Rec_input} exists in the environment and covers year
#' \code{y}, those values override the stochastic draw entirely.
#'
#' Recruitment deviation sharing follows the same population/region logic as
#' \code{\link{generate_initial_age_structure}}: deviations are drawn once
#' per population (\code{n_pop > 1}) or once per region (\code{n_pop = 1},
#' local density-dependence). Populations with \code{R0 = 0} receive zero
#' deviations. \code{sigma_idx} selects the natal region's \code{ln_sigmaR}
#' for the bias-correction term.
#'
#' @param y Integer. Year index for which recruitment is generated.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by
#'   \code{\link{Setup_sim_env}}. Modified in place:
#'   \code{$ln_RecDevs[p, r, y, sim]}, \code{$Rec[p, r, y, sim]},
#'   \code{$NAA[p, r, y, 1, 1, s, sim]}, and
#'   \code{$NAA0[p, r, y, 1, 1, s, sim]} are updated.
#'
#' @return \code{invisible(NULL)}. All modifications are made by reference
#'   within \code{sim_env}.
#'
#'
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
                                       rec_region_prop = array(t(apply(R0[,,y,sim, drop = FALSE], c(1), function(x) x / sum(x))), dim = c(n_pop, n_regions)), # get R0 proportion
                                       rec_seas_prop = array(rec_seas_prop[,,sim], dim = c(n_pop, n_seas)),
                                       h = array(h[,,y,sim], dim = c(n_pop, n_regions)),
                                       n_pop = n_pop,
                                       n_regions = n_regions,
                                       n_ages = n_ages,
                                       natal_region = natal_region,

                                       # Note: Using first year and female quantities to compute unfished SSB0
                                       sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,1,1,sim], dim = c(n_pop, n_regions)),
                                       WAA = array(WAA[,,1,,,1,sim], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                       MatAA = array(MatAA[,,1,,,1,sim], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                       natmort = array(natmort[,,1,,1,sim], dim = c(n_pop, n_regions, n_ages)),
                                       stray_rate = array(stray_rate[,1,sim], dim = c(n_pop)),
                                       Movement = array(Movement[,,,1,,,1,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                       sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,1,,1,sim], dim = c(n_pop, n_regions, n_regions, n_ages)),
                                       SSB_vals = array(SSB[,,,sim], dim = c(n_pop, n_regions, n_yrs)),
                                       n_fish_fleets = n_fish_fleets,
                                       t_spawn = t_spawn,
                                       n_seas = n_seas,
                                       spawn_seas = spawn_seas,
                                       seasdur = seasdur,
                                       do_recruits_move = do_recruits_move,
                                       init_F = init_F, # initial F applied
                                       dmr = array(dmr[,1,,,sim], dim = c(n_regions, n_seas, n_fish_fleets)), # discard mortality rate
                                       fish_sel = array(fish_sel[,,1,,,1,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # total fishery selectivity in first year
                                       ret_sel = array(ret_sel[,,1,,,1,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)) # retained fishery selectivity in first year

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
          tmp_total_rec <- Rec_input[p,r,y,sim]
        } else {

          # get rec devs
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])
          tmp_ln_rec_devs <- stats::rnorm(1, 0, exp(ln_sigmaR[2, p, sigma_idx]))

          if(R0[p,r,y,sim] != 0) {
            sim_env$ln_RecDevs[p,r,y,sim] <- tmp_ln_rec_devs
          } else sim_env$ln_RecDevs[p,r,y,sim] <- 0

          # compute rec
          tmp_total_rec <- tmp_det_rec[p,r] * exp(sim_env$ln_RecDevs[p,r,y,sim] - exp(ln_sigmaR[2,p,sigma_idx])^2/2)
        }

        # input recruitment into first season
        for(s in 1:n_sexes) sim_env$NAA[p,r,y,1,1,s,sim] <- tmp_total_rec * rec_seas_prop[p,1,sim] * sexratio[p,r,y,s,sim]

        sim_env$Rec[p,r,y,sim] <- tmp_total_rec # Save annual recruitment estimates
        sim_env$NAA0[p,r,y,1,1,,sim] = NAA[p,r,y,1,1,,sim] # populate unfished NAA

      } # end r loop
    } # end p loop
  })
}

#' Apply population dynamics within a simulation year
#'
#' Executes the full within-year population dynamics loop for year \code{y}:
#' seasonal recruitment apportionment (seasons 2+), movement, Baranov
#' catch-equation mortality, age advancement into the following year, and
#' spawning-season biomass calculations (total biomass, SSB, dynamic \eqn{B_0},
#' and effective SSB for multi-population natal homing). Both fished
#' (\code{NAA}) and unfished (\code{NAA0}) trajectories are tracked in
#' parallel. For single-season multi-population models,
#' \code{sgl_seas_spawning_movement} is applied to \code{NAA} and
#' \code{NAA0} prior to computing spawning biomass quantities. Single-sex
#' models have SSB and \eqn{B_0} multiplied by 0.5 to obtain female-only
#' spawning biomass.
#'
#' Pre- and post-movement snapshots are stored in \code{NAA_bef} and
#' \code{NAA_aft} respectively. Movement is only applied when
#' \code{n_regions > 1}; recruits (\code{a = 1}) are excluded from movement
#' when \code{do_recruits_move = 0}.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by
#'   \code{\link{Setup_sim_env}}. Modified in place: \code{$ZAA},
#'   \code{$NAA}, \code{$NAA0}, \code{$NAA_bef}, \code{$NAA_aft},
#'   \code{$Total_Biom}, \code{$SSB}, \code{$Dynamic_SSB0}, and
#'   \code{$eff_SSB} are updated.
#'
#' @return \code{invisible(NULL)}. All modifications are made by reference
#'   within \code{sim_env}.
#'
#'
#' @keywords internal
apply_pop_dy <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {

    for(seas in 1:n_seas) {

      # apportion recruitment across seasons
      if(seas > 1) {
        for(p in 1:n_pop) {
          for(r in 1:n_regions) {
            for(s in 1:n_sexes) {

              # accumulate recruits - fished
              sim_env$NAA[p,r,y,seas,1,s,sim] <- NAA[p,r,y,seas,1,s,sim] +
                Rec[p,r,y,sim] * rec_seas_prop[p,seas,sim] * sexratio[p,r,y,s,sim]

              # accumulate recruits - unfished
              sim_env$NAA0[p,r,y,seas,1,s,sim] <- NAA0[p,r,y,seas,1,s,sim] +
                Rec[p,r,y,sim] * rec_seas_prop[p,seas,sim] * sexratio[p,r,y,s,sim]

            } # end s loop
          } # end r loop
        } # end p loop
      } # end if seas > 1

      # Mortality and Ageing
      tmp_Fmort <- array(Fmort[,y,seas,,sim], dim = c(n_regions, n_fish_fleets))
      tmp_dmr <- array(dmr[,y,seas,,sim], dim = c(n_regions, n_fish_fleets))
      tmp_natmort <- array((natmort[,,y,,,sim] * seasdur[seas]), dim = c(n_pop, n_regions, n_ages, n_sexes))

      for(p in 1:n_pop) {

        # Get retained catch selectivity
        tmp_fish_sel <- array(fish_sel[p,,y,seas,,,,sim], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)) # total selectivtiy
        tmp_ret_sel <- array(ret_sel[p,,y,seas,,,,sim], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)) # retained selectivity
        tmp_ret_FAA <- apply(sweep(tmp_fish_sel * tmp_ret_sel, c(1,4), tmp_Fmort, "*"), c(1,2,3), sum) # apply Frate to retained selectivity
        tmp_ret_FAA <- array(tmp_ret_FAA, dim = c(n_regions, n_ages, n_sexes)) # reshape

        # Get discarded catch selectivity
        tmp_disc_FAA <- apply(sweep(tmp_fish_sel * (1 - tmp_ret_sel), c(1,4), tmp_Fmort * tmp_dmr, "*"), c(1,2,3), sum) # apply Frate and dmr to discarded selectivity

        # Get natural mortality
        tmp_nm  <- array(tmp_natmort[p,,,,drop=FALSE], dim = c(n_regions, n_ages, n_sexes)) # reshape natural mortality

        # Get total mortality
        sim_env$ZAA[p,,y,seas,,,sim] <- tmp_nm + tmp_ret_FAA + tmp_disc_FAA
      }

      # Movement
      # Record values prior to movement
      NAA_bef[,,y,seas,,,sim] = NAA[,,y,seas,,,sim]

      if(n_regions > 1) {
        for(p in 1:n_pop) {

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
        } # end p loop

        # Record values after movement
        NAA_aft[,,y,seas,,,sim] = NAA[,,y,seas,,,sim]

      } # only compute if spatial

      if(seas < n_seas) { # Within year seasonal mortality
        sim_env$NAA[,,y,seas + 1,1:n_ages,,sim] = NAA[,,y,seas,1:n_ages,,sim] * exp(-ZAA[,,y,seas,1:n_ages,,sim]) # fished
        sim_env$NAA0[,,y,seas + 1,1:n_ages,,sim] <- NAA0[,,y,seas,1:n_ages,,sim] * exp(-(tmp_natmort[,,1:n_ages,])) # unfished
      } else {
        # Advance into the next year, season 1
        sim_env$NAA[,,y+1,1,2:n_ages,,sim] = NAA[,,y,seas,1:(n_ages-1),,sim] * exp(-ZAA[,,y,seas,1:(n_ages-1),,sim]) # fished
        sim_env$NAA[,,y+1,1,n_ages,,sim] = NAA[,,y+1,1,n_ages,,sim] + NAA[,,y,seas,n_ages,,sim] * exp(-ZAA[,,y,seas,n_ages,,sim]) # Acuumulate plus group (fished)
        sim_env$NAA0[,,y+1,1,2:n_ages,,sim] = NAA0[,,y,seas,1:(n_ages-1),,sim] * exp(-(tmp_natmort[,,1:(n_ages - 1),])) # fished
        sim_env$NAA0[,,y+1,1,n_ages,,sim] = NAA0[,,y+1,1,n_ages,,sim] + NAA0[,,y,seas,n_ages,,sim] * exp(-(tmp_natmort[,,n_ages,])) # Acuumulate plus group (unfished)
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
        sim_env$SSB[,, y, sim] <- apply(tmp_NAA_spawn[,, , , , 1, 1,drop = FALSE] *
                                          WAA[,, y, spawn_seas, , 1, sim,drop = FALSE] *
                                          MatAA[,, y, spawn_seas, , 1, sim,drop = FALSE] *
                                          exp(-ZAA[,, y, spawn_seas, , 1, sim,drop = FALSE] * t_spawn), c(1,2), sum)

        # Get dynamic B0
        SSB0_array <- tmp_NAA0_spawn[,, , , , 1, 1,drop = FALSE] *  WAA[,,  y, spawn_seas, , 1, sim, drop = FALSE] * MatAA[,,y, spawn_seas, , 1, sim, drop = FALSE]
        mort_spawn <- exp(-natmort[,, y, , 1, sim, drop = FALSE] * t_spawn * seasdur[spawn_seas])
        mort_spawn <- array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
        sim_env$Dynamic_SSB0[,,y,sim] <- apply(SSB0_array * mort_spawn, c(1,2), sum) # Dynamic B0

        if(n_sexes == 1) { # If single sex model, multiply SSB calculations by 0.5
          sim_env$SSB[,,y,sim] <- SSB[,,y,sim] * 0.5
          sim_env$Dynamic_SSB0[,,y,sim] <- Dynamic_SSB0[,,y,sim] * 0.5
        }

        # Accumulate effective SSB at each population's natal region
        # across all source populations (captures stray contributions)
        # then inside y/sim loops:
        if(n_pop > 1) {
          n_pop_in_region = array(0, dim = n_regions)
          for(p in 1:n_pop) n_pop_in_region[natal_region[p]] = n_pop_in_region[natal_region[p]] + 1
          for(p2 in 1:n_pop) {
            for(p in 1:n_pop) {
              if(p == p2) {
                sim_env$eff_SSB[p2, y, sim] = sim_env$eff_SSB[p2, y, sim] + SSB[p, natal_region[p2], y, sim]
              } else {
                n_receivers = n_pop_in_region[natal_region[p2]]
                sim_env$eff_SSB[p2, y, sim] = sim_env$eff_SSB[p2, y, sim] + (stray_rate[p,y,sim] / n_receivers) * SSB[p, natal_region[p2], y, sim]
              }
            }
          }
        } else sim_env$eff_SSB[1, y, sim] = sum(SSB[1,,y,sim])

      } # if season = spawning season
    } # end seas loop
  })
}

#' Generate fishery catches, compositions, and indices in simulation
#'
#' Applies Baranov's catch equation to compute retained catch-at-age
#' (\code{CAA}) and dead discard catch-at-age (\code{DAA}) for all
#' populations, regions, seasons, and fleets, derives catch-at-length
#' (\code{CAL} and \code{DAL}) when a size-age transition matrix is available,
#' and generates observed catch and discard indices (with lognormal error),
#' fishery abundance or biomass indices, and age and length composition
#' samples for both retained and discarded catch. Composition sampling calls
#' \code{\link{simulate_comps}} and respects the likelihood type
#' (\code{comp_fishage_like}, \code{comp_fishlen_like}) and aggregation
#' structure (\code{FishAgeComps_Type}, \code{FishLenComps_Type}) specified in
#' \code{sim_env}.
#'
#' Composition draws are skipped for fleet-season cells with zero fishing
#' mortality (\code{Fmort = 0}). Discard composition draws are additionally
#' skipped when retention selectivity is fully 1 for the fleet-region-year-
#' season cell (i.e., no discarding occurs). Discard indices support four
#' unit types: abundance (\code{discard_units = 0}), biomass (\code{1}),
#' abundance fraction (\code{2}), and biomass fraction (\code{3}).
#'
#' When \code{ISS_FishAgeComps_fill = "F_pattern"} and feedback is active,
#' sample sizes for retained and discard compositions in the current and
#' prior years are updated via \code{\link{predict_sim_fish_iss_fmort}}
#' (scaled by fishing mortality) before sampling.
#'
#' @param y Integer. Year index.
#'
#' @param sim Integer. Simulation replicate index.
#'
#' @param sim_env Simulation environment created by \code{\link{Setup_sim_env}}.
#'   Modified in place. The following elements are updated:
#'   \describe{
#'     \item{\code{CAA}, \code{DAA}}{Retained and dead discard catch-at-age for all populations, regions, seasons, and fleets.}
#'     \item{\code{CAL} and \code{DAL}}{Retained and dead discard catch-at-length if \code{SizeAgeTrans} is present.}
#'     \item{\code{TrueCatch}, \code{ObsCatch}}{Regional retained catch indices (abundance or biomass).}
#'     \item{\code{TrueCatch_pop}, \code{ObsCatch_pop}}{Population-specific retained catch indices.}
#'     \item{\code{TrueDiscard}, \code{ObsDiscard}}{Regional discard indices (abundance, biomass, or fraction).}
#'     \item{\code{TrueDiscard_pop}, \code{ObsDiscard_pop}}{Population-specific discard indices.}
#'     \item{\code{TrueFishIdx}, \code{ObsFishIdx}}{Regional fishery indices (abundance or biomass).}
#'     \item{\code{TrueFishIdx_pop}, \code{ObsFishIdx_pop}}{Population-specific fishery indices.}
#'     \item{\code{ObsFishAgeComps}, \code{ObsFishAgeComps_pop}}{Observed retained fishery age compositions.}
#'     \item{\code{ObsFishLenComps}, \code{ObsFishLenComps_pop}}{Observed retained fishery length compositions if \code{SizeAgeTrans} is available.}
#'     \item{\code{ObsFishAgeComps_discard}, \code{ObsFishAgeComps_discard_pop}}{Observed discard fishery age compositions.}
#'     \item{\code{ObsFishLenComps_discard}, \code{ObsFishLenComps_discard_pop}}{Observed discard fishery length compositions if \code{SizeAgeTrans} is available.}
#'     \item{\code{ISS_FishAgeComps}, \code{ISS_FishAgeComps_pop}, \code{ISS_FishLenComps}, \code{ISS_FishLenComps_pop}, \code{ISS_FishAgeComps_discard}, \code{ISS_FishAgeComps_discard_pop}, \code{ISS_FishLenComps_discard}, \code{ISS_FishLenComps_discard_pop}}{Effective sample sizes for retained and discard age and length compositions.}
#'   }
#'
#' @details For each combination of season, region, and fleet, the function:
#' \enumerate{
#'   \item Applies Baranov's catch equation to compute retained and dead discard catch-at-age.
#'   \item Converts catch-at-age to catch-at-length if \code{SizeAgeTrans} is available.
#'   \item Calculates true regional and population-specific catch, discard, and fishery indices.
#'   \item Applies lognormal observation error to generate observed indices.
#'   \item Simulates retained age and length compositions using \code{\link{simulate_comps}}, skipping fleet-season cells with zero fishing mortality.
#'   \item Simulates discard age and length compositions when retention selectivity is not fully 1.
#' }
#'
#' @return \code{invisible(NULL)}. All modifications are made by reference
#'   within \code{sim_env}.
#'
#' @keywords internal
generate_fishery_catch_comp_idx <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {
    for(seas in 1:n_seas) {
      for(r in 1:n_regions) {
        for(f in 1:n_fish_fleets) {

          for(p in 1:n_pop) {

            # Baranov's catch equation (retained catch-at-age)
            sim_env$CAA[p,r,y,seas,,,f,sim] <- (Fmort[r,y,seas,f,sim] * fish_sel[p,r,y,seas,,,f,sim] * ret_sel[p,r,y,seas,,,f,sim]) / ZAA[p,r,y,seas,,,sim] *
              NAA[p,r,y,seas,,,sim] * (1 - exp(-ZAA[p,r,y,seas,,,sim]))

            # Baranov's catch equation (dead discard catch-at-age)
            sim_env$DAA[p,r,y,seas,,,f,sim] <- (Fmort[r,y,seas,f,sim] * fish_sel[p,r,y,seas,,,f,sim] * (1 - ret_sel[p,r,y,seas,,,f,sim]) * dmr[r,y,seas,f,sim]) / ZAA[p,r,y,seas,,,sim] *
              NAA[p,r,y,seas,,,sim] * (1 - exp(-ZAA[p,r,y,seas,,,sim]))

            # Catch-at-length
            if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) for(s in 1:n_sexes) sim_env$CAL[p,r,y,seas,,s,f,sim] <- SizeAgeTrans[p,r,y,seas,,,s,sim] %*% CAA[p,r,y,seas,,s,f,sim] # Retained Catch at length
            if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) for(s in 1:n_sexes) sim_env$DAL[p,r,y,seas,,s,f,sim] <- SizeAgeTrans[p,r,y,seas,,,s,sim] %*% DAA[p,r,y,seas,,s,f,sim] # Discarded Catch at length

          } # end p loop


          # Regional Retained Catch
          if(catch_units[f] == 0) sim_env$TrueCatch[r,y,seas,f,sim] <- sum(CAA[,r,y,seas,,,f,sim]) # abundance
          if(catch_units[f] == 1) sim_env$TrueCatch[r,y,seas,f,sim] <- sum(CAA[,r,y,seas,,,f,sim] * WAA_fish[,r,y,seas,,,f,sim]) # biomass
          sim_env$ObsCatch[r,y,seas,f,sim] <- TrueCatch[r,y,seas,f,sim] * exp(stats::rnorm(1, 0, exp(ln_sigmaC[r,y,seas,f]))) # Observed Catch w/ lognormal deviations

          # Population Specific Catch
          if(catch_units[f] == 0) sim_env$TrueCatch_pop[,r,y,seas,f,sim] <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)  # abundance
          if(catch_units[f] == 1) sim_env$TrueCatch_pop[,r,y,seas,f,sim] <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE] * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)  # biomass
          sim_env$ObsCatch_pop[,r,y,seas,f,sim] <- sim_env$TrueCatch_pop[,r,y,seas,f,sim] * exp(stats::rnorm(n_pop, 0, exp(ln_sigmaC_pop[,r,y,seas,f])))

          # Regional Discards
          if(discard_units[f] == 0) sim_env$TrueDiscard[r,y,seas,f,sim] <- sum(DAA[,r,y,seas,,,f,sim]  / dmr[r,y,seas,f,sim]) # abd
          if(discard_units[f] == 1) sim_env$TrueDiscard[r,y,seas,f,sim] <- sum((DAA[,r,y,seas,,,f,sim]  / dmr[r,y,seas,f,sim]) * WAA_fish[,r,y,seas,,,f,sim]) # biom
          if(discard_units[f] == 2) { # abd frac
            total_catch <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            sim_env$TrueDiscard[r,y,seas,f,sim] <- 1 - sum(CAA[,r,y,seas,,,f,sim]) / sum(total_catch)
          }
          if(discard_units[f] == 3) { # biom frac
            total_catch <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            sim_env$TrueDiscard[r,y,seas,f,sim] <- 1 - sum(CAA[,r,y,seas,,,f,sim] * WAA_fish[,r,y,seas,,,f,sim]) / sum(total_catch * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE])
          }

          # lognormal
          sim_env$ObsDiscard[r,y,seas,f,sim] <- TrueDiscard[r,y,seas,f,sim] * exp(stats::rnorm(1, 0, exp(ln_sigmaD[r,y,seas,f])))

          # Population Specific Discards
          if(discard_units[f] == 0) sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- apply(DAA[,r,y,seas,,,f,sim, drop = FALSE]  / dmr[r,y,seas,f,sim], 1, sum) #abd
          if(discard_units[f] == 1) sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- apply((DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]) * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum) # biom
          if(discard_units[f] == 2) { # abd frac
            total_catch_pop <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            tmp_c <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)
            tmp_total <- apply(total_catch_pop, 1, sum)
            sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- 1 - tmp_c / tmp_total
          }
          if(discard_units[f] == 3) { # biom frac
            total_catch_pop <- CAA[,r,y,seas,,,f,sim, drop = FALSE] + DAA[,r,y,seas,,,f,sim, drop = FALSE] / dmr[r,y,seas,f,sim]
            tmp_c <- apply(CAA[,r,y,seas,,,f,sim, drop = FALSE] * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)
            tmp_total <- apply(total_catch_pop * WAA_fish[,r,y,seas,,,f,sim, drop = FALSE], 1, sum)
            sim_env$TrueDiscard_pop[,r,y,seas,f,sim] <- 1 - tmp_c / tmp_total
          }

          # Lognormal
          sim_env$ObsDiscard_pop[,r,y,seas,f,sim] <- sim_env$TrueDiscard_pop[,r,y,seas,f,sim] * exp(stats::rnorm(n_pop, 0, exp(ln_sigmaD_pop[,r,y,seas,f])))

          # Fishery Index
          tmp_expl_abd <- sweep(NAA[,r,y,seas,,,sim, drop = F], c(1,5,6), fish_sel[,r,y,seas,,,f,sim, drop = F] * ret_sel[,r,y,seas,,,f,sim, drop = F], "*")
          tmp_expl_biom <- sweep(tmp_expl_abd, c(1,5,6), WAA_fish[,r,y,seas,,,f,sim, drop = F], "*") # get exploitable abundance
          if(fish_idx_type[f] == 0) sim_env$TrueFishIdx[r,y,seas,f,sim] <- fish_q[r,y,f,sim] * sum(tmp_expl_abd) # True Fishery Index (abundance)
          if(fish_idx_type[f] == 1) sim_env$TrueFishIdx[r,y,seas,f,sim] <- fish_q[r,y,f,sim] * sum(tmp_expl_biom) # True Fishery Index (biomass)
          sim_env$ObsFishIdx[r,y,seas,f,sim] <- TrueFishIdx[r,y,seas,f,sim] * exp(stats::rnorm(1, 0, ObsFishIdx_SE[r,y,seas,f])) # Observed Fishery index w/ lognormal deviations

          # Population-specific Fishery Index
          if(fish_idx_type[f] == 0) sim_env$TrueFishIdx_pop[,r,y,seas,f,sim] <- fish_q[r,y,f,sim] * apply(tmp_expl_abd[,1,1,1,,,1, drop = FALSE], 1, sum)  # abundance
          if(fish_idx_type[f] == 1) sim_env$TrueFishIdx_pop[,r,y,seas,f,sim] <- fish_q[r,y,f,sim] * apply(tmp_expl_biom[,1,1,1,,,1, drop = FALSE], 1, sum)  # biomass
          sim_env$ObsFishIdx_pop[,r,y,seas,f,sim] <- sim_env$TrueFishIdx_pop[,r,y,seas,f,sim] * exp(stats::rnorm(n_pop, 0, ObsFishIdx_pop_SE[,r,y,seas,f]))

          # Fishery Compositions
          if(Fmort[r,y,seas,f,sim] > 0) { # only simulate if Fishing Mortality > 0

            # Retained Compositions
            # Age Compositions (Dynamic ISS based on feedback fishing mortality)
            if(exists("ISS_FishAgeComps_fill") && isTRUE(ISS_FishAgeComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              sim_env$ISS_FishAgeComps[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishAgeComps, Fmort = Fmort, y = y, sim = sim, seas = seas)
            }
            if(exists("ISS_FishAgeComps_pop_fill") && isTRUE(ISS_FishAgeComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              for(p in 1:n_pop) sim_env$ISS_FishAgeComps_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = array(ISS_FishAgeComps_pop[p,,1:y,,,,sim],
                                                                                                                                    dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                                                                                                              Fmort = Fmort, y = y, sim = sim, seas = seas)
            }

            # Length Compositions (Dynamic ISS based on feedback fishing mortality)
            if(exists("ISS_FishLenComps_fill") && isTRUE(ISS_FishLenComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              sim_env$ISS_FishLenComps[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishLenComps, Fmort = Fmort, y = y, sim = sim, seas = seas)
            }
            if(exists("ISS_FishLenComps_pop_fill") && isTRUE(ISS_FishLenComps_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
              for(p in 1:n_pop) sim_env$ISS_FishLenComps_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = array(ISS_FishLenComps_pop[p,,1:y,,,,sim],
                                                                                                                                    dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                                                                                                              Fmort = Fmort, y = y, sim = sim, seas = seas)
            }

            # Sample fishery ages (non-population specific, retained compositions)
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
                                                      pop_specific = FALSE,
                                                      age_or_len = 0)

            # Sample fishery ages (population specific, retained compositions)
            sim_env$ObsFishAgeComps_pop <- simulate_comps(r = r,
                                                          y = y,
                                                          seas = seas,
                                                          f = f,
                                                          sim = sim,
                                                          Exp = CAA,
                                                          ISS_pop = ISS_FishAgeComps_pop,
                                                          AgeingError = AgeingError,
                                                          pop_comp_like = comp_fishage_pop_like,
                                                          ln_pop_theta = ln_FishAge_pop_theta,
                                                          ln_pop_theta_agg = ln_FishAge_pop_theta_agg,
                                                          pop_corr_pars = FishAge_pop_corr_pars,
                                                          pop_corr_pars_agg = FishAge_pop_corr_pars_agg,
                                                          pop_comp_type = FishAgeComps_pop_Type,
                                                          n_sexes = n_sexes,
                                                          n_regions = n_regions,
                                                          n_pop = n_pop,
                                                          n_cat = n_ages,
                                                          Obs = ObsFishAgeComps_pop,
                                                          pop_specific = TRUE,
                                                          age_or_len = 0)

            # Sample fishery lengths (retained compositions)
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
                                                        pop_specific = FALSE,
                                                        age_or_len = 1)

              # Sample fishery lengths (population specific)
              sim_env$ObsFishLenComps_pop <- simulate_comps(r = r,
                                                            y = y,
                                                            seas = seas,
                                                            f = f,
                                                            sim = sim,
                                                            Exp = CAL,
                                                            ISS_pop = ISS_FishLenComps_pop,
                                                            AgeingError = NULL,
                                                            pop_comp_like = comp_fishlen_pop_like,
                                                            ln_pop_theta = ln_FishLen_pop_theta,
                                                            ln_pop_theta_agg = ln_FishLen_pop_theta_agg,
                                                            pop_corr_pars = FishLen_pop_corr_pars,
                                                            pop_corr_pars_agg = FishLen_pop_corr_pars_agg,
                                                            pop_comp_type = FishLenComps_pop_Type,
                                                            n_sexes = n_sexes,
                                                            n_regions = n_regions,
                                                            n_pop = n_pop,
                                                            n_cat = n_lens,
                                                            Obs = ObsFishLenComps_pop,
                                                            pop_specific = TRUE,
                                                            age_or_len = 1)

            } # end if size age transition if availiable

            # if there is discarding occuring
            if(!all(DAA[,r,y,seas,,,f,sim] == 0)) {

              # Discarded Compositions (Dynamic ISS based on feedback fishing mortality)
              if(exists("ISS_FishAgeComps_discard_fill") && isTRUE(ISS_FishAgeComps_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                sim_env$ISS_FishAgeComps_discard[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishAgeComps_discard, Fmort = Fmort, y = y, sim = sim, seas = seas)
              }
              if(exists("ISS_FishAgeComps_pop_discard_fill") && isTRUE(ISS_FishAgeComps_pop_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                for(p in 1:n_pop) sim_env$ISS_FishAgeComps_discard_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = array(ISS_FishAgeComps_discard_pop[p,,1:y,,,,sim],
                                                                                                                                              dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                                                                                                                        Fmort = Fmort, y = y, sim = sim, seas = seas)
              }

              # Length Compositions (Dynamic ISS based on feedback fishing mortality)
              if(exists("ISS_FishLenComps_discard_fill") && isTRUE(ISS_FishLenComps_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                sim_env$ISS_FishLenComps_discard[,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = ISS_FishLenComps_discard, Fmort = Fmort, y = y, sim = sim, seas = seas)
              }
              if(exists("ISS_FishLenComps_pop_discard_fill") && isTRUE(ISS_FishLenComps_pop_discard_fill == "F_pattern") && isTRUE(run_feedback) && y >= feedback_start_yr + 1 && r == 1 && f == 1) {
                for(p in 1:n_pop) sim_env$ISS_FishLenComps_discard_pop[p,,1:y,seas,,,sim] <- predict_sim_fish_iss_fmort(ISS_FishComps = array(ISS_FishLenComps_discard_pop[p,,1:y,,,,sim],
                                                                                                                                              dim = c(n_regions, length(1:y), n_seas, n_sexes, n_fish_fleets, n_sims)),
                                                                                                                        Fmort = Fmort, y = y, sim = sim, seas = seas)
              }

              # Sample fishery lengths (non-population specific, discard compositions)
              sim_env$ObsFishAgeComps_discard <- simulate_comps(r = r,
                                                                y = y,
                                                                seas = seas,
                                                                f = f,
                                                                sim = sim,
                                                                Exp = DAA,
                                                                ISS = ISS_FishAgeComps_discard,
                                                                AgeingError = AgeingError,
                                                                comp_like = comp_fishage_discard_like,
                                                                ln_theta = ln_FishAge_discard_theta,
                                                                ln_theta_agg = ln_FishAge_discard_theta_agg,
                                                                corr_pars = FishAge_discard_corr_pars,
                                                                corr_pars_agg = FishAge_discard_corr_pars_agg,
                                                                comp_type = FishAgeComps_discard_Type,
                                                                n_sexes = n_sexes,
                                                                n_regions = n_regions,
                                                                n_cat = n_ages,
                                                                Obs = ObsFishAgeComps_discard,
                                                                pop_specific = FALSE,
                                                                age_or_len = 0)

              # Sample fishery ages (population specific, discard compositions)
              sim_env$ObsFishAgeComps_discard_pop <- simulate_comps(r = r,
                                                                    y = y,
                                                                    seas = seas,
                                                                    f = f,
                                                                    sim = sim,
                                                                    Exp = DAA,
                                                                    ISS_pop = ISS_FishAgeComps_discard_pop,
                                                                    AgeingError = AgeingError,
                                                                    pop_comp_like = comp_fishage_discard_pop_like,
                                                                    ln_pop_theta = ln_FishAge_discard_pop_theta,
                                                                    ln_pop_theta_agg = ln_FishAge_discard_pop_theta_agg,
                                                                    pop_corr_pars = FishAge_discard_pop_corr_pars,
                                                                    pop_corr_pars_agg = FishAge_discard_pop_corr_pars_agg,
                                                                    pop_comp_type = FishAgeComps_discard_pop_Type,
                                                                    n_sexes = n_sexes,
                                                                    n_regions = n_regions,
                                                                    n_pop = n_pop,
                                                                    n_cat = n_ages,
                                                                    Obs = ObsFishAgeComps_discard_pop,
                                                                    pop_specific = TRUE,
                                                                    age_or_len = 0)

              # Sample fishery lengths (retained compositions)
              if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) {
                # Sample fishery ages (non-population specific, discard compositions)
                sim_env$ObsFishLenComps_discard <- simulate_comps(r = r,
                                                                  y = y,
                                                                  seas = seas,
                                                                  f = f,
                                                                  sim = sim,
                                                                  Exp = DAL,
                                                                  ISS = ISS_FishLenComps_discard,
                                                                  AgeingError = NULL,
                                                                  comp_like = comp_fishlen_discard_like,
                                                                  ln_theta = ln_FishLen_discard_theta,
                                                                  ln_theta_agg = ln_FishLen_discard_theta_agg,
                                                                  corr_pars = FishLen_discard_corr_pars,
                                                                  corr_pars_agg = FishLen_discard_corr_pars_agg,
                                                                  comp_type = FishLenComps_discard_Type,
                                                                  n_sexes = n_sexes,
                                                                  n_regions = n_regions,
                                                                  n_cat = n_lens,
                                                                  Obs = ObsFishLenComps_discard,
                                                                  pop_specific = FALSE,
                                                                  age_or_len = 1)

                # Sample fishery lengths (population specific, discard compositions)
                sim_env$ObsFishLenComps_discard_pop <- simulate_comps(r = r,
                                                                      y = y,
                                                                      seas = seas,
                                                                      f = f,
                                                                      sim = sim,
                                                                      Exp = DAL,
                                                                      ISS_pop = ISS_FishLenComps_discard_pop,
                                                                      AgeingError = NULL,
                                                                      pop_comp_like = comp_fishlen_discard_pop_like,
                                                                      ln_pop_theta = ln_FishLen_discard_pop_theta,
                                                                      ln_pop_theta_agg = ln_FishLen_discard_pop_theta_agg,
                                                                      pop_corr_pars = FishLen_discard_pop_corr_pars,
                                                                      pop_corr_pars_agg = FishLen_discard_pop_corr_pars_agg,
                                                                      pop_comp_type = FishLenComps_discard_pop_Type,
                                                                      n_sexes = n_sexes,
                                                                      n_regions = n_regions,
                                                                      n_pop = n_pop,
                                                                      n_cat = n_lens,
                                                                      Obs = ObsFishLenComps_discard_pop,
                                                                      pop_specific = TRUE,
                                                                      age_or_len = 1)

              } # end if size age transition if availiable

            } # end if dmr > 0
          } # end if Fmort > 0

        } # end f loop
      } # end r loop
    } # end seas loop
  })

}


#' Generate survey indices and compositions in a simulation
#'
#' Computes survey index-at-age (\code{SrvIAA}) for all populations using
#' the mid-survey abundance formula \eqn{N \cdot s \cdot e^{-t_{\text{srv}} Z}},
#' derives index-at-length (\code{SrvIAL}) when a size-age transition matrix
#' is available, generates observed survey indices (with lognormal error) as
#' abundance or biomass depending on \code{srv_idx_type}, and draws age and
#' length composition samples via \code{\link{simulate_comps}}.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by \code{\link{Setup_sim_env}}.
#'   Modified in place. The following elements are updated:
#'   \describe{
#'     \item{\code{SrvIAA}}{Survey index-at-age for all populations.}
#'     \item{\code{SrvIAL}}{Survey index-at-length if \code{SizeAgeTrans} is present.}
#'     \item{\code{TrueSrvIdx}, \code{ObsSrvIdx}}{Aggregated survey index values.}
#'     \item{\code{TrueSrvIdx_pop}, \code{ObsSrvIdx_pop}}{Population-specific survey index values.}
#'     \item{\code{ObsSrvAgeComps}, \code{ObsSrvAgeComps_pop}}{Observed survey age compositions.}
#'     \item{\code{ObsSrvLenComps}, \code{ObsSrvLenComps_pop}}{Observed survey length compositions if \code{SizeAgeTrans} is available.}
#'   }
#'
#' @details This function loops over seasons, regions, and survey fleets for all
#' populations and replicates. It computes mid-period abundance, applies survey
#' selectivity, calculates true survey indices (abundance or biomass), applies
#' lognormal observation error, and simulates age and length composition samples.
#' Population-specific compositions are also generated if requested.
#'
#' @return \code{invisible(NULL)}. All modifications are performed by reference
#'   within \code{sim_env}.
#'
#' @keywords internal
#' @seealso \code{\link{simulate_comps}}, \code{\link{Setup_sim_env}}
generate_survey_comp_idx <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {

    for(seas in 1:n_seas) {
      for(r in 1:n_regions) {
        for(sf in 1:n_srv_fleets) {

          for(p in 1:n_pop) {
            # Survey Ages Indexed (midpoint year)
            sim_env$SrvIAA[p,r,y,seas,,,sf,sim] <- NAA[p,r,y,seas,,,sim] * srv_sel[p,r,y,seas,,,sf,sim] * exp(-t_srv[r,seas,sf] * ZAA[p,r,y,seas,,,sim])
            if(exists("SizeAgeTrans") && !is.null(SizeAgeTrans)) for(s in 1:n_sexes) sim_env$SrvIAL[p,r,y,seas,,s,sf,sim] <- SizeAgeTrans[p,r,y,seas,,,s,sim] %*% SrvIAA[p,r,y,seas,,s,sf,sim] # Survey index at length
          } # end p loop

          # Survey Index - Regional
          if(srv_idx_type[sf] == 0) sim_env$TrueSrvIdx[r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[,r,y,seas,,,sf,sim]) # True Survey Index (abundance)
          if(srv_idx_type[sf] == 1) sim_env$TrueSrvIdx[r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * sum(SrvIAA[,r,y,seas,,,sf,sim] * WAA_srv[,r,y,seas,,,sf,sim]) # True Survey Index (biomass)
          sim_env$ObsSrvIdx[r,y,seas,sf,sim] <- TrueSrvIdx[r,y,seas,sf,sim] * exp(stats::rnorm(1, 0, ObsSrvIdx_SE[r,y,seas,sf])) # Observed survey index w/ lognormal deviations

          # Survey Index - Population-Specific
          if(srv_idx_type[sf] == 0) sim_env$TrueSrvIdx_pop[,r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * apply(SrvIAA[,r,y,seas,,,sf,sim, drop = FALSE], 1, sum) # True Survey Index (abundance)
          if(srv_idx_type[sf] == 1) sim_env$TrueSrvIdx_pop[,r,y,seas,sf,sim] <- srv_q[r,y,sf,sim] * apply(SrvIAA[,r,y,seas,,,sf,sim, drop = FALSE] * WAA_srv[,r,y,seas,,,sf,sim, drop = FALSE], 1, sum) # True Survey Index (biomass)
          sim_env$ObsSrvIdx_pop[,r,y,seas,sf,sim] <- TrueSrvIdx_pop[,r,y,seas,sf,sim] * exp(stats::rnorm(n_pop, 0, ObsSrvIdx_pop_SE[,r,y,seas,sf])) # Observed survey index w/ lognormal deviations

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
                                                   age_or_len = 0)

          # Sample survey ages (population specific)
          sim_env$ObsSrvAgeComps_pop <- simulate_comps(r = r,
                                                       y = y,
                                                       seas = seas,
                                                       f = sf,
                                                       sim = sim,
                                                       Exp = SrvIAA,
                                                       ISS_pop = ISS_SrvAgeComps_pop,
                                                       AgeingError = AgeingError,
                                                       pop_comp_like = comp_srvage_pop_like,
                                                       ln_pop_theta = ln_SrvAge_pop_theta,
                                                       ln_pop_theta_agg = ln_SrvAge_pop_theta_agg,
                                                       pop_corr_pars = SrvAge_pop_corr_pars,
                                                       pop_corr_pars_agg = SrvAge_pop_corr_pars_agg,
                                                       pop_comp_type = SrvAgeComps_pop_Type,
                                                       n_sexes = n_sexes,
                                                       n_regions = n_regions,
                                                       n_pop = n_pop,
                                                       n_cat = n_ages,
                                                       Obs = ObsSrvAgeComps_pop,
                                                       pop_specific = TRUE,
                                                       age_or_len = 0)

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
                                                     age_or_len = 1)

            # Sample survey lengths (population specific)
            sim_env$ObsSrvLenComps_pop <- simulate_comps(r = r,
                                                         y = y,
                                                         seas = seas,
                                                         f = sf,
                                                         sim = sim,
                                                         Exp = SrvIAL,
                                                         ISS_pop = ISS_SrvLenComps_pop,
                                                         AgeingError = NULL,
                                                         pop_comp_like = comp_srvlen_pop_like,
                                                         ln_pop_theta = ln_SrvLen_pop_theta,
                                                         ln_pop_theta_agg = ln_SrvLen_pop_theta_agg,
                                                         pop_corr_pars = SrvLen_pop_corr_pars,
                                                         pop_corr_pars_agg = SrvLen_pop_corr_pars_agg,
                                                         pop_comp_type = SrvLenComps_pop_Type,
                                                         n_sexes = n_sexes,
                                                         n_regions = n_regions,
                                                         n_pop = n_pop,
                                                         n_cat = n_lens,
                                                         Obs = ObsSrvLenComps_pop,
                                                         pop_specific = TRUE,
                                                         age_or_len = 1)

          } # end if size age transition if availiable

        } # end sf loop
      } # end r loop
    } # end seas loop

  })
}

#' Release conventional tags in the simulation
#'
#' For each tag cohort scheduled for release in year \code{y}, distributes
#' \code{n_tags} (or \code{n_tags_rel_input} if provided) across populations,
#' ages, and sexes proportional to the selectivity-weighted abundance
#' (\code{NAA_bef}) of the release platform (survey, fishery, or population).
#' Tagged fish counts are rounded to integers. The attended attribute string
#' \code{conv_fish_tag_attr} is then applied via
#' \code{\link{marginalize_conv_fish_tags}} to produce the observation-level
#' release array \code{conv_tagged_fish_attr}, which is consistent with the
#' dimension resolution of the recapture likelihood.
#'
#' For survey and fishery platforms, total tags in the release region are
#' scaled relative to the selectivity-weighted global abundance to allocate
#' region-specific cohort sizes when \code{n_tags_rel_input} is not provided.
#' For the population platform, scaling is proportional to the region's share
#' of total abundance.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by
#'   \code{\link{Setup_sim_env}}. Modified in place:
#'   \code{$conv_tagged_fish[tc, , , , sim]} and
#'   \code{$conv_tagged_fish_attr[tc, , , , sim]} for each cohort released
#'   in year \code{y}.
#'
#' @return \code{invisible(NULL)}. All modifications are made by reference
#'   within \code{sim_env}.
#'
#'
#' @keywords internal
release_conv_tags <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {
    for(seas in 1:n_seas) {
      for(r in 1:n_regions) {

        # Get indices for tag cohorts in the current year and region
        tag_rel <- which(conv_tag_release_indicator[,1] == r & conv_tag_release_indicator[,2] == y & conv_tag_release_indicator[,3] == seas) # Get tag cohort (release event)

        # Release Tags if any events
        if(length(tag_rel) != 0) {

          # Tag Indexing
          tr <- conv_tag_release_indicator[tag_rel,1] # tag release region
          ty <- conv_tag_release_indicator[tag_rel,2] # tag release year
          tseas <- conv_tag_release_indicator[tag_rel,3] # tag release season
          tplat <- conv_tag_release_platform[tag_rel, 1] # get tagging platform information
          tplat_f <- as.numeric(conv_tag_release_platform[tag_rel, 2]) # get tagging fleet

          # Survey
          # Survey
          if(tplat[1] == 'survey') {
            NAA_slice  <- NAA_bef[, tr, ty, tseas, , , sim, drop = FALSE]
            dim(NAA_slice) <- c(n_pop, n_ages, n_sexes)
            NAA_sel <- array(0, dim = c(n_pop, n_ages, n_sexes))
            for(p in 1:n_pop) {
              sel_slice <- array(srv_sel[p, tr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
              NAA_sel[p,,] <- NAA_slice[p,,] * sel_slice
            }
            if(!exists("n_tags_rel_input")) {
              denom <- 0
              for(p in 1:n_pop) {
                NAA_denom_p <- array(NAA_bef[p, , ty, tseas, , , sim], dim = c(n_regions, n_ages, n_sexes))
                for(rr in 1:n_regions) {
                  sel_d <- array(srv_sel[p, rr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
                  denom <- denom + sum(NAA_denom_p[rr,,] * sel_d)
                }
              }
              n_tags_rel <- round(sum(NAA_sel) / denom * n_tags)
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel]
            }
            tmp_props <- NAA_sel / sum(NAA_sel)
          }

          # distribute tags for fishery
          if(tplat[1] == 'fishery') {
            NAA_slice  <- NAA_bef[, tr, ty, tseas, , , sim, drop = FALSE]
            dim(NAA_slice) <- c(n_pop, n_ages, n_sexes)
            NAA_sel <- array(0, dim = c(n_pop, n_ages, n_sexes))
            for(p in 1:n_pop) {
              sel_slice <- array(fish_sel[p, tr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
              NAA_sel[p,,] <- NAA_slice[p,,] * sel_slice
            }
            if(!exists("n_tags_rel_input")) {
              denom <- 0
              for(p in 1:n_pop) {
                NAA_denom_p <- array(NAA_bef[p, , ty, tseas, , , sim], dim = c(n_regions, n_ages, n_sexes))
                for(rr in 1:n_regions) {
                  sel_d <- array(fish_sel[p, rr, ty, tseas, , , tplat_f, sim], dim = c(n_ages, n_sexes))
                  denom <- denom + sum(NAA_denom_p[rr,,] * sel_d)
                }
              }
              n_tags_rel <- round(sum(NAA_sel) / denom * n_tags)
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel]
            }
            tmp_props <- NAA_sel / sum(NAA_sel)
          }

          # distribute tags by population
          if(tplat[1] == 'population') {
            if(!exists("n_tags_rel_input")) {
              n_tags_rel <- round(sum(NAA_bef[,tr,ty,tseas,,,sim]) / sum(NAA_bef[,,ty,tseas,,,sim]) * n_tags) # get region specific tags
            } else {
              n_tags_rel <- n_tags_rel_input[tag_rel] # use input tags by cohort if availiable
            }
            tmp_props <- NAA_bef[, tr, ty, tseas, , , sim] / sum(NAA_bef[, tr, ty, tseas, , , sim]) # get proportions by population, age, and sex
          }

          # multiply by tags in each region and distributa cross ages and sexes
          sim_env$conv_tagged_fish[tag_rel, , , , sim] <- array(round(tmp_props * n_tags_rel), dim = c(n_pop, n_ages, n_sexes))

          # marginalize across appropriate dimensions of what recaptures attributes are and report out
          sim_env$conv_tagged_fish_attr[tag_rel, , , , sim] <- marginalize_conv_fish_tags(conv_tagged_fish[tag_rel, , , , sim],
                                                                                              conv_fish_tag_attr, n_pop, n_ages, n_sexes)

        } # end if no tag releases

      } # end r loop
    } # end seas loop
  })
}

#' Generate conventional tag recaptures from fisheries in simulation
#'
#' For each tag cohort (\code{tc}) and recovery season (\code{rseas}) in year
#' \code{y}, advances available tagged fish through movement and mortality,
#' applies Baranov's equation to compute predicted recaptures
#' (\code{pred_conv_tag_fish_recap}), and draws observed recaptures
#' (\code{obs_conv_tag_fish_recap}) via
#' \code{\link{simulate_conv_tag_fish_recaptures}}. Cohorts not yet released,
#' already at \code{conv_tag_max_liberty}, or with release year in the future
#' are silently skipped.
#'
#' Total fishing mortality entering Z is decomposed into retained
#' (\eqn{F \cdot s_{\text{fish}} \cdot s_{\text{ret}}}) and dead discard
#' (\eqn{F \cdot s_{\text{fish}} \cdot (1 - s_{\text{ret}}) \cdot \text{dmr}})
#' components, consistent with \code{\link{apply_pop_dy}}. Predicted
#' recaptures use only the retained component in the Baranov numerator,
#' reflecting that tags are recovered from retained catch only.
#'
#' At initial release (\code{ry = 1}, \code{rseas = tseas}), tags are
#' placed into \code{conv_tag_fish_avail[1, rseas, tc, ...]} after discounting
#' for initial tag-induced mortality (\code{ln_init_conv_tag_mort}). When
#' \code{conv_tag_t_tagging < 1}, total mortality is scaled by the fraction of
#' the season remaining at release for that cell only. Chronic shedding
#' (\code{ln_conv_tag_shed}) enters the total mortality rate alongside natural
#' and fishing mortality. At the end of each season, survivors advance to the
#' next season or the next year's first season with plus-group accumulation.
#' Tag reporting rates from \code{conv_tag_fish_reporting} are applied fleet-
#' and region-specifically.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by \code{\link{Setup_sim_env}}.
#'   Modified in place. The following elements are updated:
#'   \describe{
#'     \item{\code{conv_tag_fish_avail}}{Available tagged fish at age/region/season/fleet for each cohort.}
#'     \item{\code{pred_conv_tag_fish_recap}}{Predicted conventional tag recaptures by fleet, region, season, and cohort.}
#'     \item{\code{obs_conv_tag_fish_recap}}{Observed conventional tag recaptures after sampling error.}
#'     \item{\code{conv_tag_fish_surv}}{Surviving tagged fish after mortality and movement.}
#'     \item{\code{conv_tag_fish_reported}}{Reporting-adjusted recapture counts by fleet and region.}
#'   }
#'
#' @return \code{invisible(NULL)}. All modifications are made by reference
#'   within \code{sim_env}.
#'
#' @details Tagged fish dynamics follow the same seasonal progression logic
#' as the population projection, including natural mortality, fishing mortality,
#' movement, and tag shedding. Recaptures are computed only from the retained
#' catch component, consistent with tag return processes. Cohorts exceeding
#' \code{conv_tag_max_liberty} are removed from the active tracking pool.
#'
#' @keywords internal
generate_fishery_conv_tags_recap <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env,{

      for(rseas in 1:n_seas) {
        for(tc in 1:n_tag_rel_events) {

          # get indexing
          tr <- conv_tag_release_indicator[tc,1] # tag release region
          ty <- conv_tag_release_indicator[tc,2] # tag release year
          tseas <- conv_tag_release_indicator[tc,3] # tag release seasons

          # Skipping stuff if hasn't occurred yet, or if max liberty
          if(y < ty || (y == ty && rseas < tseas)) next
          ry <- y - ty + 1 # get tag liberty
          if(ry > conv_tag_max_liberty) next # skip if max liberty

          # get fishing mortality
          tmp_FAA = array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
          tmp_ret_FAA = array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
          tmp_disc_DAA = array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
          for(p in 1:n_pop) for(f in 1:n_fish_fleets) {
            if(use_conv_fish_tagging[f] == 1) {
              tmp_ret_FAA[p,,1,,,f] = Fmort[, y, rseas, f, sim] * fish_sel[p,,y,rseas,,,f,sim] * ret_sel[p,,y,rseas,,,f,sim]  # Retained fishing mortality
              tmp_disc_DAA[p,,1,,,f] = Fmort[, y, rseas, f, sim] * fish_sel[p,,y,rseas,,,f,sim] * (1 - ret_sel[p,,y,rseas,,,f,sim]) * dmr[,y,rseas,f,sim] # Dead discard fishing mortality
              tmp_FAA[p,,1,,,f] = tmp_ret_FAA[p,,1,,,f] + tmp_disc_DAA[p,,1,,,f] # Total fishing mortality
            } # end if
          } # end p loop

          # get total mortality
          tmp_natmort = array(natmort[,,y,,,sim], dim = c(n_pop, n_regions, 1, n_ages, n_sexes))
          tmp_ZAA = (tmp_natmort * seasdur[rseas]) + apply(tmp_FAA, 1:5, sum) + (exp(ln_conv_tag_shed) * seasdur[rseas])

          # Discount with tagging time (conv_tag_t_tagging) if it doesn't happen at the start of the season / year
          if(ry == 1 && rseas == tseas) {
            if(conv_tag_t_tagging != 1) tmp_ZAA <- tmp_ZAA * conv_tag_t_tagging
            # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
            sim_env$conv_tag_fish_avail[1, rseas, tc, , tr, , , sim] <- array(conv_tagged_fish[tc, , , , sim] * exp(-exp(ln_init_conv_tag_mort)), dim = c(n_pop, n_ages, n_sexes))
          }

          # get temporary survival value
          tmp_SAA <- exp(-tmp_ZAA)

          # Move tagged fish around (skip only in first release year + tagging season when tagging occurs mid-season)
          if(conv_tag_t_tagging == 1 || ry != 1 || rseas != tseas) {
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
              tmp_SAA[,,1,,]

          } else {

            # End of year mortality and age advancement (end of season)
            sim_env$conv_tag_fish_avail[ry + 1, 1, tc, , , 2:n_ages, , sim] <-
              conv_tag_fish_avail[ry, n_seas, tc, , , 1:(n_ages-1), , sim] *
              tmp_SAA[,,1,1:(n_ages - 1),]

            # Accumulate plus group
            sim_env$conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, , sim] <-
              conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, , sim] +
              conv_tag_fish_avail[ry, n_seas, tc, , , n_ages, , sim] *
              tmp_SAA[,,1,n_ages,]
          }

          # # Apply Baranov's to get predicted recaptures
          for(f in 1:n_fish_fleets) {
            for(p in 1:n_pop) {
              sim_env$pred_conv_tag_fish_recap[ry,rseas,tc,p,,,,f,sim] <- conv_tag_fish_reporting[,y,f,sim] *
                (tmp_ret_FAA[p,,1,,,f] / tmp_ZAA[p,,1,,]) *
                conv_tag_fish_avail[ry,rseas,tc,p,,,,sim] *
                (1 - tmp_SAA[p,,1,,])
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


#' Run the annual cycle for a single simulation year
#'
#' Orchestrates the complete annual sequence of operating model processes for
#' year \code{y} and simulation replicate \code{sim}: initialises age
#' structure and generates first-year recruitment at \code{y = 1}; applies
#' population dynamics (movement, mortality, biomass); generates fishery
#' catches, indices, and compositions; generates survey indices and
#' compositions; releases conventional tags; generates fishery tag
#' recaptures (when any \code{use_conv_fish_tagging = 1}); and generates
#' recruitment for the following year (\code{y + 1}) when \code{y < n_yrs}.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment created by
#'   \code{\link{Setup_sim_env}} and passed by reference. All annual-cycle
#'   helper functions modify this environment in place.
#'
#' @return \code{invisible(NULL)}.
#'
#'
#' @importFrom stats rnorm rmultinom
#' @export run_annual_cycle
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

  if(any(sim_env$use_conv_fish_tagging == 1)) {
    release_conv_tags(y, sim, sim_env) # Release conventional tags
    generate_fishery_conv_tags_recap(y, sim, sim_env) # Generate fishery conventional tag recaptures
  }

  if(y < sim_env$n_yrs) generate_recruitment(y = y + 1, sim, sim_env) # Get recruitment in the following year


  return(invisible(NULL))

}

#' Simulate a static (open-loop) spatial age- and sex-structured population
#'
#' Runs a complete multi-replicate operating model simulation with no
#' feedback between the population and the harvest control rule (i.e.,
#' fishing mortality is fixed as supplied in \code{sim_list}). Calls
#' \code{\link{Setup_sim_env}} to create an isolated execution environment
#' and then iterates \code{\link{run_annual_cycle}} over all years and
#' simulation replicates. All simulation outputs are collected from the
#' environment and returned as a named list. Optionally writes the output
#' to an RDS file.
#'
#' @param sim_list Simulation list returned by the last upstream setup
#'   function (typically \code{\link{Setup_Sim_Rec}} or
#'   \code{\link{Setup_Sim_Tagging}}).
#' @param output_path Character string. File path for saving the output list
#'   as an RDS file via \code{saveRDS}. If \code{NULL} (default), no file is
#'   written.
#'
#' @return A named list containing all simulation outputs, including (among
#'   others): \code{NAA}, \code{NAA0}, \code{SSB}, \code{Dynamic_SSB0},
#'   \code{eff_SSB}, \code{Rec}, \code{ln_RecDevs}, \code{ln_InitDevs},
#'   \code{ZAA}, \code{TrueCatch}, \code{ObsCatch}, \code{TrueCatch_pop},
#'   \code{ObsCatch_pop}, \code{CAA}, \code{CAL},
#'   \code{ObsFishAgeComps}, \code{ObsFishAgeComps_pop},
#'   \code{ObsFishLenComps}, \code{ObsFishLenComps_pop},
#'   \code{ObsFishIdx}, \code{TrueFishIdx},
#'   \code{ObsFishIdx_pop}, \code{TrueFishIdx_pop},
#'   \code{SrvIAA}, \code{SrvIAL},
#'   \code{ObsSrvAgeComps}, \code{ObsSrvAgeComps_pop},
#'   \code{ObsSrvLenComps}, \code{ObsSrvLenComps_pop},
#'   \code{ObsSrvIdx}, \code{TrueSrvIdx},
#'   \code{ObsSrvIdx_pop}, \code{TrueSrvIdx_pop},
#'   \code{conv_tagged_fish}, \code{conv_tagged_fish_attr},
#'   \code{conv_tag_fish_avail}, \code{pred_conv_tag_fish_recap},
#'   \code{obs_conv_tag_fish_recap}, and key dimension scalars
#'   (\code{n_regions}, \code{n_pop}, \code{n_yrs}, \code{n_ages}, etc.).
#'   Note that \code{n_years} and \code{n_yrs} are both present for backwards
#'   compatibility.
#'
#'
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
                  dmr = sim_env$dmr,
                  ln_sigmaC = sim_env$ln_sigmaC,
                  ln_sigmaC_pop = sim_env$ln_sigmaC_pop,
                  ln_sigmaD = sim_env$ln_sigmaD,
                  ln_sigmaD_pop = sim_env$ln_sigmaD_pop,
                  fish_sel = sim_env$fish_sel,
                  ret_sel = sim_env$ret_sel,
                  fish_q = sim_env$fish_q,
                  ln_RecDevs = sim_env$ln_RecDevs,
                  ln_InitDevs = sim_env$ln_InitDevs,
                  natmort = sim_env$natmort,
                  ZAA = sim_env$ZAA,
                  sexratio = sim_env$sexratio,
                  R0 = sim_env$R0,
                  Rec = sim_env$Rec,
                  natal_region = sim_env$natal_region,
                  WAA = sim_env$WAA,
                  rec_seas_prop = sim_env$rec_seas_prop,
                  WAA_fish = sim_env$WAA_fish,
                  WAA_srv = sim_env$WAA_srv,
                  MatAA = sim_env$MatAA,
                  h = sim_env$h,
                  do_recruits_move = sim_env$do_recruits_move,
                  ln_sigmaR = sim_env$ln_sigmaR,
                  Movement = sim_env$Movement,
                  sgl_seas_spawning_movement = sim_env$sgl_seas_spawning_movement,
                  NAA = sim_env$NAA,
                  NAA_bef = sim_env$NAA_bef,
                  NAA_aft = sim_env$NAA_aft,
                  NAA0 = sim_env$NAA0,
                  Dynamic_SSB0 = sim_env$Dynamic_SSB0,
                  SSB = sim_env$SSB,
                  eff_SSB = sim_env$eff_SSB,
                  stray_rate = sim_env$stray_rate,
                  t_spawn = sim_env$t_spawn,
                  Total_Biom = sim_env$Total_Biom,

                  # Aggregated fishery obs
                  TrueCatch = sim_env$TrueCatch,
                  ObsCatch = sim_env$ObsCatch,
                  ObsFishIdx = sim_env$ObsFishIdx,
                  TrueFishIdx = sim_env$TrueFishIdx,
                  ObsFishIdx_SE = sim_env$ObsFishIdx_SE,
                  ObsFishAgeComps = sim_env$ObsFishAgeComps,
                  ObsFishLenComps = sim_env$ObsFishLenComps,

                  # Aggregated fishery discards
                  TrueDiscard = sim_env$TrueDiscard,
                  ObsDiscard = sim_env$ObsDiscard,
                  ObsFishAgeComps_discard = sim_env$ObsFishAgeComps_discard,
                  ObsFishLenComps_discard = sim_env$ObsFishLenComps_discard,

                  # Population-specific fishery obs
                  TrueCatch_pop = sim_env$TrueCatch_pop,
                  ObsCatch_pop = sim_env$ObsCatch_pop,
                  ObsFishIdx_pop = sim_env$ObsFishIdx_pop,
                  TrueFishIdx_pop = sim_env$TrueFishIdx_pop,
                  ObsFishAgeComps_pop = sim_env$ObsFishAgeComps_pop,
                  ObsFishIdx_pop_SE = sim_env$ObsFishIdx_pop_SE,
                  ObsFishLenComps_pop = sim_env$ObsFishLenComps_pop,

                  # Population-specific fishery discards
                  TrueDiscard_pop = sim_env$TrueDiscard_pop,
                  ObsDiscard_pop = sim_env$ObsDiscard_pop,
                  ObsFishAgeComps_discard_pop = sim_env$ObsFishAgeComps_discard_pop,
                  ObsFishLenComps_discard_pop = sim_env$ObsFishLenComps_discard_pop,

                  # True fishery compositions
                  CAA = sim_env$CAA,
                  CAL = sim_env$CAL,

                  # True Discards
                  DAA = sim_env$DAA,
                  DAL = sim_env$DAL,

                  # Aggregated survey obs
                  ObsSrvIdx = sim_env$ObsSrvIdx,
                  TrueSrvIdx = sim_env$TrueSrvIdx,
                  ObsSrvIdx_SE = sim_env$ObsSrvIdx_SE,
                  SrvIAA = sim_env$SrvIAA,
                  SrvIAL = sim_env$SrvIAL,
                  srv_sel = sim_env$srv_sel,
                  srv_q = sim_env$srv_q,
                  ObsSrvAgeComps = sim_env$ObsSrvAgeComps,
                  ObsSrvLenComps = sim_env$ObsSrvLenComps,

                  # Population-specific survey obs
                  ObsSrvIdx_pop = sim_env$ObsSrvIdx_pop,
                  TrueSrvIdx_pop = sim_env$TrueSrvIdx_pop,
                  ObsSrvIdx_pop_SE = sim_env$ObsSrvIdx_pop_SE,
                  ObsSrvAgeComps_pop = sim_env$ObsSrvAgeComps_pop,
                  ObsSrvLenComps_pop = sim_env$ObsSrvLenComps_pop,

                  # Tagging
                  conv_tag_release_indicator = as.matrix(sim_env$conv_tag_release_indicator),
                  conv_tag_fish_reporting = sim_env$conv_tag_fish_reporting,
                  conv_tagged_fish = sim_env$conv_tagged_fish,
                  conv_tagged_fish_attr = sim_env$conv_tagged_fish_attr,
                  ln_init_conv_tag_mort = sim_env$ln_init_conv_tag_mort,
                  ln_conv_tag_shed = sim_env$ln_conv_tag_shed,
                  conv_tag_fish_avail = sim_env$conv_tag_fish_avail,
                  use_conv_fish_tagging = sim_env$use_conv_fish_tagging,
                  pred_conv_tag_fish_recap = sim_env$pred_conv_tag_fish_recap,
                  obs_conv_tag_fish_recap = sim_env$obs_conv_tag_fish_recap,

                  # Composition infrastructure
                  SizeAgeTrans = if(!is.null(sim_env$SizeAgeTrans)) sim_env$SizeAgeTrans else NULL,
                  AgeingError = sim_env$AgeingError,
                  ISS_FishAgeComps = sim_env$ISS_FishAgeComps,
                  ISS_FishLenComps = sim_env$ISS_FishLenComps,
                  ISS_SrvAgeComps = sim_env$ISS_SrvAgeComps,
                  ISS_SrvLenComps = sim_env$ISS_SrvLenComps,
                  ISS_FishAgeComps_pop = sim_env$ISS_FishAgeComps_pop,
                  ISS_FishLenComps_pop = sim_env$ISS_FishLenComps_pop,
                  ISS_SrvAgeComps_pop = sim_env$ISS_SrvAgeComps_pop,
                  ISS_SrvLenComps_pop = sim_env$ISS_SrvLenComps_pop,

                  # Discard composition ISS
                  ISS_FishAgeComps_discard = sim_env$ISS_FishAgeComps_discard,
                  ISS_FishLenComps_discard = sim_env$ISS_FishLenComps_discard,
                  ISS_FishAgeComps_discard_pop = sim_env$ISS_FishAgeComps_discard_pop,
                  ISS_FishLenComps_discard_pop = sim_env$ISS_FishLenComps_discard_pop,

                  # Dimensions
                  n_sims = sim_env$n_sims,
                  n_regions = sim_env$n_regions,
                  n_pop = sim_env$n_pop,
                  n_years = sim_env$n_yrs,  # duplicated to ensure backwards compatibility
                  n_yrs = sim_env$n_yrs,    # duplicated to ensure backwards compatibility
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

} # end f


#' Run a simulation self-test of a fitted RTMB estimation model
#'
#' Validates model performance by: (1) generating \code{n_sims} new datasets
#' from the fitted model parameters using \code{\link{Simulate_Pop_Static}},
#' (2) re-fitting the estimation model to each simulated dataset, and (3)
#' storing user-specified report quantities for comparison against true values.
#' Supports sequential or parallel execution via
#' \code{future}/\code{future.apply}. Likelihood weights from the original fit
#' are propagated into the simulation (e.g., ISS scaled by
#' \code{Wt_FishAgeComps}; \code{ObsSrvIdx_SE} divided by
#' \code{sqrt(Wt_SrvIdx)}); all weights are reset to 1 when re-fitting.
#' Failed replicates are silently stored as \code{NA}.
#'
#' @param data Named list of model data from a fitted RTMB object
#'   (\code{$data}).
#' @param parameters Named list of fitted parameter values (\code{$par} or
#'   equivalent).
#' @param mapping Named list of parameter factor maps (\code{$map}).
#' @param random Character vector of random effect names passed to RTMB.
#' @param rep Named list of report values from the fitted model
#'   (\code{obj$rep}).
#' @param sd_rep \code{sdreport} object from the fitted model, used to
#'   extract optimised parameter values in list format via
#'   \code{get_optim_param_list}.
#' @param n_sims Integer. Number of simulation replicates.
#' @param newton_loops Integer. Number of Newton refinement steps applied
#'   during re-fitting. Default \code{3}.
#' @param do_sdrep Logical. Whether to compute \code{sdreport} for each
#'   fitted replicate. Results stored as \code{$sd_rep} in the output list;
#'   failed \code{sdreport} calls stored as \code{NA}. Default \code{FALSE}.
#' @param do_par Logical. Whether to run replicates in parallel via
#'   \code{future::multisession}. Default \code{FALSE}.
#' @param n_cores Integer. Number of parallel workers. If \code{NULL}
#'   (default), \code{parallel::detectCores() - 1} is used.
#' @param output_path Character string. Path to save the simulated dataset
#'   RDS file. Passed to \code{\link{Simulate_Pop_Static}}. Default
#'   \code{NULL}.
#' @param what Character vector. Names of report elements (keys of
#'   \code{rep}) to extract and store from each replicate. An error is raised
#'   if any name is not found in \code{rep}. Default \code{c("SSB", "Rec")}.
#'
#' @return Named list with one element per entry in \code{what}, each an
#'   array with the last dimension indexing simulation replicates (via
#'   \code{simplify2array}). If \code{do_sdrep = TRUE}, an additional element
#'   \code{"sd_rep"} contains a list of \code{sdreport} objects (or \code{NA}
#'   for failed replicates).
#'
#'
#' @export
#' @family Simulation Setup
#'
#' @examples
#' \dontrun{
#' res <- simulation_self_test(
#'   data = obj$data, parameters = obj$par, mapping = obj$map,
#'   random = obj$random, rep = obj$rep, sd_rep = obj$sd_rep,
#'   n_sims = 100, what = c("SSB", "Rec", "Fmort")
#' )
#' str(res$SSB)
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
  if(any(is.na(data$Wt_Catch_pop))) data$Wt_Catch_pop[is.na(data$Wt_Catch_pop)] <- 0
  if(any(is.na(data$Wt_FishAgeComps))) data$Wt_FishAgeComps[is.na(data$Wt_FishAgeComps)] <- 0
  if(any(is.na(data$Wt_FishLenComps))) data$Wt_FishLenComps[is.na(data$Wt_FishLenComps)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_discard))) data$Wt_FishAgeComps_discard[is.na(data$Wt_FishAgeComps_discard)] <- 0
  if(any(is.na(data$Wt_FishLenComps_discard))) data$Wt_FishLenComps_discard[is.na(data$Wt_FishLenComps_discard)] <- 0
  if(any(is.na(data$Wt_FishIdx))) data$Wt_FishIdx[is.na(data$Wt_FishIdx)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_pop))) data$Wt_FishAgeComps_pop[is.na(data$Wt_FishAgeComps_pop)] <- 0
  if(any(is.na(data$Wt_FishLenComps_pop))) data$Wt_FishLenComps_pop[is.na(data$Wt_FishLenComps_pop)] <- 0
  if(any(is.na(data$Wt_FishAgeComps_discard_pop))) data$Wt_FishAgeComps_pop[is.na(data$Wt_FishAgeComps_discard_pop)] <- 0
  if(any(is.na(data$Wt_FishLenComps_discard_pop))) data$Wt_FishLenComps_pop[is.na(data$Wt_FishLenComps_discard_pop)] <- 0
  if(any(is.na(data$Wt_FishIdx_pop))) data$Wt_FishIdx_pop[is.na(data$Wt_FishIdx_pop)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps))) data$Wt_SrvAgeComps[is.na(data$Wt_SrvAgeComps)] <- 0
  if(any(is.na(data$Wt_SrvLenComps))) data$Wt_SrvLenComps[is.na(data$Wt_SrvLenComps)] <- 0
  if(any(is.na(data$Wt_SrvIdx))) data$Wt_SrvIdx[is.na(data$Wt_SrvIdx)] <- 0
  if(any(is.na(data$Wt_SrvAgeComps_pop))) data$Wt_SrvAgeComps_pop[is.na(data$Wt_SrvAgeComps_pop)] <- 0
  if(any(is.na(data$Wt_SrvLenComps_pop))) data$Wt_SrvLenComps_pop[is.na(data$Wt_SrvLenComps_pop)] <- 0
  if(any(is.na(data$Wt_SrvIdx_pop))) data$Wt_SrvIdx_pop[is.na(data$Wt_SrvIdx_pop)] <- 0
  if(any(is.na(data$Wt_Tagging))) data$Wt_Tagging[is.na(data$Wt_Tagging)] <- 0

  # Setup Model Dimensions --------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims, # number of simulations
                            n_yrs = length(data$years), # number of years
                            n_regions = data$n_regions,  # number of regions
                            n_ages = length(data$ages), # number of ages
                            # Use fishery or survey observed ages depending on what is availiable
                            n_obs_ages = if(any(data$UseFishAgeComps == 1)) {
                              dim(data$ObsFishAgeComps)[4]
                            } else if(any(data$UseFishAgeComps_pop == 1)) {
                              dim(data$ObsFishAgeComps_pop)[5]
                            } else if(any(data$UseSrvAgeComps == 1)) {
                              dim(data$ObsSrvAgeComps)[4]
                            } else if(any(data$UseSrvAgeComps_pop == 1)) {
                              dim(data$ObsSrvAgeComps_pop)[5]
                            },
                            n_lens = length(data$lens), # number of lengths
                            n_sexes = data$n_sexes, # number of sexes
                            n_fish_fleets = data$n_fish_fleets, # number of fishery fleets
                            n_srv_fleets = data$n_srv_fleets, # number of survey fleets
                            # Seasonal stuff
                            n_seas = data$n_seas,
                            seasdur = data$seasdur,
                            # Population stuff
                            n_pop = data$n_pop,
                            natal_region = data$natal_region
  )

  # Setup Simulation Containers ---------------------------------------------
  sim_list <- Setup_Sim_Containers(sim_list)

  # Setup Fishing Processes -------------------------------------------------

  # Region-specific sigmaC
  ln_sigmaC <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch)) ln_sigmaC[r,,,f] <- log(exp(optim_parameters_list$ln_sigmaC[r,,,f]) / sqrt(data$Wt_Catch[r,,,f]))
    else ln_sigmaC[r,,,f] <- log(exp(optim_parameters_list$ln_sigmaC[r,,,f]) / sqrt(data$Wt_Catch))
  }

  # Population-specific sigmaC
  ln_sigmaC_pop <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaC container
  # Loop through to populate ln_sigmaC with associated weights
  for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Catch_pop)) ln_sigmaC_pop[p,r,,,f] <- log(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]) / sqrt(data$Wt_Catch_pop[p,r,,,f]))
    else ln_sigmaC_pop[p,r,,,f] <- log(exp(optim_parameters_list$ln_sigmaC_pop[p,r,,,f]) / sqrt(data$Wt_Catch_pop))
  }

  # Region-specific sigmaD
  ln_sigmaD <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaD container
  # Loop through to populate ln_sigmaD with associated weights
  for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Discard)) ln_sigmaD[r,,,f] <- log(exp(optim_parameters_list$ln_sigmaD[r,,,f]) / sqrt(data$Wt_Discard[r,,,f]))
    else ln_sigmaD[r,,,f] <- log(exp(optim_parameters_list$ln_sigmaD[r,,,f]) / sqrt(data$Wt_Discard))
  }

  # Population-specific sigmaD
  ln_sigmaD_pop <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)) # setup sigmaD container
  # Loop through to populate ln_sigmaD with associated weights
  for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) for(f in 1:sim_list$n_fish_fleets) {
    if(!is.vector(data$Wt_Discard_pop)) ln_sigmaD_pop[p,r,,,f] <- log(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]) / sqrt(data$Wt_Discard_pop[p,r,,,f]))
    else ln_sigmaD_pop[p,r,,,f] <- log(exp(optim_parameters_list$ln_sigmaD_pop[p,r,,,f]) / sqrt(data$Wt_Discard_pop))
  }

  # setup fishery simulation processes
  sim_list <- Setup_Sim_Fishing(sim_list = sim_list,
                                ln_sigmaC = ln_sigmaC,
                                ln_sigmaC_pop = ln_sigmaC_pop,
                                ln_sigmaD = ln_sigmaD,
                                ln_sigmaD_pop = ln_sigmaD_pop,
                                catch_units = data$catch_units,
                                discard_units = data$discard_units,
                                Fmort_input = replicate(n = sim_list$n_sims, rep$Fmort[,1:length(data$years),,,drop = FALSE]),
                                dmr_input = replicate(n = sim_list$n_sims, rep$dmr[,1:length(data$years),,,drop = FALSE]),
                                fish_sel_input = replicate(n = sim_list$n_sims, rep$fish_sel[,,1:length(data$years),,,,,drop = FALSE]),
                                ret_sel_input = replicate(n = sim_list$n_sims, rep$ret_sel[,,1:length(data$years),,,,,drop = FALSE]),
                                fish_q_input = replicate(n = sim_list$n_sims, rep$fish_q[,1:length(data$years),,drop = FALSE]),
                                ObsFishIdx_SE = data$ObsFishIdx_SE / sqrt(data$Wt_FishIdx),
                                ObsFishIdx_pop_SE = if(any(data$UseFishIdx_pop == 1)) {
                                  data$ObsFishIdx_pop_SE / sqrt(data$Wt_FishIdx_pop)
                                } else {
                                  array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
                                },
                                fish_idx_type = data$fish_idx_type,
                                init_F_val = rep$init_F,

                                # fishery age composition specifications
                                comp_fishage_like = data$FishAgeComps_LikeType,
                                FishAgeComps_Type = data$FishAgeComps_Type,
                                ISS_FishAgeComps = replicate(sim_list$n_sims, data$ISS_FishAgeComps[,,,,,drop = F] * data$Wt_FishAgeComps),
                                ln_FishAge_theta = optim_parameters_list$ln_FishAge_theta[,,,drop = F],
                                ln_FishAge_theta_agg = optim_parameters_list$ln_FishAge_theta_agg,
                                FishAge_corr_pars_agg = optim_parameters_list$FishAge_corr_pars_agg,
                                FishAge_corr_pars = optim_parameters_list$FishAge_corr_pars[,,,,drop = F],

                                # fishery length composition specifications
                                comp_fishlen_like = data$FishLenComps_LikeType,
                                FishLenComps_Type = data$FishLenComps_Type,
                                ISS_FishLenComps = replicate(sim_list$n_sims, data$ISS_FishLenComps[,,,,,drop = F] * data$Wt_FishLenComps),
                                ln_FishLen_theta = optim_parameters_list$ln_FishLen_theta[,,,drop = F],
                                ln_FishLen_theta_agg = optim_parameters_list$ln_FishLen_theta_agg,
                                FishLen_corr_pars_agg = optim_parameters_list$FishLen_corr_pars_agg,
                                FishLen_corr_pars = optim_parameters_list$FishLen_corr_pars[,,,,drop = F],

                                # population-specific age composition specifications
                                comp_fishage_pop_like = data$pop_FishAgeComps_LikeType,
                                FishAgeComps_pop_Type = data$FishAgeComps_pop_Type,
                                ISS_FishAgeComps_pop = if(any(data$UseFishAgeComps_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishAgeComps_pop[,,,,,,drop = F] * data$Wt_FishAgeComps_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishAge_pop_theta = optim_parameters_list$ln_FishAge_pop_theta[,,,,drop = F],
                                ln_FishAge_pop_theta_agg = optim_parameters_list$ln_FishAge_pop_theta_agg,
                                FishAge_pop_corr_pars_agg = optim_parameters_list$FishAge_pop_corr_pars_agg,
                                FishAge_pop_corr_pars = optim_parameters_list$FishAge_pop_corr_pars[,,,,,drop = F],

                                # population-specific length composition specifications
                                comp_fishlen_pop_like = data$FishLenComps_pop_LikeType,
                                FishLenComps_pop_Type = data$FishLenComps_pop_Type,
                                ISS_FishLenComps_pop = if(any(data$UseFishLenComps_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishLenComps_pop[,,,,,,drop = F] * data$Wt_FishLenComps_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishLen_pop_theta = optim_parameters_list$ln_FishLen_pop_theta[,,,,drop = F],
                                ln_FishLen_pop_theta_agg = optim_parameters_list$ln_FishLen_pop_theta_agg,
                                FishLen_pop_corr_pars_agg = optim_parameters_list$FishLen_pop_corr_pars_agg,
                                FishLen_pop_corr_pars = optim_parameters_list$FishLen_pop_corr_pars[,,,,,drop = F],

                                # discarded fishery age composition specifications
                                comp_fishage_discard_like = data$FishAgeComps_discard_LikeType,
                                FishAgeComps_discard_Type = data$FishAgeComps_discard_Type,
                                ISS_FishAgeComps_discard = replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard[,,,,,drop = F] * data$Wt_FishAgeComps_discard),
                                ln_FishAge_discard_theta = optim_parameters_list$ln_FishAge_discard_theta[,,,drop = F],
                                ln_FishAge_discard_theta_agg = optim_parameters_list$ln_FishAge_discard_theta_agg,
                                FishAge_discard_corr_pars_agg = optim_parameters_list$FishAge_discard_corr_pars_agg,
                                FishAge_discard_corr_pars = optim_parameters_list$FishAge_discard_corr_pars[,,,,drop = F],

                                # discarded fishery length composition specifications
                                comp_fishlen_discard_like = data$FishLenComps_discard_LikeType,
                                FishLenComps_discard_Type = data$FishLenComps_discard_Type,
                                ISS_FishLenComps_discard = replicate(sim_list$n_sims, data$ISS_FishLenComps_discard[,,,,,drop = F] * data$Wt_FishLenComps_discard),
                                ln_FishLen_discard_theta = optim_parameters_list$ln_FishLen_discard_theta[,,,drop = F],
                                ln_FishLen_discard_theta_agg = optim_parameters_list$ln_FishLen_discard_theta_agg,
                                FishLen_discard_corr_pars_agg = optim_parameters_list$FishLen_discard_corr_pars_agg,
                                FishLen_discard_corr_pars = optim_parameters_list$FishLen_discard_corr_pars[,,,,drop = F],

                                # discarded population-specific age composition specifications
                                comp_fishage_discard_pop_like = data$FishAgeComps_discard_pop_LikeType,
                                FishAgeComps_discard_pop_Type = data$FishAgeComps_discard_pop_Type,
                                ISS_FishAgeComps_discard_pop = if(any(data$UseFishAgeComps_discard_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishAgeComps_discard_pop[,,,,,,drop = F] * data$Wt_FishAgeComps_discard_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishAge_discard_pop_theta = optim_parameters_list$ln_FishAge_discard_pop_theta[,,,,drop = F],
                                ln_FishAge_discard_pop_theta_agg = optim_parameters_list$ln_FishAge_discard_pop_theta_agg,
                                FishAge_discard_pop_corr_pars_agg = optim_parameters_list$FishAge_discard_pop_corr_pars_agg,
                                FishAge_discard_pop_corr_pars = optim_parameters_list$FishAge_discard_pop_corr_pars[,,,,,drop = F],

                                # discarded population-specific length composition specifications
                                comp_fishlen_discard_pop_like = data$FishLenComps_discard_pop_LikeType,
                                FishLenComps_discard_pop_Type = data$FishLenComps_discard_pop_Type,
                                ISS_FishLenComps_discard_pop = if(any(data$UseFishLenComps_discard_pop == 1)) {
                                  replicate(sim_list$n_sims, data$ISS_FishLenComps_discard_pop[,,,,,,drop = F] * data$Wt_FishLenComps_discard_pop)
                                } else {
                                  array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims))
                                },
                                ln_FishLen_discard_pop_theta = optim_parameters_list$ln_FishLen_discard_pop_theta[,,,,drop = F],
                                ln_FishLen_discard_pop_theta_agg = optim_parameters_list$ln_FishLen_discard_pop_theta_agg,
                                FishLen_discard_pop_corr_pars_agg = optim_parameters_list$FishLen_discard_pop_corr_pars_agg,
                                FishLen_discard_pop_corr_pars = optim_parameters_list$FishLen_discard_pop_corr_pars[,,,,,drop = F]
  )

  # Setup Survey Processes --------------------------------------------------
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(n = sim_list$n_sims, rep$srv_sel[,,1:length(data$years),,,,,drop = FALSE]),
    srv_q_input = replicate(n = sim_list$n_sims, rep$srv_q[,1:length(data$years),,drop = FALSE]),
    ObsSrvIdx_SE = data$ObsSrvIdx_SE / sqrt(data$Wt_SrvIdx),
    ObsSrvIdx_pop_SE = if(any(data$UseSrvIdx_pop == 1)) {
      data$ObsSrvIdx_pop_SE / sqrt(data$Wt_SrvIdx_pop)
    } else {
      array(0.2, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))
    },
    srv_idx_type = data$srv_idx_type,
    t_srv = data$t_srv,

    # survey age composition specifications
    comp_srvage_like = data$SrvAgeComps_LikeType,
    SrvAgeComps_Type = data$SrvAgeComps_Type,
    ISS_SrvAgeComps = replicate(sim_list$n_sims, data$ISS_SrvAgeComps[,,,,,drop = F] * data$Wt_SrvAgeComps),
    ln_SrvAge_theta = optim_parameters_list$ln_SrvAge_theta[,,,drop = F],
    ln_SrvAge_theta_agg = optim_parameters_list$ln_SrvAge_theta_agg,
    SrvAge_corr_pars_agg = optim_parameters_list$SrvAge_corr_pars_agg,
    SrvAge_corr_pars = optim_parameters_list$SrvAge_corr_pars[,,,,drop = F],

    # survey length composition specifications
    comp_srvlen_like = data$SrvLenComps_LikeType,
    SrvLenComps_Type = data$SrvLenComps_Type,
    ISS_SrvLenComps = replicate(sim_list$n_sims, data$ISS_SrvLenComps[,,,,,drop = F] * data$Wt_SrvLenComps),
    ln_SrvLen_theta = optim_parameters_list$ln_SrvLen_theta[,,,drop = F],
    ln_SrvLen_theta_agg = optim_parameters_list$ln_SrvLen_theta_agg,
    SrvLen_corr_pars_agg = optim_parameters_list$SrvLen_corr_pars_agg,
    SrvLen_corr_pars = optim_parameters_list$SrvLen_corr_pars[,,,,drop = F],

    # population-specific age composition specifications
    comp_srvage_pop_like = data$SrvAgeComps_pop_LikeType,
    SrvAgeComps_pop_Type = data$SrvAgeComps_pop_Type,
    ISS_SrvAgeComps_pop = if(any(data$UseSrvAgeComps_pop == 1)) {
      replicate(sim_list$n_sims, data$ISS_SrvAgeComps_pop[,,,,,,drop = F] * data$Wt_SrvAgeComps_pop)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    },
    ln_SrvAge_pop_theta = optim_parameters_list$ln_SrvAge_pop_theta[,,,,drop = F],
    ln_SrvAge_pop_theta_agg = optim_parameters_list$ln_SrvAge_pop_theta_agg,
    SrvAge_pop_corr_pars_agg = optim_parameters_list$SrvAge_pop_corr_pars_agg,
    SrvAge_pop_corr_pars = optim_parameters_list$SrvAge_pop_corr_pars[,,,,,drop = F],

    # population-specific length composition specifications
    comp_srvlen_pop_like = data$SrvLenComps_pop_LikeType,
    SrvLenComps_pop_Type = data$SrvLenComps_pop_Type,
    ISS_SrvLenComps_pop = if(any(data$UseSrvLenComps_pop == 1)) {
      replicate(sim_list$n_sims, data$ISS_SrvLenComps_pop[,,,,,,drop = F] * data$Wt_SrvLenComps_pop)
    } else {
      array(100, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims))
    },
    ln_SrvLen_pop_theta = optim_parameters_list$ln_SrvLen_pop_theta[,,,,drop = F],
    ln_SrvLen_pop_theta_agg = optim_parameters_list$ln_SrvLen_pop_theta_agg,
    SrvLen_pop_corr_pars_agg = optim_parameters_list$SrvLen_pop_corr_pars_agg,
    SrvLen_pop_corr_pars = optim_parameters_list$SrvLen_pop_corr_pars[,,,,,drop = F]
  )

  # Setup Biological Dynamics -----------------------------------------------
  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list, # simualtion list
    natmort_input = replicate(n = sim_list$n_sims, rep$natmort[,,1:length(data$years),,,drop = FALSE]), # natuyral mortality
    WAA_input = replicate(n = sim_list$n_sims, data$WAA[,,1:length(data$years),,,,drop = FALSE]), # weight at age
    WAA_fish_input = replicate(n = sim_list$n_sims, data$WAA_fish[,,1:length(data$years),,,,,drop = FALSE]), # fishery weight at age
    WAA_srv_input = replicate(n = sim_list$n_sims, data$WAA_srv[,,1:length(data$years),,,,,drop = FALSE]), # survey weight at age
    MatAA_input = replicate(n = sim_list$n_sims, data$MatAA[,,1:length(data$years),,,,drop = FALSE]), # maturity at age
    AgeingError_input = replicate(n = sim_list$n_sims, data$AgeingError[1:length(data$years),,,drop = FALSE]), # ageing error
    SizeAgeTrans_input = if(data$fit_lengths == 0) NULL else replicate(n = sim_list$n_sims, data$SizeAgeTrans[,,1:length(data$years),,,,,drop = FALSE]) # size age transition matrix
  )

  # Movement
  sim_list$Movement <- replicate(n = sim_list$n_sims, rep$Movement[,,,1:length(data$years),,,,drop = FALSE])
  sim_list$sgl_seas_spawning_movement <- replicate(n = sim_list$n_sims, rep$sgl_seas_spawning_movement[,,,1:length(data$years),,,drop = FALSE])

  # Setup Recruitment Processes ---------------------------------------------
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    spawn_seas = data$spawn_seas, # spawning season
    do_recruits_move = data$do_recruits_move, # whether recruits move
    t_spawn = data$t_spawn, # spawn timing
    init_age_strc = data$init_age_strc, # initilaizing age structure
    h_input = replicate(n = sim_list$n_sims, array(rep$h_trans, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs))), # steepness
    R0_input = {
      tmp = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
      for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) tmp[p,r,,] = rep$R0[p] * rep$rec_region_prop[p,r]
      tmp
    },
    rinit_input = {
      tmp = array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_sims))
      for(p in 1:sim_list$n_pop) for(r in 1:sim_list$n_regions) tmp[p,r,] = rep$rinit[p] * rep$rec_region_prop[p,r]
      tmp
    },
    use_rinit = data$use_rinit,
    sexratio_input = replicate(n = sim_list$n_sims, expr = rep$sexratio[,,1:length(data$years),,drop = FALSE]), # sex ratio
    ln_sigmaR = optim_parameters_list$ln_sigmaR / sqrt(data$Wt_Rec), # ln_sigmaR
    Rec_input = replicate(n = sim_list$n_sims, expr = rep$Rec[,,1:length(data$years),drop = FALSE]), # recruitment time series
    ln_InitDevs_input = replicate(sim_list$n_sims, optim_parameters_list$ln_InitDevs),  # init devs
    stray_rate_input = replicate(sim_list$n_sims, data$stray_rate[,1:length(data$years), drop = FALSE]),
    rec_seas_prop_input = array(
      rep(rep$rec_seas_prop, times = sim_list$n_sims),
      dim = c(data$n_pop, data$n_seas, sim_list$n_sims)), # seasonal recruitment apportionment

    # Not needed; already specified in Rec_input and ln_InitDevs_input
    recruitment_opt = data$rec_model,
    rec_dd = data$rec_dd,
    init_dd = data$rec_dd,
    rec_lag = data$rec_lag
  )

  # Setup Tagging -----------------------------------------------------------
  if(!is.na(sum(data$conv_tagged_fish))) n_tags_rel_input <- apply(data$conv_tagged_fish, 1, sum) else n_tags_rel_input <- NA
  if(exists("conv_tag_release_indicator", data)) conv_tag_release_indicator <- data$conv_tag_release_indicator  else conv_tag_release_indicator <- NA
  conv_tag_fish_reporting_input <- if(!is.null(rep$conv_tag_fish_reporting)) replicate(n = sim_list$n_sims, rep$conv_tag_fish_reporting) else NULL

  sim_list <- Setup_Sim_Tagging(
    sim_list = sim_list, # simulation list
    conv_tag_max_liberty = data$conv_tag_max_liberty, # maximum tag liberty
    conv_tag_t_tagging = data$conv_tag_t_tagging, # time of tagging
    n_tags_rel_input = n_tags_rel_input * data$Wt_Tagging,  # number of tags to release per event
    conv_tag_release_indicator = conv_tag_release_indicator,  # tag release indicator
    ln_init_conv_tag_mort = optim_parameters_list$ln_init_conv_tag_mort,  # inital tagging mortality
    ln_conv_tag_shed = optim_parameters_list$ln_conv_tag_shed, # chronic tag shedding
    conv_tag_fish_reporting_input = conv_tag_fish_reporting_input, # tag reporting rates
    use_conv_fish_tagging = data$use_conv_fish_tagging, # whether or not tagging is used / simulated
    conv_fish_tag_like = data$conv_fish_tag_like, # tag likelihood
    ln_conv_fish_tag_theta = parameters$ln_conv_fish_tag_theta, # tag overdispersion
    conv_tag_release_platform = data$conv_tag_release_platform,  # tag release platform
    conv_fish_tag_attr = data$conv_fish_tag_attr # tag attributes
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
        tmp_data <- data
        tmp_pars <- parameters
        tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx))
        tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx))
        tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch))
        tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps))
        tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps))
        if(tmp_data$fit_lengths != 0) {
          tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps))
          tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps))
        }

        # population-specific observations
        if(any(tmp_data$UseFishIdx_pop == 1)) {
          tmp_data$ObsFishIdx_pop <- array(sim_obj$ObsFishIdx_pop[,,,,,i], dim = dim(tmp_data$ObsFishIdx_pop))
        }
        if(any(tmp_data$UseSrvIdx_pop == 1)) {
          tmp_data$ObsSrvIdx_pop <- array(sim_obj$ObsSrvIdx_pop[,,,,,i], dim = dim(tmp_data$ObsSrvIdx_pop))
        }
        if(any(tmp_data$UseFishAgeComps_pop == 1)) {
          tmp_data$ObsFishAgeComps_pop <- array(sim_obj$ObsFishAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_pop))
        }
        if(any(tmp_data$UseSrvAgeComps_pop == 1)) {
          tmp_data$ObsSrvAgeComps_pop <- array(sim_obj$ObsSrvAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps_pop))
        }
        if(tmp_data$fit_lengths != 0) {
          if(any(tmp_data$UseFishLenComps_pop == 1)) {
            tmp_data$ObsFishLenComps_pop <- array(sim_obj$ObsFishLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_pop))
          }
          if(any(tmp_data$UseSrvLenComps_pop == 1)) {
            tmp_data$ObsSrvLenComps_pop <- array(sim_obj$ObsSrvLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps_pop))
          }
        }
        if(any(tmp_data$UseCatch_pop == 1)) {
          tmp_data$ObsCatch_pop <- array(sim_obj$ObsCatch_pop[,,,,,i], dim = dim(tmp_data$ObsCatch_pop))
        }

        # set up discarding stuff
        if(any(tmp_data$UseDiscard == 1)) tmp_data$ObsDiscard <- array(sim_obj$ObsDiscard[,,,,i], dim = dim(tmp_data$ObsDiscard))
        if(any(tmp_data$UseDiscard_pop == 1)) tmp_data$ObsDiscard_pop <- array(sim_obj$ObsDiscard_pop[,,,,i], dim = dim(tmp_data$ObsDiscard_pop))
        if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ObsFishAgeComps_discard <- array(sim_obj$ObsFishAgeComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard))
        if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ObsFishLenComps_discard <- array(sim_obj$ObsFishLenComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard))
        if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ObsFishAgeComps_discard_pop <- array(sim_obj$ObsFishAgeComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard_pop))
        if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ObsFishLenComps_discard_pop <- array(sim_obj$ObsFishLenComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard_pop))

        # setup tagging data stuff if tagging is done
        if(any(tmp_data$use_conv_fish_tagging == 1)) {
          tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,,i], dim = dim(tmp_data$conv_tagged_fish))
          tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap))
          tmp_data$conv_tag_release_indicator <- sim_obj$conv_tag_release_indicator
        }

        # reset weights
        tmp_data$Wt_Rec <- 1
        tmp_data$Wt_D <- 1
        tmp_data$Wt_Tagging <- 1
        tmp_data$Wt_Catch[] <- 1
        tmp_data$Wt_Discard[] <- 1
        tmp_data$Wt_FishAgeComps[] <- 1
        tmp_data$Wt_FishAgeComps_discard[] <- 1
        tmp_data$Wt_FishIdx <- 1
        tmp_data$Wt_FishLenComps[] <- 1
        tmp_data$Wt_FishLenComps_discard[] <- 1
        tmp_data$Wt_SrvAgeComps[] <- 1
        tmp_data$Wt_SrvIdx <- 1
        tmp_data$Wt_SrvLenComps[] <- 1
        tmp_data$Wt_Catch_pop[] <- 1
        tmp_data$Wt_Discard_pop[] <- 1
        tmp_data$Wt_FishIdx_pop[] <- 1
        tmp_data$Wt_SrvIdx_pop[] <- 1
        tmp_data$Wt_FishAgeComps_pop[] <- 1
        tmp_data$Wt_FishAgeComps_discard_pop[] <- 1
        tmp_data$Wt_SrvAgeComps_pop[] <- 1
        tmp_data$Wt_FishLenComps_pop[] <- 1
        tmp_data$Wt_FishLenComps_discard_pop[] <- 1
        tmp_data$Wt_SrvLenComps_pop[] <- 1

        # input simulated uncertainty
        tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
        tmp_pars$ln_sigmaC_pop[] <- sim_list$ln_sigmaC_pop
        tmp_pars$ln_sigmaD[] <- sim_list$ln_sigmaD
        tmp_pars$ln_sigmaD_pop[] <- sim_list$ln_sigmaD_pop
        tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
        tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
        tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,,i]
        tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,,i]
        if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ISS_FishAgeComps_discard[] <- sim_list$ISS_FishAgeComps_discard[,,,,,i]
        if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ISS_FishLenComps_discard[] <- sim_list$ISS_FishLenComps_discard[,,,,,i]
        tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,,i]
        tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,,i]
        if(any(tmp_data$UseFishIdx_pop == 1)) tmp_data$ObsFishIdx_pop_SE[] <- sim_list$ObsFishIdx_pop_SE
        if(any(tmp_data$UseSrvIdx_pop == 1)) tmp_data$ObsSrvIdx_pop_SE[] <- sim_list$ObsSrvIdx_pop_SE
        if(any(tmp_data$UseFishAgeComps_pop == 1)) tmp_data$ISS_FishAgeComps_pop[] <- sim_list$ISS_FishAgeComps_pop[,,,,,,i]
        if(any(tmp_data$UseFishLenComps_pop == 1)) tmp_data$ISS_FishLenComps_pop[] <- sim_list$ISS_FishLenComps_pop[,,,,,,i]
        if(any(tmp_data$UseSrvAgeComps_pop == 1)) tmp_data$ISS_SrvAgeComps_pop[] <- sim_list$ISS_SrvAgeComps_pop[,,,,,,i]
        if(any(tmp_data$UseSrvLenComps_pop == 1)) tmp_data$ISS_SrvLenComps_pop[] <- sim_list$ISS_SrvLenComps_pop[,,,,,,i]
        if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ISS_FishAgeComps_discard_pop[] <- sim_list$ISS_FishAgeComps_discard_pop[,,,,,i]
        if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ISS_FishLenComps_discard_pop[] <- sim_list$ISS_FishLenComps_discard_pop[,,,,,i]

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
          tmp_data <- data
          tmp_pars <- parameters
          tmp_data$ObsFishIdx <- array(sim_obj$ObsFishIdx[,,,,i], dim = dim(tmp_data$ObsFishIdx))
          tmp_data$ObsSrvIdx <- array(sim_obj$ObsSrvIdx[,,,,i], dim = dim(tmp_data$ObsSrvIdx))
          tmp_data$ObsCatch <- array(sim_obj$ObsCatch[,,,,i], dim = dim(tmp_data$ObsCatch))
          tmp_data$ObsFishAgeComps <- array(sim_obj$ObsFishAgeComps[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps))
          tmp_data$ObsSrvAgeComps  <- array(sim_obj$ObsSrvAgeComps[,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps))
          if(tmp_data$fit_lengths != 0) {
            tmp_data$ObsFishLenComps <- array(sim_obj$ObsFishLenComps[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps))
            tmp_data$ObsSrvLenComps  <- array(sim_obj$ObsSrvLenComps[,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps))
          }

          # population-specific observations
          if(any(tmp_data$UseFishIdx_pop == 1)) {
            tmp_data$ObsFishIdx_pop <- array(sim_obj$ObsFishIdx_pop[,,,,,i], dim = dim(tmp_data$ObsFishIdx_pop))
          }
          if(any(tmp_data$UseSrvIdx_pop == 1)) {
            tmp_data$ObsSrvIdx_pop <- array(sim_obj$ObsSrvIdx_pop[,,,,,i], dim = dim(tmp_data$ObsSrvIdx_pop))
          }
          if(any(tmp_data$UseFishAgeComps_pop == 1)) {
            tmp_data$ObsFishAgeComps_pop <- array(sim_obj$ObsFishAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_pop))
          }
          if(any(tmp_data$UseSrvAgeComps_pop == 1)) {
            tmp_data$ObsSrvAgeComps_pop <- array(sim_obj$ObsSrvAgeComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvAgeComps_pop))
          }
          if(tmp_data$fit_lengths != 0) {
            if(any(tmp_data$UseFishLenComps_pop == 1)) {
              tmp_data$ObsFishLenComps_pop <- array(sim_obj$ObsFishLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_pop))
            }
            if(any(tmp_data$UseSrvLenComps_pop == 1)) {
              tmp_data$ObsSrvLenComps_pop <- array(sim_obj$ObsSrvLenComps_pop[,,,,,,,i], dim = dim(tmp_data$ObsSrvLenComps_pop))
            }
          }
          if(any(tmp_data$UseCatch_pop == 1)) {
            tmp_data$ObsCatch_pop <- array(sim_obj$ObsCatch_pop[,,,,,i], dim = dim(tmp_data$ObsCatch_pop))
          }

          # set up discarding stuff
          if(any(tmp_data$UseDiscard == 1)) tmp_data$ObsDiscard <- array(sim_obj$ObsDiscard[,,,,i], dim = dim(tmp_data$ObsDiscard))
          if(any(tmp_data$UseDiscard_pop == 1)) tmp_data$ObsDiscard_pop <- array(sim_obj$ObsDiscard_pop[,,,,i], dim = dim(tmp_data$ObsDiscard_pop))
          if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ObsFishAgeComps_discard <- array(sim_obj$ObsFishAgeComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard))
          if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ObsFishLenComps_discard <- array(sim_obj$ObsFishLenComps_discard[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard))
          if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ObsFishAgeComps_discard_pop <- array(sim_obj$ObsFishAgeComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishAgeComps_discard_pop))
          if(tmp_data$fit_lengths != 0) if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ObsFishLenComps_discard_pop <- array(sim_obj$ObsFishLenComps_discard_pop[,,,,,,i], dim = dim(tmp_data$ObsFishLenComps_discard_pop))

          # setup tagging data stuff if tagging is done
          if(any(tmp_data$use_conv_fish_tagging == 1)) {
            tmp_data$conv_tagged_fish <- array(sim_obj$conv_tagged_fish[,,,,i], dim = dim(tmp_data$conv_tagged_fish))
            tmp_data$obs_conv_tag_fish_recap <- array(sim_obj$obs_conv_tag_fish_recap[,,,,,,,,i], dim = dim(tmp_data$obs_conv_tag_fish_recap))
            tmp_data$conv_tag_release_indicator <- sim_obj$conv_tag_release_indicator
          }

          # reset weights
          tmp_data$Wt_Rec <- 1
          tmp_data$Wt_D <- 1
          tmp_data$Wt_Tagging <- 1
          tmp_data$Wt_Catch[] <- 1
          tmp_data$Wt_Discard[] <- 1
          tmp_data$Wt_FishAgeComps[] <- 1
          tmp_data$Wt_FishAgeComps_discard[] <- 1
          tmp_data$Wt_FishIdx <- 1
          tmp_data$Wt_FishLenComps[] <- 1
          tmp_data$Wt_FishLenComps_discard[] <- 1
          tmp_data$Wt_SrvAgeComps[] <- 1
          tmp_data$Wt_SrvIdx <- 1
          tmp_data$Wt_SrvLenComps[] <- 1
          tmp_data$Wt_Catch_pop[] <- 1
          tmp_data$Wt_Discard_pop[] <- 1
          tmp_data$Wt_FishIdx_pop[] <- 1
          tmp_data$Wt_SrvIdx_pop[] <- 1
          tmp_data$Wt_FishAgeComps_pop[] <- 1
          tmp_data$Wt_FishAgeComps_discard_pop[] <- 1
          tmp_data$Wt_SrvAgeComps_pop[] <- 1
          tmp_data$Wt_FishLenComps_pop[] <- 1
          tmp_data$Wt_FishLenComps_discard_pop[] <- 1
          tmp_data$Wt_SrvLenComps_pop[] <- 1

          # input simulated uncertainty
          tmp_pars$ln_sigmaC[] <- sim_list$ln_sigmaC
          tmp_pars$ln_sigmaC_pop[] <- sim_list$ln_sigmaC_pop
          tmp_pars$ln_sigmaD[] <- sim_list$ln_sigmaD
          tmp_pars$ln_sigmaD_pop[] <- sim_list$ln_sigmaD_pop
          tmp_data$ObsFishIdx_SE[] <- sim_list$ObsFishIdx_SE
          tmp_data$ObsSrvIdx_SE[] <- sim_list$ObsSrvIdx_SE
          tmp_data$ISS_FishAgeComps[] <- sim_list$ISS_FishAgeComps[,,,,,i]
          tmp_data$ISS_FishLenComps[] <- sim_list$ISS_FishLenComps[,,,,,i]
          if(any(tmp_data$UseFishAgeComps_discard == 1)) tmp_data$ISS_FishAgeComps_discard[] <- sim_list$ISS_FishAgeComps_discard[,,,,,i]
          if(any(tmp_data$UseFishLenComps_discard == 1)) tmp_data$ISS_FishLenComps_discard[] <- sim_list$ISS_FishLenComps_discard[,,,,,i]
          tmp_data$ISS_SrvAgeComps[] <- sim_list$ISS_SrvAgeComps[,,,,,i]
          tmp_data$ISS_SrvLenComps[] <- sim_list$ISS_SrvLenComps[,,,,,i]
          if(any(tmp_data$UseFishIdx_pop == 1)) tmp_data$ObsFishIdx_pop_SE[] <- sim_list$ObsFishIdx_pop_SE
          if(any(tmp_data$UseSrvIdx_pop == 1)) tmp_data$ObsSrvIdx_pop_SE[] <- sim_list$ObsSrvIdx_pop_SE
          if(any(tmp_data$UseFishAgeComps_pop == 1)) tmp_data$ISS_FishAgeComps_pop[] <- sim_list$ISS_FishAgeComps_pop[,,,,,,i]
          if(any(tmp_data$UseFishLenComps_pop == 1)) tmp_data$ISS_FishLenComps_pop[] <- sim_list$ISS_FishLenComps_pop[,,,,,,i]
          if(any(tmp_data$UseSrvAgeComps_pop == 1)) tmp_data$ISS_SrvAgeComps_pop[] <- sim_list$ISS_SrvAgeComps_pop[,,,,,,i]
          if(any(tmp_data$UseSrvLenComps_pop == 1)) tmp_data$ISS_SrvLenComps_pop[] <- sim_list$ISS_SrvLenComps_pop[,,,,,,i]
          if(any(tmp_data$UseFishAgeComps_discard_pop == 1)) tmp_data$ISS_FishAgeComps_discard_pop[] <- sim_list$ISS_FishAgeComps_discard_pop[,,,,,i]
          if(any(tmp_data$UseFishLenComps_discard_pop == 1)) tmp_data$ISS_FishLenComps_discard_pop[] <- sim_list$ISS_FishLenComps_discard_pop[,,,,,i]


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

#' Extract simulation outputs into SPoRC estimation model format
#'
#' Subsets and reshapes biological, tagging, fishery, and survey arrays from a
#' simulation environment or output list to cover years \code{1:y} and
#' simulation replicate \code{sim}, producing a named list ready for direct
#' use in \code{\link{Setup_Mod_Biologicals}}, \code{\link{Setup_Mod_Catch_and_F}},
#' \code{\link{Setup_Mod_SrvIdx_and_Comps}}, and \code{\link{Setup_Mod_Tagging}}.
#' Binary \code{Use*} indicator arrays are derived automatically from the
#' extracted observation arrays (non-NA, positive values set to 1).
#'
#' Population-specific arrays (\code{ObsCatch_pop}, \code{ObsFishIdx_pop},
#' \code{ObsFishAgeComps_pop}, \code{ObsFishLenComps_pop}, \code{ObsSrvIdx_pop},
#' \code{ObsSrvAgeComps_pop}, \code{ObsSrvLenComps_pop}) are extracted when
#' present in \code{sim_env}; corresponding \code{Use*_PopSpec} flags are
#' derived automatically.
#'
#' Input sample sizes for age and length compositions are extracted for both
#' aggregate and population-specific data streams, including retained
#' (\code{ISS_FishAgeComps}, \code{ISS_FishLenComps}, \code{ISS_FishAgeComps_pop},
#' \code{ISS_FishLenComps_pop}, \code{ISS_SrvAgeComps}, \code{ISS_SrvLenComps},
#' \code{ISS_SrvAgeComps_pop}, \code{ISS_SrvLenComps_pop}) and discard
#' (\code{ISS_FishAgeComps_discard}, \code{ISS_FishLenComps_discard},
#' \code{ISS_FishAgeComps_discard_pop}, \code{ISS_FishLenComps_discard_pop})
#' data streams.
#'
#' Length composition outputs (\code{ObsFishLenComps}, \code{ObsSrvLenComps},
#' and their population-specific and discard variants) and \code{SizeAgeTrans}
#' are \code{NULL} when no size-age transition matrix is present in
#' \code{sim_env}. Tagging outputs are \code{NULL} when
#' \code{use_conv_fish_tagging = 0}; otherwise, only cohorts with release
#' years in \code{1:y} are retained.
#'
#' @param sim_env Simulation environment or list (e.g., output from
#'   \code{\link{Simulate_Pop_Static}} or a \code{\link{Setup_sim_env}}
#'   environment) containing all operating model arrays.
#' @param y Integer. Last year to include; years \code{1:y} are retained.
#' @param sim Integer. Simulation replicate index to extract.
#'
#' @return Named list with the following elements (all arrays have \code{y}
#'   in the year dimension unless noted):
#'   \code{WAA} \code{[n_pop x n_regions x y x n_seas x n_ages x n_sexes]},
#'   \code{WAA_fish} \code{[... x n_fish_fleets]},
#'   \code{WAA_srv} \code{[... x n_srv_fleets]},
#'   \code{MatAA}, \code{SizeAgeTrans} (or \code{NULL}),
#'   \code{AgeingError} \code{[y x n_obs_ages x n_ages]},
#'   \code{use_conv_fish_tagging},
#'   \code{conv_tag_release_indicator}, \code{obs_conv_tag_fish_recap},
#'   \code{conv_tagged_fish}, \code{conv_tagged_fish_attr},
#'   \code{n_tag_cohorts} (all \code{NULL} when tagging inactive),
#'   \code{ObsCatch}, \code{ln_sigmaC}, \code{UseCatch},
#'   \code{ObsCatch_pop}, \code{ln_sigmaC_pop}, \code{UseCatch_pop},
#'   \code{ObsDiscard}, \code{ln_sigmaD}, \code{UseDiscard},
#'   \code{ObsDiscard_pop}, \code{ln_sigmaD_pop}, \code{UseDiscard_pop},
#'   \code{ObsFishIdx}, \code{ObsFishIdx_SE}, \code{UseFishIdx},
#'   \code{ObsFishIdx_pop}, \code{ObsFishIdx_pop_SE}, \code{UseFishIdx_pop},
#'   \code{ObsFishAgeComps}, \code{ISS_FishAgeComps}, \code{UseFishAgeComps},
#'   \code{ObsFishAgeComps_pop}, \code{ISS_FishAgeComps_pop},
#'   \code{UseFishAgeComps_pop},
#'   \code{ObsFishLenComps} (or \code{NULL}), \code{ISS_FishLenComps} (or \code{NULL}),
#'   \code{UseFishLenComps},
#'   \code{ObsFishLenComps_pop} (or \code{NULL}), \code{ISS_FishLenComps_pop} (or \code{NULL}),
#'   \code{UseFishLenComps_pop},
#'   \code{ObsFishAgeComps_discard}, \code{ISS_FishAgeComps_discard}, \code{UseFishAgeComps_discard},
#'   \code{ObsFishAgeComps_discard_pop}, \code{ISS_FishAgeComps_discard_pop},
#'   \code{UseFishAgeComps_discard_pop},
#'   \code{ObsFishLenComps_discard} (or \code{NULL}), \code{ISS_FishLenComps_discard} (or \code{NULL}),
#'   \code{UseFishLenComps_discard},
#'   \code{ObsFishLenComps_discard_pop} (or \code{NULL}), \code{ISS_FishLenComps_discard_pop} (or \code{NULL}),
#'   \code{UseFishLenComps_discard_pop},
#'   \code{ObsSrvIdx}, \code{ObsSrvIdx_SE}, \code{UseSrvIdx},
#'   \code{ObsSrvIdx_pop}, \code{ObsSrvIdx_pop_SE}, \code{UseSrvIdx_pop},
#'   \code{ObsSrvAgeComps}, \code{ISS_SrvAgeComps}, \code{UseSrvAgeComps},
#'   \code{ObsSrvAgeComps_pop}, \code{ISS_SrvAgeComps_pop},
#'   \code{UseSrvAgeComps_pop},
#'   \code{ObsSrvLenComps} (or \code{NULL}), \code{ISS_SrvLenComps} (or \code{NULL}),
#'   \code{UseSrvLenComps},
#'   \code{ObsSrvLenComps_pop} (or \code{NULL}), \code{ISS_SrvLenComps_pop} (or \code{NULL}),
#'   \code{UseSrvLenComps_pop}.
#'
#' @seealso \code{\link{Setup_Mod_Biologicals}}, \code{\link{Setup_Mod_Catch_and_F}},
#'   \code{\link{Setup_Mod_SrvIdx_and_Comps}}, \code{\link{Setup_Mod_Tagging}},
#'   \code{\link{Simulate_Pop_Static}}, \code{\link{Setup_sim_env}}
#'
#' @export simulation_data_to_SPoRC
#' @family Simulation Utilities
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
  AgeingError <- array(sim_env$AgeingError[1:y,,,sim, drop = FALSE],
                       dim = c(length(1:y), dim(sim_env$AgeingError)[2], dim(sim_env$AgeingError)[3]))

  # Tagging
  if(sim_env$use_conv_fish_tagging == 1) {
    keep_tag_cohorts <- which(sim_env$conv_tag_release_indicator[,2] %in% 1:y)
    conv_tag_release_indicator <- sim_env$conv_tag_release_indicator[keep_tag_cohorts,,drop = FALSE]
    obs_conv_tag_fish_recap <- array(sim_env$obs_conv_tag_fish_recap[,,keep_tag_cohorts,,,,,,sim], dim = dim(sim_env$obs_conv_tag_fish_recap)[-length(dim(sim_env$obs_conv_tag_fish_recap))])
    conv_tagged_fish <- array(sim_env$conv_tagged_fish[keep_tag_cohorts,,,,sim], dim = c(dim(sim_env$conv_tagged_fish)[-length(dim(sim_env$conv_tagged_fish))]))
    conv_tagged_fish_attr <- array(sim_env$conv_tagged_fish_attr[keep_tag_cohorts,,,,sim], dim = c(dim(sim_env$conv_tagged_fish_attr)[-length(dim(sim_env$conv_tagged_fish_attr))]))
    n_tag_cohorts <- nrow(conv_tag_release_indicator)
  } else {
    conv_tag_release_indicator = obs_conv_tag_fish_recap = conv_tagged_fish = conv_tagged_fish_attr = n_tag_cohorts = NULL
  }

  # Fishery Landed Catches
  ObsCatch <- array(sim_env$ObsCatch[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaC <- array(sim_env$ln_sigmaC[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseCatch <- array(0, dim = dim(ObsCatch))
  UseCatch[!is.na(ObsCatch) & ObsCatch > 0] <- 1

  # Population-specific catches
  ObsCatch_pop <- array(sim_env$ObsCatch_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseCatch_pop <- array(0, dim = dim(ObsCatch_pop))
  UseCatch_pop[!is.na(ObsCatch_pop) & ObsCatch_pop > 0] <- 1
  ln_sigmaC_pop <- array(sim_env$ln_sigmaC_pop[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Discards
  ObsDiscard <- array(sim_env$ObsDiscard[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaD <- array(sim_env$ln_sigmaD[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseDiscard <- array(0, dim = dim(ObsDiscard))
  UseDiscard[!is.na(ObsDiscard) & ObsDiscard > 0] <- 1

  # Population-specific discards
  ObsDiscard_pop <- array(sim_env$ObsDiscard_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ln_sigmaD_pop <- array(sim_env$ln_sigmaD_pop[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseDiscard_pop <- array(0, dim = dim(ObsDiscard_pop))
  UseDiscard_pop[!is.na(ObsDiscard_pop) & ObsDiscard_pop > 0] <- 1

  # Fishery Indices
  ObsFishIdx <- array(sim_env$ObsFishIdx[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  ObsFishIdx_SE <- array(sim_env$ObsFishIdx_SE[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseFishIdx <- array(0, dim = dim(ObsFishIdx))
  UseFishIdx[!is.na(ObsFishIdx) & ObsFishIdx > 0] <- 1

  # Population-specific fishery indices
  ObsFishIdx_pop <- array(sim_env$ObsFishIdx_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))
  UseFishIdx_pop <- array(0, dim = dim(ObsFishIdx_pop))
  UseFishIdx_pop[!is.na(ObsFishIdx_pop) & ObsFishIdx_pop > 0] <- 1
  ObsFishIdx_pop_SE <- array(sim_env$ObsFishIdx_pop_SE[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Retained Fishery Compositions
  ObsFishAgeComps <- array(sim_env$ObsFishAgeComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps <- array(sim_env$ISS_FishAgeComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps <- apply(ObsFishAgeComps, c(1,2,3,6), sum)
  UseFishAgeComps[!is.na(UseFishAgeComps) & UseFishAgeComps > 0] <- 1

  ObsFishLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps <- apply(ObsFishLenComps, c(1,2,3,6), sum)
    UseFishLenComps[!is.na(UseFishLenComps) & UseFishLenComps > 0] <- 1
  } else UseFishLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Population-specific retained fishery compositions
  ObsFishAgeComps_pop <- array(sim_env$ObsFishAgeComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_pop <- array(sim_env$ISS_FishAgeComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_pop <- apply(ObsFishAgeComps_pop, c(1,2,3,4,7), sum)
  UseFishAgeComps_pop[!is.na(UseFishAgeComps_pop) & UseFishAgeComps_pop > 0] <- 1

  ObsFishLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_pop <- apply(ObsFishLenComps_pop, c(1,2,3,4,7), sum)
    UseFishLenComps_pop[!is.na(UseFishLenComps_pop) & UseFishLenComps_pop > 0] <- 1
  } else UseFishLenComps_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Discarded Fishery Compositions
  ObsFishAgeComps_discard <- array(sim_env$ObsFishAgeComps_discard[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_discard <- array(sim_env$ISS_FishAgeComps_discard[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_discard <- apply(ObsFishAgeComps_discard, c(1,2,3,6), sum)
  UseFishAgeComps_discard[!is.na(UseFishAgeComps_discard) & UseFishAgeComps_discard > 0] <- 1

  ObsFishLenComps_discard <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_discard[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_discard <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_discard[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_discard <- apply(ObsFishLenComps_discard, c(1,2,3,6), sum)
    UseFishLenComps_discard[!is.na(UseFishLenComps_discard) & UseFishLenComps_discard > 0] <- 1
  } else UseFishLenComps_discard <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Population-specific discarded fishery compositions
  ObsFishAgeComps_discard_pop <- array(sim_env$ObsFishAgeComps_discard_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_fish_fleets))
  ISS_FishAgeComps_discard_pop <- array(sim_env$ISS_FishAgeComps_discard_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  UseFishAgeComps_discard_pop <- apply(ObsFishAgeComps_discard_pop, c(1,2,3,4,7), sum)
  UseFishAgeComps_discard_pop[!is.na(UseFishAgeComps_discard_pop) & UseFishAgeComps_discard_pop > 0] <- 1

  ObsFishLenComps_discard_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsFishLenComps_discard_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  ISS_FishLenComps_discard_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_FishLenComps_discard_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_fish_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseFishLenComps_discard_pop <- apply(ObsFishLenComps_discard_pop, c(1,2,3,4,7), sum)
    UseFishLenComps_discard_pop[!is.na(UseFishLenComps_discard_pop) & UseFishLenComps_discard_pop > 0] <- 1
  } else UseFishLenComps_discard_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_fish_fleets))

  # Survey Indices
  ObsSrvIdx <- array(sim_env$ObsSrvIdx[,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  ObsSrvIdx_SE <- array(sim_env$ObsSrvIdx_SE[,1:y,,, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  UseSrvIdx <- array(0, dim = dim(ObsSrvIdx))
  UseSrvIdx[!is.na(ObsSrvIdx) & ObsSrvIdx > 0] <- 1

  # Population-specific survey indices
  ObsSrvIdx_pop <- array(sim_env$ObsSrvIdx_pop[,,1:y,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))
  UseSrvIdx_pop <- array(0, dim = dim(ObsSrvIdx_pop))
  UseSrvIdx_pop[!is.na(ObsSrvIdx_pop) & ObsSrvIdx_pop > 0] <- 1
  ObsSrvIdx_pop_SE <- array(sim_env$ObsSrvIdx_pop_SE[,,1:y,,, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Survey Compositions
  ObsSrvAgeComps <- array(sim_env$ObsSrvAgeComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
  ISS_SrvAgeComps <- array(sim_env$ISS_SrvAgeComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  UseSrvAgeComps <- apply(ObsSrvAgeComps, c(1,2,3,6), sum)
  UseSrvAgeComps[!is.na(UseSrvAgeComps) & UseSrvAgeComps > 0] <- 1

  ObsSrvLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsSrvLenComps[,1:y,,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  ISS_SrvLenComps <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_SrvLenComps[,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseSrvLenComps <- apply(ObsSrvLenComps, c(1,2,3,6), sum)
    UseSrvLenComps[!is.na(UseSrvLenComps) & UseSrvLenComps > 0] <- 1
  } else UseSrvLenComps <- array(0, dim = c(sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Population-specific survey compositions
  ObsSrvAgeComps_pop <- array(sim_env$ObsSrvAgeComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, dim(sim_env$AgeingError)[3], sim_env$n_sexes, sim_env$n_srv_fleets))
  ISS_SrvAgeComps_pop <- array(sim_env$ISS_SrvAgeComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  UseSrvAgeComps_pop <- apply(ObsSrvAgeComps_pop, c(1,2,3,4,7), sum)
  UseSrvAgeComps_pop[!is.na(UseSrvAgeComps_pop) & UseSrvAgeComps_pop > 0] <- 1

  ObsSrvLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ObsSrvLenComps_pop[,,1:y,,,,,sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_lens, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  ISS_SrvLenComps_pop <- if(!is.null(sim_env$n_lens)) {
    array(sim_env$ISS_SrvLenComps_pop[,,1:y,,,, sim, drop = FALSE], dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_sexes, sim_env$n_srv_fleets))
  } else NULL
  if(!is.null(sim_env$n_lens)) {
    UseSrvLenComps_pop <- apply(ObsSrvLenComps_pop, c(1,2,3,4,7), sum)
    UseSrvLenComps_pop[!is.na(UseSrvLenComps_pop) & UseSrvLenComps_pop > 0] <- 1
  } else UseSrvLenComps_pop <- array(0, dim = c(sim_env$n_pop, sim_env$n_regions, length(1:y), sim_env$n_seas, sim_env$n_srv_fleets))

  # Return
  return(list(
    # Biologicals
    WAA = WAA,
    WAA_fish = WAA_fish,
    WAA_srv = WAA_srv,
    MatAA = MatAA,
    SizeAgeTrans = SizeAgeTrans,
    AgeingError = AgeingError,

    # Tagging
    use_conv_fish_tagging = sim_env$use_conv_fish_tagging,
    conv_tag_release_indicator = conv_tag_release_indicator,
    obs_conv_tag_fish_recap = obs_conv_tag_fish_recap,
    conv_tagged_fish = conv_tagged_fish,
    conv_tagged_fish_attr = conv_tagged_fish_attr,
    n_tag_cohorts = n_tag_cohorts,

    # Aggregated catches
    ObsCatch = ObsCatch,
    ln_sigmaC = ln_sigmaC,
    UseCatch = UseCatch,
    ObsDiscard = ObsDiscard,
    ln_sigmaD = ln_sigmaD,
    UseDiscard = UseDiscard,

    # Population-specific catches
    ObsCatch_pop = ObsCatch_pop,
    ln_sigmaC_pop = ln_sigmaC_pop,
    UseCatch_pop = UseCatch_pop,
    ObsDiscard_pop = ObsDiscard_pop,
    ln_sigmaD_pop = ln_sigmaD_pop,
    UseDiscard_pop = UseDiscard_pop,

    # Aggregated fishery indices
    ObsFishIdx = ObsFishIdx,
    ObsFishIdx_SE = ObsFishIdx_SE,
    UseFishIdx = UseFishIdx,

    # Population-specific fishery indices
    ObsFishIdx_pop = ObsFishIdx_pop,
    ObsFishIdx_pop_SE = ObsFishIdx_pop_SE,
    UseFishIdx_pop = UseFishIdx_pop,

    # Aggregated retained fishery compositions
    ObsFishAgeComps = ObsFishAgeComps,
    ISS_FishAgeComps = ISS_FishAgeComps,
    UseFishAgeComps = UseFishAgeComps,
    ObsFishLenComps = ObsFishLenComps,
    ISS_FishLenComps = ISS_FishLenComps,
    UseFishLenComps = UseFishLenComps,

    # Population-specific retained fishery compositions
    ObsFishAgeComps_pop = ObsFishAgeComps_pop,
    ISS_FishAgeComps_pop = ISS_FishAgeComps_pop,
    UseFishAgeComps_pop = UseFishAgeComps_pop,
    ObsFishLenComps_pop = ObsFishLenComps_pop,
    ISS_FishLenComps_pop = ISS_FishLenComps_pop,
    UseFishLenComps_pop = UseFishLenComps_pop,

    # Aggregated discarded fishery compositions
    ObsFishAgeComps_discard = ObsFishAgeComps_discard,
    ISS_FishAgeComps_discard = ISS_FishAgeComps_discard,
    UseFishAgeComps_discard = UseFishAgeComps_discard,
    ObsFishLenComps_discard = ObsFishLenComps_discard,
    ISS_FishLenComps_discard = ISS_FishLenComps_discard,
    UseFishLenComps_discard = UseFishLenComps_discard,

    # Population-specific discarded fishery compositions
    ObsFishAgeComps_discard_pop = ObsFishAgeComps_discard_pop,
    ISS_FishAgeComps_discard_pop = ISS_FishAgeComps_discard_pop,
    UseFishAgeComps_discard_pop = UseFishAgeComps_discard_pop,
    ObsFishLenComps_discard_pop = ObsFishLenComps_discard_pop,
    ISS_FishLenComps_discard_pop = ISS_FishLenComps_discard_pop,
    UseFishLenComps_discard_pop = UseFishLenComps_discard_pop,

    # Aggregated survey indices
    ObsSrvIdx = ObsSrvIdx,
    ObsSrvIdx_SE = ObsSrvIdx_SE,
    UseSrvIdx = UseSrvIdx,

    # Population-specific survey indices
    ObsSrvIdx_pop = ObsSrvIdx_pop,
    ObsSrvIdx_pop_SE = ObsSrvIdx_pop_SE,
    UseSrvIdx_pop = UseSrvIdx_pop,

    # Aggregated survey compositions
    ObsSrvAgeComps = ObsSrvAgeComps,
    ISS_SrvAgeComps = ISS_SrvAgeComps,
    UseSrvAgeComps = UseSrvAgeComps,
    ObsSrvLenComps = ObsSrvLenComps,
    ISS_SrvLenComps = ISS_SrvLenComps,
    UseSrvLenComps = UseSrvLenComps,

    # Population-specific survey compositions
    ObsSrvAgeComps_pop = ObsSrvAgeComps_pop,
    ISS_SrvAgeComps_pop = ISS_SrvAgeComps_pop,
    UseSrvAgeComps_pop = UseSrvAgeComps_pop,
    ObsSrvLenComps_pop = ObsSrvLenComps_pop,
    ISS_SrvLenComps_pop = ISS_SrvLenComps_pop,
    UseSrvLenComps_pop = UseSrvLenComps_pop
  ))

}

#' Predict fishery and discarded ISS under projected fishing mortality
#'
#' Scales fishery input sample sizes for the projection year \code{y} based
#' on the relationship between fishing mortality and historical ISS values.
#' For each region–sex–fleet cell, the minimum and maximum ISS from the
#' conditioning period (\code{1:(y-1)}) are identified from years with
#' positive, non-NA values, and the projected ISS is obtained by linear
#' interpolation between those bounds using the ratio of projected
#' \eqn{F_y} to the historical maximum \eqn{F} (capped at 1). If no valid
#' historical observations exist for a cell, ISS is set to zero. If
#' conditions for scaling are not met (e.g., maximum historical \eqn{F = 0}),
#' the mean historical ISS is used as a fallback. All prior years
#' (\code{1:(y-1)}) are carried over unchanged from \code{ISS_FishComps}.
#'
#' @param ISS_FishComps Array of fishery ISS values
#'   \code{[n_regions × n_yrs × n_seas × n_sexes × n_fish_fleets × n_sims]}.
#' @param Fmort Array of fishing mortality rates
#'   \code{[n_regions × n_yrs × n_seas × n_fish_fleets × n_sims]}.
#' @param y Integer. Projection year index for which ISS is predicted.
#' @param sim Integer. Simulation replicate index.
#' @param seas Integer. Season index.
#'
#' @return Array \code{[n_regions × y × 1 × n_sexes × n_fish_fleets]} with
#'   historical ISS values filled in for years \code{1:(y-1)} and the
#'   predicted ISS in year \code{y}.
#'
#'
#' @keywords internal
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
  n_fish_fleets <- dims[5]

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
