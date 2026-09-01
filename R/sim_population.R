# Operating model
#
# The operating model's true population dynamics: initial age structure,
# recruitment, biomass, and the annual cycle that applies mortality and movement.
# Simulate_Pop_Static is the entry point and drives the year loop, calling into
# sim_observations.R once the true state for a year exists.

#' Initialize age structure for a simulation replicate
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
#' density-dependence). Across sexes the draw follows
#' \code{InitDevs_sex_spec}: \code{"est_shared_s"} (the default) draws one
#' curve and gives it to every sex, \code{"est_all"} draws each sex its own.
#' If \code{ln_InitDevs_input} exists in the simulation environment, those
#' values are used directly rather than simulating new draws. Populations with
#' \code{R0 = 0} receive zero deviations. The equilibrium solver uses
#' \code{init_iter = n_ages × 5} iterations.
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
          tmp_ln_init_devs <- array(ln_InitDevs_input[p,r,,,sim], dim = c(n_ages - 1, n_sexes))
        } else { # simulate new initial age devs otherwise

          # get init devs devs
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])
          # Draw the deviations. Under est_shared_s one curve is drawn and every
          # sex reads it; under est_all each sex draws its own
          init_sex_spec <- if(exists("InitDevs_sex_spec")) InitDevs_sex_spec else "est_shared_s"
          if(is.null(tmp_ln_init_devs)) {
            n_dev_draws <- if(init_sex_spec == "est_all") n_sexes else 1
            init_draws <- stats::rnorm(n_dev_draws * (n_ages - 1), -exp(ln_sigmaR[1,p,sigma_idx])^2/2, exp(ln_sigmaR[1,p,sigma_idx]))
            tmp_ln_init_devs <- array(init_draws, dim = c(n_ages - 1, n_sexes)) # recycled across sexes when one curve was drawn
          }
        }

        # input age deviations
        if(R0[p,r,1,sim] != 0) {
          sim_env$ln_InitDevs[p,r,,,sim] <- tmp_ln_init_devs
        } else sim_env$ln_InitDevs[p,r,,,sim] <- 0

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
      ln_InitDevs = array(ln_InitDevs[,,,,sim], dim = c(n_pop, n_regions, n_ages - 1, n_sexes)), # initial deviations
      # Movement / mortality sequencing must match the estimation model, otherwise the
      # operating model's initial age structure is built under a different set of dynamics
      Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # rates in first year
      move_timing = move_timing
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
      ln_InitDevs = array(ln_InitDevs[,,,,sim], dim = c(n_pop, n_regions, n_ages - 1, n_sexes)), # initial deviations
      Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # rates in first year
      move_timing = move_timing
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
#' Unfished NAA (\code{NAA0}) is synchronized to match fished NAA at
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
#'   \code{$NAA[p, r, y, seas, 1, s, sim]}, and
#'   \code{$NAA0[p, r, y, seas, 1, s, sim]} are updated.
#' @param seas Integer. Season this recruitment first enters the population
#'   in, using \code{rec_seas_prop[p, seas, sim]}. Default \code{1}, matching
#'   the classic \code{rec_lag >= 1} case where recruitment for the whole year
#'   is already known before season 1 starts. \code{rec_lag = 0} (age-0
#'   recruitment) instead calls this with \code{seas = spawn_seas}, since that
#'   is the earliest season this year's own SSB (and hence recruitment) is
#'   knowable, see \code{\link{apply_pop_dy}}.
#'
#' @return \code{invisible(NULL)}. All modifications are made by reference
#'   within \code{sim_env}.
#'
#'
#' @keywords internal
generate_recruitment <- function(y,
                                 sim,
                                 sim_env,
                                 seas = 1) {

  sim_env$y    <- y
  sim_env$sim  <- sim
  sim_env$seas <- seas

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
                                       ret_sel = array(ret_sel[,,1,,,1,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # retained fishery selectivity in first year
                                       # Sequencing must match the estimation model; SSB0 and the SPR
                                       # machinery behind Beverton-Holt depend on it
                                       Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,1,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                       move_timing = move_timing

    ,
    expm_nsub = expm_nsub)


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

          # keep the deviation container consistent with the recruitment being
          # used, so an operating model conditioned on a fit reports the
          # deviations it is running on rather than zeros. Anything reading them,
          # a recruitment deviation index among them, then sees the same series
          # the population does.
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])
          if(tmp_det_rec[p,r] > 0 && tmp_total_rec > 0) {
            sim_env$ln_RecDevs[p,r,y,sim] <- log(tmp_total_rec / tmp_det_rec[p,r]) + exp(ln_sigmaR[2,p,sigma_idx])^2/2
          } else sim_env$ln_RecDevs[p,r,y,sim] <- 0
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

        # input recruitment into the season it first enters the population
        for(s in 1:n_sexes) sim_env$NAA[p,r,y,seas,1,s,sim] <- tmp_total_rec * rec_seas_prop[p,seas,sim] * sexratio[p,r,y,s,sim]

        sim_env$Rec[p,r,y,sim] <- tmp_total_rec # Save annual recruitment estimates
        sim_env$NAA0[p,r,y,seas,1,,sim] = NAA[p,r,y,seas,1,,sim] # populate unfished NAA

      } # end r loop
    } # end p loop
  })
}

#' Compute spawning-time biomass quantities for one simulation year/season
#'
#' Computes Total_Biom, SSB, Dynamic_SSB0, and eff_SSB for year \code{y} at
#' season \code{seas} (always called with \code{seas == spawn_seas}) from the
#' current \code{NAA}/\code{NAA0} state in \code{sim_env}. Factored out of
#' \code{\link{apply_pop_dy}} so it can be evaluated either before or after
#' that season's mortality/ageing step depending on \code{rec_lag}, without
#' duplicating the underlying math. Pure/read-only: returns a list rather
#' than modifying \code{sim_env}.
#'
#' @param y Year integer
#' @param seas Season integer
#' @param sim Simulation integer
#' @param sim_env Simulation environment
#'
#' @keywords internal
compute_biom_y_sim <- function(y, seas, sim, sim_env) {

  NAA <- sim_env$NAA; NAA0 <- sim_env$NAA0
  WAA <- sim_env$WAA; MatAA <- sim_env$MatAA; ZAA <- sim_env$ZAA
  natmort <- sim_env$natmort; t_spawn <- sim_env$t_spawn; seasdur <- sim_env$seasdur
  n_pop <- sim_env$n_pop; n_regions <- sim_env$n_regions; n_seas <- sim_env$n_seas
  n_ages <- sim_env$n_ages; n_sexes <- sim_env$n_sexes
  natal_region <- sim_env$natal_region; stray_rate <- sim_env$stray_rate
  sgl_seas_spawning_movement <- sim_env$sgl_seas_spawning_movement
  Movement <- sim_env$Movement; Mrate <- sim_env$Mrate
  move_timing <- if(is.null(sim_env$move_timing)) 0 else sim_env$move_timing
  expm_nsub <- if(is.null(sim_env$expm_nsub)) 0 else sim_env$expm_nsub
  do_recruits_move <- sim_env$do_recruits_move

  tmp_NAA_spawn <- NAA[,,y,seas,,,sim, drop = FALSE]
  tmp_NAA0_spawn <- NAA0[,,y,seas,,,sim, drop = FALSE]

  # Propagate to the spawning point
  if(move_timing != 0 && n_regions > 1) {
    for(p in 1:n_pop) {
      for(a in 1:n_ages) {
        moves <- (do_recruits_move == 1 || a > 1)
        for(s in 1:n_sexes) {
          Mv <- if(moves) Movement[p,,,y,seas,a,s,sim] else diag(n_regions)
          Qv <- if(moves) Mrate[p,,,y,seas,a,s,sim] else matrix(0, n_regions, n_regions)
          tmp_NAA_spawn[p,,1,1,a,s,1] <- spawn_state(tmp_NAA_spawn[p,,1,1,a,s,1], Mv,
                                                     ZAA[p,,y,seas,a,s,sim], Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
          tmp_NAA0_spawn[p,,1,1,a,s,1] <- spawn_state(tmp_NAA0_spawn[p,,1,1,a,s,1], Mv,
                                                      natmort[p,,y,a,s,sim] * seasdur[seas], Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
        } # end s loop
      } # end a loop
    } # end p loop
  }

  # If we we are natal homing with 1 season
  if(n_seas == 1 && n_pop > 1) {
    for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
      tmp_NAA_spawn[p,,1,1,a,s,1] <- tmp_NAA_spawn[p,,1,1,a,s,1] %*% sgl_seas_spawning_movement[p,,,y,a,s,sim]
      tmp_NAA0_spawn[p,,1,1,a,s,1] <- tmp_NAA0_spawn[p,,1,1,a,s,1] %*% sgl_seas_spawning_movement[p,,,y,a,s,sim]
    } # end s loop
  }

  # Mortality discount up to spawning. Already folded into spawn_state above when
  # move_timing != 0, so it collapses to 1 in that case to avoid applying it twice.
  spawn_disc_Z <- if(move_timing == 0 || n_regions == 1) exp(-ZAA[,,y,seas,,,sim,drop = FALSE] * t_spawn) else 1
  spawn_disc_Z_f <- if(move_timing == 0 || n_regions == 1) exp(-ZAA[,, y, seas, , 1, sim,drop = FALSE] * t_spawn) else 1

  # Total Biomass
  Total_Biom_y <- apply(tmp_NAA_spawn *
                          WAA[,, y, seas, , , sim,drop = FALSE] *
                          spawn_disc_Z, c(1,2), sum)

  # Spawning Stock Biomass
  SSB_y <- apply(tmp_NAA_spawn[,, , , , 1, 1,drop = FALSE] *
                   WAA[,, y, seas, , 1, sim,drop = FALSE] *
                   MatAA[,, y, seas, , 1, sim,drop = FALSE] *
                   spawn_disc_Z_f, c(1,2), sum)

  # Get dynamic B0
  SSB0_array <- tmp_NAA0_spawn[,, , , , 1, 1,drop = FALSE] *  WAA[,,  y, seas, , 1, sim, drop = FALSE] * MatAA[,,y, seas, , 1, sim, drop = FALSE]
  if(move_timing == 0 || n_regions == 1) {
    mort_spawn <- exp(-natmort[,, y, , 1, sim, drop = FALSE] * t_spawn * seasdur[seas])
    mort_spawn <- array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
  } else mort_spawn <- 1
  Dynamic_SSB0_y <- apply(SSB0_array * mort_spawn, c(1,2), sum) # Dynamic B0

  if(n_sexes == 1) { # If single sex model, multiply SSB calculations by 0.5
    SSB_y <- SSB_y * 0.5
    Dynamic_SSB0_y <- Dynamic_SSB0_y * 0.5
  }

  # Accumulate effective SSB at each population's natal region
  # across all source populations (captures stray contributions)
  eff_SSB_y <- array(0, dim = n_pop)
  if(n_pop > 1) {
    n_pop_in_region = array(0, dim = n_regions)
    for(p in 1:n_pop) n_pop_in_region[natal_region[p]] = n_pop_in_region[natal_region[p]] + 1
    for(p2 in 1:n_pop) {
      for(p in 1:n_pop) {
        if(p == p2) {
          eff_SSB_y[p2] = eff_SSB_y[p2] + SSB_y[p, natal_region[p2]]
        } else {
          n_receivers = n_pop_in_region[natal_region[p2]]
          eff_SSB_y[p2] = eff_SSB_y[p2] + (stray_rate[p,y,sim] / n_receivers) * SSB_y[p, natal_region[p2]]
        }
      }
    }
  } else eff_SSB_y[1] = sum(SSB_y[1,])

  list(Total_Biom_y = Total_Biom_y, SSB_y = SSB_y, Dynamic_SSB0_y = Dynamic_SSB0_y, eff_SSB_y = eff_SSB_y)
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
#' When \code{rec_lag == 0} (age-0 recruitment), this year's recruitment
#' can't be known until \code{spawn_seas} is reached (it depends on this
#' year's own SSB), so \code{\link{generate_recruitment}} is called from
#' inside this function at \code{seas == spawn_seas} instead of beforehand -
#' see the "rec_lag == 0" block below, which mirrors the equivalent
#' restructuring in the estimation model (\code{SPoRC_rtmb.R}).
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

      # apportion recruitment across seasons already known from earlier this
      # year:
      # - rec_lag != 0: the year's recruitment is already known (computed by
      #   generate_recruitment() before this function ran), so any season
      #   past the first gets its share here, as before.
      # - rec_lag == 0: recruitment isn't known until spawn_seas is reached
      #   (below), so only seasons strictly after spawn_seas are handled here;
      #   spawn_seas itself generates and inserts its own share.
      if(if(rec_lag != 0) seas > 1 else seas > spawn_seas) {
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
      }

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

      # Movement is applied at the start of the season only under move_timing == 0; under
      # timings 1 and 2 it is folded into the mortality/ageing step at the end of the season.
      if(n_regions > 1 && move_timing == 0) {
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

      # Compute Biomass Quantities + Recruitment (rec_lag == 0 only) --------
      if(rec_lag == 0 && seas == spawn_seas) {

        # SSB from survivors only, used to feed generate_recruitment() below.
        spawn_biom <- compute_biom_y_sim(y, seas, sim, sim_env)
        sim_env$SSB[,, y, sim] <- spawn_biom$SSB_y

        generate_recruitment(y, sim, sim_env, seas = spawn_seas)

        # Recruits just generated above missed this season's movement step so apply movement if needed.
        # Only under move_timing == 0; timings 1 and 2 pick them up in the end-of-season step.
        if(do_recruits_move == 1 && n_regions > 1 && move_timing == 0) {
          for(p in 1:n_pop) {
            for(s in 1:n_sexes) {
              sim_env$NAA[p,,y,seas,1,s,sim] <- t(NAA[p,,y,seas,1,s,sim]) %*% Movement[p,,,y,seas,1,s,sim]
              sim_env$NAA0[p,,y,seas,1,s,sim] <- t(NAA0[p,,y,seas,1,s,sim]) %*% Movement[p,,,y,seas,1,s,sim]
            } # end s loop
          } # end p loop
          NAA_aft[,,y,seas,1,,sim] <- NAA[,,y,seas,1,,sim]
        }

        # Recompute now that this year's recruits are included
        spawn_biom <- compute_biom_y_sim(y, seas, sim, sim_env)
        sim_env$Total_Biom[,, y, sim] <- spawn_biom$Total_Biom_y
        sim_env$SSB[,, y, sim] <- spawn_biom$SSB_y
        sim_env$Dynamic_SSB0[,,y,sim] <- spawn_biom$Dynamic_SSB0_y
        sim_env$eff_SSB[,y,sim] <- spawn_biom$eff_SSB_y

      } # end if rec_lag == 0 && seas == spawn_seas

      # Post-season state at every age, before the ageing shift. Under move_timing == 0
      # movement was applied above, so this reduces to the original elementwise survival.
      if(move_timing == 0 || n_regions == 1) {
        # array() guards against R dropping length-1 pop/region/sex dimensions
        sstep_NAA <- array(NAA[,,y,seas,1:n_ages,,sim] * exp(-ZAA[,,y,seas,1:n_ages,,sim]),
                           dim = c(n_pop, n_regions, n_ages, n_sexes))
        sstep_NAA0 <- array(NAA0[,,y,seas,1:n_ages,,sim] * exp(-(tmp_natmort[,,1:n_ages,])),
                            dim = c(n_pop, n_regions, n_ages, n_sexes))
      } else {
        sstep_NAA <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        sstep_NAA0 <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        for(p in 1:n_pop) {
          for(a in 1:n_ages) {
            moves <- (do_recruits_move == 1 || a > 1)
            for(s in 1:n_sexes) {
              Mv <- if(moves) Movement[p,,,y,seas,a,s,sim] else diag(n_regions)
              Qv <- if(moves) Mrate[p,,,y,seas,a,s,sim] else matrix(0, n_regions, n_regions)
              sstep_NAA[p,,a,s] <- advance_seas(NAA[p,,y,seas,a,s,sim], Mv, ZAA[p,,y,seas,a,s,sim],
                                                Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
              sstep_NAA0[p,,a,s] <- advance_seas(NAA0[p,,y,seas,a,s,sim], Mv, tmp_natmort[p,,a,s],
                                                 Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
            } # end s loop
          } # end a loop
        } # end p loop
        # Movement happens at season end under these timings, so record it here
        NAA_aft[,,y,seas,,,sim] <- sstep_NAA
      }

      if(seas < n_seas) { # Within year seasonal mortality
        sim_env$NAA[,,y,seas + 1,1:n_ages,,sim] = sstep_NAA # fished
        sim_env$NAA0[,,y,seas + 1,1:n_ages,,sim] <- sstep_NAA0 # unfished
      } else {
        # Advance into the next year, season 1
        sim_env$NAA[,,y+1,1,2:n_ages,,sim] = sstep_NAA[,,1:(n_ages-1),] # fished
        sim_env$NAA[,,y+1,1,n_ages,,sim] = NAA[,,y+1,1,n_ages,,sim] + sstep_NAA[,,n_ages,] # Acuumulate plus group (fished)
        sim_env$NAA0[,,y+1,1,2:n_ages,,sim] = sstep_NAA0[,,1:(n_ages-1),] # fished
        sim_env$NAA0[,,y+1,1,n_ages,,sim] = NAA0[,,y+1,1,n_ages,,sim] + sstep_NAA0[,,n_ages,] # Acuumulate plus group (unfished)

        # State-space numbers at age, applied where the estimation model applies it: after the plus
        # group accumulates, at the year boundary, with the unfished numbers taking the same factor
        if(sim_env$NAA_re > 0 && (y + 1) <= n_yrs) {
          sim_env$NAA_pred[,,y+1,,,sim] = sim_env$NAA[,,y+1,1,,,sim]
          fac <- exp(sim_env$naa_eta[,,y+1,,])
          sim_env$NAA[,,y+1,1,,,sim] = sim_env$NAA[,,y+1,1,,,sim] * fac
          sim_env$NAA0[,,y+1,1,,,sim] = sim_env$NAA0[,,y+1,1,,,sim] * fac
        }
      }

      # Compute Biomass Quantities (rec_lag != 0: unchanged original timing)
      if(rec_lag != 0 && seas == spawn_seas) {
        spawn_biom <- compute_biom_y_sim(y, seas, sim, sim_env)
        sim_env$Total_Biom[,, y, sim] <- spawn_biom$Total_Biom_y
        sim_env$SSB[,, y, sim] <- spawn_biom$SSB_y
        sim_env$Dynamic_SSB0[,,y,sim] <- spawn_biom$Dynamic_SSB0_y
        sim_env$eff_SSB[,y,sim] <- spawn_biom$eff_SSB_y
      } # if season = spawning season
    } # end seas loop
  })
}


#' Run the annual cycle for a single simulation year
#'
#' Orchestrates the complete annual sequence of operating model processes for
#' year \code{y} and simulation replicate \code{sim}: initializes age
#' structure and generates first-year recruitment at \code{y = 1}; applies
#' population dynamics (movement, mortality, biomass); generates fishery
#' catches, indices, and compositions; generates survey indices and
#' compositions; releases conventional tags; generates fishery tag
#' recaptures (when any \code{use_conv_fish_tagging = 1}); and generates
#' recruitment for the following year (\code{y + 1}) when \code{y < n_yrs}.
#'
#' The two standalone \code{generate_recruitment()} calls described above
#' (at \code{y = 1} and for \code{y + 1}) only run when \code{rec_lag != 0}.
#' For \code{rec_lag = 0} (age-0 recruitment), recruitment for year \code{y}
#' depends on year \code{y}'s own SSB, which isn't known until
#' \code{\link{apply_pop_dy}} reaches \code{spawn_seas} within that year -
#' \code{generate_recruitment()} is called from inside \code{apply_pop_dy()}
#' instead, once that SSB is available.
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
    # The whole replicate's innovations are drawn up front: the matrix is conditional on all
    # previous states and formed rectangular over age and year
    if(isTRUE(sim_env$NAA_re > 0)) {
      sim_env$naa_eta <- draw_naa_innovations(sim_env)
      sim_env$naa_eta_all[,,,,,sim] <- sim_env$naa_eta
    }
    generate_initial_age_structure(y = 1, sim, sim_env) # Initialize age structure
    if(sim_env$rec_lag != 0) generate_recruitment(y = 1, sim, sim_env) # Get recruitment in the first year
  }

  apply_pop_dy(y, sim, sim_env) # Apply population dynamics (movement, mortality, and biomass calculations)
  generate_fishery_catch_comp_idx(y, sim, sim_env) # Get Fishery Catches, Compositions, and Indices
  generate_survey_comp_idx(y, sim, sim_env) # Get Fishery Catches, Compositions, and Indices

  if(any(sim_env$use_conv_fish_tagging == 1)) {
    release_conv_tags(y, sim, sim_env) # Release conventional tags
    generate_fishery_conv_tags_recap(y, sim, sim_env) # Generate fishery conventional tag recaptures
  }

  if(y < sim_env$n_yrs && sim_env$rec_lag != 0) generate_recruitment(y = y + 1, sim, sim_env) # Get recruitment in the following year


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
                  naa_eta = sim_env$naa_eta_all,
                  NAA_pred = sim_env$NAA_pred,
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
                  Mrate = sim_env$Mrate,
                  move_timing = sim_env$move_timing,
                  expm_nsub = if(is.null(sim_env$expm_nsub)) 0 else sim_env$expm_nsub,
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
                  TrueCatchAA = sim_env$TrueCatchAA, ObsCatchAA = sim_env$ObsCatchAA,
                  TrueSrvIdxAA = sim_env$TrueSrvIdxAA, ObsSrvIdxAA = sim_env$ObsSrvIdxAA,
                  ObsCatch = sim_env$ObsCatch,
                  ObsFishIdx = sim_env$ObsFishIdx,
                  TrueFishIdx = sim_env$TrueFishIdx,
                  ObsFishIdx_SE = sim_env$ObsFishIdx_SE,
                  ObsFishAgeComps = sim_env$ObsFishAgeComps,
                  ObsFishLenComps = sim_env$ObsFishLenComps,
                  ObsFish_caal = sim_env$ObsFish_caal,
                  ObsSrv_caal = sim_env$ObsSrv_caal,
                  ISS_Fish_caal = sim_env$ISS_Fish_caal,
                  ISS_Srv_caal = sim_env$ISS_Srv_caal,
                  do_fish_caal = sim_env$do_fish_caal,
                  do_srv_caal = sim_env$do_srv_caal,

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
                  SizeAgeTrans_fish = if(!is.null(sim_env$SizeAgeTrans_fish)) sim_env$SizeAgeTrans_fish else NULL,
                  SizeAgeTrans_srv = if(!is.null(sim_env$SizeAgeTrans_srv)) sim_env$SizeAgeTrans_srv else NULL,
                  AgeingError = sim_env$AgeingError,
                  AgeingError_fish = if(!is.null(sim_env$AgeingError_fish)) sim_env$AgeingError_fish else NULL,
                  AgeingError_srv = if(!is.null(sim_env$AgeingError_srv)) sim_env$AgeingError_srv else NULL,
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
                  # Same quantity under both spellings. n_yrs is what the simulation
                  # code reads; n_years is the name used throughout the dimension
                  # documentation, so both are exposed for downstream scripts.
                  n_years = sim_env$n_yrs,
                  n_yrs = sim_env$n_yrs,
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
