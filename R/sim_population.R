# Operating model
#
# The operating model's true population dynamics. Simulate_Pop_Static drives the year loop and calls
# into sim_observations.R once the true state for a year exists.

# Initial Age Structure and Recruitment -------------------------------------

#' Initialize age structure for a simulation replicate
#'
#' Draws or reads the initial age deviations and calls
#' \code{\link{Get_Init_NAA}} for the fished and unfished equilibrium numbers at
#' age in year 1, season 1, writing them into \code{NAA} and \code{NAA0}. Called
#' once per replicate at \code{y = 1} by \code{\link{run_annual_cycle}}.
#'
#' Sharing follows the estimation model: one draw per population when
#' \code{n_pop > 1}, or one per region when \code{n_pop = 1} and
#' \code{init_dd = 0}. Across sexes it follows \code{InitDevs_sex_spec}, with
#' \code{"est_shared_s"} (default) drawing one curve for every sex and
#' \code{"est_all"} drawing each its own. An \code{ln_InitDevs_input} in the
#' environment is used directly rather than drawn. Populations with \code{R0 = 0}
#' get zero deviations, and the equilibrium solver runs \code{n_ages × 5}
#' iterations.
#'
#' @param y Integer. Year index, which must be \code{1}.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}},
#'   modified in place: \code{$ln_InitDevs}, \code{$NAA[,,1,1,,,sim]} and
#'   \code{$NAA0[,,1,1,,,sim]}.
#'
#' @return \code{invisible(NULL)}; everything is modified by reference within
#'   \code{sim_env}.
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

          # get init devs
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])
          # Draw the deviations. Under est_shared_s one curve is drawn and every
          # sex reads it; under est_all each sex draws its own
          init_sex_spec <- if(exists("InitDevs_sex_spec")) InitDevs_sex_spec else "est_shared_s"
          if(is.null(tmp_ln_init_devs)) {
            n_dev_draws <- if(init_sex_spec == "est_all") n_sexes else 1
            init_center <- if(isTRUE(rec_bias_correct == 0)) 0 else -exp(ln_sigmaR[1,p,sigma_idx])^2 / 2 # whether to do bias correction
            init_draws <- stats::rnorm(n_dev_draws * (n_ages - 1), init_center, exp(ln_sigmaR[1,p,sigma_idx]))
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
      natmort = array(natmort[,,1,,,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)), # natural mortality in first year
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
      natmort = array(natmort[,,1,,,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)), # natural mortality in first year
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
#' Takes deterministic recruitment from \code{\link{Get_Det_Recruitment}},
#' multiplies it by lognormal deviations, apportions it across sexes and seasons,
#' and writes it into the age-one slot of \code{sim_env$NAA}, with \code{NAA0}
#' synchronized to match. A \code{Rec_input} covering year \code{y} overrides the
#' draw entirely.
#'
#' Deviation sharing follows \code{\link{generate_initial_age_structure}}: one
#' draw per population when \code{n_pop > 1}, or one per region when
#' \code{n_pop = 1} under local density dependence. Populations with
#' \code{R0 = 0} get zero deviations, and \code{sigma_idx} picks the natal
#' region's \code{ln_sigmaR} for the bias correction. \code{RecDevs_model} sets
#' what the draw is centered on: zero for independent deviations, the previous
#' year's for a random walk, and \code{RecDevs_rho} times it for an AR1. Only the
#' independent draws are bias corrected, a walk's deviation not being mean zero.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}},
#'   modified in place: \code{$ln_RecDevs}, \code{$Rec}, \code{$NAA} and
#'   \code{$NAA0}.
#' @param seas Integer. Season this recruitment first enters in, through
#'   \code{rec_seas_prop[p, seas, sim]}. Default \code{1}, the classic
#'   \code{rec_lag >= 1} case where the whole year's recruitment is known before
#'   season one. \code{rec_lag = 0} instead calls this with
#'   \code{seas = spawn_seas}, the earliest season this year's own SSB is
#'   knowable. See \code{\link{apply_pop_dy}}.
#'
#' @return \code{invisible(NULL)}; everything is modified by reference within
#'   \code{sim_env}.
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

                                       # Note: every input to unfished SSB0 is taken at SR_ref_yr, female quantities only
                                       sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,SR_ref_yr,1,sim], dim = c(n_pop, n_regions)),
                                       WAA = array(WAA[,,SR_ref_yr,,,1,sim], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                       MatAA = array(MatAA[,,SR_ref_yr,,,1,sim], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                       natmort = array(natmort[,,SR_ref_yr,,,1,sim], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                       stray_rate = array(stray_rate[,SR_ref_yr,sim], dim = c(n_pop)),
                                       Movement = array(Movement[,,,SR_ref_yr,,,1,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                       sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,SR_ref_yr,,1,sim], dim = c(n_pop, n_regions, n_regions, n_ages)),
                                       SSB_vals = array(SSB[,,,sim], dim = c(n_pop, n_regions, n_yrs)),
                                       n_fish_fleets = n_fish_fleets,
                                       t_spawn = t_spawn,
                                       n_seas = n_seas,
                                       spawn_seas = spawn_seas,
                                       seasdur = seasdur,
                                       do_recruits_move = do_recruits_move,
                                       init_F = init_F, # initial F applied
                                       dmr = array(dmr[,SR_ref_yr,,,sim], dim = c(n_regions, n_seas, n_fish_fleets)), # discard mortality rate
                                       fish_sel = array(fish_sel[,,SR_ref_yr,,,1,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # total fishery selectivity at SR_ref_yr
                                       ret_sel = array(ret_sel[,,SR_ref_yr,,,1,,sim], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # retained fishery selectivity at SR_ref_yr
                                       Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,SR_ref_yr,,,1,sim], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                       move_timing = move_timing,
                                       expm_nsub = expm_nsub)


    # if Rec_input exists and year index is within bounds
    use_rec_input <- exists("Rec_input") && (y <= dim(Rec_input)[3])
    tmp_ln_rec_devs <- NULL

    for(p in 1:n_pop) {

      # reset deviations for each population (draws for each popn)
      if(n_pop > 1) tmp_ln_rec_devs <- NULL

      for(r in 1:n_regions) {

        # if local DD and n_pop = 1, reset deviations for each region (draws for each region)
        if(n_pop == 1 && rec_dd == 0) tmp_ln_rec_devs <- NULL

        if(use_rec_input) { # if jut using recruitment input

          tmp_total_rec <- Rec_input[p,r,y,sim]
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])

          if(isTRUE(tmp_det_rec[p,r] > 0) && tmp_total_rec > 0) { # back out the true ln Rec Devs from a conditioned fit if needed
            sim_env$ln_RecDevs[p,r,y,sim] <- log(tmp_total_rec / tmp_det_rec[p,r]) + exp(ln_sigmaR[2,p,sigma_idx])^2 / 2
          } else sim_env$ln_RecDevs[p,r,y,sim] <- 0

        } else {

          # get rec devs
          sigma_idx <- ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])

          # doing random walk or AR 1
          dev_mu <- 0
          dev_sd <- exp(ln_sigmaR[2, p, sigma_idx])
          if(RecDevs_model != 1 && y > 1) {
            prev_dev <- sim_env$ln_RecDevs[p,r,y - 1,sim]
            dev_mu <- if(RecDevs_model == 2) prev_dev else RecDevs_rho[p,r] * prev_dev
          }

          # get staionary marginal for ar1
          if(RecDevs_model == 3 && y == 1) dev_sd <- dev_sd / sqrt(1 - RecDevs_rho[p,r]^2)

          # if using a dsem, rec devs already drawn, otherwise, draw here
          dsem_cell <- dsem_drawn$ln_RecDevs[p,r,y]
          tmp_ln_rec_devs <- if(dsem_cell) sim_env$ln_RecDevs[p,r,y,sim] else stats::rnorm(1, dev_mu, dev_sd)

          # input devs here
          if(R0[p,r,y,sim] != 0) {
            sim_env$ln_RecDevs[p,r,y,sim] <- tmp_ln_rec_devs
          } else sim_env$ln_RecDevs[p,r,y,sim] <- 0

          # doing bias correction
          bias_corr <- if(dsem_cell || isTRUE(rec_bias_correct == 0)) 0 else if(RecDevs_model == 1) exp(ln_sigmaR[2,p,sigma_idx])^2 / 2 else 0
          tmp_total_rec <- tmp_det_rec[p,r] * exp(sim_env$ln_RecDevs[p,r,y,sim] - bias_corr)

        }

        # input recruitment into the season it first enters the population
        for(s in 1:n_sexes) sim_env$NAA[p,r,y,seas,1,s,sim] <- tmp_total_rec * rec_seas_prop[p,seas,sim] * sexratio[p,r,y,s,sim]

        sim_env$Rec[p,r,y,sim] <- tmp_total_rec # Save annual recruitment estimates
        sim_env$NAA0[p,r,y,seas,1,,sim] = NAA[p,r,y,seas,1,,sim] # populate unfished NAA

      } # end r loop
    } # end p loop
  })
}

# Biomass -------------------------------------------------------------------

#' Compute spawning-time biomass quantities for one simulation year/season
#'
#' The operating model's state at year \code{y} and season \code{seas}, always the spawning
#' season, sliced at replicate \code{sim} and given to \code{\link{biom_at_spawn}}, which the
#' estimation model and the forward projection also run. 
#'
#' @param y Year integer
#' @param seas Season integer
#' @param sim Simulation integer
#' @param sim_env Simulation environment
#'
#' @keywords internal
compute_biom_y_sim <- function(y, seas, sim, sim_env) {

  n_pop <- sim_env$n_pop
  n_regions <- sim_env$n_regions
  n_seas <- sim_env$n_seas
  move_timing <- if(is.null(sim_env$move_timing)) 0 else sim_env$move_timing

  # drop sim dimension
  drop_sim <- function(x) {
    dim(x) <- dim(x)[-length(dim(x))]
    x
  }

  propagates <- move_timing != 0 && n_regions > 1 # whether movement happens
  homes <- n_seas == 1 && n_pop > 1 # whether homing happens

  biom_at_spawn(
    NAA_s = drop_sim(sim_env$NAA[,,y,seas,,,sim, drop = FALSE]),
    NAA0_s = drop_sim(sim_env$NAA0[,,y,seas,,,sim, drop = FALSE]),
    WAA_s = drop_sim(sim_env$WAA[,,y,seas,,,sim, drop = FALSE]),
    MatAA_s = drop_sim(sim_env$MatAA[,,y,seas,,,sim, drop = FALSE]),
    ZAA_s = drop_sim(sim_env$ZAA[,,y,seas,,,sim, drop = FALSE]),
    natmort_s = drop_sim(sim_env$natmort[,,y,seas,,,sim, drop = FALSE]),
    spawn_move_s = if(homes) drop_sim(sim_env$sgl_seas_spawning_movement[,,,y,,,sim, drop = FALSE]) else NULL,
    stray_rate_y = if(n_pop > 1) sim_env$stray_rate[,y,sim] else NULL,
    t_spawn = sim_env$t_spawn,
    seasdur_seas = sim_env$seasdur[seas],
    n_seas = n_seas,
    n_pop = n_pop,
    n_regions = n_regions,
    n_ages = sim_env$n_ages,
    n_sexes = sim_env$n_sexes,
    natal_region = sim_env$natal_region,
    Movement_s = if(propagates) drop_sim(sim_env$Movement[,,,y,seas,,,sim, drop = FALSE]) else NULL,
    Mrate_s = if(propagates) drop_sim(sim_env$Mrate[,,,y,seas,,,sim, drop = FALSE]) else NULL,
    move_timing = move_timing,
    do_recruits_move = sim_env$do_recruits_move,
    expm_nsub = if(is.null(sim_env$expm_nsub)) 0 else sim_env$expm_nsub
  )
}

# Annual Cycle --------------------------------------------------------------

#' Apply population dynamics within a simulation year
#'
#' Runs the within-year loop for year \code{y}: seasonal recruitment
#' apportionment from season two on, movement, Baranov mortality, age advancement
#' into the following year, and the spawning-season biomass quantities (total
#' biomass, SSB, dynamic \eqn{B_0} and effective SSB under natal homing). The
#' fished and unfished trajectories are tracked together, with the snapshots
#' before and after movement stored in \code{NAA_bef} and \code{NAA_aft}. Movement
#' runs only when \code{n_regions > 1}, and recruits are left out of it when
#' \code{do_recruits_move = 0}. With one season and several populations,
#' \code{sgl_seas_spawning_movement} is applied to both trajectories before the
#' biomass quantities; with one sex, SSB and \eqn{B_0} are halved to give
#' female-only spawning biomass.
#'
#' Under \code{rec_lag == 0} this year's recruitment is not knowable until
#' \code{spawn_seas}, since it depends on this year's own SSB, so
#' \code{\link{generate_recruitment}} is called from inside this function at
#' \code{seas == spawn_seas}, mirroring the estimation model.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}},
#'   modified in place: \code{$ZAA}, \code{$NAA}, \code{$NAA0}, \code{$NAA_bef},
#'   \code{$NAA_aft}, \code{$Total_Biom}, \code{$SSB}, \code{$Dynamic_SSB0} and
#'   \code{$eff_SSB}.
#'
#' @return \code{invisible(NULL)}; everything is modified by reference within
#'   \code{sim_env}.
#'
#' @keywords internal
apply_pop_dy <- function(y, sim, sim_env) {

  sim_env$y   <- y
  sim_env$sim <- sim

  with(sim_env, {

    for(seas in 1:n_seas) {

      # apportion recruitment across seasons already known from earlier this year.
      # under rec_lag == 0 spawn_seas generates and inserts its own share below
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
      tmp_natmort <- array((natmort[,,y,seas,,,sim] * seasdur[seas]), dim = c(n_pop, n_regions, n_ages, n_sexes))

      for(p in 1:n_pop) {

        # Get retained catch selectivity
        tmp_fish_sel <- array(fish_sel[p,,y,seas,,,,sim], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)) # total selectivtiy
        tmp_ret_sel <- array(ret_sel[p,,y,seas,,,,sim], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)) # retained selectivity
        tmp_ret_FAA <- apply(sweep(tmp_fish_sel * tmp_ret_sel, c(1,4), tmp_Fmort, "*"), c(1,2,3), sum) # apply Frate to retained selectivity
        tmp_ret_FAA <- array(tmp_ret_FAA, dim = c(n_regions, n_ages, n_sexes)) # reshape

        # Get discarded catch selectivity
        tmp_disc_FAA <- apply(sweep(tmp_fish_sel * (1 - tmp_ret_sel), c(1,4), tmp_Fmort * tmp_dmr, "*"), c(1,2,3), sum) # apply Frate and dmr to discarded selectivity

        # Get natural mortality
        tmp_MAA <- array(tmp_natmort[p,,,,drop = FALSE], dim = c(n_regions, n_ages, n_sexes)) # reshape natural mortality

        # Get total mortality
        sim_env$ZAA[p,,y,seas,,,sim] <- tmp_MAA + tmp_ret_FAA + tmp_disc_FAA
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

        # State-space numbers at age at a within-year boundary, where the estimation model
        # applies it: on the survival and movement step, with no ageing shift
        if(sim_env$NAA_re > 0) {
          sim_env$NAA_pred[,,y,seas + 1,,,sim] = sim_env$NAA[,,y,seas + 1,,,sim]
          fac <- exp(sim_env$naa_eta[,,y,seas + 1,,])
          sim_env$NAA[,,y,seas + 1,,,sim] = sim_env$NAA[,,y,seas + 1,,,sim] * fac
          sim_env$NAA0[,,y,seas + 1,,,sim] = sim_env$NAA0[,,y,seas + 1,,,sim] * fac
        }
      } else {
        # Advance into the next year, season 1
        sim_env$NAA[,,y + 1,1,2:n_ages,,sim] = sstep_NAA[,,1:(n_ages - 1),] # fished
        sim_env$NAA[,,y + 1,1,n_ages,,sim] = NAA[,,y + 1,1,n_ages,,sim] + sstep_NAA[,,n_ages,] # Acuumulate plus group (fished)
        sim_env$NAA0[,,y + 1,1,2:n_ages,,sim] = sstep_NAA0[,,1:(n_ages - 1),] # fished
        sim_env$NAA0[,,y + 1,1,n_ages,,sim] = NAA0[,,y + 1,1,n_ages,,sim] + sstep_NAA0[,,n_ages,] # Acuumulate plus group (unfished)

        # State-space numbers at age, applied where the estimation model applies it: after the plus
        # group accumulates, at the year boundary, with the unfished numbers taking the same factor
        if(sim_env$NAA_re > 0 && (y + 1) <= n_yrs) {
          sim_env$NAA_pred[,,y + 1,1,,,sim] = sim_env$NAA[,,y + 1,1,,,sim]
          fac <- exp(sim_env$naa_eta[,,y + 1,1,,])
          sim_env$NAA[,,y + 1,1,,,sim] = sim_env$NAA[,,y + 1,1,,,sim] * fac
          sim_env$NAA0[,,y + 1,1,,,sim] = sim_env$NAA0[,,y + 1,1,,,sim] * fac
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
#' Runs the operating model's processes for year \code{y} and replicate
#' \code{sim} in order: initial age structure and first-year recruitment at
#' \code{y = 1}, population dynamics, fishery catches, indices and compositions,
#' survey indices and compositions, tag releases, fishery tag recaptures when any
#' \code{use_conv_fish_tagging = 1}, and next year's recruitment when
#' \code{y < n_yrs}.
#'
#' Those two standalone \code{generate_recruitment()} calls only run when
#' \code{rec_lag != 0}. Under \code{rec_lag = 0} recruitment depends on year
#' \code{y}'s own SSB, which is not known until \code{\link{apply_pop_dy}} reaches
#' \code{spawn_seas}, so it is called from inside \code{apply_pop_dy()} instead.
#'
#' @param y Integer. Year index.
#' @param sim Integer. Simulation replicate index.
#' @param sim_env Simulation environment from \code{\link{Setup_sim_env}}, passed
#'   by reference and modified in place by every annual-cycle helper.
#'
#' @return \code{invisible(NULL)}.
#'
#' @importFrom stats rnorm rmultinom
#' @export run_annual_cycle
#' @family Simulation Setup
run_annual_cycle <- function(y,
                             sim,
                             sim_env) {

  if(y == 1) {

    # note that some innovations are drawn before hand (in y = 1 here)
    if(isTRUE(sim_env$NAA_re > 0)) {

      sim_env$naa_eta <- draw_naa_innovations(sim_env) # get naa PE
      if(!is.null(sim_env$naa_eta_input)) { # if provided input eta
        n_cond <- dim(sim_env$naa_eta_input)[3]
        d <- dim(sim_env$naa_eta)
        sim_env$naa_eta[,,seq_len(n_cond),,,] <- array(sim_env$naa_eta_input[,,seq_len(n_cond),,,,sim], dim = c(d[1], d[2], n_cond, d[4], d[5], d[6]))
      }

      # if drawing own eta
      drawn <- sim_env$dsem_drawn$naa_eta_all
      if(!is.null(drawn) && any(drawn)) sim_env$naa_eta[drawn] <- array(sim_env$naa_eta_all[,,,,,,sim], dim = dim(drawn))[drawn]
      sim_env$naa_eta_all[,,,,,,sim] <- sim_env$naa_eta
    }

    draw_sim_q_devs(sim, sim_env)  # get catchability deviations
    generate_initial_age_structure(y = 1, sim, sim_env) # Initialize age structure
    if(sim_env$rec_lag != 0) generate_recruitment(y = 1, sim, sim_env) # Get recruitment in the first year
  }

  # growth kept cohort by cohort takes this year's start of year numbers, as the fit's population loop does
  if(!is.null(sim_env$growth_state) && y >= sim_env$dsem_growth_args$growth_cohort_styr) advance_sim_growth_year(y, sim, sim_env)

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

# Entry Point ---------------------------------------------------------------

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
                  srv_q = sim_env$srv_q, # note already scaled by catchability deviations
                  ln_fish_q_devs = sim_env$ln_fish_q_devs,
                  ln_srv_q_devs = sim_env$ln_srv_q_devs,
                  ln_RecDevs = sim_env$ln_RecDevs,
                  dsem_x_sim = sim_env$dsem_x_sim, # the dsem grid every replicate was drawn on, NULL without one
                  dsem_cov_obs_sim = sim_env$dsem_cov_obs_sim, # dsem covariate observations a refit reads
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
                  TrueDiscardAA = sim_env$TrueDiscardAA, ObsDiscardAA = sim_env$ObsDiscardAA,
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
