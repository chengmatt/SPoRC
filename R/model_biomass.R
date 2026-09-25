# Stage 2 of 3: objective function
#
# Calculate spawning and total biomass derived from an abundance array. biom_at_spawn does the calculation and
# the three wrappers slice their own state into it.

#' Spawning-time biomass for one year and season
#'
#' The calculation doing the annual cycle, the forward projection and the operating model.
#'
#' @param NAA_s Array \code{[pop, region, 1, 1, age, sex]} of abundance at the start of the season.
#' @param NAA0_s Array on the same dims of unfished abundance.
#' @param WAA_s Array on the same dims of weight at age.
#' @param MatAA_s Array on the same dims of maturity at age.
#' @param ZAA_s Array on the same dims of total mortality over the season.
#' @param natmort_s Array on the same dims of the natural mortality rate.
#' @param Movement_s Array \code{[pop, region, region, 1, 1, age, sex]} of movement fractions, or
#'   \code{NULL} when movement does not run at spawning.
#' @param Mrate_s Array on the same dims of the movement generator, or \code{NULL}.
#' @param spawn_move_s Array \code{[pop, region, region, 1, age, sex]} of natal homing movement,
#'   read only by a single-season multi-population model, \code{NULL} otherwise.
#' @param stray_rate_y Numeric \code{[pop]} of this year's stray rates, read only by a
#'   multi-population model, \code{NULL} otherwise.
#' @param t_spawn Fraction of the season elapsed at spawning.
#' @param seasdur_seas Duration of this season as a fraction of the year.
#' @param n_seas,n_pop,n_regions,n_ages,n_sexes Model dimensions.
#' @param natal_region Integer vector of each population's natal region.
#' @param move_timing \code{0} movement before mortality, \code{1} and \code{2} movement within the
#'   season, which folds the mortality discount into the propagation.
#' @param do_recruits_move Whether age one moves.
#' @param expm_nsub Substeps for the implicit matrix exponential.
#'
#' @return List of \code{Total_Biom_y}, \code{SSB_y}, \code{Dynamic_SSB0_y} and \code{eff_SSB_y}.
#'
#' @keywords internal
#' @import RTMB
biom_at_spawn = function(
  NAA_s,
  NAA0_s,
  WAA_s,
  MatAA_s,
  ZAA_s,
  natmort_s,
  spawn_move_s,
  stray_rate_y,
  t_spawn,
  seasdur_seas,
  n_seas,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  natal_region,
  Movement_s = NULL,
  Mrate_s = NULL,
  move_timing = 0,
  do_recruits_move = 1,
  expm_nsub = 0
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  tmp_NAA_spawn = NAA_s
  tmp_NAA0_spawn = NAA0_s

  # propagate to the spawning point. a no-op under move_timing == 0, where movement already ran;
  # under 1 and 2 spawners advance through part of the season and the discount below is skipped
  if(move_timing != 0 && n_regions > 1) {
    for(p in 1:n_pop) {
      for(a in 1:n_ages) {
        moves = (do_recruits_move == 1 || a > 1)
        for(s in 1:n_sexes) {
          Mv = if(moves) Movement_s[p,,,1,1,a,s] else diag(n_regions)
          Qv = if(moves) Mrate_s[p,,,1,1,a,s] else matrix(0, n_regions, n_regions)
          tmp_NAA_spawn[p,,1,1,a,s] = spawn_state(tmp_NAA_spawn[p,,1,1,a,s], Mv,
                                                  ZAA_s[p,,1,1,a,s], Qv, seasdur_seas, t_spawn, move_timing, expm_nsub = expm_nsub)
          tmp_NAA0_spawn[p,,1,1,a,s] = spawn_state(tmp_NAA0_spawn[p,,1,1,a,s], Mv,
                                                   natmort_s[p,,1,1,a,s] * seasdur_seas, Qv, seasdur_seas, t_spawn, move_timing, expm_nsub = expm_nsub)
        } # end s loop
      } # end a loop
    } # end p loop
  }

  # If we are natal homing with 1 season
  if(n_seas == 1 && n_pop > 1) {
    # Get NAA during spawning
    for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
      tmp_NAA_spawn[p,,1,1,a,s] = tmp_NAA_spawn[p,,1,1,a,s] %*% spawn_move_s[p,,,1,a,s]
      tmp_NAA0_spawn[p,,1,1,a,s] = tmp_NAA0_spawn[p,,1,1,a,s] %*% spawn_move_s[p,,,1,a,s]
    } # end s loop
  }

  # Mortality discount up to spawning. Already folded into spawn_state above when
  # move_timing != 0, so it collapses to 1 in that case to avoid applying it twice.
  spawn_disc_Z = if(move_timing == 0 || n_regions == 1) exp(-ZAA_s * t_spawn) else 1
  spawn_disc_Z_f = if(move_timing == 0 || n_regions == 1) exp(-ZAA_s[,, 1, 1, , 1,drop = FALSE] * t_spawn) else 1

  # Total Biomass
  Total_Biom_y = apply(tmp_NAA_spawn *
                         WAA_s *
                         spawn_disc_Z, c(1,2), sum)

  # Spawning Stock Biomass
  SSB_y = apply(tmp_NAA_spawn[,, 1, 1, , 1,drop = FALSE] *
                  WAA_s[,, 1, 1, , 1,drop = FALSE] *
                  MatAA_s[,, 1, 1, , 1,drop = FALSE] *
                  spawn_disc_Z_f, c(1,2), sum)

  # Get dynamic B0
  SSB0_array = tmp_NAA0_spawn[,, 1, 1, , 1,drop = FALSE] * WAA_s[,, 1, 1, , 1, drop = FALSE] * MatAA_s[,, 1, 1, , 1, drop = FALSE]
  if(move_timing == 0 || n_regions == 1) {
    mort_spawn = exp(-natmort_s[,, 1, 1, , 1, drop = FALSE] * t_spawn * seasdur_seas)
    mort_spawn = array(mort_spawn, dim = dim(SSB0_array)) # coerce array
  } else mort_spawn = 1
  Dynamic_SSB0_y = apply(SSB0_array * mort_spawn, c(1,2), sum) # Dynamic B0

  if(n_sexes == 1) { # If single sex model, multiply SSB calculations by 0.5
    SSB_y = SSB_y * 0.5
    Dynamic_SSB0_y = Dynamic_SSB0_y * 0.5
  }

  # Accumulate effective SSB at each population's natal region
  # across all source populations (captures stray contributions)
  eff_SSB_y = array(0, dim = n_pop)
  if(n_pop > 1) {

    # get number of pops in a given region
    n_pop_in_region = array(0, dim = n_regions)
    for(p in 1:n_pop) n_pop_in_region[natal_region[p]] = n_pop_in_region[natal_region[p]] + 1

    for(p2 in 1:n_pop) {
      for(p in 1:n_pop) {
        if(p == p2) {
          eff_SSB_y[p2] = eff_SSB_y[p2] + SSB_y[p, natal_region[p2]]
        } else {
          n_receivers = n_pop_in_region[natal_region[p2]]
          eff_SSB_y[p2] = eff_SSB_y[p2] + (stray_rate_y[p] / n_receivers) * SSB_y[p, natal_region[p2]]
        }
      }
    }
  } else eff_SSB_y[1] = sum(SSB_y[1,])

  list(
    Total_Biom_y = Total_Biom_y,
    SSB_y = SSB_y,
    Dynamic_SSB0_y = Dynamic_SSB0_y,
    eff_SSB_y = eff_SSB_y
  )
}

#' Compute Biomass
#'
#' Spawning-time biomass quantities for year \code{y} from the annual cycle's current state at
#' season \code{seas}, always the spawning season. 
#'
#' @param y Year integer
#' @param seas Season integer
#' @param NAA,NAA0 Abundance arrays \code{[pop, region, year, season, age, sex]}, fished and unfished.
#' @param WAA,MatAA,ZAA,natmort Weight, maturity, total mortality and the natural mortality rate on
#'   the same dims.
#' @param t_spawn Fraction of the season elapsed at spawning.
#' @param seasdur Numeric vector of season durations.
#' @param n_seas,n_pop,n_regions,n_ages,n_sexes Model dimensions.
#' @param sgl_seas_spawning_movement Natal homing movement \code{[pop, region, region, year, age, sex]}.
#' @param natal_region Integer vector of each population's natal region.
#' @param stray_rate Matrix \code{[pop, year]} of stray rates.
#' @param Movement,Mrate Movement fractions and generator
#'   \code{[pop, region, region, year, season, age, sex]}, or \code{NULL}.
#' @param move_timing,do_recruits_move,expm_nsub Movement timing, whether age one moves, and
#'   substeps for the implicit matrix exponential.
#'
#' @return List of \code{Total_Biom_y}, \code{SSB_y}, \code{Dynamic_SSB0_y} and \code{eff_SSB_y}.
#' @keywords internal
compute_biom_y = function(
  y,
  seas,
  NAA,
  NAA0,
  WAA,
  MatAA,
  ZAA,
  natmort,
  t_spawn,
  seasdur,
  n_seas,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  sgl_seas_spawning_movement,
  natal_region,
  stray_rate,
  Movement = NULL,
  Mrate = NULL,
  move_timing = 0,
  do_recruits_move = 1,
  expm_nsub = 0
) {

  propagates = move_timing != 0 && n_regions > 1 # if movement happening
  homes = n_seas == 1 && n_pop > 1 # if there is natal homing here

  biom_at_spawn(
    NAA_s = NAA[,,y,seas,,, drop = FALSE],
    NAA0_s = NAA0[,,y,seas,,, drop = FALSE],
    WAA_s = WAA[,,y,seas,,, drop = FALSE],
    MatAA_s = MatAA[,,y,seas,,, drop = FALSE],
    ZAA_s = ZAA[,,y,seas,,, drop = FALSE],
    natmort_s = natmort[,,y,seas,,, drop = FALSE],
    spawn_move_s = if(homes) sgl_seas_spawning_movement[,,,y,,, drop = FALSE] else NULL,
    stray_rate_y = if(n_pop > 1) stray_rate[,y] else NULL,
    t_spawn = t_spawn,
    seasdur_seas = seasdur[seas],
    n_seas = n_seas,
    n_pop = n_pop,
    n_regions = n_regions,
    n_ages = n_ages,
    n_sexes = n_sexes,
    natal_region = natal_region,
    Movement_s = if(propagates) Movement[,,,y,seas,,, drop = FALSE] else NULL,
    Mrate_s = if(propagates) Mrate[,,,y,seas,,, drop = FALSE] else NULL,
    move_timing = move_timing,
    do_recruits_move = do_recruits_move,
    expm_nsub = expm_nsub
  )
}

#' Compute Biomass for Population Projections
#'
#' Spawning-time biomass for projection year \code{y} at the spawning season, from the current
#' projected state. Plain R, so it can run either side of mortality depending on \code{rec_lag}.
#'
#' @param y Projection year integer
#' @param seas Season integer (always spawn_seas)
#' @param proj_NAA,proj_NAA0 Projected abundance arrays, fished and unfished.
#' @param proj_ZAA Projected total mortality.
#' @inheritParams compute_biom_y
#'
#' @return List of \code{Total_Biom_y}, \code{SSB_y}, \code{Dynamic_SSB0_y} and \code{eff_SSB_y}.
#' @keywords internal
derive_proj_biom = function(
  y,
  seas,
  proj_NAA,
  proj_NAA0,
  WAA,
  MatAA,
  proj_ZAA,
  natmort,
  t_spawn,
  seasdur,
  n_seas,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  sgl_seas_spawning_movement,
  natal_region,
  stray_rate,
  Movement = NULL,
  Mrate = NULL,
  move_timing = 0,
  do_recruits_move = 1,
  expm_nsub = 0
) {

  compute_biom_y(
    y = y,
    seas = seas,
    NAA = proj_NAA,
    NAA0 = proj_NAA0,
    WAA = WAA,
    MatAA = MatAA,
    ZAA = proj_ZAA,
    natmort = natmort,
    t_spawn = t_spawn,
    seasdur = seasdur,
    n_seas = n_seas,
    n_pop = n_pop,
    n_regions = n_regions,
    n_ages = n_ages,
    n_sexes = n_sexes,
    sgl_seas_spawning_movement = sgl_seas_spawning_movement,
    natal_region = natal_region,
    stray_rate = stray_rate,
    Movement = Movement,
    Mrate = Mrate,
    move_timing = move_timing,
    do_recruits_move = do_recruits_move,
    expm_nsub = expm_nsub
  )
}
