# Stage 2 of 3: objective function
#
# Spawning and total biomass derived from an abundance array. compute_biom_y
# serves the annual cycle inside the objective, derive_proj_biom serves the
# forward projection.

#' Compute Biomass
#'
# Computes spawning-time biomass quantities (Total_Biom, SSB, Dynamic_SSB0,
# eff_SSB) for year y using the current NAA/NAA0 state at season "seas"
# (always called with seas == spawn_seas). Takes all inputs explicitly
# (rather than relying on lexical scoping) since it's called from the main
# RTMB model function's local frame but defined at the top level of a
# different file.
#'
#'
#' @param y Year integer
#' @param seas Season integer
#' @keywords internal
compute_biom_y = function(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, seasdur,
                          n_seas, n_pop, n_regions, n_ages, n_sexes,
                          sgl_seas_spawning_movement, natal_region, stray_rate,
                          Movement = NULL, Mrate = NULL, move_timing = 0,
                          do_recruits_move = 1) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Get NAA for spawning
  tmp_NAA_spawn = NAA[,,y,seas,,, drop = FALSE]
  tmp_NAA0_spawn = NAA0[,,y,seas,,, drop = FALSE]

  # Propagate to the spawning point. Under move_timing == 0 movement has already been
  # applied at the top of the season, so the t_spawn mortality discount below is applied
  # as-is and this block is a no-op. Under move_timing 1 and 2 movement has not yet
  # happened, so spawners must be advanced through the appropriate fraction of the
  # season here, and the discount below is skipped (spawn_state already applies it).
  if(move_timing != 0 && n_regions > 1) {
    for(p in 1:n_pop) {
      for(a in 1:n_ages) {
        moves = (do_recruits_move == 1 || a > 1)
        Mv = if(moves) Movement[p,,,y,seas,a,1] else diag(n_regions)
        Qv = if(moves) Mrate[p,,,y,seas,a,1] else matrix(0, n_regions, n_regions)
        for(s in 1:n_sexes) {
          if(moves) { Mv = Movement[p,,,y,seas,a,s]; Qv = Mrate[p,,,y,seas,a,s] }
          tmp_NAA_spawn[p,,1,1,a,s] = spawn_state(tmp_NAA_spawn[p,,1,1,a,s], Mv,
                                                  ZAA[p,,y,seas,a,s], Qv, seasdur[seas], t_spawn, move_timing)
          tmp_NAA0_spawn[p,,1,1,a,s] = spawn_state(tmp_NAA0_spawn[p,,1,1,a,s], Mv,
                                                   natmort[p,,y,a,s] * seasdur[seas], Qv, seasdur[seas], t_spawn, move_timing)
        } # end s loop
      } # end a loop
    } # end p loop
  }

  # If we we are natal homing with 1 season
  if(n_seas == 1 && n_pop > 1) {
    # Get NAA during spawning
    for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
      tmp_NAA_spawn[p,,1,1,a,s] = tmp_NAA_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
      tmp_NAA0_spawn[p,,1,1,a,s] = tmp_NAA0_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
    } # end s loop
  }

  # Mortality discount up to spawning. Already folded into spawn_state above when
  # move_timing != 0, so it collapses to 1 in that case to avoid applying it twice.
  spawn_disc_Z = if(move_timing == 0 || n_regions == 1) exp(-ZAA[,,y,seas,,,drop = FALSE] * t_spawn) else 1
  spawn_disc_Z_f = if(move_timing == 0 || n_regions == 1) exp(-ZAA[,, y, seas, , 1,drop = FALSE] * t_spawn) else 1

  # Total Biomass
  Total_Biom_y = apply(tmp_NAA_spawn *
                         WAA[,, y, seas, , ,drop = FALSE] *
                         spawn_disc_Z, c(1,2), sum)

  # Spawning Stock Biomass
  SSB_y = apply(tmp_NAA_spawn[,, 1, 1, , 1,drop = FALSE] *
                  WAA[,, y, seas, , 1,drop = FALSE] *
                  MatAA[,, y, seas, , 1,drop = FALSE] *
                  spawn_disc_Z_f, c(1,2), sum)

  # Get dynamic B0
  SSB0_array = tmp_NAA0_spawn[,, 1, 1, , 1,drop = FALSE] *  WAA[,,  y, seas, , 1, drop = FALSE] * MatAA[,,y, seas, , 1, drop = FALSE]
  if(move_timing == 0 || n_regions == 1) {
    mort_spawn = exp(-natmort[,, y, , 1, drop = FALSE] * t_spawn * seasdur[seas])
    mort_spawn = array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
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
          eff_SSB_y[p2] = eff_SSB_y[p2] + (stray_rate[p,y] / n_receivers) * SSB_y[p, natal_region[p2]]
        }
      }
    }
  } else eff_SSB_y[1] = sum(SSB_y[1,])

  list(Total_Biom_y = Total_Biom_y, SSB_y = SSB_y, Dynamic_SSB0_y = Dynamic_SSB0_y, eff_SSB_y = eff_SSB_y)
}

#' Compute Biomass for Population Projections
#'
# Computes SSB / Total_Biom / Dynamic_SSB0 / eff_SSB for projection year y at season seas
# (always called with seas == spawn_seas) from the current proj_NAA/proj_NAA0
# state in Do_Population_Projection(). Factored out (plain R, no RTMB/AD
# concerns since Do_Population_Projection is never used inside an AD tape) so
# it can run either before or after that season's mortality/ageing step
# depending on whether rec_lag == 0 (age0_bh), without duplicating the math.
# Mirrors compute_biom_y() above, which serves the same role for the RTMB
# estimation model.
#'
#' @param y Projection year integer
#' @param seas Season integer (always spawn_seas)
#' @keywords internal
derive_proj_biom = function(y, seas, proj_NAA, proj_NAA0, WAA, MatAA, proj_ZAA, natmort, t_spawn, seasdur,
                            n_seas, n_pop, n_regions, n_ages, n_sexes,
                            sgl_seas_spawning_movement, natal_region, stray_rate,
                            Movement = NULL, Mrate = NULL, move_timing = 0,
                            do_recruits_move = 1) {

  tmp_NAA_spawn = proj_NAA[,,y,seas,,, drop = FALSE]
  tmp_NAA0_spawn = proj_NAA0[,,y,seas,,, drop = FALSE]

  # Propagate to the spawning point, as compute_biom_y() does for the estimation model.
  # Under move_timing 1 and 2 movement has not been applied yet at this point in the
  # season, so spawners are advanced here and the t_spawn discount below is dropped.
  if(move_timing != 0 && n_regions > 1) {
    for(p in 1:n_pop) {
      for(a in 1:n_ages) {
        moves = (do_recruits_move == 1 || a > 1)
        for(s in 1:n_sexes) {
          Mv = if(moves) Movement[p,,,y,seas,a,s] else diag(n_regions)
          Qv = if(moves) Mrate[p,,,y,seas,a,s] else matrix(0, n_regions, n_regions)
          tmp_NAA_spawn[p,,1,1,a,s] = spawn_state(tmp_NAA_spawn[p,,1,1,a,s], Mv,
                                                  proj_ZAA[p,,y,seas,a,s], Qv, seasdur[seas], t_spawn, move_timing)
          tmp_NAA0_spawn[p,,1,1,a,s] = spawn_state(tmp_NAA0_spawn[p,,1,1,a,s], Mv,
                                                   natmort[p,,y,a,s] * seasdur[seas], Qv, seasdur[seas], t_spawn, move_timing)
        } # end s loop
      } # end a loop
    } # end p loop
  }

  # If we we are natal homing with 1 season
  if(n_seas == 1 && n_pop > 1) {
    for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
      tmp_NAA_spawn[p,,1,1,a,s] = tmp_NAA_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
      tmp_NAA0_spawn[p,,1,1,a,s] = tmp_NAA0_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
    } # end s loop
  }

  # Mortality discount up to spawning; already folded into spawn_state when move_timing != 0
  spawn_disc_Z_f = if(move_timing == 0 || n_regions == 1) exp(-proj_ZAA[,, y, seas, , 1,drop = FALSE] * t_spawn) else 1
  spawn_disc_Z = if(move_timing == 0 || n_regions == 1) exp(-proj_ZAA[,, y, seas, , ,drop = FALSE] * t_spawn) else 1

  # get total biomass
  Total_Biom_y = apply(tmp_NAA_spawn *
                         WAA[,, y, seas, , ,drop = FALSE] *
                         spawn_disc_Z, c(1,2), sum)

  # get SSB
  SSB_y = apply(tmp_NAA_spawn[,, 1, 1, , 1,drop = FALSE] *
                  WAA[,, y, seas, , 1,drop = FALSE] *
                  MatAA[,, y, seas, , 1,drop = FALSE] *
                  spawn_disc_Z_f, c(1,2), sum)

  # Get dynamic B0
  SSB0_array = tmp_NAA0_spawn[,, 1, 1, , 1,drop = FALSE] *  WAA[,,  y, seas, , 1, drop = FALSE] * MatAA[,,y, seas, , 1, drop = FALSE]
  if(move_timing == 0 || n_regions == 1) {
    mort_spawn = exp(-natmort[,, y, , 1, drop = FALSE] * t_spawn * seasdur[seas])
    mort_spawn = array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
  } else mort_spawn = 1
  Dynamic_SSB0_y = apply(SSB0_array * mort_spawn, c(1,2), sum) # Dynamic B0

  if(n_sexes == 1) { # If single sex model, multiply SSB calculations by 0.5
    SSB_y = SSB_y * 0.5
    Dynamic_SSB0_y = Dynamic_SSB0_y * 0.5
  }

  # Accumulate effective SSB at each population's natal region across all source populations to capture stray contributions
  eff_SSB_y <- rep(0, n_pop)
  if(n_pop > 1) {
    n_pop_in_region <- rep(0, n_regions)
    for(p in 1:n_pop) n_pop_in_region[natal_region[p]] <- n_pop_in_region[natal_region[p]] + 1
    for(p2 in 1:n_pop) {
      for(p in 1:n_pop) {
        if(p == p2) {
          # Own population contribution - no stray scaling
          eff_SSB_y[p2] = eff_SSB_y[p2] + SSB_y[p, natal_region[p2]]
        } else {
          # Cross-population contribution scaled by stray_rate
          eff_SSB_y[p2] = eff_SSB_y[p2] + (stray_rate[p, y] / n_pop_in_region[natal_region[p2]]) * SSB_y[p, natal_region[p2]]
        }
      } # end p loop
    } # end p2 loop
  } else eff_SSB_y[1] = sum(SSB_y[1,])

  list(SSB_y = SSB_y, Total_Biom_y = Total_Biom_y, Dynamic_SSB0_y = Dynamic_SSB0_y, eff_SSB_y = eff_SSB_y)
}
