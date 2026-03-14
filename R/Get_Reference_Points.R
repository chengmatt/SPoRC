#' Get SPR reference points (Single Region)
#'
#' @param pars Parameter List
#' @param data Data List
#'
#' @keywords internal
#' @import RTMB
single_region_SPR <- function(pars,
                              data
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # Exponentiate reference point
  F_x = exp(log_F_x)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_ages)) # 2 slots in rows, for unfished, and fished at F_x

  # Set up the initial recruits
  for(p in 1:n_pop) Nspr[,p,1] = sex_ratio_f[p] * rec_seas_prop[p,1]

  ## Loop through ages
  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        tmp_unfished = Nspr[1,p, j - 1]
        tmp_fished   = Nspr[2,p, j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          tmp_unfished = tmp_unfished + rec_seas_prop[p,seas] * sex_ratio_f[p]
          tmp_fished   = tmp_fished   + rec_seas_prop[p,seas] * sex_ratio_f[p]
        }

        ## Spawning biomass
        if (seas == spawn_seas) {
          SB_age[1, p, j - 1] = tmp_unfished * WAA[p, spawn_seas, j - 1] * MatAA[p, spawn_seas, j - 1] * exp(-t_spawn * natmort[p, j - 1] * seasdur[seas])
          SB_age[2, p, j - 1] = tmp_fished * WAA[p, spawn_seas, j - 1] * MatAA[p, spawn_seas, j - 1] *
            exp(-t_spawn * (natmort[p, j - 1] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[j-1,]) ))
        }

        ## Mortality and ageing
        if (seas < n_seas) { # Within season mortality
          Nspr[1, p, j - 1] = tmp_unfished * exp(-natmort[p, j - 1] * seasdur[seas])
          Nspr[2, p, j - 1] = tmp_fished * exp(-(natmort[p, j - 1] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[j-1,]) ))
        } else {
          # Ageing
          Nspr[1, p,  j] = tmp_unfished * exp(-natmort[p, j - 1] * seasdur[seas])
          Nspr[2, p,  j] = tmp_fished * exp(-(natmort[p, j - 1] * seasdur[seas] +  sum(F_fract_flt[seas,] * F_x * fish_sel[j-1,]) ))
        }

      } # end seas loop
    } # end j loop
  } # end p loop

  # Age n_ages-1 is now at start of year after the loop
  tmp_unfished = Nspr[1,,n_ages-1]
  tmp_fished = Nspr[2,,n_ages-1]

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal mortality
        tmp_unfished[p] = tmp_unfished[p] * exp(-(natmort[p,n_ages-1] * seasdur[seas]))
        tmp_fished[p] = tmp_fished[p] * exp(-(natmort[p,n_ages-1] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[n_ages-1,]) ))

      } # end seas loop
    } # end p loop
  }

  SB_age[1, , n_ages - 1] = tmp_unfished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] * exp(-t_spawn * natmort[,n_ages - 1] * seasdur[spawn_seas])
  SB_age[2, , n_ages - 1] = tmp_fished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[,n_ages - 1] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[n_ages-1,]) ))


  # Plus group
  M_penult = natmort[,n_ages - 1]
  Z_penult = natmort[,n_ages - 1] + F_x * sum(colSums(F_fract_flt) * fish_sel[n_ages - 1,])

  Z_plus = natmort[,n_ages] + F_x * sum(colSums(F_fract_flt) * fish_sel[n_ages,])
  M_plus = natmort[,n_ages]

  Nspr[1,,n_ages] = Nspr[1,,n_ages-1] * exp(-M_penult) / (1 - exp(-M_plus))
  Nspr[2,,n_ages] = Nspr[2,,n_ages-1] * exp(-Z_penult) / (1 - exp(-Z_plus))

  tmp_unfished = Nspr[1,,n_ages]
  tmp_fished = Nspr[2,,n_ages]

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal mortality
        tmp_unfished[p] = tmp_unfished[p] * exp(-(natmort[p,n_ages] * seasdur[seas]))
        tmp_fished[p] = tmp_fished[p] * exp(-(natmort[p,n_ages] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[n_ages,]) ))

      } # end seas loop
    } # end p loop
  }

  ## Plus group spawning biomass
  SB_age[1,,n_ages] = tmp_unfished * WAA[,spawn_seas, n_ages] * MatAA[,spawn_seas, n_ages] * exp(-t_spawn * natmort[,n_ages] * seasdur[spawn_seas])
  SB_age[2,,n_ages] = tmp_fished * WAA[,spawn_seas, n_ages] * MatAA[,spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[,n_ages] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[n_ages,])))

  # Get effective SB after straying
  SB = apply(SB_age[2,,,drop = FALSE], 2, sum)
  SB0 = apply(SB_age[1,,,drop = FALSE], 2, sum)

  effective_SB = array(0, n_pop)
  effective_SB0 = array(0, n_pop)

  # compute effective SSB after straying
  for(p2 in 1:n_pop) {
    for(p in 1:n_pop) {
      if(p == p2) {
        # Own population contribution - no stray scaling
        effective_SB[p2] = effective_SB[p2] + SB[p]
        effective_SB0[p2]  = effective_SB0[p2]  + SB0[p]
      } else {
        effective_SB[p2] = effective_SB[p2] + stray_rate[p] * SB[p]
        effective_SB0[p2]  = effective_SB0[p2]  + stray_rate[p] * SB0[p]
      }
    }
  }

  # Get spawning biomass per recruit to get spawning potential ratio
  SPR = sum(effective_SB) / sum(effective_SB0)

  # compute objective function to get F_x
  sprpen = 100 * (SPR - SPR_x)^2

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(effective_SB0)
  RTMB::REPORT(effective_SB)
  RTMB::REPORT(F_x)
  RTMB::REPORT(SB)
  RTMB::REPORT(SB0)

  return(sprpen)
}


























































#' Title Get Global SPR Reference Points (Spatial)
#'
#' @param pars Parameter List from RTMB
#' @param data Data List from RTMB
#' @keywords internal
#' @import RTMB
global_SPR <- function(pars,
                       data
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # Exponentiate reference points
  F_x = exp(log_F_x)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy

  # Set up the initial recruits
  for(p in 1:n_pop) Nspr[2,p,,1] = Nspr[1,p,,1] = rec_region_prop[p,] * sex_ratio_f[p,] * rec_seas_prop[p,1]

  ## Loop through ages
  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        tmp_unfished = Nspr[1,p,, j - 1]
        tmp_fished   = Nspr[2,p,, j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          tmp_unfished = tmp_unfished + rec_seas_prop[p,seas] * sex_ratio_f[p,] * rec_region_prop[p,]
          tmp_fished   = tmp_fished  + rec_seas_prop[p,seas] * sex_ratio_f[p,] * rec_region_prop[p,]
        }

        ## Movement
        if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
          tmp_unfished = as.vector(tmp_unfished %*% Movement[p,,,seas, j - 1])
          tmp_fished   = as.vector(tmp_fished   %*% Movement[p,,,seas, j - 1])
        }

        ## Spawning biomass
        if (seas == spawn_seas) {

          # Extract temporary variables out
          tmp_unfished_spawn = tmp_unfished
          tmp_fished_spawn = tmp_fished

          # If single season natal homing population
          if(n_pop > 1 && n_seas == 1) {
            # Get NAA during spawning in single season case
            tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
            tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
          }

          SB_age[1,p,, j - 1] = tmp_unfished_spawn * WAA[p,, spawn_seas, j - 1] * MatAA[p,, spawn_seas, j - 1] * exp(-t_spawn * natmort[p,, j - 1] * seasdur[seas])
          SB_age[2,p,, j - 1] = tmp_fished_spawn * WAA[p,, spawn_seas, j - 1] * MatAA[p,, spawn_seas, j - 1] * exp(-t_spawn *
                                                                                                               (natmort[p,, j - 1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,j-1,,drop = F]) ))
        }

        ## Mortality and ageing
        if (seas < n_seas) { # Within season mortality
          Nspr[1,p,, j - 1] = tmp_unfished * exp(-natmort[p,, j - 1] * seasdur[seas])
          Nspr[2,p,, j - 1] = tmp_fished * exp(-(natmort[p,, j - 1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,j-1,,drop = F], 1, sum) ))
        } else {
          # Ageing
          Nspr[1,p,, j] = tmp_unfished * exp(-natmort[p,, j - 1] * seasdur[seas])
          Nspr[2,p,, j] = tmp_fished * exp(-(natmort[p,, j - 1] * seasdur[seas] +  apply(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,j-1,,drop = F], 1, sum) ))
        }
      }
    }
  }

  # Age n_ages-1 is now at start of year after the loop
  tmp_unfished = array(Nspr[1,,,n_ages-1], dim = c(n_pop, n_regions))
  tmp_fished = array(Nspr[2,,,n_ages-1], dim = c(n_pop, n_regions))

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal movement
        tmp_unfished[p,] = tmp_unfished[p,] %*% Movement[p,,,seas,n_ages-1]
        tmp_fished[p,] = tmp_fished[p,] %*% Movement[p,,,seas,n_ages-1]

        # Apply seasonal mortality
        tmp_unfished[p,] = tmp_unfished[p,] * exp(-(natmort[p,,n_ages-1] * seasdur[seas]))
        tmp_fished[p,] = tmp_fished[p,] * exp(-(natmort[p,,n_ages-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages - 1,,drop = F], 1, sum) ))

      } # end seas loop
    }
  }

  ## Penultimate age spawning biomass
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn = tmp_fished

  for(p in 1:n_pop) {
    tmp_unfished_spawn[p,] = as.vector(tmp_unfished_spawn[p,] %*% Movement[p,,, spawn_seas, n_ages - 1])
    tmp_fished_spawn[p,]   = as.vector(tmp_fished_spawn[p,] %*% Movement[p,,, spawn_seas, n_ages - 1])

    # If single season natal homing population
    if(n_pop > 1 && n_seas == 1) {
      # Get NAA during spawning in single season case
      tmp_unfished_spawn[p,] = as.vector(tmp_unfished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages-1])
      tmp_fished_spawn[p,] = as.vector(tmp_fished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages-1])
    }

  }


  SB_age[1,,, n_ages - 1] = tmp_unfished_spawn * WAA[,, spawn_seas, n_ages - 1] * MatAA[,, spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[,, n_ages - 1] * seasdur[spawn_seas])
  SB_age[2,,, n_ages - 1] = tmp_fished_spawn * WAA[,, spawn_seas, n_ages - 1] * MatAA[,, spawn_seas, n_ages - 1] * exp(-t_spawn *
                                                                                                                   (natmort[,, n_ages - 1] * seasdur[spawn_seas] + apply(F_fract_flt[,spawn_seas,,drop = F] * F_x * fish_sel[,n_ages - 1,,drop = F], 1, sum) ))

  ## Plus group analytical solution
  for(p in 1:n_pop) {

    T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

    for (seas in 1:n_seas) {
      # Get survival
      S_penult_unfished = diag(exp(-natmort[p,, n_ages - 1] * seasdur[seas]), n_regions)
      S_plus_unfished = diag(exp(-natmort[p,, n_ages] * seasdur[seas]), n_regions)
      S_penult_fished = diag(exp(-(natmort[p,, n_ages - 1] * seasdur[seas] +  apply(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages-1,,drop = F], 1, sum) )), n_regions)
      S_plus_fished = diag(exp(-(natmort[p,, n_ages] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages,,drop = F], 1, sum) )), n_regions)

      # Get transition matrices
      T_penult_unfished = S_penult_unfished %*% t(Movement[p,,,seas, n_ages - 1]) %*% T_penult_unfished
      T_plus_unfished = S_plus_unfished %*% t(Movement[p,,,seas, n_ages]) %*% T_plus_unfished
      T_penult_fished = S_penult_fished %*% t(Movement[p,,,seas, n_ages - 1]) %*% T_penult_fished
      T_plus_fished = S_plus_fished %*% t(Movement[p,,,seas, n_ages]) %*% T_plus_fished
    }

    source_unfished = T_penult_unfished %*% Nspr[1,p,, n_ages - 1]
    source_fished   = T_penult_fished %*% Nspr[2,p,, n_ages - 1]

    Nspr[1,p,, n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
    Nspr[2,p,, n_ages] = solve(diag(n_regions) - T_plus_fished, source_fished)
  }

  tmp_unfished = array(Nspr[1,,,n_ages], dim = c(n_pop, n_regions))
  tmp_fished = array(Nspr[2,,,n_ages], dim = c(n_pop, n_regions))

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal movement
        tmp_unfished[p,] = tmp_unfished[p,] %*% Movement[p,,,seas,n_ages]
        tmp_fished[p,] = tmp_fished[p,] %*% Movement[p,,,seas,n_ages]

        # Apply seasonal mortality
        tmp_unfished[p,] = tmp_unfished[p,] * exp(-(natmort[p,,n_ages] * seasdur[seas]))
        tmp_fished[p,] = tmp_fished[p,] * exp(-(natmort[p,,n_ages] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages,,drop = F], 1, sum) ))

      } # end seas loop
    }
  }

  # Extract temporary variables out
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn = tmp_fished

  ## Plus group spawning biomass
  for(p in 1:n_pop) {
    tmp_unfished_spawn[p,] = as.vector(tmp_unfished_spawn[p,] %*% Movement[p,,, spawn_seas, n_ages])
    tmp_fished_spawn[p,]   = as.vector(tmp_fished_spawn[p,] %*% Movement[p,,, spawn_seas, n_ages])

    # If single season natal homing population
    if(n_pop > 1 && n_seas == 1) {
      # Get NAA during spawning in single season case
      tmp_unfished_spawn[p,] = as.vector(tmp_unfished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages])
      tmp_fished_spawn[p,] = as.vector(tmp_fished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages])
    }
  }

  SB_age[1,,,n_ages] = tmp_unfished_spawn * WAA[,,spawn_seas,n_ages] * MatAA[,,spawn_seas,n_ages] * exp(-t_spawn * natmort[,,n_ages] * seasdur[spawn_seas])
  SB_age[2,,,n_ages] = tmp_fished_spawn  * WAA[,,spawn_seas,n_ages] * MatAA[,,spawn_seas,n_ages] * exp(-t_spawn * (natmort[,,n_ages] * seasdur[spawn_seas] + apply(F_fract_flt[,spawn_seas,,drop=F] * F_x * fish_sel[,n_ages,,drop=F], 1, sum)))

  # Get effective SB after straying
  SB = apply(SB_age[2,,,,drop = FALSE], 2:3, sum)
  SB0 = apply(SB_age[1,,,,drop = FALSE], 2:3, sum)
  effective_SB = array(0, dim = n_pop)
  effective_SB0 = array(0, dim = n_pop)

  # Accumulate effective SSB at each population's natal region
  # across all source populations (captures stray contributions)
  for(p2 in 1:n_pop) {
    for(p in 1:n_pop) {
      if(p == p2) {
        effective_SB[p2] = effective_SB[p2] + SB[p, natal_region[p2]]
        effective_SB0[p2] = effective_SB0[p2] + SB0[p, natal_region[p2]]
      } else {
        effective_SB[p2] = effective_SB[p2] + stray_rate[p] * SB[p, natal_region[p2]]
        effective_SB0[p2] = effective_SB0[p2] + stray_rate[p] * SB0[p, natal_region[p2]]
      }
    } # end p loop
  } # end p2 loop

  # Get spawning biomass per recruit to get spawning potential ratio
  if(n_pop == 1) {
    # Single pop - SPR is global sum across all regions
    SPR = sum(SB) / sum(SB0)
  } else {
    # Multi-pop natal homing - weighted by stray rates at natal regions
    SPR = sum(effective_SB) / sum(effective_SB0)
  }

  # compute objective function to get F_x
  sprpen = 100 * (SPR - SPR_x)^2

  RTMB::REPORT(SPR)
  RTMB::REPORT(effective_SB0)
  RTMB::REPORT(effective_SB)
  RTMB::REPORT(F_x)
  RTMB::REPORT(SB)
  RTMB::REPORT(SB0)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB_age)

  return(sprpen)
}

#'
#' Title Get singl region BH MSY Reference Points
#'
#' @param pars Parameter List from RTMB
#' @param data Data List from RTMB
#' @keywords internal
#' @import RTMB
single_region_BH_Fmsy <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # Exponentiate Fmsy
  Fmsy = exp(log_Fmsy)

  SB_age = Nspr = array(0, dim = c(2, n_pop, n_ages))
  CAA = array(0, c(n_pop, n_seas, n_ages))

  for(p in 1:n_pop) Nspr[,p,1] = sex_ratio_f[p] * rec_seas_prop[p,1]

  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        # Extract out numbers
        tmp_unfished = Nspr[1, p, j - 1]
        tmp_fished   = Nspr[2, p,  j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          tmp_unfished = tmp_unfished + rec_seas_prop[p,seas] * sex_ratio_f[p]
          tmp_fished   = tmp_fished   + rec_seas_prop[p,seas] * sex_ratio_f[p]
        }

        # Compute F and Z for this age/season
        F_a_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[j-1,])
        Z_a_seas = natmort[p, j - 1] * seasdur[seas] + F_a_seas

        # Spawning biomass
        if (seas == spawn_seas) {
          SB_age[1, p, j - 1] = tmp_unfished * WAA[p, spawn_seas, j - 1] * MatAA[p, spawn_seas, j - 1] *
            exp(-t_spawn * natmort[p, j - 1] * seasdur[seas])
          SB_age[2, p, j - 1] = tmp_fished * WAA[p, spawn_seas, j - 1] * MatAA[p, spawn_seas, j - 1] *
            exp(-t_spawn * (natmort[p, j - 1] * seasdur[seas] + sum(F_fract_flt[seas,] * Fmsy * fish_sel[j-1,])))
        }

        # Catch-at-age (Baranov)
        CAA[p, seas, j - 1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))

        # Mortality and ageing
        if (seas < n_seas) {
          Nspr[1, p, j - 1] = tmp_unfished * exp(-natmort[p, j - 1] * seasdur[seas])
          Nspr[2, p, j - 1] = tmp_fished * exp(-Z_a_seas)
        } else {
          Nspr[1, p, j] = tmp_unfished * exp(-natmort[p, j - 1] * seasdur[seas])
          Nspr[2, p, j] = tmp_fished * exp(-Z_a_seas)
        }

      } # end seas loop
    } # end j loop
  } # end p loop

  # Age n_ages-1 at start of year after loop
  tmp_unfished = Nspr[1,,n_ages - 1]
  tmp_fished = Nspr[2,, n_ages - 1]
  tmp_fished_caa = tmp_fished

  # Catch-at-age for penultimate age
  for (seas in 1:n_seas) {
    F_a_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages-1,])
    Z_a_seas = natmort[, n_ages-1] * seasdur[seas] + F_a_seas
    CAA[,seas, n_ages-1] = tmp_fished_caa * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
    tmp_fished_caa = tmp_fished_caa * exp(-Z_a_seas)  # advance
  }

  if (spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {
        # Exponential mortality
        tmp_unfished[p] = tmp_unfished[p] * exp(-natmort[p,n_ages - 1] * seasdur[seas])
        tmp_fished[p] = tmp_fished[p] * exp(-(sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages-1,]) + natmort[p,n_ages - 1] * seasdur[seas]))
      }
    }
  }

  # Get spawning biomass after mortality decrement
  SB_age[1,,n_ages - 1] = tmp_unfished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[,n_ages - 1] * seasdur[spawn_seas])
  SB_age[2,, n_ages - 1] = tmp_fished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[,n_ages - 1] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[n_ages-1,])))

  # Plus group
  Z_plus = natmort[,n_ages] + sum(colSums(F_fract_flt) * Fmsy * fish_sel[n_ages,])
  M_plus = natmort[,n_ages]
  Nspr[1, ,n_ages] = Nspr[1, ,n_ages - 1] * exp(-M_plus) / (1 - exp(-M_plus))
  Nspr[2, ,n_ages] = Nspr[2, ,n_ages - 1] * exp(-Z_plus) / (1 - exp(-Z_plus))

  # Plus group catch across seasons
  tmp_unfished = Nspr[1,, n_ages]
  tmp_fished = Nspr[2,, n_ages]
  tmp_fished_caa = tmp_fished

  # Catch-at-age for plus group
  for (seas in 1:n_seas) {
    F_a_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages,])
    Z_a_seas = natmort[,n_ages] * seasdur[seas] + F_a_seas
    CAA[, seas, n_ages] = tmp_fished_caa * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
    tmp_fished_caa = tmp_fished_caa * exp(-Z_a_seas)  # advance
  }

  if (spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {
        tmp_unfished[p] = tmp_unfished[p] * exp(-natmort[p,n_ages] * seasdur[seas])
        tmp_fished[p] = tmp_fished[p] * exp(-(natmort[p,n_ages] * seasdur[seas] + sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages,])))
      }
    }
  }

  SB_age[1,, n_ages] = tmp_unfished * WAA[,spawn_seas, n_ages] * MatAA[,spawn_seas, n_ages] *
    exp(-t_spawn * natmort[,n_ages] * seasdur[spawn_seas])
  SB_age[2,, n_ages] = tmp_fished * WAA[,spawn_seas, n_ages] * MatAA[,spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[,n_ages] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[n_ages,])))

  # Get effective SB after straying
  SB = apply(SB_age[2,,,drop = FALSE], 2, sum)
  SB0 = apply(SB_age[1,,,drop = FALSE], 2, sum)

  effective_SB = array(0, n_pop)
  effective_SB0 = array(0, n_pop)

  # compute effective SSB after straying
  for(p2 in 1:n_pop) {
    for(p in 1:n_pop) {
      if(p == p2) {
        # Own population contribution - no stray scaling
        effective_SB[p2] = effective_SB[p2] + SB[p]
        effective_SB0[p2]  = effective_SB0[p2]  + SB0[p]

      } else {
        effective_SB[p2] = effective_SB[p2] + stray_rate[p] * SB[p]
        effective_SB0[p2]  = effective_SB0[p2]  + stray_rate[p] * SB0[p]
      }
    }
  }

  Req = R0 * ((4 * h * effective_SB) - (1 - h) * effective_SB0) / ((5 * h - 1) * effective_SB)

  Yield = array(0, dim = c(n_pop))
  for(p in 1:n_pop) Yield[p] = sum(CAA[p,,] * WAA[p,,]) * Req[p]

  obj_fun = -sum(Yield)

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(effective_SB)
  RTMB::REPORT(effective_SB0)
  RTMB::REPORT(SB)
  RTMB::REPORT(SB0)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Yield)
  RTMB::REPORT(Req)
  RTMB::REPORT(CAA)

  return(obj_fun)
}

#' Title Get Global FMSY from a Beverton-Holt (Spatial)
#'
#' @param pars Parameter List
#' @param data Data List
#' @keywords internal
#' @import RTMB
global_BH_Fmsy <- function(pars,
                           data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # exponentitate reference points to "estimate"
  Fmsy = exp(log_Fmsy)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_seas, n_ages)) # catch at age

  # Set up the initial recruits
  Nspr[1,,1] = rec_region_prop * sex_ratio_f * rec_seas_prop[1]
  Nspr[2,,1] = rec_region_prop * sex_ratio_f * rec_seas_prop[1]

  ## Loop through ages
  for (j in 2:(n_ages - 1)) {
    for (seas in 1:n_seas) {

      tmp_unfished = Nspr[1,, j - 1]
      tmp_fished   = Nspr[2,, j - 1]

      # add in seasonal recruits
      if(seas > 1 && j - 1 == 1) {
        tmp_unfished = tmp_unfished + rec_seas_prop[seas] * sex_ratio_f
        tmp_fished   = tmp_fished   + rec_seas_prop[seas] * sex_ratio_f
      }

      ## Movement
      if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
        tmp_unfished = as.vector(tmp_unfished %*% Movement[,,seas, j - 1])
        tmp_fished   = as.vector(tmp_fished   %*% Movement[,,seas, j - 1])
      }

      ## Spawning biomass
      if (seas == spawn_seas) {
        SB_age[1,, j - 1] = tmp_unfished * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] * exp(-t_spawn * natmort[, j - 1] * seasdur[seas])
        SB_age[2,, j - 1] = tmp_fished * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] * exp(-t_spawn *
                                                                                                       (natmort[, j - 1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F]) ))
      }

      # Compute F and Z for this age/season
      F_a_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F])
      Z_a_seas = natmort[,j - 1] * seasdur[seas] + F_a_seas

      # Catch-at-age (Baranov)
      CAA[,seas, j - 1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))

      ## Mortality and ageing
      if (seas < n_seas) { # Within season mortality
        Nspr[1,, j - 1] = tmp_unfished * exp(-natmort[, j - 1] * seasdur[seas])
        Nspr[2,, j - 1] =
          tmp_fished * exp(-(natmort[, j - 1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F]) ))
      } else {
        # Ageing
        Nspr[1,, j] = tmp_unfished * exp(-natmort[, j - 1] * seasdur[seas])
        Nspr[2,, j] = tmp_fished * exp(-(natmort[, j - 1] * seasdur[seas] +  rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F]) ))
      }
    }
  }

  # Age n_ages-1 is now at start of year after the loop
  tmp_unfished = Nspr[1,,n_ages-1]
  tmp_fished = Nspr[2,,n_ages-1]

  # Catch-at-age for penultimate age
  for (seas in 1:n_seas) {
    # Compute F and Z for this age/season
    F_a_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F])
    Z_a_seas = natmort[,n_ages - 1] * seasdur[seas] + F_a_seas

    # Catch-at-age (Baranov)
    CAA[,seas, n_ages - 1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
  }

  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # Apply seasonal movement
      tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages-1]
      tmp_fished = tmp_fished %*% Movement[,,seas,n_ages-1]

      # Apply seasonal mortality
      tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages-1] * seasdur[seas]))
      tmp_fished = tmp_fished * exp(-(natmort[,n_ages-1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages - 1,,drop = F]) ))

    } # end seas loop
  }

  ## Penultimate age spawning biomass
  tmp_unfished = as.vector(tmp_unfished %*% Movement[,, spawn_seas, n_ages - 1])
  tmp_fished   = as.vector(tmp_fished %*% Movement[,, spawn_seas, n_ages - 1])
  SB_age[1,, n_ages - 1] = tmp_unfished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[, n_ages - 1] * seasdur[spawn_seas])
  SB_age[2,, n_ages - 1] = tmp_fished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] * exp(-t_spawn *
                                                                                                                (natmort[, n_ages - 1] * seasdur[spawn_seas] + rowSums(F_fract_flt[,spawn_seas,,drop = F] * Fmsy * fish_sel[,n_ages - 1,,drop = F]) ))

  ## Plus group analytical solution
  T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

  for (seas in 1:n_seas) {
    # Get survival
    S_penult_unfished = diag(exp(-natmort[, n_ages - 1] * seasdur[seas]), n_regions)
    S_plus_unfished = diag(exp(-natmort[, n_ages] * seasdur[seas]), n_regions)
    S_penult_fished = diag(exp(-(natmort[, n_ages - 1] * seasdur[seas] +  rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F]))), n_regions)
    S_plus_fished = diag(exp(-(natmort[, n_ages] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F]) )), n_regions)

    # Get transition matrices
    T_penult_unfished = S_penult_unfished %*% t(Movement[,,seas, n_ages - 1]) %*% T_penult_unfished
    T_plus_unfished = S_plus_unfished %*% t(Movement[,,seas, n_ages]) %*% T_plus_unfished
    T_penult_fished = S_penult_fished %*% t(Movement[,,seas, n_ages - 1]) %*% T_penult_fished
    T_plus_fished = S_plus_fished %*% t(Movement[,,seas, n_ages]) %*% T_plus_fished
  }

  source_unfished = T_penult_unfished %*% Nspr[1,, n_ages - 1]
  source_fished   = T_penult_fished %*% Nspr[2,, n_ages - 1]

  Nspr[1,, n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
  Nspr[2,, n_ages] = solve(diag(n_regions) - T_plus_fished, source_fished)

  tmp_unfished = Nspr[1,,n_ages]
  tmp_fished = Nspr[2,,n_ages]

  # Catch-at-age for penultimate age
  for (seas in 1:n_seas) {
    # Compute F and Z for this age/season
    F_a_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F])
    Z_a_seas = natmort[,n_ages] * seasdur[seas] + F_a_seas

    # Catch-at-age (Baranov)
    CAA[,seas, n_ages] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
  }


  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # Apply seasonal movement
      tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages]
      tmp_fished = tmp_fished %*% Movement[,,seas,n_ages]

      # Apply seasonal mortality
      tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages] * seasdur[seas]))
      tmp_fished = tmp_fished * exp(-(natmort[,n_ages] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F]) ))

    } # end seas loop
  }

  ## Plus group spawning biomass
  tmp_unfished = as.vector(tmp_unfished %*% Movement[,, spawn_seas, n_ages])
  tmp_fished   = as.vector(tmp_fished %*% Movement[,, spawn_seas, n_ages])
  SB_age[1,, n_ages] = tmp_unfished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] * exp(-t_spawn * natmort[, n_ages] * seasdur[spawn_seas])
  SB_age[2,, n_ages] = tmp_fished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[, n_ages] * seasdur[spawn_seas] + rowSums(F_fract_flt[,spawn_seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F]) ))

  # Get spawning biomass per recruit to get spawning potential ratio
  SBPR_0 = sum(SB_age[1,,])
  SBPR_F = sum(SB_age[2,,])
  SPR = SBPR_F / SBPR_0

  # Get equilibrium recruitment
  Req = R0 * ( (4*h*SBPR_F) - (1 - h) * SBPR_0) / ((5 * h -1) * SBPR_F)

  # Get yield
  Yield = sum(CAA * WAA) * Req
  Yield_r = rowSums(CAA * WAA) * Req

  # Get Bmsy
  # Bmsy = SBPR_F * Req
  # B0 = SBPR_0 * R0

  # compute objective function to get Fmsy
  obj_fun = -Yield

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SBPR_0)
  RTMB::REPORT(SBPR_F)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Yield)
  RTMB::REPORT(Yield_r)
  # RTMB::REPORT(Bmsy)
  # RTMB::REPORT(B0)
  RTMB::REPORT(Req)
  RTMB::REPORT(SPR)

  return(obj_fun)
}

#' Title Get Local FMSY from a Beverton-Holt (Spatial)
#'
#' @param pars Parameter List
#' @param data Data List
#' @keywords internal
#' @import RTMB
local_BH_Fmsy_sglpop <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c( n_regions, n_regions, n_seas, n_ages)) # catch at age
  Yield_r = array(0, dim = n_regions) # yield by region
  SB_unfished_mat = matrix(0, n_regions, n_regions)  # unfished spawning biomass per recruit
  SB_fished_mat = matrix(0,  n_regions, n_regions) # fished spawning biomass per recruit

  # exponentitate reference points to "solve"
  Fmsy = exp(log_Fmsy)

  # Set up the initial recruits (1 recruit per area)
  for(o in 1:n_regions) {
    for(d in 1:n_regions) {
      if(o == d) Nspr[1,o,d,1] = Nspr[2,o,d,1] = 1 * sex_ratio_f[o]
      else Nspr[1,o,d,1] = Nspr[2,o,d,1] = 0
    } # end d loop
  } # end o loop

  # Loop through, apply movement first, then decrement recruit
  # Loop through ages, projecting each cohort through the full annual cycle
  for(j in 2:(n_ages-1)){

    # Project age j-1 through all seasons to become age j
    for(seas in 1:n_seas) {

      for(o in 1:n_regions) {
        # Get temporary values from origin region
        tmp_unfished = Nspr[1,o,,j-1]
        tmp_fished = Nspr[2,o,,j-1]

        F_a_seas = apply(F_fract_flt[,seas,,drop=F] * Fmsy * fish_sel[,j-1,,drop=F], 1, sum)
        Z_a_seas = natmort[,j - 1] * seasdur[seas] + F_a_seas

        # Apply movement
        if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
          tmp_unfished = tmp_unfished %*% Movement[,,seas,j-1]
          tmp_fished = tmp_fished %*% Movement[,,seas,j-1]
        }

        # Calculate spawning biomass if this is the spawning season
        if(seas == spawn_seas) {
          for(d in 1:n_regions) {
            SB_age[1,o,d,j-1] = tmp_unfished[d] * WAA[d,spawn_seas,j-1] * MatAA[d,spawn_seas,j-1] *
              exp(-(t_spawn * natmort[d,j-1] * seasdur[spawn_seas]))
            SB_age[2,o,d,j-1] =
              tmp_fished[d] * WAA[d,spawn_seas,j-1] * MatAA[d,spawn_seas,j-1] *
              exp(-t_spawn * (natmort[d,j-1] * seasdur[spawn_seas] + F_a_seas[d]))
          }
        }

        # Compute F and Z for this age/season
        CAA[o,,seas,j-1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))

        # Apply mortality
        if(seas < n_seas) {
          # Within-season mortality, no ageing yet
          Nspr[1,o,,j-1] = tmp_unfished * exp(-(natmort[,j-1] * seasdur[seas]))
          Nspr[2,o,,j-1] = tmp_fished * exp(-(natmort[,j-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F], 1, sum) ))
        } else {
          # Last season: mortality + ageing
          Nspr[1,o,,j] = tmp_unfished * exp(-(natmort[,j-1] * seasdur[seas]))
          Nspr[2,o,,j] = tmp_fished * exp(-(natmort[,j-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F], 1, sum) ))
        }

      } # end o loop
    } # end seas loop
  } # end j loop

  # Now calculate spawning biomass for penultimate age (n_ages-1)
  for(o in 1:n_regions) {
    # Age n_ages-1 is now at start of year after the loop
    tmp_unfished = Nspr[1,o,,n_ages-1]
    tmp_fished = Nspr[2,o,,n_ages-1]

    for(seas in 1:n_seas) {
      # Compute F and Z for this age/season
      F_a_seas = apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F], 1, sum)
      Z_a_seas = natmort[,n_ages - 1] * seasdur[seas] + F_a_seas
      CAA[o,,seas,n_ages-1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
    }

    if(spawn_seas > 1) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal movement
        tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages-1]
        tmp_fished = tmp_fished %*% Movement[,,seas,n_ages-1]

        # Apply seasonal mortality
        tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages-1] * seasdur[seas]))
        tmp_fished = tmp_fished * exp(-(natmort[,n_ages-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F], 1, sum) ))

      } # end seas loop
    }

    # Apply movement for spawning season
    tmp_unfished = tmp_unfished %*% Movement[,,spawn_seas,n_ages-1]
    tmp_fished = tmp_fished %*% Movement[,,spawn_seas,n_ages-1]

    # Calculate spawning biomass
    for(d in 1:n_regions) {
      SB_age[1,o,d,n_ages-1] = tmp_unfished[d] * WAA[d,spawn_seas,n_ages-1] * MatAA[d,spawn_seas,n_ages-1] *
        exp(-(t_spawn * natmort[d,n_ages-1] * seasdur[spawn_seas]))
      SB_age[2,o,d,n_ages-1] = tmp_fished[d] * WAA[d,spawn_seas,n_ages-1] * MatAA[d,spawn_seas,n_ages-1] *
        exp(-t_spawn * ((natmort[d,n_ages-1] * seasdur[spawn_seas]) + sum(F_fract_flt[d,spawn_seas,] * Fmsy[d] * fish_sel[d,n_ages-1,])))
    }
  }

  # Set up analytical solution for plus group
  for(o in 1:n_regions) {

    # Build FULL annual transition for penultimate and plus ages
    T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

    # Loop through ALL seasons to build annual transition matrix
    for(seas in 1:n_seas) {
      # Unfished
      S_penult_unfished = diag(exp(-(natmort[,n_ages-1] * seasdur[seas])), n_regions)
      S_plus_unfished = diag(exp(-(natmort[,n_ages] * seasdur[seas])), n_regions)
      T_penult_unfished = S_penult_unfished %*% t(Movement[,,seas,n_ages-1]) %*% T_penult_unfished
      T_plus_unfished = S_plus_unfished %*% t(Movement[,,seas,n_ages]) %*% T_plus_unfished

      # Fished
      S_penult_fished = diag(exp(-(natmort[,n_ages-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F], 1, sum))), n_regions)
      S_plus_fished = diag(exp(-(natmort[,n_ages] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F], 1, sum))), n_regions)
      T_penult_fished = S_penult_fished %*% t(Movement[,,seas,n_ages-1]) %*% T_penult_fished
      T_plus_fished = S_plus_fished %*% t(Movement[,,seas,n_ages]) %*% T_plus_fished
    } # end seas loop

    # Solve for equilibrium plus group (at start of year)
    source_unfished = T_penult_unfished %*% Nspr[1,o,,n_ages-1]
    Nspr[1,o,,n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
    source_fished = T_penult_fished %*% Nspr[2,o,,n_ages-1]
    Nspr[2,o,,n_ages] = solve(diag(n_regions) - T_plus_fished, source_fished)

  } # end o loop

  # Now calculate spawning biomass for penultimate age (n_ages-1)
  for(o in 1:n_regions) {
    # Age n_ages-1 is now at start of year after the loop
    tmp_unfished = Nspr[1,o,,n_ages]
    tmp_fished = Nspr[2,o,,n_ages]

    for(seas in 1:n_seas) {
      # Compute F and Z for this age/season
      F_a_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F])
      Z_a_seas = natmort[,n_ages ] * seasdur[seas] + F_a_seas
      CAA[o,,seas,n_ages] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
    }

    if(spawn_seas > 1) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal movement
        tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages]
        tmp_fished = tmp_fished %*% Movement[,,seas,n_ages]

        # Apply seasonal mortality
        tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages] * seasdur[seas]))
        tmp_fished = tmp_fished * exp(-(natmort[,n_ages] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F], 1, sum)))

      } # end seas loop
    }

    # Apply movement for spawning season
    tmp_unfished = tmp_unfished %*% Movement[,,spawn_seas,n_ages]
    tmp_fished = tmp_fished %*% Movement[,,spawn_seas,n_ages]

    # Calculate spawning biomass
    for(d in 1:n_regions) {
      SB_age[1,o,d,n_ages] = tmp_unfished[d] * WAA[d,spawn_seas,n_ages] * MatAA[d,spawn_seas,n_ages] *
        exp(-(t_spawn * natmort[d,n_ages] * seasdur[spawn_seas]))
      SB_age[2,o,d,n_ages] = tmp_fished[d] * WAA[d,spawn_seas,n_ages] * MatAA[d,spawn_seas,n_ages] *
        exp(-t_spawn * ((natmort[d,n_ages] * seasdur[spawn_seas]) + sum(F_fract_flt[d,spawn_seas,] * Fmsy[d] * fish_sel[d,n_ages,])))
    } # end d loop

  } # end o loop

  # Determine equilibrium recruitment for destination region
  # parse out and compute unfished and fished spawning biomass per recruit
  for(o in 1:n_regions) {
    for(d in 1:n_regions) {
      SB_unfished_mat[o, d] = sum(SB_age[1, o, d, ])  # unfished
      SB_fished_mat[o, d] = sum(SB_age[2, o, d, ])  # fished at Fmsy
    } # end o loop
  } # end d loop

  A = 4 * h * rec_region_prop * R0 # define first part of the numerator of BH recruitment
  B = rep(0, n_regions) # define first part of the denominator of BH recruitment
  for(d in 1:n_regions) B[d] = (1 - h[d]) * sum(SB_unfished_mat[,d] * rec_region_prop * R0)
  C = 5 * h - 1 # define second part of the denominator for BH recruitment

  # define initial guess to solve for equilibrium recruitment from origin region
  Req_o = R0 * rec_region_prop

  for(nit in 1:newton_steps) {
    # compute equilibrium spawning biomass (SSBR * Req) in destination region
    x_vec = as.numeric(t(SB_fished_mat) %*% Req_o)  # function of equilibrium recruitment in origin region (what we are solving for)
    numer_vec = A * x_vec # compute numerator of BH
    denom_vec = B + (C * x_vec) # compute denominator of BH
    g_vec = numer_vec / denom_vec # equilibrium recruitment in destination region

    # define root and deffine Jacobian
    iter_vec = Req_o - g_vec # find values of origin recruitment that are consisitent w/ destination recruitment such that pop'n is in equilibrium

    # construct Jacobian for root
    # we need J = df (iter_vec)/dReq = dReq/dReq (or I) - dg/dReq
    # we basically want to know dg / dReq (how does destination equil rec change as orign equil rec change)
    # dg / dReq = (dg / dxk) * (dxk / dReq)
    # to get (dg / dxk), use quotient rule of (BH recruitment)
    # note that dxk / dReq S_2mat * Req = S_2mat
    dg_dxk = (A * B) / (denom_vec^2)
    dg_dReq = matrix(0, n_regions, n_regions)
    for(d in 1:n_regions) dg_dReq[d, ] = dg_dxk[d] * SB_fished_mat[, d]# now compute to see how destination equilibrium rec changes, as origin equil rec changes

    # compute jacobian
    J = diag(1, n_regions) - dg_dReq
    delta = solve(J, iter_vec) # get step to move towards solution
    Req_o = Req_o - delta # newton raphson update
  }

  # get destination reigon yield
  for(d in 1:n_regions) {
    tmp = 0 # define temp variable
    for(seas in 1:n_seas) for(o in 1:n_regions) tmp = tmp + sum(CAA[o, d, seas, ] * WAA[d, seas, ]) * Req_o[o] # get yield to destination
    Yield_r[d] = tmp
  } # end d loop

  SB  <- matrix(0, n_pop, n_regions)
  SB0 <- matrix(0, n_pop, n_regions)
  for(p in 1:n_pop) for(r in 1:n_regions) {
    SB[p, r]  <- SB_fished_mat[p, r]    # already weighted by rec_region_prop
    SB0[p, r] <- SB_unfished_mat[p, r]
  }


  # maximize total yield
  Yield_total = sum(Yield_r)
  obj_fun = -Yield_total

  sum_SB_unfished_mat = sum(SB_unfished_mat)

  # RTMB::REPORT(eqrec_prop)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Req_o)
  RTMB::REPORT(Yield_r)
  RTMB::REPORT(Yield_total)
  RTMB::REPORT(dg_dReq)
  RTMB::REPORT(iter_vec)
  RTMB::REPORT(SB_fished_mat)
  RTMB::REPORT(SB_unfished_mat)
  RTMB::REPORT(sum_SB_unfished_mat)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB_age)
  RTMB::REPORT(SB)
  RTMB::REPORT(SB0)

  return(obj_fun)
}

#' Title Get Local FMSY from a Beverton-Holt (Spatial)
#'
#' @param pars Parameter List
#' @param data Data List
#' @keywords internal
#' @import RTMB
local_BH_Fmsy_multipop <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_regions, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_pop, n_regions, n_regions, n_seas, n_ages)) # catch at age
  Yield_r = array(0, dim = n_regions) # yield by region
  SB_unfished_mat = array(0, dim = c(n_pop, n_regions, n_regions))  # unfished spawning biomass per recruit
  SB_fished_mat = array(0, dim = c(n_pop, n_regions, n_regions)) # fished spawning biomass per recruit

  # exponentitate reference points to "solve"
  Fmsy = exp(log_Fmsy)

  # Set up the initial recruits (1 recruit per natal region / area)
  for(p in 1:n_pop) {
    for(o in 1:n_regions) {
      for(d in 1:n_regions) {

        if(o == d) Nspr[1,p,o,d,1] = Nspr[2,p,o,d,1] = sex_ratio_f[p,o] * rec_seas_prop[p,1] * rec_region_prop[p,o]
        else Nspr[1,p,o,d,1] = Nspr[2,p,o,d,1] = 0

      } # end d loop
    } # end o loop
  } # end p loop

  # Loop through, apply movement first, then decrement recruit
  # Loop through ages, projecting each cohort through the full annual cycle
  for(j in 2:(n_ages-1)){

    # Project age j-1 through all seasons to become age j
    for(seas in 1:n_seas) {

      for(p in 1:n_pop) {
        for(o in 1:n_regions) {

          F_a_seas = apply(F_fract_flt[,seas,,drop=F] * Fmsy * fish_sel[,j-1,,drop=F], 1, sum)
          Z_a_seas = natmort[p,,j - 1] * seasdur[seas] + F_a_seas

          # Get temporary values from origin region
          tmp_unfished = Nspr[1,p,o,,j-1]
          tmp_fished = Nspr[2,p,o,,j-1]

          # add in seasonal recruits
          if(seas > 1 && j - 1 == 1) {
            tmp_unfished[o] = tmp_unfished[o] + rec_seas_prop[p,seas] * sex_ratio_f[p,o] * rec_region_prop[p,o]
            tmp_fished[o]   = tmp_fished[o]   + rec_seas_prop[p,seas] * sex_ratio_f[p,o] * rec_region_prop[p,o]
          }

          # Apply movement
          if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
            tmp_unfished = tmp_unfished %*% Movement[p,,,seas,j-1]
            tmp_fished = tmp_fished %*% Movement[p,,,seas,j-1]
          }

          # Calculate spawning biomass if this is the spawning season
          if(seas == spawn_seas) {

            # Extract temporary variables out
            tmp_unfished_spawn = tmp_unfished
            tmp_fished_spawn = tmp_fished

            # If single season natal homing population
            if(n_pop > 1 && n_seas == 1) {
              # Get NAA during spawning in single season case
              tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
              tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
            }

            for(d in 1:n_regions) {
              SB_age[1,p,o,d,j-1] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1] *
                exp(-(t_spawn * natmort[p,d,j-1] * seasdur[spawn_seas]))
              SB_age[2,p,o,d,j-1] =
                tmp_fished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1] *
                exp(-t_spawn * (natmort[p,d,j-1] * seasdur[spawn_seas] + F_a_seas[d]))
            }
          }

          # Compute F and Z for this age/season
          CAA[p,o,,seas,j-1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))

          # Apply mortality
          if(seas < n_seas) {
            # Within-season mortality, no ageing yet
            Nspr[1,p,o,,j-1] = tmp_unfished * exp(-(natmort[p,,j-1] * seasdur[seas]))
            Nspr[2,p,o,,j-1] = tmp_fished * exp(-(natmort[p,,j-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F], 1, sum) ))
          } else {
            # Last season: mortality + ageing
            Nspr[1,p,o,,j] = tmp_unfished * exp(-(natmort[p,,j-1] * seasdur[seas]))
            Nspr[2,p,o,,j] = tmp_fished * exp(-(natmort[p,,j-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F], 1, sum) ))
          }

        } # end o loop
      } # end p loop

    } # end seas loop
  } # end j loop

  # Now calculate spawning biomass for penultimate age (n_ages-1)
  for(p in 1:n_pop) {
    for(o in 1:n_regions) {

      # Age n_ages-1 is now at start of year after the loop
      tmp_unfished = Nspr[1,p,o,,n_ages-1]
      tmp_fished = Nspr[2,p,o,,n_ages-1]

      for(seas in 1:n_seas) {
        # Compute F and Z for this age/season
        F_a_seas = apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F], 1, sum)
        Z_a_seas = natmort[p,,n_ages - 1] * seasdur[seas] + F_a_seas
        CAA[p,o,,seas,n_ages-1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
      }

      if(spawn_seas > 1) {
        for (seas in 1:(spawn_seas - 1)) {

          # Apply seasonal movement
          tmp_unfished = tmp_unfished %*% Movement[p,,,seas,n_ages-1]
          tmp_fished = tmp_fished %*% Movement[p,,,seas,n_ages-1]

          # Apply seasonal mortality
          tmp_unfished = tmp_unfished * exp(-(natmort[p,,n_ages-1] * seasdur[seas]))
          tmp_fished = tmp_fished * exp(-(natmort[p,,n_ages-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F], 1, sum) ))

        } # end seas loop
      }

      # Extract temporary variables out
      tmp_unfished_spawn = tmp_unfished
      tmp_fished_spawn = tmp_fished

      tmp_unfished_spawn = as.vector(tmp_unfished_spawn %*% Movement[p,,,spawn_seas,n_ages-1])
      tmp_fished_spawn   = as.vector(tmp_fished_spawn   %*% Movement[p,,,spawn_seas,n_ages-1])

      # If single season natal homing population
      if(n_pop > 1 && n_seas == 1) {
        # Get NAA during spawning in single season case
        tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages-1]
        tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages-1]
      }

      F_spawn = apply(F_fract_flt[,spawn_seas,,drop=F] * Fmsy * fish_sel[,n_ages-1,,drop=F], 1, sum)
      for(d in 1:n_regions) {
        SB_age[1,p,o,d,n_ages-1] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,n_ages-1] * MatAA[p,d,spawn_seas,n_ages-1] *
          exp(-t_spawn * natmort[p,d,n_ages-1] * seasdur[spawn_seas])
        SB_age[2,p,o,d,n_ages-1] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,n_ages-1] * MatAA[p,d,spawn_seas,n_ages-1] *
          exp(-t_spawn * (natmort[p,d,n_ages-1] * seasdur[spawn_seas] + F_spawn[d]))
      } # end

    } # end o
  } # end p

  # Set up analytical solution for plus group
  for(p in 1:n_pop) {

    T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

    for(seas in 1:n_seas) {
      S_penult_unfished = diag(exp(-(natmort[p,,n_ages-1] * seasdur[seas])), n_regions)
      S_plus_unfished   = diag(exp(-(natmort[p,,n_ages]   * seasdur[seas])), n_regions)
      S_penult_fished   = diag(exp(-(natmort[p,,n_ages-1] * seasdur[seas] + apply(F_fract_flt[,seas,,drop=F] * Fmsy * fish_sel[,n_ages-1,,drop=F], 1, sum))), n_regions)
      S_plus_fished     = diag(exp(-(natmort[p,,n_ages]   * seasdur[seas] + apply(F_fract_flt[,seas,,drop=F] * Fmsy * fish_sel[,n_ages,,drop=F],   1, sum))), n_regions)
      T_penult_unfished = S_penult_unfished %*% t(Movement[p,,,seas,n_ages-1]) %*% T_penult_unfished
      T_plus_unfished   = S_plus_unfished   %*% t(Movement[p,,,seas,n_ages])   %*% T_plus_unfished
      T_penult_fished   = S_penult_fished   %*% t(Movement[p,,,seas,n_ages-1]) %*% T_penult_fished
      T_plus_fished     = S_plus_fished     %*% t(Movement[p,,,seas,n_ages])   %*% T_plus_fished
    }

    for(o in 1:n_regions) {
      source_unfished = T_penult_unfished %*% Nspr[1,p,o,,n_ages-1]
      Nspr[1,p,o,,n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
      source_fished   = T_penult_fished   %*% Nspr[2,p,o,,n_ages-1]
      Nspr[2,p,o,,n_ages] = solve(diag(n_regions) - T_plus_fished,   source_fished)
    }

  } # end p loop

  # Now calculate spawning biomass for penultimate age (n_ages-1)
  for(p in 1:n_pop) {
    for(o in 1:n_regions) {

      # Age n_ages-1 is now at start of year after the loop
      tmp_unfished = Nspr[1,p,o,,n_ages]
      tmp_fished = Nspr[2,p,o,,n_ages]

      for(seas in 1:n_seas) {
        # Compute F and Z for this age/season
        F_a_seas = apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F], 1, sum)
        Z_a_seas = natmort[p,,n_ages ] * seasdur[seas] + F_a_seas
        CAA[p,o,,seas,n_ages] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
      }

      if(spawn_seas > 1) {
        for (seas in 1:(spawn_seas - 1)) {

          # Apply seasonal movement
          tmp_unfished = tmp_unfished %*% Movement[p,,,seas,n_ages]
          tmp_fished = tmp_fished %*% Movement[p,,,seas,n_ages]

          # Apply seasonal mortality
          tmp_unfished = tmp_unfished * exp(-(natmort[p,,n_ages] * seasdur[seas]))
          tmp_fished = tmp_fished * exp(-(natmort[p,,n_ages] * seasdur[seas] + apply(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F], 1, sum)))

        } # end seas loop
      }

      # Extract temporary variables out
      tmp_unfished_spawn = tmp_unfished
      tmp_fished_spawn = tmp_fished

      tmp_unfished_spawn = as.vector(tmp_unfished_spawn %*% Movement[p,,,spawn_seas,n_ages])
      tmp_fished_spawn   = as.vector(tmp_fished_spawn   %*% Movement[p,,,spawn_seas,n_ages])

      # If single season natal homing population
      if(n_pop > 1 && n_seas == 1) {
        # Get NAA during spawning in single season case
        tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages]
        tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages]
      }

      # Calculate spawning biomass
      F_spawn = apply(F_fract_flt[,spawn_seas,,drop=F] * Fmsy * fish_sel[,n_ages,,drop=F], 1, sum)
      for(d in 1:n_regions) {
        SB_age[1,p,o,d,n_ages] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,n_ages] * MatAA[p,d,spawn_seas,n_ages] *
          exp(-(t_spawn * natmort[p,d,n_ages] * seasdur[spawn_seas]))
        SB_age[2,p,o,d,n_ages] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,n_ages] * MatAA[p,d,spawn_seas,n_ages] *
          exp(-t_spawn * ((natmort[p,d,n_ages] * seasdur[spawn_seas]) + F_spawn[d]))
      } # end d loop

    } # end o loop
  } # end p loop

  # Determine equilibrium recruitment for destination region
  # parse out and compute unfished and fished spawning biomass per recruit
  for(p in 1:n_pop) {
    for(o in 1:n_regions) {
      for(d in 1:n_regions) {
        SB_unfished_mat[p, o, d] = sum(SB_age[1, p, o, d, ])  # unfished
        SB_fished_mat[p, o, d] = sum(SB_age[2, p, o, d, ])  # fished at Fmsy
      } # end o loop
    } # end d loop
  }

  SBPR_fished   <- matrix(0, n_pop, n_regions)  # [p, d]
  SBPR_unfished <- matrix(0, n_pop, n_regions)  # [p, d]
  for(p in 1:n_pop) for(d in 1:n_regions) {
    SBPR_fished[p, d]   <- sum(SB_fished_mat[p,,d])
    SBPR_unfished[p, d] <- sum(SB_unfished_mat[p,,d])
  }

  # Virgin effective SSB at each pop's natal region
  eff_SSB0_virgin = rep(0, n_pop)
  for(p2 in 1:n_pop) {
    r = natal_region[p2]
    for(p in 1:n_pop) {
      sc = if(p == p2) 1 else stray_rate[p]
      eff_SSB0_virgin[p2] = eff_SSB0_virgin[p2] + sc * SBPR_unfished[p, r] * R0[p]
    }
  }

  A_p <- 4 * h[cbind(1:n_pop, natal_region)] * R0   # define first part of the numerator of BH recruitment
  B_p <- (1 - h[cbind(1:n_pop, natal_region)]) * eff_SSB0_virgin # define first part of the denominator of BH recruitment
  C_p <- 5 * h[cbind(1:n_pop, natal_region)] - 1 # define second part of the denominator for BH recruitment

  # define initial guess to solve for equilibrium recruitment from origin region
  Req_o <- R0

  for(nit in 1:newton_steps) {

    # Effective fished SSB at each pop's natal region
    eff_SSB_fished <- rep(0, n_pop)
    for(p2 in 1:n_pop) {
      r <- natal_region[p2]
      for(p in 1:n_pop) {
        sc <- if(p == p2) 1 else stray_rate[p]
        eff_SSB_fished[p2] <- eff_SSB_fished[p2] + sc * SBPR_fished[p, r] * Req_o[p]
      }
    }

    # BH equilibrium total recruitment per pop
    g_vec <- A_p * eff_SSB_fished / (B_p + C_p * eff_SSB_fished)

    # define root and deffine Jacobian
    iter_vec = Req_o - g_vec # find values of origin recruitment that are consisitent w/ destination recruitment such that pop'n is in equilibrium

    # Jacobian
    dg_deff <- (A_p * B_p) / (B_p + C_p * eff_SSB_fished)^2
    J <- diag(n_pop)
    for(p2 in 1:n_pop) {
      r <- natal_region[p2]
      for(p in 1:n_pop) {
        sc <- if(p == p2) 1 else stray_rate[p]
        J[p2, p] <- J[p2, p] - dg_deff[p2] * sc * SBPR_fished[p, r]
      }
    }
    Req_o <- Req_o - solve(J, iter_vec)

  }


  # get destination reigon yield
  for(d in 1:n_regions) {
    tmp <- 0
    for(p in 1:n_pop) for(seas in 1:n_seas) for(o in 1:n_regions) {
      tmp <- tmp + sum(CAA[p,o,d,seas,] * WAA[p,d,seas,]) * rec_region_prop[p,o] * Req_o[p]
    }
    Yield_r[d] <- tmp
  }

  SB  <- matrix(0, n_pop, n_regions)
  SB0 <- matrix(0, n_pop, n_regions)
  for(p in 1:n_pop) for(r in 1:n_regions) {
    SB[p, r]  <- SBPR_fished[p, r]    # already weighted by rec_region_prop
    SB0[p, r] <- SBPR_unfished[p, r]
  }


  # maximize total yield
  Yield_total = sum(Yield_r)
  obj_fun = -Yield_total

  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Req_o)
  RTMB::REPORT(Yield_r)
  RTMB::REPORT(Yield_total)
  RTMB::REPORT(iter_vec)
  RTMB::REPORT(SB_fished_mat)
  RTMB::REPORT(SB_unfished_mat)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB_age)
  RTMB::REPORT(SB)
  RTMB::REPORT(SB0)

  return(obj_fun)
}

#' Wrapper function to get reference points
#'
#' Wrapper function to compute fishing and biological reference points given data and report
#' objects from an assessment or simulation. Supports both single-region and multi-region
#' calculations with options for SPR or Beverton–Holt MSY reference points.
#'
#' @param data List. Data object containing agesears, weight-at-age, maturity, natural mortality, and other simulation/assessment info.
#' @param rep List. Report object from RTMB containing estimated parameters like Fmort, selectivity, recruitment, steepness.
#' @param SPR_x Numeric. Target Spawning Potential Ratio fraction. Required for SPR-based reference points.
#' @param t_spawn Numeric. Mortality time until spawning.
#' @param sex_ratio_f Numeric vector. Female sex ratio by region.
#' @param calc_rec_st_yr Integer. First year used to compute mean recruitment.
#' @param rec_age Integer. Age at recruitment.
#' @param type Character. "single_region" or "multi_region".
#' @param what Character. Type of reference point:
#'   \describe{
#'     \item{SPR}{Single-region SPR reference point}
#'     \item{independent_SPR}{Multi-region SPR without movement}
#'     \item{global_SPR}{Multi-region SPR with movement}
#'     \item{BH_MSY}{Single-region Beverton–Holt MSY}
#'     \item{independent_BH_MSY}{Multi-region BH-MSY without movement}
#'     \item{global_BH_MSY}{Multi-region global BH-MSY with movement}
#'     \item{local_BH_MSY}{Multi-region local BH-MSY with movement}
#'   }
#' @param n_avg_yrs Integer. Number of years to average demographic rates when calculating reference points.
#' @param local_bh_msy_newton_steps Number of newton steps to take to solve for equilibrium recruitment in the origin region
#' when local_BH_MSY is assumed.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{f_ref_pt}{Vector of fishing reference points for each region.}
#'     \item{b_ref_pt}{Vector of biological reference points for each region.}
#'     \item{virgin_b_ref_pt}{Vector of virgin biomass reference points for each region.}
#'   }
#'
#' @importFrom stats nlminb
#' @import RTMB
#' @export Get_Reference_Points
#' @family Reference Points and Projections
Get_Reference_Points <- function(data,
                                 rep,
                                 SPR_x = NULL,
                                 t_spawn = 0,
                                 sex_ratio_f = array(0.5, dim = c(data$n_pop, data$n_regions)),
                                 calc_rec_st_yr = 1,
                                 rec_age = 1,
                                 type,
                                 what,
                                 n_avg_yrs = 1,
                                 local_bh_msy_newton_steps = 6
                                 ) {

  # Dimensions
  n_years <- length(data$years) # number of years
  n_ages <- length(data$ages) # number of ages
  n_pop <- data$n_pop # number of populations
  n_seas <- data$n_seas # number of populations
  n_regions <- data$n_regions # number of regions
  n_fish_fleets <- data$n_fish_fleets # number of fleets

  f_ref_pt <- vector() # set up storage
  virgin_pop_b_ref_pt <- pop_b_ref_pt <- virgin_b_ref_pt <- b_ref_pt <- array(0, dim = c(n_pop, n_regions)) # set up storage

  # determine years to average over demogrphaics
  n_yrs <- length(data$years)
  avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs

  if(type == "single_region") {

    data_list <- list() # set up data list

    # Seasonal stuff
    data_list$t_spawn <- t_spawn # specified mortality time up until spawning
    data_list$n_seas <- data$n_seas # number of seasons
    data_list$seasdur <- data$seasdur # seasonal duration
    data_list$spawn_seas <- data$spawn_seas # spawning season

    # Dimensions
    data_list$n_pop <- n_pop
    data_list$n_ages <- n_ages

    if(!what %in% c("SPR", "BH_MSY")) stop("what is not correctly specified! Should be SPR, BH_MSY for type = single_region")

    # fishing mortality fraction
    data_list$F_fract_flt <- array(rep$Fmort[1,n_years,,] / sum(rep$Fmort[1,n_years,,]),
                                   dim = c(data$n_seas, data$n_fish_fleets)) # get fleet F fraction to derive population level selectivity

    # fishery selectivity
    fish_sel_avg <- apply(rep$fish_sel[1,avg_yrs,,1,,drop = FALSE], c(3,5), mean)
    data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

    # natural mortality
    natmort_avg <- apply(rep$natmort[,1,avg_yrs,,1,drop = FALSE], c(1,4), mean)
    data_list$natmort <- natmort_avg # get female natural mortality

    # weight at age
    WAA_avg <- apply(data$WAA[,1,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
    data_list$WAA <- array(WAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # weight at age for females

    # maturity at age
    MatAA_avg <- apply(data$MatAA[,1,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
    data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # maturity at age for females

    data_list$rec_seas_prop <- array(rep$rec_seas_prop[,], dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
    data_list$sex_ratio_f <- sex_ratio_f # recritment sex ratio
    data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate

    if(what == 'SPR') {
      data_list$SPR_x <- SPR_x # SPR fraction
      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # Make adfun object
      tmp_obj <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
      tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt[1] <- tmp_obj$rep$F_x

      # Compute population specific reference points, by using stray rates
      mean_rec <- apply(rep$Rec[,1,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], 1, mean)
      for(p2 in 1:n_pop) {
        pop_b_ref_pt[p2,1]        <- tmp_obj$rep$SB[p2]  * mean_rec[p2]
        virgin_pop_b_ref_pt[p2,1] <- tmp_obj$rep$SB0[p2] * mean_rec[p2]
        for(p in 1:n_pop) {
          if(p != p2) {
            pop_b_ref_pt[p2,1]        <- pop_b_ref_pt[p2,1]        + stray_rate[p] * tmp_obj$rep$SB[p]  * mean_rec[p]
            virgin_pop_b_ref_pt[p2,1] <- virgin_pop_b_ref_pt[p2,1] + stray_rate[p] * tmp_obj$rep$SB0[p] * mean_rec[p]
          } # end if
        } # end p loop
      } # end p2 loop

      # Compute global reference points (sum across populations)
      b_ref_pt[,1] <- tmp_obj$rep$SB * apply(rep$Rec[,1,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)
      virgin_b_ref_pt[,1] <- tmp_obj$rep$SB0 * apply(rep$Rec[,1,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)

    } # end SPR reference points

    if(what == 'BH_MSY') {

      # extract out beverton-holt parameters
      data_list$h <- array(rep$h_trans[,1], dim = n_pop) # steepness
      data_list$R0 <- array(rep$R0, dim = n_pop) # unfished recruitment

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # make adfun ect
      tmp_obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
      tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt[1] <- tmp_obj$rep$Fmsy

      for(p2 in 1:n_pop) {
        pop_b_ref_pt[p2,1]        <- tmp_obj$rep$SB[p2]  * tmp_obj$rep$Req[p2]
        virgin_pop_b_ref_pt[p2,1] <- tmp_obj$rep$SB0[p2] * rep$R0[p2]
        for(p in 1:n_pop) {
          if(p != p2) {
            pop_b_ref_pt[p2,1]        <- pop_b_ref_pt[p2,1]        + stray_rate[p] * tmp_obj$rep$SB[p]  * tmp_obj$rep$Req[p]
            virgin_pop_b_ref_pt[p2,1] <- virgin_pop_b_ref_pt[p2,1] + stray_rate[p] * tmp_obj$rep$SB0[p] * data_list$R0[p]
          }
        }
      }

      b_ref_pt[,1] <- tmp_obj$rep$SB * tmp_obj$rep$Req
      virgin_b_ref_pt[,1] <- tmp_obj$rep$SB0 * data_list$R0

    }
  }

  if(type == 'multi_region') {

    if(!what %in% c("independent_SPR", "independent_BH_MSY", "global_SPR", "global_BH_MSY", "local_BH_MSY"))
      stop("what is not correctly specified! Should be independent_SPR, independent_BH_MSY, global_SPR, global_BH_MSY, local_BH_MSY for type = multi_region")

    data_list <- list() # set up data list

    # Seasonal stuff
    data_list$t_spawn <- t_spawn # specified mortality time up until spawning
    data_list$n_seas <- data$n_seas # number of seasons
    data_list$seasdur <- data$seasdur # seasonal duration
    data_list$spawn_seas <- data$spawn_seas # spawning season

    # Dimensions
    data_list$n_pop <- n_pop
    data_list$n_ages <- n_ages
    data_list$n_regions <- n_regions

    if(what == "independent_SPR") {

      for(r in 1:data$n_regions) {

        data_list$F_fract_flt <- array(rep$Fmort[r,n_years,,] / sum(rep$Fmort[r,n_years,,]),
                                       dim = c(data$n_seas, data$n_fish_fleets)) # get fleet F fraction to derive population level selectivity

        # fishery selectivity
        fish_sel_avg <- apply(rep$fish_sel[r,avg_yrs,,1,,drop = FALSE], c(1, 3, 4, 5), mean)
        data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

        # natural mortality
        natmort_avg <- apply(rep$natmort[,r,avg_yrs,,1,drop = FALSE], c(1,4), mean)
        data_list$natmort <- array(natmort_avg, dim = c(n_pop, n_ages)) # get female natural mortality

        # weight at age
        WAA_avg <- apply(data$WAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$WAA <- array(WAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # weight at age for females

        # maturity at age
        MatAA_avg <- apply(data$MatAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # maturity at age for females

        data_list$SPR_x <- SPR_x # SPR fraction
        data_list$sex_ratio_f <- array(sex_ratio_f[,r], dim = n_pop) # recritment sex ratio
        data_list$rec_seas_prop <- array(rep$rec_seas_prop[,], dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
        data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate

        par_list <- list() # set up parameter list
        par_list$log_F_x <- log(0.1) # F_x starting value

        # Make adfun object
        tmp_obj <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
        tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
        tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

        # Output reference points
        f_ref_pt[r] <- tmp_obj$rep$F_x

        # Compute population specific reference points, by using stray rates
        mean_rec <- apply(rep$Rec[,r,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], 1, mean)
        if(n_pop > 1) {
          for(p2 in 1:n_pop) {
            pop_b_ref_pt[p2,r]        <- tmp_obj$rep$SB[p2]  * mean_rec[p2]
            virgin_pop_b_ref_pt[p2,r] <- tmp_obj$rep$SB0[p2] * mean_rec[p2]
            for(p in 1:n_pop) {
              if(p != p2) {
                pop_b_ref_pt[p2,r]        <- pop_b_ref_pt[p2,r]        + stray_rate[p] * tmp_obj$rep$SB[p]  * mean_rec[p]
                virgin_pop_b_ref_pt[p2,r] <- virgin_pop_b_ref_pt[p2,r] + stray_rate[p] * tmp_obj$rep$SB0[p] * mean_rec[p]
              } # end if
            } # end p loop
          } # end p2 loop
        }

        # Compute global reference points (sum across populations)
        b_ref_pt[,r] <- tmp_obj$rep$SB * apply(rep$Rec[,r,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)
        virgin_b_ref_pt[,r] <- tmp_obj$rep$SB0 * apply(rep$Rec[,r,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)

      } # end r loop

      if(n_pop == 1) {
        pop_b_ref_pt[1,1] = sum(b_ref_pt)
        virgin_pop_b_ref_pt[1,1] = sum(virgin_b_ref_pt)
      }


    } # end independent_SPR

    if(what == "independent_BH_MSY") {

      for(r in 1:data$n_regions) {

        data_list$F_fract_flt <- array(rep$Fmort[r,n_years,,] / sum(rep$Fmort[r,n_years,,]),
                                       dim = c(data$n_seas, data$n_fish_fleets)) # get fleet F fraction to derive population level selectivity

        # fishery selectivity
        fish_sel_avg <- apply(rep$fish_sel[r,avg_yrs,,1,,drop = FALSE], c(1, 3, 4, 5), mean)
        data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

        # natural mortality
        natmort_avg <- apply(rep$natmort[,r,avg_yrs,,1,drop = FALSE], c(1,4), mean)
        data_list$natmort <- array(natmort_avg, dim = c(n_pop, n_ages)) # get female natural mortality

        # weight at age
        WAA_avg <- apply(data$WAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$WAA <- array(WAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # weight at age for females

        # maturity at age
        MatAA_avg <- apply(data$MatAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # maturity at age for females

        # Beverton Holt parameters
        data_list$h <- array(rep$h_trans[,r], dim = n_pop) # steepness
        data_list$R0 <- array(rep$R0 * rep$rec_region_prop[,r], dim = n_pop) # unfished recruitment by region
        data_list$sex_ratio_f <- array(sex_ratio_f[,r], dim = n_pop) # recritment sex ratio
        data_list$rec_seas_prop <- array(rep$rec_seas_prop, dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
        data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate

        par_list <- list() # set up parameter list
        par_list$log_Fmsy <- log(0.1) # Fmsy starting value

        # Make adfun object
        tmp_obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
        tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
        tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

        # Output reference points
        f_ref_pt[r] <- tmp_obj$rep$Fmsy
        if(n_pop > 1) {
          for(p2 in 1:n_pop) {
            pop_b_ref_pt[p2,r]        <- tmp_obj$rep$SB[p2]  * tmp_obj$rep$Req[p2]
            virgin_pop_b_ref_pt[p2,r] <- tmp_obj$rep$SB0[p2] * rep$R0[p2]
            for(p in 1:n_pop) {
              if(p != p2) {
                pop_b_ref_pt[p2,r]        <- pop_b_ref_pt[p2,r]        + stray_rate[p] * tmp_obj$rep$SB[p]  * tmp_obj$rep$Req[p]
                virgin_pop_b_ref_pt[p2,r] <- virgin_pop_b_ref_pt[p2,r] + stray_rate[p] * tmp_obj$rep$SB0[p] * data_list$R0[p]
              }
            }
          }
        }

        b_ref_pt[,r] <-  tmp_obj$rep$SB * tmp_obj$rep$Req
        virgin_b_ref_pt[,r] <-  tmp_obj$rep$SB0 * data_list$R0
      } # end r loop

      if(n_pop == 1) {
        pop_b_ref_pt[1,1] = sum(b_ref_pt)
        virgin_pop_b_ref_pt[1,1] = sum(virgin_b_ref_pt)
      }


    } # end independent_SPR

    if(what == 'global_SPR') {

      # Fleet fraction F
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio

      # fishery selectivity
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

      # natural mortality
      natmort_avg <- apply(rep$natmort[,,avg_yrs,,1,drop = FALSE], c(1,2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_pop, n_regions, n_ages)) # get female natural mortality

      # weight at age
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_pop, n_regions, data$n_seas, n_ages)) # weight at age for females

      # maturity at age
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, n_regions, data$n_seas, n_ages)) # maturity at age for females

      # Movement
      Movement_avg <- apply(rep$Movement[,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)) # Movement

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$sex_ratio_f <- sex_ratio_f # recritment sex ratio
      data_list$rec_seas_prop <- array(rep$rec_seas_prop[,], dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
      data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate
      data_list$SPR_x <- SPR_x # SPR fraction
      sgl_seas_spawning_movement_avg <- apply(rep$sgl_seas_spawning_movement[,,,avg_yrs,,1,drop = FALSE], c(1,2,3,5), mean)
      data_list$sgl_seas_spawning_movement <- array(sgl_seas_spawning_movement_avg, dim = c(n_pop, n_regions, n_regions, n_ages)) # Movement
      data_list$natal_region <- data$natal_region

      mean_rec <- apply(rep$Rec[,,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], c(1,2), mean) # [n_pop, n_regions]
      total_mean_rec <- apply(mean_rec, 1, sum) # [n_pop] - total recruitment across regions
      data_list$rec_region_prop <- mean_rec / total_mean_rec # recruitment proportions

      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # make adfn object
      tmp_obj <- RTMB::MakeADFun(cmb(global_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
      tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

      # output reference points
      f_ref_pt <- rep(tmp_obj$rep$F_x, n_regions)

      # Region-specific physical SSB
      for(r in 1:n_regions) {
        b_ref_pt[,r]        <- tmp_obj$rep$SB[,r]  * total_mean_rec
        virgin_b_ref_pt[,r] <- tmp_obj$rep$SB0[,r] * total_mean_rec
      }

      # Population-specific effective SSB at natal region
      if(n_pop > 1) {
        for(p2 in 1:n_pop) {
          r_natal <- natal_region[p2]
          pop_b_ref_pt[p2, r_natal]        <- tmp_obj$rep$SB[p2, r_natal]  * total_mean_rec[p2]
          virgin_pop_b_ref_pt[p2, r_natal] <- tmp_obj$rep$SB0[p2, r_natal] * total_mean_rec[p2]
          for(p in 1:n_pop) {
            if(p != p2) {
              pop_b_ref_pt[p2, r_natal]        <- pop_b_ref_pt[p2, r_natal]        + stray_rate[p] * tmp_obj$rep$SB[p, r_natal]  * total_mean_rec[p]
              virgin_pop_b_ref_pt[p2, r_natal] <- virgin_pop_b_ref_pt[p2, r_natal] + stray_rate[p] * tmp_obj$rep$SB0[p, r_natal] * total_mean_rec[p]
            }
          }
        }
      } else {
        pop_b_ref_pt[1,1] = sum(b_ref_pt)
        virgin_pop_b_ref_pt[1,1] = sum(virgin_b_ref_pt)
      }

    } # end global SPR

    if(what == 'global_BH_MSY') {

      # Error out if invalid recruitment density dependent option
      if(n_pop > 1) stop("Invalid reference point option! When n_pop > 1 reference points must either be independent_SPR, independent_BH_MSY, global_SPR, or local_BH_MSY.")

      # Fleet fraction F
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio

      # fishery selectivity
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

      # natural mortality
      natmort_avg <- apply(rep$natmort[1,,avg_yrs,,1,drop = FALSE], c(2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_regions, n_ages)) # get female natural mortality

      # weight at age
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # weight at age for females

      # maturity at age
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # maturity at age for females

      # Movement
      Movement_avg <- apply(rep$Movement[1,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_regions, n_regions, n_seas, n_ages)) # Movement

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$rec_region_prop <- rep$rec_region_prop[1,] # recruitment proportions
      data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio to use
      data_list$rec_seas_prop <- rep$rec_seas_prop[1,] # seasonal recruitment
      data_list$h <- mean(rep$h_trans[1,]) # steepness
      data_list$R0 <- rep$R0[1]  # unfished recruitment

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # Make adfun object
      tmp_obj <- RTMB::MakeADFun(cmb(global_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
      tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt <- rep(tmp_obj$rep$Fmsy, n_regions)
      b_ref_pt[1,] <- apply(tmp_obj$rep$SB_age[2,,,drop = F], 2, sum) * tmp_obj$rep$Req
      pop_b_ref_pt[1,1] <- sum(b_ref_pt)
      virgin_b_ref_pt[1,] <- apply(tmp_obj$rep$SB_age[1,,,drop = F], 2, sum) * data_list$R0
      virgin_pop_b_ref_pt[1,1]  <- sum(virgin_b_ref_pt)
    }

    if(what == 'local_BH_MSY') {

      # Fleet fraction F
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio

      # fishery selectivity
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

      # natural mortality
      natmort_avg <- apply(rep$natmort[,,avg_yrs,,1,drop = FALSE], c(1,2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, n_ages)) # get female natural mortality

      # weight at age
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, data$n_seas, n_ages)) # weight at age for females

      # maturity at age
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, data$n_seas, n_ages)) # maturity at age for females

      # Movement
      Movement_avg <- apply(rep$Movement[,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean)
      data_list$Movement <- array(Movement_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, n_regions, n_seas, n_ages)) # Movement

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$rec_region_prop <- array(rep$rec_region_prop, dim = c(if(n_pop > 1) n_pop else NULL, n_regions)) # recruitment proportions
      data_list$sex_ratio_f <- array(sex_ratio_f, dim = c(if(n_pop > 1) n_pop else NULL, n_regions)) # recruitment sex ratio to use
      data_list$rec_seas_prop <- array(rep$rec_seas_prop, dim = c(if(n_pop > 1) n_pop else NULL, data$n_seas)) # seasonal recruitment
      data_list$h <- array(rep$h_trans, dim = c(if(n_pop > 1) n_pop else NULL, n_regions)) # steepness
      data_list$R0 <- rep$R0  # unfished recruitment
      data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = data$n_pop) # stray rate
      data_list$newton_steps <- local_bh_msy_newton_steps # number of newton steps to take
      data_list$natal_region <- data$natal_region
      sgl_seas_spawning_movement_avg <- apply(rep$sgl_seas_spawning_movement[,,,avg_yrs,,1,drop = FALSE], c(1,2,3,5), mean)
      data_list$sgl_seas_spawning_movement <- array(sgl_seas_spawning_movement_avg, dim = c(n_pop, n_regions, n_regions, n_ages)) # Movement


      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- rep(log(0.1), n_regions) # Fmsy starting value

      # Make adfun object
      if(n_pop == 1) tmp_obj <- RTMB::MakeADFun(cmb(local_BH_Fmsy_sglpop, data_list), parameters = par_list, map = NULL, silent = TRUE)
      if(n_pop > 1) tmp_obj <- RTMB::MakeADFun(cmb(local_BH_Fmsy_multipop, data_list), parameters = par_list, map = NULL, silent = TRUE)
      tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt <- tmp_obj$rep$Fmsy


      if(n_pop > 1) {
        for(r in 1:n_regions) {
          b_ref_pt[,r]        <- tmp_obj$rep$SB[,r]  * tmp_obj$rep$Req_o
          virgin_b_ref_pt[,r] <- tmp_obj$rep$SB0[,r] * data_list$R0
        }

        for(p2 in 1:n_pop) {
          r_natal <- data$natal_region[p2]
          pop_b_ref_pt[p2, r_natal]        <- tmp_obj$rep$SB[p2, r_natal]  * tmp_obj$rep$Req_o[p2]
          virgin_pop_b_ref_pt[p2, r_natal] <- tmp_obj$rep$SB0[p2, r_natal] * data_list$R0[p2]
          for(p in 1:n_pop) {
            if(p != p2) {
              pop_b_ref_pt[p2, r_natal]        <- pop_b_ref_pt[p2, r_natal]        + stray_rate[p] * tmp_obj$rep$SB[p, r_natal]  * tmp_obj$rep$Req_o[p]
              virgin_pop_b_ref_pt[p2, r_natal] <- virgin_pop_b_ref_pt[p2, r_natal] + stray_rate[p] * tmp_obj$rep$SB0[p, r_natal] * data_list$R0[p]
            }
          }
        }
      } else{
        for(r in 1:n_regions) {
          b_ref_pt[1, r]        <- tmp_obj$rep$SB[r] * tmp_obj$rep$Req_o[r]
          virgin_b_ref_pt[1, r] <- tmp_obj$rep$SB0[r] * data_list$R0 * data_list$rec_region_prop[r]
        }
        pop_b_ref_pt[1, 1]        <- sum(b_ref_pt)
        virgin_pop_b_ref_pt[1, 1] <- sum(virgin_b_ref_pt)
      }

      # see if Newton Raphson calcs for equil rec converged
      # if(sum(tmp_obj$rep$iter_vec) > 1e-10) warning("Calculations for equilibrium recruits from origin regions might not have converged! Try increasing local_bh_msy_newton_steps or be wary of these values!")
      # if(sum(tmp_obj$rep$Fmsy) == sum(exp(par_list$log_Fmsy))) warning("It is unlikely this converged. Starting values of log Fmsy have not changed (specified at log (0.1).")
    }

  } # end multi region

  return(list(f_ref_pt = f_ref_pt,
              b_ref_pt = b_ref_pt,
              virgin_b_ref_pt = virgin_b_ref_pt,
              pop_b_ref_pt = pop_b_ref_pt,
              virgin_pop_b_ref_pt = virgin_pop_b_ref_pt))

}


