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

  n_ages = dim(fish_sel)[1] # number of ages

  # exponentitate reference points to "estimate"
  F_x = exp(log_F_x)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_ages)) # 2 slots in rows, for unfished, and fished at F_x

  # Set up the initial recruits
  Nspr[,1] = 1 * sex_ratio_f

  ## Loop through ages
  for (j in 2:(n_ages - 1)) {
    for (seas in 1:n_seas) {

      tmp_unfished = Nspr[1,j - 1]
      tmp_fished   = Nspr[2,j - 1]

      ## Spawning biomass
      if (seas == spawn_seas) {
        SB_age[1, j - 1] = tmp_unfished * WAA[spawn_seas, j - 1] * MatAA[spawn_seas, j - 1] * exp(-t_spawn * natmort[j - 1] * seasdur[seas])
        SB_age[2, j - 1] = tmp_fished * WAA[spawn_seas, j - 1] * MatAA[spawn_seas, j - 1] *
          exp(-t_spawn * (natmort[j - 1] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[j-1,]) ))
      }

      ## Mortality and ageing
      if (seas < n_seas) { # Within season mortality
        Nspr[1, j - 1] = tmp_unfished * exp(-natmort[j - 1] * seasdur[seas])
        Nspr[2, j - 1] = tmp_fished * exp(-(natmort[j - 1] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[j-1,]) ))
      } else {
        # Ageing
        Nspr[1, j] = tmp_unfished * exp(-natmort[j - 1] * seasdur[seas])
        Nspr[2, j] = tmp_fished * exp(-(natmort[j - 1] * seasdur[seas] +  sum(F_fract_flt[seas,] * F_x * fish_sel[j-1,]) ))
      }
    }
  }

  # Age n_ages-1 is now at start of year after the loop
  tmp_unfished = Nspr[1,n_ages-1]
  tmp_fished = Nspr[2,n_ages-1]

  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # Apply seasonal mortality
      tmp_unfished = tmp_unfished * exp(-(natmort[n_ages-1] * seasdur[seas]))
      tmp_fished = tmp_fished * exp(-(natmort[n_ages-1] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[n_ages-1,]) ))

    } # end seas loop
  }

  SB_age[1, n_ages - 1] = tmp_unfished * WAA[spawn_seas, n_ages - 1] * MatAA[spawn_seas, n_ages - 1] * exp(-t_spawn * natmort[n_ages - 1] * seasdur[spawn_seas])
  SB_age[2, n_ages - 1] = tmp_fished * WAA[spawn_seas, n_ages - 1] * MatAA[spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[n_ages - 1] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[n_ages-1,]) ))


  # Plus group
  M_penult = natmort[n_ages - 1]
  Z_penult = natmort[n_ages - 1] + F_x * sum(colSums(F_fract_flt) * fish_sel[n_ages - 1,])

  Z_plus = natmort[n_ages] + F_x * sum(colSums(F_fract_flt) * fish_sel[n_ages,])
  M_plus = natmort[n_ages]

  Nspr[1,n_ages] = Nspr[1,n_ages-1] * exp(-M_penult) / (1 - exp(-M_plus))
  Nspr[2,n_ages] = Nspr[2,n_ages-1] * exp(-Z_penult) / (1 - exp(-Z_plus))

  tmp_unfished = Nspr[1,n_ages]
  tmp_fished = Nspr[2,n_ages]

  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # Apply seasonal mortality
      tmp_unfished = tmp_unfished * exp(-(natmort[n_ages] * seasdur[seas]))
      tmp_fished = tmp_fished * exp(-(natmort[n_ages] * seasdur[seas] + sum(F_fract_flt[seas,] * F_x * fish_sel[n_ages,]) ))

    } # end seas loop
  }

  ## Plus group spawning biomass
  SB_age[1, n_ages] = tmp_unfished * WAA[spawn_seas, n_ages] * MatAA[spawn_seas, n_ages] * exp(-t_spawn * natmort[n_ages] * seasdur[spawn_seas])
  SB_age[2, n_ages] = tmp_fished * WAA[spawn_seas, n_ages] * MatAA[spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[n_ages] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[n_ages,])))

  # Get spawning biomass per recruit to get spawning potential ratio
  SB0 = sum(SB_age[1,])
  SB_F_x = sum(SB_age[2,])
  SPR = SB_F_x / SB0

  # compute objective function to get F_x
  sprpen = 100 * (SPR - SPR_x)^2

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB0)
  RTMB::REPORT(SB_F_x)
  RTMB::REPORT(F_x)

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

  n_regions = dim(fish_sel)[1] # number of regions
  n_ages = dim(fish_sel)[2] # number of model ages

  # exponentitate reference points to "estimate"
  F_x = exp(log_F_x)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy

  # Set up the initial recruits
  Nspr[1,,1] = Rec_Prop * sex_ratio_f
  Nspr[2,,1] = Rec_Prop * sex_ratio_f

  ## Loop through ages
  for (j in 2:(n_ages - 1)) {
    for (seas in 1:n_seas) {

      tmp_unfished = Nspr[1,, j - 1]
      tmp_fished   = Nspr[2,, j - 1]

      ## Movement
      if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
        tmp_unfished = as.vector(tmp_unfished %*% Movement[,,seas, j - 1])
        tmp_fished   = as.vector(tmp_fished   %*% Movement[,,seas, j - 1])
      }

      ## Spawning biomass
      if (seas == spawn_seas) {
        SB_age[1,, j - 1] = tmp_unfished * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] * exp(-t_spawn * natmort[, j - 1] * seasdur[seas])
        SB_age[2,, j - 1] = tmp_fished * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] * exp(-t_spawn *
                                                                                                       (natmort[, j - 1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,j-1,,drop = F]) ))
      }

      ## Mortality and ageing
      if (seas < n_seas) { # Within season mortality
        Nspr[1,, j - 1] = tmp_unfished * exp(-natmort[, j - 1] * seasdur[seas])
        Nspr[2,, j - 1] =
          tmp_fished * exp(-(natmort[, j - 1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,j-1,,drop = F]) ))
      } else {
        # Ageing
        Nspr[1,, j] = tmp_unfished * exp(-natmort[, j - 1] * seasdur[seas])
        Nspr[2,, j] = tmp_fished * exp(-(natmort[, j - 1] * seasdur[seas] +  rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,j-1,,drop = F]) ))
      }
    }
  }
  # Age n_ages-1 is now at start of year after the loop
  tmp_unfished = Nspr[1,,n_ages-1]
  tmp_fished = Nspr[2,,n_ages-1]

  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # Apply seasonal movement
      tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages-1]
      tmp_fished = tmp_fished %*% Movement[,,seas,n_ages-1]

      # Apply seasonal mortality
      tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages-1] * seasdur[seas]))
      tmp_fished = tmp_fished * exp(-(natmort[,n_ages-1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages - 1,,drop = F]) ))

    } # end seas loop
  }

  ## Penultimate age spawning biomass
  tmp_unfished = as.vector(tmp_unfished %*% Movement[,, spawn_seas, n_ages - 1])
  tmp_fished   = as.vector(tmp_fished %*% Movement[,, spawn_seas, n_ages - 1])
  SB_age[1,, n_ages - 1] = tmp_unfished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[, n_ages - 1] * seasdur[spawn_seas])
  SB_age[2,, n_ages - 1] = tmp_fished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] * exp(-t_spawn *
                                                                                                                (natmort[, n_ages - 1] * seasdur[spawn_seas] + rowSums(F_fract_flt[,spawn_seas,,drop = F] * F_x * fish_sel[,n_ages - 1,,drop = F]) ))

  ## Plus group analytical solution
  T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

  for (seas in 1:n_seas) {
    # Get survival
    S_penult_unfished = diag(exp(-natmort[, n_ages - 1] * seasdur[seas]), n_regions)
    S_plus_unfished = diag(exp(-natmort[, n_ages] * seasdur[seas]), n_regions)
    S_penult_fished = diag(exp(-(natmort[, n_ages - 1] * seasdur[seas] +  rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages-1,,drop = F]))), n_regions)
    S_plus_fished = diag(exp(-(natmort[, n_ages] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages,,drop = F]) )), n_regions)

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

  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # Apply seasonal movement
      tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages]
      tmp_fished = tmp_fished %*% Movement[,,seas,n_ages]

      # Apply seasonal mortality
      tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages] * seasdur[seas]))
      tmp_fished = tmp_fished * exp(-(natmort[,n_ages] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * F_x * fish_sel[,n_ages,,drop = F]) ))

    } # end seas loop
  }

  ## Plus group spawning biomass
  SB_age[1,, n_ages] = tmp_unfished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] * exp(-t_spawn * natmort[, n_ages] * seasdur[spawn_seas])
  SB_age[2,, n_ages] = tmp_fished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] * exp(-t_spawn *
                                                                                                    (natmort[, n_ages] * seasdur[spawn_seas] + rowSums(F_fract_flt[,spawn_seas,,drop = F] * F_x * fish_sel[,n_ages,,drop = F]) ))


  # Get spawning biomass per recruit to get spawning potential ratio
  SB0 = sum(SB_age[1,,])
  SB_F_x = sum(SB_age[2,,])
  SPR = SB_F_x / SB0

  # compute objective function to get F_x
  sprpen = 100 * (SPR - SPR_x)^2

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB0)
  RTMB::REPORT(SB_F_x)
  RTMB::REPORT(F_x)

  return(sprpen)
}

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

  n_ages = dim(fish_sel)[1]

  Fmsy = exp(log_Fmsy)

  SB_age = Nspr = array(0, dim = c(2, n_ages))
  CAA = array(0, c(n_seas, n_ages))

  Nspr[,1] = 1 * sex_ratio_f

  for (j in 2:(n_ages - 1)) {
    for (seas in 1:n_seas) {

      # Extract out numbers
      tmp_unfished = Nspr[1, j - 1]
      tmp_fished   = Nspr[2, j - 1]

      # Compute F and Z for this age/season
      F_a_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[j-1,])
      Z_a_seas = natmort[j - 1] * seasdur[seas] + F_a_seas

      # Spawning biomass
      if (seas == spawn_seas) {
        SB_age[1, j - 1] = tmp_unfished * WAA[spawn_seas, j - 1] * MatAA[spawn_seas, j - 1] *
          exp(-t_spawn * natmort[j - 1] * seasdur[seas])
        SB_age[2, j - 1] = tmp_fished * WAA[spawn_seas, j - 1] * MatAA[spawn_seas, j - 1] *
          exp(-t_spawn * (natmort[j - 1] * seasdur[seas] + sum(F_fract_flt[seas,] * Fmsy * fish_sel[j-1,])))
      }

      # Catch-at-age (Baranov)
      CAA[seas, j - 1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))

      # Mortality and ageing
      if (seas < n_seas) {
        Nspr[1, j - 1] = tmp_unfished * exp(-natmort[j - 1] * seasdur[seas])
        Nspr[2, j - 1] = tmp_fished * exp(-Z_a_seas)
      } else {
        Nspr[1, j] = tmp_unfished * exp(-natmort[j - 1] * seasdur[seas])
        Nspr[2, j] = tmp_fished * exp(-Z_a_seas)
      }
    }
  }

  # Age n_ages-1 at start of year after loop
  tmp_unfished = Nspr[1, n_ages - 1]
  tmp_fished = Nspr[2, n_ages - 1]

  # Catch-at-age for penultimate age
  for (seas in 1:n_seas) {
    F_a_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages-1,])
    Z_a_seas = natmort[n_ages-1] * seasdur[seas] + F_a_seas
    CAA[seas, n_ages-1] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
  }

  if (spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {
      # Exponential mortality
      tmp_unfished = tmp_unfished * exp(-natmort[n_ages - 1] * seasdur[seas])
      tmp_fished = tmp_fished * exp(-(sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages-1,]) + natmort[n_ages - 1]))
    }
  }

  # Get spawning biomass after mortality decrement
  SB_age[1, n_ages - 1] = tmp_unfished * WAA[spawn_seas, n_ages - 1] * MatAA[spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[n_ages - 1] * seasdur[spawn_seas])
  SB_age[2, n_ages - 1] = tmp_fished * WAA[spawn_seas, n_ages - 1] * MatAA[spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[n_ages - 1] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[n_ages-1,])))

  # Plus group
  Z_plus = natmort[n_ages] + sum(Fmsy * fish_sel[n_ages,])
  M_plus = natmort[n_ages]
  Nspr[1, n_ages] = Nspr[1, n_ages - 1] * exp(-M_plus) / (1 - exp(-M_plus))
  Nspr[2, n_ages] = Nspr[2, n_ages - 1] * exp(-Z_plus) / (1 - exp(-Z_plus))

  # Plus group catch across seasons
  tmp_unfished = Nspr[1, n_ages]
  tmp_fished = Nspr[2, n_ages]

  # Catch-at-age for plus group
  for (seas in 1:n_seas) {
    F_a_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages,])
    Z_a_seas = natmort[n_ages] * seasdur[seas] + F_a_seas
    CAA[seas, n_ages] = tmp_fished * (F_a_seas / Z_a_seas) * (1 - exp(-Z_a_seas))
  }

  if (spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {
      tmp_unfished = tmp_unfished * exp(-natmort[n_ages] * seasdur[seas])
      tmp_fished = tmp_fished * exp(-(natmort[n_ages] * seasdur[seas] + sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages,])))
    }
  }

  SB_age[1, n_ages] = tmp_unfished * WAA[spawn_seas, n_ages] * MatAA[spawn_seas, n_ages] *
    exp(-t_spawn * natmort[n_ages] * seasdur[spawn_seas])
  SB_age[2, n_ages] = tmp_fished * WAA[spawn_seas, n_ages] * MatAA[spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[n_ages] * seasdur[spawn_seas] + sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[n_ages,])))

  SBPR_0 = sum(SB_age[1,])
  SBPR_F = sum(SB_age[2,])

  Req = R0 * ((4 * h * SBPR_F) - (1 - h) * SBPR_0) / ((5 * h - 1) * SBPR_F)

  Yield = sum(CAA * WAA) * Req

  Bmsy = SBPR_F * Req
  B0 = SBPR_0 * R0

  obj_fun = -Yield

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SBPR_0)
  RTMB::REPORT(SBPR_F)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Yield)
  RTMB::REPORT(Bmsy)
  RTMB::REPORT(B0)
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

  n_regions = dim(fish_sel)[1] # number of regions
  n_ages = dim(fish_sel)[2] # number of model ages

  # exponentitate reference points to "estimate"
  Fmsy = exp(log_Fmsy)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_seas, n_ages)) # catch at age

  # Set up the initial recruits
  Nspr[1,,1] = Rec_Prop * sex_ratio_f
  Nspr[2,,1] = Rec_Prop * sex_ratio_f

  ## Loop through ages
  for (j in 2:(n_ages - 1)) {
    for (seas in 1:n_seas) {

      tmp_unfished = Nspr[1,, j - 1]
      tmp_fished   = Nspr[2,, j - 1]

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
  SB_age[1,, n_ages] = tmp_unfished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] * exp(-t_spawn * natmort[, n_ages] * seasdur[spawn_seas])
  SB_age[2,, n_ages] = tmp_fished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] * exp(-t_spawn *
                                                                                                    (natmort[, n_ages] * seasdur[spawn_seas] + rowSums(F_fract_flt[,spawn_seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F]) ))



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
  Bmsy = SBPR_F * Req
  B0 = SBPR_0 * R0

  # compute objective function to get Fmsy
  obj_fun = -Yield

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SBPR_0)
  RTMB::REPORT(SBPR_F)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Yield)
  RTMB::REPORT(Yield_r)
  RTMB::REPORT(Bmsy)
  RTMB::REPORT(B0)
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
local_BH_Fmsy <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  n_regions = dim(fish_sel)[1] # number of regions
  n_ages = dim(fish_sel)[2] # number of model ages

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_regions, n_seas, n_ages)) # catch at age
  Yield_r = array(0, dim = n_regions) # yield by region
  SB_unfished_mat = matrix(0, n_regions, n_regions)  # unfished spawning biomass per recruit
  SB_fished_mat = matrix(0, n_regions, n_regions) # fished spawning biomass per recruit
  Bmsy_r = array(0, dim = n_regions) # BMSY
  B0_r = array(0, dim = n_regions) # unfished B0
  SPR_r = array(0, dim = n_regions) # spawning potential ratio

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

        F_a_seas = rowSums(F_fract_flt[,seas,,drop=F] * Fmsy * fish_sel[,j-1,,drop=F])
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
          Nspr[2,o,,j-1] = tmp_fished * exp(-(natmort[,j-1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F]) ))
        } else {
          # Last season: mortality + ageing
          Nspr[1,o,,j] = tmp_unfished * exp(-(natmort[,j-1] * seasdur[seas]))
          Nspr[2,o,,j] = tmp_fished * exp(-(natmort[,j-1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F]) ))
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
      F_a_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F])
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
        tmp_fished = tmp_fished * exp(-(natmort[,n_ages-1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F]) ))

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
        exp(-t_spawn * ((natmort[d,n_ages-1] * seasdur[spawn_seas]) + sum(F_fract_flt[d,spawn_seas,] * Fmsy * fish_sel[d,n_ages-1,])))
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
      S_penult_fished = diag(exp(-(natmort[,n_ages-1] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F]))), n_regions)
      S_plus_fished = diag(exp(-(natmort[,n_ages] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F]))), n_regions)
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
        tmp_fished = tmp_fished * exp(-(natmort[,n_ages] * seasdur[seas] + rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F])))

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
        exp(-t_spawn * ((natmort[d,n_ages] * seasdur[spawn_seas]) + sum(F_fract_flt[d,spawn_seas,] * Fmsy * fish_sel[d,n_ages,])))
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

  A = 4 * h * Rec_Prop * R0 # define first part of the numerator of BH recruitment
  B = rep(0, n_regions) # define first part of the denominator of BH recruitment
  for(d in 1:n_regions) B[d] = (1 - h[d]) * sum(SB_unfished_mat[,d] * Rec_Prop * R0)
  C = 5 * h - 1 # define second part of the denominator for BH recruitment

  # define initial guess to solve for equilibrium recruitment from origin region
  Req_o = R0 * Rec_Prop

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

  # get other derived quantities
  for(d in 1:n_regions) {
    Bmsy_r[d] = sum(SB_fished_mat[,d] * Req_o)
    B0_r[d] = sum(SB_unfished_mat[,d] * R0 * Rec_Prop)
    SPR_r[d] = Bmsy_r[d] / B0_r[d]
  }

  # maximize total yield
  Yield_total = sum(Yield_r)
  obj_fun = -Yield_total

  sum_SB_unfished_mat = sum(SB_unfished_mat)

  # RTMB::REPORT(eqrec_prop)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Req_o)
  RTMB::REPORT(Bmsy_r)
  RTMB::REPORT(SPR_r)
  RTMB::REPORT(Yield_r)
  RTMB::REPORT(Yield_total)
  RTMB::REPORT(dg_dReq)
  RTMB::REPORT(B0_r)
  RTMB::REPORT(SPR_r)
  RTMB::REPORT(iter_vec)
  RTMB::REPORT(SB_fished_mat)
  RTMB::REPORT(SB_unfished_mat)
  RTMB::REPORT(sum_SB_unfished_mat)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SB_age)

  return(obj_fun)
}

#' Wrapper function to get reference points
#'
#' Wrapper function to compute fishing and biological reference points given data and report
#' objects from an assessment or simulation. Supports both single-region and multi-region
#' calculations with options for SPR or Beverton–Holt MSY reference points.
#'
#' @param data List. Data object containing ages, years, weight-at-age, maturity, natural mortality, and other simulation/assessment info.
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
                                 sex_ratio_f = rep(0.5, data$n_regions),
                                 calc_rec_st_yr = 1,
                                 rec_age = 1,
                                 type,
                                 what,
                                 n_avg_yrs = 1,
                                 local_bh_msy_newton_steps = 6
                                 ) {

  f_ref_pt <- vector() # set up storage
  b_ref_pt <- vector() # set up storage
  virgin_b_ref_pt <- vector() # set up storage

  # determine years to average over demogrphaics
  n_yrs <- length(data$years)
  avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs

  if(type == "single_region") {

    if(!what %in% c("SPR", "BH_MSY")) stop("what is not correctly specified! Should be SPR, BH_MSY for type = single_region")

    data_list <- list() # set up data list
    # Extract out relevant elements
    n_ages <- length(data$ages) # number of ages
    n_years <- length(data$years) # number of years

    # Seasonal stuff
    data_list$t_spawn <- t_spawn # specified mortality time up until spawning
    data_list$n_seas <- data$n_seas # number of seasons
    data_list$seasdur <- data$seasdur # seasonal duration
    data_list$spawn_seas <- data$spawn_seas # spawning season

    # fishing mortality fraction
    data_list$F_fract_flt <- array(rep$Fmort[1,n_years,,] / sum(rep$Fmort[1,n_years,,]),
                                   dim = c(data$n_seas, data$n_fish_fleets)) # get fleet F fraction to derive population level selectivity

    # fishery selectivity
    fish_sel_avg <- apply(rep$fish_sel[1,avg_yrs,,1,,drop = FALSE], c(3,5), mean)
    data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

    # natural mortality
    natmort_avg <- apply(rep$natmort[1,avg_yrs,,1,drop = FALSE], c(3), mean)
    data_list$natmort <- as.vector(natmort_avg) # get female natural mortality

    # weight at age
    WAA_avg <- apply(data$WAA[1,avg_yrs,,,1,drop = FALSE], c(3,4), mean)
    data_list$WAA <- array(WAA_avg, dim = c(data$n_seas, n_ages)) # weight at age for females

    # maturity at age
    MatAA_avg <- apply(data$MatAA[1,avg_yrs,,,1,drop = FALSE], c(3,4), mean)
    data_list$MatAA <- array(MatAA_avg, dim = c(data$n_seas, n_ages)) # maturity at age for females
    data_list$sex_ratio_f <- sex_ratio_f # recritment sex ratio

    if(what == 'SPR') {
      data_list$SPR_x <- SPR_x # SPR fraction
      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # Make adfun object
      obj <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt[1] <- obj$rep$F_x
      b_ref_pt[1] <- obj$rep$SB_F_x * mean(rep$Rec[1,calc_rec_st_yr:(n_years - rec_age)])
      virgin_b_ref_pt[1] <- obj$rep$SB0 * mean(rep$Rec[1,calc_rec_st_yr:(n_years - rec_age)])
    } # end SPR reference points

    if(what == 'BH_MSY') {

      # extract out beverton-holt parameters
      data_list$h <- rep$h_trans # steepness
      data_list$R0 <- rep$R0 # unfished recruitment

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # make adfun ect
      obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt[1] <- obj$rep$Fmsy
      b_ref_pt[1] <- obj$rep$Bmsy
      virgin_b_ref_pt[1] <- obj$rep$B0
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

    if(what == "independent_SPR") {
      for(r in 1:data$n_regions) {

        # Extract out relevant elements for a given region
        n_years <- length(data$years) # number of years
        n_ages <- length(data$ages) # number of ages
        data_list$F_fract_flt <- array(rep$Fmort[r,n_years,,] / sum(rep$Fmort[r,n_years,,]),
                                       dim = c(data$n_seas, data$n_fish_fleets)) # get fleet F fraction to derive population level selectivity

        # fishery selectivity
        fish_sel_avg <- apply(rep$fish_sel[r,avg_yrs,,1,,drop = FALSE], c(1, 3, 4, 5), mean)
        data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

        # natural mortality
        natmort_avg <- apply(rep$natmort[r,avg_yrs,,1,drop = FALSE], 3, mean)
        data_list$natmort <- as.vector(natmort_avg) # get female natural mortality

        # weight at age
        WAA_avg <- apply(data$WAA[r,avg_yrs,,,1,drop = FALSE], c(3,4), mean)
        data_list$WAA <- array(WAA_avg, dim = c(data$n_seas, n_ages)) # weight at age for females

        # maturity at age
        MatAA_avg <- apply(data$MatAA[r,avg_yrs,,,1,drop = FALSE], c(3,4), mean)
        data_list$MatAA <- array(MatAA_avg, dim = c(data$n_seas, n_ages)) # maturity at age for females
        data_list$SPR_x <- SPR_x # SPR fraction
        data_list$sex_ratio_f <- sex_ratio_f[r] # recritment sex ratio

        par_list <- list() # set up parameter list
        par_list$log_F_x <- log(0.1) # F_x starting value

        # Make adfun object
        tmp_obj <- RTMB::MakeADFun(cmb(single_region_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
        tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
        tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

        # Output reference points
        f_ref_pt[r] <- tmp_obj$rep$F_x
        b_ref_pt[r] <- tmp_obj$rep$SB_F_x * mean(rep$Rec[r,calc_rec_st_yr:(n_years - rec_age)])
        virgin_b_ref_pt[r] <- tmp_obj$rep$SB0 * mean(rep$Rec[r,calc_rec_st_yr:(n_years - rec_age)])

      } # end r loop
    } # end independent_SPR

    if(what == "independent_BH_MSY") {
      for(r in 1:data$n_regions) {

        # Extract out relevant elements for a given region
        n_years <- length(data$years) # number of years
        n_ages <- length(data$ages) # number of ages
        data_list$F_fract_flt <- rep$Fmort[r,n_years,] / sum(rep$Fmort[r,n_years,]) # get fleet F fraction to derive population level selectivity

        # fishery selectivity
        fish_sel_avg <- apply(rep$fish_sel[r,avg_yrs,,1,,drop = FALSE], c(1, 3, 4, 5), mean)
        data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

        # natural mortality
        natmort_avg <- apply(rep$natmort[r,avg_yrs,,1,drop = FALSE], c(1, 3, 4), mean)
        data_list$natmort <- as.vector(natmort_avg) # get female natural mortality

        # weight at age
        WAA_avg <- apply(data$WAA[r,avg_yrs,,1,drop = FALSE], c(1, 3, 4), mean)
        data_list$WAA <- WAA_avg # weight at age for females

        # maturity at age
        MatAA_avg <- apply(data$MatAA[r,avg_yrs,,1,drop = FALSE], c(1, 3, 4), mean)
        data_list$MatAA <- MatAA_avg # maturity at age for females

        # Beverton Holt parameters
        data_list$h <- rep$h_trans[r] # steepness
        data_list$R0 <- rep$R0 * rep$Rec_trans_prop[r] # unfished recruitment by region
        data_list$sex_ratio_f <- sex_ratio_f[r] # recritment sex ratio

        par_list <- list() # set up parameter list
        par_list$log_Fmsy <- log(0.1) # Fmsy starting value

        # Make adfun object
        tmp_obj <- RTMB::MakeADFun(cmb(single_region_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
        tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
        tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report

        # Output reference points
        f_ref_pt[r] <- tmp_obj$rep$Fmsy
        b_ref_pt[r] <- tmp_obj$rep$Bmsy
        virgin_b_ref_pt[r] <- tmp_obj$rep$B0

      } # end r loop
    } # end independent_SPR

    if(what == 'global_SPR') {

      # Extract out relevant elements for a given region
      n_ages <- length(data$ages) # number of ages to iterate through
      n_years <- length(data$years) # number of years
      n_regions <- data$n_regions # number of regions

      # Fleet fraction F
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio

      # fishery selectivity
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

      # natural mortality
      natmort_avg <- apply(rep$natmort[,avg_yrs,,1,drop = FALSE], c(1,3), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_regions, n_ages)) # get female natural mortality

      # weight at age
      WAA_avg <- apply(data$WAA[,avg_yrs,,,1,drop = FALSE], c(1, 3, 4), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # weight at age for females

      # maturity at age
      MatAA_avg <- apply(data$MatAA[,avg_yrs,,,1,drop = FALSE], c(1, 3, 4), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # maturity at age for females

      # Movement
      Movement_avg <- apply(rep$Movement[,,avg_yrs,,,1,drop = FALSE], c(1,2,4,5), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_regions, n_regions, data$n_seas,n_ages)) # Movement

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$Rec_Prop <- rep$Rec_trans_prop # recruitment proportions
      data_list$sex_ratio_f <- sex_ratio_f # recritment sex ratio

      data_list$SPR_x <- SPR_x # SPR fraction

      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # make adfn object
      obj <- RTMB::MakeADFun(cmb(global_SPR, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # output reference points
      f_ref_pt <- rep(obj$rep$F_x, n_regions)
      b_ref_pt <- obj$rep$SB_F_x * rowMeans(rep$Rec[,calc_rec_st_yr:(n_years - rec_age)])
      virgin_b_ref_pt <- obj$rep$SB0 * rowMeans(rep$Rec[,calc_rec_st_yr:(n_years - rec_age)])

    } # end global SPR

    if(what == 'global_BH_MSY') {

      # Extract out relevant elements for a given region
      n_ages <- length(data$ages) # number of ages to iterate through
      n_years <- length(data$years) # number of years
      n_regions <- data$n_regions # number of regions

      # Fleet fraction F
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio

      # fishery selectivity
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

      # natural mortality
      natmort_avg <- apply(rep$natmort[,avg_yrs,,1,drop = FALSE], c(1,3), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_regions, n_ages)) # get female natural mortality

      # weight at age
      WAA_avg <- apply(data$WAA[,avg_yrs,,,1,drop = FALSE], c(1, 3, 4), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # weight at age for females

      # maturity at age
      MatAA_avg <- apply(data$MatAA[,avg_yrs,,,1,drop = FALSE], c(1, 3, 4), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # maturity at age for females

      # Movement
      Movement_avg <- apply(rep$Movement[,,avg_yrs,,,1,drop = FALSE], c(1,2,4,5), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_regions, n_regions, data$n_seas,n_ages)) # Movement

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$Rec_Prop <- rep$Rec_trans_prop # recruitment proportions
      data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio to use
      data_list$h <- mean(rep$h_trans) # steepness
      data_list$R0 <- rep$R0  # unfished recruitment

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- log(0.1) # Fmsy starting value

      # Make adfun object
      obj <- RTMB::MakeADFun(cmb(global_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt <- rep(obj$rep$Fmsy, n_regions)
      b_ref_pt <- obj$rep$Bmsy * rep$Rec_trans_prop
      virgin_b_ref_pt <- obj$rep$B0 * rep$Rec_trans_prop
    }

    if(what == 'local_BH_MSY') {

      # Extract out relevant elements for a given region
      n_ages <- length(data$ages) # number of ages to iterate through
      n_years <- length(data$years) # number of years
      n_regions <- data$n_regions # number of regions

      # Fleet fraction F
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio

      # fishery selectivity
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets

      # natural mortality
      natmort_avg <- apply(rep$natmort[,avg_yrs,,1,drop = FALSE], c(1,3), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_regions, n_ages)) # get female natural mortality

      # weight at age
      WAA_avg <- apply(data$WAA[,avg_yrs,,,1,drop = FALSE], c(1, 3, 4), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # weight at age for females

      # maturity at age
      MatAA_avg <- apply(data$MatAA[,avg_yrs,,,1,drop = FALSE], c(1, 3, 4), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # maturity at age for females

      # Movement
      Movement_avg <- apply(rep$Movement[,,avg_yrs,,,1,drop = FALSE], c(1,2,4,5), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_regions, n_regions, data$n_seas,n_ages)) # Movement

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$Rec_Prop <- rep$Rec_trans_prop # recruitment proportions
      data_list$h <- rep$h_trans # steepness
      data_list$R0 <- rep$R0  # unfished recruitment
      data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio to use
      data_list$newton_steps <- local_bh_msy_newton_steps # number of newton steps to take

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- rep(log(0.1), n_regions) # Fmsy starting value

      # Make adfun object
      obj <- RTMB::MakeADFun(cmb(local_BH_Fmsy, data_list), parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
      obj$rep <- obj$report(obj$env$last.par.best) # get report

      # Output reference points
      f_ref_pt <- obj$rep$Fmsy
      b_ref_pt <- obj$rep$Bmsy_r
      virgin_b_ref_pt <- obj$rep$B0_r

      # see if Newton Raphson calcs for equil rec converged
      # if(sum(obj$rep$iter_vec) > 1e-10) warning("Calculations for equilibrium recruits from origin regions might not have converged! Try increasing local_bh_msy_newton_steps or be wary of these values!")
      # if(sum(obj$rep$Fmsy) == sum(exp(par_list$log_Fmsy))) warning("It is unlikely this converged. Starting values of log Fmsy have not changed (specified at log (0.1).")
    }

  } # end multi region

  return(list(f_ref_pt = f_ref_pt,
              b_ref_pt = b_ref_pt,
              virgin_b_ref_pt = virgin_b_ref_pt))

}


