# Stage 3 of 3: post fit
#
# Spawning biomass per recruit solvers: single_region_SPR for one region at a
# time, global_SPR for a single fishing mortality shared across regions.

#' Compute SPR reference point for a single-region or non-spatial model
#'
#' Calculates the spawning potential ratio (SPR) as a function of a trial
#' fishing mortality \eqn{F_x}, then returns a squared penalty
#' \eqn{100 (SPR - SPR_x)^2} that is minimised by the outer optimizer to find
#' \eqn{F_{SPR_x}}. Supports multiple populations via stray rates but does not
#' include spatial movement.
#'
#' @details
#' **Fishing mortality decomposition**
#'
#' Fishing mortality at age is split into:
#'
#' - retained fishing mortality
#'   \code{F_ret = F * selectivity * retention}
#'
#' - discard fishing mortality (dead discards only)
#'   \code{F_disc = F * selectivity * (1 - retention) * dmr}
#'
#' where \code{dmr} is the discard mortality rate (fraction of discarded fish
#' that die). Only the dead fraction contributes to total instantaneous
#' mortality \code{Z}.
#'
#' The total mortality used for survival is:
#'
#' \code{Z = M + F_ret + F_disc}
#'
#' This formulation assumes:
#'
#' - retained fish always die,
#' - only a fraction \code{dmr} of discarded fish die,
#' - the surviving fraction \code{(1 - dmr)} of discards remains in the
#'   population and continues aging and contributing to spawning biomass.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_F_x}}{Log-scale trial fishing mortality.}
#'   }
#'
#' @param data Named list of RTMB data. Must contain:
#'   \describe{
#'     \item{\code{n_pop}}{Integer. Number of populations.}
#'     \item{\code{n_ages}}{Integer. Number of age classes.}
#'     \item{\code{n_seas}}{Integer. Number of seasons.}
#'     \item{\code{seasdur}}{Numeric vector \code{[n_seas]}. Season durations.}
#'     \item{\code{spawn_seas}}{Integer. Index of the spawning season.}
#'     \item{\code{t_spawn}}{Numeric. Fraction of spawning season elapsed
#'       before spawning (mid-season mortality correction).}
#'
#'     \item{\code{F_fract_flt}}{Numeric array
#'       \code{[n_seas, n_fish_fleets]}. Fleet F fractions.}
#'
#'     \item{\code{fish_sel}}{Numeric array
#'       \code{[n_pop, n_seas, n_ages, n_fish_fleets]}.
#'       Fishery selectivity at age for females.}
#'
#'     \item{\code{ret_sel}}{Numeric array
#'       \code{[n_pop, n_seas, n_ages, n_fish_fleets]}.
#'       Retention selectivity (fraction of selected fish retained).}
#'
#'     \item{\code{dmr}}{Numeric array
#'       \code{[n_seas, n_fish_fleets]}. Discard mortality rate (fraction of
#'       discarded fish that die).}
#'
#'     \item{\code{natmort}}{Numeric array
#'       \code{[n_pop, n_ages]}. Female natural mortality at age.}
#'
#'     \item{\code{WAA}}{Numeric array
#'       \code{[n_pop, n_seas, n_ages]}. Female weight at age.}
#'
#'     \item{\code{MatAA}}{Numeric array
#'       \code{[n_pop, n_seas, n_ages]}. Maturity at age.}
#'
#'     \item{\code{sex_ratio_f}}{Numeric vector
#'       \code{[n_pop]}. Female sex ratio at recruitment.}
#'
#'     \item{\code{rec_seas_prop}}{Numeric array
#'       \code{[n_pop, n_seas]}. Proportion of annual recruitment entering
#'       in each season.}
#'
#'     \item{\code{stray_rate}}{Numeric vector
#'       \code{[n_pop]}. Per-population stray rate used to compute effective
#'       SSB across populations.}
#'
#'     \item{\code{natal_region}}{Integer vector
#'       \code{[n_pop]}. Natal region index for each population.}
#'
#'     \item{\code{n_pop_in_region}}{Integer vector
#'       \code{[n_regions]}. Number of populations per natal region.}
#'
#'     \item{\code{SPR_x}}{Numeric. Target SPR fraction (e.g. 0.4).}
#'   }
#'
#' @return Numeric scalar. Squared penalty \eqn{(SPR - SPR_x)^2}.
#'   Minimised to zero at \eqn{F = F_{SPR_x}}.
#'
#' @seealso \code{\link{Get_Reference_Points}}
#' @keywords internal
#' @import RTMB

single_region_SPR <- function(pars,
                              data
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  F_x = exp(log_F_x) # Exponentiate reference point
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_ages)) # 2 slots in rows, for unfished, and fished at F_x

  # Set up the initial recruits
  for(p in 1:n_pop) Nspr[, p, 1] <- sex_ratio_f[p] * rec_seas_prop[p, 1]

  ## Loop through ages
  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        # compute mortality
        F_seas = sum(F_fract_flt[seas, ] * F_x * fish_sel[p,seas,j - 1, ] * ret_sel[p,seas,j - 1,]) + # retained F
          sum(F_fract_flt[seas, ] * F_x * fish_sel[p,seas,j - 1, ] * (1 - ret_sel[p,seas,j - 1,]) * dmr[seas,]) # discard F
        M_seas = natmort[p, j - 1] * seasdur[seas]
        Z_seas = F_seas + M_seas

        # get recruits out
        tmp_unfished = Nspr[1,p, j - 1]; tmp_fished   = Nspr[2,p, j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          add_rec = rec_seas_prop[p,seas] * sex_ratio_f[p]
          tmp_unfished = tmp_unfished + add_rec
          tmp_fished = tmp_fished   + add_rec
        }

        ## Spawning biomass
        if (seas == spawn_seas) {
          SB_age[1, p, j - 1] = tmp_unfished * WAA[p, spawn_seas, j - 1] *
            MatAA[p, spawn_seas, j - 1] * exp(-t_spawn * M_seas)
          SB_age[2, p, j - 1] = tmp_fished * WAA[p, spawn_seas, j - 1] *
            MatAA[p, spawn_seas, j - 1] * exp(-t_spawn * Z_seas)
        }

        ## Mortality and ageing
        if (seas < n_seas) { # Within season mortality
          Nspr[1, p, j - 1] = tmp_unfished * exp(-M_seas)
          Nspr[2, p, j - 1] = tmp_fished * exp(-Z_seas)
        } else {
          # Ageing
          Nspr[1, p,  j] = tmp_unfished * exp(-M_seas)
          Nspr[2, p,  j] = tmp_fished * exp(-Z_seas)
        }

      } # end seas loop
    } # end j loop
  } # end p loop

  # Advance penultimate age to spawning season
  tmp_unfished = Nspr[1,,n_ages-1]
  tmp_fished = Nspr[2,,n_ages-1]

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal mortality
        tmp_unfished[p] = tmp_unfished[p] * exp(-(natmort[p,n_ages-1] * seasdur[seas]))
        tmp_fished[p] = tmp_fished[p] * exp(-(natmort[p,n_ages-1] * seasdur[seas] +
                                                sum(F_fract_flt[seas,] * F_x * fish_sel[p,seas,n_ages-1,] * ret_sel[p,seas,n_ages-1, ]) +
                                                sum(F_fract_flt[seas,] * F_x * fish_sel[p,seas,n_ages-1,] * (1 - ret_sel[p,seas,n_ages-1, ]) * dmr[seas,] )))
      } # end seas loop
    } # end p loop
  }

  # Get spawning biomass of penultimate age
  for(p in 1:n_pop) {
    SB_age[1, p, n_ages - 1] = tmp_unfished[p] * WAA[p,spawn_seas, n_ages - 1] * MatAA[p,spawn_seas, n_ages - 1] *
      exp(-t_spawn * natmort[p,n_ages - 1] * seasdur[spawn_seas])
    SB_age[2, p, n_ages - 1] = tmp_fished[p] * WAA[p,spawn_seas, n_ages - 1] * MatAA[p,spawn_seas, n_ages - 1] *
      exp(-t_spawn * (natmort[p,n_ages - 1] * seasdur[spawn_seas] +
                        sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[p,spawn_seas,n_ages-1,] * ret_sel[p,spawn_seas,n_ages-1, ]) +
                        sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[p,spawn_seas,n_ages-1,] * (1 - ret_sel[p,spawn_seas,n_ages-1, ]) * dmr[spawn_seas,] )))
    }

  # Plus group (scalar, no movement)
  for(p in 1:n_pop) {
    F_annual_penult = 0; F_annual_plus = 0
    for(seas in 1:n_seas) {
      F_annual_penult = F_annual_penult + sum(F_fract_flt[seas,] * F_x * fish_sel[p,seas,n_ages-1,] * ret_sel[p,seas,n_ages-1,]) +
        sum(F_fract_flt[seas, ] * F_x * fish_sel[p,seas,n_ages - 1, ] * (1 - ret_sel[p,seas,n_ages - 1,]) * dmr[seas,])
      F_annual_plus = F_annual_plus + sum(F_fract_flt[seas,] * F_x * fish_sel[p,seas,n_ages,] * ret_sel[p,seas,n_ages,]) +
        sum(F_fract_flt[seas, ] * F_x * fish_sel[p,seas,n_ages, ] * (1 - ret_sel[p,seas,n_ages,]) * dmr[seas,])
    }
    M_penult = natmort[p,n_ages - 1]
    M_plus = natmort[p,n_ages]
    Nspr[1,p,n_ages] = Nspr[1,p,n_ages-1] * exp(-M_penult) / (1 - exp(-M_plus))
    Nspr[2,p,n_ages] = Nspr[2,p,n_ages-1] * exp(-(M_penult + F_annual_penult)) / (1 - exp(-(M_plus + F_annual_plus)))
  }

  # Advance plus group to spawning season
  tmp_unfished = Nspr[1,,n_ages]
  tmp_fished = Nspr[2,,n_ages]

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {

        # Apply seasonal mortality
        tmp_unfished[p] = tmp_unfished[p] * exp(-(natmort[p,n_ages] * seasdur[seas]))
        tmp_fished[p] = tmp_fished[p] * exp(-(natmort[p,n_ages] * seasdur[seas] +
                                                sum(F_fract_flt[seas,] * F_x * fish_sel[p,seas,n_ages,] * ret_sel[p,seas,n_ages, ]) +
                                                sum(F_fract_flt[seas,] * F_x * fish_sel[p,seas,n_ages,] * (1 - ret_sel[p,seas,n_ages, ]) * dmr[seas,] )))

      } # end seas loop
    } # end p loop
  }

  ## Plus group spawning biomass
  for(p in 1:n_pop) {
    SB_age[1,p,n_ages] = tmp_unfished[p] * WAA[p,spawn_seas, n_ages] * MatAA[p,spawn_seas, n_ages] *
      exp(-t_spawn * natmort[p,n_ages] * seasdur[spawn_seas])
    SB_age[2,p,n_ages] = tmp_fished[p] * WAA[p,spawn_seas, n_ages] * MatAA[p,spawn_seas, n_ages] *
      exp(-t_spawn * (natmort[p,n_ages] * seasdur[spawn_seas] +
                        sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[p,spawn_seas,n_ages,] * ret_sel[p,spawn_seas,n_ages, ]) +
                        sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[p,spawn_seas,n_ages,] * (1 - ret_sel[p,spawn_seas,n_ages, ]) * dmr[spawn_seas,] )))
  }

  # Get effective SB after straying
  effective_SB = array(0, n_pop)
  effective_SB0 = array(0, n_pop)

  # Get spawning biomass per recruit by population
  SB = apply(SB_age[2,,,drop = FALSE], 2, sum)
  SB0 = apply(SB_age[1,,,drop = FALSE], 2, sum)

  # compute effective SSB after straying (divide by n_pop_in_region to preserve mass balance)
  for(p2 in 1:n_pop) {
    r <- natal_region[p2]
    for(p in 1:n_pop) {
      sc <- if(p == p2) 1 else stray_rate[p] / n_pop_in_region[r]
      effective_SB[p2]  <- effective_SB[p2]  + sc * SB[p]
      effective_SB0[p2] <- effective_SB0[p2] + sc * SB0[p]
    }
  }

  # Get spawning biomass per recruit to get spawning potential ratio
  SPR = sum(effective_SB) / sum(effective_SB0)
  sprpen =  (SPR - SPR_x)^2

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(effective_SB0)
  RTMB::REPORT(effective_SB)
  RTMB::REPORT(F_x)
  RTMB::REPORT(SB)
  RTMB::REPORT(SB0)

  return(sprpen)
}

#' Compute global SPR reference point for a spatially explicit model
#'
#' Calculates a single, spatially-integrated SPR using a per-recruit cohort
#' tracked across all regions and seasons under movement. A single scalar
#' \eqn{F_x} is applied uniformly across regions (scaled by region-specific
#' fleet fractions and selectivity). Returns the squared penalty
#' \eqn{100 (SPR - SPR_x)^2} for optimisation.
#'
#' Supports single- and multi-population models. When \code{n_pop > 1},
#' effective SSB at each population's natal region accumulates straying
#' contributions from other populations. When \code{n_seas = 1} and
#' \code{n_pop > 1}, \code{sgl_seas_spawning_movement} redistributes fish to
#' natal grounds before SSB is computed.
#'
#' The plus-group is solved analytically using
#' \code{\link{build_plus_group_T}} and \code{\link{solve_plus_group}}.
#'
#' @details
#' **Fishing mortality decomposition**
#'
#' Fishing mortality at age is split into:
#'
#' - retained fishing mortality
#'   \code{F_ret = F * selectivity * retention}
#'
#' - discard fishing mortality (dead discards only)
#'   \code{F_disc = F * selectivity * (1 - retention) * dmr}
#'
#' where \code{dmr} is the discard mortality rate (fraction of discarded fish
#' that die). Only the dead fraction contributes to total instantaneous
#' mortality \code{Z}.
#'
#' The total mortality used for survival is:
#'
#' \code{Z = M + F_ret + F_disc}
#'
#' This formulation assumes:
#'
#' - retained fish always die,
#' - only a fraction \code{dmr} of discarded fish die,
#' - the surviving fraction \code{(1 - dmr)} of discards remains in the
#'   population and continues aging, moving, and contributing to spawning
#'   biomass.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_F_x}}{Log-scale trial fishing mortality.}
#'   }
#'
#' @param data Named list of RTMB data. Must contain:
#'   \describe{
#'     \item{\code{n_pop}}{Integer. Number of populations.}
#'     \item{\code{n_regions}}{Integer. Number of spatial regions.}
#'     \item{\code{n_ages}}{Integer. Number of age classes.}
#'     \item{\code{n_seas}}{Integer. Number of seasons.}
#'     \item{\code{seasdur}}{Numeric vector \code{[n_seas]}. Season durations.}
#'     \item{\code{spawn_seas}}{Integer. Index of the spawning season.}
#'     \item{\code{t_spawn}}{Numeric. Mid-season spawning timing correction.}
#'
#'     \item{\code{F_fract_flt}}{Numeric array
#'       \code{[n_regions, n_seas, n_fish_fleets]}. Fleet F fractions by region.}
#'
#'     \item{\code{fish_sel}}{Numeric array
#'       \code{[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]}.
#'       Female fishery selectivity.}
#'
#'     \item{\code{ret_sel}}{Numeric array
#'       \code{[n_pop, n_regions, n_seas, n_ages, n_fish_fleets]}.
#'       Retention selectivity (fraction of selected fish that are retained).}
#'
#'     \item{\code{dmr}}{Numeric array
#'       \code{[n_regions, n_seas, n_fish_fleets]}. Discard mortality rate
#'       (fraction of discarded fish that die).}
#'
#'     \item{\code{natmort}}{Numeric array
#'       \code{[n_pop, n_regions, n_ages]}. Female natural mortality.}
#'
#'     \item{\code{WAA}}{Numeric array
#'       \code{[n_pop, n_regions, n_seas, n_ages]}. Female weight at age.}
#'
#'     \item{\code{MatAA}}{Numeric array
#'       \code{[n_pop, n_regions, n_seas, n_ages]}. Maturity at age.}
#'
#'     \item{\code{Movement}}{Numeric array
#'       \code{[n_pop, n_regions, n_regions, n_seas, n_ages]}.
#'       Seasonal movement transition matrices.}
#'
#'     \item{\code{sgl_seas_spawning_movement}}{Numeric array
#'       \code{[n_pop, n_regions, n_regions, n_ages]}. Spawning movement for
#'       single-season natal homing models.}
#'
#'     \item{\code{do_recruits_move}}{Integer (0/1). Whether age-1 recruits
#'       are subject to movement.}
#'
#'     \item{\code{rec_region_prop}}{Numeric array
#'       \code{[n_pop, n_regions]}. Proportion of recruitment entering each region.}
#'
#'     \item{\code{sex_ratio_f}}{Numeric array
#'       \code{[n_pop, n_regions]}. Female sex ratio at recruitment.}
#'
#'     \item{\code{rec_seas_prop}}{Numeric array
#'       \code{[n_pop, n_seas]}. Seasonal recruitment proportions.}
#'
#'     \item{\code{stray_rate}}{Numeric vector \code{[n_pop]}. Per-population
#'       stray rate.}
#'
#'     \item{\code{natal_region}}{Integer vector \code{[n_pop]}. Natal region
#'       index for each population.}
#'
#'     \item{\code{n_pop_in_region}}{Integer vector \code{[n_regions]}. Number
#'       of populations per natal region.}
#'
#'     \item{\code{SPR_x}}{Numeric. Target SPR fraction.}
#'   }
#'
#' @return Numeric scalar. Squared penalty \eqn{(SPR - SPR_x)^2}.
#'
#' @keywords internal
#' @import RTMB
global_SPR <- function(pars,
                       data
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # As with rec_model, expm_nsub is only present on data lists built since the
  # implicit matrix exponential option, so fall back to the exact exponential.
  if(!exists("expm_nsub", inherits = FALSE)) expm_nsub <- 0

  F_x = exp(log_F_x) # Exponentiate reference points

  # helper function to get retained f by region for a given age and season
  ret_F_by_region <- function(p, age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + (F_fract_flt[, seas, f] * F_x * fish_sel[p, , seas, age, f] * ret_sel[p,,seas,age,f])
    }
    Fr
  }

  # helper function to get discarded f by region for a given age and season
  disc_F_by_region <- function(p, age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + (F_fract_flt[, seas, f] * F_x * fish_sel[p, , seas, age, f] * (1 - ret_sel[p,,seas,age,f]) * dmr[, seas, f])
    }
    Fr
  }

  SB_age = Nspr = array(0, dim = c(2, n_pop, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy

  # Set up the initial recruits
  for(p in 1:n_pop) Nspr[2,p,,1] = Nspr[1,p,,1] = rec_region_prop[p,] * sex_ratio_f[p,] * rec_seas_prop[p,1]

  ## Loop through ages
  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        # get mortality
        ret_F_seas = ret_F_by_region(p,j - 1, seas)
        disc_F_seas = disc_F_by_region(p,j - 1, seas)
        M_seas = natmort[p,, j - 1] * seasdur[seas]
        Z_seas = ret_F_seas + M_seas + disc_F_seas

        # extract out quantities
        tmp_unfished = Nspr[1,p,, j - 1]; tmp_fished = Nspr[2,p,, j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          add = rec_seas_prop[p, seas] * sex_ratio_f[p, ] * rec_region_prop[p, ]
          tmp_unfished = tmp_unfished + add
          tmp_fished = tmp_fished + add
        }

        ## Movement operators for this age. Ages that do not move get an identity
        ## transition and a zero generator, which leaves survival unchanged under
        ## every move_timing.
        if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
          Mv = Movement[p,,,seas, j - 1]
          Qv = Mrate[p,,,seas, j - 1]
        } else {
          Mv = diag(n_regions)
          Qv = matrix(0, n_regions, n_regions)
        }

        ## Spawning biomass
        if (seas == spawn_seas) {

          # Propagate to the spawning point; movement and t_spawn mortality are applied
          # together so the spawning location is consistent with move_timing.
          tmp_unfished_spawn = spawn_state(tmp_unfished, Mv, M_seas, Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
          tmp_fished_spawn   = spawn_state(tmp_fished,   Mv, Z_seas, Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)

          # If single season natal homing population
          if(n_pop > 1 && n_seas == 1) {
            # Get NAA during spawning in single season case
            tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
            tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
          }

          # Get spawning biomass per recruit
          SB_age[1,p,, j - 1] = tmp_unfished_spawn * WAA[p,, spawn_seas, j - 1] * MatAA[p,, spawn_seas, j - 1]
          SB_age[2,p,, j - 1] = tmp_fished_spawn * WAA[p,, spawn_seas, j - 1] * MatAA[p,, spawn_seas, j - 1]

        }

        ## Movement, mortality and ageing
        adv_unfished = advance_seas(tmp_unfished, Mv, M_seas, Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
        adv_fished   = advance_seas(tmp_fished,   Mv, Z_seas, Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)

        if (seas < n_seas) { # Within season mortality
          Nspr[1,p,, j - 1] = adv_unfished
          Nspr[2,p,, j - 1] = adv_fished
        } else {
          # Ageing
          Nspr[1,p,, j] = adv_unfished
          Nspr[2,p,, j] = adv_fished
        }
      }
    }
  }

  # Advance penultimate age into spawning season
  tmp_unfished = array(Nspr[1,,,n_ages-1], dim = c(n_pop, n_regions))
  tmp_fished = array(Nspr[2,,,n_ages-1], dim = c(n_pop, n_regions))

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {

        # Get mortality
        M_seas = natmort[p,,n_ages-1] * seasdur[seas]
        ret_F_seas = ret_F_by_region(p, n_ages - 1, seas)
        disc_F_seas = disc_F_by_region(p, n_ages - 1, seas)

        # Apply seasonal movement and mortality together, per move_timing
        tmp_unfished[p,] = advance_seas(tmp_unfished[p,], Movement[p,,,seas,n_ages-1],
                                        M_seas, Mrate[p,,,seas,n_ages-1], seasdur[seas], move_timing, expm_nsub = expm_nsub)
        tmp_fished[p,]   = advance_seas(tmp_fished[p,], Movement[p,,,seas,n_ages-1],
                                        M_seas + ret_F_seas + disc_F_seas,
                                        Mrate[p,,,seas,n_ages-1], seasdur[seas], move_timing, expm_nsub = expm_nsub)

      } # end seas loop
    }
  }

  ## Advance penultimate age into spawning biomass / season
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn = tmp_fished
  for(p in 1:n_pop) {
    # Propagate to the spawning point; t_spawn mortality is folded in here rather than
    # applied to SB_age below, so that continuous movement can redistribute spawners.
    Zu_spawn = natmort[p,, n_ages - 1] * seasdur[spawn_seas]
    Zf_spawn = Zu_spawn + ret_F_by_region(p, n_ages - 1, spawn_seas) + disc_F_by_region(p, n_ages - 1, spawn_seas)
    tmp_unfished_spawn[p,] = spawn_state(tmp_unfished_spawn[p,], Movement[p,,, spawn_seas, n_ages - 1],
                                         Zu_spawn, Mrate[p,,, spawn_seas, n_ages - 1], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
    tmp_fished_spawn[p,]   = spawn_state(tmp_fished_spawn[p,], Movement[p,,, spawn_seas, n_ages - 1],
                                         Zf_spawn, Mrate[p,,, spawn_seas, n_ages - 1], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
    # If single season natal homing population
    if(n_pop > 1 && n_seas == 1) {
      # Get NAA during spawning in single season case
      tmp_unfished_spawn[p,] = as.vector(tmp_unfished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages-1])
      tmp_fished_spawn[p,] = as.vector(tmp_fished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages-1])
    }

  }

  for(p in 1:n_pop) {
    SB_age[1,p,, n_ages - 1] = tmp_unfished_spawn[p,] * WAA[p,, spawn_seas, n_ages - 1] * MatAA[p,, spawn_seas, n_ages - 1]
    SB_age[2,p,, n_ages - 1] = tmp_fished_spawn[p,] * WAA[p,, spawn_seas, n_ages - 1] * MatAA[p,, spawn_seas, n_ages - 1]
  }

  ## Plus group analytical solution
  for(p in 1:n_pop) {

    F_penult = F_plus = matrix(0, n_regions, n_seas)
    for(seas in 1:n_seas) {
      F_penult[,seas] = ret_F_by_region(p, n_ages - 1, seas) + disc_F_by_region(p, n_ages - 1, seas)
      F_plus[,seas] = ret_F_by_region(p, n_ages, seas) + disc_F_by_region(p, n_ages, seas)
    }

    Ts <- build_plus_group_T(
      M_penult   = natmort[p,, n_ages - 1],
      M_plus     = natmort[p,, n_ages],
      F_penult   = F_penult,
      F_plus     = F_plus,
      Mov_penult = array(Movement[p,,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
      Mov_plus   = array(Movement[p,,,, n_ages], dim = c(n_regions, n_regions, n_seas)),
      n_regions  = n_regions,
      n_seas = n_seas,
      seasdur = seasdur,
      Mrate_penult = array(Mrate[p,,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
      Mrate_plus   = array(Mrate[p,,,, n_ages], dim = c(n_regions, n_regions, n_seas)),
      move_timing = move_timing,
      expm_nsub = expm_nsub)

    pg <- solve_plus_group(Ts, Nspr[1, p,, n_ages - 1], Nspr[2, p,, n_ages - 1], n_regions)
    Nspr[1, p,, n_ages] <- pg$unfished
    Nspr[2, p,, n_ages] <- pg$fished
  }


  ## Plus group: advance to spawning season
  tmp_unfished = array(Nspr[1,,,n_ages], dim = c(n_pop, n_regions))
  tmp_fished = array(Nspr[2,,,n_ages], dim = c(n_pop, n_regions))

  if(spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {
        F_seas_plus = ret_F_by_region(p, n_ages, seas) + disc_F_by_region(p, n_ages, seas)
        M_seas_plus = natmort[p,,n_ages] * seasdur[seas]
        tmp_unfished[p,] = advance_seas(tmp_unfished[p,], Movement[p,,,seas,n_ages],
                                        M_seas_plus, Mrate[p,,,seas,n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
        tmp_fished[p,]   = advance_seas(tmp_fished[p,], Movement[p,,,seas,n_ages],
                                        M_seas_plus + F_seas_plus,
                                        Mrate[p,,,seas,n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
      }
    }
  }

  # Extract temporary variables out
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn = tmp_fished

  ## Plus group spawning biomass
  for(p in 1:n_pop) {
    # As for the penultimate age, t_spawn mortality is folded into the propagation
    # so that continuous movement can redistribute spawners before spawning.
    Zu_spawn_plus = natmort[p,,n_ages] * seasdur[spawn_seas]
    Zf_spawn_plus = Zu_spawn_plus + ret_F_by_region(p, n_ages, spawn_seas) + disc_F_by_region(p, n_ages, spawn_seas)
    tmp_unfished_spawn[p,] = spawn_state(tmp_unfished_spawn[p,], Movement[p,,, spawn_seas, n_ages],
                                         Zu_spawn_plus, Mrate[p,,, spawn_seas, n_ages], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
    tmp_fished_spawn[p,]   = spawn_state(tmp_fished_spawn[p,], Movement[p,,, spawn_seas, n_ages],
                                         Zf_spawn_plus, Mrate[p,,, spawn_seas, n_ages], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)

    # If single season natal homing population
    if(n_pop > 1 && n_seas == 1) {
      # Get NAA during spawning in single season case
      tmp_unfished_spawn[p,] = as.vector(tmp_unfished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages])
      tmp_fished_spawn[p,] = as.vector(tmp_fished_spawn[p,] %*% sgl_seas_spawning_movement[p,,,n_ages])
    }
  }

  # Get spawning biomass per recruit
  for(p in 1:n_pop) {
    SB_age[1,p,,n_ages] = tmp_unfished_spawn[p,] * WAA[p,,spawn_seas,n_ages] * MatAA[p,,spawn_seas,n_ages]
    SB_age[2,p,,n_ages] = tmp_fished_spawn[p,] * WAA[p,,spawn_seas,n_ages] * MatAA[p,,spawn_seas,n_ages]
  }

  # Get effective SB after straying
  effective_SB = array(0, dim = n_pop)
  effective_SB0 = array(0, dim = n_pop)
  SB = apply(SB_age[2,,,,drop = FALSE], 2:3, sum)
  SB0 = apply(SB_age[1,,,,drop = FALSE], 2:3, sum)

  # Accumulate effective SSB at each population's natal region
  # across all source populations (captures stray contributions)
  # divide by n_pop_in_region to preserve mass balance when multiple pops share a natal region
  for(p2 in 1:n_pop) {
    r <- natal_region[p2]
    for(p in 1:n_pop) {
      sc <- if(p == p2) 1 else stray_rate[p] / n_pop_in_region[r]
      effective_SB[p2]  <- effective_SB[p2]  + sc * SB[p,  r]
      effective_SB0[p2] <- effective_SB0[p2] + sc * SB0[p, r]
    }
  }

  # Get spawning biomass per recruit to get spawning potential ratio
  SPR = if(n_pop == 1) sum(SB) / sum(SB0) else sum(effective_SB) / sum(effective_SB0)
  sprpen = (SPR - SPR_x)^2

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
