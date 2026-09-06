# Stage 3 of 3: post fit
#
# Stock-recruit Fmsy solvers, differing in what is kept common across space: single region, one global
# F, or a local F per region. rec_model 0 has no curve at all and is rejected rather than folded in.

# Equilibrium Recruitment ---------------------------------------------------

#' Equilibrium recruitment from spawning biomass per recruit
#'
#' Solves R = f(R * phi) for R, where f is the stock-recruit curve. Both forms
#' pass through (S0, R0) and diverge away from it, so a model fitted with one and
#' given reference points from the other is silently wrong rather than noisy.
#'
#' @param phi Spawning biomass per recruit under fishing.
#' @param phi0 Unfished spawning biomass per recruit.
#' @param R0 Unfished recruitment.
#' @param h Steepness.
#' @param rec_model 1 for Beverton-Holt, 2 for Ricker.
#' @keywords internal
equil_rec_phi <- function(phi, phi0, R0, h, rec_model) {
  if(rec_model == 2) {
    a <- log(4 * h / (1 - h))
    R0 * (phi0 / phi) * (1 - log(phi0 / phi) / a)
  } else {
    R0 * ((4 * h * phi) - (1 - h) * phi0) / ((5 * h - 1) * phi)
  }
}

#' Equilibrium recruitment from spawning biomass
#'
#' The same two curves evaluated at a spawning biomass rather than a per-recruit
#' quantity, for the spatial solvers that iterate on recruitment directly.
#'
#' @param S Equilibrium spawning biomass.
#' @param S0 Unfished spawning biomass.
#' @param R0 Unfished recruitment.
#' @param h Steepness.
#' @param rec_model 1 for Beverton-Holt, 2 for Ricker.
#' @keywords internal
equil_rec_ssb <- function(S, S0, R0, h, rec_model) {
  if(rec_model == 2) {
    a <- log(4 * h / (1 - h))
    R0 * (S / S0) * exp(a * (1 - S / S0))
  } else {
    (4 * h * R0 * S) / ((1 - h) * S0 + (5 * h - 1) * S)
  }
}

#' Derivative of equilibrium recruitment with respect to spawning biomass
#'
#' Supplies the Jacobian term for the Newton solves in the spatial cases.
#'
#' @inheritParams equil_rec_ssb
#' @keywords internal
equil_rec_ssb_deriv <- function(S, S0, R0, h, rec_model) {
  if(rec_model == 2) {
    a <- log(4 * h / (1 - h))
    (R0 / S0) * exp(a * (1 - S / S0)) * (1 - a * S / S0)
  } else {
    (4 * h * R0 * (1 - h) * S0) / ((1 - h) * S0 + (5 * h - 1) * S)^2
  }
}

# Single Region MSY ---------------------------------------------------------

#' Compute Beverton-Holt Fmsy for a single-region or non-spatial model
#'
#' Finds \eqn{F_{MSY}} by maximizing equilibrium yield under a
#' Beverton-Holt stock-recruit relationship. Yield is computed from
#' spawning biomass per recruit (\eqn{\phi_F}), the BH equilibrium
#' recruitment formula, and catch-at-age integrated across all seasons.
#' Supports multiple populations via stray rates but does not include
#' spatial movement. Yield includes only landings from fleets where
#' \code{is_discard_fleet == 0}; discard-only fleets contribute to total
#' mortality but not to the yield being maximized.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_Fmsy}}{Log-scale trial \eqn{F_{MSY}}.}
#'   }
#' @param data Named list of RTMB data. Must contain all fields required by
#'   \code{\link{single_region_SPR}} plus:
#'   \describe{
#'     \item{\code{h}}{Numeric vector \code{[n_pop]}. Beverton-Holt steepness.}
#'     \item{\code{R0}}{Numeric vector \code{[n_pop]}. Unfished equilibrium
#'       recruitment.}
#'     \item{\code{is_discard_fleet}}{Integer vector \code{[n_fish_fleets]}.
#'       Indicator for fleets whose catch is excluded from landed yield
#'       (0 = landing fleet, 1 = discard-only fleet). These fleets still
#'       contribute to total fishing mortality \code{Z} and affect population
#'       dynamics and spawning biomass.}
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield (minimized to
#'   find \eqn{F_{MSY}}).
#'
#' @keywords internal
#' @import RTMB
single_region_Fmsy <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # getAll defines rec_model whenever the caller supplied it, so hand-built data lists predating
  # the Ricker fall back here rather than being assigned before getAll, which would collide
  if(!exists("rec_model", inherits = FALSE)) rec_model <- 1

  # rec_model 0 is mean recruitment, which has no curve to maximize yield over. the helpers below
  # branch on Ricker against everything else, so a mean recruitment fit would report Beverton-Holt
  if(rec_model == 0) stop("rec_model = 0 (mean recruitment) has no stock-recruit curve, so Fmsy is undefined. ",
                          "MSY reference points need rec_model 1 (Beverton-Holt) or 2 (Ricker); use SPR reference points instead.",
                          call. = FALSE)

  # annual total, read by the plus group series which steps a whole year
  natmort_annual = collapse_natmort_annual(natmort, seasdur, seas_dim = 2)

  # Exponentiate Fmsy
  Fmsy = exp(log_Fmsy)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_ages))
  CAA = array(0, c(n_pop, n_seas, n_ages))
  landed_fleets = which(is_discard_fleet == 0) # figure out discard fleets

  # initialize recruits
  for(p in 1:n_pop) Nspr[,p,1] = sex_ratio_f[p] * rec_seas_prop[p,1]

  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        # get mortality
        # compute mortality
        F_seas = sum(F_fract_flt[seas, ] * Fmsy * fish_sel[p,seas,j - 1, ] * ret_sel[p,seas,j - 1,]) + # retained F
          sum(F_fract_flt[seas, ] * Fmsy * fish_sel[p,seas,j - 1, ] * (1 - ret_sel[p,seas,j - 1,]) * dmr[seas,]) # discard F
        landed_F_seas = sum(F_fract_flt[seas, landed_fleets] * Fmsy * fish_sel[p,seas,j - 1, landed_fleets] * ret_sel[p,seas,j - 1,landed_fleets]) # landed F
        M_a_seas = natmort[p,seas, j - 1] * seasdur[seas]
        Z_seas = F_seas + M_a_seas

        # Extract out numbers
        tmp_unfished = Nspr[1, p, j - 1]
        tmp_fished   = Nspr[2, p,  j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          add_rec = rec_seas_prop[p,seas] * sex_ratio_f[p]
          tmp_unfished = tmp_unfished + add_rec
          tmp_fished   = tmp_fished   + add_rec
        }

        # Spawning biomass
        if (seas == spawn_seas) {
          SB_age[1, p, j - 1] = tmp_unfished * WAA[p, spawn_seas, j - 1] * MatAA[p, spawn_seas, j - 1] * exp(-t_spawn * M_a_seas)
          SB_age[2, p, j - 1] = tmp_fished * WAA[p, spawn_seas, j - 1] * MatAA[p, spawn_seas, j - 1] * exp(-t_spawn * Z_seas)
        }

        # Catch-at-age (Baranov)
        CAA[p, seas, j - 1] = tmp_fished * (landed_F_seas / Z_seas) * (1 - exp(-Z_seas))

        # Mortality and ageing
        if (seas < n_seas) {
          Nspr[1, p, j - 1] = tmp_unfished * exp(-M_a_seas)
          Nspr[2, p, j - 1] = tmp_fished * exp(-Z_seas)
        } else {
          Nspr[1, p, j] = tmp_unfished * exp(-M_a_seas)
          Nspr[2, p, j] = tmp_fished * exp(-Z_seas)
        }

      } # end seas loop
    } # end j loop
  } # end p loop

  # Penultimate ages, catch, and advance to spawning season
  tmp_unfished = Nspr[1,,n_ages - 1]
  tmp_fished = Nspr[2,, n_ages - 1]
  tmp_fished_caa = tmp_fished

  # Catch-at-age for penultimate age
  for(p in 1:n_pop) {
    tmp_fished_caa_p = tmp_fished[p]
    for (seas in 1:n_seas) {
      F_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages-1,] * ret_sel[p, seas, n_ages-1,]) +
        sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages-1,] * (1 - ret_sel[p, seas, n_ages-1,]) * dmr[seas,])
      landed_F_seas = sum(F_fract_flt[seas, landed_fleets] * Fmsy * fish_sel[p,seas,n_ages - 1, landed_fleets] * ret_sel[p,seas,n_ages - 1,landed_fleets]) # landed F
      Z_seas = natmort[p,seas, n_ages-1] * seasdur[seas] + F_seas
      CAA[p, seas, n_ages-1] = tmp_fished_caa_p * (landed_F_seas / Z_seas) * (1 - exp(-Z_seas))
      tmp_fished_caa_p = tmp_fished_caa_p * exp(-Z_seas)
    }
  }

  # Spawning season stuff
  if (spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {
        # Exponential mortality
        M_seas = natmort[p,seas,n_ages - 1] * seasdur[seas]
        tmp_unfished[p] = tmp_unfished[p] * exp(-M_seas)
        tmp_fished[p] = tmp_fished[p] * exp(-(sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages-1,] * ret_sel[p, seas, n_ages-1, ]) +
                                              sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages-1,]  * (1 - ret_sel[p, seas, n_ages-1, ]) * dmr[seas,]) + M_seas))
      }
    }
  }

  # Get spawning biomass after mortality decrement
  for(p in 1:n_pop) {
    SB_age[1, p, n_ages - 1] = tmp_unfished[p] * WAA[p,spawn_seas, n_ages - 1] * MatAA[p,spawn_seas, n_ages - 1] *
      exp(-t_spawn * natmort[p,spawn_seas,n_ages - 1] * seasdur[spawn_seas])
    SB_age[2, p, n_ages - 1] = tmp_fished[p] * WAA[p,spawn_seas, n_ages - 1] * MatAA[p,spawn_seas, n_ages - 1] *
      exp(-t_spawn * (natmort[p,spawn_seas,n_ages - 1] * seasdur[spawn_seas] +
                        sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[p, spawn_seas, n_ages-1,] * ret_sel[p, spawn_seas, n_ages-1, ]) +
                        sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[p, spawn_seas, n_ages-1,] * (1 - ret_sel[p, spawn_seas, n_ages-1, ]) * dmr[spawn_seas, ]) ))
  }

  # Plus group (scalar, no movement)
  for(p in 1:n_pop) {
    F_annual_penult = 0; F_annual_plus = 0
    for(seas in 1:n_seas) {
      F_annual_penult = F_annual_penult + sum(F_fract_flt[seas,] * Fmsy * fish_sel[p,seas,n_ages-1,] * ret_sel[p,seas,n_ages-1,]) +
        sum(F_fract_flt[seas, ] * Fmsy * fish_sel[p,seas,n_ages - 1, ] * (1 - ret_sel[p,seas,n_ages - 1,]) * dmr[seas,])
      F_annual_plus = F_annual_plus + sum(F_fract_flt[seas,] * Fmsy * fish_sel[p,seas,n_ages,] * ret_sel[p,seas,n_ages,]) +
        sum(F_fract_flt[seas, ] * Fmsy * fish_sel[p,seas,n_ages, ] * (1 - ret_sel[p,seas,n_ages,]) * dmr[seas,])
    }
    M_penult = natmort_annual[p, n_ages - 1]
    M_plus = natmort_annual[p, n_ages]
    Nspr[1, p, n_ages] = Nspr[1, p, n_ages-1] * exp(-M_penult) / (1 - exp(-M_plus))
    Nspr[2, p, n_ages] = Nspr[2, p, n_ages-1] * exp(-(M_penult + F_annual_penult)) / (1 - exp(-(M_plus + F_annual_plus)))
  }

  # Plus group catch
  tmp_unfished = Nspr[1,, n_ages]
  tmp_fished = Nspr[2,, n_ages]
  tmp_fished_caa = tmp_fished

  # Catch-at-age for plus group
  for(p in 1:n_pop) {
    tmp_fished_caa_p = Nspr[2, p, n_ages]
    for (seas in 1:n_seas) {
      F_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages,] * ret_sel[p, seas, n_ages,]) +
        sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages,] * (1 - ret_sel[p, seas, n_ages,]) * dmr[seas,])
      landed_F_seas = sum(F_fract_flt[seas, landed_fleets] * Fmsy * fish_sel[p,seas,n_ages, landed_fleets] * ret_sel[p,seas,n_ages,landed_fleets]) # landed F
      Z_seas = natmort[p,seas, n_ages] * seasdur[seas] + F_seas
      CAA[p, seas, n_ages] = tmp_fished_caa_p * (landed_F_seas / Z_seas) * (1 - exp(-Z_seas))
      tmp_fished_caa_p = tmp_fished_caa_p * exp(-Z_seas)
    }
  }

  # Advance plus group to spawning season
  if (spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {
        tmp_unfished[p] = tmp_unfished[p] * exp(-natmort[p,seas,n_ages] * seasdur[seas])
        tmp_fished[p] = tmp_fished[p] * exp(-(natmort[p,seas,n_ages] * seasdur[seas] +
                                                sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages,] * ret_sel[p, seas, n_ages,]) +
                                                sum(F_fract_flt[seas,] * Fmsy * fish_sel[p, seas, n_ages,] * (1 - ret_sel[p, seas, n_ages,]) * dmr[seas,])
                                              ))
      }
    }
  }

  for(p in 1:n_pop) {
    SB_age[1, p, n_ages] = tmp_unfished[p] * WAA[p,spawn_seas, n_ages] * MatAA[p,spawn_seas, n_ages] *
      exp(-t_spawn * natmort[p,spawn_seas,n_ages] * seasdur[spawn_seas])
    SB_age[2, p, n_ages] = tmp_fished[p] * WAA[p,spawn_seas, n_ages] * MatAA[p,spawn_seas, n_ages] *
      exp(-t_spawn * (natmort[p,spawn_seas,n_ages] * seasdur[spawn_seas] +
                        sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[p, spawn_seas, n_ages,] * ret_sel[p, spawn_seas, n_ages, ]) +
                        sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[p, spawn_seas, n_ages,] * (1 - ret_sel[p, spawn_seas, n_ages, ]) * dmr[spawn_seas, ])
                      ))
  }

  # Get effective SB after straying
  effective_SB = array(0, n_pop)
  effective_SB0 = array(0, n_pop)
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

  # Equilibrium recruitment
  Req = equil_rec_phi(effective_SB, effective_SB0, R0, h, rec_model)

  # Yield calculations
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

# Global MSY ----------------------------------------------------------------

#' Compute Beverton-Holt Fmsy for a spatially explicit model
#'
#' Calculates \eqn{F_{MSY}} by maximizing equilibrium yield under a
#' Beverton-Holt stock-recruit relationship. Yield is computed from
#' spawning biomass per recruit (\eqn{\phi_F}), the BH equilibrium
#' recruitment formula, and catch-at-age integrated across all regions,
#' seasons, and movement transitions. Yield includes only landings from
#' fleets where \code{is_discard_fleet == 0}; discard-only fleets
#' contribute to total mortality but not to the yield being maximized.
#'
#' Supports multi-region, single-population models with seasonal movement.
#' Straying is not included here (use \code{single_region_Fmsy} for
#' multi-population non-spatial models).
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
#' Landed yield used in the objective function is computed via the Baranov
#' catch equation using only the landed fraction of fishing mortality
#' (excluding fleets where \code{is_discard_fleet == 1}). The discard
#' fleet's F remains in the Z denominator, so the partitioning correctly
#' accounts for competition between landing and discard mortality sources.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_Fmsy}}{Log-scale trial \eqn{F_{MSY}}.}
#'   }
#'
#' @param data Named list of RTMB data. Must contain:
#'   \describe{
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
#'       \code{[1, n_regions, n_seas, n_ages, n_fish_fleets]}.
#'       Fishery selectivity at age for females.}
#'
#'     \item{\code{ret_sel}}{Numeric array
#'       \code{[1, n_regions, n_seas, n_ages, n_fish_fleets]}.
#'       Retention selectivity (fraction of selected fish retained).}
#'
#'     \item{\code{dmr}}{Numeric array
#'       \code{[n_regions, n_seas, n_fish_fleets]}. Discard mortality rate
#'       (fraction of discarded fish that die).}
#'
#'     \item{\code{natmort}}{Numeric array
#'       \code{[n_regions, n_ages]}. Female natural mortality at age.}
#'
#'     \item{\code{WAA}}{Numeric array
#'       \code{[n_regions, n_seas, n_ages]}. Female weight at age.}
#'
#'     \item{\code{MatAA}}{Numeric array
#'       \code{[n_regions, n_seas, n_ages]}. Maturity at age.}
#'
#'     \item{\code{Movement}}{Numeric array
#'       \code{[n_regions, n_regions, n_seas, n_ages]}. Seasonal movement
#'       transition matrices.}
#'
#'     \item{\code{rec_region_prop}}{Numeric vector
#'       \code{[n_regions]}. Proportion of recruitment entering each region.}
#'
#'     \item{\code{sex_ratio_f}}{Numeric vector
#'       \code{[n_regions]}. Female sex ratio at recruitment.}
#'
#'     \item{\code{rec_seas_prop}}{Numeric vector
#'       \code{[n_seas]}. Seasonal recruitment proportions.}
#'
#'     \item{\code{h}}{Numeric. Beverton-Holt steepness.}
#'
#'     \item{\code{R0}}{Numeric. Unfished equilibrium recruitment.}
#'
#'     \item{\code{is_discard_fleet}}{Integer vector
#'       \code{[n_fish_fleets]}. Indicator for fleets whose catch is excluded
#'       from landed yield (0 = landing fleet, 1 = discard-only fleet). These
#'       fleets still contribute to total fishing mortality \code{Z} and
#'       affect population dynamics and spawning biomass.}
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield (minimized to
#'   find \eqn{F_{MSY}}).
#'
#' @keywords internal
#' @import RTMB
global_Fmsy <- function(pars,
                           data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # getAll defines rec_model whenever the caller supplied it, so hand-built data lists predating
  # the Ricker fall back here rather than being assigned before getAll, which would collide
  if(!exists("rec_model", inherits = FALSE)) rec_model <- 1

  # As with rec_model, expm_nsub is only present on data lists built since the
  # implicit matrix exponential option, so fall back to the exact exponential.
  if(!exists("expm_nsub", inherits = FALSE)) expm_nsub <- 0

  # rec_model 0 is mean recruitment, which has no curve to maximize yield over. the helpers below
  # branch on Ricker against everything else, so a mean recruitment fit would report Beverton-Holt
  if(rec_model == 0) stop("rec_model = 0 (mean recruitment) has no stock-recruit curve, so Fmsy is undefined. ",
                          "MSY reference points need rec_model 1 (Beverton-Holt) or 2 (Ricker); use SPR reference points instead.",
                          call. = FALSE)

  # exponentitate reference points to "estimate"
  Fmsy = exp(log_Fmsy)

  # Helper to get retained F by region for a given (age, season)
  ret_F_by_region <- function(age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[1, , seas, age, f] * ret_sel[1 ,, seas, age, f]
    }
    Fr
  }

  # Helper to get discarded F by region for a given (age, season)
  disc_F_by_region <- function(age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[1, , seas, age, f] * (1 - ret_sel[1 ,, seas, age, f]) * dmr[,seas,f]
    }
    Fr
  }

  # get landed F
  landed_F_by_region <- function(age, seas, is_discard_fleet) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      if (is_discard_fleet[f] == 0)  Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[1, , seas, age, f] * ret_sel[1,,seas, age, f]
    }
    Fr
  }

  SB_age = Nspr = array(0, dim = c(2, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_seas, n_ages)) # catch at age

  # Set up the initial recruits
  Nspr[1,,1] = rec_region_prop * sex_ratio_f * rec_seas_prop[1]
  Nspr[2,,1] = rec_region_prop * sex_ratio_f * rec_seas_prop[1]

  ## Loop through ages
  for (j in 2:(n_ages - 1)) {
    for (seas in 1:n_seas) {

      # get mortality
      ret_F_seas = ret_F_by_region(j - 1, seas)
      disc_F_seas = disc_F_by_region(j - 1, seas)
      landed_F_seas = landed_F_by_region(j - 1, seas, is_discard_fleet)
      M_a_seas = natmort[,seas,j - 1] * seasdur[seas]
      Z_seas = M_a_seas + ret_F_seas + disc_F_seas

      # extract quantities
      tmp_unfished = Nspr[1,, j - 1]
      tmp_fished   = Nspr[2,, j - 1]

      # add in seasonal recruits, apportioned across regions by rec_region_prop to match the initial
      # seeding; dropping it would spread seasons after the first uniformly across regions
      if(seas > 1 && j - 1 == 1) {
        add_rec = rec_seas_prop[seas] * sex_ratio_f * rec_region_prop
        tmp_unfished = tmp_unfished + add_rec
        tmp_fished   = tmp_fished   + add_rec
      }

      ## movement operators for this age. ages that do not move get an identity transition and a
      ## zero generator, which leaves survival unchanged under every move_timing
      if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
        Mv = Movement[,,seas, j - 1]
        Qv = Mrate[,,seas, j - 1]
      } else {
        Mv = diag(n_regions)
        Qv = matrix(0, n_regions, n_regions)
      }

      ## Spawning biomass
      if (seas == spawn_seas) {
        # Propagate to the spawning point; movement and t_spawn mortality applied together
        tmp_unfished_spawn = spawn_state(tmp_unfished, Mv, M_a_seas, Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
        tmp_fished_spawn   = spawn_state(tmp_fished,   Mv, Z_seas,   Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
        SB_age[1,, j - 1] = tmp_unfished_spawn * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1]
        SB_age[2,, j - 1] = tmp_fished_spawn * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1]
      }

      # catch-at-age, taken where the fish are during the season: post-movement under move_timing 0,
      # pre-movement under 1. under continuous movement the spatial Baranov integral is used instead
      CAA[,seas, j - 1] = catch_at_age(tmp_fished, Mv, Z_seas, Qv, seasdur[seas],
                                       landed_F_seas, move_timing, expm_nsub = expm_nsub)

      ## Movement, mortality and ageing
      adv_unfished = advance_seas(tmp_unfished, Mv, M_a_seas, Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
      adv_fished   = advance_seas(tmp_fished,   Mv, Z_seas,   Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)

      if (seas < n_seas) { # Within season mortality
        Nspr[1,, j - 1] = adv_unfished
        Nspr[2,, j - 1] = adv_fished
      } else {
        # Ageing
        Nspr[1,, j] = adv_unfished
        Nspr[2,, j] = adv_fished
      }
    }
  }

  # Age n_ages-1 is now at start of year after the loop
  tmp_unfished = Nspr[1,,n_ages-1]
  tmp_fished = Nspr[2,,n_ages-1]

  # catch-at-age for the penultimate age. tmp_caa holds the survivors season to season, so each
  # season's catch comes from the fish still alive to be caught in it
  tmp_caa = tmp_fished
  for (seas in 1:n_seas) {
    # Compute F and Z for this age/season
    ret_F_seas = ret_F_by_region(n_ages - 1, seas)
    disc_F_seas = disc_F_by_region(n_ages - 1, seas)
    landed_F_seas = landed_F_by_region(n_ages - 1, seas, is_discard_fleet)
    Z_seas = natmort[,seas,n_ages - 1] * seasdur[seas] + ret_F_seas + disc_F_seas
    CAA[,seas, n_ages - 1] = catch_at_age(tmp_caa, Movement[,,seas,n_ages-1], Z_seas,
                                          Mrate[,,seas,n_ages-1], seasdur[seas],
                                          landed_F_seas, move_timing, expm_nsub = expm_nsub)
    tmp_caa = advance_seas(tmp_caa, Movement[,,seas,n_ages-1], Z_seas,
                           Mrate[,,seas,n_ages-1], seasdur[seas], move_timing, expm_nsub = expm_nsub)
  }

  # advance into spawning season
  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      ret_F_seas = ret_F_by_region(n_ages - 1, seas)
      disc_F_seas = disc_F_by_region(n_ages - 1, seas)
      M_seas = natmort[,seas,n_ages-1] * seasdur[seas]
      Z_seas = M_seas + ret_F_seas + disc_F_seas

      # Apply seasonal movement and mortality together, per move_timing
      tmp_unfished = advance_seas(tmp_unfished, Movement[,,seas,n_ages-1], M_seas,
                                  Mrate[,,seas,n_ages-1], seasdur[seas], move_timing, expm_nsub = expm_nsub)
      tmp_fished = advance_seas(tmp_fished, Movement[,,seas,n_ages-1], Z_seas,
                                Mrate[,,seas,n_ages-1], seasdur[seas], move_timing, expm_nsub = expm_nsub)

    } # end seas loop
  }

  ## Penultimate age spawning biomass; t_spawn mortality folded into spawn_state
  Zu_spawn = natmort[,spawn_seas, n_ages - 1] * seasdur[spawn_seas]
  Zf_spawn = Zu_spawn + ret_F_by_region(n_ages - 1, spawn_seas) + disc_F_by_region(n_ages - 1, spawn_seas)
  tmp_unfished = spawn_state(tmp_unfished, Movement[,, spawn_seas, n_ages - 1], Zu_spawn,
                             Mrate[,, spawn_seas, n_ages - 1], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
  tmp_fished   = spawn_state(tmp_fished, Movement[,, spawn_seas, n_ages - 1], Zf_spawn,
                             Mrate[,, spawn_seas, n_ages - 1], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
  SB_age[1,, n_ages - 1] = tmp_unfished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1]
  SB_age[2,, n_ages - 1] = tmp_fished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1]

  # Get fishing mortality for plus group calculations
  F_penult = F_plus = matrix(0, n_regions, n_seas)
  for(seas in 1:n_seas) {
    F_penult[, seas] = ret_F_by_region(n_ages - 1, seas) + disc_F_by_region(n_ages - 1, seas)
    F_plus[, seas] = ret_F_by_region(n_ages, seas) + disc_F_by_region(n_ages, seas)
  }

  ## Plus group analytical solution
  Ts = build_plus_group_T(
    M_penult = array(natmort[,, n_ages - 1], dim = c(n_regions, n_seas)),
    M_plus = array(natmort[,, n_ages], dim = c(n_regions, n_seas)),
    F_penult   = F_penult,
    F_plus     = F_plus,
    Mov_penult = array(Movement[,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
    Mov_plus   = array(Movement[,,, n_ages], dim = c(n_regions, n_regions, n_seas)),
    n_regions  = n_regions,
    n_seas = n_seas,
    seasdur = seasdur,
    Mrate_penult = array(Mrate[,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
    Mrate_plus   = array(Mrate[,,, n_ages], dim = c(n_regions, n_regions, n_seas)),
    move_timing = move_timing,
    expm_nsub = expm_nsub
  )

  pg = solve_plus_group(Ts, Nspr[1,, n_ages - 1], Nspr[2,, n_ages - 1], n_regions)
  Nspr[1,, n_ages] = pg$unfished
  Nspr[2,, n_ages] = pg$fished

  ## Plus group catch, then advance to spawning season
  tmp_unfished = Nspr[1,, n_ages]
  tmp_fished = Nspr[2,, n_ages]

  # Catch-at-age for plus group, taken from the survivors in each season
  tmp_caa = tmp_fished
  for (seas in 1:n_seas) {
    ret_F_seas = ret_F_by_region(n_ages, seas)
    disc_F_seas = disc_F_by_region(n_ages, seas)
    landed_F_seas = landed_F_by_region(n_ages, seas, is_discard_fleet)
    Z_seas = natmort[,seas,n_ages] * seasdur[seas] + ret_F_seas + disc_F_seas
    CAA[,seas, n_ages] = catch_at_age(tmp_caa, Movement[,,seas,n_ages], Z_seas,
                                      Mrate[,,seas,n_ages], seasdur[seas],
                                      landed_F_seas, move_timing, expm_nsub = expm_nsub)
    tmp_caa = advance_seas(tmp_caa, Movement[,,seas,n_ages], Z_seas,
                           Mrate[,,seas,n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
  }


  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # get mortality
      M_seas = natmort[,seas,n_ages] * seasdur[seas]
      ret_F_seas = ret_F_by_region(n_ages, seas)
      disc_F_seas = disc_F_by_region(n_ages, seas)
      Z_seas = ret_F_seas + M_seas + disc_F_seas

      # Apply seasonal movement and mortality together, per move_timing
      tmp_unfished = advance_seas(tmp_unfished, Movement[,,seas,n_ages], M_seas,
                                  Mrate[,,seas,n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
      tmp_fished = advance_seas(tmp_fished, Movement[,,seas,n_ages], Z_seas,
                                Mrate[,,seas,n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)

    } # end seas loop
  }

  ## Plus group spawning biomass; t_spawn mortality folded into spawn_state
  Zu_spawn_plus = natmort[,spawn_seas, n_ages] * seasdur[spawn_seas]
  Zf_spawn_plus = Zu_spawn_plus + ret_F_by_region(n_ages, spawn_seas) + disc_F_by_region(n_ages, spawn_seas)
  tmp_unfished = spawn_state(tmp_unfished, Movement[,, spawn_seas, n_ages], Zu_spawn_plus,
                             Mrate[,, spawn_seas, n_ages], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
  tmp_fished   = spawn_state(tmp_fished, Movement[,, spawn_seas, n_ages], Zf_spawn_plus,
                             Mrate[,, spawn_seas, n_ages], seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
  SB_age[1,, n_ages] = tmp_unfished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages]
  SB_age[2,, n_ages] = tmp_fished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages]

  # Get spawning biomass per recruit to get spawning potential ratio
  SBPR_0 = sum(SB_age[1,,])
  SBPR_F = sum(SB_age[2,,])
  SPR = SBPR_F / SBPR_0

  # Get equilibrium recruitment
  Req = equil_rec_phi(SBPR_F, SBPR_0, R0, h, rec_model)

  # Get yield
  Yield = sum(CAA * WAA) * Req
  Yield_r = rowSums(CAA * WAA) * Req

  # compute objective function to get Fmsy
  obj_fun = -Yield

  RTMB::REPORT(SB_age)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SBPR_0)
  RTMB::REPORT(SBPR_F)
  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Yield)
  RTMB::REPORT(Yield_r)
  RTMB::REPORT(Req)
  RTMB::REPORT(SPR)

  return(obj_fun)
}

# Local MSY, Single Population ----------------------------------------------

#' Compute region-specific Beverton-Holt Fmsy for a spatially explicit
#' single-population model
#'
#' Computes the vector of regional \eqn{F_{MSY}} values that jointly maximize
#' total equilibrium yield across all regions under a spatially explicit
#' Beverton-Holt stock-recruit relationship. Unlike \code{\link{global_Fmsy}},
#' which constrains all regions to share a single fishing mortality, this
#' function allows each region to have its own optimal \eqn{F}.
#'
#' Cohorts originating in each region are tracked separately through seasonal
#' movement, mortality, and ageing using an \code{[origin, destination]}
#' per-recruit accounting framework. Spawning biomass per recruit is accumulated
#' by origin and destination region, and the plus group is solved analytically
#' using \code{\link{build_plus_group_T}} and \code{\link{solve_plus_group}}.
#'
#' Equilibrium recruitment by origin region \eqn{R_{eq,o}} is solved using a
#' Newton-Raphson algorithm applied to the fixed-point condition that
#' recruitment produced at each destination region (via the BH relationship
#' applied to effective SSB) equals the recruitment attributed to that origin.
#' The Jacobian is derived analytically using the quotient rule and the chain
#' rule through the spatial redistribution of spawning biomass.
#'
#' Yield is computed using only the landed fraction of fishing mortality,
#' excluding fleets flagged as discard-only via \code{is_discard_fleet}.
#' Discard-only fleets still contribute to total mortality \code{Z} and
#' affect population dynamics and spawning biomass.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_Fmsy}}{Numeric vector \code{[n_regions]}. Log-scale
#'       trial \eqn{F_{MSY}} values, one per region.}
#'   }
#'
#' @param data Named list of RTMB data. Must contain all spatial fields required
#'   by \code{\link{global_SPR}} (excluding \code{SPR_x}, \code{stray_rate},
#'   and \code{natal_region}) plus:
#'   \describe{
#'     \item{\code{h}}{Numeric vector \code{[n_regions]}. Beverton-Holt
#'       steepness by region.}
#'     \item{\code{R0}}{Numeric scalar. Total unfished equilibrium recruitment.}
#'     \item{\code{rec_region_prop}}{Numeric vector \code{[n_regions]}.
#'       Proportion of annual recruitment entering each region.}
#'     \item{\code{newton_steps}}{Integer. Number of Newton-Raphson iterations
#'       used to solve for equilibrium recruitment by origin region.}
#'     \item{\code{is_discard_fleet}}{Integer vector \code{[n_fish_fleets]}.
#'       Indicator for fleets whose catch is excluded from landed yield
#'       (0 = landing fleet, 1 = discard-only fleet). These fleets still
#'       contribute to total fishing mortality \code{Z} and affect population
#'       dynamics and spawning biomass.}
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield across all regions.
#'   This is minimized to obtain the vector of regional \eqn{F_{MSY}} values.
#'
#' @details
#' Fishing mortality is decomposed into retained and discarded components:
#'
#' \itemize{
#'   \item Retained fishing mortality:
#'     \deqn{F^{\mathrm{ret}}_{r,a,s,f} = F_{MSY,r} \, F_{\mathrm{fract},r,s,f} \,
#'           \mathrm{sel}_{r,a,s,f} \, \mathrm{ret}_{r,a,s,f}}
#'
#'   \item Discard fishing mortality (dead discards only):
#'     \deqn{F^{\mathrm{disc}}_{r,a,s,f} =
#'           F_{MSY,r} \, F_{\mathrm{fract},r,s,f} \,
#'           \mathrm{sel}_{r,a,s,f} \, (1 - \mathrm{ret}_{r,a,s,f}) \,
#'           \mathrm{dmr}_{r,s,f}}
#'
#'   \item Total instantaneous mortality:
#'     \deqn{Z_{r,a,s} = M_{r,a} \, \mathrm{seasdur}_s +
#'           F^{\mathrm{ret}}_{r,a,s} + F^{\mathrm{disc}}_{r,a,s}}
#' }
#'
#' Landed yield used in the objective function excludes catch from fleets
#' where \code{is_discard_fleet == 1}. The Baranov catch equation partitions
#' landed F out of total Z, so the discard fleet's contribution to mortality
#' is properly accounted for in the denominator.
#'
#' Seasonal movement is applied using the \code{Movement[origin, dest, seas, age]}
#' array. Recruitment may move immediately or only after age-1 depending on
#' \code{do_recruits_move}. Spawning biomass is accumulated at
#' \code{spawn_seas} with fractional mortality \code{t_spawn}.
#'
#' The plus group is solved analytically using the transition matrices produced
#' by \code{\link{build_plus_group_T}} and the solver \code{\link{solve_plus_group}}.
#'
#' @keywords internal
#' @import RTMB
local_Fmsy_sglpop <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # getAll defines rec_model whenever the caller supplied it, so hand-built data lists predating
  # the Ricker fall back here rather than being assigned before getAll, which would collide
  if(!exists("rec_model", inherits = FALSE)) rec_model <- 1

  # As with rec_model, expm_nsub is only present on data lists built since the
  # implicit matrix exponential option, so fall back to the exact exponential.
  if(!exists("expm_nsub", inherits = FALSE)) expm_nsub <- 0

  # rec_model 0 is mean recruitment, which has no curve to maximize yield over. the helpers below
  # branch on Ricker against everything else, so a mean recruitment fit would report Beverton-Holt
  if(rec_model == 0) stop("rec_model = 0 (mean recruitment) has no stock-recruit curve, so Fmsy is undefined. ",
                          "MSY reference points need rec_model 1 (Beverton-Holt) or 2 (Ricker); use SPR reference points instead.",
                          call. = FALSE)

  # exponentitate reference points
  Fmsy = exp(log_Fmsy)

  # helper function to get f by region for a given age and season (retained)
  ret_F_by_region <- function(age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[1, , seas, age, f] * ret_sel[1,,seas, age, f]
    }
    Fr
  }

  # helper function to get f by region (discarded)
  disc_F_by_region <- function(age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[1, , seas, age, f] * (1 - ret_sel[1,,seas, age, f]) * dmr[,seas, f]
    }
    Fr
  }

  # get landed F
  landed_F_by_region <- function(age, seas, is_discard_fleet) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      if (is_discard_fleet[f] == 0)  Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[1, , seas, age, f] * ret_sel[1,,seas, age, f]
    }
    Fr
  }

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_regions, n_seas, n_ages)) # catch at age
  Yield_r = array(0, dim = n_regions) # yield by region
  SB_unfished_mat = matrix(0, n_regions, n_regions)  # unfished spawning biomass per recruit
  SB_fished_mat = matrix(0, n_regions, n_regions) # fished spawning biomass per recruit

  # Set up the initial recruits (1 recruit per origin area on the diagonal)
  for(o in 1:n_regions) {
    for(d in 1:n_regions) {
      if(o == d) Nspr[1,o,d,1] = Nspr[2,o,d,1] = sex_ratio_f[o] * rec_seas_prop[1]
      else Nspr[1,o,d,1] = Nspr[2,o,d,1] = 0
    } # end d loop
  } # end o loop

  ## Loop through ages
  for(j in 2:(n_ages - 1)) {
    for(seas in 1:n_seas) {
      for(o in 1:n_regions) {

        # get mortality
        ret_F_seas = ret_F_by_region(j - 1, seas)
        disc_F_seas = disc_F_by_region(j - 1, seas)
        landed_F_seas = landed_F_by_region(j - 1, seas, is_discard_fleet)
        M_seas = natmort[,seas, j - 1] * seasdur[seas]
        Z_seas = ret_F_seas + M_seas + disc_F_seas

        # extract out quantities
        tmp_unfished = Nspr[1, o,, j - 1]; tmp_fished = Nspr[2, o,, j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          tmp_unfished[o] = tmp_unfished[o] + rec_seas_prop[seas] * sex_ratio_f[o]
          tmp_fished[o]   = tmp_fished[o]   + rec_seas_prop[seas] * sex_ratio_f[o]
        }

        ## Movement operators for this age. Ages that do not move get an identity
        ## transition and a zero generator, leaving survival unchanged under every timing.
        if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
          Mv = Movement[,,seas, j - 1]
          Qv = Mrate[,,seas, j - 1]
        } else {
          Mv = diag(n_regions)
          Qv = matrix(0, n_regions, n_regions)
        }

        ## Spawning biomass; movement and t_spawn mortality applied together
        if(seas == spawn_seas) {
          tmp_unfished_spawn = spawn_state(tmp_unfished, Mv, M_seas, Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
          tmp_fished_spawn   = spawn_state(tmp_fished,   Mv, Z_seas, Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
          SB_age[1, o,, j - 1] = tmp_unfished_spawn * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1]
          SB_age[2, o,, j - 1] = tmp_fished_spawn   * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1]
        }

        # Catch-at-age, taken where the fish are during the season per move_timing
        CAA[o,, seas, j - 1] = catch_at_age(tmp_fished, Mv, Z_seas, Qv, seasdur[seas],
                                            landed_F_seas, move_timing, expm_nsub = expm_nsub)

        ## Movement, mortality and ageing
        adv_unfished = advance_seas(tmp_unfished, Mv, M_seas, Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
        adv_fished   = advance_seas(tmp_fished,   Mv, Z_seas, Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)

        if(seas < n_seas) { # Within season mortality
          Nspr[1, o,, j - 1] = adv_unfished
          Nspr[2, o,, j - 1] = adv_fished
        } else {
          # Ageing
          Nspr[1, o,, j] = adv_unfished
          Nspr[2, o,, j] = adv_fished
        }

      } # end o loop
    } # end seas loop
  } # end j loop

  # Advance penultimate age into spawning season
  tmp_unfished = array(Nspr[1,,, n_ages - 1], dim = c(n_regions, n_regions))
  tmp_fished   = array(Nspr[2,,, n_ages - 1], dim = c(n_regions, n_regions))

  for(o in 1:n_regions) {
    # catch-at-age for the penultimate age. tmp_caa holds this origin cohort's survivors season to
    # season, so each season's catch comes from the fish still alive to be caught in it
    tmp_caa = tmp_fished[o,]
    for(seas in 1:n_seas) {
      ret_F_seas = ret_F_by_region(n_ages - 1, seas)
      disc_F_seas = disc_F_by_region(n_ages - 1, seas)
      landed_F_seas = landed_F_by_region(n_ages - 1, seas, is_discard_fleet)
      Z_seas = natmort[,seas, n_ages - 1] * seasdur[seas] + ret_F_seas + disc_F_seas
      CAA[o,, seas, n_ages - 1] = catch_at_age(tmp_caa, Movement[,,seas, n_ages - 1], Z_seas,
                                               Mrate[,,seas, n_ages - 1], seasdur[seas],
                                               landed_F_seas, move_timing, expm_nsub = expm_nsub)
      tmp_caa = advance_seas(tmp_caa, Movement[,,seas, n_ages - 1], Z_seas,
                             Mrate[,,seas, n_ages - 1], seasdur[seas], move_timing, expm_nsub = expm_nsub)
    }
  }

  if(spawn_seas > 1) {
    for(o in 1:n_regions) {
      for(seas in 1:(spawn_seas - 1)) {

        # get mortality
        M_seas = natmort[,seas, n_ages - 1] * seasdur[seas]
        ret_F_seas = ret_F_by_region(n_ages - 1, seas)
        disc_F_seas = disc_F_by_region(n_ages - 1, seas)

        # Apply seasonal movement and mortality together, per move_timing
        tmp_unfished[o,] = advance_seas(tmp_unfished[o,], Movement[,,seas, n_ages - 1], M_seas,
                                        Mrate[,,seas, n_ages - 1], seasdur[seas], move_timing, expm_nsub = expm_nsub)
        tmp_fished[o,]   = advance_seas(tmp_fished[o,], Movement[,,seas, n_ages - 1],
                                        M_seas + ret_F_seas + disc_F_seas,
                                        Mrate[,,seas, n_ages - 1], seasdur[seas], move_timing, expm_nsub = expm_nsub)

      } # end seas loop
    } # end o loop
  }

  ## Advance penultimate age into spawning biomass / season; t_spawn folded in here
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn   = tmp_fished
  Zu_spawn = natmort[,spawn_seas, n_ages - 1] * seasdur[spawn_seas]
  Zf_spawn = Zu_spawn + ret_F_by_region(n_ages - 1, spawn_seas) + disc_F_by_region(n_ages - 1, spawn_seas)
  for(o in 1:n_regions) {
    tmp_unfished_spawn[o,] = spawn_state(tmp_unfished_spawn[o,], Movement[,,spawn_seas, n_ages - 1],
                                         Zu_spawn, Mrate[,,spawn_seas, n_ages - 1],
                                         seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
    tmp_fished_spawn[o,]   = spawn_state(tmp_fished_spawn[o,], Movement[,,spawn_seas, n_ages - 1],
                                         Zf_spawn, Mrate[,,spawn_seas, n_ages - 1],
                                         seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
  }

  SB_age[1,,, n_ages - 1] = tmp_unfished_spawn * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1]
  SB_age[2,,, n_ages - 1] = tmp_fished_spawn * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1]

  ## Plus group analytical solution
  # Get fishing mortality for penultimate and plus ages
  F_penult = F_plus = matrix(0, n_regions, n_seas)
  for(seas in 1:n_seas) {
    F_penult[, seas] = ret_F_by_region(n_ages - 1, seas) + disc_F_by_region(n_ages - 1, seas)
    F_plus[, seas]   = ret_F_by_region(n_ages, seas) + disc_F_by_region(n_ages, seas)
  }

  # Build transition matrices
  Ts <- build_plus_group_T(
    M_penult = array(natmort[,, n_ages - 1], dim = c(n_regions, n_seas)),
    M_plus = array(natmort[,, n_ages], dim = c(n_regions, n_seas)),
    F_penult   = F_penult,
    F_plus     = F_plus,
    Mov_penult = array(Movement[,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
    Mov_plus   = array(Movement[,,, n_ages],     dim = c(n_regions, n_regions, n_seas)),
    n_regions  = n_regions,
    n_seas     = n_seas,
    seasdur    = seasdur,
    Mrate_penult = array(Mrate[,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
    Mrate_plus   = array(Mrate[,,, n_ages],     dim = c(n_regions, n_regions, n_seas)),
    move_timing  = move_timing,
    expm_nsub = expm_nsub)

  for(o in 1:n_regions) {
    # Get and input plus group
    pg <- solve_plus_group(Ts, Nspr[1, o,, n_ages - 1], Nspr[2, o,, n_ages - 1], n_regions)
    Nspr[1, o,, n_ages] <- pg$unfished
    Nspr[2, o,, n_ages] <- pg$fished
  }

  ## Plus group: catch, then advance to spawning season
  tmp_unfished = array(Nspr[1,,, n_ages], dim = c(n_regions, n_regions))
  tmp_fished   = array(Nspr[2,,, n_ages], dim = c(n_regions, n_regions))

  for(o in 1:n_regions) {
    # Catch-at-age for plus group, taken from this origin cohort's survivors in each season
    tmp_caa = tmp_fished[o,]
    for(seas in 1:n_seas) {
      ret_F_seas = ret_F_by_region(n_ages, seas)
      disc_F_seas = disc_F_by_region(n_ages, seas)
      landed_F_seas = landed_F_by_region(n_ages, seas, is_discard_fleet)
      Z_seas = natmort[,seas, n_ages] * seasdur[seas] + ret_F_seas + disc_F_seas
      CAA[o,, seas, n_ages] = catch_at_age(tmp_caa, Movement[,,seas, n_ages], Z_seas,
                                           Mrate[,,seas, n_ages], seasdur[seas],
                                           landed_F_seas, move_timing, expm_nsub = expm_nsub)
      tmp_caa = advance_seas(tmp_caa, Movement[,,seas, n_ages], Z_seas,
                             Mrate[,,seas, n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
    }
  }

  if(spawn_seas > 1) {
    for(o in 1:n_regions) {
      for(seas in 1:(spawn_seas - 1)) {

        # get mortality
        M_seas = natmort[,seas, n_ages] * seasdur[seas]
        F_seas = F_plus[, seas]

        # Apply seasonal movement and mortality together, per move_timing
        tmp_unfished[o,] = advance_seas(tmp_unfished[o,], Movement[,,seas, n_ages], M_seas,
                                        Mrate[,,seas, n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
        tmp_fished[o,]   = advance_seas(tmp_fished[o,], Movement[,,seas, n_ages], M_seas + F_seas,
                                        Mrate[,,seas, n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)

      } # end seas loop
    } # end o loop
  }

  ## Plus group spawning biomass; t_spawn mortality folded into spawn_state
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn   = tmp_fished
  Zu_spawn_plus = natmort[,spawn_seas, n_ages] * seasdur[spawn_seas]
  Zf_spawn_plus = Zu_spawn_plus + F_plus[, spawn_seas]
  for(o in 1:n_regions) {
    tmp_unfished_spawn[o,] = spawn_state(tmp_unfished_spawn[o,], Movement[,,spawn_seas, n_ages],
                                         Zu_spawn_plus, Mrate[,,spawn_seas, n_ages],
                                         seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
    tmp_fished_spawn[o,]   = spawn_state(tmp_fished_spawn[o,], Movement[,,spawn_seas, n_ages],
                                         Zf_spawn_plus, Mrate[,,spawn_seas, n_ages],
                                         seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
  }

  SB_age[1,,, n_ages] = tmp_unfished_spawn * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages]
  SB_age[2,,, n_ages] = tmp_fished_spawn   * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages]

  # Determine equilibrium recruitment for destination region
  # parse out and compute unfished and fished spawning biomass per recruit
  for(o in 1:n_regions) {
    for(d in 1:n_regions) {
      SB_unfished_mat[o, d] = sum(SB_age[1, o, d, ])  # unfished
      SB_fished_mat[o, d]   = sum(SB_age[2, o, d, ])  # fished at Fmsy
    } # end d loop
  } # end o loop

  # Unfished spawning biomass and recruitment by destination region: the two
  # quantities both stock-recruit curves are parameterized by.
  S0_d = rep(0, n_regions)
  for(d in 1:n_regions) S0_d[d] = sum(SB_unfished_mat[, d] * R0 * rec_region_prop)
  R0_d = R0 * rec_region_prop

  # define initial guess to solve for equilibrium recruitment from origin region
  Req_o = R0 * rec_region_prop

  for(nit in 1:newton_steps) {
    # compute equilibrium spawning biomass (SSBR * Req) in destination region
    x_vec = as.numeric(t(SB_fished_mat) %*% Req_o)  # function of equilibrium recruitment in origin region (what we are solving for)
    g_vec = equil_rec_ssb(x_vec, S0_d, R0_d, h, rec_model) # equilibrium recruitment in destination region

    # define root and define Jacobian
    iter_vec = Req_o - g_vec # find values of origin recruitment that are consistent w/ destination recruitment such that pop'n is in equilibrium

    # construct Jacobian for root
    # we need J = df (iter_vec)/dReq = dReq/dReq (or I) - dg/dReq
    # we basically want to know dg / dReq (how does destination equil rec change as origin equil rec change)
    # dg / dReq = (dg / dxk) * (dxk / dReq)
    # to get (dg / dxk), use quotient rule of (BH recruitment)
    # note that dxk / dReq S_2mat * Req = S_2mat
    dg_dxk = equil_rec_ssb_deriv(x_vec, S0_d, R0_d, h, rec_model)
    dg_dReq = matrix(0, n_regions, n_regions)
    for(d in 1:n_regions) dg_dReq[d, ] = dg_dxk[d] * SB_fished_mat[, d] # now compute to see how destination equilibrium rec changes, as origin equil rec changes

    # compute jacobian
    J = diag(1, n_regions) - dg_dReq
    delta = solve(J, as.vector(iter_vec)) # get step to move towards solution
    Req_o = Req_o - delta # newton raphson update
  }

  # get destination region yield
  for(d in 1:n_regions) {
    tmp = 0 # define temp variable
    for(seas in 1:n_seas) for(o in 1:n_regions) tmp = tmp + sum(CAA[o, d, seas, ] * WAA[d, seas, ]) * Req_o[o] # get yield to destination
    Yield_r[d] = tmp
  } # end d loop

  # Get spawning biomass per recruit summed across origins weighted by rec_region_prop
  SB  <- matrix(0, n_pop, n_regions)
  SB0 <- matrix(0, n_pop, n_regions)
  for(r in 1:n_regions) {
    SB0[1, r] <- sum(SB_unfished_mat[, r] * rec_region_prop)
    SB[1, r]  <- sum(SB_fished_mat[, r]   * rec_region_prop)
  }

  # maximize total yield
  Yield_total = sum(Yield_r)
  obj_fun = -Yield_total

  sum_SB_unfished_mat = sum(SB_unfished_mat)

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

# Local MSY, Multiple Populations -------------------------------------------

#' Compute local Beverton-Holt Fmsy for a spatially explicit multi-population model
#'
#' Multi-population extension of \code{\link{local_Fmsy_sglpop}}. Estimates a
#' vector of region-specific \eqn{F_{MSY}} values that jointly maximize total
#' equilibrium yield when multiple populations, each with distinct natal regions,
#' movement schedules, and Beverton-Holt parameters, co-occupy a shared spatial
#' domain.
#'
#' Cohorts are tracked using a per-recruit framework indexed by
#' \code{[population x origin x destination x age x season]}. Recruitment is
#' distributed across seasons (\code{rec_seas_prop}), regions
#' (\code{rec_region_prop}), and sex (\code{sex_ratio_f}). Initial recruits are
#' assigned in the first season at age-1, with additional seasonal recruitment
#' contributions added within the first age class prior to movement and mortality.
#'
#' Movement is applied at each seasonal step using region- and age-specific
#' transition matrices. Fishing mortality is decomposed into retained and
#' discarded components, and total mortality is applied continuously within each
#' season. Catch-at-age is accumulated across fleets, seasons, and regions using
#' only the landed fraction of fishing mortality (excluding fleets flagged as
#' discard-only via \code{is_discard_fleet}), while total mortality \code{Z}
#' includes all fleets.
#'
#' Spawning biomass per recruit (SBPR) is computed at the spawning season after
#' applying movement and partial mortality up to the spawning time
#' (\code{t_spawn}). For single-season models with multiple populations,
#' \code{sgl_seas_spawning_movement} redistributes individuals to natal spawning
#' regions prior to SSB calculation.
#'
#' The plus group is solved analytically using
#' \code{\link{build_plus_group_T}} and \code{\link{solve_plus_group}}, ensuring
#' a consistent equilibrium solution for the terminal age class.
#'
#' Effective spawning biomass at each population's natal region includes
#' contributions from all populations via straying. Stray contributions are
#' scaled by \code{stray_rate} and normalized by
#' \code{n_pop_in_region} to preserve mass balance.
#'
#' Equilibrium recruitment by population is obtained via a Newton-Raphson
#' algorithm that solves the coupled Beverton-Holt system. The Jacobian accounts
#' for cross-population dependence of spawning biomass induced by straying.
#'
#' Total equilibrium yield is computed by integrating catch-at-age over all
#' populations, regions, and seasons, scaled by equilibrium recruitment and
#' origin-region proportions.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_Fmsy}}{Log-scale trial \eqn{F_{MSY}} values, one per
#'       region (length \code{n_regions}).}
#'   }
#' @param data Named list of RTMB data. Must contain all fields required by
#'   \code{\link{global_SPR}} plus:
#'   \describe{
#'     \item{\code{h}}{Numeric array \code{[n_pop, n_regions]}. Beverton-Holt
#'       steepness evaluated at each population's natal region.}
#'     \item{\code{R0}}{Numeric vector \code{[n_pop]}. Unfished equilibrium
#'       recruitment per population.}
#'     \item{\code{stray_rate}}{Numeric vector \code{[n_pop]}. Fraction of
#'       individuals contributing to non-natal spawning regions.}
#'     \item{\code{natal_region}}{Integer vector \code{[n_pop]}. Natal region
#'       index for each population.}
#'     \item{\code{n_pop_in_region}}{Integer vector \code{[n_regions]}. Number
#'       of populations sharing each natal region (used to normalize straying).}
#'     \item{\code{newton_steps}}{Integer. Number of Newton-Raphson iterations
#'       used to solve for equilibrium recruitment.}
#'     \item{\code{is_discard_fleet}}{Integer vector \code{[n_fish_fleets]}.
#'       Indicator for fleets whose catch is excluded from landed yield
#'       (0 = landing fleet, 1 = discard-only fleet). These fleets still
#'       contribute to total fishing mortality \code{Z} and affect population
#'       dynamics and spawning biomass.}
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield across all regions.
#'   This objective is minimized to obtain the vector of regional
#'   \eqn{F_{MSY}} values.
#'
#' @details
#' Recruitment is implemented as a per-recruit process and later scaled by
#' equilibrium recruitment. Seasonal recruitment proportions are applied within
#' the first age class, allowing intra-annual timing of recruitment before
#' movement and mortality are applied.
#'
#' The objective function is the negative of total yield, summed across all
#' populations and regions. Yield includes only landings from fleets where
#' \code{is_discard_fleet == 0}; discard-only fleets contribute to mortality
#' but not to the yield being maximized.
#'
#' @keywords internal
#' @import RTMB
local_Fmsy_multipop <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # getAll defines rec_model whenever the caller supplied it, so hand-built data lists predating
  # the Ricker fall back here rather than being assigned before getAll, which would collide
  if(!exists("rec_model", inherits = FALSE)) rec_model <- 1

  # As with rec_model, expm_nsub is only present on data lists built since the
  # implicit matrix exponential option, so fall back to the exact exponential.
  if(!exists("expm_nsub", inherits = FALSE)) expm_nsub <- 0

  # rec_model 0 is mean recruitment, which has no curve to maximize yield over. the helpers below
  # branch on Ricker against everything else, so a mean recruitment fit would report Beverton-Holt
  if(rec_model == 0) stop("rec_model = 0 (mean recruitment) has no stock-recruit curve, so Fmsy is undefined. ",
                          "MSY reference points need rec_model 1 (Beverton-Holt) or 2 (Ricker); use SPR reference points instead.",
                          call. = FALSE)

  # exponentitate reference points
  Fmsy = exp(log_Fmsy)

  # helper function to get f by region for a given age and season
  ret_F_by_region <- function(p, age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[p, , seas, age, f] * ret_sel[p,, seas, age, f]
    }
    Fr
  }

  disc_F_by_region <- function(p, age, seas) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[p, , seas, age, f] * (1 - ret_sel[p,, seas, age, f]) * dmr[, seas, f]
    }
    Fr
  }

  # get landed F
  landed_F_by_region <- function(p, age, seas, is_discard_fleet) {
    Fr <- numeric(n_regions)
    for (f in seq_len(dim(fish_sel)[length(dim(fish_sel))])) {
      if (is_discard_fleet[f] == 0) {
        Fr <- Fr + F_fract_flt[, seas, f] * Fmsy * fish_sel[p, , seas, age, f] * ret_sel[p,, seas, age, f]
      }
    }
    Fr
  }

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_regions, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_pop, n_regions, n_regions, n_seas, n_ages)) # catch at age
  Yield_pr = array(0, dim = c(n_pop, n_regions)) # yield by region
  SB_unfished_mat = array(0, dim = c(n_pop, n_regions, n_regions))  # unfished spawning biomass per recruit
  SB_fished_mat   = array(0, dim = c(n_pop, n_regions, n_regions))  # fished spawning biomass per recruit

  # Set up the initial recruits (1 recruit per natal region / area)
  for(p in 1:n_pop) {
    for(o in 1:n_regions) {
      Nspr[1,p,o,o,1] = Nspr[2,p,o,o,1] = sex_ratio_f[p,o] * rec_seas_prop[p,1] * rec_region_prop[p,o]
    }
  }

  # Loop through, apply movement first, then decrement recruit
  # Loop through ages, projecting each cohort through the full annual cycle
  for(j in 2:(n_ages-1)){

    # Project age j-1 through all seasons to become age j
    for(seas in 1:n_seas) {

      for(p in 1:n_pop) {
        for(o in 1:n_regions) {

          ret_F_a_seas = ret_F_by_region(p, j -1, seas)
          disc_F_a_seas = disc_F_by_region(p, j -1, seas)
          landed_F_a_seas = landed_F_by_region(p, j - 1, seas, is_discard_fleet)
          Z_a_seas = natmort[p,,seas,j - 1] * seasdur[seas] + ret_F_a_seas + disc_F_a_seas

          # Get temporary values from origin region
          tmp_unfished = Nspr[1,p,o,,j-1]
          tmp_fished = Nspr[2,p,o,,j-1]

          # add in seasonal recruits
          if(seas > 1 && j - 1 == 1) {
            tmp_unfished[o] = tmp_unfished[o] + rec_seas_prop[p,seas] * sex_ratio_f[p,o] * rec_region_prop[p,o]
            tmp_fished[o]   = tmp_fished[o]   + rec_seas_prop[p,seas] * sex_ratio_f[p,o] * rec_region_prop[p,o]
          }

          # Movement operators for this age; non-moving ages get an identity transition
          # and a zero generator, leaving survival unchanged under every move_timing
          if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
            Mv = Movement[p,,,seas,j-1]
            Qv = Mrate[p,,,seas,j-1]
          } else {
            Mv = diag(n_regions)
            Qv = matrix(0, n_regions, n_regions)
          }
          Mu_a_seas = natmort[p,,seas,j-1] * seasdur[seas]

          # Calculate spawning biomass if this is the spawning season
          if(seas == spawn_seas) {

            # Propagate to the spawning point; movement and t_spawn mortality together
            tmp_unfished_spawn = spawn_state(tmp_unfished, Mv, Mu_a_seas, Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)
            tmp_fished_spawn   = spawn_state(tmp_fished, Mv, Z_a_seas, Qv, seasdur[seas], t_spawn, move_timing, expm_nsub = expm_nsub)

            # If single season natal homing population
            if(n_pop > 1 && n_seas == 1) {
              # Get NAA during spawning in single season case
              tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
              tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
            }

            for(d in 1:n_regions) {
              SB_age[1,p,o,d,j-1] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1]
              SB_age[2,p,o,d,j-1] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1]
            }
          }

          # Catch-at-age, taken where the fish are during the season per move_timing
          CAA[p,o,,seas,j-1] = catch_at_age(tmp_fished, Mv, Z_a_seas, Qv, seasdur[seas],
                                            landed_F_a_seas, move_timing, expm_nsub = expm_nsub)

          # Apply movement and mortality together, per move_timing
          adv_unfished = advance_seas(tmp_unfished, Mv, Mu_a_seas, Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
          adv_fished   = advance_seas(tmp_fished, Mv, Z_a_seas, Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)

          if(seas < n_seas) {
            # Within-season mortality, no ageing yet
            Nspr[1,p,o,,j-1] = adv_unfished
            Nspr[2,p,o,,j-1] = adv_fished
          } else {
            # Last season: mortality + ageing
            Nspr[1,p,o,,j] = adv_unfished
            Nspr[2,p,o,,j] = adv_fished
          }

        } # end o loop
      } # end p loop

    } # end seas loop
  } # end j loop

  # Advance penultimate age into spawning season
  for(p in 1:n_pop) {

    tmp_unfished = array(Nspr[1, p,,, n_ages - 1], dim = c(n_regions, n_regions))
    tmp_fished   = array(Nspr[2, p,,, n_ages - 1], dim = c(n_regions, n_regions))

    for(o in 1:n_regions) {
      # catch-at-age for the penultimate age. tmp_caa holds this origin cohort's survivors season to
      # season, so each season's catch comes from the fish still alive to be caught in it
      tmp_caa = tmp_fished[o,]
      for(seas in 1:n_seas) {
        ret_F_seas = ret_F_by_region(p, n_ages - 1, seas)
        disc_F_seas = disc_F_by_region(p, n_ages - 1, seas)
        landed_F_seas = landed_F_by_region(p, n_ages - 1, seas, is_discard_fleet)
        Z_seas = natmort[p,,seas, n_ages - 1] * seasdur[seas] + ret_F_seas + disc_F_seas
        CAA[p, o,, seas, n_ages - 1] = catch_at_age(tmp_caa, Movement[p,,,seas, n_ages - 1], Z_seas,
                                                    Mrate[p,,,seas, n_ages - 1], seasdur[seas],
                                                    landed_F_seas, move_timing, expm_nsub = expm_nsub)
        tmp_caa = advance_seas(tmp_caa, Movement[p,,,seas, n_ages - 1], Z_seas,
                               Mrate[p,,,seas, n_ages - 1], seasdur[seas], move_timing, expm_nsub = expm_nsub)
      }
    }

    if(spawn_seas > 1) {
      for(o in 1:n_regions) {
        for(seas in 1:(spawn_seas - 1)) {

          # get mortality
          M_seas = natmort[p,,seas, n_ages - 1] * seasdur[seas]
          ret_F_seas = ret_F_by_region(p, n_ages - 1, seas)
          disc_F_seas = disc_F_by_region(p, n_ages - 1, seas)

          # Apply seasonal movement and mortality together, per move_timing
          tmp_unfished[o,] = advance_seas(tmp_unfished[o,], Movement[p,,,seas, n_ages - 1], M_seas,
                                          Mrate[p,,,seas, n_ages - 1], seasdur[seas], move_timing, expm_nsub = expm_nsub)
          tmp_fished[o,]   = advance_seas(tmp_fished[o,], Movement[p,,,seas, n_ages - 1],
                                          M_seas + ret_F_seas + disc_F_seas,
                                          Mrate[p,,,seas, n_ages - 1], seasdur[seas], move_timing, expm_nsub = expm_nsub)

        } # end seas loop
      } # end o loop
    }

    ## Advance penultimate age into spawning biomass / season; t_spawn folded in here
    tmp_unfished_spawn = tmp_unfished; tmp_fished_spawn = tmp_fished
    Zu_spawn = natmort[p,,spawn_seas, n_ages - 1] * seasdur[spawn_seas]
    Zf_spawn = Zu_spawn + ret_F_by_region(p, n_ages - 1, spawn_seas) + disc_F_by_region(p, n_ages - 1, spawn_seas)
    for(o in 1:n_regions) {
      tmp_unfished_spawn[o,] = spawn_state(tmp_unfished_spawn[o,], Movement[p,,, spawn_seas, n_ages - 1],
                                           Zu_spawn, Mrate[p,,, spawn_seas, n_ages - 1],
                                           seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
      tmp_fished_spawn[o,]   = spawn_state(tmp_fished_spawn[o,], Movement[p,,, spawn_seas, n_ages - 1],
                                           Zf_spawn, Mrate[p,,, spawn_seas, n_ages - 1],
                                           seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)

      # If single season natal homing population
      if(n_pop > 1 && n_seas == 1) {
        # Get NAA during spawning in single season case
        tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% sgl_seas_spawning_movement[p,,, n_ages - 1])
        tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% sgl_seas_spawning_movement[p,,, n_ages - 1])
      }
    }

    # get spawning biomass
    SB_age[1, p,,, n_ages - 1] = t(t(tmp_unfished_spawn) * WAA[p,, spawn_seas, n_ages - 1] *
                                     MatAA[p,, spawn_seas, n_ages - 1])
    SB_age[2, p,,, n_ages - 1] = t(t(tmp_fished_spawn) * WAA[p,, spawn_seas, n_ages - 1] *
                                     MatAA[p,, spawn_seas, n_ages - 1])

    ## Plus group analytical solution
    # Get fishing mortality for penultimate and plus ages
    F_penult = F_plus = matrix(0, n_regions, n_seas)
    for(seas in 1:n_seas) {
      F_penult[, seas] = ret_F_by_region(p, n_ages - 1, seas) + disc_F_by_region(p, n_ages - 1, seas)
      F_plus[, seas]   = ret_F_by_region(p, n_ages,     seas) + disc_F_by_region(p, n_ages,     seas)
    }

    # Build transition matrices
    Ts <- build_plus_group_T(
      M_penult = array(natmort[p,,, n_ages - 1], dim = c(n_regions, n_seas)),
      M_plus = array(natmort[p,,, n_ages], dim = c(n_regions, n_seas)),
      F_penult   = F_penult,
      F_plus     = F_plus,
      Mov_penult = array(Movement[p,,, , n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
      Mov_plus   = array(Movement[p,,, , n_ages],     dim = c(n_regions, n_regions, n_seas)),
      n_regions  = n_regions,
      n_seas     = n_seas,
      seasdur    = seasdur,
      Mrate_penult = array(Mrate[p,,, , n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
      Mrate_plus   = array(Mrate[p,,, , n_ages],     dim = c(n_regions, n_regions, n_seas)),
      move_timing  = move_timing,
      expm_nsub = expm_nsub)

    for(o in 1:n_regions) {
      # Get and input plus group
      pg <- solve_plus_group(Ts, Nspr[1, p, o,, n_ages - 1], Nspr[2, p, o,, n_ages - 1], n_regions)
      Nspr[1, p, o,, n_ages] <- pg$unfished
      Nspr[2, p, o,, n_ages] <- pg$fished
    }

    ## Plus group: catch, then advance to spawning season
    tmp_unfished = array(Nspr[1, p,,, n_ages], dim = c(n_regions, n_regions))
    tmp_fished   = array(Nspr[2, p,,, n_ages], dim = c(n_regions, n_regions))

    for(o in 1:n_regions) {
      # Catch-at-age for plus group, taken from this origin cohort's survivors in each season
      tmp_caa = tmp_fished[o,]
      for(seas in 1:n_seas) {
        ret_F_seas = ret_F_by_region(p, n_ages,     seas)
        disc_F_seas = disc_F_by_region(p, n_ages,     seas)
        landed_F_seas = landed_F_by_region(p, n_ages, seas, is_discard_fleet)
        Z_seas = natmort[p,,seas, n_ages] * seasdur[seas] + ret_F_seas + disc_F_seas
        CAA[p, o,, seas, n_ages] = catch_at_age(tmp_caa, Movement[p,,,seas, n_ages], Z_seas,
                                                Mrate[p,,,seas, n_ages], seasdur[seas],
                                                landed_F_seas, move_timing, expm_nsub = expm_nsub)
        tmp_caa = advance_seas(tmp_caa, Movement[p,,,seas, n_ages], Z_seas,
                               Mrate[p,,,seas, n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
      }
    }

    if(spawn_seas > 1) {
      for(o in 1:n_regions) {
        for(seas in 1:(spawn_seas - 1)) {

          # get mortality
          M_seas = natmort[p,,seas, n_ages] * seasdur[seas]
          F_seas = F_plus[, seas]

          # Apply seasonal movement and mortality together, per move_timing
          tmp_unfished[o,] = advance_seas(tmp_unfished[o,], Movement[p,,,seas, n_ages], M_seas,
                                          Mrate[p,,,seas, n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)
          tmp_fished[o,]   = advance_seas(tmp_fished[o,], Movement[p,,,seas, n_ages], M_seas + F_seas,
                                          Mrate[p,,,seas, n_ages], seasdur[seas], move_timing, expm_nsub = expm_nsub)

        } # end seas loop
      } # end o loop
    }

    ## Plus group spawning biomass; t_spawn mortality folded into spawn_state
    tmp_unfished_spawn = tmp_unfished; tmp_fished_spawn = tmp_fished
    Zu_spawn_plus = natmort[p,,spawn_seas, n_ages] * seasdur[spawn_seas]
    Zf_spawn_plus = Zu_spawn_plus + F_plus[, spawn_seas]
    for(o in 1:n_regions) {
      tmp_unfished_spawn[o,] = spawn_state(tmp_unfished_spawn[o,], Movement[p,,, spawn_seas, n_ages],
                                           Zu_spawn_plus, Mrate[p,,, spawn_seas, n_ages],
                                           seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)
      tmp_fished_spawn[o,]   = spawn_state(tmp_fished_spawn[o,], Movement[p,,, spawn_seas, n_ages],
                                           Zf_spawn_plus, Mrate[p,,, spawn_seas, n_ages],
                                           seasdur[spawn_seas], t_spawn, move_timing, expm_nsub = expm_nsub)

      # If single season natal homing population
      if(n_pop > 1 && n_seas == 1) {
        # Get NAA during spawning in single season case
        tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% sgl_seas_spawning_movement[p,,, n_ages])
        tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% sgl_seas_spawning_movement[p,,, n_ages])
      }
    }

    # get spawning biomass
    SB_age[1, p,,, n_ages] = t(t(tmp_unfished_spawn) * WAA[p,, spawn_seas, n_ages] *
                                 MatAA[p,, spawn_seas, n_ages])
    SB_age[2, p,,, n_ages] = t(t(tmp_fished_spawn) * WAA[p,, spawn_seas, n_ages] *
                                 MatAA[p,, spawn_seas, n_ages])
  } # end p loop

  # Determine equilibrium recruitment for destination region
  # parse out and compute unfished and fished spawning biomass per recruit
  for(p in 1:n_pop) {
    for(o in 1:n_regions) {
      for(d in 1:n_regions) {
        SB_unfished_mat[p, o, d] = sum(SB_age[1, p, o, d, ])  # unfished
        SB_fished_mat[p, o, d]   = sum(SB_age[2, p, o, d, ])  # fished at Fmsy
      } # end d loop
    } # end o loop
  } # end p loop

  # Sum over origins to get total SBPR at each destination by population
  SBPR_fished   <- matrix(0, n_pop, n_regions)  # [p, d]
  SBPR_unfished <- matrix(0, n_pop, n_regions)  # [p, d]
  for(p in 1:n_pop) for(d in 1:n_regions) {
    SBPR_fished[p, d]   <- sum(SB_fished_mat[p,, d])
    SBPR_unfished[p, d] <- sum(SB_unfished_mat[p,, d])
  }

  # Virgin effective SSB at each pop's natal region
  # divide stray contributions by n_pop_in_region to preserve mass balance
  eff_SSB0_virgin = rep(0, n_pop)
  for(p2 in 1:n_pop) {
    r = natal_region[p2]
    for(p in 1:n_pop) {
      sc = if(p == p2) 1 else stray_rate[p] / n_pop_in_region[r]
      eff_SSB0_virgin[p2] = eff_SSB0_virgin[p2] + sc * SBPR_unfished[p, r] * R0[p]
    }
  }

  h_p <- h[cbind(1:n_pop, natal_region)]   # steepness at each population's natal region

  # define initial guess to solve for equilibrium recruitment from origin region
  Req_o <- R0

  for(nit in 1:newton_steps) {

    # Effective fished SSB at each pop's natal region
    # divide stray contributions by n_pop_in_region to preserve mass balance
    eff_SSB_fished <- rep(0, n_pop)
    for(p2 in 1:n_pop) {
      r <- natal_region[p2]
      for(p in 1:n_pop) {
        sc <- if(p == p2) 1 else stray_rate[p] / n_pop_in_region[r]
        eff_SSB_fished[p2] <- eff_SSB_fished[p2] + sc * SBPR_fished[p, r] * Req_o[p]
      }
    }

    # BH equilibrium total recruitment per pop
    g_vec <- equil_rec_ssb(eff_SSB_fished, eff_SSB0_virgin, R0, h_p, rec_model)

    # find values of origin recruitment consistent w/ destination recruitment such that pop'n is in equilibrium
    iter_vec = Req_o - g_vec

    # Jacobian: dg/deff = (A*B) / (B + C*x)^2, then chain rule through eff_SSB_fished
    dg_deff <- equil_rec_ssb_deriv(eff_SSB_fished, eff_SSB0_virgin, R0, h_p, rec_model)
    J <- diag(n_pop)
    for(p2 in 1:n_pop) {
      r <- natal_region[p2]
      for(p in 1:n_pop) {
        sc <- if(p == p2) 1 else stray_rate[p] / n_pop_in_region[r]
        J[p2, p] <- J[p2, p] - dg_deff[p2] * sc * SBPR_fished[p, r]
      }
    }
    Req_o <- Req_o - solve(J, iter_vec) # newton raphson update

  } # end newton steps

  # get destination region yield
  for(d in 1:n_regions) {
    for(p in 1:n_pop) {
      tmp <- 0
      for(seas in 1:n_seas) {
        for(o in 1:n_regions) {
          tmp <- tmp + sum(CAA[p, o, d, seas, ] * WAA[p, d, seas, ]) * rec_region_prop[p, o] * Req_o[p]
        } # end o loop
      } # end seas loop
      Yield_pr[p,d] <- tmp
    } # end p loop
  } # end d loop

  # Get spawning biomass per recruit summed across origins (SBPR already sums origins in SBPR_fished/unfished)
  SB  <- matrix(0, n_pop, n_regions)
  SB0 <- matrix(0, n_pop, n_regions)
  for(p in 1:n_pop) for(r in 1:n_regions) {
    SB[p, r]  <- SBPR_fished[p, r]
    SB0[p, r] <- SBPR_unfished[p, r]
  }

  # maximize total yield
  Yield_total = sum(Yield_pr)
  obj_fun = -Yield_total

  RTMB::REPORT(Fmsy)
  RTMB::REPORT(Req_o)
  RTMB::REPORT(Yield_pr)
  RTMB::REPORT(Yield_total)
  RTMB::REPORT(iter_vec)
  RTMB::REPORT(SB_fished_mat)
  RTMB::REPORT(SB_unfished_mat)
  RTMB::REPORT(Nspr)
  RTMB::REPORT(SBPR_unfished)
  RTMB::REPORT(SBPR_fished)
  RTMB::REPORT(SB_age)
  RTMB::REPORT(SB)
  RTMB::REPORT(SB0)

  return(obj_fun)
}
