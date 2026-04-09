#' Optimize Reference Point Models
#'
#' Constructs and optimizes an RTMB automatic differentiation objective function
#' for estimating fisheries reference points (e.g., SPR-based or Fmsy-based
#' biological reference points). After optimization, retrieves the model report
#' and standard deviation report from the best parameter estimates.
#'
#' @param model_name Function. An RTMB-compatible model function (e.g.,
#'   \code{global_SPR}, \code{global_BH_Fmsy}, \code{local_BH_Fmsy}) that
#'   defines the objective function for the reference point calculation.
#' @param data_list List. A named list of data inputs passed to \code{model_name}
#'   via \code{\link[RTMB]{MakeADFun}}. Should contain all quantities treated as
#'   fixed data within the reference point model (e.g., biological parameters,
#'   selectivity, movement matrices).
#' @param pars_list List. A named list of initial parameter values passed to
#'   \code{\link[RTMB]{MakeADFun}}. These are the parameters over which the
#'   objective function is optimized (e.g., \code{ln_Fmsy}, \code{ln_F_spr}).
#'
#' @returns An RTMB AD function object (list) with the following additional
#'   elements appended after optimization:
#'   \describe{
#'     \item{\code{$optim}}{Output from \code{\link[stats]{nlminb}}, including
#'       convergence code, final objective value, and optimized parameter
#'       estimates.}
#'     \item{\code{$rep}}{Named list of reported quantities from the model
#'       (e.g., equilibrium SSB, yield, reference point values), evaluated at
#'       the best parameter estimates via \code{$env$last.par.best}.}
#'     \item{\code{$sd_rep}}{Output from \code{\link[RTMB]{sdreport}}, containing
#'       standard errors and summary statistics for all estimated and derived
#'       quantities.}
#'   }
#'
#'
#' @keywords internal
optim_ref_pts <- function(model_name, data_list, pars_list) {
  tmp_obj <- RTMB::MakeADFun(cmb(model_name, data_list), parameters = pars_list, random = NULL, silent = TRUE)
  tmp_obj$optim <- stats::nlminb(tmp_obj$par, tmp_obj$fn, tmp_obj$gr, control = list(iter.max = 1e6, eval.max = 1e6, rel.tol = 1e-15))
  tmp_obj$rep <- tmp_obj$report(tmp_obj$env$last.par.best) # get report
  tmp_obj$sd_rep <- sdreport(tmp_obj)
  return(tmp_obj)
}

#' Build annual transition matrices for the plus-group analytical solution
#'
#' Constructs the four annual transition matrices needed to solve for the
#' equilibrium plus-group abundance analytically. Each matrix accumulates
#' survival and movement across all seasons for either the penultimate age or
#' the plus-group age, under either unfished or fished conditions.
#'
#' The equilibrium plus-group vector \eqn{N_+} satisfies
#' \deqn{N_+ = T_+ N_+ + T_{n-1} N_{n-1}}
#' which rearranges to \eqn{(I - T_+) N_+ = T_{n-1} N_{n-1}}, solved in
#' \code{\link{solve_plus_group}}.
#'
#' All arguments are sliced by the caller to remove the population dimension,
#' so this helper works identically for the single-population spatial case
#' (\code{global_BH_Fmsy}, \code{local_BH_Fmsy_sglpop}) and the
#' multi-population case (\code{global_SPR}, \code{local_BH_Fmsy_multipop}).
#'
#' @param M_penult Numeric vector \code{[n_regions]}. Natural mortality for
#'   the penultimate age class, used as an annual rate (scaled by
#'   \code{seasdur} internally).
#' @param M_plus Numeric vector \code{[n_regions]}. Natural mortality for the
#'   plus-group age class.
#' @param F_penult Numeric matrix \code{[n_regions, n_seas]}. Total fishing
#'   mortality per season for the penultimate age, already summed across
#'   fleets.
#' @param F_plus Numeric matrix \code{[n_regions, n_seas]}. Total fishing
#'   mortality per season for the plus-group age.
#' @param Mov_penult Numeric array \code{[n_regions, n_regions, n_seas]}.
#'   Movement transition matrices for the penultimate age. Entry
#'   \code{[r1, r2, s]} is the probability of moving from region \code{r1}
#'   to region \code{r2} in season \code{s}.
#' @param Mov_plus Numeric array \code{[n_regions, n_regions, n_seas]}.
#'   Movement transition matrices for the plus-group age.
#' @param n_regions Integer. Number of spatial regions.
#' @param n_seas Integer. Number of seasons.
#' @param seasdur Numeric vector \code{[n_seas]}. Fractional duration of each
#'   season (must sum to one).
#'
#' @return A named list with four transition matrices, each of dimension
#'   \code{[n_regions, n_regions]}:
#'   \describe{
#'     \item{\code{T_penult_unfished}}{Annual transition for the penultimate
#'       age under unfished conditions.}
#'     \item{\code{T_plus_unfished}}{Annual transition for the plus-group age
#'       under unfished conditions.}
#'     \item{\code{T_penult_fished}}{Annual transition for the penultimate age
#'       under fished conditions.}
#'     \item{\code{T_plus_fished}}{Annual transition for the plus-group age
#'       under fished conditions.}
#'   }
#'
#' @keywords internal
build_plus_group_T <- function(M_penult, M_plus, F_penult, F_plus,
                               Mov_penult, Mov_plus,
                               n_regions, n_seas, seasdur) {

  # initialize transition matrices
  T_pu <- T_lu <- T_pf <- T_lf <- diag(n_regions)

  for (seas in seq_len(n_seas)) {
    Su_p <- diag(exp(-M_penult * seasdur[seas]),  n_regions) # unfished survival
    Su_l <- diag(exp(-M_plus * seasdur[seas]),  n_regions) # unfished survival
    Sf_p <- diag(exp(-(M_penult * seasdur[seas] + F_penult[,seas])), n_regions) # fished
    Sf_l <- diag(exp(-(M_plus   * seasdur[seas] + F_plus[,seas])), n_regions) # fished
    Mp   <- Mov_penult[,, seas] # movement
    Ml   <- Mov_plus[,,   seas] # movement
    # transition matrices
    T_pu <- Su_p %*% t(Mp) %*% T_pu
    T_lu <- Su_l %*% t(Ml) %*% T_lu
    T_pf <- Sf_p %*% t(Mp) %*% T_pf
    T_lf <- Sf_l %*% t(Ml) %*% T_lf
  } # end seas loop


  list(T_penult_unfished = T_pu,
       T_plus_unfished = T_lu,
       T_penult_fished = T_pf,
       T_plus_fished = T_lf)
}


#' Solve for equilibrium plus-group numbers given transition matrices
#'
#' Given the annual transition matrices produced by
#' \code{\link{build_plus_group_T}} and the penultimate-age abundance vector,
#' solves the linear system \eqn{(I - T_+) N_+ = T_{n-1} N_{n-1}} for the
#' equilibrium plus-group abundance under both unfished and fished conditions.
#'
#' @param Ts Named list returned by \code{\link{build_plus_group_T}}.
#' @param N_penult_u Numeric vector \code{[n_regions]}. Unfished
#'   penultimate-age abundance (per-recruit) at the start of the year.
#' @param N_penult_f Numeric vector \code{[n_regions]}. Fished
#'   penultimate-age abundance (per-recruit) at the start of the year.
#' @param n_regions Integer. Number of spatial regions.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{unfished}}{Numeric vector \code{[n_regions]}. Equilibrium
#'       plus-group abundance per recruit under unfished conditions.}
#'     \item{\code{fished}}{Numeric vector \code{[n_regions]}. Equilibrium
#'       plus-group abundance per recruit under fished conditions.}
#'   }
#'
#' @keywords internal
solve_plus_group <- function(Ts, N_penult_u, N_penult_f, n_regions) {
  I <- diag(n_regions)
  list(
    unfished = solve(I - Ts$T_plus_unfished, Ts$T_penult_unfished %*% N_penult_u),
    fished   = solve(I - Ts$T_plus_fished,   Ts$T_penult_fished   %*% N_penult_f)
  )
}

#' Compute SPR reference point for a single-region or non-spatial model
#'
#' Calculates the spawning potential ratio (SPR) as a function of a trial
#' fishing mortality \eqn{F_x}, then returns a squared penalty
#' \eqn{100 (SPR - SPR_x)^2} that is minimised by the outer optimizer to find
#' \eqn{F_{SPR_x}}. Supports multiple populations via stray rates but does not
#' include spatial movement.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_F_x}}{Log-scale trial fishing mortality.}
#'   }
#' @param data Named list of RTMB data. Must contain:
#'   \describe{
#'     \item{\code{n_pop}}{Integer. Number of populations.}
#'     \item{\code{n_ages}}{Integer. Number of age classes.}
#'     \item{\code{n_seas}}{Integer. Number of seasons.}
#'     \item{\code{seasdur}}{Numeric vector \code{[n_seas]}. Season durations.}
#'     \item{\code{spawn_seas}}{Integer. Index of the spawning season.}
#'     \item{\code{t_spawn}}{Numeric. Fraction of spawning season elapsed
#'       before spawning (used for mid-season mortality correction).}
#'     \item{\code{F_fract_flt}}{Numeric array \code{[n_seas, n_fish_fleets]}.
#'       Fleet F fractions (see \code{\link{make_fratio}}).}
#'     \item{\code{fish_sel}}{Numeric array \code{[n_ages, n_fish_fleets]}.
#'       Fishery selectivity at age for females.}
#'     \item{\code{natmort}}{Numeric array \code{[n_pop, n_ages]}. Female
#'       natural mortality at age.}
#'     \item{\code{WAA}}{Numeric array \code{[n_pop, n_seas, n_ages]}. Female
#'       weight at age.}
#'     \item{\code{MatAA}}{Numeric array \code{[n_pop, n_seas, n_ages]}.
#'       Maturity at age.}
#'     \item{\code{sex_ratio_f}}{Numeric vector \code{[n_pop]}. Female sex
#'       ratio at recruitment.}
#'     \item{\code{rec_seas_prop}}{Numeric array \code{[n_pop, n_seas]}.
#'       Proportion of annual recruitment entering in each season.}
#'     \item{\code{stray_rate}}{Numeric vector \code{[n_pop]}. Per-population
#'       stray rate used to compute effective SSB across populations.}
#'     \item{\code{natal_region}}{Integer vector \code{[n_pop]}. Natal region
#'       index for each population.}
#'     \item{\code{n_pop_in_region}}{Integer vector \code{[n_regions]}. Number
#'       of populations per natal region.}
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
        F_seas = sum(F_fract_flt[seas, ] * F_x * fish_sel[j - 1, ])
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
                                                sum(F_fract_flt[seas,] * F_x * fish_sel[n_ages-1,]) ))

      } # end seas loop
    } # end p loop
  }

  # Get spawning biomass of penultimate age
  SB_age[1, , n_ages - 1] = tmp_unfished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[,n_ages - 1] * seasdur[spawn_seas])
  SB_age[2, , n_ages - 1] = tmp_fished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[,n_ages - 1] * seasdur[spawn_seas] +
                      sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[n_ages-1,]) ))


  # Plus group (scalar, no movement)
  M_penult = natmort[,n_ages - 1]
  Z_penult = natmort[,n_ages - 1] + F_x * sum(colSums(F_fract_flt) * fish_sel[n_ages - 1,])
  Z_plus = natmort[,n_ages] + F_x * sum(colSums(F_fract_flt) * fish_sel[n_ages,])
  M_plus = natmort[,n_ages]

  # Geometric series plus group
  Nspr[1,,n_ages] = Nspr[1,,n_ages-1] * exp(-M_penult) / (1 - exp(-M_plus))
  Nspr[2,,n_ages] = Nspr[2,,n_ages-1] * exp(-Z_penult) / (1 - exp(-Z_plus))

  # Advance plus group to spawning season
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
  SB_age[1,,n_ages] = tmp_unfished * WAA[,spawn_seas, n_ages] * MatAA[,spawn_seas, n_ages] *
    exp(-t_spawn * natmort[,n_ages] * seasdur[spawn_seas])
  SB_age[2,,n_ages] = tmp_fished * WAA[,spawn_seas, n_ages] * MatAA[,spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[,n_ages] * seasdur[spawn_seas] +
                      sum(F_fract_flt[spawn_seas,] * F_x * fish_sel[n_ages,])))

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
#' that is tracked across all regions and seasons under movement. A single
#' scalar \eqn{F_x} is applied uniformly across regions (scaled by
#' region-specific fleet fractions and selectivity). Returns the squared
#' penalty \eqn{100(SPR - SPR_x)^2} for optimisation.
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
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_F_x}}{Log-scale trial fishing mortality.}
#'   }
#' @param data Named list of RTMB data. Must contain:
#'   \describe{
#'     \item{\code{n_pop}}{Integer. Number of populations.}
#'     \item{\code{n_regions}}{Integer. Number of spatial regions.}
#'     \item{\code{n_ages}}{Integer. Number of age classes.}
#'     \item{\code{n_seas}}{Integer. Number of seasons.}
#'     \item{\code{seasdur}}{Numeric vector \code{[n_seas]}. Season durations.}
#'     \item{\code{spawn_seas}}{Integer. Index of the spawning season.}
#'     \item{\code{t_spawn}}{Numeric. Mid-season spawning timing correction.}
#'     \item{\code{F_fract_flt}}{Numeric array \code{[n_regions, n_seas, n_fish_fleets]}.
#'       Fleet F fractions by region.}
#'     \item{\code{fish_sel}}{Numeric array \code{[n_regions, n_ages, n_fish_fleets]}.
#'       Female fishery selectivity.}
#'     \item{\code{natmort}}{Numeric array \code{[n_pop, n_regions, n_ages]}.
#'       Female natural mortality.}
#'     \item{\code{WAA}}{Numeric array \code{[n_pop, n_regions, n_seas, n_ages]}.
#'       Female weight at age.}
#'     \item{\code{MatAA}}{Numeric array \code{[n_pop, n_regions, n_seas, n_ages]}.
#'       Maturity at age.}
#'     \item{\code{Movement}}{Numeric array \code{[n_pop, n_regions, n_regions, n_seas, n_ages]}.
#'       Seasonal movement transition matrices.}
#'     \item{\code{sgl_seas_spawning_movement}}{Numeric array
#'       \code{[n_pop, n_regions, n_regions, n_ages]}. Spawning movement for
#'       single-season natal homing models.}
#'     \item{\code{do_recruits_move}}{Integer (0/1). Whether age-1 recruits
#'       are subject to movement.}
#'     \item{\code{rec_region_prop}}{Numeric array \code{[n_pop, n_regions]}.
#'       Proportion of recruitment entering each region.}
#'     \item{\code{sex_ratio_f}}{Numeric array \code{[n_pop, n_regions]}.
#'       Female sex ratio at recruitment.}
#'     \item{\code{rec_seas_prop}}{Numeric array \code{[n_pop, n_seas]}.
#'       Seasonal recruitment proportions.}
#'     \item{\code{stray_rate}}{Numeric vector \code{[n_pop]}. Per-population
#'       stray rate.}
#'     \item{\code{natal_region}}{Integer vector \code{[n_pop]}. Natal region
#'       index for each population.}
#'     \item{\code{n_pop_in_region}}{Integer vector \code{[n_regions]}. Number
#'       of populations per natal region.}
#'     \item{\code{SPR_x}}{Numeric. Target SPR fraction.}
#'   }
#'
#' @return Numeric scalar. Squared penalty \eqn{(SPR - SPR_x)^2}.
#'
#'
#' @keywords internal
#' @import RTMB
global_SPR <- function(pars,
                       data
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # helper function to get f by region for a given age and season
  F_by_region <- function(age, seas) apply(F_fract_flt[, seas,, drop = FALSE] * F_x * fish_sel[, age,, drop = FALSE], 1, sum)

  F_x = exp(log_F_x) # Exponentiate reference points
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy

  # Set up the initial recruits
  for(p in 1:n_pop) Nspr[2,p,,1] = Nspr[1,p,,1] = rec_region_prop[p,] * sex_ratio_f[p,] * rec_seas_prop[p,1]

  ## Loop through ages
  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        # get mortality
        F_seas = F_by_region(j - 1, seas)
        M_seas = natmort[p,, j - 1] * seasdur[seas]
        Z_seas = F_seas + M_seas

        # extract out quantities
        tmp_unfished = Nspr[1,p,, j - 1]; tmp_fished = Nspr[2,p,, j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          add = rec_seas_prop[p, seas] * sex_ratio_f[p, ] * rec_region_prop[p, ]
          tmp_unfished = tmp_unfished + add
          tmp_fished = tmp_fished + add
        }

        ## Movement
        if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
          tmp_unfished = as.vector(tmp_unfished %*% Movement[p,,,seas, j - 1])
          tmp_fished   = as.vector(tmp_fished   %*% Movement[p,,,seas, j - 1])
        }

        ## Spawning biomass
        if (seas == spawn_seas) {

          # Extract temporary variables out
          tmp_unfished_spawn = tmp_unfished; tmp_fished_spawn = tmp_fished

          # If single season natal homing population
          if(n_pop > 1 && n_seas == 1) {
            # Get NAA during spawning in single season case
            tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
            tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
          }

          # Get spawning biomass per recruit
          SB_age[1,p,, j - 1] = tmp_unfished_spawn * WAA[p,, spawn_seas, j - 1] * MatAA[p,, spawn_seas, j - 1] *
            exp(-t_spawn * M_seas)
          SB_age[2,p,, j - 1] = tmp_fished_spawn * WAA[p,, spawn_seas, j - 1] * MatAA[p,, spawn_seas, j - 1] *
            exp(-t_spawn * Z_seas )

        }

        ## Mortality and ageing
        if (seas < n_seas) { # Within season mortality
          Nspr[1,p,, j - 1] = tmp_unfished * exp(-M_seas)
          Nspr[2,p,, j - 1] = tmp_fished * exp(-Z_seas)
        } else {
          # Ageing
          Nspr[1,p,, j] = tmp_unfished * exp(-M_seas)
          Nspr[2,p,, j] = tmp_fished * exp(-Z_seas)
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
        F_seas = F_by_region(n_ages - 1, seas)

        # Apply seasonal movement
        tmp_unfished[p,] = tmp_unfished[p,] %*% Movement[p,,,seas,n_ages-1]
        tmp_fished[p,] = tmp_fished[p,] %*% Movement[p,,,seas,n_ages-1]

        # Apply seasonal mortality
        tmp_unfished[p,] = tmp_unfished[p,] * exp(-M_seas)
        tmp_fished[p,] = tmp_fished[p,] * exp(-(M_seas + F_seas))

      } # end seas loop
    }
  }

  ## Advance penultimate age into spawning biomass / season
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
  SB_age[2,,, n_ages - 1] = tmp_fished_spawn * WAA[,, spawn_seas, n_ages - 1] * MatAA[,, spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[,, n_ages - 1] * seasdur[spawn_seas] +
                      apply(F_fract_flt[,spawn_seas,,drop=F] * F_x * fish_sel[,n_ages-1,,drop=F], 1, sum)))

  ## Plus group analytical solution
  # Get fishing mortality
  F_penult = F_plus = matrix(0, n_regions, n_seas)
  for(seas in 1:n_seas) {
    F_penult[,seas] = F_by_region(n_ages - 1, seas)
    F_plus[,seas] = F_by_region(n_ages, seas)
  }

  for(p in 1:n_pop) {

    # Build transition matrices
    Ts <- build_plus_group_T(
      M_penult   = natmort[p,, n_ages - 1],
      M_plus     = natmort[p,, n_ages],
      F_penult   = F_penult,
      F_plus     = F_plus,
      Mov_penult = array(Movement[p,,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
      Mov_plus   = array(Movement[p,,,, n_ages], dim = c(n_regions, n_regions, n_seas)),
      n_regions  = n_regions,
      n_seas = n_seas,
      seasdur = seasdur
    )

    # Get and input plus group
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

        # Apply seasonal movement
        tmp_unfished[p,] = tmp_unfished[p,] %*% Movement[p,,,seas,n_ages]
        tmp_fished[p,] = tmp_fished[p,] %*% Movement[p,,,seas,n_ages]

        # Apply seasonal mortality
        tmp_unfished[p,] = tmp_unfished[p,] * exp(-(natmort[p,,n_ages] * seasdur[seas]))
        tmp_fished[p,] = tmp_fished[p,] * exp(-(natmort[p,,n_ages] * seasdur[seas] + F_plus[,seas] ))

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

  # Get spawning biomass per recruit
  SB_age[1,,,n_ages] = tmp_unfished_spawn * WAA[,,spawn_seas,n_ages] * MatAA[,,spawn_seas,n_ages] *
    exp(-t_spawn * natmort[,,n_ages] * seasdur[spawn_seas])
  SB_age[2,,,n_ages] = tmp_fished_spawn  * WAA[,,spawn_seas,n_ages] * MatAA[,,spawn_seas,n_ages] *
    exp(-t_spawn * (natmort[,,n_ages] * seasdur[spawn_seas] + F_plus[,spawn_seas]))

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


#' Compute Beverton-Holt Fmsy for a single-region or non-spatial model
#'
#' Finds \eqn{F_{MSY}} by maximising equilibrium yield under a
#' Beverton-Holt stock-recruit relationship. Yield is computed from
#' spawning biomass per recruit (\eqn{\phi_F}), the BH equilibrium
#' recruitment formula, and catch-at-age integrated across all seasons.
#' Supports multiple populations via stray rates but does not include
#' spatial movement.
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
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield (minimised to
#'   find \eqn{F_{MSY}}).
#'
#'
#' @keywords internal
#' @import RTMB
single_region_BH_Fmsy <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # Exponentiate Fmsy
  Fmsy = exp(log_Fmsy)

  # set up containers
  SB_age = Nspr = array(0, dim = c(2, n_pop, n_ages))
  CAA = array(0, c(n_pop, n_seas, n_ages))

  # initialize recruits
  for(p in 1:n_pop) Nspr[,p,1] = sex_ratio_f[p] * rec_seas_prop[p,1]

  for(p in 1:n_pop) {
    for (j in 2:(n_ages - 1)) {
      for (seas in 1:n_seas) {

        # get mortality
        F_seas = sum(F_fract_flt[seas, ] * Fmsy * fish_sel[j - 1, ])
        M_a_seas = natmort[p, j - 1] * seasdur[seas]
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
        CAA[p, seas, j - 1] = tmp_fished * (F_seas / Z_seas) * (1 - exp(-Z_seas))

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
  for (seas in 1:n_seas) {
    F_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages-1,])
    Z_seas = natmort[, n_ages-1] * seasdur[seas] + F_seas
    CAA[,seas, n_ages-1] = tmp_fished_caa * (F_seas / Z_seas) * (1 - exp(-Z_seas))
    tmp_fished_caa = tmp_fished_caa * exp(-Z_seas)  # advance
  }

  # Spawning season stuff
  if (spawn_seas > 1) {
    for(p in 1:n_pop) {
      for (seas in 1:(spawn_seas - 1)) {
        # Exponential mortality
        M_seas = natmort[p,n_ages - 1] * seasdur[seas]
        tmp_unfished[p] = tmp_unfished[p] * exp(-M_seas)
        tmp_fished[p] = tmp_fished[p] * exp(-(sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages-1,]) + M_seas ))
      }
    }
  }

  # Get spawning biomass after mortality decrement
  SB_age[1,,n_ages - 1] = tmp_unfished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[,n_ages - 1] * seasdur[spawn_seas])
  SB_age[2,, n_ages - 1] = tmp_fished * WAA[,spawn_seas, n_ages - 1] * MatAA[,spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[,n_ages - 1] * seasdur[spawn_seas] +
                      sum(F_fract_flt[spawn_seas,] * Fmsy * fish_sel[n_ages-1,])))

  # Plus group (scalar, no movement)
  Z_plus = natmort[,n_ages] + sum(colSums(F_fract_flt) * Fmsy * fish_sel[n_ages,])
  M_plus = natmort[,n_ages]
  M_penult = natmort[, n_ages - 1]
  Z_penult = natmort[, n_ages - 1] + sum(colSums(F_fract_flt) * Fmsy * fish_sel[n_ages - 1, ])
  Nspr[1, ,n_ages] = Nspr[1, ,n_ages - 1] * exp(-M_penult) / (1 - exp(-M_plus))
  Nspr[2, ,n_ages] = Nspr[2, ,n_ages - 1] * exp(-Z_penult) / (1 - exp(-Z_plus))

  # Plus group catch
  tmp_unfished = Nspr[1,, n_ages]
  tmp_fished = Nspr[2,, n_ages]
  tmp_fished_caa = tmp_fished

  # Catch-at-age for plus group
  for (seas in 1:n_seas) {
    F_seas = sum(F_fract_flt[seas,] * Fmsy * fish_sel[n_ages,])
    Z_seas = natmort[,n_ages] * seasdur[seas] + F_seas
    CAA[, seas, n_ages] = tmp_fished_caa * (F_seas / Z_seas) * (1 - exp(-Z_seas))
    tmp_fished_caa = tmp_fished_caa * exp(-Z_seas)  # advance
  }

  # Advance plus group to spawning season
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
  Req = R0 * ((4 * h * effective_SB) - (1 - h) * effective_SB0) / ((5 * h - 1) * effective_SB)

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

#' Compute global Beverton-Holt Fmsy for a spatially explicit single-population model
#'
#' Finds a single \eqn{F_{MSY}} that is applied uniformly across all regions
#' by maximising total equilibrium yield under a Beverton-Holt stock-recruit
#' relationship. Cohorts are tracked spatially across regions and seasons under
#' movement. Only valid for single-population models (\code{n_pop = 1});
#' multi-population models should use \code{\link{local_BH_Fmsy_multipop}}.
#'
#' The plus-group is solved analytically using
#' \code{\link{build_plus_group_T}} and \code{\link{solve_plus_group}}.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_Fmsy}}{Log-scale trial \eqn{F_{MSY}}.}
#'   }
#' @param data Named list of RTMB data. Must contain all spatial fields
#'   required by \code{\link{global_SPR}} (excluding \code{SPR_x} and
#'   \code{stray_rate}) plus:
#'   \describe{
#'     \item{\code{h}}{Numeric scalar. Beverton-Holt steepness.}
#'     \item{\code{R0}}{Numeric scalar. Unfished equilibrium recruitment.}
#'     \item{\code{rec_region_prop}}{Numeric vector \code{[n_regions]}.
#'       Proportion of annual recruitment entering each region.}
#'     \item{\code{rec_seas_prop}}{Numeric vector \code{[n_seas]}. Seasonal
#'       recruitment proportions (no population dimension for single-pop).}
#'     \item{\code{sex_ratio_f}}{Numeric scalar or vector. Female sex ratio.}
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield (minimised to
#'   find \eqn{F_{MSY}}).
#'
#' @keywords internal
#' @import RTMB
global_BH_Fmsy <- function(pars,
                           data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # Helper to get total F by region for a given (age, season)
  F_by_region = function(age, seas) rowSums(F_fract_flt[, seas,, drop = FALSE] * Fmsy * fish_sel[, age,, drop = FALSE])

  # exponentitate reference points to "estimate"
  Fmsy = exp(log_Fmsy)
  SB_age = Nspr = array(0, dim = c(2, n_regions, n_ages)) # 2 slots in rows, for unfished, and fished at Fmsy
  CAA = array(0, c(n_regions, n_seas, n_ages)) # catch at age

  # Set up the initial recruits
  Nspr[1,,1] = rec_region_prop * sex_ratio_f * rec_seas_prop[1]
  Nspr[2,,1] = rec_region_prop * sex_ratio_f * rec_seas_prop[1]

  ## Loop through ages
  for (j in 2:(n_ages - 1)) {
    for (seas in 1:n_seas) {

      # get mortality
      F_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,j-1,,drop = F])
      M_a_seas = natmort[,j - 1] * seasdur[seas]
      Z_seas = M_a_seas + F_seas

      # extract quantities
      tmp_unfished = Nspr[1,, j - 1]
      tmp_fished   = Nspr[2,, j - 1]

      # add in seasonal recruits
      if(seas > 1 && j - 1 == 1) {
        add_rec = rec_seas_prop[seas] * sex_ratio_f
        tmp_unfished = tmp_unfished + add_rec
        tmp_fished   = tmp_fished   + add_rec
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

      # Catch-at-age (Baranov)
      CAA[,seas, j - 1] = tmp_fished * (F_seas / Z_seas) * (1 - exp(-Z_seas))

      ## Mortality and ageing
      if (seas < n_seas) { # Within season mortality
        Nspr[1,, j - 1] = tmp_unfished * exp(-M_a_seas)
        Nspr[2,, j - 1] = tmp_fished * exp(-Z_seas)
      } else {
        # Ageing
        Nspr[1,, j] = tmp_unfished * exp(-M_a_seas)
        Nspr[2,, j] = tmp_fished * exp(-Z_seas)
      }
    }
  }

  # Age n_ages-1 is now at start of year after the loop
  tmp_unfished = Nspr[1,,n_ages-1]
  tmp_fished = Nspr[2,,n_ages-1]

  # Catch-at-age for penultimate age
  for (seas in 1:n_seas) {
    # Compute F and Z for this age/season
    F_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages-1,,drop = F])
    Z_seas = natmort[,n_ages - 1] * seasdur[seas] + F_seas
    CAA[,seas, n_ages - 1] = tmp_fished * (F_seas / Z_seas) * (1 - exp(-Z_seas))
  }

  # advance into spawning season
  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      F_seas = F_by_region(n_ages - 1, seas)
      M_seas = natmort[,n_ages-1] * seasdur[seas]
      Z_seas = M_seas + F_seas

      # Apply seasonal movement and mortality
      tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages-1] * exp(-M_seas)
      tmp_fished = tmp_fished %*% Movement[,,seas,n_ages-1] * exp(-Z_seas)

    } # end seas loop
  }

  ## Penultimate age spawning biomass
  tmp_unfished = as.vector(tmp_unfished %*% Movement[,, spawn_seas, n_ages - 1])
  tmp_fished   = as.vector(tmp_fished %*% Movement[,, spawn_seas, n_ages - 1])
  SB_age[1,, n_ages - 1] = tmp_unfished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[, n_ages - 1] * seasdur[spawn_seas])
  SB_age[2,, n_ages - 1] = tmp_fished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[, n_ages - 1] * seasdur[spawn_seas] +
                      rowSums(F_fract_flt[,spawn_seas,,drop = F] * Fmsy * fish_sel[,n_ages - 1,,drop = F]) ))

  # Get fishing mortality for plus group calculations
  F_penult = F_plus = matrix(0, n_regions, n_seas)
  for(seas in 1:n_seas) {
    F_penult[, seas] = F_by_region(n_ages - 1, seas)
    F_plus[, seas] = F_by_region(n_ages,     seas)
  }

  ## Plus group analytical solution
  Ts = build_plus_group_T(
    M_penult   = natmort[, n_ages - 1],
    M_plus     = natmort[, n_ages],
    F_penult   = F_penult,
    F_plus     = F_plus,
    Mov_penult = array(Movement[,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
    Mov_plus   = array(Movement[,,, n_ages], dim = c(n_regions, n_regions, n_seas)),
    n_regions  = n_regions,
    n_seas = n_seas,
    seasdur = seasdur
  )

  pg = solve_plus_group(Ts, Nspr[1,, n_ages - 1], Nspr[2,, n_ages - 1], n_regions)
  Nspr[1,, n_ages] = pg$unfished
  Nspr[2,, n_ages] = pg$fished

  ## Plus group catch, then advance to spawning season
  tmp_unfished = Nspr[1,, n_ages]
  tmp_fished = Nspr[2,, n_ages]

  # Catch-at-age for plus group
  for (seas in 1:n_seas) {
    F_seas = F_plus[seas, ]
    Z_seas = natmort[,n_ages] * seasdur[seas] + F_seas
    CAA[,seas, n_ages] = tmp_fished * (F_seas / Z_seas) * (1 - exp(-Z_seas))
  }


  if(spawn_seas > 1) {
    for (seas in 1:(spawn_seas - 1)) {

      # get mortality
      M_seas = natmort[,n_ages] * seasdur[seas]
      F_seas = rowSums(F_fract_flt[,seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F])
      Z_seas = F_seas + M_seas

      # Apply seasonal movement and mortality
      tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages] * exp(-M_seas)
      tmp_fished = tmp_fished %*% Movement[,,seas,n_ages] * exp(-Z_seas)

    } # end seas loop
  }

  ## Plus group spawning biomass
  tmp_unfished = as.vector(tmp_unfished %*% Movement[,, spawn_seas, n_ages])
  tmp_fished   = as.vector(tmp_fished %*% Movement[,, spawn_seas, n_ages])
  SB_age[1,, n_ages] = tmp_unfished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] *
    exp(-t_spawn * natmort[, n_ages] * seasdur[spawn_seas])
  SB_age[2,, n_ages] = tmp_fished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[, n_ages] * seasdur[spawn_seas] +
                      rowSums(F_fract_flt[,spawn_seas,,drop = F] * Fmsy * fish_sel[,n_ages,,drop = F]) ))

  # Get spawning biomass per recruit to get spawning potential ratio
  SBPR_0 = sum(SB_age[1,,])
  SBPR_F = sum(SB_age[2,,])
  SPR = SBPR_F / SBPR_0

  # Get equilibrium recruitment
  Req = R0 * ( (4*h*SBPR_F) - (1 - h) * SBPR_0) / ((5 * h -1) * SBPR_F)

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

#' Compute local Beverton-Holt Fmsy for a spatially explicit single-population model
#'
#' Finds region-specific \eqn{F_{MSY}} values that jointly maximise total
#' equilibrium yield across all regions under a spatially explicit
#' Beverton-Holt stock-recruit relationship. Unlike \code{\link{global_BH_Fmsy}},
#' which constrains all regions to a single F, this function allows each region
#' to have its own optimal fishing mortality.
#'
#' Cohorts originating from each region are tracked separately through
#' movement, mortality, and ageing using a \code{[origin x destination]}
#' per-recruit accounting framework. The equilibrium recruitment vector
#' \eqn{R_{eq,o}} (recruits by origin region) is solved iteratively via a
#' Newton-Raphson algorithm that accounts for the full spatial redistribution
#' of spawning biomass. The plus-group is solved analytically using
#' \code{\link{build_plus_group_T}} and \code{\link{solve_plus_group}}.
#'
#' @param pars Named list of RTMB parameters. Must contain:
#'   \describe{
#'     \item{\code{log_Fmsy}}{Log-scale trial \eqn{F_{MSY}} values, one per
#'       region (length \code{n_regions}).}
#'   }
#' @param data Named list of RTMB data. Must contain all spatial fields
#'   required by \code{\link{global_SPR}} (excluding \code{SPR_x},
#'   \code{stray_rate}, and \code{natal_region}) plus:
#'   \describe{
#'     \item{\code{h}}{Numeric vector \code{[n_regions]}. Beverton-Holt
#'       steepness by region.}
#'     \item{\code{R0}}{Numeric scalar. Total unfished equilibrium
#'       recruitment.}
#'     \item{\code{rec_region_prop}}{Numeric vector \code{[n_regions]}.
#'       Proportion of annual recruitment entering each region.}
#'     \item{\code{newton_steps}}{Integer. Number of Newton-Raphson iterations
#'       used to solve for equilibrium recruitment by origin region.}
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield across all regions
#'   (minimised to find the vector of regional \eqn{F_{MSY}} values).
#'
#' @details
#' The Newton-Raphson solver finds the origin-region recruitment vector
#' \eqn{R_{eq,o}} satisfying the fixed-point condition that recruitment
#' produced at each destination region (via the BH relationship applied to
#' effective SSB) equals the recruitment attributed to that origin region.
#' The Jacobian is derived analytically from the quotient rule applied to the
#' BH formula and the chain rule through the spatial redistribution of SSB.
#'
#'
#' @keywords internal
#' @import RTMB
local_BH_Fmsy_sglpop <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # exponentitate reference points
  Fmsy = exp(log_Fmsy)

  # helper function to get f by region for a given age and season
  F_by_region <- function(age, seas) apply(F_fract_flt[, seas,, drop = FALSE] * Fmsy * fish_sel[, age,, drop = FALSE], 1, sum)

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
        F_seas = F_by_region(j - 1, seas)
        M_seas = natmort[, j - 1] * seasdur[seas]
        Z_seas = F_seas + M_seas

        # extract out quantities
        tmp_unfished = Nspr[1, o,, j - 1]; tmp_fished = Nspr[2, o,, j - 1]

        # add in seasonal recruits
        if(seas > 1 && j - 1 == 1) {
          tmp_unfished[o] = tmp_unfished[o] + rec_seas_prop[seas] * sex_ratio_f[o]
          tmp_fished[o]   = tmp_fished[o]   + rec_seas_prop[seas] * sex_ratio_f[o]
        }

        ## Movement
        if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
          tmp_unfished = as.vector(tmp_unfished %*% Movement[,,seas, j - 1])
          tmp_fished   = as.vector(tmp_fished   %*% Movement[,,seas, j - 1])
        }

        ## Spawning biomass
        if(seas == spawn_seas) {
          SB_age[1, o,, j - 1] = tmp_unfished * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] *
            exp(-t_spawn * M_seas)
          SB_age[2, o,, j - 1] = tmp_fished   * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] *
            exp(-t_spawn * Z_seas)
        }

        # Catch-at-age (Baranov)
        CAA[o,, seas, j - 1] = tmp_fished * (F_seas / Z_seas) * (1 - exp(-Z_seas))

        ## Mortality and ageing
        if(seas < n_seas) { # Within season mortality
          Nspr[1, o,, j - 1] = tmp_unfished * exp(-M_seas)
          Nspr[2, o,, j - 1] = tmp_fished   * exp(-Z_seas)
        } else {
          # Ageing
          Nspr[1, o,, j] = tmp_unfished * exp(-M_seas)
          Nspr[2, o,, j] = tmp_fished   * exp(-Z_seas)
        }

      } # end o loop
    } # end seas loop
  } # end j loop

  # Advance penultimate age into spawning season
  tmp_unfished = array(Nspr[1,,, n_ages - 1], dim = c(n_regions, n_regions))
  tmp_fished   = array(Nspr[2,,, n_ages - 1], dim = c(n_regions, n_regions))

  for(o in 1:n_regions) {
    # Catch-at-age for penultimate age
    for(seas in 1:n_seas) {
      F_seas = F_by_region(n_ages - 1, seas)
      Z_seas = natmort[, n_ages - 1] * seasdur[seas] + F_seas
      CAA[o,, seas, n_ages - 1] = tmp_fished[o,] * (F_seas / Z_seas) * (1 - exp(-Z_seas))
    }
  }

  if(spawn_seas > 1) {
    for(o in 1:n_regions) {
      for(seas in 1:(spawn_seas - 1)) {

        # get mortality
        M_seas = natmort[, n_ages - 1] * seasdur[seas]
        F_seas = F_by_region(n_ages - 1, seas)

        # Apply seasonal movement
        tmp_unfished[o,] = tmp_unfished[o,] %*% Movement[,,seas, n_ages - 1]
        tmp_fished[o,]   = tmp_fished[o,]   %*% Movement[,,seas, n_ages - 1]

        # Apply seasonal mortality
        tmp_unfished[o,] = tmp_unfished[o,] * exp(-M_seas)
        tmp_fished[o,]   = tmp_fished[o,]   * exp(-(M_seas + F_seas))

      } # end seas loop
    } # end o loop
  }

  ## Advance penultimate age into spawning biomass / season
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn   = tmp_fished
  for(o in 1:n_regions) {
    tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% Movement[,,spawn_seas, n_ages - 1])
    tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% Movement[,,spawn_seas, n_ages - 1])
  }

  SB_age[1,,, n_ages - 1] = tmp_unfished_spawn * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] *
    exp(-t_spawn * natmort[, n_ages - 1] * seasdur[spawn_seas])
  SB_age[2,,, n_ages - 1] = tmp_fished_spawn   * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] *
    exp(-t_spawn * (natmort[, n_ages - 1] * seasdur[spawn_seas] +
                      apply(F_fract_flt[, spawn_seas,, drop = F] * Fmsy * fish_sel[, n_ages - 1,, drop = F], 1, sum)))

  ## Plus group analytical solution
  # Get fishing mortality for penultimate and plus ages
  F_penult = F_plus = matrix(0, n_regions, n_seas)
  for(seas in 1:n_seas) {
    F_penult[, seas] = F_by_region(n_ages - 1, seas)
    F_plus[, seas]   = F_by_region(n_ages,     seas)
  }

  # Build transition matrices
  Ts <- build_plus_group_T(
    M_penult   = natmort[, n_ages - 1],
    M_plus     = natmort[, n_ages],
    F_penult   = F_penult,
    F_plus     = F_plus,
    Mov_penult = array(Movement[,,, n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
    Mov_plus   = array(Movement[,,, n_ages],     dim = c(n_regions, n_regions, n_seas)),
    n_regions  = n_regions,
    n_seas     = n_seas,
    seasdur    = seasdur
  )

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
    # Catch-at-age for plus group
    for(seas in 1:n_seas) {
      F_seas = F_plus[, seas]
      Z_seas = natmort[, n_ages] * seasdur[seas] + F_seas
      CAA[o,, seas, n_ages] = tmp_fished[o,] * (F_seas / Z_seas) * (1 - exp(-Z_seas))
    }
  }

  if(spawn_seas > 1) {
    for(o in 1:n_regions) {
      for(seas in 1:(spawn_seas - 1)) {

        # get mortality
        M_seas = natmort[, n_ages] * seasdur[seas]
        F_seas = F_plus[, seas]

        # Apply seasonal movement
        tmp_unfished[o,] = tmp_unfished[o,] %*% Movement[,,seas, n_ages]
        tmp_fished[o,]   = tmp_fished[o,]   %*% Movement[,,seas, n_ages]

        # Apply seasonal mortality
        tmp_unfished[o,] = tmp_unfished[o,] * exp(-M_seas)
        tmp_fished[o,]   = tmp_fished[o,]   * exp(-(M_seas + F_seas))

      } # end seas loop
    } # end o loop
  }

  ## Plus group spawning biomass
  tmp_unfished_spawn = tmp_unfished
  tmp_fished_spawn   = tmp_fished
  for(o in 1:n_regions) {
    tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% Movement[,,spawn_seas, n_ages])
    tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% Movement[,,spawn_seas, n_ages])
  }

  SB_age[1,,, n_ages] = tmp_unfished_spawn * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] *
    exp(-t_spawn * natmort[, n_ages] * seasdur[spawn_seas])
  SB_age[2,,, n_ages] = tmp_fished_spawn   * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] *
    exp(-t_spawn * (natmort[, n_ages] * seasdur[spawn_seas] + F_plus[, spawn_seas]))

  # Determine equilibrium recruitment for destination region
  # parse out and compute unfished and fished spawning biomass per recruit
  for(o in 1:n_regions) {
    for(d in 1:n_regions) {
      SB_unfished_mat[o, d] = sum(SB_age[1, o, d, ])  # unfished
      SB_fished_mat[o, d]   = sum(SB_age[2, o, d, ])  # fished at Fmsy
    } # end d loop
  } # end o loop

  A = 4 * h * rec_region_prop * R0 # define first part of the numerator of BH recruitment
  B = rep(0, n_regions) # define first part of the denominator of BH recruitment
  for(d in 1:n_regions) B[d] = (1 - h[d]) * sum(SB_unfished_mat[, d] * R0 * rec_region_prop)
  C = 5 * h - 1 # define second part of the denominator for BH recruitment

  # define initial guess to solve for equilibrium recruitment from origin region
  Req_o = R0 * rec_region_prop

  for(nit in 1:newton_steps) {
    # compute equilibrium spawning biomass (SSBR * Req) in destination region
    x_vec = as.numeric(t(SB_fished_mat) %*% Req_o)  # function of equilibrium recruitment in origin region (what we are solving for)
    numer_vec = A * x_vec # compute numerator of BH
    denom_vec = B + (C * x_vec) # compute denominator of BH
    g_vec = numer_vec / denom_vec # equilibrium recruitment in destination region

    # define root and define Jacobian
    iter_vec = Req_o - g_vec # find values of origin recruitment that are consistent w/ destination recruitment such that pop'n is in equilibrium

    # construct Jacobian for root
    # we need J = df (iter_vec)/dReq = dReq/dReq (or I) - dg/dReq
    # we basically want to know dg / dReq (how does destination equil rec change as origin equil rec change)
    # dg / dReq = (dg / dxk) * (dxk / dReq)
    # to get (dg / dxk), use quotient rule of (BH recruitment)
    # note that dxk / dReq S_2mat * Req = S_2mat
    dg_dxk = (A * B) / (denom_vec^2)
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

#' Compute local Beverton-Holt Fmsy for a spatially explicit multi-population model
#'
#' Multi-population extension of \code{\link{local_BH_Fmsy_sglpop}}. Finds
#' region-specific \eqn{F_{MSY}} values that jointly maximise total
#' equilibrium yield when multiple populations with distinct natal regions,
#' movement schedules, and Beverton-Holt parameters co-occupy the same spatial
#' domain.
#'
#' Each population's cohorts are tracked separately using a
#' \code{[population x origin x destination]} per-recruit accounting
#' framework. Effective SSB at each population's natal region accumulates
#' contributions from all source populations scaled by stray rates.
#' Equilibrium recruitment by population is solved via a Newton-Raphson
#' algorithm whose Jacobian accounts for cross-population SSB coupling through
#' straying. The plus-group is solved analytically using
#' \code{\link{build_plus_group_T}} and \code{\link{solve_plus_group}}.
#'
#' When \code{n_seas = 1}, \code{sgl_seas_spawning_movement} redistributes
#' fish to natal spawning grounds before SSB is computed.
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
#'       steepness; indexed at the natal region via \code{natal_region}.}
#'     \item{\code{R0}}{Numeric vector \code{[n_pop]}. Unfished equilibrium
#'       recruitment per population.}
#'     \item{\code{stray_rate}}{Numeric vector \code{[n_pop]}. Per-population
#'       stray rate controlling cross-population SSB contributions.}
#'     \item{\code{natal_region}}{Integer vector \code{[n_pop]}. Natal region
#'       index for each population.}
#'     \item{\code{n_pop_in_region}}{Integer vector \code{[n_regions]}. Number
#'       of populations per natal region.}
#'     \item{\code{newton_steps}}{Integer. Number of Newton-Raphson iterations
#'       used to solve for equilibrium recruitment by population.}
#'   }
#'
#' @return Numeric scalar. Negative total equilibrium yield across all regions
#'   (minimised to find the vector of regional \eqn{F_{MSY}} values).
#'
#' @keywords internal
#' @import RTMB
local_BH_Fmsy_multipop <- function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data) # get parameters and data

  # exponentitate reference points
  Fmsy = exp(log_Fmsy)

  # helper function to get f by region for a given age and season
  F_by_region <- function(age, seas) apply(F_fract_flt[, seas,, drop = FALSE] * Fmsy * fish_sel[, age,, drop = FALSE], 1, sum)

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

  # Advance penultimate age into spawning season
  for(p in 1:n_pop) {

    tmp_unfished = array(Nspr[1, p,,, n_ages - 1], dim = c(n_regions, n_regions))
    tmp_fished   = array(Nspr[2, p,,, n_ages - 1], dim = c(n_regions, n_regions))

    for(o in 1:n_regions) {
      # Catch-at-age for penultimate age
      for(seas in 1:n_seas) {
        F_seas = F_by_region(n_ages - 1, seas)
        Z_seas = natmort[p,, n_ages - 1] * seasdur[seas] + F_seas
        CAA[p, o,, seas, n_ages - 1] = tmp_fished[o,] * (F_seas / Z_seas) * (1 - exp(-Z_seas))
      }
    }

    if(spawn_seas > 1) {
      for(o in 1:n_regions) {
        for(seas in 1:(spawn_seas - 1)) {

          # get mortality
          M_seas = natmort[p,, n_ages - 1] * seasdur[seas]
          F_seas = F_by_region(n_ages - 1, seas)

          # Apply seasonal movement
          tmp_unfished[o,] = tmp_unfished[o,] %*% Movement[p,,,seas, n_ages - 1]
          tmp_fished[o,]   = tmp_fished[o,]   %*% Movement[p,,,seas, n_ages - 1]

          # Apply seasonal mortality
          tmp_unfished[o,] = tmp_unfished[o,] * exp(-M_seas)
          tmp_fished[o,]   = tmp_fished[o,]   * exp(-(M_seas + F_seas))

        } # end seas loop
      } # end o loop
    }

    ## Advance penultimate age into spawning biomass / season
    tmp_unfished_spawn = tmp_unfished; tmp_fished_spawn = tmp_fished
    for(o in 1:n_regions) {
      tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% Movement[p,,, spawn_seas, n_ages - 1])
      tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% Movement[p,,, spawn_seas, n_ages - 1])

      # If single season natal homing population
      if(n_pop > 1 && n_seas == 1) {
        # Get NAA during spawning in single season case
        tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% sgl_seas_spawning_movement[p,,, n_ages - 1])
        tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% sgl_seas_spawning_movement[p,,, n_ages - 1])
      }
    }

    # get spawning biomass
    SB_age[1, p,,, n_ages - 1] = t(t(tmp_unfished_spawn) * WAA[p,, spawn_seas, n_ages - 1] *
                                     MatAA[p,, spawn_seas, n_ages - 1] *
                                     exp(-t_spawn * natmort[p,, n_ages - 1] * seasdur[spawn_seas]))
    SB_age[2, p,,, n_ages - 1] = t(t(tmp_fished_spawn) * WAA[p,, spawn_seas, n_ages - 1] *
                                     MatAA[p,, spawn_seas, n_ages - 1] *
                                     exp(-t_spawn * (natmort[p,, n_ages - 1] * seasdur[spawn_seas] +
                                                       apply(F_fract_flt[, spawn_seas,, drop = F] * Fmsy * fish_sel[, n_ages - 1,, drop = F], 1, sum))))

    ## Plus group analytical solution
    # Get fishing mortality for penultimate and plus ages
    F_penult = F_plus = matrix(0, n_regions, n_seas)
    for(seas in 1:n_seas) {
      F_penult[, seas] = F_by_region(n_ages - 1, seas)
      F_plus[, seas]   = F_by_region(n_ages,     seas)
    }

    # Build transition matrices
    Ts <- build_plus_group_T(
      M_penult   = natmort[p,, n_ages - 1],
      M_plus     = natmort[p,, n_ages],
      F_penult   = F_penult,
      F_plus     = F_plus,
      Mov_penult = array(Movement[p,,, , n_ages - 1], dim = c(n_regions, n_regions, n_seas)),
      Mov_plus   = array(Movement[p,,, , n_ages],     dim = c(n_regions, n_regions, n_seas)),
      n_regions  = n_regions,
      n_seas     = n_seas,
      seasdur    = seasdur
    )

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
      # Catch-at-age for plus group
      for(seas in 1:n_seas) {
        F_seas = F_plus[, seas]
        Z_seas = natmort[p,, n_ages] * seasdur[seas] + F_seas
        CAA[p, o,, seas, n_ages] = tmp_fished[o,] * (F_seas / Z_seas) * (1 - exp(-Z_seas))
      }
    }

    if(spawn_seas > 1) {
      for(o in 1:n_regions) {
        for(seas in 1:(spawn_seas - 1)) {

          # get mortality
          M_seas = natmort[p,, n_ages] * seasdur[seas]
          F_seas = F_plus[, seas]

          # Apply seasonal movement
          tmp_unfished[o,] = tmp_unfished[o,] %*% Movement[p,,,seas, n_ages]
          tmp_fished[o,]   = tmp_fished[o,]   %*% Movement[p,,,seas, n_ages]

          # Apply seasonal mortality
          tmp_unfished[o,] = tmp_unfished[o,] * exp(-M_seas)
          tmp_fished[o,]   = tmp_fished[o,]   * exp(-(M_seas + F_seas))

        } # end seas loop
      } # end o loop
    }

    ## Plus group spawning biomass
    tmp_unfished_spawn = tmp_unfished; tmp_fished_spawn = tmp_fished
    for(o in 1:n_regions) {
      tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% Movement[p,,, spawn_seas, n_ages])
      tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% Movement[p,,, spawn_seas, n_ages])

      # If single season natal homing population
      if(n_pop > 1 && n_seas == 1) {
        # Get NAA during spawning in single season case
        tmp_unfished_spawn[o,] = as.vector(tmp_unfished_spawn[o,] %*% sgl_seas_spawning_movement[p,,, n_ages])
        tmp_fished_spawn[o,]   = as.vector(tmp_fished_spawn[o,]   %*% sgl_seas_spawning_movement[p,,, n_ages])
      }
    }

    # get spawning biomass
    SB_age[1, p,,, n_ages] = t(t(tmp_unfished_spawn) * WAA[p,, spawn_seas, n_ages] *
                                 MatAA[p,, spawn_seas, n_ages] *
                                 exp(-t_spawn * natmort[p,, n_ages] * seasdur[spawn_seas]))
    SB_age[2, p,,, n_ages] = t(t(tmp_fished_spawn) * WAA[p,, spawn_seas, n_ages] *
                                 MatAA[p,, spawn_seas, n_ages] *
                                 exp(-t_spawn * (natmort[p,, n_ages] * seasdur[spawn_seas] + F_plus[, spawn_seas])))
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

  A_p <- 4 * h[cbind(1:n_pop, natal_region)] * R0                   # define first part of the numerator of BH recruitment
  B_p <- (1 - h[cbind(1:n_pop, natal_region)]) * eff_SSB0_virgin    # define first part of the denominator of BH recruitment
  C_p <- 5 * h[cbind(1:n_pop, natal_region)] - 1                    # define second part of the denominator for BH recruitment

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
    g_vec <- A_p * eff_SSB_fished / (B_p + C_p * eff_SSB_fished)

    # find values of origin recruitment consistent w/ destination recruitment such that pop'n is in equilibrium
    iter_vec = Req_o - g_vec

    # Jacobian: dg/deff = (A*B) / (B + C*x)^2, then chain rule through eff_SSB_fished
    dg_deff <- (A_p * B_p) / (B_p + C_p * eff_SSB_fished)^2
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


#' Compute fishing and biological reference points from an assessment or simulation
#'
#' Wrapper that constructs the appropriate data list, calls the relevant
#' inner objective function via RTMB, and returns fishing and biological
#' reference points for use in projections or harvest control rules. Supports
#' single-region and spatially explicit multi-region models, with options for
#' SPR-based or Beverton-Holt MSY-based reference points.
#'
#' @param data List. SPoRC data object containing age structure, weight-at-age,
#'   maturity, natural mortality, seasons, and spatial configuration.
#' @param rep List. SPoRC report object from RTMB containing estimated or
#'   simulated quantities including \code{Fmort}, \code{fish_sel},
#'   \code{natmort}, \code{Rec}, \code{SSB}, \code{h_trans}, \code{R0},
#'   \code{rec_region_prop}, \code{rec_seas_prop}, \code{Movement}, and
#'   \code{stray_rate}.
#' @param SPR_x Numeric. Target spawning potential ratio fraction (e.g. 0.4).
#'   Required when \code{what} is \code{"SPR"}, \code{"independent_SPR"}, or
#'   \code{"global_SPR"}.
#' @param t_spawn Numeric. Fraction of the spawning season elapsed before
#'   spawning, used for the mid-season mortality correction. Default = 0.
#' @param sex_ratio_f Numeric array \code{[n_pop, n_regions]}. Female sex
#'   ratio at recruitment. Default = 0.5 everywhere.
#' @param calc_rec_st_yr Integer. First year included when computing mean
#'   historical recruitment for biological reference point scaling. Default = 1.
#' @param rec_age Integer. Recruitment lag in years, used to exclude the most
#'   recent years from the mean recruitment calculation. Default = 1.
#' @param type Character. Spatial structure of the model:
#'   \describe{
#'     \item{\code{"single_region"}}{No spatial movement; supports
#'       \code{"SPR"} and \code{"BH_MSY"}.}
#'     \item{\code{"multi_region"}}{Spatially explicit; supports
#'       \code{"independent_SPR"}, \code{"independent_BH_MSY"},
#'       \code{"global_SPR"}, \code{"global_BH_MSY"}, and
#'       \code{"local_BH_MSY"}.}
#'   }
#' @param what Character. Reference point method:
#'   \describe{
#'     \item{\code{"SPR"}}{Single-region \eqn{F_{SPR_x}}.}
#'     \item{\code{"BH_MSY"}}{Single-region Beverton-Holt \eqn{F_{MSY}}.}
#'     \item{\code{"independent_SPR"}}{Per-region \eqn{F_{SPR_x}} computed
#'       independently for each region without movement.}
#'     \item{\code{"independent_BH_MSY"}}{Per-region \eqn{F_{MSY}} computed
#'       independently for each region without movement.}
#'     \item{\code{"global_SPR"}}{Single shared \eqn{F_{SPR_x}} with
#'       movement, integrated across all regions.}
#'     \item{\code{"global_BH_MSY"}}{Single shared \eqn{F_{MSY}} with
#'       movement. Valid for single-population models only.}
#'     \item{\code{"local_BH_MSY"}}{Region-specific \eqn{F_{MSY}} values
#'       that jointly maximise total yield with movement. Valid for both
#'       single- and multi-population models.}
#'   }
#' @param n_avg_yrs Integer. Number of terminal years over which demographic
#'   rates (selectivity, natural mortality, weight, maturity, movement) are
#'   averaged before computing reference points. Default = 1.
#' @param local_bh_msy_newton_steps Integer. Number of Newton-Raphson
#'   iterations used to solve for equilibrium recruitment by origin region
#'   when \code{what = "local_BH_MSY"}. Increase if convergence is suspect.
#'   Default = 6.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{f_ref_pt}}{Numeric vector \code{[n_regions]}. Fishing
#'       mortality reference point by region. All regions share the same value
#'       for global methods; regions have independent values for local or
#'       independent methods.}
#'     \item{\code{b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Equilibrium spawning biomass at the reference point by population and
#'       region (\eqn{SBPR_F \times R_{eq}} or \eqn{SBPR_F \times \bar{R}}).}
#'     \item{\code{virgin_b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Virgin (unfished) spawning biomass by population and region
#'       (\eqn{SBPR_0 \times R_0} or \eqn{SBPR_0 \times \bar{R}}).}
#'     \item{\code{pop_b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Population-specific effective spawning biomass at the reference point,
#'       evaluated at each population's natal region and incorporating stray
#'       contributions from other populations.}
#'     \item{\code{virgin_pop_b_ref_pt}}{Numeric array \code{[n_pop, n_regions]}.
#'       Population-specific effective virgin spawning biomass, evaluated at
#'       each population's natal region.}
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

  f_ref_pt <- vector()
  virgin_pop_b_ref_pt <- pop_b_ref_pt <- virgin_b_ref_pt <- b_ref_pt <- array(0, dim = c(n_pop, n_regions))

  # determine years to average over demographics
  n_yrs <- length(data$years)
  avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs

  # precompute number of populations per natal region (used throughout for mass-balance stray scaling)
  n_pop_in_region <- rep(0, n_regions)
  for(p in 1:n_pop) n_pop_in_region[data$natal_region[p]] <- n_pop_in_region[data$natal_region[p]] + 1

  if(type == "single_region") {

    if(n_regions > 1) stop("Single region reference points specified, but n_regions > 1!")

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

    # setup shared data lists
    data_list$F_fract_flt <- array(rep$Fmort[1,n_years,,] / sum(rep$Fmort[1,n_years,,]), dim = c(data$n_seas, data$n_fish_fleets)) # fishing mortality fraction
    fish_sel_avg <- apply(rep$fish_sel[1,avg_yrs,,1,,drop = FALSE], c(3,5), mean)
    data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
    natmort_avg <- apply(rep$natmort[,1,avg_yrs,,1,drop = FALSE], c(1,4), mean)
    data_list$natmort <- natmort_avg # get female natural mortality
    WAA_avg <- apply(data$WAA[,1,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
    data_list$WAA <- array(WAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # weight at age for females
    MatAA_avg <- apply(data$MatAA[,1,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
    data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # maturity at age for females

    # Other recruitment stuff
    data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio
    data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate
    data_list$rec_seas_prop <- array(rep$rec_seas_prop[,], dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
    data_list$natal_region <- data$natal_region
    data_list$n_pop_in_region <- n_pop_in_region

    if(what == 'SPR') {

      data_list$SPR_x <- SPR_x # SPR fraction
      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # Make adfun object
      tmp_obj <- optim_ref_pts(single_region_SPR, data_list, par_list)

      # Output reference points
      f_ref_pt[1] <- tmp_obj$rep$F_x

      # Compute population specific reference points, by using stray rates
      mean_rec <- apply(rep$Rec[,1,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], 1, mean)

      for(p2 in 1:n_pop) {
        r <- data$natal_region[p2]
        pop_b_ref_pt[p2,1]        <- tmp_obj$rep$SB[p2]  * mean_rec[p2]
        virgin_pop_b_ref_pt[p2,1] <- tmp_obj$rep$SB0[p2] * mean_rec[p2]
        for(p in 1:n_pop) {
          if(p != p2) {
            sc <- data_list$stray_rate[p] / n_pop_in_region[r]
            pop_b_ref_pt[p2,1]        <- pop_b_ref_pt[p2,1]        + sc * tmp_obj$rep$SB[p]  * mean_rec[p]
            virgin_pop_b_ref_pt[p2,1] <- virgin_pop_b_ref_pt[p2,1] + sc * tmp_obj$rep$SB0[p] * mean_rec[p]
          }
        }
      }

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

      # make adfun etc
      tmp_obj <- optim_ref_pts(single_region_BH_Fmsy, data_list, par_list)

      # Output reference points
      f_ref_pt[1] <- tmp_obj$rep$Fmsy

      # Accumulate biomass reference points
      for(p2 in 1:n_pop) {
        r <- data$natal_region[p2]
        pop_b_ref_pt[p2,1]        <- tmp_obj$rep$SB[p2]  * tmp_obj$rep$Req[p2]
        virgin_pop_b_ref_pt[p2,1] <- tmp_obj$rep$SB0[p2] * rep$R0[p2]
        for(p in 1:n_pop) {
          if(p != p2) {
            sc <- data_list$stray_rate[p] / n_pop_in_region[r]
            pop_b_ref_pt[p2,1]        <- pop_b_ref_pt[p2,1]        + sc * tmp_obj$rep$SB[p]  * tmp_obj$rep$Req[p]
            virgin_pop_b_ref_pt[p2,1] <- virgin_pop_b_ref_pt[p2,1] + sc * tmp_obj$rep$SB0[p] * data_list$R0[p]
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

    if(what %in% c("independent_SPR", 'independent_BH_MSY')) {

      tmp_obj <- list() # save optimized object as a list

      for(r in 1:data$n_regions) {

        # create shared data lists
        data_list$F_fract_flt <- array(rep$Fmort[r,n_years,,] / sum(rep$Fmort[r,n_years,,]), dim = c(data$n_seas, data$n_fish_fleets)) # get fleet F fraction to derive population level selectivity
        fish_sel_avg <- apply(rep$fish_sel[r,avg_yrs,,1,,drop = FALSE], c(1, 3, 4, 5), mean)
        data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
        natmort_avg <- apply(rep$natmort[,r,avg_yrs,,1,drop = FALSE], c(1,4), mean)
        data_list$natmort <- array(natmort_avg, dim = c(n_pop, n_ages)) # get female natural mortality
        WAA_avg <- apply(data$WAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$WAA <- array(WAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # weight at age for females
        MatAA_avg <- apply(data$MatAA[,r,avg_yrs,,,1,drop = FALSE], c(1,4,5), mean)
        data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, data$n_seas, n_ages)) # maturity at age for females
        data_list$sex_ratio_f <- array(sex_ratio_f[,r], dim = n_pop) # recruitment sex ratio
        data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate
        data_list$rec_seas_prop <- array(rep$rec_seas_prop, dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
        data_list$natal_region <- data$natal_region
        data_list$n_pop_in_region <- n_pop_in_region

        if(what == 'independent_SPR') {

          data_list$SPR_x <- SPR_x # SPR fraction

          par_list <- list() # set up parameter list
          par_list$log_F_x <- log(0.1) # F_x starting value

          # Make adfun object
          tmp_obj[[r]] <- optim_ref_pts(single_region_SPR, data_list, par_list)

          # Output reference points
          f_ref_pt[r] <- tmp_obj[[r]]$rep$F_x

          # Compute population specific reference points, by using stray rates
          mean_rec <- apply(rep$Rec[,r,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], 1, mean)
          if(n_pop > 1) {
            for(p2 in 1:n_pop) {
              rn <- data$natal_region[p2]
              pop_b_ref_pt[p2,r]        <- tmp_obj[[r]]$rep$SB[p2]  * mean_rec[p2]
              virgin_pop_b_ref_pt[p2,r] <- tmp_obj[[r]]$rep$SB0[p2] * mean_rec[p2]
              for(p in 1:n_pop) {
                if(p != p2) {
                  sc <- data_list$stray_rate[p] / n_pop_in_region[rn]
                  pop_b_ref_pt[p2,r]        <- pop_b_ref_pt[p2,r]        + sc * tmp_obj[[r]]$rep$SB[p]  * mean_rec[p]
                  virgin_pop_b_ref_pt[p2,r] <- virgin_pop_b_ref_pt[p2,r] + sc * tmp_obj[[r]]$rep$SB0[p] * mean_rec[p]
                } # end if
              } # end p loop
            } # end p2 loop
          }

          # Compute global reference points (sum across populations)
          b_ref_pt[,r] <- tmp_obj[[r]]$rep$SB * apply(rep$Rec[,r,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)
          virgin_b_ref_pt[,r] <- tmp_obj[[r]]$rep$SB0 * apply(rep$Rec[,r,calc_rec_st_yr:(n_years - rec_age), drop = F], 1, mean)

        } # independent SPR

        if(what == 'independent_BH_MSY') {

          # Beverton Holt parameters
          data_list$h <- array(rep$h_trans[,r], dim = n_pop) # steepness
          data_list$R0 <- array(rep$R0 * rep$rec_region_prop[,r], dim = n_pop) # unfished recruitment by region

          par_list <- list() # set up parameter list
          par_list$log_Fmsy <- log(0.1) # Fmsy starting value

          # optimize model
          tmp_obj[[r]] <- optim_ref_pts(single_region_BH_Fmsy, data_list, par_list)
          f_ref_pt[r] <- tmp_obj[[r]]$rep$Fmsy

          # get and accumulate biomass reference points
          if(n_pop > 1) {
            for(p2 in 1:n_pop) {
              rn <- data$natal_region[p2]
              pop_b_ref_pt[p2,r]        <- tmp_obj[[r]]$rep$SB[p2]  * tmp_obj[[r]]$rep$Req[p2]
              virgin_pop_b_ref_pt[p2,r] <- tmp_obj[[r]]$rep$SB0[p2] * rep$R0[p2]
              for(p in 1:n_pop) {
                if(p != p2) {
                  sc <- data_list$stray_rate[p] / n_pop_in_region[rn]
                  pop_b_ref_pt[p2,r]        <- pop_b_ref_pt[p2,r]        + sc * tmp_obj[[r]]$rep$SB[p]  * tmp_obj[[r]]$rep$Req[p]
                  virgin_pop_b_ref_pt[p2,r] <- virgin_pop_b_ref_pt[p2,r] + sc * tmp_obj[[r]]$rep$SB0[p] * data_list$R0[p]
                }
              }
            }
          }

          b_ref_pt[,r] <-  tmp_obj[[r]]$rep$SB * tmp_obj[[r]]$rep$Req
          virgin_b_ref_pt[,r] <-  tmp_obj[[r]]$rep$SB0 * data_list$R0

        } # independent BH MSY

      } # end r loop

      # sum up biomass reference points to pop-specific quantities
      if(n_pop == 1) {
        pop_b_ref_pt[1,1] = sum(b_ref_pt)
        virgin_pop_b_ref_pt[1,1] = sum(virgin_b_ref_pt)
      }

    } # end independent methods

    # Global SPR
    if(what == 'global_SPR') {

      # create data lists
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
      natmort_avg <- apply(rep$natmort[,,avg_yrs,,1,drop = FALSE], c(1,2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_pop, n_regions, n_ages)) # get female natural mortality
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_pop, n_regions, data$n_seas, n_ages)) # weight at age for females
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_pop, n_regions, data$n_seas, n_ages)) # maturity at age for females
      Movement_avg <- apply(rep$Movement[,,,avg_yrs,,,1,drop = FALSE], c(1,2,3,5,6), mean)
      data_list$Movement <- array(Movement_avg, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)) # Movement

      # Recruitment options
      data_list$do_recruits_move <- data$do_recruits_move # whether recruits move
      data_list$sex_ratio_f <- sex_ratio_f # recruitment sex ratio
      data_list$rec_seas_prop <- array(rep$rec_seas_prop[,], dim = c(n_pop, data$n_seas)) # recruitment seasonal proportion
      data_list$stray_rate <- array(apply(rep$stray_rate[,avg_yrs, drop = FALSE], 1, mean), dim = n_pop) # stray rate
      sgl_seas_spawning_movement_avg <- apply(rep$sgl_seas_spawning_movement[,,,avg_yrs,,1,drop = FALSE], c(1,2,3,5), mean)
      data_list$sgl_seas_spawning_movement <- array(sgl_seas_spawning_movement_avg, dim = c(n_pop, n_regions, n_regions, n_ages)) # Movement
      data_list$natal_region <- data$natal_region
      data_list$n_pop_in_region <- n_pop_in_region

      data_list$SPR_x <- SPR_x # SPR fraction
      mean_rec <- apply(rep$Rec[,,calc_rec_st_yr:(n_years-rec_age),drop=FALSE], c(1,2), mean) # [n_pop, n_regions]
      total_mean_rec <- apply(mean_rec, 1, sum) # [n_pop] - total recruitment across regions
      data_list$rec_region_prop <- mean_rec / total_mean_rec # recruitment proportions

      par_list <- list() # set up parameter list
      par_list$log_F_x <- log(0.1) # F_x starting value

      # make adfn object
      tmp_obj <- optim_ref_pts(global_SPR, data_list, par_list)

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
          r_natal <- data$natal_region[p2]
          pop_b_ref_pt[p2, r_natal]        <- tmp_obj$rep$SB[p2, r_natal]  * total_mean_rec[p2]
          virgin_pop_b_ref_pt[p2, r_natal] <- tmp_obj$rep$SB0[p2, r_natal] * total_mean_rec[p2]
          for(p in 1:n_pop) {
            if(p != p2) {
              sc <- data_list$stray_rate[p] / n_pop_in_region[r_natal]
              pop_b_ref_pt[p2, r_natal]        <- pop_b_ref_pt[p2, r_natal]        + sc * tmp_obj$rep$SB[p, r_natal]  * total_mean_rec[p]
              virgin_pop_b_ref_pt[p2, r_natal] <- virgin_pop_b_ref_pt[p2, r_natal] + sc * tmp_obj$rep$SB0[p, r_natal] * total_mean_rec[p]
            }
          }
        }
      } else {
        pop_b_ref_pt[1,1] = sum(b_ref_pt)
        virgin_pop_b_ref_pt[1,1] = sum(virgin_b_ref_pt)
      }

    } # end global SPR

    # Global BH MSY
    if(what == 'global_BH_MSY') {

      # Error out if invalid recruitment density dependent option
      if(n_pop > 1) stop("Invalid reference point option! When n_pop > 1 reference points must either be independent_SPR, independent_BH_MSY, global_SPR, or local_BH_MSY.")

      # create a data list
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
      natmort_avg <- apply(rep$natmort[1,,avg_yrs,,1,drop = FALSE], c(2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(n_regions, n_ages)) # get female natural mortality
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # weight at age for females
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(n_regions, data$n_seas, n_ages)) # maturity at age for females
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

      # make adfn object
      tmp_obj <- optim_ref_pts(global_BH_Fmsy, data_list, par_list)

      # Output reference points
      f_ref_pt <- rep(tmp_obj$rep$Fmsy, n_regions)
      b_ref_pt[1,] <- apply(tmp_obj$rep$SB_age[2,,,drop = F], 2, sum) * tmp_obj$rep$Req
      pop_b_ref_pt[1,1] <- sum(b_ref_pt)
      virgin_b_ref_pt[1,] <- apply(tmp_obj$rep$SB_age[1,,,drop = F], 2, sum) * data_list$R0
      virgin_pop_b_ref_pt[1,1]  <- sum(virgin_b_ref_pt)

    } # end global BH MSY

    if(what == 'local_BH_MSY') {

      # create data list
      fratio <- array(0, dim = c(n_regions, data$n_seas, data$n_fish_fleets))
      terminal_F <- array(rep$Fmort[,n_years,,], dim = dim(fratio))
      for(r in 1:n_regions) for(seas in 1:data$n_seas) for(f in 1:data$n_fish_fleets) fratio[r,seas,f] <- terminal_F[r,seas,f] / sum(terminal_F[r,,])
      data_list$F_fract_flt <- fratio
      fish_sel_avg <- apply(rep$fish_sel[,avg_yrs,,1,,drop = FALSE], c(1,3,5), mean)
      data_list$fish_sel <- array(fish_sel_avg, dim = c(n_regions, n_ages, data$n_fish_fleets)) # get female selectivity for all fleets
      natmort_avg <- apply(rep$natmort[,,avg_yrs,,1,drop = FALSE], c(1,2,4), mean)
      data_list$natmort <- array(natmort_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, n_ages)) # get female natural mortality
      WAA_avg <- apply(data$WAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$WAA <- array(WAA_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, data$n_seas, n_ages)) # weight at age for females
      MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,1,drop = FALSE], c(1, 2, 4, 5), mean)
      data_list$MatAA <- array(MatAA_avg, dim = c(if(n_pop > 1) n_pop else NULL, n_regions, data$n_seas, n_ages)) # maturity at age for females
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
      data_list$n_pop_in_region <- n_pop_in_region
      sgl_seas_spawning_movement_avg <- apply(rep$sgl_seas_spawning_movement[,,,avg_yrs,,1,drop = FALSE], c(1,2,3,5), mean)
      data_list$sgl_seas_spawning_movement <- array(sgl_seas_spawning_movement_avg, dim = c(n_pop, n_regions, n_regions, n_ages)) # Movement

      par_list <- list() # set up parameter list
      par_list$log_Fmsy <- rep(log(0.1), n_regions) # Fmsy starting value

      # Make adfun object
      tmp_obj <- optim_ref_pts(if(n_pop == 1) local_BH_Fmsy_sglpop else local_BH_Fmsy_multipop, data_list, par_list)

      # Output reference points
      f_ref_pt <- tmp_obj$rep$Fmsy

      # multi-population reference points
      if(n_pop > 1) {
        for(r in 1:n_regions) {
          b_ref_pt[,r]        <- tmp_obj$rep$SB[,r]  * tmp_obj$rep$Req_o
          virgin_b_ref_pt[,r] <- tmp_obj$rep$SB0[,r] * data_list$R0
        }

        # accumulate stray rates (divided by n_pop_in_region for mass balance)
        for(p2 in 1:n_pop) {
          r_natal <- data$natal_region[p2]
          pop_b_ref_pt[p2, r_natal]        <- tmp_obj$rep$SB[p2, r_natal]  * tmp_obj$rep$Req_o[p2]
          virgin_pop_b_ref_pt[p2, r_natal] <- tmp_obj$rep$SB0[p2, r_natal] * data_list$R0[p2]

          for(p in 1:n_pop) {
            if(p != p2) {
              sc <- data_list$stray_rate[p] / n_pop_in_region[r_natal]
              pop_b_ref_pt[p2, r_natal]        <- pop_b_ref_pt[p2, r_natal]        + sc * tmp_obj$rep$SB[p, r_natal]  * tmp_obj$rep$Req_o[p]
              virgin_pop_b_ref_pt[p2, r_natal] <- virgin_pop_b_ref_pt[p2, r_natal] + sc * tmp_obj$rep$SB0[p, r_natal] * data_list$R0[p]
            }
          }
        }
      } else { # single population reference points
        for(r in 1:n_regions) {
          b_ref_pt[1, r]        <- sum(tmp_obj$rep$SB_fished_mat[, r]   * tmp_obj$rep$Req_o)
          virgin_b_ref_pt[1, r] <- sum(tmp_obj$rep$SB_unfished_mat[, r] * data_list$R0 * data_list$rec_region_prop)
        }
        pop_b_ref_pt[1, 1]        <- sum(b_ref_pt)
        virgin_pop_b_ref_pt[1, 1] <- sum(virgin_b_ref_pt)
      }
      # see if Newton Raphson calcs for equil rec converged
      if(sum(tmp_obj$rep$iter_vec) > 1e-10) warning("Calculations for equilibrium recruits from origin regions might not have converged! Try increasing local_bh_msy_newton_steps or be wary of these values!")
      if(sum(tmp_obj$rep$Fmsy) == sum(exp(par_list$log_Fmsy))) warning("It is unlikely this converged. Starting values of log Fmsy have not changed (specified at log (0.1).")

    }

  } # end multi region

  return(list(f_ref_pt = f_ref_pt,
              b_ref_pt = b_ref_pt,
              virgin_b_ref_pt = virgin_b_ref_pt,
              pop_b_ref_pt = pop_b_ref_pt,
              virgin_pop_b_ref_pt = virgin_pop_b_ref_pt,
              obj = tmp_obj))

}
