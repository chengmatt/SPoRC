# Stage 2 of 3: objective function
#
# Deterministic recruitment: mean recruitment or Beverton Holt, apportioned
# across regions and seasons. When Beverton Holt is used, unfished spawning
# biomass per recruit is solved here, which is why this file steps through the
# seasonal operators in model_transition.R.

#' Deterministic Recruitment
#'
#' Computes deterministic recruitment by population and region using either
#' a mean recruitment model or a Beverton-Holt stock-recruitment relationship.
#'
#' Recruitment is distributed spatially using regional recruitment proportions
#' and seasonal recruitment timing. When Beverton-Holt recruitment is used,
#' unfished spawning biomass per recruit (\eqn{S_0}) is calculated internally
#' by projecting a single recruit through the full seasonal population dynamics,
#' including movement and mortality.
#'
#' @param recruitment_model Integer flag specifying the recruitment model:
#'   \itemize{
#'   \item \code{0} Mean recruitment
#'   \item \code{1} Beverton-Holt recruitment with steepness
#'   \item \code{2} Ricker recruitment with steepness
#'   }
#'
#' @param rec_dd Integer flag specifying the density dependence structure:
#'   \itemize{
#'   \item \code{0} Local density dependence (population or region specific)
#'   \item \code{1} Global density dependence (shared across regions;
#'   only valid when \code{n_pop = 1})
#'   }
#'
#' @param y Current model year index.
#' @param rec_lag Recruitment lag (in seasons) between spawning and
#'   recruitment. \code{1} is the classic lagged case: recruitment uses
#'   \code{SSB_vals} from \code{rec_lag} seasons prior. \code{0} is age-0
#'   recruitment: recruitment uses the SAME year's SSB
#'   (\code{SSB_vals[,,y]}). The caller is responsible for supplying that
#'   value already computed from survivors only (i.e. before this year's
#'   recruits exist) when \code{rec_lag = 0}, see \code{SPoRC_rtmb.R},
#'   \code{Simulate_Population.R}, and \code{Do_Population_Projection.R} for
#'   how each population-dynamics loop does this.
#' @param R0 Numeric vector (\code{n_pop}) of unfished recruitment by population.
#' @param rec_region_prop Matrix (\code{n_pop × n_regions}) giving the proportion
#'   of recruitment allocated to each region.
#' @param rec_seas_prop Matrix (\code{n_pop × n_seas}) giving seasonal recruitment
#'   proportions. When \code{rec_lag = 0}, must be zero for every season
#'   before \code{spawn_seas} (age-0 recruits can't predate the spawning
#'   event that produced them), validated at setup by
#'   \code{Setup_Mod_Rec}/\code{Setup_Sim_Rec}.
#' @param h Matrix (\code{n_pop × n_regions}) of Beverton-Holt steepness values.
#' @param n_pop Number of populations.
#' @param n_regions Number of spatial regions.
#' @param n_ages Number of age classes (including the plus group).
#' @param WAA Array (\code{n_pop × n_regions × n_seas × n_ages}) of weight-at-age.
#' @param MatAA Array (\code{n_pop × n_regions × n_seas × n_ages}) of maturity-at-age.
#' @param natmort Array (\code{n_pop × n_regions × n_ages}) of natural mortality.
#' @param SSB_vals Array (\code{n_pop × n_regions × n_years}) of spawning biomass.
#' @param Movement Array
#'   (\code{n_pop × origin × destination × n_seas × n_ages}) giving seasonal
#'   movement probabilities.
#' @param sgl_seas_spawning_movement Array
#'   (\code{n_pop × origin × destination × n_ages}) describing spawning movement
#'   when a single season is used and \code{n_pop > 1}.
#' @param stray_rate Numeric vector of stray rates by population.
#' @param do_recruits_move Indicator for whether recruits move in their first year.
#' @param t_spawn Fraction of the spawning season that occurs before spawning.
#' @param init_F Array (\code{n_regions × n_seas × n_fish_fleets}) of initial fishing mortality.
#' @param fish_sel Array (\code{n_pop x n_regions × n_seas x n_ages x n_fish_fleets}) of total fishery selectivity.
#' @param n_seas Number of seasons per year.
#' @param spawn_seas Season index in which spawning occurs.
#' @param natal_region Integer vector (\code{n_pop}) mapping each population
#'   to its natal region.
#' @param seasdur Numeric vector (\code{n_seas}) giving seasonal durations
#'   as fractions of a year.
#' @param sexratio_f Matrix (\code{n_pop × n_regions}) giving female recruitment
#'   proportions.
#' @param n_fish_fleets Integer. Number of fishery fleets.
#' @param dmr Array (\code{n_regions × n_seas × n_fish_fleets}) of initial (first year) discard mortality.
#' @param ret_sel Array (\code{n_pop x n_regions × n_seas x n_ages x n_fish_fleets}) of retained fishery selectivity.
#'
#' @details
#'
#' Two recruitment formulations are supported.
#'
#' **Mean recruitment**
#'
#' When \code{recruitment_model = 0}, recruitment is constant:
#'
#' \deqn{R_{p,r} = R_{0,p} \times RecProp_{p,r}}
#'
#' where recruitment is distributed spatially according to
#' \code{rec_region_prop}.
#'
#' **Beverton-Holt recruitment**
#'
#' When \code{recruitment_model = 1}, recruitment follows the
#' Beverton-Holt relationship:
#'
#' \deqn{
#' R = \frac{4hR_0SSB}{(1-h)S_0 + (5h-1)SSB}
#' }
#'
#' where:
#' \itemize{
#' \item \eqn{SSB} is spawning biomass lagged by \code{rec_lag} seasons
#'   (or, when \code{rec_lag = 0}, the current year's own spawning biomass,
#'   see the \code{rec_lag} parameter above)
#' \item \eqn{S_0} is unfished spawning biomass per recruit
#' \item \eqn{h} is steepness
#' }
#'
#' **Ricker recruitment**
#'
#' When \code{recruitment_model = 2}, recruitment follows the Ricker
#' relationship, written in the depletion form used by the EBS pollock
#' assessment (2024, \code{SrType = 1}):
#'
#' \deqn{
#' R = R_0 \frac{SSB}{S_0} \exp\left(\alpha \left(1 - \frac{SSB}{S_0}\right)\right),
#' \quad \alpha = \log\left(\frac{4h}{1-h}\right)
#' }
#'
#' The curve passes through \eqn{(S_0, R_0)} by construction. Note that
#' \eqn{\alpha} is set so the Ricker carries the same compensation ratio as a
#' Beverton-Holt at the same \eqn{h}, rather than by the textbook definition
#' \eqn{R(0.2 S_0) = h R_0}. Steepness is therefore not interchangeable between
#' the two curves: the Ricker here gives \eqn{R(0.2S_0)/R_0 = 0.2(4h/(1-h))^{0.8}},
#' which exceeds \eqn{h} and is not bounded by 1.
#'
#' \eqn{S_0} (and the age-composition of spawning biomass per recruit more
#' generally) does not depend on \code{rec_lag}; it is a pure per-recruit,
#' equilibrium quantity. The recruit age class (the first age) is always
#' included in the sum; when \code{rec_lag = 0}, maturity at that age is
#' required to be exactly zero (validated at setup by
#' \code{Setup_Mod_Biologicals}/\code{Setup_Sim_Biologicals}), so it
#' contributes nothing regardless.
#'
#' Spawning biomass per recruit (\eqn{S_0}) is computed internally by
#' projecting a single recruit through all ages and seasons under both
#' unfished and fished conditions. The algorithm:
#'
#' \enumerate{
#' \item Allocates a recruit across regions and seasons.
#' \item Applies seasonal movement.
#' \item Applies natural and fishing mortality, where fishing mortality
#'   is decomposed into retained
#'   (\eqn{F \cdot sel \cdot ret}) and dead discard
#'   (\eqn{F \cdot sel \cdot (1 - ret) \cdot dmr}) components.
#' \item Computes spawning biomass during the spawning season.
#' \item Solves the plus group analytically using annual transition matrices.
#' }
#'
#' When multiple populations are modeled, recruitment for each population
#' depends on spawning biomass in its natal region. Contributions from
#' other populations are scaled by the specified stray rates.
#'
#'
#' @keywords internal
Get_Det_Recruitment <- function(recruitment_model,
                                rec_dd,
                                y,
                                rec_lag,
                                R0,
                                rec_region_prop,
                                rec_seas_prop,
                                h,
                                n_pop,
                                n_regions,
                                n_ages,
                                n_fish_fleets,
                                WAA,
                                MatAA,
                                natmort,
                                SSB_vals,
                                Movement,
                                sgl_seas_spawning_movement,
                                stray_rate,
                                do_recruits_move,
                                t_spawn,
                                init_F,
                                dmr,
                                fish_sel,
                                ret_sel,
                                n_seas,
                                spawn_seas,
                                natal_region,
                                seasdur,
                                sexratio_f,
                                Mrate = NULL,
                                move_timing = 0
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(recruitment_model == 0) {
    rec = array(0, dim = c(n_pop, n_regions))
    for(p in 1:n_pop) rec[p,] = R0[p] * rec_region_prop[p,] # mean recruitment apportioned across n_pop and n_regions
  }

  # Beverton-Holt (recruitment_model == 1) and Ricker (recruitment_model == 2).
  # Both are density-dependent forms driven by unfished spawning biomass, so they
  # share the whole spawning-biomass-per-recruit calculation and differ only in the
  # curve evaluated at the end.
  if(recruitment_model %in% c(1, 2)) {

    # Storage for recruitment, S0, and SF (equilibrium fished)
    rec = S0 = SF = array(0, dim = c(n_pop, n_regions))

    # Ricker log-slope, derived from steepness as log(4h / (1 - h)) - Dorn parameterization
    ricker_alpha = function(hh) log(4 * hh / (1 - hh))

    # local density dependence
    if(rec_dd == 0) {

      # Calculate unexploited naa per recruit by origin area and destination area
      SB_fished_age = Nspr_fished = SB_age = Nspr = array(0, dim = c(n_pop, n_regions, n_regions, n_ages))
      SB_fished_mat = SB_unfished_mat = array(0, c(n_pop, n_regions, n_regions))

      # Set up the initial recruits (1 recruit per area)
      for(p in 1:n_pop) {
        for(o in 1:n_regions) {
          for(d in 1:n_regions) {

            if(o == d) Nspr_fished[p,o,d,1] = Nspr[p,o,d,1] = sexratio_f[p,o] * rec_region_prop[p,o] * rec_seas_prop[p,1] # apportion recruits in the first season
            else Nspr_fished[p,o,d,1] = Nspr[p,o,d,1] = 0

          } # end d loop
        } # end o loop
      } # end p loop

      # Loop through ages, projecting each cohort through the full annual cycle
      for(j in 2:(n_ages-1)){

        # Project age j-1 through all seasons to become age j
        for(seas in 1:n_seas) {

          for(p in 1:n_pop) {
            for(o in 1:n_regions) {

              # Get temporary values from origin region
              tmp_unfished = Nspr[p,o,,j-1]
              tmp_fished = Nspr_fished[p,o,,j-1]

              # add in seasonal recruits
              if(seas > 1 && j - 1 == 1) {
                tmp_unfished[o] = tmp_unfished[o] + rec_seas_prop[p,seas] * sexratio_f[p,o] * rec_region_prop[p,o]
                tmp_fished[o]   = tmp_fished[o]   + rec_seas_prop[p,seas] * sexratio_f[p,o] * rec_region_prop[p,o]
              }

              # Movement operators for this age; non-moving ages get an identity transition
              # and a zero generator, which leaves survival unchanged under every move_timing
              if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
                Mv = Movement[p,,,seas,j-1]
                Qv = Mrate[p,,,seas,j-1]
              } else {
                Mv = diag(n_regions)
                Qv = matrix(0, n_regions, n_regions)
              }

              # Total mortality for this season, by region
              Zu_seas = natmort[p,,j-1] * seasdur[seas]
              Zf_seas = Zu_seas + rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,j-1,] * ret_sel[p,,seas,j-1,] +
                                                                    fish_sel[p,,seas,j-1,] * (1 - ret_sel[p,,seas,j-1,]) * dmr[,seas,]),
                                                dim = c(n_regions, n_fish_fleets)))

              # Calculate spawning biomass if this is the spawning season
              if(seas == spawn_seas) {

                # Propagate to the spawning point; movement and t_spawn mortality applied together
                tmp_unfished_spawn = spawn_state(tmp_unfished, Mv, Zu_seas, Qv, seasdur[seas], t_spawn, move_timing)
                tmp_fished_spawn = spawn_state(tmp_fished, Mv, Zf_seas, Qv, seasdur[seas], t_spawn, move_timing)

                # If single season natal homing population
                if(n_pop > 1 && n_seas == 1) {
                  # Get NAA during spawning in single season case
                  tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
                  tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,j-1]
                }

                # Get spawning biomass per recruit by age
                for(d in 1:n_regions) {
                  SB_age[p,o,d,j-1] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1]
                  SB_fished_age[p,o,d,j-1] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1]
                  }
              }

              # Apply movement and mortality together, per move_timing
              adv_unfished = advance_seas(tmp_unfished, Mv, Zu_seas, Qv, seasdur[seas], move_timing)
              adv_fished   = advance_seas(tmp_fished, Mv, Zf_seas, Qv, seasdur[seas], move_timing)

              if(seas < n_seas) {
                # Within-season mortality, no ageing yet
                Nspr[p,o,,j-1] = adv_unfished
                Nspr_fished[p,o,,j-1] = adv_fished
                } else {
                # Last season: mortality + ageing
                Nspr[p,o,,j] = adv_unfished
                Nspr_fished[p,o,,j] = adv_fished
                }

            } # end o loop
          } # end p loop

        } # end seas loop
      } # end j loop

      # Now calculate spawning biomass for penultimate age (n_ages-1)
      for(p in 1:n_pop) {
        for(o in 1:n_regions) {

          # Age n_ages-1 is now at start of year after the loop
          tmp_unfished = Nspr[p,o,,n_ages-1]
          tmp_fished = Nspr_fished[p,o,,n_ages-1]

          if(spawn_seas > 1) {
            for (seas in 1:(spawn_seas - 1)) {

              # Apply seasonal movement and mortality together, per move_timing
              Zu_seas = natmort[p,,n_ages-1] * seasdur[seas]
              Zf_seas = Zu_seas + rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages-1,] * ret_sel[p,,seas,n_ages-1,] +
                                                                    fish_sel[p,,seas,n_ages-1,] * (1 - ret_sel[p,,seas,n_ages-1,]) * dmr[,seas,]),
                                                dim = c(n_regions, n_fish_fleets)))
              tmp_unfished = advance_seas(tmp_unfished, Movement[p,,,seas,n_ages-1], Zu_seas,
                                          Mrate[p,,,seas,n_ages-1], seasdur[seas], move_timing)
              tmp_fished = advance_seas(tmp_fished, Movement[p,,,seas,n_ages-1], Zf_seas,
                                        Mrate[p,,,seas,n_ages-1], seasdur[seas], move_timing)
            } # end seas loop
          }

          # Propagate into the spawning season; t_spawn mortality folded in below
          Zu_spawn = natmort[p,,n_ages-1] * seasdur[spawn_seas]
          Zf_spawn = Zu_spawn + rowSums(array(init_F[,spawn_seas,] * (fish_sel[p,,spawn_seas,n_ages-1,] * ret_sel[p,,spawn_seas,n_ages-1,] +
                                                                        fish_sel[p,,spawn_seas,n_ages-1,] * (1 - ret_sel[p,,spawn_seas,n_ages-1,]) * dmr[,spawn_seas,]),
                                              dim = c(n_regions, n_fish_fleets)))
          tmp_unfished_spawn = spawn_state(tmp_unfished, Movement[p,,,spawn_seas,n_ages-1], Zu_spawn,
                                           Mrate[p,,,spawn_seas,n_ages-1], seasdur[spawn_seas], t_spawn, move_timing)
          tmp_fished_spawn = spawn_state(tmp_fished, Movement[p,,,spawn_seas,n_ages-1], Zf_spawn,
                                         Mrate[p,,,spawn_seas,n_ages-1], seasdur[spawn_seas], t_spawn, move_timing)

          # If single season natal homing population
          if(n_pop > 1 && n_seas == 1) {
            # Get NAA during spawning in single season case
            tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages-1]
            tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages-1]
          }

          # Calculate spawning biomass
          for(d in 1:n_regions) {
            SB_age[p,o,d,n_ages-1] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,n_ages-1] * MatAA[p,d,spawn_seas,n_ages-1]
            SB_fished_age[p,o,d,n_ages-1] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,n_ages-1] * MatAA[p,d,spawn_seas,n_ages-1]
            }
        }
      }

      # Set up analytical solution for plus group
      for(p in 1:n_pop) {

        # Build FULL annual transition for penultimate and plus ages
        T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

        # Loop through ALL seasons to build annual transition matrix
        for(seas in 1:n_seas) {
          # Fished mortality components
          F_penult = rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages-1,] * ret_sel[p,,seas,n_ages-1,] +
                             fish_sel[p,,seas,n_ages-1,] * (1 - ret_sel[p,,seas,n_ages-1,]) * dmr[,seas,]),
                                   dim = c(n_regions, n_fish_fleets)))
          F_plus = rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages,] * ret_sel[p,,seas,n_ages,] +
                           fish_sel[p,,seas,n_ages,] * (1 - ret_sel[p,,seas,n_ages,]) * dmr[,seas,]),
                                 dim = c(n_regions, n_fish_fleets)))

          # Column-convention seasonal operators (build_seas_operator returns row convention),
          # composed left-to-right so that season 1 is applied first
          op <- function(age, Z) t(build_seas_operator(Movement[p,,,seas,age], Z, Mrate[p,,,seas,age], seasdur[seas], move_timing))
          T_penult_unfished = op(n_ages-1, natmort[p,,n_ages-1] * seasdur[seas]) %*% T_penult_unfished
          T_plus_unfished   = op(n_ages,   natmort[p,,n_ages]   * seasdur[seas]) %*% T_plus_unfished
          T_penult_fished   = op(n_ages-1, natmort[p,,n_ages-1] * seasdur[seas] + F_penult) %*% T_penult_fished
          T_plus_fished     = op(n_ages,   natmort[p,,n_ages]   * seasdur[seas] + F_plus) %*% T_plus_fished
        } # end seas loop

        for(o in 1:n_regions) {

          # Solve for equilibrium plus group (at start of year)
          source_unfished = T_penult_unfished %*% Nspr[p,o,,n_ages-1]
          Nspr[p,o,,n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
          source_fished = T_penult_fished %*% Nspr_fished[p,o,,n_ages-1]
          Nspr_fished[p,o,,n_ages] = solve(diag(n_regions) - T_plus_fished, source_fished)

        } # end o loop
      } # end p loop

      # Now calculate spawning biomass for plus age (n_ages)
      for(p in 1:n_pop) {
        for(o in 1:n_regions) {

          # Age n_ages-1 is now at start of year after the loop
          tmp_unfished = Nspr[p,o,,n_ages]
          tmp_fished = Nspr_fished[p,o,,n_ages]

          if(spawn_seas > 1) {
            for (seas in 1:(spawn_seas - 1)) {

              # Apply seasonal movement and mortality together, per move_timing
              Zu_seas = natmort[p,,n_ages] * seasdur[seas]
              Zf_seas = Zu_seas + rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages,] * ret_sel[p,,seas,n_ages,] +
                                                                    fish_sel[p,,seas,n_ages,] * (1 - ret_sel[p,,seas,n_ages,]) * dmr[,seas,]),
                                                dim = c(n_regions, n_fish_fleets)))
              tmp_unfished = advance_seas(tmp_unfished, Movement[p,,,seas,n_ages], Zu_seas,
                                          Mrate[p,,,seas,n_ages], seasdur[seas], move_timing)
              tmp_fished = advance_seas(tmp_fished, Movement[p,,,seas,n_ages], Zf_seas,
                                        Mrate[p,,,seas,n_ages], seasdur[seas], move_timing)

            } # end seas loop
          }

          # Propagate into the spawning season; t_spawn mortality folded in below
          Zu_spawn = natmort[p,,n_ages] * seasdur[spawn_seas]
          Zf_spawn = Zu_spawn + rowSums(array(init_F[,spawn_seas,] * (fish_sel[p,,spawn_seas,n_ages,] * ret_sel[p,,spawn_seas,n_ages,] +
                                                                        fish_sel[p,,spawn_seas,n_ages,] * (1 - ret_sel[p,,spawn_seas,n_ages,]) * dmr[,spawn_seas,]),
                                              dim = c(n_regions, n_fish_fleets)))
          tmp_unfished_spawn = spawn_state(tmp_unfished, Movement[p,,,spawn_seas,n_ages], Zu_spawn,
                                           Mrate[p,,,spawn_seas,n_ages], seasdur[spawn_seas], t_spawn, move_timing)
          tmp_fished_spawn = spawn_state(tmp_fished, Movement[p,,,spawn_seas,n_ages], Zf_spawn,
                                         Mrate[p,,,spawn_seas,n_ages], seasdur[spawn_seas], t_spawn, move_timing)

          # If single season natal homing population
          if(n_pop > 1 && n_seas == 1) {
            # Get NAA during spawning in single season case
            tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages]
            tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages]
          }

          # Calculate spawning biomass
          for(d in 1:n_regions) {
            SB_age[p,o,d,n_ages] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,n_ages] * MatAA[p,d,spawn_seas,n_ages]
            SB_fished_age[p,o,d,n_ages] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,n_ages] * MatAA[p,d,spawn_seas,n_ages]
            } # end d loop

        } # end o loop
      } # end p loop

      # Remove the old spawning biomass calculation loop entirely
      # parse out and compute unfished spawning biomass per recruit
      for(p in 1:n_pop) {
        for(o in 1:n_regions) {
          for(d in 1:n_regions) {
            # unfished
            SB_unfished_mat[p, o, d] = sum(SB_age[p, o, d, 1:n_ages])
            # fished
            SB_fished_mat[p, o, d] = sum(SB_fished_age[p, o, d, 1:n_ages])
          } # end d
        } # end o
      } # end p

      for(p in 1:n_pop) {
        for(d in 1:n_regions) {
          S0[p,d] = sum(SB_unfished_mat[p,,d] * R0[p]) # unfished
          SF[p,d] = sum(SB_fished_mat[p,,d] * R0[p]) # fished
        } # end d
      } # end p loop

    } # end if rec_dd == 0

    # global density dependence
    if (rec_dd == 1) {

      # Error out if invalid recruitment density dependent option
      if(n_pop > 1) stop("Invalid recruitment density-dependence option! When n_pop > 1 rec_dd must be local (0).")

      # Setup containers
      SB_fished_age = SB_age = Nspr_fished = Nspr = array(0, dim = c(n_regions, n_ages))

      # Initial recruits: 1 recruit globally, split by rec_region_prop and seasonal recruitment
      Nspr[,1] = rec_region_prop[1,] * sexratio_f[1,] * rec_seas_prop[1,1]
      Nspr_fished[,1] = rec_region_prop[1,] * sexratio_f[1,] * rec_seas_prop[1,1]

      ## Loop through ages
      for (j in 2:(n_ages - 1)) {
        for (seas in 1:n_seas) {

          tmp_unfished = Nspr[, j - 1]
          tmp_fished   = Nspr_fished[, j - 1]

          # apportion seasonal recruits
          if(seas > 1 && j - 1 == 1) {
            tmp_unfished = tmp_unfished + rec_seas_prop[1,seas] * sexratio_f[1,] * rec_region_prop[1,]
            tmp_fished = tmp_fished + rec_seas_prop[1,seas] * sexratio_f[1,] * rec_region_prop[1,]
          }

          ## Movement operators for this age; non-moving ages get an identity transition
          ## and a zero generator, leaving survival unchanged under every move_timing
          if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
            Mv = Movement[1,,,seas, j - 1]
            Qv = Mrate[1,,,seas, j - 1]
          } else {
            Mv = diag(n_regions)
            Qv = matrix(0, n_regions, n_regions)
          }

          # Total mortality for this season, by region
          Zu_seas = natmort[1,, j - 1] * seasdur[seas]
          Zf_seas = Zu_seas + rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,j-1,] * ret_sel[1,,seas,j-1,] +
                                                                fish_sel[1,,seas,j-1,] * (1 - ret_sel[1,,seas,j-1,]) * dmr[,seas,]),
                                            dim = c(n_regions, n_fish_fleets)))

          ## Spawning biomass
          if (seas == spawn_seas) {
            tmp_unfished_spawn = spawn_state(tmp_unfished, Mv, Zu_seas, Qv, seasdur[seas], t_spawn, move_timing)
            tmp_fished_spawn = spawn_state(tmp_fished, Mv, Zf_seas, Qv, seasdur[seas], t_spawn, move_timing)
            SB_age[, j - 1] = tmp_unfished_spawn * WAA[1,, spawn_seas, j - 1] * MatAA[1,, spawn_seas, j - 1]
            SB_fished_age[, j - 1] = tmp_fished_spawn * WAA[1,, spawn_seas, j - 1] * MatAA[1,, spawn_seas, j - 1]
            }

          ## Movement, mortality and ageing
          adv_unfished = advance_seas(tmp_unfished, Mv, Zu_seas, Qv, seasdur[seas], move_timing)
          adv_fished   = advance_seas(tmp_fished, Mv, Zf_seas, Qv, seasdur[seas], move_timing)

          if (seas < n_seas) { # Within season mortality
            Nspr[, j - 1] = adv_unfished
            Nspr_fished[, j - 1] = adv_fished
            } else {
            # Ageing
            Nspr[, j] = adv_unfished
            Nspr_fished[, j] = adv_fished
          }
        } # end seas loop
      } # end j loop

      # Age n_ages-1 is now at start of year after the loop
      tmp_unfished = Nspr[,n_ages-1]
      tmp_fished = Nspr_fished[,n_ages-1]

      if(spawn_seas > 1) {
        for (seas in 1:(spawn_seas - 1)) {

          # Apply seasonal movement and mortality together, per move_timing
          Zu_seas = natmort[1,,n_ages-1] * seasdur[seas]
          Zf_seas = Zu_seas + rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages-1,] * ret_sel[1,,seas,n_ages-1,] +
                                                                fish_sel[1,,seas,n_ages-1,] * (1 - ret_sel[1,,seas,n_ages-1,]) * dmr[,seas,]),
                                            dim = c(n_regions, n_fish_fleets)))
          tmp_unfished = advance_seas(tmp_unfished, Movement[1,,,seas,n_ages-1], Zu_seas,
                                      Mrate[1,,,seas,n_ages-1], seasdur[seas], move_timing)
          tmp_fished = advance_seas(tmp_fished, Movement[1,,,seas,n_ages-1], Zf_seas,
                                    Mrate[1,,,seas,n_ages-1], seasdur[seas], move_timing)

        } # end seas loop
      } # end if spawn_seas > 1

      ## Penultimate age spawning biomass; t_spawn mortality folded into spawn_state
      Zu_spawn = natmort[1,, n_ages - 1] * seasdur[spawn_seas]
      Zf_spawn = Zu_spawn + rowSums(array(init_F[,spawn_seas,] * (fish_sel[1,,spawn_seas,n_ages-1,] * ret_sel[1,,spawn_seas,n_ages-1,] +
                                                                    fish_sel[1,,spawn_seas,n_ages-1,] * (1 - ret_sel[1,,spawn_seas,n_ages-1,]) * dmr[,spawn_seas,]),
                                          dim = c(n_regions, n_fish_fleets)))
      tmp_unfished = spawn_state(tmp_unfished, Movement[1,,, spawn_seas, n_ages - 1], Zu_spawn,
                                 Mrate[1,,, spawn_seas, n_ages - 1], seasdur[spawn_seas], t_spawn, move_timing)
      tmp_fished = spawn_state(tmp_fished, Movement[1,,, spawn_seas, n_ages - 1], Zf_spawn,
                               Mrate[1,,, spawn_seas, n_ages - 1], seasdur[spawn_seas], t_spawn, move_timing)
      SB_age[, n_ages - 1] = tmp_unfished * WAA[1,, spawn_seas, n_ages - 1] * MatAA[1,, spawn_seas, n_ages - 1]
      SB_fished_age[, n_ages - 1] = tmp_fished * WAA[1,, spawn_seas, n_ages - 1] * MatAA[1,, spawn_seas, n_ages - 1]

      ## Plus group analytical solution
      T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

      for (seas in 1:n_seas) {
        # Total mortality by region for each age, built directly rather than recovered
        # from the survival matrices
        Zg_penult_u = natmort[1,, n_ages - 1] * seasdur[seas]
        Zg_plus_u   = natmort[1,, n_ages] * seasdur[seas]
        Zg_penult_f = Zg_penult_u + rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages-1,] * ret_sel[1,,seas,n_ages-1,] +
                                                                      fish_sel[1,,seas,n_ages-1,] * (1 - ret_sel[1,,seas,n_ages-1,]) * dmr[,seas,]),
                                                  dim = c(n_regions, n_fish_fleets)))
        Zg_plus_f = Zg_plus_u + rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages,] * ret_sel[1,,seas,n_ages,] +
                                                                  fish_sel[1,,seas,n_ages,] * (1 - ret_sel[1,,seas,n_ages,]) * dmr[,seas,]),
                                              dim = c(n_regions, n_fish_fleets)))

        # Get transition matrices. build_seas_operator returns row convention, so transpose
        # into the column convention this recursion uses; composed so season 1 applies first.
        opg <- function(age, Z) t(build_seas_operator(Movement[1,,,seas,age], Z, Mrate[1,,,seas,age], seasdur[seas], move_timing))
        T_penult_unfished = opg(n_ages - 1, Zg_penult_u) %*% T_penult_unfished
        T_plus_unfished   = opg(n_ages,     Zg_plus_u) %*% T_plus_unfished
        T_penult_fished   = opg(n_ages - 1, Zg_penult_f) %*% T_penult_fished
        T_plus_fished     = opg(n_ages,     Zg_plus_f) %*% T_plus_fished
      }

      source_unfished = T_penult_unfished %*% Nspr[, n_ages - 1]
      source_fished   = T_penult_fished %*% Nspr_fished[, n_ages - 1]

      Nspr[, n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
      Nspr_fished[, n_ages] = solve(diag(n_regions) - T_plus_fished, source_fished)

      tmp_unfished = Nspr[,n_ages]
      tmp_fished = Nspr_fished[,n_ages]

      if(spawn_seas > 1) {
        for (seas in 1:(spawn_seas - 1)) {

          # Apply seasonal movement and mortality together, per move_timing
          Zu_seas = natmort[1,,n_ages] * seasdur[seas]
          Zf_seas = Zu_seas + rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages,] * ret_sel[1,,seas,n_ages,] +
                                                                fish_sel[1,,seas,n_ages,] * (1 - ret_sel[1,,seas,n_ages,]) * dmr[,seas,]),
                                            dim = c(n_regions, n_fish_fleets)))
          tmp_unfished = advance_seas(tmp_unfished, Movement[1,,,seas,n_ages], Zu_seas,
                                      Mrate[1,,,seas,n_ages], seasdur[seas], move_timing)
          tmp_fished = advance_seas(tmp_fished, Movement[1,,,seas,n_ages], Zf_seas,
                                    Mrate[1,,,seas,n_ages], seasdur[seas], move_timing)
          } # end seas loop
      }

      ## Plus group spawning biomass. The plus group is carried to spawning the same way
      ## as every other age: spawn_state applies both the spawning season movement step
      ## and the t_spawn mortality discount, ordered by move_timing. Applying the discount
      ## alone would fix the ordering at "mortality, then no movement", which is wrong for
      ## any multi-region model at move_timing = 0.
      Zu_spawn_plus = natmort[1,, n_ages] * seasdur[spawn_seas]
      Zf_spawn_plus = Zu_spawn_plus +
        rowSums(array(init_F[,spawn_seas,] * (fish_sel[1,,spawn_seas,n_ages,] * ret_sel[1,,spawn_seas,n_ages,] +
                                                fish_sel[1,,spawn_seas,n_ages,] * (1 - ret_sel[1,,spawn_seas,n_ages,]) * dmr[,spawn_seas,]),
                      dim = c(n_regions, n_fish_fleets)))
      tmp_unfished_spawn = spawn_state(tmp_unfished, Movement[1,,, spawn_seas, n_ages], Zu_spawn_plus,
                                       Mrate[1,,, spawn_seas, n_ages], seasdur[spawn_seas], t_spawn, move_timing)
      tmp_fished_spawn   = spawn_state(tmp_fished, Movement[1,,, spawn_seas, n_ages], Zf_spawn_plus,
                                       Mrate[1,,, spawn_seas, n_ages], seasdur[spawn_seas], t_spawn, move_timing)
      SB_age[, n_ages] = tmp_unfished_spawn * WAA[1,, spawn_seas, n_ages] * MatAA[1,, spawn_seas, n_ages]
      SB_fished_age[, n_ages] = tmp_fished_spawn * WAA[1,, spawn_seas, n_ages] * MatAA[1,, spawn_seas, n_ages]

      # Get global spawning biomass per recruit (scalar)
      S0 = sum(SB_age[,1:n_ages]) * R0[1]
      SF = sum(SB_fished_age[,1:n_ages]) * R0[1]
    }

    # get SSB to use to predict recruitment. When rec_lag == 0 (age-0
    # recruitment), y <= rec_lag is never true for y >= 1, so this always
    # takes the SSB_vals[,,y-rec_lag] = SSB_vals[,,y] branch, the CURRENT
    # year's SSB. The caller (SPoRC_rtmb.R) is responsible for computing
    # SSB_vals[,,y] from survivors only (before this year's recruits exist)
    # and passing it in before calling this function in that case.
    if(y <= rec_lag) SSB = SF else SSB = array(SSB_vals[,,y-rec_lag], dim = c(n_pop, n_regions))

    # Get recruitment based on SSB and R0
    # Single population recruitment
    if(n_pop == 1) {
      for(r in 1:n_regions) {

        # Local Density Dependence (using h[1,r] b/c steepness is region-specific)
        if(rec_dd == 0) {
          local_R0 = R0[1] * rec_region_prop[1,r] # get local R0 based on recruitment proportions
          if(recruitment_model == 1) rec[1,r] = (4 * h[1,r] * local_R0 * SSB[1,r] ) / ( (1 - h[1,r] ) * S0[1,r] + (5 * h[1,r] - 1) * SSB[1,r])
          if(recruitment_model == 2) {
            rel_SSB = SSB[1,r] / S0[1,r] # depletion drives the Ricker, so R0 and S0 stay local
            rec[1,r] = local_R0 * rel_SSB * exp(ricker_alpha(h[1,r]) * (1 - rel_SSB))
          }
        }

        # Global Density Dependence (using h[1,1] b/c steepness is global )
        if(rec_dd == 1) {
          if(recruitment_model == 1) rec[1,r] = (4* h[1,1] * R0[1] * sum(SSB) ) / ((1 - h[1,1] ) * S0 + (5 * h[1,1] - 1) * sum(SSB) ) * rec_region_prop[1,r]
          if(recruitment_model == 2) {
            rel_SSB = sum(SSB) / S0 # S0 is a scalar here b/c density dependence is global
            rec[1,r] = R0[1] * rel_SSB * exp(ricker_alpha(h[1,1]) * (1 - rel_SSB)) * rec_region_prop[1,r]
          }
        }

      } # end r loop
    }

    # Local Density Dependence w/ more than 1 population (using h[p,p] since a given population has the same steepness)
    if(rec_dd == 0 && n_pop > 1) {

      # count number of populations in a given region
      n_pop_in_region = array(0, dim = n_regions)
      for(p in 1:n_pop) n_pop_in_region[natal_region[p]] = n_pop_in_region[natal_region[p]] + 1

      effective_SSB = array(0, dim = n_pop)
      effective_S0  = array(0, dim = n_pop)

      for(p2 in 1:n_pop) {
        for(p in 1:n_pop) {
          if(p == p2) {
            # Own population contribution - no stray scaling
            effective_SSB[p2] = effective_SSB[p2] + SSB[p, natal_region[p2]]
            effective_S0[p2]  = effective_S0[p2]  + S0[p, natal_region[p2]]
          } else {
            n_receivers = n_pop_in_region[natal_region[p2]] # get number of populations in a givenr egion
            # Cross-population contribution scaled by stray_rate
            # SSB[p, natal_region[p2]] already reflects skip spawning via spawning migration
            # stray_rate[p] controls what fraction of those actually contribute here
            effective_SSB[p2] = effective_SSB[p2] + (stray_rate[p] / n_receivers) * SSB[p, natal_region[p2]]
            effective_S0[p2]  = effective_S0[p2]  + (stray_rate[p] / n_receivers) * S0[p, natal_region[p2]]
          }
        }
      }

      for(p in 1:n_pop) {
        if(recruitment_model == 1) rec[p,] = (4 * h[p, natal_region[p]] * R0[p] * effective_SSB[p]) /
          ((1 - h[p, natal_region[p]]) * effective_S0[p] + (5 * h[p, natal_region[p]] - 1) * effective_SSB[p]) * rec_region_prop[p,]
        if(recruitment_model == 2) {
          rel_SSB = effective_SSB[p] / effective_S0[p]
          rec[p,] = R0[p] * rel_SSB * exp(ricker_alpha(h[p, natal_region[p]]) * (1 - rel_SSB)) * rec_region_prop[p,]
        }
      } # end p loop
    }

  } # end Beverton-Holt and Ricker

  # resampling
  if(recruitment_model == 999) rec = NA

  # coerce into array at the end
  rec = array(rec, dim = c(n_pop, n_regions))

  return(rec)
}
