#' Deterministic Recruitment
#'
#' Computes deterministic recruitment by population and region using either
#' a mean recruitment model or a Beverton–Holt stock–recruitment relationship.
#'
#' Recruitment is distributed spatially using regional recruitment proportions
#' and seasonal recruitment timing. When Beverton–Holt recruitment is used,
#' unfished spawning biomass per recruit (\eqn{S_0}) is calculated internally
#' by projecting a single recruit through the full seasonal population dynamics,
#' including movement and mortality.
#'
#' @param recruitment_model Integer flag specifying the recruitment model:
#'   \itemize{
#'   \item \code{0} Mean recruitment
#'   \item \code{1} Beverton–Holt recruitment with steepness
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
#'   recruits exist) when \code{rec_lag = 0} -- see \code{SPoRC_rtmb.R},
#'   \code{Simulate_Population.R}, and \code{Do_Population_Projection.R} for
#'   how each population-dynamics loop does this.
#' @param R0 Numeric vector (\code{n_pop}) of unfished recruitment by population.
#' @param rec_region_prop Matrix (\code{n_pop × n_regions}) giving the proportion
#'   of recruitment allocated to each region.
#' @param rec_seas_prop Matrix (\code{n_pop × n_seas}) giving seasonal recruitment
#'   proportions. When \code{rec_lag = 0}, must be zero for every season
#'   before \code{spawn_seas} (age-0 recruits can't predate the spawning
#'   event that produced them) -- validated at setup by
#'   \code{Setup_Mod_Rec}/\code{Setup_Sim_Rec}.
#' @param h Matrix (\code{n_pop × n_regions}) of Beverton–Holt steepness values.
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
#' **Beverton–Holt recruitment**
#'
#' When \code{recruitment_model = 1}, recruitment follows the
#' Beverton–Holt relationship:
#'
#' \deqn{
#' R = \frac{4hR_0SSB}{(1-h)S_0 + (5h-1)SSB}
#' }
#'
#' where:
#' \itemize{
#' \item \eqn{SSB} is spawning biomass lagged by \code{rec_lag} seasons
#'   (or, when \code{rec_lag = 0}, the current year's own spawning biomass --
#'   see the \code{rec_lag} parameter above)
#' \item \eqn{S_0} is unfished spawning biomass per recruit
#' \item \eqn{h} is steepness
#' }
#'
#' \eqn{S_0} (and the age-composition of spawning biomass per recruit more
#' generally) does not depend on \code{rec_lag} -- it is a pure per-recruit,
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
                                sexratio_f
) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(recruitment_model == 0) {
    rec = array(0, dim = c(n_pop, n_regions))
    for(p in 1:n_pop) rec[p,] = R0[p] * rec_region_prop[p,] # mean recruitment apportioned across n_pop and n_regions
  }

  # Beverton-Holt
  if(recruitment_model == 1) {

    # Storage for recruitment, S0, and SF (equilibrium fished)
    rec = S0 = SF = array(0, dim = c(n_pop, n_regions))

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

                # Get spawning biomass per recruit by age
                for(d in 1:n_regions) {
                  SB_age[p,o,d,j-1] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1] *
                    exp(-(t_spawn * natmort[p,d,j-1] * seasdur[spawn_seas]))
                  SB_fished_age[p,o,d,j-1] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1] *
                    exp(-t_spawn * ((natmort[p,d,j-1] * seasdur[spawn_seas]) +
                                      sum(init_F[d,spawn_seas,] * (fish_sel[p,d,spawn_seas,j-1,] * ret_sel[p,d,spawn_seas,j-1,] +
                                                                     fish_sel[p,d,spawn_seas,j-1,] * (1 - ret_sel[p,d,spawn_seas,j-1,]) * dmr[d,spawn_seas,])) ))
                  }
              }

              # Apply mortality
              if(seas < n_seas) {
                # Within-season mortality, no ageing yet
                Nspr[p,o,,j-1] = tmp_unfished * exp(-(natmort[p,,j-1] * seasdur[seas]))
                Nspr_fished[p,o,,j-1] = tmp_fished * exp(-(natmort[p,,j-1] * seasdur[seas] +
                                                             rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,j-1,] * ret_sel[p,,seas,j-1,] +
                                                                      fish_sel[p,,seas,j-1,] * (1 - ret_sel[p,,seas,j-1,]) * dmr[,seas,]),
                                                                           dim = c(n_regions, n_fish_fleets))) ))
                } else {
                # Last season: mortality + ageing
                Nspr[p,o,,j] = tmp_unfished * exp(-(natmort[p,,j-1] * seasdur[seas]))
                Nspr_fished[p,o,,j] = tmp_fished * exp(-(natmort[p,,j-1] * seasdur[seas] +
                                                           rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,j-1,] * ret_sel[p,,seas,j-1,] +
                                                                                             fish_sel[p,,seas,j-1,] * (1 - ret_sel[p,,seas,j-1,]) * dmr[,seas,]),
                                                                         dim = c(n_regions, n_fish_fleets))) ))
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

              # Apply seasonal movement
              tmp_unfished = tmp_unfished %*% Movement[p,,,seas,n_ages-1]
              tmp_fished = tmp_fished %*% Movement[p,,,seas,n_ages-1]

              # Apply seasonal mortality
              tmp_unfished = tmp_unfished * exp(-(natmort[p,,n_ages-1] * seasdur[seas]))
              tmp_fished = tmp_fished * exp(-(natmort[p,,n_ages-1] * seasdur[seas] +
                                                rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages-1,] * ret_sel[p,,seas,n_ages-1,] +
                                                              fish_sel[p,,seas,n_ages-1,] * (1 - ret_sel[p,,seas,n_ages-1,]) * dmr[,seas,]),
                                                              dim = c(n_regions, n_fish_fleets))) ))
            } # end seas loop
          }

          # Extract temporary variables out
          tmp_unfished_spawn = tmp_unfished
          tmp_fished_spawn = tmp_fished

          # Apply movement for spawning season
          tmp_unfished_spawn = tmp_unfished_spawn %*% Movement[p,,,spawn_seas,n_ages-1]
          tmp_fished_spawn = tmp_fished_spawn %*% Movement[p,,,spawn_seas,n_ages-1]

          # If single season natal homing population
          if(n_pop > 1 && n_seas == 1) {
            # Get NAA during spawning in single season case
            tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages-1]
            tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages-1]
          }

          # Calculate spawning biomass
          for(d in 1:n_regions) {
            SB_age[p,o,d,n_ages-1] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,n_ages-1] * MatAA[p,d,spawn_seas,n_ages-1] *
              exp(-(t_spawn * natmort[p,d,n_ages-1] * seasdur[spawn_seas]))
            SB_fished_age[p,o,d,n_ages-1] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,n_ages-1] * MatAA[p,d,spawn_seas,n_ages-1] *
              exp(-t_spawn * ((natmort[p,d,n_ages-1] * seasdur[spawn_seas]) +
                                sum(init_F[d,spawn_seas,] * (fish_sel[p,d,spawn_seas,n_ages-1,] * ret_sel[p,d,spawn_seas,n_ages-1,] +
                                                               fish_sel[p,d,spawn_seas,n_ages-1,] * (1 - ret_sel[p,d,spawn_seas,n_ages-1,]) * dmr[d,spawn_seas,])) ))
            }
        }
      }

      # Set up analytical solution for plus group
      for(p in 1:n_pop) {

        # Build FULL annual transition for penultimate and plus ages
        T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

        # Loop through ALL seasons to build annual transition matrix
        for(seas in 1:n_seas) {
          # Unfished
          S_penult_unfished = diag(exp(-(natmort[p,,n_ages-1] * seasdur[seas])), n_regions)
          S_plus_unfished = diag(exp(-(natmort[p,,n_ages] * seasdur[seas])), n_regions)
          T_penult_unfished = S_penult_unfished %*% t(Movement[p,,,seas,n_ages-1]) %*% T_penult_unfished
          T_plus_unfished = S_plus_unfished %*% t(Movement[p,,,seas,n_ages]) %*% T_plus_unfished

          # Fished
          F_penult = rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages-1,] * ret_sel[p,,seas,n_ages-1,] +
                             fish_sel[p,,seas,n_ages-1,] * (1 - ret_sel[p,,seas,n_ages-1,]) * dmr[,seas,]),
                                   dim = c(n_regions, n_fish_fleets)))
          F_plus = rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages,] * ret_sel[p,,seas,n_ages,] +
                           fish_sel[p,,seas,n_ages,] * (1 - ret_sel[p,,seas,n_ages,]) * dmr[,seas,]),
                                 dim = c(n_regions, n_fish_fleets)))
          S_penult_fished = diag(exp(-(natmort[p,,n_ages-1] * seasdur[seas] + F_penult)), n_regions)
          S_plus_fished = diag(exp(-(natmort[p,,n_ages] * seasdur[seas] + F_plus)), n_regions)
          T_penult_fished = S_penult_fished %*% t(Movement[p,,,seas,n_ages-1]) %*% T_penult_fished
          T_plus_fished = S_plus_fished %*% t(Movement[p,,,seas,n_ages]) %*% T_plus_fished
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

              # Apply seasonal movement
              tmp_unfished = tmp_unfished %*% Movement[p,,,seas,n_ages]
              tmp_fished = tmp_fished %*% Movement[p,,,seas,n_ages]

              # Apply seasonal mortality
              tmp_unfished = tmp_unfished * exp(-(natmort[p,,n_ages] * seasdur[seas]))
              tmp_fished = tmp_fished * exp(-(natmort[p,,n_ages] * seasdur[seas] +
                                                rowSums(array(init_F[,seas,] * (fish_sel[p,,seas,n_ages,] * ret_sel[p,,seas,n_ages,] +
                                                        fish_sel[p,,seas,n_ages,] * (1 - ret_sel[p,,seas,n_ages,]) * dmr[,seas,]),
                                                              dim = c(n_regions, n_fish_fleets)))))

            } # end seas loop
          }

          # Extract temporary variables out
          tmp_unfished_spawn = tmp_unfished
          tmp_fished_spawn = tmp_fished

          # Apply movement for spawning season
          tmp_unfished_spawn = tmp_unfished_spawn %*% Movement[p,,,spawn_seas,n_ages]
          tmp_fished_spawn = tmp_fished_spawn %*% Movement[p,,,spawn_seas,n_ages]

          # If single season natal homing population
          if(n_pop > 1 && n_seas == 1) {
            # Get NAA during spawning in single season case
            tmp_unfished_spawn = tmp_unfished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages]
            tmp_fished_spawn = tmp_fished_spawn %*% sgl_seas_spawning_movement[p,,,n_ages]
          }

          # Calculate spawning biomass
          for(d in 1:n_regions) {
            SB_age[p,o,d,n_ages] = tmp_unfished_spawn[d] * WAA[p,d,spawn_seas,n_ages] * MatAA[p,d,spawn_seas,n_ages] *
              exp(-(t_spawn * natmort[p,d,n_ages] * seasdur[spawn_seas]))
            SB_fished_age[p,o,d,n_ages] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,n_ages] * MatAA[p,d,spawn_seas,n_ages] *
              exp(-t_spawn * ((natmort[p,d,n_ages] * seasdur[spawn_seas]) +
                                sum(init_F[d,spawn_seas,] * (fish_sel[p,d,spawn_seas,n_ages,] * ret_sel[p,d,spawn_seas,n_ages,] +
                                                               fish_sel[p,d,spawn_seas,n_ages,] * (1 - ret_sel[p,d,spawn_seas,n_ages,]) * dmr[d,spawn_seas,]))))
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

          ## Movement
          if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
            tmp_unfished = as.vector(tmp_unfished %*% Movement[1,,,seas, j - 1])
            tmp_fished   = as.vector(tmp_fished   %*% Movement[1,,,seas, j - 1])
          }

          ## Spawning biomass
          if (seas == spawn_seas) {
            SB_age[, j - 1] = tmp_unfished * WAA[1,, spawn_seas, j - 1] * MatAA[1,, spawn_seas, j - 1] * exp(-t_spawn * natmort[1,, j - 1] * seasdur[seas])
            SB_fished_age[, j - 1] = tmp_fished * WAA[1,, spawn_seas, j - 1] * MatAA[1,, spawn_seas, j - 1] *
              exp(-t_spawn * (natmort[1,, j - 1] * seasdur[seas] +
                                rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,j-1,] * ret_sel[1,,seas,j-1,] +
                                                                  fish_sel[1,,seas,j-1,] * (1 - ret_sel[1,,seas,j-1,]) * dmr[,seas,]),
                                              dim = c(n_regions, n_fish_fleets))) ))
            }

          ## Mortality and ageing
          if (seas < n_seas) { # Within season mortality
            Nspr[, j - 1] = tmp_unfished * exp(-natmort[1,, j - 1] * seasdur[seas])
            Nspr_fished[, j - 1] = tmp_fished * exp(-(natmort[1,, j - 1] * seasdur[seas] +
                                                        rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,j-1,] * ret_sel[1,,seas,j-1,] +
                                                                      fish_sel[1,,seas,j-1,] * (1 - ret_sel[1,,seas,j-1,]) * dmr[,seas,]),
                                                                      dim = c(n_regions, n_fish_fleets))) ))
            } else {
            # Ageing
            Nspr[, j] = tmp_unfished * exp(-natmort[1,, j - 1] * seasdur[seas])
            Nspr_fished[, j] = tmp_fished * exp(-(natmort[1,, j - 1] * seasdur[seas] +
                                                    rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,j-1,] * ret_sel[1,,seas,j-1,] +
                                                                                      fish_sel[1,,seas,j-1,] * (1 - ret_sel[1,,seas,j-1,]) * dmr[,seas,]),
                                                                  dim = c(n_regions, n_fish_fleets))) ))
          }
        } # end seas loop
      } # end j loop

      # Age n_ages-1 is now at start of year after the loop
      tmp_unfished = Nspr[,n_ages-1]
      tmp_fished = Nspr_fished[,n_ages-1]

      if(spawn_seas > 1) {
        for (seas in 1:(spawn_seas - 1)) {

          # Apply seasonal movement
          tmp_unfished = tmp_unfished %*% Movement[1,,,seas,n_ages-1]
          tmp_fished = tmp_fished %*% Movement[1,,,seas,n_ages-1]

          # Apply seasonal mortality
          tmp_unfished = tmp_unfished * exp(-(natmort[1,,n_ages-1] * seasdur[seas]))
          tmp_fished = tmp_fished * exp(-(natmort[1,,n_ages-1] * seasdur[seas] +
                                            rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages-1,] * ret_sel[1,,seas,n_ages-1,] +
                                                                              fish_sel[1,,seas,n_ages-1,] * (1 - ret_sel[1,,seas,n_ages-1,]) * dmr[,seas,]),
                                                          dim = c(n_regions, n_fish_fleets)))))

        } # end seas loop
      } # end if spawn_seas > 1

      ## Penultimate age spawning biomass
      tmp_unfished = as.vector(tmp_unfished %*% Movement[1,,, spawn_seas, n_ages - 1])
      tmp_fished   = as.vector(tmp_fished %*% Movement[1,,, spawn_seas, n_ages - 1])
      SB_age[, n_ages - 1] = tmp_unfished * WAA[1,, spawn_seas, n_ages - 1] * MatAA[1,, spawn_seas, n_ages - 1] *
        exp(-t_spawn * natmort[1,, n_ages - 1] * seasdur[spawn_seas])
      SB_fished_age[, n_ages - 1] = tmp_fished * WAA[1,, spawn_seas, n_ages - 1] * MatAA[1,, spawn_seas, n_ages - 1] *
        exp(-t_spawn * (natmort[1,, n_ages - 1] * seasdur[spawn_seas] +
                          rowSums(array(init_F[,spawn_seas,] * (fish_sel[1,,spawn_seas,n_ages-1,] * ret_sel[1,,spawn_seas,n_ages-1,] +
                                                                  fish_sel[1,,spawn_seas,n_ages-1,] * (1 - ret_sel[1,,spawn_seas,n_ages-1,]) * dmr[,spawn_seas,]),
                                        dim = c(n_regions, n_fish_fleets)))))

      ## Plus group analytical solution
      T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

      for (seas in 1:n_seas) {
        # Get survival
        S_penult_unfished = diag(exp(-natmort[1,, n_ages - 1] * seasdur[seas]), n_regions)
        S_plus_unfished = diag(exp(-natmort[1,, n_ages] * seasdur[seas]), n_regions)
        S_penult_fished = diag(exp(-(natmort[1,, n_ages - 1] * seasdur[seas] +
                                       rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages-1,] * ret_sel[1,,seas,n_ages-1,] +
                                                                         fish_sel[1,,seas,n_ages-1,] * (1 - ret_sel[1,,seas,n_ages-1,]) * dmr[,seas,]),
                                                     dim = c(n_regions, n_fish_fleets))))), n_regions)
        S_plus_fished = diag(exp(-(natmort[1,, n_ages] * seasdur[seas] +
                                     rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages,] * ret_sel[1,,seas,n_ages,] +
                                                                       fish_sel[1,,seas,n_ages,] * (1 - ret_sel[1,,seas,n_ages,]) * dmr[,seas,]),
                                                   dim = c(n_regions, n_fish_fleets))))), n_regions)
        # Get transition matrices
        T_penult_unfished = S_penult_unfished %*% t(Movement[1,,,seas, n_ages - 1]) %*% T_penult_unfished
        T_plus_unfished = S_plus_unfished %*% t(Movement[1,,,seas, n_ages]) %*% T_plus_unfished
        T_penult_fished = S_penult_fished %*% t(Movement[1,,,seas, n_ages - 1]) %*% T_penult_fished
        T_plus_fished = S_plus_fished %*% t(Movement[1,,,seas, n_ages]) %*% T_plus_fished
      }

      source_unfished = T_penult_unfished %*% Nspr[, n_ages - 1]
      source_fished   = T_penult_fished %*% Nspr_fished[, n_ages - 1]

      Nspr[, n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
      Nspr_fished[, n_ages] = solve(diag(n_regions) - T_plus_fished, source_fished)

      tmp_unfished = Nspr[,n_ages]
      tmp_fished = Nspr_fished[,n_ages]

      if(spawn_seas > 1) {
        for (seas in 1:(spawn_seas - 1)) {

          # Apply seasonal movement
          tmp_unfished = tmp_unfished %*% Movement[1,,,seas,n_ages]
          tmp_fished = tmp_fished %*% Movement[1,,,seas,n_ages]

          # Apply seasonal mortality
          tmp_unfished = tmp_unfished * exp(-(natmort[1,,n_ages] * seasdur[seas]))
          tmp_fished = tmp_fished * exp(-(natmort[1,,n_ages] * seasdur[seas] +
                                            rowSums(array(init_F[,seas,] * (fish_sel[1,,seas,n_ages,] * ret_sel[1,,seas,n_ages,] +
                                                                              fish_sel[1,,seas,n_ages,] * (1 - ret_sel[1,,seas,n_ages,]) * dmr[,seas,]),
                                                          dim = c(n_regions, n_fish_fleets)))))
          } # end seas loop
      }

      ## Plus group spawning biomass
      SB_age[, n_ages] = tmp_unfished * WAA[1,, spawn_seas, n_ages] * MatAA[1,, spawn_seas, n_ages] *
        exp(-t_spawn * natmort[1,, n_ages] * seasdur[spawn_seas])
      SB_fished_age[, n_ages] = tmp_fished * WAA[1,, spawn_seas, n_ages] * MatAA[1,, spawn_seas, n_ages] *
        exp(-t_spawn * (natmort[1,, n_ages] * seasdur[spawn_seas] +
                          rowSums(array(init_F[,spawn_seas,] * (fish_sel[1,,spawn_seas,n_ages,] * ret_sel[1,,spawn_seas,n_ages,] +
                                                                  fish_sel[1,,spawn_seas,n_ages,] * (1 - ret_sel[1,,spawn_seas,n_ages,]) * dmr[,spawn_seas,]),
                                        dim = c(n_regions, n_fish_fleets)))))

      # Get global spawning biomass per recruit (scalar)
      S0 = sum(SB_age[,1:n_ages]) * R0[1]
      SF = sum(SB_fished_age[,1:n_ages]) * R0[1]
    }

    # get SSB to use to predict recruitment. When rec_lag == 0 (age-0
    # recruitment), y <= rec_lag is never true for y >= 1, so this always
    # takes the SSB_vals[,,y-rec_lag] = SSB_vals[,,y] branch -- the CURRENT
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
          rec[1,r] = (4 * h[1,r] * local_R0 * SSB[1,r] ) / ( (1 - h[1,r] ) * S0[1,r] + (5 * h[1,r] - 1) * SSB[1,r])
        }

        # Global Density Dependence (using h[1,1] b/c steepness is global )
        if(rec_dd == 1) {
          rec[1,r] = (4* h[1,1] * R0[1] * sum(SSB) ) / ((1 - h[1,1] ) * S0 + (5 * h[1,1] - 1) * sum(SSB) ) * rec_region_prop[1,r]
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

      for(p in 1:n_pop) rec[p,] = (4 * h[p, natal_region[p]] * R0[p] * effective_SSB[p]) /
        ((1 - h[p, natal_region[p]]) * effective_S0[p] + (5 * h[p, natal_region[p]] - 1) * effective_SSB[p]) * rec_region_prop[p,]
    }

  } # end Beverton-Holt

  # resampling
  if(recruitment_model == 999) rec = NA

  # coerce into array at the end
  rec = array(rec, dim = c(n_pop, n_regions))

  return(rec)
}

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
                          sgl_seas_spawning_movement, natal_region, stray_rate) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Get NAA for spawning
  tmp_NAA_spawn = NAA[,,y,seas,,, drop = FALSE]
  tmp_NAA0_spawn = NAA0[,,y,seas,,, drop = FALSE]

  # If we we are natal homing with 1 season
  if(n_seas == 1 && n_pop > 1) {
    # Get NAA during spawning
    for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
      tmp_NAA_spawn[p,,1,1,a,s] = tmp_NAA_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
      tmp_NAA0_spawn[p,,1,1,a,s] = tmp_NAA0_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
    } # end s loop
  }

  # Total Biomass
  Total_Biom_y = apply(tmp_NAA_spawn *
                         WAA[,, y, seas, , ,drop = FALSE] *
                         exp(-ZAA[,,y,seas,,,drop = FALSE] * t_spawn), c(1,2), sum)

  # Spawning Stock Biomass
  SSB_y = apply(tmp_NAA_spawn[,, 1, 1, , 1,drop = FALSE] *
                  WAA[,, y, seas, , 1,drop = FALSE] *
                  MatAA[,, y, seas, , 1,drop = FALSE] *
                  exp(-ZAA[,, y, seas, , 1,drop = FALSE] * t_spawn), c(1,2), sum)

  # Get dynamic B0
  SSB0_array = tmp_NAA0_spawn[,, 1, 1, , 1,drop = FALSE] *  WAA[,,  y, seas, , 1, drop = FALSE] * MatAA[,,y, seas, , 1, drop = FALSE]
  mort_spawn = exp(-natmort[,, y, , 1, drop = FALSE] * t_spawn * seasdur[seas])
  mort_spawn = array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
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
# Computes SSB / Dynamic_SSB0 / eff_SSB for projection year y at season seas
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
                            sgl_seas_spawning_movement, natal_region, stray_rate) {

  tmp_NAA_spawn = proj_NAA[,,y,seas,,, drop = FALSE]
  tmp_NAA0_spawn = proj_NAA0[,,y,seas,,, drop = FALSE]

  # If we we are natal homing with 1 season
  if(n_seas == 1 && n_pop > 1) {
    for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
      tmp_NAA_spawn[p,,1,1,a,s] = tmp_NAA_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
      tmp_NAA0_spawn[p,,1,1,a,s] = tmp_NAA0_spawn[p,,1,1,a,s] %*% sgl_seas_spawning_movement[p,,,y,a,s]
    } # end s loop
  }

  # get SSB
  SSB_y = apply(tmp_NAA_spawn[,, 1, 1, , 1,drop = FALSE] *
                  WAA[,, y, seas, , 1,drop = FALSE] *
                  MatAA[,, y, seas, , 1,drop = FALSE] *
                  exp(-proj_ZAA[,, y, seas, , 1,drop = FALSE] * t_spawn), c(1,2), sum)

  # Get dynamic B0
  SSB0_array = tmp_NAA0_spawn[,, 1, 1, , 1,drop = FALSE] *  WAA[,,  y, seas, , 1, drop = FALSE] * MatAA[,,y, seas, , 1, drop = FALSE]
  mort_spawn = exp(-natmort[,, y, , 1, drop = FALSE] * t_spawn * seasdur[seas])
  mort_spawn = array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
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

  list(SSB_y = SSB_y, Dynamic_SSB0_y = Dynamic_SSB0_y, eff_SSB_y = eff_SSB_y)
}
