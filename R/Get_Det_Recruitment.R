#' Get Deterministic Recruitment
#'
#' Computes deterministic recruitment for each region based on either
#' a mean recruitment model or a Beverton–Holt stock–recruitment relationship.
#'
#' @param recruitment_model Integer flag specifying the recruitment model:
#'   \itemize{
#'     \item `0` = mean recruitment
#'     \item `1` = Beverton–Holt recruitment with steepness
#'   }
#' @param rec_dd Integer flag specifying the scale of density dependence:
#'   \itemize{
#'     \item `0` = local (region-specific)
#'     \item `1` = global (shared across regions)
#'   }
#' @param y Current model year (used for SSB lag indexing)
#' @param rec_lag Recruitment lag (number of years between spawning and recruitment)
#' @param R0 Virgin or mean recruitment by population (vector)
#' @param Rec_Prop Array of recruitment proportions by population and region (used to allocate population `R0`)
#' @param h Array of Beverton–Holt steepness values by population, region
#' @param n_regions Number of spatial regions
#' @param n_ages Number of modeled age classes
#' @param WAA Array of weight-at-age by population, region and age
#' @param MatAA Array of maturity-at-age by population, region and age
#' @param natmort Array or vector of natural mortality by population, region and age
#' @param SSB_vals Array of spawning stock biomass (SSB) by population, region and year
#' @param Movement Array of movement probabilities between regions by age (`[pop, origin, destination, seas, age]` or alternatively, `[n_pop, n_regions, n_regions, n_seas, age]`)
#' @param do_recruits_move Logical or integer flag (0/1) indicating whether recruits move during their first year
#' @param t_spawn Fraction of the year at which spawning occurs (used for survival to spawning)
#' @param init_F Scalar for initial F value to apply
#' @param fish_sel Array of fishery selectivity of dominant fleet (fleet 1) dimensioned by n_regions x n_ages
#' @param n_seas Number of seasons
#' @param spawn_seas Season in which spawning happens
#' @param seasdur Fraction of year within a given season
#' @param sexratio_f Array of sex-ratio values for females (n_pop x n_regions)
#' @param n_pop Number of populations
#' @param sgl_seas_spawning_movement Array of spawning movement probabilities between regions by age (`[pop, origin, destination, age]` or alternatively, `[n_pop, n_regions, n_regions, age]`)
#' @param natal_region Integer vector of length \code{n_pop}. Maps each population
#'   to its natal region (1-indexed).
#'
#' @details
#' The function returns region-specific deterministic recruitment estimates
#' based on the chosen recruitment model and density dependence structure.
#'
#' When `recruitment_model = 0`, recruitment is fixed at mean values (`R0 * Rec_Prop`).
#' When `recruitment_model = 1`, Beverton–Holt recruitment is applied using:
#' \deqn{R = \frac{4hR_0SSB}{(1 - h)S_0 + (5h - 1)SSB}}
#' where `S_0` is unfished spawning biomass per recruit, computed separately for each
#' region (local) or summed across all regions (global).
#'
#' @return A numeric vector of length `n_regions` containing deterministic
#' recruitment values for each region.
#'
#' @keywords internal
Get_Det_Recruitment <- function(recruitment_model,
                                rec_dd,
                                y,
                                rec_lag,
                                R0,
                                Rec_Prop,
                                h,
                                n_pop,
                                n_regions,
                                n_ages,
                                WAA,
                                MatAA,
                                natmort,
                                SSB_vals,
                                Movement,
                                sgl_seas_spawning_movement,
                                do_recruits_move,
                                t_spawn,
                                init_F,
                                fish_sel,
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
    for(p in 1:n_pop) rec[p,] = R0[p] * Rec_Prop[p,] # mean recruitment apportioned across n_pop and n_regions
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

            if(o == d) Nspr_fished[p,o,d,1] = Nspr[p,o,d,1] = 1 * sexratio_f[p,o] * Rec_Prop[p,o]
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
                    exp(-(t_spawn * natmort[p,d,j-1] * seasdur[seas]))
                  SB_fished_age[p,o,d,j-1] = tmp_fished_spawn[d] * WAA[p,d,spawn_seas,j-1] * MatAA[p,d,spawn_seas,j-1] *
                    exp(-t_spawn * ((natmort[p,d,j-1] * seasdur[seas]) + init_F[seas] * fish_sel[d,j-1]))
                }
              }

              # Apply mortality
              if(seas < n_seas) {
                # Within-season mortality, no ageing yet
                Nspr[p,o,,j-1] = tmp_unfished * exp(-(natmort[p,,j-1] * seasdur[seas]))
                Nspr_fished[p,o,,j-1] = tmp_fished * exp(-(natmort[p,,j-1] * seasdur[seas] + init_F[seas] * fish_sel[,j-1]))
              } else {
                # Last season: mortality + ageing
                Nspr[p,o,,j] = tmp_unfished * exp(-(natmort[p,,j-1] * seasdur[seas]))
                Nspr_fished[p,o,,j] = tmp_fished * exp(-(natmort[p,,j-1] * seasdur[seas] + init_F[seas] * fish_sel[,j-1]))
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
              tmp_fished = tmp_fished * exp(-(natmort[p,,n_ages-1] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages-1]))

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
              exp(-t_spawn * ((natmort[p,d,n_ages-1] * seasdur[spawn_seas]) + init_F[spawn_seas] * fish_sel[d,n_ages-1]))
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
          S_penult_fished = diag(exp(-(natmort[p,,n_ages-1] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages-1])), n_regions)
          S_plus_fished = diag(exp(-(natmort[p,,n_ages] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages])), n_regions)
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

      # Now calculate spawning biomass for penultimate age (n_ages-1)
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
              tmp_fished = tmp_fished * exp(-(natmort[p,,n_ages] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages]))

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
              exp(-t_spawn * ((natmort[p,d,n_ages] * seasdur[spawn_seas]) + init_F[spawn_seas] * fish_sel[d,n_ages]))
          } # end d loop

        } # end o loop
      } # end p loop

      # Remove the old spawning biomass calculation loop entirely
      # parse out and compute unfished spawning biomass per recruit
      for(p in 1:n_pop) {
        for(o in 1:n_regions) {
          for(d in 1:n_regions) {
            # unfished
            SB_unfished_mat[p, o, d] = sum(SB_age[p, o, d, ])
            # fished
            SB_fished_mat[p, o, d] = sum(SB_fished_age[p, o, d, ])
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

      # Initial recruits: 1 recruit globally, split by Rec_Prop
      Nspr[,1] = Rec_Prop[1,] * sexratio_f[1,]
      Nspr_fished[,1] = Rec_Prop[1,] * sexratio_f[1,]

      ## Loop through ages
      for (j in 2:(n_ages - 1)) {
        for (seas in 1:n_seas) {

          tmp_unfished = Nspr[, j - 1]
          tmp_fished   = Nspr_fished[, j - 1]

          ## Movement
          if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
            tmp_unfished = as.vector(tmp_unfished %*% Movement[1,,,seas, j - 1])
            tmp_fished   = as.vector(tmp_fished   %*% Movement[1,,,seas, j - 1])
          }

          ## Spawning biomass
          if (seas == spawn_seas) {
            SB_age[, j - 1] = tmp_unfished * WAA[1,, spawn_seas, j - 1] * MatAA[1,, spawn_seas, j - 1] * exp(-t_spawn * natmort[1,, j - 1] * seasdur[seas])
            SB_fished_age[, j - 1] = tmp_fished * WAA[1,, spawn_seas, j - 1] * MatAA[1,, spawn_seas, j - 1] *
              exp(-t_spawn * (natmort[1,, j - 1] * seasdur[seas] + init_F[seas] * fish_sel[, j - 1]))
          }

          ## Mortality and ageing
          if (seas < n_seas) { # Within season mortality
            Nspr[, j - 1] = tmp_unfished * exp(-natmort[1,, j - 1] * seasdur[seas])
            Nspr_fished[, j - 1] = tmp_fished * exp(-(natmort[1,, j - 1] * seasdur[seas] + init_F[seas] * fish_sel[, j - 1]))
          } else {
            # Ageing
            Nspr[, j] = tmp_unfished * exp(-natmort[1,, j - 1] * seasdur[seas])
            Nspr_fished[, j] = tmp_fished * exp(-(natmort[1,, j - 1] * seasdur[seas] +  init_F[seas] * fish_sel[, j - 1]))
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
          tmp_fished = tmp_fished * exp(-(natmort[1,,n_ages-1] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages-1]))

        } # end seas loop
      } # end if spawn_seas > 1

      ## Penultimate age spawning biomass
      tmp_unfished = as.vector(tmp_unfished %*% Movement[1,,, spawn_seas, n_ages - 1])
      tmp_fished   = as.vector(tmp_fished %*% Movement[1,,, spawn_seas, n_ages - 1])
      SB_age[, n_ages - 1] = tmp_unfished * WAA[1,, spawn_seas, n_ages - 1] * MatAA[1,, spawn_seas, n_ages - 1] *
        exp(-t_spawn * natmort[1,, n_ages - 1] * seasdur[spawn_seas])
      SB_fished_age[, n_ages - 1] = tmp_fished * WAA[1,, spawn_seas, n_ages - 1] * MatAA[1,, spawn_seas, n_ages - 1] *
        exp(-t_spawn * (natmort[1,, n_ages - 1] * seasdur[spawn_seas] + init_F[spawn_seas] * fish_sel[, n_ages - 1]))

      ## Plus group analytical solution
      T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

      for (seas in 1:n_seas) {
        # Get survival
        S_penult_unfished = diag(exp(-natmort[1,, n_ages - 1] * seasdur[seas]), n_regions)
        S_plus_unfished = diag(exp(-natmort[1,, n_ages] * seasdur[seas]), n_regions)
        S_penult_fished = diag(exp(-(natmort[1,, n_ages - 1] * seasdur[seas] +  init_F[seas] * fish_sel[, n_ages - 1])), n_regions)
        S_plus_fished = diag(exp(-(natmort[1,, n_ages] * seasdur[seas] + init_F[seas] * fish_sel[, n_ages])), n_regions)

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
          tmp_fished = tmp_fished * exp(-(natmort[1,,n_ages] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages]))

        } # end seas loop
      }

      ## Plus group spawning biomass
      SB_age[, n_ages] = tmp_unfished * WAA[1,, spawn_seas, n_ages] * MatAA[1,, spawn_seas, n_ages] *
        exp(-t_spawn * natmort[1,, n_ages] * seasdur[spawn_seas])
      SB_fished_age[, n_ages] = tmp_fished * WAA[1,, spawn_seas, n_ages] * MatAA[1,, spawn_seas, n_ages] *
        exp(-t_spawn * (natmort[1,, n_ages] * seasdur[spawn_seas] + init_F[spawn_seas] * fish_sel[, n_ages]))

      # Get global spawning biomass per recruit (scalar)
      S0 = sum(SB_age) * R0[1]
      SF = sum(SB_fished_age) * R0[1]
    }

    # get SSB to use to predict recruitment
    if(y <= rec_lag) SSB = SF else SSB = array(SSB_vals[,,y-rec_lag], dim = c(n_pop, n_regions))

    # Get recruitment based on SSB and R0
    # Single population recruitment
    if(n_pop == 1) {
      for(r in 1:n_regions) {

        # Local Density Dependence (using h[1,r] b/c steepness is region-specific)
        if(rec_dd == 0) {
          local_R0 = R0[1] * Rec_Prop[1,r] # get local R0 based on recruitment proportions
          rec[1,r] = (4 * h[1,r] * local_R0 * SSB[1,r] ) / ( (1 - h[1,r] ) * S0[1,r] + (5 * h[1,r] - 1) * SSB[1,r])
        }

        # Global Density Dependence (using h[1,1] b/c steepness is global )
        if(rec_dd == 1) {
          rec[1,r] = (4* h[1,1] * R0[1] * sum(SSB) ) / ((1 - h[1,1] ) * S0 + (5 * h[1,1] - 1) * sum(SSB) ) * Rec_Prop[1,r]
        }

      } # end r loop
    }

    # Local Density Dependence w/ more than 1 population (using h[p,p] since a given population has the same steepness)
    if(rec_dd == 0 && n_pop > 1) {
      for(p in 1:n_pop) rec[p,] = (4 * h[p,natal_region[p]] *  R0[p] * sum(SSB[p,]) ) /
          ( (1 - h[p,natal_region[p]] ) * sum(S0[p,]) + (5 * h[p,natal_region[p]] -1) * sum(SSB[p,])) * Rec_Prop[p,]
    }

  } # end Beverton-Holt

  # coerce into array at the end
  rec = array(rec, dim = c(n_pop, n_regions))

  return(rec)
}
