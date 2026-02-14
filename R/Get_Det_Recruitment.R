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
#' @param recruitment_dd Integer flag specifying the scale of density dependence:
#'   \itemize{
#'     \item `0` = local (region-specific)
#'     \item `1` = global (shared across regions)
#'   }
#' @param y Current model year (used for SSB lag indexing)
#' @param rec_lag Recruitment lag (number of years between spawning and recruitment)
#' @param R0 Virgin or mean recruitment (global scalar)
#' @param Rec_Prop Vector of recruitment proportions by region (used to allocate global `R0` under local density dependence)
#' @param h Vector of Beverton–Holt steepness values by region
#' @param n_regions Number of spatial regions
#' @param n_ages Number of modeled age classes
#' @param WAA Matrix of weight-at-age by region and age
#' @param MatAA Matrix of maturity-at-age by region and age
#' @param natmort Matrix or vector of natural mortality by region and age
#' @param SSB_vals Matrix of spawning stock biomass (SSB) by region and year
#' @param Movement 3D array of movement probabilities between regions by age (`[origin, destination, age]` or alternatively, `[n_regions, n_regions, age]`)
#' @param do_recruits_move Logical or integer flag (0/1) indicating whether recruits move during their first year
#' @param t_spawn Fraction of the year at which spawning occurs (used for survival to spawning)
#' @param init_F Scalar for initial F value to apply
#' @param fish_sel Array of fishery selectivity of dominant fleet (fleet 1) dimensioned by n_regions x n_ages
#' @param n_seas Number of seasons
#' @param spawn_seas Season in which spawning happens
#' @param seasdur Fraction of year within a given season
#' @param sexratio_f Vector of sex-ratio values for females (n_regions)
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
                                recruitment_dd,
                                y,
                                rec_lag,
                                R0,
                                Rec_Prop,
                                h,
                                n_regions,
                                n_ages,
                                WAA,
                                MatAA,
                                natmort,
                                SSB_vals,
                                Movement,
                                do_recruits_move,
                                t_spawn,
                                init_F,
                                fish_sel,
                                n_seas,
                                spawn_seas,
                                seasdur,
                                sexratio_f
                                ) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  if(recruitment_model == 0) rec = R0 * Rec_Prop # mean recruitment apportioned across n_regions

  # Beverton-Holt
  if(recruitment_model == 1) {

    # Storage for recruitment, S0, and SF (equilibrium fished)
    rec = S0 = SF = rep(0, n_regions)

    # local density dependence
    if(recruitment_dd == 0) {

      # Calculate unexploited naa per recruit by origin area and destination area
      SB_fished_age = Nspr_fished = SB_age = Nspr = array(0, dim = c(n_regions, n_regions, n_ages))
      SB_fished_mat = SB_unfished_mat = array(0, c(n_regions, n_regions))

      # Set up the initial recruits (1 recruit per area)
      for(o in 1:n_regions) {
        for(d in 1:n_regions) {
          if(o == d) Nspr_fished[o,d,1] = Nspr[o,d,1] = 1 * sexratio_f[o]
          else Nspr_fished[o,d,1] = Nspr[o,d,1] = 0
        } # end d loop
      } # end o loop

      # Loop through ages, projecting each cohort through the full annual cycle
      for(j in 2:(n_ages-1)){

        # Project age j-1 through all seasons to become age j
        for(seas in 1:n_seas) {

          for(o in 1:n_regions) {
            # Get temporary values from origin region
            tmp_unfished = Nspr[o,,j-1]
            tmp_fished = Nspr_fished[o,,j-1]

            # Apply movement
            if(do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
              tmp_unfished = tmp_unfished %*% Movement[,,seas,j-1]
              tmp_fished = tmp_fished %*% Movement[,,seas,j-1]
            }

            # Calculate spawning biomass if this is the spawning season
            if(seas == spawn_seas) {
              for(d in 1:n_regions) {
                SB_age[o,d,j-1] = tmp_unfished[d] * WAA[d,spawn_seas,j-1] * MatAA[d,spawn_seas,j-1] *
                  exp(-(t_spawn * natmort[d,j-1] * seasdur[seas]))
                SB_fished_age[o,d,j-1] = tmp_fished[d] * WAA[d,spawn_seas,j-1] * MatAA[d,spawn_seas,j-1] *
                  exp(-t_spawn * ((natmort[d,j-1] * seasdur[seas]) + init_F[seas] * fish_sel[d,j-1]))
              }
            }

            # Apply mortality
            if(seas < n_seas) {
              # Within-season mortality, no ageing yet
              Nspr[o,,j-1] = tmp_unfished * exp(-(natmort[,j-1] * seasdur[seas]))
              Nspr_fished[o,,j-1] = tmp_fished * exp(-(natmort[,j-1] * seasdur[seas] + init_F[seas] * fish_sel[,j-1]))
            } else {
              # Last season: mortality + ageing
              Nspr[o,,j] = tmp_unfished * exp(-(natmort[,j-1] * seasdur[seas]))
              Nspr_fished[o,,j] = tmp_fished * exp(-(natmort[,j-1] * seasdur[seas] + init_F[seas] * fish_sel[,j-1]))
            }

          } # end o loop
        } # end seas loop
      } # end j loop

      # Now calculate spawning biomass for penultimate age (n_ages-1)
      for(o in 1:n_regions) {
        # Age n_ages-1 is now at start of year after the loop
        tmp_unfished = Nspr[o,,n_ages-1]
        tmp_fished = Nspr_fished[o,,n_ages-1]

        if(spawn_seas > 1) {
          for (seas in 1:(spawn_seas - 1)) {

            # Apply seasonal movement
            tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages-1]
            tmp_fished = tmp_fished %*% Movement[,,seas,n_ages-1]

            # Apply seasonal mortality
            tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages-1] * seasdur[seas]))
            tmp_fished = tmp_fished * exp(-(natmort[,n_ages-1] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages-1]))

          } # end seas loop
        }

        # Apply movement for spawning season
        tmp_unfished = tmp_unfished %*% Movement[,,spawn_seas,n_ages-1]
        tmp_fished = tmp_fished %*% Movement[,,spawn_seas,n_ages-1]

        # Calculate spawning biomass
        for(d in 1:n_regions) {
          SB_age[o,d,n_ages-1] = tmp_unfished[d] * WAA[d,spawn_seas,n_ages-1] * MatAA[d,spawn_seas,n_ages-1] *
            exp(-(t_spawn * natmort[d,n_ages-1] * seasdur[spawn_seas]))
          SB_fished_age[o,d,n_ages-1] = tmp_fished[d] * WAA[d,spawn_seas,n_ages-1] * MatAA[d,spawn_seas,n_ages-1] *
            exp(-t_spawn * ((natmort[d,n_ages-1] * seasdur[spawn_seas]) + init_F[spawn_seas] * fish_sel[d,n_ages-1]))
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
          S_penult_fished = diag(exp(-(natmort[,n_ages-1] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages-1])), n_regions)
          S_plus_fished = diag(exp(-(natmort[,n_ages] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages])), n_regions)
          T_penult_fished = S_penult_fished %*% t(Movement[,,seas,n_ages-1]) %*% T_penult_fished
          T_plus_fished = S_plus_fished %*% t(Movement[,,seas,n_ages]) %*% T_plus_fished
        } # end seas loop

        # Solve for equilibrium plus group (at start of year)
        source_unfished = T_penult_unfished %*% Nspr[o,,n_ages-1]
        Nspr[o,,n_ages] = solve(diag(n_regions) - T_plus_unfished, source_unfished)
        source_fished = T_penult_fished %*% Nspr_fished[o,,n_ages-1]
        Nspr_fished[o,,n_ages] = solve(diag(n_regions) - T_plus_fished, source_fished)

      } # end o loop

        # Now calculate spawning biomass for penultimate age (n_ages-1)
        for(o in 1:n_regions) {
          # Age n_ages-1 is now at start of year after the loop
          tmp_unfished = Nspr[o,,n_ages]
          tmp_fished = Nspr_fished[o,,n_ages]

          if(spawn_seas > 1) {
            for (seas in 1:(spawn_seas - 1)) {

              # Apply seasonal movement
              tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages]
              tmp_fished = tmp_fished %*% Movement[,,seas,n_ages]

              # Apply seasonal mortality
              tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages] * seasdur[seas]))
              tmp_fished = tmp_fished * exp(-(natmort[,n_ages] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages]))

            } # end seas loop
          }

          # Apply movement for spawning season
          tmp_unfished = tmp_unfished %*% Movement[,,spawn_seas,n_ages]
          tmp_fished = tmp_fished %*% Movement[,,spawn_seas,n_ages]

          # Calculate spawning biomass
          for(d in 1:n_regions) {
            SB_age[o,d,n_ages] = tmp_unfished[d] * WAA[d,spawn_seas,n_ages] * MatAA[d,spawn_seas,n_ages] *
              exp(-(t_spawn * natmort[d,n_ages] * seasdur[spawn_seas]))
            SB_fished_age[o,d,n_ages] = tmp_fished[d] * WAA[d,spawn_seas,n_ages] * MatAA[d,spawn_seas,n_ages] *
              exp(-t_spawn * ((natmort[d,n_ages] * seasdur[spawn_seas]) + init_F[spawn_seas] * fish_sel[d,n_ages]))
          } # end d loop

        } # end o loop

      # Remove the old spawning biomass calculation loop entirely
      # parse out and compute unfished spawning biomass per recruit
      for(o in 1:n_regions) {
        for(d in 1:n_regions) {
          # unfished
          SB_unfished_mat[o, d] = sum(SB_age[o, d, ])
          # fished
          SB_fished_mat[o, d] = sum(SB_fished_age[o, d, ])
        } # end d
      } # end o

      for(d in 1:n_regions) {
        S0[d] = sum(SB_unfished_mat[,d] * Rec_Prop * R0) # unfished
        SF[d] = sum(SB_fished_mat[,d] * Rec_Prop * R0) # fished
      } # end d

    } # end if recruitment_dd == 0

    # global density dependence
    if (recruitment_dd == 1) {

      # Setup containers
      SB_fished_age = SB_age = Nspr_fished = Nspr = array(0, dim = c(n_regions, n_ages))

      # Initial recruits: 1 recruit globally, split by Rec_Prop
      Nspr[,1] = Rec_Prop * sexratio_f
      Nspr_fished[,1] = Rec_Prop * sexratio_f

      ## Loop through ages
      for (j in 2:(n_ages - 1)) {
        for (seas in 1:n_seas) {

          tmp_unfished = Nspr[, j - 1]
          tmp_fished   = Nspr_fished[, j - 1]

          ## Movement
          if (do_recruits_move == 1 || (do_recruits_move == 0 && j > 2)) {
            tmp_unfished = as.vector(tmp_unfished %*% Movement[,,seas, j - 1])
            tmp_fished   = as.vector(tmp_fished   %*% Movement[,,seas, j - 1])
          }

          ## Spawning biomass
          if (seas == spawn_seas) {
            SB_age[, j - 1] = tmp_unfished * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] * exp(-t_spawn * natmort[, j - 1] * seasdur[seas])
            SB_fished_age[, j - 1] = tmp_fished * WAA[, spawn_seas, j - 1] * MatAA[, spawn_seas, j - 1] * exp(-t_spawn *
                    (natmort[, j - 1] * seasdur[seas] + init_F[seas] * fish_sel[, j - 1]))
          }

          ## Mortality and ageing
          if (seas < n_seas) { # Within season mortality
            Nspr[, j - 1] = tmp_unfished * exp(-natmort[, j - 1] * seasdur[seas])
            Nspr_fished[, j - 1] =
              tmp_fished * exp(-(natmort[, j - 1] * seasdur[seas] + init_F[seas] * fish_sel[, j - 1]))
          } else {
            # Ageing
            Nspr[, j] = tmp_unfished * exp(-natmort[, j - 1] * seasdur[seas])
            Nspr_fished[, j] = tmp_fished * exp(-(natmort[, j - 1] * seasdur[seas] +  init_F[seas] * fish_sel[, j - 1]))
          }
        }
      }


      # Age n_ages-1 is now at start of year after the loop
      tmp_unfished = Nspr[,n_ages-1]
      tmp_fished = Nspr_fished[,n_ages-1]

      if(spawn_seas > 1) {
        for (seas in 1:(spawn_seas - 1)) {

          # Apply seasonal movement
          tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages-1]
          tmp_fished = tmp_fished %*% Movement[,,seas,n_ages-1]

          # Apply seasonal mortality
          tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages-1] * seasdur[seas]))
          tmp_fished = tmp_fished * exp(-(natmort[,n_ages-1] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages-1]))

        } # end seas loop
      }

      ## Penultimate age spawning biomass
      tmp_unfished = as.vector(tmp_unfished %*% Movement[,, spawn_seas, n_ages - 1])
      tmp_fished   = as.vector(tmp_fished %*% Movement[,, spawn_seas, n_ages - 1])
      SB_age[, n_ages - 1] = tmp_unfished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] *
        exp(-t_spawn * natmort[, n_ages - 1] * seasdur[spawn_seas])
      SB_fished_age[, n_ages - 1] = tmp_fished * WAA[, spawn_seas, n_ages - 1] * MatAA[, spawn_seas, n_ages - 1] * exp(-t_spawn *
              (natmort[, n_ages - 1] * seasdur[spawn_seas] + init_F[spawn_seas] * fish_sel[, n_ages - 1]))

      ## Plus group analytical solution
      T_plus_fished = T_penult_fished = T_plus_unfished = T_penult_unfished = diag(n_regions)

      for (seas in 1:n_seas) {
        # Get survival
        S_penult_unfished = diag(exp(-natmort[, n_ages - 1] * seasdur[seas]), n_regions)
        S_plus_unfished = diag(exp(-natmort[, n_ages] * seasdur[seas]), n_regions)
        S_penult_fished = diag(exp(-(natmort[, n_ages - 1] * seasdur[seas] +  init_F[seas] * fish_sel[, n_ages - 1])), n_regions)
        S_plus_fished = diag(exp(-(natmort[, n_ages] * seasdur[seas] + init_F[seas] * fish_sel[, n_ages])), n_regions)

        # Get transition matrices
        T_penult_unfished = S_penult_unfished %*% t(Movement[,,seas, n_ages - 1]) %*% T_penult_unfished
        T_plus_unfished = S_plus_unfished %*% t(Movement[,,seas, n_ages]) %*% T_plus_unfished
        T_penult_fished = S_penult_fished %*% t(Movement[,,seas, n_ages - 1]) %*% T_penult_fished
        T_plus_fished = S_plus_fished %*% t(Movement[,,seas, n_ages]) %*% T_plus_fished
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
          tmp_unfished = tmp_unfished %*% Movement[,,seas,n_ages]
          tmp_fished = tmp_fished %*% Movement[,,seas,n_ages]

          # Apply seasonal mortality
          tmp_unfished = tmp_unfished * exp(-(natmort[,n_ages] * seasdur[seas]))
          tmp_fished = tmp_fished * exp(-(natmort[,n_ages] * seasdur[seas] + init_F[seas] * fish_sel[,n_ages]))

        } # end seas loop
      }

      ## Plus group spawning biomass
      SB_age[, n_ages] = tmp_unfished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] * exp(-t_spawn * natmort[, n_ages] * seasdur[spawn_seas])

      SB_fished_age[, n_ages] = tmp_fished * WAA[, spawn_seas, n_ages] * MatAA[, spawn_seas, n_ages] * exp(-t_spawn *
              (natmort[, n_ages] * seasdur[spawn_seas] + init_F[spawn_seas] * fish_sel[, n_ages]))

      # Get global spawning biomass per recruit
      S0 = sum(SB_age) * R0
      SF = sum(SB_fished_age) * R0
    }

    # get SSB to use to predict recruitment
    if(y <= rec_lag) SSB = SF else SSB = SSB_vals[,y-rec_lag]

    # Get recruitment based on SSB and R0
    for(r in 1:n_regions) {
      # Local Density Dependence
      if(recruitment_dd == 0) {
        local_R0 = R0 * Rec_Prop[r] # get local R0 based on recruitment proportions
        rec[r] = (4*h[r]*local_R0*SSB[r]) / ((1-h[r])*S0[r] + (5*h[r]-1)*SSB[r]) # get local beverton holt
      }
      # Global Density Dependence
      if(recruitment_dd == 1) {
        rec[r] = (4*h[r]*R0*sum(SSB)) / ((1-h[r])*sum(S0) + (5*h[r]-1)*sum(SSB)) * Rec_Prop[r] # get global beverton holt and then apportion to different regions
      }

    } # end r loop
  } # end Beverton-Holt

  return(rec)
}
