# Stage 2 of 3: objective function
#
# The forward projection at the centre of the objective function. Walks the
# population through every year and season, applying recruitment, mortality and
# movement in the order set by move_timing, and records the state each
# observation model needs.

#' Population projection (numbers-at-age dynamics)
#'
#' Advances numbers-at-age forward through all modeled years and seasons:
#' inserts recruitment (timing controlled by \code{rec_lag}), applies
#' movement, computes SSB/biomass quantities via \code{compute_biom_y}, and
#' applies mortality/ageing. Called once from the "Population Projection"
#' section of \code{SPoRC_rtmb.R}. \code{ZAA} (total mortality at age) must
#' already be computed before calling this, since it is treated as an input
#' here rather than derived from \code{NAA}.
#'
#' All array arguments matching an output name (\code{NAA}, \code{NAA0},
#' \code{NAA_bef}, \code{NAA_aft}, \code{Rec}, \code{SSB}, \code{Total_Biom},
#' \code{Dynamic_SSB0}, \code{eff_SSB}) are passed in already dimensioned
#' (typically all-zero, aside from any initial-year values already inserted
#' upstream) and returned fully populated over \code{1:n_yrs}.
#'
#' @param n_pop,n_regions,n_seas,n_ages,n_sexes,n_yrs,n_fish_fleets,n_est_rec_devs
#'   Dimension sizes.
#' @param rec_lag Integer. Recruitment timing: \code{0} inserts recruitment
#'   within the spawning-season biomass computation; non-zero inserts
#'   recruitment once per year ahead of the seasonal loop.
#' @param rec_model,rec_dd,R0,rec_region_prop,rec_seas_prop,h_trans,natal_region,t_spawn,spawn_seas,seasdur,init_F
#'   Recruitment and timing arguments passed through to
#'   \code{Get_Det_Recruitment}.
#' @param n_est_rec_devs Number of estimated recruitment deviations.
#' @param ln_RecDevs Array \code{[pop, region, year]} of log recruitment
#'   deviations; applied multiplicatively to deterministic recruitment for
#'   \code{y <= n_est_rec_devs}.
#' @param sexratio Array \code{[pop, region, year, sex]} of recruitment sex
#'   ratio.
#' @param WAA,MatAA Arrays \code{[pop, region, year, season, age, sex]} of
#'   weight-at-age and maturity-at-age.
#' @param natmort Array \code{[pop, region, year, age, sex]} of natural
#'   mortality at age.
#' @param Movement Array \code{[pop, region_from, region_to, year, season,
#'   age, sex]} of movement rates.
#' @param stray_rate Array \code{[pop, year]} of stray rate.
#' @param sgl_seas_spawning_movement Array \code{[pop, region_from,
#'   region_to, year, age, sex]} of single-season-spawning movement rates.
#' @param do_recruits_move Integer (0/1) switch for whether age-1 recruits
#'   are subject to movement.
#' @param fish_sel,ret_sel Arrays \code{[pop, region, year, season, age, sex,
#'   fish_fleet]} of total/retained fishery selectivity.
#' @param dmr Array \code{[region, year, season, fish_fleet]} of discard
#'   mortality rate.
#' @param ZAA Array \code{[pop, region, year, season, age, sex]} of total
#'   mortality at age (precomputed).
#' @param NAA,NAA0 Arrays \code{[pop, region, year+1, season, age, sex]},
#'   output containers for fished/unfished numbers at age.
#' @param NAA_bef,NAA_aft Arrays \code{[pop, region, year+1, season, age,
#'   sex]}, output containers for numbers at age immediately before/after
#'   movement.
#' @param Rec Array \code{[pop, region, year]}, output container for total
#'   recruitment before seasonal apportionment.
#' @param SSB,Total_Biom,Dynamic_SSB0 Arrays \code{[pop, region, year]},
#'   output containers.
#' @param eff_SSB Array \code{[pop, year]}, output container for effective
#'   (natal-homing-adjusted) SSB.
#'
#' @return List with elements \code{NAA}, \code{NAA0}, \code{NAA_bef},
#'   \code{NAA_aft}, \code{Rec}, \code{SSB}, \code{Total_Biom},
#'   \code{Dynamic_SSB0}, \code{eff_SSB}, \code{Aggregated_SSB} (array
#'   \code{[year]}, SSB summed across pop/region),
#'   \code{Dynamic_Aggregated_SSB0} (array \code{[year]}, likewise for
#'   \code{Dynamic_SSB0}), and \code{NAA_int} (array \code{[pop, region, year,
#'   season, age, sex]}). \code{NAA_int} holds the season-integrated abundance
#'   needed by the spatial Baranov catch equation and is populated only when
#'   \code{move_timing = 2}; it is all zeros otherwise.
#'
#' @keywords internal
#' @import RTMB
get_population_projection <- function(n_pop, n_regions, n_seas, n_ages, n_sexes, n_yrs, n_fish_fleets, n_est_rec_devs,
                                       rec_lag, rec_model, rec_dd, R0, rec_region_prop, rec_seas_prop, h_trans,
                                       natal_region, t_spawn, spawn_seas, seasdur, init_F, ln_RecDevs,
                                       sexratio, WAA, MatAA, natmort, Movement, stray_rate, sgl_seas_spawning_movement,
                                       do_recruits_move, fish_sel, ret_sel, dmr, ZAA,
                                       NAA, NAA0, NAA_bef, NAA_aft, Rec, SSB, Total_Biom, Dynamic_SSB0, eff_SSB,
                                       Mrate = NULL, move_timing = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Season-integrated abundance for the spatial Baranov, filled in below only under move_timing == 2.
  NAA_int <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))

  for(y in 1:n_yrs) {

    ### Annual Recruitment (rec_lag != 0 only) -----------------------------------
    if(rec_lag != 0) {

      # Get Deterministic Recruitment
      tmp_Det_Rec <- Get_Det_Recruitment(recruitment_model = rec_model,
                                        rec_dd = rec_dd,
                                        R0 = R0,
                                        rec_region_prop = rec_region_prop,
                                        rec_seas_prop = rec_seas_prop,
                                        h = h_trans,
                                        n_pop = n_pop,
                                        n_ages = n_ages,
                                        n_regions = n_regions,
                                        # Note: Using first year and female quantities to compute unfished SSB0
                                        sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,1,1], dim = c(n_pop, n_regions)),
                                        WAA = array(WAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                        MatAA = array(MatAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                        natmort = array(natmort[,,1,,1], dim = c(n_pop, n_regions, n_ages)),
                                        Movement = array(Movement[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                        stray_rate = array(stray_rate[,1], dim = c(n_pop)),
                                        sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,1,,1], dim = c(n_pop, n_regions, n_regions, n_ages)),
                                        do_recruits_move = do_recruits_move,
                                        natal_region = natal_region,
                                        t_spawn = t_spawn,
                                        SSB_vals = SSB,
                                        y = y,
                                        n_seas = n_seas,
                                        spawn_seas = spawn_seas,
                                        seasdur = seasdur,
                                        rec_lag = rec_lag,
                                        n_fish_fleets = n_fish_fleets,
                                        init_F = init_F, # initF
                                        fish_sel = array(fish_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # total fishery selectivity
                                        ret_sel = array(ret_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # retained fishery selectivity in first year
                                        dmr = array(dmr[,1,,], dim = c(n_regions, n_seas, n_fish_fleets)),
                                        Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                        move_timing = move_timing
      )

      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            if(y <= n_est_rec_devs) tmp_total_rec <- tmp_Det_Rec[p,r] * exp(ln_RecDevs[p,r,y])
            if(y > n_est_rec_devs) tmp_total_rec <- tmp_Det_Rec[p,r]
            # season 1 fraction
            NAA[p,r,y,1,1,s] <- tmp_total_rec * rec_seas_prop[p,1] * sexratio[p,r,y,s]
          }
          Rec[p,r,y] <- tmp_total_rec  # store total before seasonal split
          NAA0[p,r,y,1,1,] <- NAA[p,r,y,1,1,]
        } # end r loop
      } # end p loop
    } # end if rec_lag != 0

    for(seas in 1:n_seas) {

      # Insert seasonal recruits
      if(if(rec_lag != 0) seas > 1 else seas > spawn_seas) {
        for(p in 1:n_pop) {
          for(r in 1:n_regions) {
            for(s in 1:n_sexes) {
              NAA[p,r,y,seas,1,s]  <- NAA[p,r,y,seas,1,s]  + Rec[p,r,y] * rec_seas_prop[p,seas] * sexratio[p,r,y,s]
              NAA0[p,r,y,seas,1,s] <- NAA0[p,r,y,seas,1,s] + Rec[p,r,y] * rec_seas_prop[p,seas] * sexratio[p,r,y,s]
            } # end s loop
          } # end r loop
        } # end p loop
      }

      ### Movement ----------------------------------------------------------------
      # Record values prior to movement
      NAA_bef[,,y,seas,,] <- NAA[,,y,seas,,]

      # Movement is applied at the start of the season only under move_timing == 0.
      # Under move_timing 1 (mortality then movement) and 2 (continuous) it is folded
      # into the mortality/ageing step at the end of the season instead.
      if(n_regions > 1 && move_timing == 0) {
        for(p in 1:n_pop) {
          # Recruits don't move
          if(do_recruits_move == 0) {
            # Apply movement after ageing processes - start movement at age 2
            for(a in 2:n_ages) {
              for(s in 1:n_sexes) {
                NAA[p,,y,seas,a,s] <- t(NAA[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # Fished
              } # end s loop
            } # end a loop
          } # end if recruits don't move

          # Recruits move here
          if(do_recruits_move == 1) {
            for(a in 1:n_ages) {
              for(s in 1:n_sexes) {
                NAA[p,,y,seas,a,s] <- t(NAA[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # Fished
              } # end s loop
            } # end a loop
          } # end if
        } # end p loop

        # Record values after movement
        NAA_aft[,,y,seas,,] <- NAA[,,y,seas,,]

      } # only compute if spatial

      ### Compute Biomass Quantities + Recruitment (rec_lag == 0 only) ------------
      if(rec_lag == 0 && seas == spawn_seas) {

        # SSB from survivors only
        biom <- compute_biom_y(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, seasdur,
                              n_seas, n_pop, n_regions, n_ages, n_sexes,
                              sgl_seas_spawning_movement, natal_region, stray_rate,
                              Movement, Mrate, move_timing, do_recruits_move)
        SSB[,, y] <- biom$SSB_y

        tmp_Det_Rec <- Get_Det_Recruitment(recruitment_model = rec_model,
                                          rec_dd = rec_dd,
                                          R0 = R0,
                                          rec_region_prop = rec_region_prop,
                                          rec_seas_prop = rec_seas_prop,
                                          h = h_trans,
                                          n_pop = n_pop,
                                          n_ages = n_ages,
                                          n_regions = n_regions,
                                          sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,1,1], dim = c(n_pop, n_regions)),
                                          WAA = array(WAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                          MatAA = array(MatAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                          natmort = array(natmort[,,1,,1], dim = c(n_pop, n_regions, n_ages)),
                                          Movement = array(Movement[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                          stray_rate = array(stray_rate[,1], dim = c(n_pop)),
                                          sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,1,,1], dim = c(n_pop, n_regions, n_regions, n_ages)),
                                          do_recruits_move = do_recruits_move,
                                          natal_region = natal_region,
                                          t_spawn = t_spawn,
                                          SSB_vals = SSB,
                                          y = y,
                                          n_seas = n_seas,
                                          spawn_seas = spawn_seas,
                                          seasdur = seasdur,
                                          rec_lag = rec_lag,
                                          n_fish_fleets = n_fish_fleets,
                                          init_F = init_F, # initF
                                          fish_sel = array(fish_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # total fishery selectivity
                                          ret_sel = array(ret_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # retained fishery selectivity in first year
                                          dmr = array(dmr[,1,,], dim = c(n_regions, n_seas, n_fish_fleets)),
                                        Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                        move_timing = move_timing
        )

        for(p in 1:n_pop) {
          for(r in 1:n_regions) {
            for(s in 1:n_sexes) {
              if(y <= n_est_rec_devs) tmp_total_rec <- tmp_Det_Rec[p,r] * exp(ln_RecDevs[p,r,y])
              if(y > n_est_rec_devs) tmp_total_rec <- tmp_Det_Rec[p,r]
              NAA[p,r,y,spawn_seas,1,s]  <- tmp_total_rec * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,s]
              NAA0[p,r,y,spawn_seas,1,s] <- tmp_total_rec * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,s]
            }
            Rec[p,r,y] <- tmp_total_rec # store total before seasonal split
          } # end r loop
        } # end p loop

        # Recruits just inserted above missed this season's movement step. Move recruits if allowed.
        # Only needed under move_timing == 0, where movement happens at the start of the
        # season; under timings 1 and 2 these recruits are picked up by the end-of-season step.
        if(do_recruits_move == 1 && n_regions > 1 && move_timing == 0) {
          for(p in 1:n_pop) {
            for(s in 1:n_sexes) {
              NAA[p,,y,seas,1,s] <- t(NAA[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
              NAA0[p,,y,seas,1,s] <- t(NAA0[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
            } # end s loop
          } # end p loop
          NAA_aft[,,y,seas,1,] <- NAA[,,y,seas,1,]
        }

        # Recompute biomass quantities now that this year's recruits are included
        biom <- compute_biom_y(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, seasdur,
                              n_seas, n_pop, n_regions, n_ages, n_sexes,
                              sgl_seas_spawning_movement, natal_region, stray_rate,
                              Movement, Mrate, move_timing, do_recruits_move)
        Total_Biom[,, y] <- biom$Total_Biom_y
        SSB[,, y] <- biom$SSB_y
        Dynamic_SSB0[,,y] <- biom$Dynamic_SSB0_y
        eff_SSB[,y] <- biom$eff_SSB_y

      } # end if rec_lag == 0 && seas == spawn_seas

      ### Movement (timing 1 and 2), Mortality and Ageing or Continuous Movement and Mortality ---------------------------
      # Post-season state at every age, before the ageing shift. Under move_timing == 0
      # movement was already applied at the top of the season, so this reduces to the
      # original elementwise survival; under move_timing 1 and 2 the seasonal transition
      # operator carries movement and mortality together.
      if(move_timing == 0 || n_regions == 1) {
        step_NAA <- array(NAA[,,y,seas,1:n_ages,] * exp(-ZAA[,,y,seas,1:n_ages,]),
                          dim = c(n_pop, n_regions, n_ages, n_sexes))
        step_NAA0 <- array(NAA0[,,y,seas,1:n_ages,] * exp(-(natmort[,,y,1:n_ages,] * seasdur[seas])),
                           dim = c(n_pop, n_regions, n_ages, n_sexes))
        if(move_timing == 2) {
          for(p in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes) {
            NAA_int[p,,y,seas,a,s] <- integrate_seas_abundance(NAA[p,,y,seas,a,s], ZAA[p,,y,seas,a,s],
                                                               Mrate[p,,,y,seas,a,s], seasdur[seas])
          }
        }
      } else {
        step_NAA <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        step_NAA0 <- array(0, dim = c(n_pop, n_regions, n_ages, n_sexes))
        for(p in 1:n_pop) {
          for(a in 1:n_ages) {
            # Recruits only move when allowed; otherwise an identity transition and a zero
            # generator leave survival unchanged under either timing.
            moves <- (do_recruits_move == 1 || a > 1)
            for(s in 1:n_sexes) {
              Mv <- if(moves) Movement[p,,,y,seas,a,s] else diag(n_regions)
              Qv <- if(moves) Mrate[p,,,y,seas,a,s] else matrix(0, n_regions, n_regions)
              if(move_timing == 2) {
                # The fished stratum needs both the transition operator and the catch
                # integral over the same generator, so take them from one exponential.
                both <- seas_operator_and_integral(ZAA[p,,y,seas,a,s], Qv, seasdur[seas])
                step_NAA[p,,a,s] <- as.vector(t(NAA[p,,y,seas,a,s]) %*% both$T)
                NAA_int[p,,y,seas,a,s] <- as.vector(both$Integral %*% NAA[p,,y,seas,a,s])
              } else {
                step_NAA[p,,a,s] <- advance_seas(NAA[p,,y,seas,a,s], Mv, ZAA[p,,y,seas,a,s],
                                                 Qv, seasdur[seas], move_timing)
              }
              step_NAA0[p,,a,s] <- advance_seas(NAA0[p,,y,seas,a,s], Mv, natmort[p,,y,a,s] * seasdur[seas],
                                                Qv, seasdur[seas], move_timing)
            } # end s loop
          } # end a loop
        } # end p loop
        # Movement happens at season end under these timings, so record it here
        NAA_aft[,,y,seas,,] <- step_NAA
      }

      if(seas < n_seas) {
        # within year / seasonal mortality
        NAA[,,y,seas+1,1:n_ages,] <- step_NAA
        NAA0[,,y,seas+1,1:n_ages,] <- step_NAA0
      } else {
        # age advancement and enter into first season of next year
        # Fished
        NAA[,,y+1,1,2:n_ages,] <- step_NAA[,,1:(n_ages-1),] # Exponential mortality for individuals not in plus group
        NAA[,,y+1,1,n_ages,] <- NAA[,,y+1,1,n_ages,] + step_NAA[,,n_ages,] # Acuumulate plus group
        # Unfished
        NAA0[,,y+1,1,2:n_ages,] <- step_NAA0[,,1:(n_ages-1),] # Exponential mortality for individuals not in plus group
        NAA0[,,y+1,1,n_ages,] <- NAA0[,,y+1,1,n_ages,] + step_NAA0[,,n_ages,] # Acuumulate plus group
      }

      ### Compute Biomass Quantities (rec_lag != 0)
      if(rec_lag != 0 && seas == spawn_seas) {
        spawn_biom <- compute_biom_y(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, seasdur,
                                    n_seas, n_pop, n_regions, n_ages, n_sexes,
                                    sgl_seas_spawning_movement, natal_region, stray_rate,
                              Movement, Mrate, move_timing, do_recruits_move)
        Total_Biom[,, y] <- spawn_biom$Total_Biom_y
        SSB[,, y] <- spawn_biom$SSB_y
        Dynamic_SSB0[,,y] <- spawn_biom$Dynamic_SSB0_y
        eff_SSB[,y] <- spawn_biom$eff_SSB_y
      }

    } # end seas loop
  } # end y loop

  # Get aggregated SSB values
  Aggregated_SSB <- apply(SSB, 3, sum)
  Dynamic_Aggregated_SSB0 <- apply(Dynamic_SSB0, 3, sum)

  return(list(NAA = NAA, NAA0 = NAA0, NAA_bef = NAA_bef, NAA_aft = NAA_aft, Rec = Rec,
              SSB = SSB, Total_Biom = Total_Biom, Dynamic_SSB0 = Dynamic_SSB0, eff_SSB = eff_SSB,
              Aggregated_SSB = Aggregated_SSB, Dynamic_Aggregated_SSB0 = Dynamic_Aggregated_SSB0,
              NAA_int = NAA_int))
}
