# Stage 2 of 3: objective function
#
# The forward projection at the center of the objective function. Walks the
# population through every year and season, applying recruitment, mortality and
# movement in the order set by move_timing, and records the state each
# observation model needs.

#' Gather the arguments \code{\link{compute_mortality_year}} needs from the model frame
#'
#' @param env Environment holding the unpacked data and parameters, i.e. the
#'   \code{SPoRC_rtmb} frame after \code{RTMB::getAll}. Defaults to the caller.
#'
#' @return A named list of arguments for \code{\link{compute_mortality_year}}.
#'
#' @keywords internal
mortality_args_from_model = function(env = parent.frame()) {
  nm = c("growth_model", "derive_waa", "fish_selex_type", "ret_selex_type", "srv_selex_type",
         "fish_waa_selected", "srv_waa_selected", "fish_sel_l", "ret_sel_l", "srv_sel_l",
         "wt_len_pars", "growth_len_mid_vals",
         "UseCatch", "UseCatch_pop", "missing_catch", "UseCatchAA", "UseCatchAA_pop", "use_catch_aa",
         "ln_F_mean", "ln_F_devs", "logit_dmr_mean", "logit_dmr_devs",
         "SizeAgeTrans", "natmort", "seasdur",
         "n_pop", "n_regions", "n_seas", "n_ages", "n_sexes", "n_fish_fleets")
  stats::setNames(lapply(nm, function(x) get(x, envir = env)), nm)
}

#' Fishing and total mortality for one year
#'
#' Also derives the selection-weighted weight at age where a fleet asks for
#' it. Both are done a year at a time because under cohort growth the key a
#' length-based selectivity acts through is only known once the population
#' loop reaches that year; every other model runs them for all years before
#' the loop starts. Under cohort growth this runs inside the population loop,
#' so the year's state is taken as an argument and handed back rather than
#' assigned into this frame.
#'
#' @param y Year index.
#' @param st Named list carrying \code{Fmort}, \code{dmr}, \code{fish_sel},
#'   \code{ret_sel}, \code{ret_FAA}, \code{disc_FAA}, \code{tot_FAA},
#'   \code{ZAA}, \code{WAA_fish}, \code{WAA_srv}, \code{SizeAgeTrans_fish} and
#'   \code{SizeAgeTrans_srv}, returned with year \code{y} updated. The last
#'   two may be \code{NULL}, in which case the shared \code{SizeAgeTrans} key
#'   is used instead.
#' @param growth_model,derive_waa Growth module switches.
#' @param fish_selex_type,ret_selex_type,srv_selex_type Integer (0/1);
#'   \code{1} for length-based selectivity.
#' @param fish_waa_selected,srv_waa_selected Integer vectors \code{[fleet]}
#'   (0/1) for which fleets get a selection-weighted weight at age.
#' @param fish_sel_l,ret_sel_l,srv_sel_l Arrays \code{[region, year, len,
#'   sex, fleet]} of selectivity at length.
#' @param wt_len_pars Array \code{[pop, region, sex, 2]} of weight-length
#'   parameters.
#' @param growth_len_mid_vals Length bin midpoints the weight-length
#'   relationship is read at.
#' @param UseCatch,UseCatch_pop Arrays flagging which cells fit an
#'   aggregate/pop-specific catch observation.
#' @param missing_catch Logical array, \code{TRUE} where the aggregate catch
#'   observation is missing (not a true recorded zero).
#' @param UseCatchAA,UseCatchAA_pop Integer arrays \code{[n_regions, n_years,
#'   n_seas, n_ages, n_sexes, n_fish_fleets]}, with a leading population
#'   dimension for the second, non-zero where a catch at age observation is fit.
#'   A fleet fitting catch at age is fished in any cell where at least one age is
#'   fit, in either stream.
#' @param use_catch_aa Integer vector \code{[n_fish_fleets]}, non-zero for fleets
#'   fitting catch at age rather than aggregated catch.
#' @param ln_F_mean,ln_F_devs Log fishing mortality mean and deviations.
#' @param logit_dmr_mean,logit_dmr_devs Logit discard mortality rate mean and
#'   deviations.
#' @param SizeAgeTrans Shared size-age key, used when no growth-derived
#'   per-fleet key is supplied in \code{st}.
#' @param natmort,seasdur Natural mortality at age and season duration.
#' @param n_pop,n_regions,n_seas,n_ages,n_sexes,n_fish_fleets Dimensions.
#'
#' @return \code{st} with year \code{y} updated.
#'
#' @keywords internal
#' @import RTMB
compute_mortality_year = function(y, st, growth_model, derive_waa, fish_selex_type, ret_selex_type, srv_selex_type,
                          fish_waa_selected, srv_waa_selected, fish_sel_l, ret_sel_l, srv_sel_l,
                          wt_len_pars, growth_len_mid_vals,
                          UseCatch, UseCatch_pop, missing_catch, UseCatchAA, UseCatchAA_pop, use_catch_aa,
                          ln_F_mean, ln_F_devs, logit_dmr_mean, logit_dmr_devs,
                          SizeAgeTrans, natmort, seasdur,
                          n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  Fmort = st$Fmort; dmr = st$dmr; fish_sel = st$fish_sel; ret_sel = st$ret_sel
  ret_FAA = st$ret_FAA; disc_FAA = st$disc_FAA; tot_FAA = st$tot_FAA; ZAA = st$ZAA
  WAA_fish = st$WAA_fish; SizeAgeTrans_fish = st$SizeAgeTrans_fish
  WAA_srv = st$WAA_srv; SizeAgeTrans_srv = st$SizeAgeTrans_srv

  if(growth_model != 0 && derive_waa == 1 && fish_selex_type == 1 && any(fish_waa_selected == 1)) {
    WAA_fish = growth_selected_waa_year(WAA_fish, SizeAgeTrans_fish, fish_sel_l, wt_len_pars,
                                        growth_len_mid_vals, fish_waa_selected, y,
                                        n_pop, n_regions, n_seas, n_sexes)
  }

  if(growth_model != 0 && derive_waa == 1 && srv_selex_type == 1 && any(srv_waa_selected == 1)) {
    WAA_srv = growth_selected_waa_year(WAA_srv, SizeAgeTrans_srv, srv_sel_l, wt_len_pars,
                                       growth_len_mid_vals, srv_waa_selected, y,
                                       n_pop, n_regions, n_seas, n_sexes)
  }

  for(r in 1:n_regions) {
    for(seas in 1:n_seas) {
      for(f in 1:n_fish_fleets) {

        # A cell is a true closure only when no catch is fit
        # A fleet fitting catch at age is fished wherever any age is fit there
        caa_open = if(use_catch_aa[f] == 1) any(UseCatchAA[r,y,seas,,,f] == 1) ||
                                            any(UseCatchAA_pop[,r,y,seas,,,f] == 1) else FALSE
        is_closed = (UseCatch[r,y,seas,f] == 0) && all(UseCatch_pop[,r,y,seas,f] == 0) &&
          !missing_catch[r,y,seas,f] && !caa_open

        if(is_closed) {
          Fmort[r,y,seas,f] = 0
          dmr[r,y,seas,f] = 0
        } else {
          Fmort[r,y,seas,f] = exp(ln_F_mean[r,seas,f] + ln_F_devs[r,y,seas,f])
          dmr[r,y,seas,f] = RTMB::plogis(logit_dmr_mean[r,seas,f] + logit_dmr_devs[r,y,seas,f])
        }

        # get fishing mortality at age
        for(p in 1:n_pop) {
          # the key is the fleet's own when the growth module built it, the shared data key otherwise
          if(fish_selex_type == 1 || ret_selex_type == 1) {
            for(s in 1:n_sexes) {
              key_f = if(is.null(SizeAgeTrans_fish)) SizeAgeTrans[p,r,y,seas,,,s] else SizeAgeTrans_fish[p,r,y,seas,,,s,f]
              if(fish_selex_type == 1) fish_sel[p,r,y,seas,,s,f] = fish_sel_l[r,y,,s,f] %*% key_f
              if(ret_selex_type == 1) ret_sel[p,r,y,seas,,s,f] = ret_sel_l[r,y,,s,f] %*% key_f
            } # end s loop

          } # length based selectivity
          ret_FAA[p,r,y,seas,,,f] = Fmort[r,y,seas,f] * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f] # Retained fishing mortality at age
          disc_FAA[p,r,y,seas,,,f] = Fmort[r,y,seas,f] * fish_sel[p,r,y,seas,,,f] * (1 - ret_sel[p,r,y,seas,,,f]) * dmr[r,y,seas,f] # Discarded fishing mortality at age
          tot_FAA[p,r,y,seas,,,f] = ret_FAA[p,r,y,seas,,,f] +  disc_FAA[p,r,y,seas,,,f]# Total fishing mortality at age

        } # end p loop

      } # end f loop

      # get total mortality
      for(p1 in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes)
        ZAA[p1,r,y,seas,a,s] = sum(ret_FAA[p1,r,y,seas,a,s,]) + sum(disc_FAA[p1,r,y,seas,a,s,]) + (natmort[p1,r,y,a,s] * seasdur[seas])
    } # end seas loop
  } # end r loop

  st$Fmort = Fmort; st$dmr = dmr; st$fish_sel = fish_sel; st$ret_sel = ret_sel
  st$ret_FAA = ret_FAA; st$disc_FAA = disc_FAA; st$tot_FAA = tot_FAA; st$ZAA = ZAA
  st$WAA_fish = WAA_fish; st$WAA_srv = WAA_srv

  return(st)

} # end compute_mortality_year

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
#' @param n_pop,n_regions,n_seas,n_ages,n_sexes,n_yrs,n_fish_fleets
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
#' @param SR_ref_yr Integer year index supplying the biological inputs, weight
#'   at age, maturity, natural mortality and movement, to unfished spawning
#'   biomass per recruit, and so to \code{S0} and the scale of the stock-recruit
#'   curve. Default \code{1}, the first model year, which is what the function
#'   used to hardcode. Set to \code{n_yrs} to condition the curve on terminal
#'   weight at age, which is what several ADMB assessments do; with time-varying
#'   weight at age the two differ and the whole curve shifts with them. It is a
#'   year INDEX, not a calendar year, so callers that truncate the year
#'   dimension (retrospectives) must clamp it.
#' @param growth_mortality_year_fn Optional function of \code{(y, NAA_y, growth_mortality_state)} called at
#'   the top of every year with the numbers at age at the start of that year,
#'   array \code{[pop, region, age, sex]}, and the state carried from the
#'   previous year. It returns a list with \code{state}, carried forward to the
#'   next call and returned to the caller, and \code{ZAA_y}, \code{WAA_y} and
#'   \code{MatAA_y}, the year's slices of total mortality, weight and maturity
#'   at age, which replace those handed in for that year. Passing the state in
#'   and out keeps the per-year step a function of its arguments.
#' @param growth_mortality_state Initial state for \code{growth_mortality_year_fn}, carried through the
#'   year loop and returned as \code{growth_mortality_state}. Ignored when
#'   \code{growth_mortality_year_fn} is \code{NULL}.
#'   This is how cohort growth, whose plus group blends by numbers, is evaluated
#'   inside the year loop. \code{NULL} (the default) uses the arrays as given.
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
                                       Mrate = NULL, move_timing = 0, SR_ref_yr = 1,
                                       sr_penalty = 0, sr_R0 = NULL, growth_mortality_year_fn = NULL, growth_mortality_state = NULL,
                                       expm_nsub = 0) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  # Season-integrated abundance for the spatial Baranov, filled in below only under move_timing == 2.
  NAA_int <- array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes))

  # Stock recruitment prediction - only used when mean recruitment, but penalize the mean recruits to an SR curve
  SR_pred <- array(1, dim = c(n_pop, n_regions, n_yrs))

  # Arguments for determinstic recruitment
  det_rec_args <- list(
    rec_dd = rec_dd,
    rec_region_prop = rec_region_prop,
    rec_seas_prop = rec_seas_prop,
    h = h_trans,
    n_pop = n_pop,
    n_ages = n_ages,
    n_regions = n_regions,
    sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,1,1], dim = c(n_pop, n_regions)),
    WAA = array(WAA[,,SR_ref_yr,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    MatAA = array(MatAA[,,SR_ref_yr,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    natmort = array(natmort[,,SR_ref_yr,,1], dim = c(n_pop, n_regions, n_ages)),
    Movement = array(Movement[,,,SR_ref_yr,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
    stray_rate = array(stray_rate[,1], dim = c(n_pop)),
    sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,SR_ref_yr,,1], dim = c(n_pop, n_regions, n_regions, n_ages)),
    do_recruits_move = do_recruits_move,
    natal_region = natal_region,
    t_spawn = t_spawn,
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

  for(y in 1:n_yrs) {

    # Growth carried cohort by cohort needs this year's start-of-year numbers
    # before anything else in the year is formed
    if(!is.null(growth_mortality_year_fn)) {
      hk <- growth_mortality_year_fn(y, array(NAA[,,y,1,,], dim = c(n_pop, n_regions, n_ages, n_sexes)), growth_mortality_state)
      growth_mortality_state <- hk$state
      ZAA[,,y,,,] <- hk$ZAA_y
      WAA[,,y,,,] <- hk$WAA_y
      MatAA[,,y,,,] <- hk$MatAA_y
    }

    ### Annual Recruitment (rec_lag != 0 only) -----------------------------------
    if(rec_lag != 0) {

      # Get Deterministic Recruitment
      tmp_Det_Rec <- do.call(Get_Det_Recruitment, c(det_rec_args, list(recruitment_model = rec_model, R0 = R0, SSB_vals = SSB, y = y)))

      # Predicts the SR curve for use as a penalty 
      if(sr_penalty > 0) {
        SR_pred[,, y] <- do.call(Get_Det_Recruitment, c(det_rec_args, list(recruitment_model = sr_penalty, R0 = sr_R0, SSB_vals = SSB, y = y)))
      }

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
                              Movement, Mrate, move_timing, do_recruits_move, expm_nsub = expm_nsub)
        SSB[,, y] <- biom$SSB_y

        # Deterministic recruitment
        tmp_Det_Rec <- do.call(Get_Det_Recruitment, c(det_rec_args, list(recruitment_model = rec_model, R0 = R0, SSB_vals = SSB, y = y)))

        # stock recruitment penalty against a curve
        if(sr_penalty > 0) {
          SR_pred[,, y] <- do.call(Get_Det_Recruitment, c(det_rec_args, list(recruitment_model = sr_penalty, R0 = sr_R0, SSB_vals = SSB, y = y)))
        }

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
                              Movement, Mrate, move_timing, do_recruits_move, expm_nsub = expm_nsub)
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
                                                               Mrate[p,,,y,seas,a,s], seasdur[seas], expm_nsub = expm_nsub)
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
                both <- seas_operator_and_integral(ZAA[p,,y,seas,a,s], Qv, seasdur[seas], expm_nsub = expm_nsub)
                step_NAA[p,,a,s] <- as.vector(t(NAA[p,,y,seas,a,s]) %*% both$T)
                NAA_int[p,,y,seas,a,s] <- as.vector(both$Integral %*% NAA[p,,y,seas,a,s])
              } else {
                step_NAA[p,,a,s] <- advance_seas(NAA[p,,y,seas,a,s], Mv, ZAA[p,,y,seas,a,s],
                                                 Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
              }
              step_NAA0[p,,a,s] <- advance_seas(NAA0[p,,y,seas,a,s], Mv, natmort[p,,y,a,s] * seasdur[seas],
                                                Qv, seasdur[seas], move_timing, expm_nsub = expm_nsub)
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
                              Movement, Mrate, move_timing, do_recruits_move, expm_nsub = expm_nsub)
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
              NAA_int = NAA_int, SR_pred = SR_pred, growth_mortality_state = growth_mortality_state))
}
