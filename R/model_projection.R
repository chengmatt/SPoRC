# Stage 2 of 3: objective function
#
#' Threshold Harvest Control Rule
#'
#' The target F above the trigger biomass \code{brp}, no fishing below \code{alpha}
#' times it, and a straight ramp between. Status is compared on its plain numbers
#' through \code{\link{ad_value}}, so the rule can run inside the model despite its
#' \code{if}. The model keeps the side of the comparison that was true when it was built,
#' and \code{\link{fit_model}} builds it at the estimates, so that side is the
#' estimated stock's.
#'
#' @param x Spawning biomass the rule reads.
#' @param frp Fishing mortality reference point.
#' @param brp Biomass reference point.
#' @param alpha Share of \code{brp} below which F is zero. Default 0.05.
#'
#' @return The F to apply. Inside the model it still depends on \code{x} and
#'   \code{frp}, so \code{sdreport} accounts for their uncertainty.
#' @export hcr_threshold
#' @family Reference Points and Projections
#' @examples
#' hcr_threshold(x = 80, frp = 0.1, brp = 100) # on the ramp: 0.1 * (0.8 - 0.05) / 0.95
hcr_threshold <- function(x, frp, brp, alpha = 0.05) {
  status <- ad_value(x) / ad_value(brp) # compare status
  if(status >= 1) return(frp) # above the trigger, the target F
  if(status <= alpha) return(0 * frp) # below the cutoff, no fishing
  return(frp * (x / brp - alpha) / (1 - alpha)) # the ramp between
}


#' Average An Array Over Some Of Its Years
#'
#' Keeps demographics at their recent average in the projection: weight, maturity, natural
#' mortality, selectivity and movement averaged over the last \code{n_avg_yrs} years. The array
#' can have any number of dims, so a year's slice is taken by building the index for \code{[} as
#' a list, every dim kept whole except the year dim, which picks one year. It uses only indexing
#' and addition, so it runs on the model's own arrays inside the objective.
#'
#' @param x Array with a year dim.
#' @param yr_dim Integer. Which dim of \code{x} is years.
#' @param yrs Integer vector. The years averaged, as positions along \code{yr_dim}.
#' @return Array shaped as \code{x}, with the year dim kept at length one so
#'   \code{repeat_proj_years} can spread it over the projection years.
#' @keywords internal
#' @noRd
avg_over_years <- function(x, yr_dim, yrs) {

  # every dim kept whole (TRUE) except the year dim, which is set to one year at a time
  year_slice_index <- rep(list(TRUE), length(dim(x)))

  # x at one year, that year kept as a dim of length one
  year_slice <- function(year) {
    year_slice_index[[yr_dim]] <- year
    do.call("[", c(list(x), year_slice_index, list(drop = FALSE)))
  }

  # add the years one slice at a time, then divide by how many there were
  year_sum <- year_slice(yrs[1])
  for(year in yrs[-1]) year_sum <- year_sum + year_slice(year)

  return(year_sum / length(yrs))
}


#' Repeat One Year's Values Over The Projection Years
#'
#' Fills every projection year with the same values, the output of \code{avg_over_years} or a
#' terminal year slice.
#'
#' @param one_yr Array whose year dim has length one.
#' @param yr_dim Integer. Which dim is years.
#' @param n_proj_yrs Integer. Projection years.
#' @return Array shaped as \code{one_yr} with \code{n_proj_yrs} years.
#' @keywords internal
#' @noRd
repeat_proj_years <- function(one_yr, yr_dim, n_proj_yrs) {

  "[<-" <- RTMB::ADoverload("[<-") # so model quantities can be written into the plain array below

  # one_yr's shape with the year dim stretched to the projection years
  out_dim <- dim(one_yr)
  out_dim[yr_dim] <- n_proj_yrs
  out <- array(0, dim = out_dim)

  # write one_yr into each projection year in turn, every other dim whole
  proj_year_index <- rep(list(TRUE), length(out_dim))
  for(y in 1:n_proj_yrs) {
    proj_year_index[[yr_dim]] <- y
    out <- do.call("[<-", c(list(out), proj_year_index, list(value = one_yr)))
  } # end y loop

  return(out)
}


#' The Terminal Year Then The Years The Fit Projected
#'
#' For a process the fit itself projected past the terminal year through \code{n_proj_yrs_devs}
#' (time-varying selectivity, movement or growth). Projection year 1 takes the terminal year,
#' and projection year \code{y} the fit's year \code{n_yrs + y - 1}.
#'
#' @param x Array whose year dim runs past the fitted years.
#' @param yr_dim Integer. Which dim is years.
#' @param n_yrs Integer. Fitted years.
#' @param n_proj_yrs Integer. Projection years, the first being the terminal year.
#' @return Array shaped as \code{x} with \code{n_proj_yrs} years.
#' @keywords internal
#' @noRd
extended_proj_years <- function(x, yr_dim, n_yrs, n_proj_yrs) {

  "[<-" <- RTMB::ADoverload("[<-") # so model quantities can be written into the plain array below
  n_yrs_past_fit <- dim(x)[yr_dim] - n_yrs
  if(n_yrs_past_fit < n_proj_yrs - 1) stop("A process projected by the fit needs n_proj_yrs_devs of at least n_proj_yrs - 1 in Setup_Mod_Dim, but this array has ", n_yrs_past_fit, " years past the fit.", call. = FALSE)

  # x's shape with the year dim cut to the projection years
  out_dim <- dim(x)
  out_dim[yr_dim] <- n_proj_yrs
  out <- array(0, dim = out_dim)

  # copy fit year n_yrs + y - 1 into projection year y, every other dim whole
  fit_year_index <- proj_year_index <- rep(list(TRUE), length(out_dim))
  for(y in 1:n_proj_yrs) {
    fit_year_index[[yr_dim]] <- n_yrs + y - 1 # the terminal year, then the years the fit projected
    proj_year_index[[yr_dim]] <- y
    fit_year <- do.call("[", c(list(x), fit_year_index, list(drop = FALSE)))
    out <- do.call("[<-", c(list(out), proj_year_index, list(value = fit_year)))
  } # end y loop

  return(out)
}


#' Run The Projection Inside The Objective
#'
#' Builds what \code{\link{Do_Population_Projection}} takes from the model's own arrays, the terminal
#' numbers and F, demographics averaged over the terminal years or read from the years the fit
#' extended, and mean recruitment over the chosen years times the fit's projected deviations, then
#' runs it under the F option, rule and constraints \code{\link{Setup_Mod_Projection}} stored.
#'
#' @param proj_n_yrs,proj_f_opt,proj_f_input,proj_catch_input,proj_f_ref_pt,proj_b_ref_pt,proj_hcr_function,proj_hcr_scope,proj_F_constraint,proj_n_avg_yrs,proj_sel_opt,proj_move_opt,proj_growth_opt,proj_calc_rec_yrs,proj_rec_opt,proj_naa_opt
#'   The settings as \code{Setup_Mod_Projection} stores them.
#' @param n_yrs,n_pop,n_regions,n_ages,n_sexes,n_seas,n_fish_fleets Dimensions.
#' @param seasdur,spawn_seas,natal_region,do_recruits_move,move_timing,expm_nsub,t_spawn Model settings.
#' @param NAA,NAA0,Fmort,Rec,natmort,WAA,WAA_fish,MatAA,fish_sel,ret_sel,dmr,Movement,Mrate,sgl_seas_spawning_movement,sexratio,stray_rate,rec_seas_prop,rec_region_prop,R0,h_trans,ln_RecDevs
#'   The model's arrays, on the tape.
#' @param ln_NAA,n_est_naa_re,naa_re_ages,naa_re_seas,naa_re_where The numbers at age state and where it
#'   runs. Cells past the fitted years hold the projected innovations, which multiply the projected numbers
#'   where the state applies.
#'
#' @return Named list: \code{proj_SSB}, \code{proj_Catch}, \code{proj_Catch_reg} (by region and year,
#'   populations, seasons and fleets summed), \code{proj_F}, \code{proj_F_rule} and \code{proj_Total_Biom}.
#' @keywords internal
#' @noRd
get_internal_projection <- function(proj_n_yrs, proj_f_opt, proj_f_input, proj_f_ref_pt, proj_b_ref_pt,
                                    proj_hcr_function, proj_hcr_scope, proj_F_constraint, proj_catch_input,
                                    proj_refpt_type, proj_refpt_what, proj_spr_x, proj_sex_ratio_f, proj_is_discard_fleet,
                                    proj_local_bh_msy_newton_steps, proj_refpt_newton_steps, proj_refpt_log_f_start,
                                    proj_n_avg_yrs, proj_sel_opt, proj_move_opt, proj_growth_opt, proj_calc_rec_yrs, proj_rec_opt, proj_naa_opt,
                                    n_yrs, n_pop, n_regions, n_ages, n_sexes, n_seas, n_fish_fleets,
                                    seasdur, spawn_seas, natal_region, do_recruits_move, move_timing, expm_nsub, t_spawn, rec_model,
                                    NAA, NAA0, Fmort, Rec, natmort, WAA, WAA_fish, MatAA, fish_sel, ret_sel, dmr,
                                    Movement, Mrate, sgl_seas_spawning_movement, sexratio, stray_rate, rec_seas_prop, rec_region_prop,
                                    R0, h_trans, ln_RecDevs, ln_NAA, n_est_naa_re, naa_re_ages, naa_re_seas, naa_re_where) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  avg_yrs <- (n_yrs - proj_n_avg_yrs + 1):n_yrs # terminal years averaged

  # reference points solved on the tape, when the rule reads them
  refpts <- NULL
  if(proj_f_opt == 4) {

    # the female share of recruits the reference points start from: the model's own unless one was given, so an
    # estimated sex ratio and its uncertainty reach them. a one-sex model halves its spawning biomass, so 0.5
    sex_ratio_f <- proj_sex_ratio_f
    if(is.null(sex_ratio_f)) {
      if(n_sexes == 2) {
        sex_ratio_f <- avg_over_years(sexratio[,,,1,drop = FALSE], 3, avg_yrs) # [p, r, 1, 1]
        dim(sex_ratio_f) <- c(n_pop, n_regions)
      } else sex_ratio_f <- array(0.5, dim = c(n_pop, n_regions))
    }

    refpts <- get_refpts_tape(refpt_type = proj_refpt_type, refpt_what = proj_refpt_what, SPR_x = proj_spr_x,
                              sex_ratio_f = sex_ratio_f, is_discard_fleet = proj_is_discard_fleet,
                              local_bh_msy_newton_steps = proj_local_bh_msy_newton_steps,
                              newton_steps = proj_refpt_newton_steps, log_F_start = proj_refpt_log_f_start,
                              avg_yrs = avg_yrs, rec_yrs = proj_calc_rec_yrs,
                              n_yrs = n_yrs, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages, n_seas = n_seas, n_fish_fleets = n_fish_fleets,
                              seasdur = seasdur, spawn_seas = spawn_seas, t_spawn = t_spawn, natal_region = natal_region,
                              do_recruits_move = do_recruits_move, move_timing = move_timing, expm_nsub = expm_nsub, rec_model = rec_model,
                              Fmort = Fmort, dmr = dmr, fish_sel = fish_sel, ret_sel = ret_sel, natmort = natmort, WAA = WAA, MatAA = MatAA,
                              Movement = Movement, Mrate = Mrate, sgl_seas_spawning_movement = sgl_seas_spawning_movement,
                              stray_rate = stray_rate, rec_seas_prop = rec_seas_prop, rec_region_prop = rec_region_prop,
                              Rec = Rec, R0 = R0, h_trans = h_trans)
  }

  # terminal state
  terminal_NAA <- NAA[,,n_yrs,,,,drop = FALSE]
  dim(terminal_NAA) <- c(n_pop, n_regions, n_seas, n_ages, n_sexes)
  terminal_NAA0 <- NAA0[,,n_yrs,,,,drop = FALSE]
  dim(terminal_NAA0) <- c(n_pop, n_regions, n_seas, n_ages, n_sexes)
  terminal_F <- Fmort[,n_yrs,,,drop = FALSE]
  dim(terminal_F) <- c(n_regions, n_seas, n_fish_fleets)
  terminal_dmr <- dmr[,n_yrs,,,drop = FALSE] # discard mortality by region, season and fleet
  dim(terminal_dmr) <- c(n_regions, n_seas, n_fish_fleets)
  recruitment <- Rec[,,proj_calc_rec_yrs,drop = FALSE] # the years mean recruitment is taken over

  # demographics over the projection years
  WAA_p <- if(proj_growth_opt == 1) repeat_proj_years(avg_over_years(WAA, 3, avg_yrs), 3, proj_n_yrs) else extended_proj_years(WAA, 3, n_yrs, proj_n_yrs)
  WAA_fish_p <- if(proj_growth_opt == 1) repeat_proj_years(avg_over_years(WAA_fish, 3, avg_yrs), 3, proj_n_yrs) else extended_proj_years(WAA_fish, 3, n_yrs, proj_n_yrs)
  MatAA_p <- repeat_proj_years(avg_over_years(MatAA, 3, avg_yrs), 3, proj_n_yrs)
  natmort_p <- repeat_proj_years(avg_over_years(natmort, 3, avg_yrs), 3, proj_n_yrs)
  sexratio_p <- repeat_proj_years(avg_over_years(sexratio, 3, avg_yrs), 3, proj_n_yrs)
  stray_rate_p <- repeat_proj_years(avg_over_years(stray_rate, 2, avg_yrs), 2, proj_n_yrs)
  sgl_seas_p <- repeat_proj_years(sgl_seas_spawning_movement[,,,n_yrs,,,drop = FALSE], 4, proj_n_yrs) # terminal year, as the fit used it
  fish_sel_p <- if(proj_sel_opt == 1) repeat_proj_years(avg_over_years(fish_sel, 3, avg_yrs), 3, proj_n_yrs) else extended_proj_years(fish_sel, 3, n_yrs, proj_n_yrs)
  ret_sel_p <- if(proj_sel_opt == 1) repeat_proj_years(avg_over_years(ret_sel, 3, avg_yrs), 3, proj_n_yrs) else extended_proj_years(ret_sel, 3, n_yrs, proj_n_yrs)
  Movement_p <- if(proj_move_opt == 1) repeat_proj_years(avg_over_years(Movement, 4, avg_yrs), 4, proj_n_yrs) else extended_proj_years(Movement, 4, n_yrs, proj_n_yrs)
  Mrate_p <- NULL
  if(move_timing == 2) Mrate_p <- if(proj_move_opt == 1) repeat_proj_years(avg_over_years(Mrate, 4, avg_yrs), 4, proj_n_yrs) else extended_proj_years(Mrate, 4, n_yrs, proj_n_yrs)

  # fishing mortality by region and year, or the reference points the rule reads
  f_proj <- array(0, dim = c(n_regions, proj_n_yrs))
  b_proj <- array(0, dim = c(n_pop, n_regions, proj_n_yrs))
  if(proj_f_opt == 1) for(r in 1:n_regions) for(seas in 1:n_seas) for(f in 1:n_fish_fleets) f_proj[r,] <- f_proj[r,] + Fmort[r,n_yrs,seas,f] # terminal total F
  if(proj_f_opt == 2) f_proj <- f_proj + proj_f_input
  if(proj_f_opt == 3) {
    f_proj <- f_proj + proj_f_ref_pt
    b_proj <- b_proj + proj_b_ref_pt
  }
  if(proj_f_opt == 4) {
    for(y in 1:proj_n_yrs) {
      f_proj[,y] <- refpts$f_ref_pt # the same reference points in every projection year
      b_proj[,,y] <- refpts$b_ref_pt
    } # end y loop
  }
  rule_opt <- if(proj_f_opt %in% c(3, 4)) (if(proj_hcr_scope == 2) "HCR_global" else "HCR") else "Input"
  fmort_opt <- if(is.null(proj_catch_input)) rule_opt else "Catch" # years with catch given are solved for F, the rest fished under rule_opt
  rule <- if(is.null(proj_hcr_function)) hcr_threshold else proj_hcr_function # input lists stored before the default was the rule itself

  # recruitment in the projection years is mean recruitment times the fit's projected deviations, alone, or
  # times the terminal year's deviation kept; year one gets no new recruits
  rec_devs <- array(1, dim = c(n_pop, n_regions, proj_n_yrs))
  if(proj_rec_opt == 1) for(k in seq_len(proj_n_yrs - 1)) rec_devs[,,k + 1] <- exp(ln_RecDevs[,,n_yrs + k])
  if(proj_rec_opt == 3) for(k in seq_len(proj_n_yrs - 1)) rec_devs[,,k + 1] <- exp(ln_RecDevs[,,n_yrs])

  # the fit's projected state innovations multiply the numbers at age where the state runs, or none
  naa_devs <- NULL
  n_state_ahead <- if(n_est_naa_re > 0) max(0, dim(ln_NAA)[3] - n_yrs) else 0
  if(proj_naa_opt == 1 && n_state_ahead > 0) {
    naa_devs <- array(1, dim = c(n_pop, n_regions, proj_n_yrs, n_seas, n_ages, n_sexes))
    for(k in seq_len(min(proj_n_yrs - 1, n_state_ahead))) {
      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          if(naa_re_where[p,r] == 0) next # a population that never reaches this region has no state here
          for(seas in naa_re_seas) for(a in naa_re_ages) naa_devs[p,r,k + 1,seas,a,] <- exp(ln_NAA[p,r,n_yrs + k,seas,a,])
        } # end r loop
      } # end p loop
    } # end k loop
  }

  out <- Do_Population_Projection(n_proj_yrs = proj_n_yrs,
                                  n_pop = n_pop,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  n_seas = n_seas,
                                  n_fish_fleets = n_fish_fleets,
                                  seasdur = seasdur,
                                  spawn_seas = spawn_seas,
                                  natal_region = natal_region,
                                  do_recruits_move = do_recruits_move,
                                  rec_seas_prop = rec_seas_prop,
                                  t_spawn = t_spawn,
                                  sexratio = sexratio_p,
                                  stray_rate = stray_rate_p,
                                  recruitment = recruitment,
                                  terminal_NAA = terminal_NAA,
                                  terminal_NAA0 = terminal_NAA0,
                                  terminal_F = terminal_F,
                                  natmort = natmort_p,
                                  WAA = WAA_p,
                                  WAA_fish = WAA_fish_p,
                                  MatAA = MatAA_p,
                                  fish_sel = fish_sel_p,
                                  ret_sel = ret_sel_p,
                                  dmr = terminal_dmr,
                                  Movement = Movement_p,
                                  Mrate = Mrate_p,
                                  sgl_seas_spawning_movement = sgl_seas_p,
                                  move_timing = move_timing,
                                  expm_nsub = expm_nsub,
                                  f_ref_pt = f_proj,
                                  b_ref_pt = b_proj,
                                  HCR_function = rule,
                                  recruitment_opt = "mean_rec",
                                  fmort_opt = fmort_opt,
                                  catch_input = proj_catch_input,
                                  catch_fallback_opt = rule_opt,
                                  rec_devs = rec_devs,
                                  naa_devs = naa_devs,
                                  F_constraint = proj_F_constraint)

  # catch by region and year, summed over populations, seasons and fleets
  proj_Catch_reg <- array(0, dim = c(n_regions, proj_n_yrs))
  for(r in 1:n_regions) for(y in 1:proj_n_yrs) proj_Catch_reg[r,y] <- sum(out$proj_Catch[,r,y,,])

  # spawning biomass over its target by region and year, populations summed, what the rule reads
  proj_status <- array(0, dim = c(n_regions, proj_n_yrs))
  if(proj_f_opt %in% c(3, 4)) for(r in 1:n_regions) for(y in 1:proj_n_yrs) proj_status[r,y] <- sum(out$proj_SSB[,r,y]) / sum(b_proj[,r,y])

  return(list(proj_SSB = out$proj_SSB,
              proj_Catch = out$proj_Catch,
              proj_Catch_reg = proj_Catch_reg,
              proj_F = out$proj_F,
              proj_F_rule = out$proj_F_rule,
              proj_Total_Biom = out$proj_Total_Biom,
              proj_NAA = out$proj_NAA,
              proj_status = proj_status,
              refpts = refpts))

} # end get_internal_projection


#' Newton Steps From A Solved Reference Point, So It Moves With The Estimates
#'
#' \code{nlminb} finds the reference point's log F on plain numbers, which leaves it a fixed
#' number to the model. These Newton steps start from that answer and are computed from the
#' model's own quantities: the gradient and Hessian of the objective come from central
#' differences, so each step is arithmetic the model can differentiate. At the minimum the
#' gradient is zero, so a step changes nothing in value, but its derivative is the derivative of
#' the minimizing log F in the model's parameters (the implicit function theorem). One step gives
#' that derivative; the second covers the rounding left by the plain solve.
#'
#' @param objective Function of \code{(log_F, data_list)} returning the scalar to minimize.
#' @param log_F Numeric vector. Where the steps start, the plain solve's answer: one F per
#'   region under local MSY, one shared F otherwise.
#' @param data_list Data the objective reads, the model's own quantities included.
#' @param n_steps Integer. Newton steps.
#' @param h Numeric. Central difference step on \code{log_F}.
#' @param log_F_floor Numeric. A component below this is at zero fishing, where the objective is
#'   flat in log F and a step means nothing, so it is left where the plain solve put it. Default
#'   \code{log(1e-6)}.
#' @return \code{log_F} after the steps.
#' @keywords internal
#' @noRd
refine_refpt_log_F <- function(objective, log_F, data_list, n_steps, h, log_F_floor = log(1e-6)) {

  "[<-" <- RTMB::ADoverload("[<-") # so model quantities can be written into log_F and the arrays below

  # only the F's with fishing are stepped; one left at zero fishing keeps the plain solve's answer
  fished <- which(ad_value(log_F) > log_F_floor)
  n_fished <- length(fished)
  if(n_fished == 0) return(log_F)

  for(step in seq_len(n_steps)) {

    # the objective at the current log F, and its slope and curvature in each fished component
    obj_center <- objective(log_F, data_list)
    gradient <- rep(0, n_fished)
    hessian <- matrix(0, n_fished, n_fished)

    for(i in 1:n_fished) {

      # log F with component i moved up, and down, by h
      log_F_up_i <- log_F
      log_F_up_i[fished[i]] <- log_F_up_i[fished[i]] + h
      log_F_down_i <- log_F
      log_F_down_i[fished[i]] <- log_F_down_i[fished[i]] - h
      obj_up_i <- objective(log_F_up_i, data_list)
      obj_down_i <- objective(log_F_down_i, data_list)

      # central differences: the slope from the two sides, the curvature from the two sides and the center
      gradient[i] <- (obj_up_i - obj_down_i) / (2 * h)
      hessian[i,i] <- (obj_up_i - 2 * obj_center + obj_down_i) / h^2

      # cross curvature with each earlier component j, from moving i and j up or down together
      for(j in seq_len(i - 1)) {
        log_F_up_i_up_j <- log_F_up_i
        log_F_up_i_up_j[fished[j]] <- log_F_up_i_up_j[fished[j]] + h
        log_F_up_i_down_j <- log_F_up_i
        log_F_up_i_down_j[fished[j]] <- log_F_up_i_down_j[fished[j]] - h
        log_F_down_i_up_j <- log_F_down_i
        log_F_down_i_up_j[fished[j]] <- log_F_down_i_up_j[fished[j]] + h
        log_F_down_i_down_j <- log_F_down_i
        log_F_down_i_down_j[fished[j]] <- log_F_down_i_down_j[fished[j]] - h
        hessian[i,j] <- (objective(log_F_up_i_up_j, data_list) - objective(log_F_up_i_down_j, data_list) - objective(log_F_down_i_up_j, data_list) + objective(log_F_down_i_down_j, data_list)) / (4 * h^2)
        hessian[j,i] <- hessian[i,j] # the hessian is symmetric
      } # end j loop

    } # end i loop

    # newton step: the gradient scaled by the inverse curvature, taken off log F
    log_F[fished] <- log_F[fished] - solve(hessian, gradient)

  } # end step loop

  return(log_F)
} # end refine_refpt_log_F


#' Solve A Reference Point Inside The Model
#'
#' Finds the minimizing log F on plain numbers with \code{nlminb}, as
#' \code{\link{Get_Reference_Points}} does, then refines it with \code{refine_refpt_log_F} on the
#' model's own quantities, so the reference point has a derivative with respect to everything in
#' \code{data_list}.
#'
#' @param quantities One of the \code{*_quantities} cores, e.g. \code{single_region_SPR_quantities}.
#' @param data_list Data the core reads, the model's own quantities included.
#' @param log_F_start Numeric vector. Where the plain solve starts.
#' @param newton_steps Integer. Newton steps after the plain solve. Default 2.
#' @param h Numeric. Central difference step on log F. Default 1e-4.
#' @return The core's list at the solution, plus \code{log_F}.
#' @keywords internal
#' @noRd
solve_refpt_tape <- function(quantities, data_list, log_F_start, newton_steps = 2, h = 1e-4) {

  # the minimizing log F on plain numbers. inside the model the core returns model quantities even from plain
  # inputs, so the objective nlminb sees is read through ad_value
  data_plain <- ad_value(data_list)
  objective_plain <- function(log_F) ad_value(quantities(log_F, data_plain)$obj_fun)
  plain_solve <- stats::nlminb(log_F_start, objective_plain, control = list(iter.max = 1e4, eval.max = 1e4, rel.tol = 1e-15))

  # refined with newton steps on the model's own quantities, which gives the reference point its derivative
  objective_model <- function(log_F, model_data) quantities(log_F, model_data)$obj_fun
  log_F <- refine_refpt_log_F(objective_model, plain_solve$par, data_list, newton_steps, h)

  # every quantity the core reports, at the refined log F
  out <- quantities(log_F, data_list)
  out$log_F <- log_F
  return(out)
} # end solve_refpt_tape


#' Reference Points From The Model's Own Arrays
#'
#' Solves the reference points inside the model, so their uncertainty reaches the projection.
#' \code{\link{Get_Reference_Points}} builds each solver's data list from a fitted report; this
#' builds the same list from the model's own arrays while the objective runs, and solves it with
#' \code{solve_refpt_tape}. Demographics are the females' values averaged over \code{avg_yrs}, and a
#' reference F is split among seasons and fleets as the terminal year's F was.
#'
#' A biomass reference point is spawning biomass per recruit at the solved F times a recruitment:
#' mean recruitment over \code{rec_yrs} for the SPR options, the stock-recruit curve's equilibrium
#' recruitment for the MSY options. The unfished (virgin) forms are the same at F = 0. The
#' population forms add, at a population's natal region, what the other populations stray there.
#'
#' @param refpt_type Integer. 1 single region, 2 multi region.
#' @param refpt_what Integer, the \code{Get_Reference_Points} option in order: 1 SPR, 2 MSY,
#'   3 independent_SPR, 4 independent_MSY, 5 global_SPR, 6 global_MSY, 7 local_MSY.
#' @param SPR_x Numeric. Target spawning potential ratio.
#' @param sex_ratio_f Array \code{[n_pop, n_regions]}. Female share of recruits.
#' @param is_discard_fleet Integer vector \code{[n_fish_fleets]}.
#' @param local_bh_msy_newton_steps Integer. Inner steps for \code{"local_MSY"}.
#' @param newton_steps,log_F_start Passed to \code{solve_refpt_tape}.
#' @param avg_yrs Integer vector. Terminal years averaged.
#' @param rec_yrs Integer vector. Years mean recruitment is taken over.
#' @param n_yrs,n_pop,n_regions,n_ages,n_seas,n_fish_fleets,seasdur,spawn_seas,t_spawn,natal_region,do_recruits_move,move_timing,expm_nsub,rec_model
#'   Dimensions and settings.
#' @param Fmort,dmr,fish_sel,ret_sel,natmort,WAA,MatAA,Movement,Mrate,sgl_seas_spawning_movement,stray_rate,rec_seas_prop,rec_region_prop,Rec,R0,h_trans
#'   Model arrays.
#' @return List with \code{f_ref_pt} \code{[n_regions]}, and \code{b_ref_pt}, \code{virgin_b_ref_pt},
#'   \code{pop_b_ref_pt}, \code{virgin_pop_b_ref_pt}, each \code{[n_pop, n_regions]}.
#' @keywords internal
#' @noRd
get_refpts_tape <- function(refpt_type, refpt_what, SPR_x, sex_ratio_f, is_discard_fleet, local_bh_msy_newton_steps,
                            newton_steps, log_F_start, avg_yrs, rec_yrs,
                            n_yrs, n_pop, n_regions, n_ages, n_seas, n_fish_fleets, seasdur, spawn_seas, t_spawn,
                            natal_region, do_recruits_move, move_timing, expm_nsub, rec_model,
                            Fmort, dmr, fish_sel, ret_sel, natmort, WAA, MatAA, Movement, Mrate, sgl_seas_spawning_movement,
                            stray_rate, rec_seas_prop, rec_region_prop, Rec, R0, h_trans) {

  "[<-" <- RTMB::ADoverload("[<-") # so model quantities can be written into the plain containers below

  # the option by name, so the branches below read as Get_Reference_Points' options
  what <- c("SPR", "MSY", "independent_SPR", "independent_MSY", "global_SPR", "global_MSY", "local_MSY")[refpt_what]

  # how many populations have each region as their natal region, which shares out the strays arriving there
  n_pop_in_region <- rep(0, n_regions)
  for(p in 1:n_pop) n_pop_in_region[natal_region[p]] <- n_pop_in_region[natal_region[p]] + 1

  # females' demographics (sex 1) averaged over the terminal years. every dim is kept, so the year and sex dims
  # of length one are dropped below, in the shape each solver reads
  fish_sel_avg <- avg_over_years(fish_sel[,,,,,1,,drop = FALSE], 3, avg_yrs) # [p, r, 1, seas, a, 1, f]
  ret_sel_avg <- avg_over_years(ret_sel[,,,,,1,,drop = FALSE], 3, avg_yrs)
  natmort_avg <- avg_over_years(natmort[,,,,,1,drop = FALSE], 3, avg_yrs) # [p, r, 1, seas, a, 1]
  WAA_avg <- avg_over_years(WAA[,,,,,1,drop = FALSE], 3, avg_yrs)
  MatAA_avg <- avg_over_years(MatAA[,,,,,1,drop = FALSE], 3, avg_yrs)
  Movement_avg <- avg_over_years(Movement[,,,,,,1,drop = FALSE], 4, avg_yrs) # [p, r, r, 1, seas, a, 1]
  Mrate_avg <- if(move_timing == 2) avg_over_years(Mrate[,,,,,,1,drop = FALSE], 4, avg_yrs) else Movement_avg * 0
  sgl_seas_avg <- avg_over_years(sgl_seas_spawning_movement[,,,,,1,drop = FALSE], 4, avg_yrs) # [p, r, r, 1, a, 1]
  stray_rate_avg <- avg_over_years(stray_rate, 2, avg_yrs) # [p, 1]
  dim(stray_rate_avg) <- n_pop

  # the terminal year's F and discard mortality by region, season and fleet, which split a reference F the same way
  terminal_F <- Fmort[,n_yrs,,,drop = FALSE] # [r, 1, seas, f]
  dim(terminal_F) <- c(n_regions, n_seas, n_fish_fleets)
  dmr_terminal <- dmr[,n_yrs,,,drop = FALSE]
  dim(dmr_terminal) <- c(n_regions, n_seas, n_fish_fleets)

  # mean recruitment by population and region, which turns spawning biomass per recruit into a biomass target
  mean_rec <- avg_over_years(Rec, 3, rec_yrs) # [p, r, 1]
  dim(mean_rec) <- c(n_pop, n_regions)

  # one region's arrays in the shape the single region solvers read: no region, year or sex dims
  single_region_data <- function(r) {
    refpt_data <- list(rec_model = rec_model, t_spawn = t_spawn, n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
                       n_pop = n_pop, n_ages = n_ages, natal_region = natal_region, n_pop_in_region = n_pop_in_region,
                       SPR_x = SPR_x, is_discard_fleet = is_discard_fleet)
    refpt_data$F_fract_flt <- terminal_F[r,,] / sum(terminal_F[r,,]) # each season and fleet's share of the region's terminal F
    dim(refpt_data$F_fract_flt) <- c(n_seas, n_fish_fleets)
    refpt_data$dmr <- dmr_terminal[r,,]
    dim(refpt_data$dmr) <- c(n_seas, n_fish_fleets)
    refpt_data$fish_sel <- fish_sel_avg[,r,1,,,1,,drop = FALSE]
    dim(refpt_data$fish_sel) <- c(n_pop, n_seas, n_ages, n_fish_fleets)
    refpt_data$ret_sel <- ret_sel_avg[,r,1,,,1,,drop = FALSE]
    dim(refpt_data$ret_sel) <- c(n_pop, n_seas, n_ages, n_fish_fleets)
    refpt_data$natmort <- natmort_avg[,r,1,,,1,drop = FALSE]
    dim(refpt_data$natmort) <- c(n_pop, n_seas, n_ages)
    refpt_data$WAA <- WAA_avg[,r,1,,,1,drop = FALSE]
    dim(refpt_data$WAA) <- c(n_pop, n_seas, n_ages)
    refpt_data$MatAA <- MatAA_avg[,r,1,,,1,drop = FALSE]
    dim(refpt_data$MatAA) <- c(n_pop, n_seas, n_ages)
    refpt_data$sex_ratio_f <- sex_ratio_f[,r]
    dim(refpt_data$sex_ratio_f) <- n_pop
    refpt_data$stray_rate <- stray_rate_avg
    refpt_data$rec_seas_prop <- rec_seas_prop
    dim(refpt_data$rec_seas_prop) <- c(n_pop, n_seas)
    refpt_data$h <- h_trans[,r]
    dim(refpt_data$h) <- n_pop
    refpt_data$R0 <- if(refpt_type == 1) R0 else R0 * rec_region_prop[,r] # unfished recruitment, this region's share when regions are solved apart
    dim(refpt_data$R0) <- n_pop
    return(refpt_data)
  }

  # every region's arrays in the shape the spatial solvers read: no year or sex dims, and under drop_pop no
  # population dim either, for the solvers written for one population
  spatial_data <- function(drop_pop) {
    pop_dim <- if(drop_pop) NULL else n_pop
    refpt_data <- list(rec_model = rec_model, t_spawn = t_spawn, n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
                       n_pop = n_pop, n_ages = n_ages, n_regions = n_regions, natal_region = natal_region, n_pop_in_region = n_pop_in_region,
                       SPR_x = SPR_x, is_discard_fleet = is_discard_fleet, move_timing = move_timing, expm_nsub = expm_nsub,
                       do_recruits_move = do_recruits_move, newton_steps = local_bh_msy_newton_steps)
    refpt_data$F_fract_flt <- terminal_F * 0
    for(r in 1:n_regions) refpt_data$F_fract_flt[r,,] <- terminal_F[r,,] / sum(terminal_F[r,,]) # each region's season and fleet shares
    refpt_data$dmr <- dmr_terminal
    refpt_data$fish_sel <- fish_sel_avg
    dim(refpt_data$fish_sel) <- c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)
    refpt_data$ret_sel <- ret_sel_avg
    dim(refpt_data$ret_sel) <- c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)
    refpt_data$natmort <- natmort_avg
    dim(refpt_data$natmort) <- c(pop_dim, n_regions, n_seas, n_ages)
    refpt_data$WAA <- WAA_avg
    dim(refpt_data$WAA) <- c(pop_dim, n_regions, n_seas, n_ages)
    refpt_data$MatAA <- MatAA_avg
    dim(refpt_data$MatAA) <- c(pop_dim, n_regions, n_seas, n_ages)
    refpt_data$Movement <- Movement_avg
    dim(refpt_data$Movement) <- c(pop_dim, n_regions, n_regions, n_seas, n_ages)
    refpt_data$Mrate <- Mrate_avg
    dim(refpt_data$Mrate) <- c(pop_dim, n_regions, n_regions, n_seas, n_ages)
    refpt_data$sgl_seas_spawning_movement <- sgl_seas_avg
    dim(refpt_data$sgl_seas_spawning_movement) <- c(n_pop, n_regions, n_regions, n_ages)
    refpt_data$sex_ratio_f <- sex_ratio_f
    dim(refpt_data$sex_ratio_f) <- c(pop_dim, n_regions)
    refpt_data$rec_seas_prop <- rec_seas_prop
    dim(refpt_data$rec_seas_prop) <- c(pop_dim, n_seas)
    refpt_data$stray_rate <- stray_rate_avg
    refpt_data$h <- h_trans
    dim(refpt_data$h) <- c(pop_dim, n_regions)
    refpt_data$R0 <- R0
    refpt_data$rec_region_prop <- rec_region_prop
    dim(refpt_data$rec_region_prop) <- c(pop_dim, n_regions)
    return(refpt_data)
  }

  # a population's effective spawning biomass at region r: its own spawning biomass per recruit times its
  # recruitment, plus what each other population strays there, its stray rate shared among the populations
  # whose natal region r is
  effective_biomass <- function(p_home, r, spawn_per_rec, rec) {
    biomass <- spawn_per_rec[p_home] * rec[p_home]
    for(p_stray in 1:n_pop) {
      if(p_stray == p_home) next
      stray_share <- stray_rate_avg[p_stray] / n_pop_in_region[r]
      biomass <- biomass + stray_share * spawn_per_rec[p_stray] * rec[p_stray]
    } # end p_stray loop
    return(biomass)
  }

  f_ref_pt <- rep(0, n_regions)
  b_ref_pt <- virgin_b_ref_pt <- pop_b_ref_pt <- virgin_pop_b_ref_pt <- array(0, dim = c(n_pop, n_regions))

  # one region, or each region on its own: the single region solvers, region by region
  if(what %in% c("SPR", "MSY", "independent_SPR", "independent_MSY")) {

    for(r in 1:n_regions) {

      refpt_data <- single_region_data(r)

      # F at the target SPR with mean recruitment, or F at maximum yield with its equilibrium recruitment
      if(what %in% c("SPR", "independent_SPR")) {
        solved <- solve_refpt_tape(single_region_SPR_quantities, refpt_data, log_F_start, newton_steps)
        f_ref_pt[r] <- solved$F_x
        rec_fished <- rec_unfished <- mean_rec[,r]
      } else {
        solved <- solve_refpt_tape(single_region_Fmsy_quantities, refpt_data, log_F_start, newton_steps)
        f_ref_pt[r] <- solved$Fmsy
        rec_fished <- solved$Req
        rec_unfished <- refpt_data$R0
      }

      # spawning biomass per recruit times recruitment, fished at the reference F and unfished
      b_ref_pt[,r] <- solved$SB * rec_fished
      virgin_b_ref_pt[,r] <- solved$SB0 * rec_unfished

      # each population's effective biomass in this region, with the others' strays
      for(p_home in 1:n_pop) {
        pop_b_ref_pt[p_home,r] <- effective_biomass(p_home, r, solved$SB, rec_fished)
        virgin_pop_b_ref_pt[p_home,r] <- effective_biomass(p_home, r, solved$SB0, rec_unfished)
      } # end p_home loop

    } # end r loop

    # one population has one natal region, so its effective biomass is the total over regions
    if(n_pop == 1) {
      pop_b_ref_pt[] <- 0
      virgin_pop_b_ref_pt[] <- 0
      pop_b_ref_pt[1,1] <- sum(b_ref_pt)
      virgin_pop_b_ref_pt[1,1] <- sum(virgin_b_ref_pt)
    }

  } # end single region solvers

  # one F shared by every region, at the target SPR
  if(what == "global_SPR") {

    refpt_data <- spatial_data(drop_pop = FALSE)
    total_mean_rec <- rep(0, n_pop)
    for(p in 1:n_pop) total_mean_rec[p] <- sum(mean_rec[p,])
    refpt_data$rec_region_prop <- mean_rec / total_mean_rec # where recruits settle, read off the fit's own recruitment
    solved <- solve_refpt_tape(global_SPR_quantities, refpt_data, log_F_start, newton_steps)

    # spawning biomass per recruit by region times each population's total mean recruitment
    f_ref_pt[] <- solved$F_x
    for(r in 1:n_regions) {
      b_ref_pt[,r] <- solved$SB[,r] * total_mean_rec
      virgin_b_ref_pt[,r] <- solved$SB0[,r] * total_mean_rec
    } # end r loop

    # each population's effective biomass at its natal region; one population's is the total over regions
    if(n_pop > 1) {
      for(p_home in 1:n_pop) {
        r_natal <- natal_region[p_home]
        pop_b_ref_pt[p_home,r_natal] <- effective_biomass(p_home, r_natal, solved$SB[,r_natal], total_mean_rec)
        virgin_pop_b_ref_pt[p_home,r_natal] <- effective_biomass(p_home, r_natal, solved$SB0[,r_natal], total_mean_rec)
      } # end p_home loop
    } else {
      pop_b_ref_pt[1,1] <- sum(b_ref_pt)
      virgin_pop_b_ref_pt[1,1] <- sum(virgin_b_ref_pt)
    }

  } # end global SPR

  # one F shared by every region, at maximum yield; one population only
  if(what == "global_MSY") {

    if(n_pop > 1) stop("global_MSY needs one population. Use independent_SPR, independent_MSY, global_SPR or local_MSY.", call. = FALSE)
    refpt_data <- spatial_data(drop_pop = TRUE)
    refpt_data$rec_region_prop <- rec_region_prop[1,]
    refpt_data$sex_ratio_f <- sex_ratio_f
    refpt_data$rec_seas_prop <- rec_seas_prop[1,]
    refpt_data$h <- mean(h_trans[1,])
    refpt_data$R0 <- R0[1]
    solved <- solve_refpt_tape(global_Fmsy_quantities, refpt_data, log_F_start, newton_steps)

    # spawning biomass by region, fished (row 2) and unfished (row 1), times equilibrium and unfished recruitment
    f_ref_pt[] <- solved$Fmsy
    for(r in 1:n_regions) {
      b_ref_pt[1,r] <- sum(solved$SB_age[2,r,]) * solved$Req
      virgin_b_ref_pt[1,r] <- sum(solved$SB_age[1,r,]) * refpt_data$R0
    } # end r loop
    pop_b_ref_pt[1,1] <- sum(b_ref_pt)
    virgin_pop_b_ref_pt[1,1] <- sum(virgin_b_ref_pt)

  } # end global MSY

  # one F per region, at maximum total yield
  if(what == "local_MSY") {

    refpt_data <- spatial_data(drop_pop = n_pop == 1)
    quantities <- if(n_pop == 1) local_Fmsy_sglpop_quantities else local_Fmsy_multipop_quantities
    solved <- solve_refpt_tape(quantities, refpt_data, rep(log_F_start, n_regions), newton_steps)
    f_ref_pt <- solved$Fmsy

    if(n_pop > 1) {

      # spawning biomass per recruit by region times each population's equilibrium and unfished recruitment
      for(r in 1:n_regions) {
        b_ref_pt[,r] <- solved$SB[,r] * solved$Req_o
        virgin_b_ref_pt[,r] <- solved$SB0[,r] * R0
      } # end r loop

      # each population's effective biomass at its natal region
      for(p_home in 1:n_pop) {
        r_natal <- natal_region[p_home]
        pop_b_ref_pt[p_home,r_natal] <- effective_biomass(p_home, r_natal, solved$SB[,r_natal], solved$Req_o)
        virgin_pop_b_ref_pt[p_home,r_natal] <- effective_biomass(p_home, r_natal, solved$SB0[,r_natal], R0)
      } # end p_home loop

    } else {

      # one population: spawning biomass by region summed over where its recruits settle
      for(r in 1:n_regions) {
        b_ref_pt[1,r] <- sum(solved$SB_fished_mat[,r] * solved$Req_o)
        virgin_b_ref_pt[1,r] <- sum(solved$SB_unfished_mat[,r] * R0 * refpt_data$rec_region_prop)
      } # end r loop
      pop_b_ref_pt[1,1] <- sum(b_ref_pt)
      virgin_pop_b_ref_pt[1,1] <- sum(virgin_b_ref_pt)

    }

  } # end local MSY

  return(list(f_ref_pt = f_ref_pt,
              b_ref_pt = b_ref_pt,
              virgin_b_ref_pt = virgin_b_ref_pt,
              pop_b_ref_pt = pop_b_ref_pt,
              virgin_pop_b_ref_pt = virgin_pop_b_ref_pt))

} # end get_refpts_tape
