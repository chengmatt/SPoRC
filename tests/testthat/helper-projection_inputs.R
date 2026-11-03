# Post hoc projection inputs built from a data list and report, the archived builder kept as a test helper so the
# internal projection can be checked against Do_Population_Projection on any packaged model. Not package code.

#' Build the inputs of a projection from a fitted model
#'
#' Reads the terminal state off a report, averages the demographics over the last
#' \code{n_avg_yrs} years, repeats them over the projection years, and places the
#' reference points as every projection year's F and biomass targets. Projection
#' year one replays the terminal year, so its spawning biomass should reproduce the
#' report's.
#'
#' @param data Data list from the fit.
#' @param rep Report list from the fit.
#' @param refpts Output of \code{\link{Get_Reference_Points}}.
#' @param n_proj_yrs Integer. Projection years, counting the replayed terminal year.
#' @param n_avg_yrs Integer. Terminal years the demographics are averaged over. Default 1.
#' @param calc_rec_st_yr Integer. First year of the mean recruitment window. Default 1.
#' @param rec_age Integer. Recruitment lag in years. Default 1.
#' @param t_spawn Numeric. Fraction of the spawning season elapsed before spawning.
#'   \code{NULL} (default) reads it from the data list.
#' @param recruitment_opt Recruitment scenario, as in
#'   \code{\link{Do_Population_Projection}}. \code{"bh_rec"} and \code{"ricker_rec"}
#'   also build \code{srr_opt} from the fitted curve.
#'
#' @return Named list of arguments for \code{\link{Do_Population_Projection}}. The
#'   harvest rule, \code{fmort_opt} and any catch targets are added by the caller.
#'
#' @family Reference Points and Projections
#'
#' @examples
#' \dontrun{
#' data("dusky_rtmb_model")
#' rp <- Get_Reference_Points(data = dusky_rtmb_model$data, rep = dusky_rtmb_model$rep,
#'                            SPR_x = 0.4, type = "single_region", what = "SPR")
#' inputs <- build_projection_inputs(data = dusky_rtmb_model$data, rep = dusky_rtmb_model$rep,
#'                                   refpts = rp, n_proj_yrs = 10)
#' out <- do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "Input")))
#' }
build_projection_inputs <- function(data,
                                    rep,
                                    refpts,
                                    n_proj_yrs,
                                    n_avg_yrs = 1,
                                    calc_rec_st_yr = 1,
                                    rec_age = 1,
                                    t_spawn = NULL,
                                    recruitment_opt = "mean_rec") {

  # dimensions
  n_pop <- data$n_pop
  n_regions <- data$n_regions
  n_ages <- length(data$ages)
  n_sexes <- data$n_sexes
  n_seas <- data$n_seas
  n_fish_fleets <- data$n_fish_fleets
  n_yrs <- length(data$years)

  if(n_proj_yrs < 1) stop("n_proj_yrs must be at least 1, the replayed terminal year.", call. = FALSE)
  if(n_avg_yrs < 1 || n_avg_yrs > n_yrs) stop("n_avg_yrs must be between 1 and the number of assessment years.", call. = FALSE)
  avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs # terminal years the demographics are averaged over

  if(is.null(t_spawn)) t_spawn <- data$t_spawn
  if(is.null(t_spawn)) stop("t_spawn is not in the data list, so pass it directly.", call. = FALSE)

  # older reports have no season dim on natural mortality, and older data lists no fishery weight. weight
  # at age comes from the growth module when the model derived it, since the data array is a placeholder then
  natmort_fit <- expand_natmort_seasons(rep$natmort, n_seas)
  WAA_fit <- if(is.null(rep$WAA)) data$WAA else rep$WAA
  WAA_fish_fit <- if(is.null(rep$WAA_fish)) data$WAA_fish else rep$WAA_fish
  if(is.null(WAA_fish_fit)) {
    WAA_fish_fit <- array(0, dim = c(dim(WAA_fit), n_fish_fleets))
    for(f in 1:n_fish_fleets) WAA_fish_fit[,,,,,,f] <- WAA_fit # spawning weight stands in for every fleet
  }

  if(length(dim(rep$fish_sel)) != 7) stop("rep$fish_sel should be population by region by year by season by age by sex by fleet. Older reports lack the population and season dims; refit or expand them first.", call. = FALSE)

  # terminal state
  terminal_NAA <- array(rep$NAA[,,n_yrs,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)) # fished numbers at age
  terminal_NAA0 <- array(rep$NAA0[,,n_yrs,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)) # unfished numbers at age
  terminal_F <- array(rep$Fmort[,n_yrs,,], dim = c(n_regions, n_seas, n_fish_fleets)) # terminal F by fleet
  rec_yrs <- calc_rec_st_yr:(n_yrs - rec_age) # years the mean recruitment is taken over
  recruitment <- array(rep$Rec[,,rec_yrs], dim = c(n_pop, n_regions, length(rec_yrs)))

  # terminal averages of the demographics
  WAA_avg <- apply(WAA_fit[,,avg_yrs,,,,drop = FALSE], c(1,2,4,5,6), mean) # spawning weight
  WAA_fish_avg <- apply(WAA_fish_fit[,,avg_yrs,,,,,drop = FALSE], c(1,2,4,5,6,7), mean) # fishery weight
  MatAA_avg <- apply(data$MatAA[,,avg_yrs,,,,drop = FALSE], c(1,2,4,5,6), mean) # maturity
  natmort_avg <- apply(natmort_fit[,,avg_yrs,,,,drop = FALSE], c(1,2,4,5,6), mean) # natural mortality
  fish_sel_avg <- apply(rep$fish_sel[,,avg_yrs,,,,,drop = FALSE], c(1,2,4,5,6,7), mean) # fishery selectivity
  Movement_avg <- apply(rep$Movement[,,,avg_yrs,,,,drop = FALSE], c(1,2,3,5,6,7), mean) # movement fractions
  sexratio_avg <- apply(rep$sexratio[,,avg_yrs,,drop = FALSE], c(1,2,4), mean) # recruitment sex ratio
  stray_rate_avg <- apply(rep$stray_rate[,avg_yrs,drop = FALSE], 1, mean) # stray rate by population

  # the generator is averaged on the rate scale under continuous movement, since mean(expm(Q)) != expm(mean(Q))
  move_timing <- if(is.null(data$move_timing)) 0 else data$move_timing
  expm_nsub <- if(is.null(data$move_expm_nsub)) 0 else data$move_expm_nsub
  Mrate_avg <- if(move_timing == 2) apply(rep$Mrate[,,,avg_yrs,,,,drop = FALSE], c(1,2,3,5,6,7), mean) else NULL

  # repeat each average over the projection years
  WAA <- MatAA <- natmort <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  WAA_fish <- fish_sel <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  Mrate <- if(move_timing == 2) Movement else NULL
  sgl_seas_spawning_movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_ages, n_sexes))
  sexratio <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes))
  stray_rate <- array(0, dim = c(n_pop, n_proj_yrs))

  for(y in 1:n_proj_yrs) {
    WAA[,,y,,,] <- WAA_avg
    WAA_fish[,,y,,,,] <- WAA_fish_avg
    MatAA[,,y,,,] <- MatAA_avg
    natmort[,,y,,,] <- natmort_avg
    fish_sel[,,y,,,,] <- fish_sel_avg
    Movement[,,,y,,,] <- Movement_avg
    if(move_timing == 2) Mrate[,,,y,,,] <- Mrate_avg
    sgl_seas_spawning_movement[,,,y,,] <- rep$sgl_seas_spawning_movement[,,,n_yrs,,] # terminal year, as the fit used it
    sexratio[,,y,] <- sexratio_avg
    stray_rate[,y] <- stray_rate_avg
  } # end y loop

  # reference points as every year's targets
  f_ref_pt <- array(refpts$f_ref_pt, dim = c(n_regions, n_proj_yrs))
  b_ref_pt <- array(refpts$b_ref_pt, dim = c(n_pop, n_regions, n_proj_yrs))

  # stock-recruit inputs, taken at the year the fit scaled unfished spawning biomass from
  srr_opt <- NULL
  if(recruitment_opt %in% c("bh_rec", "ricker_rec")) {
    sr_yr <- if(is.null(data$SR_ref_yr)) 1 else data$SR_ref_yr
    srr_opt <- list(
      do_recruits_move = data$do_recruits_move,
      rec_dd = data$rec_dd,
      rec_lag = data$rec_lag,
      R0 = rep$R0,
      rec_region_prop = rep$rec_region_prop,
      h = rep$h_trans,
      SSB = rep$SSB,
      WAA = array(WAA_fit[,,sr_yr,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)), # female spawning weight
      MatAA = array(data$MatAA[,,sr_yr,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)), # female maturity
      natmort = array(natmort_fit[,,sr_yr,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)), # female natural mortality
      sgl_seas_spawning_movement = array(rep$sgl_seas_spawning_movement[,,,sr_yr,,1], dim = c(n_pop, n_regions, n_regions, n_ages)),
      Movement = array(rep$Movement[,,,sr_yr,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
      Mrate = if(move_timing == 2) array(rep$Mrate[,,,sr_yr,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)) else NULL,
      stray_rate = array(rep$stray_rate[,1], dim = n_pop),
      sex_ratio_f = array(if(n_sexes == 1) 0.5 else rep$sexratio[,,1,1], dim = c(n_pop, n_regions)),
      # unfished spawning biomass per recruit is taken at zero f, so the fleet arrays only need the right shape
      fish_sel = array(rep$fish_sel[,,sr_yr,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
      ret_sel = array(rep$ret_sel[,,sr_yr,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
      dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
      init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets))
    )
  } # end stock-recruit inputs

  list(n_proj_yrs = n_proj_yrs,
       n_pop = n_pop,
       n_regions = n_regions,
       n_ages = n_ages,
       n_sexes = n_sexes,
       n_seas = n_seas,
       n_fish_fleets = n_fish_fleets,
       seasdur = data$seasdur,
       spawn_seas = data$spawn_seas,
       natal_region = data$natal_region,
       do_recruits_move = data$do_recruits_move,
       rec_seas_prop = rep$rec_seas_prop,
       t_spawn = t_spawn,
       sexratio = sexratio,
       stray_rate = stray_rate,
       recruitment = recruitment,
       terminal_NAA = terminal_NAA,
       terminal_NAA0 = terminal_NAA0,
       terminal_F = terminal_F,
       natmort = natmort,
       WAA = WAA,
       WAA_fish = WAA_fish,
       MatAA = MatAA,
       fish_sel = fish_sel,
       Movement = Movement,
       Mrate = Mrate,
       sgl_seas_spawning_movement = sgl_seas_spawning_movement,
       move_timing = move_timing,
       expm_nsub = expm_nsub,
       f_ref_pt = f_ref_pt,
       b_ref_pt = b_ref_pt,
       recruitment_opt = recruitment_opt,
       srr_opt = srr_opt)

}
