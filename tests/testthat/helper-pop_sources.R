# A model with the population-specific data sources on: catch, fishery index, fishery comps,
# survey index and survey comps. Nothing else in the suite turns those five on together.

#' Build a model with the population-specific data sources fit
#'
#' Layers the \code{_pop} arrays onto \code{\link{sweep_input}}, which already
#' fits the regional counterparts, so a likelihood that reads the wrong one of
#' the pair moves a number rather than going quiet.
#'
#' @param n_pop,n_regions,n_yrs,n_ages,n_seas,n_sexes,n_fish,n_srv Dimensions.
#' @param like Index likelihood, one of \code{"lognormal"}, \code{"normal"} or
#'   \code{"mvn"}. The regional and the population-specific index blocks split
#'   the likelihoods between them differently, so both routes need covering.
#'
#' @return An input list ready for \code{fit_model}.
#'
#' @keywords internal
pop_sources_input <- function(n_pop = 2,
                              n_regions = 2,
                              n_yrs = 10,
                              n_ages = 6,
                              n_seas = 1,
                              n_sexes = 1,
                              n_fish = 1,
                              n_srv = 1,
                              like = "lognormal") {

  # the population dim leads every population-specific array, so each shape is the
  # regional one with n_pop in front
  agg_fish <- c(n_pop, n_regions, n_yrs, n_seas, n_fish)                     # one value per fleet
  agg_srv <- c(n_pop, n_regions, n_yrs, n_seas, n_srv)
  comp_fish <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish)   # at age
  comp_srv <- c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv)
  iss_fish <- c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish)            # one per sex
  iss_srv <- c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_srv)

  # a multivariate normal reads a fleet's whole series at once, so its covariance is
  # sized to the cells that fleet fits rather than to any one dim
  n_cells <- n_regions * n_yrs * n_seas
  idx_cov <- if(like == "mvn") {
    corr <- outer(rep(0.7, n_cells), rep(0.7, n_cells)); diag(corr) <- 1
    outer(rep(0.2, n_cells), rep(0.2, n_cells)) * corr
  } else NULL

  sweep_input(

    dims = list(
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_ages = n_ages,
      n_seas = n_seas,
      n_sexes = n_sexes,
      n_fish_fleets = n_fish,
      n_srv_fleets = n_srv,
      natal_region = seq_len(n_pop)
    ),

    # more than one population recruits under local density dependence
    rec = list(rec_dd = "local", ln_global_R0 = rep(log(1e6), n_pop)),

    # catch by population, alongside the regional catch the sweep already fits
    catch = list(
      ObsCatch_pop = array(5e3, dim = agg_fish),
      UseCatch_pop = array(1, dim = agg_fish),
      sigmaC_pop_spec = "fix"
    ),

    # fishery index and age compositions by population
    fishidx = list(
      ObsFishIdx_pop = array(5e4, dim = agg_fish),
      ObsFishIdx_pop_SE = array(0.2, dim = agg_fish),
      UseFishIdx_pop = array(1, dim = agg_fish),
      ObsFishAgeComps_pop = array(1 / n_ages, dim = comp_fish),
      UseFishAgeComps_pop = array(1, dim = agg_fish),
      ISS_FishAgeComps_pop = array(100, dim = iss_fish),
      FishAgeComps_pop_LikeType = rep("Multinomial", n_fish),
      FishAgeComps_pop_Type = paste0("spltRspltS_Year_1-terminal_Fleet_", seq_len(n_fish)),
      FishIdx_LikeType = rep(like, n_fish),
      FishIdx_Cov = if(like == "mvn") rep(list(idx_cov), n_fish) else NULL
    ),

    # survey index and age compositions by population
    srvidx = list(
      ObsSrvIdx_pop = array(5e4, dim = agg_srv),
      ObsSrvIdx_pop_SE = array(0.2, dim = agg_srv),
      UseSrvIdx_pop = array(1, dim = agg_srv),
      ObsSrvAgeComps_pop = array(1 / n_ages, dim = comp_srv),
      UseSrvAgeComps_pop = array(1, dim = agg_srv),
      ISS_SrvAgeComps_pop = array(100, dim = iss_srv),
      SrvAgeComps_pop_LikeType = rep("Multinomial", n_srv),
      SrvAgeComps_pop_Type = paste0("spltRspltS_Year_1-terminal_Fleet_", seq_len(n_srv)),
      SrvIdx_LikeType = rep(like, n_srv),
      SrvIdx_Cov = if(like == "mvn") rep(list(idx_cov), n_srv) else NULL
    ),

    wt = list(
      Wt_Catch_pop = 1,
      Wt_FishIdx_pop = 1,
      Wt_SrvIdx_pop = 1,
      Wt_FishAgeComps_pop = array(1, dim = iss_fish),
      Wt_SrvAgeComps_pop = array(1, dim = iss_srv)
    )
  )
}

#' The population-specific likelihood terms this model produces
#'
#' @param input_list From \code{\link{pop_sources_input}}.
#'
#' @return Named list of the five population-specific likelihood totals, plus the
#'   two regional index totals, which are filled by the other of the two index
#'   routes and are otherwise unchecked.
#'
#' @keywords internal
pop_sources_nLL <- function(input_list) {

  rep <- fit_model(input_list$data, input_list$par, input_list$map,
                   do_optim = FALSE, silent = TRUE)$rep

  list(
    Catch_pop = sum(rep$Catch_pop_nLL),
    FishIdx_pop = sum(rep$FishIdx_pop_nLL),
    FishAgeComps_pop = sum(rep$FishAgeComps_pop_nLL),
    SrvIdx_pop = sum(rep$SrvIdx_pop_nLL),
    SrvAgeComps_pop = sum(rep$SrvAgeComps_pop_nLL),
    FishIdx = sum(rep$FishIdx_nLL),
    SrvIdx = sum(rep$SrvIdx_nLL)
  )
}
