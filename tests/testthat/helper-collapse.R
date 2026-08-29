# Fixtures for the collapse tests: the same population described at two different
# resolutions, so the finer description must reproduce the coarser one.
#
# The observation streams are switched off here. A finer model carries more
# observations than a coarser one and its likelihood is a different number by
# construction, so what has to agree is the population the two describe and the
# quantities predicted from it, not the joint negative log likelihood.

collapse_cfg <- list(n_yrs = 10, n_ages = 6)

#' A population described at a chosen resolution
#'
#' The catch is split evenly across regions and fleets so the total removals are
#' the same however finely the fishery is described.
#'
#' @param nr,nx,nf,ns Regions, sexes, fishery fleets, seasons.
#' @param catch_tot Total catch, divided evenly across region and fleet cells.
#' @param f_scale Multiplier on the mean fishing mortality. At fixed parameters
#'   every fleet carries its own F, so a model split into k fleets fishes k times
#'   as hard unless each fleet is set to 1/k of the rate.
#'
#' @keywords internal
collapse_input <- function(nr = 1, nx = 1, nf = 1, ns = 1, catch_tot = 1e4, f_scale = 1) {
  NY <- collapse_cfg$n_yrs; NAG <- collapse_cfg$n_ages
  off_f <- array(0, dim = c(nr, NY, ns, nf))
  off_s <- array(0, dim = c(nr, NY, ns, 1))

  il <- sweep_input(
    dims = list(n_regions = nr, n_sexes = nx, n_fish_fleets = nf, n_srv_fleets = 1,
                n_seas = ns, n_yrs = NY, n_ages = NAG),
    catch = list(ObsCatch = array(catch_tot / (nr * nf * ns), dim = c(nr, NY, ns, nf)),
                 UseCatch = array(1, dim = c(nr, NY, ns, nf))),
    fishidx = list(
      ObsFishIdx = array(1e5, dim = c(nr, NY, ns, nf)),
      ObsFishIdx_SE = array(0.2, dim = c(nr, NY, ns, nf)),
      UseFishIdx = off_f, fish_idx_type = rep("none", nf),
      ObsFishAgeComps = array(1 / NAG, dim = c(nr, NY, ns, NAG, nx, nf)),
      UseFishAgeComps = off_f, ISS_FishAgeComps = array(100, dim = c(nr, NY, ns, nx, nf)),
      FishAgeComps_LikeType = rep("none", nf),
      ObsFishLenComps = array(0, dim = c(nr, NY, ns, 1, nx, nf)),
      UseFishLenComps = off_f, ISS_FishLenComps = array(0, dim = c(nr, NY, ns, nx, nf))),
    srvidx = list(
      ObsSrvIdx = array(1e5, dim = c(nr, NY, ns, 1)),
      ObsSrvIdx_SE = array(0.2, dim = c(nr, NY, ns, 1)),
      UseSrvIdx = off_s, srv_idx_type = "none",
      ObsSrvAgeComps = array(1 / NAG, dim = c(nr, NY, ns, NAG, nx, 1)),
      UseSrvAgeComps = off_s, ISS_SrvAgeComps = array(100, dim = c(nr, NY, ns, nx, 1)),
      SrvAgeComps_LikeType = "none",
      ObsSrvLenComps = array(0, dim = c(nr, NY, ns, 1, nx, 1)),
      UseSrvLenComps = off_s, ISS_SrvLenComps = array(0, dim = c(nr, NY, ns, nx, 1))),
    wt = list(Wt_FishAgeComps = array(0, dim = c(nr, NY, ns, nx, nf)),
              Wt_FishLenComps = array(0, dim = c(nr, NY, ns, nx, nf)),
              Wt_SrvAgeComps = array(0, dim = c(nr, NY, ns, nx, 1)),
              Wt_SrvLenComps = array(0, dim = c(nr, NY, ns, nx, 1))))

  if(f_scale != 1) il$par$ln_F_mean <- il$par$ln_F_mean + log(f_scale)
  il
}

#' Report from a collapse fixture, evaluated rather than fitted
#'
#' @keywords internal
collapse_rep <- function(...) {
  il <- collapse_input(...)
  fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep
}

#' Assert two reports describe the same population
#'
#' Every reported array is laid out population by region by year, so summing over
#' everything but the year margin gives the total the two resolutions must agree
#' on.
#'
#' @param coarse,fine Reports from \code{collapse_rep}.
#' @param label Name of the relation, used in failure messages.
#' @param what Reported quantities to compare.
#' @param tolerance Relative tolerance. The default is loose enough for the
#'   reordered floating point sums a different array shape produces and far
#'   tighter than any real discrepancy in the dynamics.
#'
#' @keywords internal
expect_collapses <- function(coarse, fine, label,
                             what = c("NAA", "SSB", "Total_Biom", "Rec", "CAA"),
                             tolerance = 1e-10) {
  for(nm in what) {
    if(is.null(coarse[[nm]]) || is.null(fine[[nm]])) next
    a <- apply(coarse[[nm]], 3, sum)
    b <- apply(fine[[nm]], 3, sum)
    testthat::expect_equal(as.numeric(b), as.numeric(a), tolerance = tolerance,
                           label = sprintf("%s: %s", label, nm))
  }
}


# ---------------------------------------------------------------------------
# A small fitted model, for the tests that need parameters, a mapping and an
# sdreport rather than just a report. Only one fitted model ships with the
# package and it is single-region, single-sex, so anything checking the spatial
# or sexed paths has to fit its own.
# ---------------------------------------------------------------------------

#' Fit a small model at chosen dimensions
#'
#' @param nr,nx,nf,np Regions, sexes, fishery fleets, populations.
#' @param n_yrs,n_ages Time series and age range. Small: these are fitted, and
#'   what is being tested is structure rather than estimation quality.
#'
#' @return List with the input list, the fit, and its sdreport.
#'
#' @keywords internal
fitted_small_model <- local({
  cache <- list()
  function(nr = 1, nx = 1, nf = 1, n_yrs = 15, n_ages = 6, np = 1) {
    key <- paste(nr, nx, nf, n_yrs, n_ages, np, sep = "-")
    if(!is.null(cache[[key]])) return(cache[[key]])

    il <- sweep_input(dims = list(n_regions = nr, n_sexes = nx, n_fish_fleets = nf,
                                  n_srv_fleets = 1, n_yrs = n_yrs, n_ages = n_ages,
                                  n_pop = np, natal_region = if(np > 1) rep(1, np) else NA),
                      # populations recruit under local density dependence
                      rec = if(np > 1) list(rec_dd = "local") else list())
    fit <- fit_model(il$data, il$par, il$map, do_optim = TRUE, silent = TRUE, newton_loops = 0)
    sdrep <- RTMB::sdreport(fit, getJointPrecision = FALSE)
    cache[[key]] <<- list(il = il, fit = fit, sdrep = sdrep)
    cache[[key]]
  }
})


#' Fit a small model carrying at-age catch observations
#'
#' The at-age streams replace the aggregated catch for a fleet rather than
#' joining it, so switching them on means switching \code{UseCatch} off. The type
#' is sex-split because the default sums over sexes, which a two-sex model of
#' observations cannot do.
#'
#' @inheritParams fitted_small_model
#'
#' @keywords internal
fitted_at_age_model <- local({
  cache <- list()
  function(nr = 1, nx = 1, n_yrs = 12, n_ages = 6) {
    key <- paste(nr, nx, n_yrs, n_ages, sep = "-")
    if(!is.null(cache[[key]])) return(cache[[key]])

    aa <- c(nr, n_yrs, 1, n_ages, nx, 1)
    il <- sweep_input(
      dims = list(n_regions = nr, n_sexes = nx, n_fish_fleets = 1, n_srv_fleets = 1,
                  n_yrs = n_yrs, n_ages = n_ages),
      catch = list(ObsCatchAA = array(100, dim = aa), UseCatchAA = array(1, dim = aa),
                   CatchAA_Type = "spltRspltS",
                   UseCatch = array(0, dim = c(nr, n_yrs, 1, 1))))
    fit <- fit_model(il$data, il$par, il$map, do_optim = TRUE, silent = TRUE, newton_loops = 0)
    sdrep <- RTMB::sdreport(fit, getJointPrecision = FALSE)
    cache[[key]] <<- list(il = il, fit = fit, sdrep = sdrep, n_yrs = n_yrs, n_ages = n_ages)
    cache[[key]]
  }
})
