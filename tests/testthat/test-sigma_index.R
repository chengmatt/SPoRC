# Estimated index observation error. An index carries a standard error from its
# own survey design, and an assessment may additionally estimate a component
# covering everything that design does not. These tests pin the three ways the
# two are combined, the separation of the aggregated and population-specific
# streams, and the two guard rails: a spec that cannot be identified, and a spec
# on a likelihood that ignores standard deviations entirely.

library(SPoRC)
library(testthat)

# A minimal single region, single sex, single season model. Small enough to fit
# in a test, large enough that an index sigma is identified.
build_toy <- function(n_yrs = 20, n_srv = 1, ...) {
  yrs <- seq_len(n_yrs); ages <- 1:8; n_ages <- length(ages)
  d1 <- c(1, 1, n_yrs, 1, n_ages, 1)

  il <- Setup_Mod_Dim(n_pop = 1, years = yrs, ages = ages, lens = NA,
                      n_regions = 1, n_sexes = 1, n_seas = 1,
                      n_fish_fleets = 1, n_srv_fleets = n_srv, verbose = FALSE)

  il <- Setup_Mod_Rec(il, rec_model = "mean_rec", sigmaR_spec = "fix",
                      do_rec_bias_ramp = 0, init_age_strc = 1, ln_global_R0 = log(1e6))

  il <- suppressWarnings(Setup_Mod_Biologicals(
    il, WAA = array(1, dim = d1),
    WAA_fish = array(1, dim = c(d1, 1)), WAA_srv = array(1, dim = c(d1, n_srv)),
    MatAA = array(1, dim = d1), fit_lengths = 0, M_spec = "fix",
    Fixed_natmort = array(0.2, dim = c(1, 1, n_yrs, n_ages, 1))))

  il <- Setup_Mod_Movement(il, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  il <- Setup_Mod_Tagging(il, use_conv_fish_tagging = 0)

  il <- suppressWarnings(Setup_Mod_Catch_and_F(
    il, ObsCatch = array(1e4, dim = c(1, n_yrs, 1, 1)),
    UseCatch = array(1, dim = c(1, n_yrs, 1, 1)),
    sigmaC_spec = "fix", sigmaF_spec = "fix"))

  il <- Setup_Mod_FishIdx_and_Comps(
    il, ObsFishIdx = array(NA, dim = c(1, n_yrs, 1, 1)),
    ObsFishIdx_SE = array(NA, dim = c(1, n_yrs, 1, 1)),
    UseFishIdx = array(0, dim = c(1, n_yrs, 1, 1)),
    ObsFishAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
    UseFishAgeComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_FishAgeComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    ObsFishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1, 1)),
    UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    fish_idx_type = "none",
    FishAgeComps_LikeType = "none", FishLenComps_LikeType = "none",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "agg_Year_1-terminal_Fleet_1")

  fl <- function(tag) paste0(tag, "_Fleet_", seq_len(n_srv))
  il <- Setup_Mod_SrvIdx_and_Comps(
    il, ObsSrvIdx = array(1e5, dim = c(1, n_yrs, 1, n_srv)),
    ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, n_srv)),
    UseSrvIdx = array(1, dim = c(1, n_yrs, 1, n_srv)),
    ObsSrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, n_srv)),
    UseSrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
    ISS_SrvAgeComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
    ObsSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1, n_srv)),
    UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_srv)),
    ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, n_srv)),
    srv_idx_type = rep("abd", n_srv),
    SrvAgeComps_LikeType = rep("none", n_srv), SrvLenComps_LikeType = rep("none", n_srv),
    SrvAgeComps_Type = fl("agg_Year_1-terminal"), SrvLenComps_Type = fl("agg_Year_1-terminal"),
    ...)

  il <- Setup_Mod_Fishsel_and_Q(
    il, cont_tv_fish_sel = "none_Fleet_1", fish_sel_blocks = "none_Fleet_1",
    fish_sel_model = "logist1_Fleet_1", fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "fix")

  il <- Setup_Mod_Srvsel_and_Q(
    il, cont_tv_srv_sel = fl("none"), srv_sel_blocks = fl("none"),
    srv_sel_model = fl("logist1"), srv_q_blocks = fl("none"),
    srv_fixed_sel_pars_spec = rep("est_all", n_srv),
    srv_q_spec = rep("est_all", n_srv))

  cd <- c(1, n_yrs, 1, 1, 1); sd <- c(1, n_yrs, 1, 1, n_srv)
  Setup_Mod_Weighting(il, Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1,
                      Wt_F = 1, Wt_Tagging = 0,
                      Wt_FishAgeComps = array(0, dim = cd), Wt_FishLenComps = array(0, dim = cd),
                      Wt_SrvAgeComps = array(0, dim = sd), Wt_SrvLenComps = array(0, dim = sd))
}

sd_of <- function(il) {
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  obj$rep$SrvIdx_SD
}

test_that("combine_idx_sd implements the three forms it documents", {
  se <- c(0.2, 0.3); extra <- 0.4
  expect_equal(combine_idx_sd(se, extra, 0), se)
  expect_equal(combine_idx_sd(se, extra, 1), se + extra)
  expect_equal(combine_idx_sd(se, extra, 2), sqrt(se^2 + extra^2))
  expect_equal(combine_idx_sd(se, extra, 3), rep(extra, length(se)))
})

test_that("build_idx_sd applies the estimated component fleet by fleet", {
  se <- array(0.2, dim = c(1, 3, 1, 2))
  se[,,,2] <- 0.5
  out <- build_idx_sd(se, log(c(0.1, 0.3)), 1)
  expect_equal(as.numeric(out[,,,1]), rep(0.3, 3))
  expect_equal(as.numeric(out[,,,2]), rep(0.8, 3))
  # and a fixed form is a passthrough, whatever the parameter holds
  expect_equal(build_idx_sd(se, log(c(9, 9)), 0), se)
})

test_that("the default leaves the reported standard errors untouched", {
  il <- build_toy()
  expect_equal(il$data$sigmaSrvIdx_form, 0)
  expect_true(all(is.na(as.integer(il$map$ln_sigmaSrvIdx))))
  expect_equal(as.numeric(sd_of(il)), as.numeric(il$data$ObsSrvIdx_SE))
})

test_that("each estimated form reaches the likelihood with the right total", {
  extra <- 0.35
  for(spec in c("est_additive", "est_quadrature", "est_replace")) {
    il <- build_toy(sigmaSrvIdx_spec = spec, ln_sigmaSrvIdx = log(extra))
    form <- match(spec, c("est_additive", "est_quadrature", "est_replace"))
    expect_equal(il$data$sigmaSrvIdx_form, form, info = spec)
    expect_equal(as.numeric(sd_of(il)),
                 as.numeric(combine_idx_sd(il$data$ObsSrvIdx_SE, extra, form)),
                 info = spec)
  }
})

test_that("the fleet map shares and holds fleets as asked", {
  il <- build_toy(n_srv = 3, sigmaSrvIdx_spec = "est_additive",
                  sigmaSrvIdx_map = c(NA, 1, 1))
  m <- as.integer(il$map$ln_sigmaSrvIdx)
  expect_true(is.na(m[1]))          # held at its starting value
  expect_equal(m[2], m[3])          # and these two share one parameter
  expect_equal(length(unique(stats::na.omit(m))), 1)
})

test_that("the aggregated and population streams carry separate parameters", {
  il <- build_toy(sigmaSrvIdx_spec = "est_additive", sigmaSrvIdx_pop_spec = "est_replace")
  expect_equal(il$data$sigmaSrvIdx_form, 1)
  expect_equal(il$data$sigmaSrvIdx_pop_form, 3)
  expect_true("ln_sigmaSrvIdx" %in% names(il$par))
  expect_true("ln_sigmaSrvIdx_pop" %in% names(il$par))
  # they are distinct entries, so one can be held while the other is estimated
  il2 <- build_toy(sigmaSrvIdx_spec = "est_additive", sigmaSrvIdx_pop_spec = "fix")
  expect_false(all(is.na(as.integer(il2$map$ln_sigmaSrvIdx))))
  expect_true(all(is.na(as.integer(il2$map$ln_sigmaSrvIdx_pop))))
})

test_that("an estimated sigma is refused on a multivariate normal index", {
  n_yrs <- 20
  cov1 <- list(diag(n_yrs))
  expect_error(
    build_toy(SrvIdx_LikeType = "mvn", SrvIdx_Cov = cov1,
              sigmaSrvIdx_spec = "est_additive"),
    "multivariate normal", ignore.case = TRUE)
})

test_that("an unrecognised spec is rejected by name", {
  expect_error(build_toy(sigmaSrvIdx_spec = "estimate_it"),
               "sigmaSrvIdx_spec is .estimate_it., which is not recognized")
})

test_that("a catch sigma with one observation per parameter is refused", {
  # ln_sigmaC has exactly the dimensions of ObsCatch, so any spec leaving the
  # year dimension free puts one variance on one observation. That likelihood is
  # unbounded rather than merely poorly determined, so it must not be reachable.
  dims <- c(region = 1, year = 20, season = 1, fleet = 1)
  ab <- c(r = "region", y = "year", seas = "season", f = "fleet")
  use <- array(1, dim = unname(dims))

  expect_error(build_shared_spec_map(dims, "est_all", ab, use = use, what = "sigmaC"),
               "cannot be estimated")
  # sharing over every dimension is one parameter over every observation, fine
  expect_silent(build_shared_spec_map(dims, "est_shared_r_y_seas_f", ab, use = use, what = "sigmaC"))
  # and without a use array the check is skipped, so old callers are unaffected
  expect_silent(build_shared_spec_map(dims, "est_all", ab))
})

test_that("a thinly informed catch sigma warns without erroring", {
  dims <- c(region = 1, year = 4, season = 1, fleet = 2)
  ab <- c(r = "region", y = "year", seas = "season", f = "fleet")
  use <- array(1, dim = unname(dims))
  # sharing over years leaves one parameter per fleet, four observations each
  expect_warning(build_shared_spec_map(dims, "est_shared_y", ab, use = use, what = "sigmaC"),
                 "poorly determined")
})

test_that("estimating an index sigma under a likelihood weight warns", {
  expect_warning(build_toy(sigmaSrvIdx_spec = "est_additive"),
                 NA) # unit weight is fine
  il <- build_toy(sigmaSrvIdx_spec = "est_additive")
  cd <- c(1, 20, 1, 1, 1)
  expect_warning(
    Setup_Mod_Weighting(il, Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 5, Wt_Rec = 1,
                        Wt_F = 1, Wt_Tagging = 0,
                        Wt_FishAgeComps = array(0, dim = cd), Wt_FishLenComps = array(0, dim = cd),
                        Wt_SrvAgeComps = array(0, dim = cd), Wt_SrvLenComps = array(0, dim = cd)),
    "confounded")
})

test_that("an input list built before the feature existed still runs", {
  il <- build_toy()
  # strip every trace of the feature, as an older saved list would be
  il$data$sigmaSrvIdx_form <- NULL; il$data$sigmaFishIdx_form <- NULL
  il$data$sigmaSrvIdx_pop_form <- NULL; il$data$sigmaFishIdx_pop_form <- NULL
  il$par$ln_sigmaSrvIdx <- NULL; il$par$ln_sigmaFishIdx <- NULL
  il$par$ln_sigmaSrvIdx_pop <- NULL; il$par$ln_sigmaFishIdx_pop <- NULL
  il$map$ln_sigmaSrvIdx <- NULL; il$map$ln_sigmaFishIdx <- NULL
  il$map$ln_sigmaSrvIdx_pop <- NULL; il$map$ln_sigmaFishIdx_pop <- NULL
  il$data$ObsSrvIdx_pop_SE <- NULL; il$data$ObsFishIdx_pop_SE <- NULL

  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$rep$jnLL))
  expect_equal(as.numeric(obj$rep$SrvIdx_SD), as.numeric(il$data$ObsSrvIdx_SE))
})
