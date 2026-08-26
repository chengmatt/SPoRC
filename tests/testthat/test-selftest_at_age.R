library(SPoRC)
library(testthat)

# Parameter recovery for the age-disaggregated streams. The operating model
# draws catch at age through the same lognormal error the estimation model fits,
# so refitting the simulated data should return the truth without median bias.
#
# This exercises the whole simulation path, not just the likelihood: the draw in
# sim_observations.R, the arguments and stores in setup_sim_fleets.R, the
# containers, and the self test's reset blocks. A draw with no parameter behind
# it passes every structural test and fails here.

at_age_cfg <- list(n_yrs = 40, n_ages = 8, n_sims = 20, sigmaCAA = 0.25, waa = 1, mat = 1)

# operating model fitting catch at age on a single fishery fleet
at_age_om <- function(seed = 321) {

  n_yrs <- at_age_cfg$n_yrs; n_ages <- at_age_cfg$n_ages

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages,
                            n_lens = NULL, n_sexes = 1, n_fish_fleets = 1,
                            n_srv_fleets = 1, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)

  curve <- function(slope, infl) array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
                                       dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))

  # every age fit, with one standard deviation shared across ages
  use_aa <- array(1, dim = c(1, n_yrs, 1, n_ages, 1))

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, curve(3, 2)),
    ret_sel_input = replicate(1, curve(3, 2)),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(rep(seq(0.05, 0.25, length.out = n_yrs), each = 1), dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 1, 1)),
    ln_sigmaCAA = array(log(at_age_cfg$sigmaCAA), dim = c(n_ages, 1)),
    UseCatchAA = use_aa,
    use_catch_aa = 1
  )

  # the survey draws both an aggregated index and an index at age, so each test
  # can take whichever stream it fits
  sim_list <- Setup_Sim_Survey(sim_list = sim_list, srv_sel_input = replicate(1, curve(1, 3)),
                               ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)),
                               ISS_SrvAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 1, 1)),
                               ln_sigmaSrvIdxAA = array(log(0.2), dim = c(n_ages, 1)),
                               UseSrvIdxAA = use_aa, use_srv_idx_aa = 1)

  biol <- function(val) array(rep(val, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, biol(at_age_cfg$waa)),
    WAA_fish_input = replicate(1, array(at_age_cfg$waa, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(at_age_cfg$waa, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol(at_age_cfg$mat))
  ))

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  sim_list <- Setup_Sim_Rec(sim_list = sim_list,
                            R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
                            ln_sigmaR = array(log(0.3), dim = c(2, 1, 1)),
                            recruitment_opt = "mean_rec", init_age_strc = 1)

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}


test_that("the operating model draws catch at age with the error it was given", {

  om <- at_age_om()

  # log residual of the draw against the truth recovers the standard deviation
  true_caa <- om$TrueCatchAA[1, , 1, , 1, 1]
  obs_caa <- om$ObsCatchAA[1, , 1, , 1, 1]
  fit <- true_caa > 0 & obs_caa > 0

  expect_gt(sum(fit), 100)  # the draw actually happened
  expect_equal(stats::sd(log(obs_caa[fit] / true_caa[fit])), at_age_cfg$sigmaCAA, tolerance = 0.12)
  expect_lt(abs(stats::median(log(obs_caa[fit] / true_caa[fit]))), 0.05)
})

# Estimation model over one at-age stream. Catchability and R0 are pinned at the
# operating model's values so the comparison is about the dynamics rather than
# the abundance/catchability ridge, and init devs are held because the operating
# model starts from deterministic equilibrium.
at_age_em <- function(om, stream = "catch", extra_disc = NULL, extra_fidx = NULL) {

  n_yrs <- at_age_cfg$n_yrs; n_ages <- at_age_cfg$n_ages
  yrs <- seq_len(n_yrs); ages <- seq_len(n_ages); d1 <- c(1, 1, n_yrs, 1, n_ages, 1)
  zero4 <- array(0, dim = c(1, n_yrs, 1, 1))
  aa_dim <- c(1, n_yrs, 1, n_ages, 1)

  obs_aa <- array(om$TrueCatchAA[1, , 1, , 1, 1], dim = aa_dim)
  use_aa <- array(as.numeric(obs_aa > 0), dim = aa_dim)

  il <- Setup_Mod_Dim(n_pop = 1, years = yrs, ages = ages, lens = NA, n_regions = 1,
                      n_sexes = 1, n_seas = 1, n_fish_fleets = 1, n_srv_fleets = 1, verbose = FALSE)
  il <- Setup_Mod_Rec(il, rec_model = "mean_rec", sigmaR_spec = "fix",
                      do_rec_bias_ramp = 0, init_age_strc = 1, ln_global_R0 = log(5))
  il$map$ln_global_R0 <- factor(NA)
  il$map$ln_InitDevs <- factor(rep(NA, length(il$par$ln_InitDevs)))

  il <- suppressWarnings(Setup_Mod_Biologicals(
    il, WAA = array(at_age_cfg$waa, dim = d1), WAA_fish = array(at_age_cfg$waa, dim = c(d1, 1)),
    WAA_srv = array(at_age_cfg$waa, dim = c(d1, 1)), MatAA = array(at_age_cfg$mat, dim = d1),
    fit_lengths = 0, M_spec = "fix", Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))))
  il <- Setup_Mod_Movement(il, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  il <- Setup_Mod_Tagging(il, use_conv_fish_tagging = 0)

  # catch at age is always fit: without it fishing mortality has nothing to
  # constrain it and drifts whatever the survey says
  catch_args <- list(input_list = il, ObsCatch = zero4, UseCatch = zero4,
                     catch_units = array("abd", dim = 1), Use_F_pen = 1,
                     sigmaC_spec = "fix", sigmaF_spec = "fix",
                     ObsCatchAA = obs_aa, UseCatchAA = use_aa,
                     sigmaCAA_key = array(1L, dim = c(n_ages, 1)), sigmaCAA_spec = "fix",
                     ln_sigmaCAA = array(log(at_age_cfg$sigmaCAA), dim = c(n_ages, 1)))
  if(stream == "discard") {
    # the operating model discards nothing, so this checks the discard stream
    # reaches the likelihood and leaves the fit intact
    catch_args <- c(catch_args, list(
      ObsDiscardAA = array(1e-6, dim = aa_dim), UseDiscardAA = array(0, dim = aa_dim),
      sigmaDAA_spec = "fix"))
  }
  if(!is.null(extra_disc)) {
    catch_args <- c(catch_args, list(
      ObsDiscardAA = extra_disc, UseDiscardAA = array(1, dim = aa_dim),
      sigmaDAA_key = array(1L, dim = c(n_ages, 1)), sigmaDAA_spec = "fix"))
  }
  il <- suppressWarnings(do.call(Setup_Mod_Catch_and_F, catch_args))

  il <- Setup_Mod_FishIdx_and_Comps(
    il, ObsFishIdx = array(NA, dim = c(1, n_yrs, 1, 1)), ObsFishIdx_SE = array(NA, dim = c(1, n_yrs, 1, 1)),
    UseFishIdx = zero4, ObsFishAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
    UseFishAgeComps = zero4, ISS_FishAgeComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    ObsFishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1, 1)), UseFishLenComps = zero4,
    ISS_FishLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    fish_idx_type = "none", FishAgeComps_LikeType = "none", FishLenComps_LikeType = "none",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1", FishLenComps_Type = "agg_Year_1-terminal_Fleet_1",
    ObsFishIdxAA = extra_fidx,
    UseFishIdxAA = if(is.null(extra_fidx)) NULL else array(1, dim = aa_dim),
    sigmaFishIdxAA_key = if(is.null(extra_fidx)) NULL else array(1L, dim = c(n_ages, 1)),
    sigmaFishIdxAA_spec = "fix")

  srv_args <- list(input_list = il,
    ObsSrvIdx = array(om$TrueSrvIdx[1, , 1, 1, 1], dim = c(1, n_yrs, 1, 1)),
    ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)), UseSrvIdx = array(1, dim = c(1, n_yrs, 1, 1)),
    ObsSrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1)), UseSrvAgeComps = zero4,
    ISS_SrvAgeComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    ObsSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1, 1)), UseSrvLenComps = zero4,
    ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    srv_idx_type = "abd", SrvAgeComps_LikeType = "none", SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1", SrvLenComps_Type = "agg_Year_1-terminal_Fleet_1")
  if(stream == "srv_index") {
    # survey index at age instead of the aggregate
    srv_aa <- array(om$TrueSrvIdxAA[1, , 1, , 1, 1], dim = aa_dim)
    use_srv_aa <- array(as.numeric(srv_aa > 0), dim = aa_dim)
    srv_args$UseSrvIdx <- zero4
    srv_args$ObsSrvIdxAA <- srv_aa
    srv_args$UseSrvIdxAA <- use_srv_aa
    srv_args$sigmaSrvIdxAA_key <- array(1L, dim = c(n_ages, 1))
    srv_args$sigmaSrvIdxAA_spec <- "fix"
    srv_args$ln_sigmaSrvIdxAA <- array(log(0.2), dim = c(n_ages, 1))
  }
  il <- do.call(Setup_Mod_SrvIdx_and_Comps, srv_args)

  il <- Setup_Mod_Fishsel_and_Q(il, cont_tv_fish_sel = "none_Fleet_1", fish_sel_blocks = "none_Fleet_1",
    fish_sel_model = "logist1_Fleet_1", fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "fix")
  # an index fit age by age carries its age shape in selectivity, which is what
  # the "nonparfree" form is for: one free value per age, no standardization, so
  # the values hold the height of the curve as well as its shape
  aa_srv <- stream == "srv_index"
  il <- Setup_Mod_Srvsel_and_Q(il, cont_tv_srv_sel = "none_Fleet_1", srv_sel_blocks = "none_Fleet_1",
    srv_sel_model = if(aa_srv) "nonparfree_Fleet_1" else "logist1_Fleet_1",
    srv_sel_nonpar_est_bins = if(aa_srv) list(list(as.list(seq_len(n_ages)))) else NULL,
    srv_q_blocks = "none_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "fix", ln_srv_q = array(log(1), dim = c(1, 1, 1)))

  cd <- c(1, n_yrs, 1, 1, 1)
  Setup_Mod_Weighting(il, Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1,
    Wt_Tagging = 0, Wt_FishAgeComps = array(0, dim = cd), Wt_FishLenComps = array(0, dim = cd),
    Wt_SrvAgeComps = array(0, dim = cd), Wt_SrvLenComps = array(0, dim = cd))
}

at_age_fit <- function(il) {
  fit <- fit_model(il$data, il$par, il$map, random = NULL, silent = TRUE)
  for(i in 1:3) fit$optim <- stats::nlminb(fit$optim$par, fit$fn, fit$gr,
                    control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
  fit
}


test_that("selectivity is estimated for a fleet fitting catch at age", {

  # the mapping gates on whether a fleet has data, and a catch-at-age fleet has
  # none in the aggregated array; holding selectivity fixed here took the
  # recovery below from 1.3% error to 25%
  il <- at_age_em(at_age_om())
  free <- length(unique(stats::na.omit(as.integer(il$map$fish_fixed_sel_pars))))
  expect_gt(free, 0)
})


test_that("catch at age recovers the operating model's spawning biomass", {

  om <- at_age_om()
  fit <- at_age_fit(at_age_em(om, "catch"))
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-3)

  truth <- as.numeric(om$SSB[1, 1, , 1])
  est <- as.numeric(fit$report(fit$env$last.par.best)$SSB)
  expect_gt(stats::cor(log(truth), log(est)), 0.95)
  expect_lt(stats::median(abs((est - truth) / truth)), 0.05)
})


test_that("a survey index at age recovers the operating model's spawning biomass", {

  om <- at_age_om()
  fit <- at_age_fit(at_age_em(om, "srv_index"))
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-3)

  truth <- as.numeric(om$SSB[1, 1, , 1])
  est <- as.numeric(fit$report(fit$env$last.par.best)$SSB)
  expect_gt(stats::cor(log(truth), log(est)), 0.95)
})


test_that("a discard-at-age stream declared but unused leaves the fit unchanged", {

  om <- at_age_om()
  base <- at_age_fit(at_age_em(om, "catch"))
  with_disc <- at_age_fit(at_age_em(om, "discard"))

  expect_equal(with_disc$optim$objective, base$optim$objective, tolerance = 1e-8)
  expect_equal(sum(with_disc$report(with_disc$env$last.par.best)$DiscardAA_nLL), 0)
})

test_that("one-step-ahead residuals work on every at-age stream", {

  # OSA on a stream needs a converged fit, so each stream is seeded from the
  # model's own prediction and given the same standard deviation the fit assumes.
  om <- at_age_om()
  il <- at_age_em(om, "catch")
  r0 <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep

  n_yrs <- at_age_cfg$n_yrs; n_ages <- at_age_cfg$n_ages
  aa_dim <- c(1, n_yrs, 1, n_ages, 1)
  # the operating model retains everything, so it produces no discards to fit;
  # the discard stream is covered by the declared-but-unused test above
  # seeded near, not at, the prediction: an exactly zero residual degenerates
  # the numerical integration one-step-ahead residuals rely on
  set.seed(4242)
  fidx <- array(0, dim = aa_dim)
  for(y in 1:n_yrs) {
    for(a in 1:n_ages) {
      fidx[1, y, 1, a, 1] <- sum(r0$FishIAA[, 1, y, 1, a, , 1]) * exp(stats::rnorm(1, 0, 0.2))
    } # end a loop
  } # end y loop

  il2 <- at_age_em(om, "catch", extra_fidx = fidx)
  fit <- at_age_fit(il2)
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-3)

  for(stream_name in c("CatchAA", "FishIdxAA")) {
    # the observations are Gaussian on the log scale, so the Gaussian method is
    # exact here and avoids the numerical integration the generic method uses
    osa <- suppressWarnings(get_osa(model = fit, data = il2$data, index_source = stream_name,
                                    osa_method = "oneStepGaussian"))
    expect_true(nrow(osa$res) > 0, info = stream_name)
    expect_true("age" %in% names(osa$res), info = stream_name)
    expect_true(all(is.finite(osa$res$resid)), info = stream_name)
    expect_silent(invisible(plot_resids(osa)))
  } # end stream_name loop
})
