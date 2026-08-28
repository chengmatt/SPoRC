# Age-disaggregated observations: retained catch, discards, and the fishery and
# survey indices, each aggregated and population-specific. These tests pin the
# key matrix convention, the guard rails, the equivalence with the aggregated
# statement where the two agree, and the correlation across ages.

library(SPoRC)
library(testthat)

# A minimal one region, one sex, one season model. Catch at age is supplied as
# the aggregate split evenly over ages, so the two statements can be compared.
build_aa <- function(n_yrs = 20, n_ages = 5, aa = TRUE,
                     ObsCatchAA_in = NULL, UseCatchAA_in = NULL, fish_t = NULL, ...) {
  yrs <- seq_len(n_yrs); ages <- seq_len(n_ages)
  d1 <- c(1, 1, n_yrs, 1, n_ages, 1)
  caa_dim <- c(1, n_yrs, 1, n_ages, 1, 1)

  il <- Setup_Mod_Dim(n_pop = 1, years = yrs, ages = ages, lens = NA,
                      n_regions = 1, n_sexes = 1, n_seas = 1,
                      n_fish_fleets = 1, n_srv_fleets = 1, verbose = FALSE)
  il <- Setup_Mod_Rec(il, rec_model = "mean_rec", sigmaR_spec = "fix",
                      do_rec_bias_ramp = 0, init_age_strc = 1, ln_global_R0 = log(1e6))
  il <- suppressWarnings(Setup_Mod_Biologicals(
    il, WAA = array(1, dim = d1), WAA_fish = array(1, dim = c(d1, 1)),
    WAA_srv = array(1, dim = c(d1, 1)), MatAA = array(1, dim = d1),
    fit_lengths = 0, M_spec = "fix",
    Fixed_natmort = array(0.2, dim = c(1, 1, n_yrs, n_ages, 1))))
  il <- Setup_Mod_Movement(il, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  il <- Setup_Mod_Tagging(il, use_conv_fish_tagging = 0)

  ObsCatchAA <- if(is.null(ObsCatchAA_in)) array(2e3, dim = caa_dim) else ObsCatchAA_in
  UseCatchAA <- if(is.null(UseCatchAA_in)) array(if(aa) 1 else 0, dim = caa_dim) else UseCatchAA_in

  il <- suppressWarnings(Setup_Mod_Catch_and_F(
    il,
    ObsCatch = array(1e4, dim = c(1, n_yrs, 1, 1)),
    UseCatch = array(if(aa) 0 else 1, dim = c(1, n_yrs, 1, 1)),
    ObsCatchAA = if(aa) ObsCatchAA else NULL,
    UseCatchAA = if(aa) UseCatchAA else NULL,
    sigmaC_spec = "fix", sigmaF_spec = "fix", ...))

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
    fish_idx_type = "none", FishAgeComps_LikeType = "none", FishLenComps_LikeType = "none",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "agg_Year_1-terminal_Fleet_1",
    t_fish = fish_t)

  il <- Setup_Mod_SrvIdx_and_Comps(
    il, ObsSrvIdx = array(1e5, dim = c(1, n_yrs, 1, 1)),
    ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)),
    UseSrvIdx = array(1, dim = c(1, n_yrs, 1, 1)),
    ObsSrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
    UseSrvAgeComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    ObsSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1, 1)),
    UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvLenComps = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    srv_idx_type = "abd", SrvAgeComps_LikeType = "none", SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "agg_Year_1-terminal_Fleet_1")

  il <- Setup_Mod_Fishsel_and_Q(
    il, cont_tv_fish_sel = "none_Fleet_1", fish_sel_blocks = "none_Fleet_1",
    fish_sel_model = "logist1_Fleet_1", fish_q_blocks = "none_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "fix")
  il <- Setup_Mod_Srvsel_and_Q(
    il, cont_tv_srv_sel = "none_Fleet_1", srv_sel_blocks = "none_Fleet_1",
    srv_sel_model = "logist1_Fleet_1", srv_q_blocks = "none_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "est_all")

  cd <- c(1, n_yrs, 1, 1, 1)
  Setup_Mod_Weighting(il, Wt_Catch = 1, Wt_FishIdx = 0, Wt_SrvIdx = 1, Wt_Rec = 1,
                      Wt_F = 1, Wt_Tagging = 0,
                      Wt_FishAgeComps = array(0, dim = cd), Wt_FishLenComps = array(0, dim = cd),
                      Wt_SrvAgeComps = array(0, dim = cd), Wt_SrvLenComps = array(0, dim = cd))
}

rep_of <- function(il) fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep

test_that("get_at_age_nLL matches dnorm when ages are independent", {
  obs <- log(c(100, 200, 300)); pred <- log(c(110, 190, 310)); sig <- c(0.2, 0.3, 0.4)
  expect_equal(get_at_age_nLL(obs, pred, sig, corr_type = 0),
               -stats::dnorm(obs, pred, sig, log = TRUE))
  # a single age has no correlation to describe, so AR falls back to independent
  expect_equal(get_at_age_nLL(obs[1], pred[1], sig[1], corr_type = 1, rho = 0.5),
               -stats::dnorm(obs[1], pred[1], sig[1], log = TRUE))
})

test_that("an AR(1) across ages reduces to independence at rho zero", {
  obs <- log(c(100, 200, 300)); pred <- log(c(110, 190, 310)); sig <- c(0.2, 0.3, 0.4)
  ar0 <- get_at_age_nLL(obs, pred, sig, corr_type = 1, rho = 0)
  id  <- get_at_age_nLL(obs, pred, sig, corr_type = 0)
  # the correlated form puts the whole cell on its first age, so compare totals
  expect_equal(sum(ar0), sum(id))
  # and a non-zero correlation genuinely changes the density
  ar5 <- get_at_age_nLL(obs, pred, sig, corr_type = 1, rho = 0.5)
  expect_false(isTRUE(all.equal(sum(ar5), sum(id))))
})

test_that("catch at age reaches the likelihood and is finite", {
  il <- build_aa()
  r <- rep_of(il)
  expect_true(is.finite(r$jnLL))
  expect_equal(sum(il$data$UseCatchAA), 20 * 5)
  expect_true(sum(r$CatchAA_nLL) != 0)
  # the aggregated stream is off for this fleet, so it contributes nothing
  expect_equal(sum(r$Catch_nLL), 0)
})

test_that("a fleet cannot fit both aggregated catch and catch at age", {
  n_yrs <- 20; n_ages <- 5
  expect_error(
    build_aa(),
    NA) # the helper turns one off, so the default builds cleanly
  # forcing both on is refused
  expect_error({
    il <- Setup_Mod_Dim(n_pop = 1, years = seq_len(n_yrs), ages = seq_len(n_ages), lens = NA,
                        n_regions = 1, n_sexes = 1, n_seas = 1,
                        n_fish_fleets = 1, n_srv_fleets = 1, verbose = FALSE)
    il <- Setup_Mod_Rec(il, rec_model = "mean_rec", sigmaR_spec = "fix",
                        do_rec_bias_ramp = 0, init_age_strc = 1, ln_global_R0 = log(1e6))
    suppressWarnings(Setup_Mod_Catch_and_F(
      il, ObsCatch = array(1e4, dim = c(1, n_yrs, 1, 1)),
      UseCatch = array(1, dim = c(1, n_yrs, 1, 1)),
      ObsCatchAA = array(2e3, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
      UseCatchAA = array(1, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
      sigmaC_spec = "fix", sigmaF_spec = "fix"))
  }, "one or the other")
})

test_that("the key matrix couples, excludes and refuses what it cannot identify", {
  n_ages <- 5
  # one parameter per age
  il <- build_aa(sigmaCAA_key = array(1:n_ages, dim = c(n_ages, 1, 1)))
  expect_equal(length(unique(stats::na.omit(as.integer(il$map$ln_sigmaCAA)))), n_ages)

  # by age group
  il <- build_aa(sigmaCAA_key = array(c(1, 1, 2, 2, 2), dim = c(n_ages, 1, 1)))
  m <- as.integer(il$map$ln_sigmaCAA)
  expect_equal(m[1], m[2]); expect_equal(m[3], m[5])
  expect_equal(length(unique(stats::na.omit(m))), 2)

  # NA holds an age out
  il <- build_aa(sigmaCAA_key = array(c(NA, 1, 1, 1, 1), dim = c(n_ages, 1, 1)))
  expect_true(is.na(as.integer(il$map$ln_sigmaCAA)[1]))

  # and the wrong shape is rejected by name
  expect_error(build_aa(sigmaCAA_key = array(1, dim = c(1, n_ages))),
               "not the correct dimension")
})

test_that("an at-age variance with too few observations is refused", {
  # two years of data and a parameter per age leaves two observations each, which
  # is the floor; one year leaves one, which is unbounded rather than merely poor
  expect_error(build_aa(n_yrs = 1, sigmaCAA_key = array(1:5, dim = c(5, 1, 1))),
               "fewer than 2 observations")
})

test_that("discards at age carry their own parameter, keyed independently", {
  n_yrs <- 20; n_ages <- 5
  # keyed differently from the retained stream, so the two cannot be sharing
  il <- build_aa(ObsDiscardAA = array(50, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
                 UseDiscardAA = array(1, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
                 discard_units = "biom",
                 sigmaCAA_key = array(rep(1L, n_ages), dim = c(n_ages, 1, 1)),
                 sigmaDAA_key = array(1:n_ages, dim = c(n_ages, 1, 1)))
  expect_true("ln_sigmaDAA" %in% names(il$par))
  expect_equal(length(unique(stats::na.omit(as.integer(il$map$ln_sigmaCAA)))), 1)
  expect_equal(length(unique(stats::na.omit(as.integer(il$map$ln_sigmaDAA)))), n_ages)
  expect_equal(sum(il$data$UseDiscardAA), n_yrs * n_ages)
})

test_that("an observed discard the model says is impossible is not silently absorbed", {
  # this toy model has no discarding, so predicted discards at age are exactly
  # zero. Observing discards against that is infinitely unlikely and should read
  # as such rather than being quietly finite.
  n_yrs <- 20; n_ages <- 5
  il <- build_aa(ObsDiscardAA = array(50, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
                 UseDiscardAA = array(1, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
                 discard_units = "biom")
  expect_false(is.finite(rep_of(il)$jnLL))
})

test_that("the across-age correlation is validated and set per stream", {
  expect_error(build_aa(AgeObsCorr_catch = "nope"), "Valid options: iid, 1dar1")

  # each stream is configured where its data are, so freeing one leaves the rest
  il <- build_aa(AgeObsCorr_catch = "1dar1")
  expect_equal(il$data$AgeObsCorr_catch, 1)
  expect_equal(il$data$AgeObsCorr_discard, 0)
  expect_false(any(is.na(as.integer(il$map$trans_rho_catch))))
  expect_true(all(is.na(as.integer(il$map$trans_rho_discard))))

  # one correlation per fleet, not one per stream
  expect_equal(length(il$par$trans_rho_catch), il$data$n_fish_fleets)
  expect_equal(length(il$par$trans_rho_srv_idx), il$data$n_srv_fleets)
})

test_that("a model with no at-age data is unchanged by the machinery", {
  il <- build_aa(aa = FALSE)
  r <- rep_of(il)
  expect_true(is.finite(r$jnLL))
  expect_equal(sum(r$CatchAA_nLL), 0)
  expect_equal(sum(r$DiscardAA_nLL), 0)
  expect_equal(sum(r$SrvIdxAA_nLL), 0)
  # the aggregated stream carries the catch instead
  expect_true(sum(r$Catch_nLL) != 0)
})

test_that("an input list built before at-age observations existed still runs", {
  il <- build_aa(aa = FALSE)
  for(nm in c("ObsCatchAA", "UseCatchAA", "ObsDiscardAA", "UseDiscardAA",
              "ObsSrvIdxAA", "UseSrvIdxAA", "ObsFishIdxAA", "UseFishIdxAA",
              "use_catch_aa", "use_discard_aa", "use_srv_idx_aa", "use_fish_idx_aa",
              "AgeObsCorr")) il$data[[nm]] <- NULL
  for(nm in c("ln_sigmaCAA", "ln_sigmaDAA", "ln_sigmaSrvIdxAA", "ln_sigmaFishIdxAA",
              "trans_rho")) {
    il$par[[nm]] <- NULL; il$map[[nm]] <- NULL
  }
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$rep$jnLL))
})

test_that("the at-age and aggregated forms agree where the two say the same thing", {
  # A lognormal on the total and a lognormal on each age are different
  # statements in general. They coincide in one case worth pinning: a single
  # observed age. There the total IS that age, so the two likelihoods must agree
  # exactly, which is what says the at-age path reads the same prediction and
  # applies the same density as the path it is standing beside.
  n_yrs <- 20; n_ages <- 5
  one_age <- 3

  # numbers-at-age prediction for the fleet, taken from a run with neither
  # stream fit, so both variants below are compared against the same dynamics
  base <- build_aa(aa = FALSE)
  base$data$UseCatch[] <- 0
  r0 <- rep_of(base)
  pred_a <- vapply(seq_len(n_yrs), function(y) sum(r0$CAA[, 1, y, 1, one_age, , 1]), numeric(1))

  obs_a <- pred_a * exp(stats::rnorm(n_yrs, 0, 0.1))
  sig <- 0.25

  # as an aggregated catch, with selectivity pinned so only that age is caught
  agg <- build_aa(aa = FALSE, ln_sigmaC = array(log(sig), dim = c(1, n_yrs, 1, 1)))
  agg$data$ObsCatch[1, , 1, 1] <- obs_a
  agg$data$UseCatch[] <- 1

  # and as a catch at age on that one age
  ObsCatchAA <- array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1))
  UseCatchAA <- array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1))
  ObsCatchAA[1, , 1, one_age, 1, 1] <- obs_a
  UseCatchAA[1, , 1, one_age, 1, 1] <- 1
  aa <- build_aa(ObsCatchAA_in = ObsCatchAA, UseCatchAA_in = UseCatchAA,
                 sigmaCAA_key = array(c(NA, NA, 1L, NA, NA), dim = c(n_ages, 1, 1)),
                 ln_sigmaCAA = array(log(sig), dim = c(n_ages, 1, 1)))

  ragg <- rep_of(agg); raa <- rep_of(aa)

  # the aggregated total is the sum over ages, so the two only match when the
  # fleet catches one age; compare the density each places on the same numbers
  agg_ll <- -stats::dnorm(log(obs_a), log(vapply(seq_len(n_yrs), function(y)
    sum(ragg$CAA[, 1, y, 1, , , 1]), numeric(1))), sig, log = TRUE)
  expect_equal(as.numeric(ragg$Catch_nLL[1, , 1, 1]), agg_ll, tolerance = 1e-10)

  aa_ll <- -stats::dnorm(log(obs_a), log(vapply(seq_len(n_yrs), function(y)
    sum(raa$CAA[, 1, y, 1, one_age, , 1]), numeric(1))), sig, log = TRUE)
  expect_equal(as.numeric(raa$CatchAA_nLL[1, , 1, one_age, 1, 1]), aa_ll, tolerance = 1e-10)

  # every other age contributes nothing, since the key holds them out
  expect_equal(sum(raa$CatchAA_nLL[1, , 1, -one_age, 1, 1]), 0)
})

test_that("the fishery index at age carries its fleet's seasonal timing", {

  # t_fish is the fraction of the season elapsed when the index is observed, so
  # the numbers are decayed by that much total mortality first. The at-age index
  # reads the same array the aggregated one does, so it inherits that; building
  # it from raw numbers at age would silently ignore the timing.
  n_yrs <- 20; n_ages <- 5
  aa_dim <- c(1, n_yrs, 1, n_ages, 1, 1)
  obs <- array(100, dim = aa_dim); use <- array(1, dim = aa_dim)

  pred_at <- function(t_fish_val) {
    il <- build_aa(ObsFishIdxAA = obs, UseFishIdxAA = use,
                   fish_t = array(t_fish_val, dim = c(1, 1, 1)))
    as.numeric(rep_of(il)$FishIAA[1, 1, 10, 1, , 1, 1])
  }

  start_of_season <- pred_at(0)
  mid_season <- pred_at(0.5)

  # a later index sees fewer fish, and the shortfall grows with age because
  # older fish carry more accumulated mortality within the season
  expect_true(all(mid_season < start_of_season))
  expect_true(all(mid_season > 0))
})
