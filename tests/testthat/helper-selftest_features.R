# Shared machinery for the feature self-tests (normal/mvn index likelihoods,
# concentrated catchability, Ricker recruitment, and the Ricker reference point
# consistency check). Each test file fits a small single-region model with one
# feature switched on, hands the fit to simulation_self_test so the operating
# model generates data under the same error structure, and asserts median
# relative error. dev/scratch/selftest_new_features.R runs the same blocks
# through this helper as a standalone driver.
#
# The index blocks keep sigmaR and the observation errors at realistic values
# because observation error is the machinery under test there; the Ricker block
# passes near-zero values for both so recovery is held to near-exact tolerances,
# where any residual bias is structural rather than statistical.

selftest_cfg <- list(
  n_yrs = 40, n_ages = 6, n_sims = 50,
  idx_se = 0.1, comp_iss = 300,
  waa = 5 / (1 + exp(-3 * ((1:6) - 3))),
  mat = 1 / (1 + exp(-3 * ((1:6) - 3))),
  # two-way trip so the stock-recruit curve sees SSB contrast; the deep trough
  # is what identifies steepness, whose recovery is skewed under weaker depletion
  f_ramp = c(seq(0.05, 0.5, length.out = 24), seq(0.5, 0.12, length.out = 16))
)

# Operating model for the baseline dataset each estimation model is fitted to.
selftest_make_om <- function(SrvIdx_LikeType = NULL, SrvIdx_Cov = NULL, UseSrvIdx = NULL,
                             recruitment_opt = "mean_rec", sigmaR = 0.3,
                             idx_se_om = selftest_cfg$idx_se, iss_om = selftest_cfg$comp_iss,
                             seed = 55) {

  n_yrs <- selftest_cfg$n_yrs; n_ages <- selftest_cfg$n_ages

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages,
                            n_lens = NULL, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)

  curve7 <- function(slope, infl, scale = 1) {
    array(rep(scale / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
          dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
  }

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, curve7(3, 2)),
    ret_sel_input = replicate(1, curve7(3, 2)),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(rep(selftest_cfg$f_ramp, each = 1), dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(iss_om, dim = c(1, n_yrs, 1, 1, 1, 1))
  )

  srv_args <- list(sim_list = sim_list, srv_sel_input = replicate(1, curve7(1, 3)),
                   ObsSrvIdx_SE = array(idx_se_om, dim = c(1, n_yrs, 1, 1)),
                   ISS_SrvAgeComps = array(iss_om, dim = c(1, n_yrs, 1, 1, 1, 1)))
  if(!is.null(SrvIdx_LikeType)) srv_args$SrvIdx_LikeType <- SrvIdx_LikeType
  if(!is.null(SrvIdx_Cov)) srv_args$SrvIdx_Cov <- SrvIdx_Cov
  if(!is.null(UseSrvIdx)) srv_args$UseSrvIdx <- UseSrvIdx
  sim_list <- do.call(Setup_Sim_Survey, srv_args)

  biol6 <- function(val) array(rep(val, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, biol6(selftest_cfg$waa)),
    WAA_fish_input = replicate(1, array(rep(selftest_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(rep(selftest_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol6(selftest_cfg$mat))
  ))

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  rec_args <- list(
    sim_list = sim_list,
    R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(sigmaR), dim = c(2, 1, 1)),
    recruitment_opt = recruitment_opt,
    init_age_strc = 1
  )
  if(recruitment_opt == "ricker_rec") rec_args$h_input <- replicate(1, array(0.7, dim = c(1, 1, n_yrs)))
  sim_list <- do.call(Setup_Sim_Rec, rec_args)

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}

# Estimation model over the baseline data, with the feature under test switched on.
selftest_build_input <- function(sim_data, SrvIdx_LikeType = "lognormal", SrvIdx_Cov = NULL,
                                 srv_q_type = "est", rec_model = "mean_rec", sigmaR = 0.3) {

  n_yrs <- selftest_cfg$n_yrs; n_ages <- selftest_cfg$n_ages

  input_list <- Setup_Mod_Dim(years = 1:n_yrs, ages = 1:n_ages, lens = NULL,
                              n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                              n_pop = 1, natal_region = 1, verbose = FALSE)

  rec_args <- list(input_list = input_list, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                   ln_sigmaR = array(log(sigmaR), c(2, 1, 1)), rec_model = rec_model,
                   sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2,
                   ln_global_R0 = log(5))
  if(rec_model == "ricker_rec") {
    rec_args$rec_lag <- 1
    rec_args$steepness_h <- array(0.7, dim = c(1, 1))
    rec_args$h_spec <- "est_shared_pop_r"
  }
  input_list <- do.call(Setup_Mod_Rec, rec_args)

  input_list <- Setup_Mod_Biologicals(input_list = input_list, WAA = sim_data$WAA,
                                      MatAA = sim_data$MatAA, WAA_fish = sim_data$WAA_fish,
                                      WAA_srv = sim_data$WAA_srv, fit_lengths = 0,
                                      AgeingError = sim_data$AgeingError, M_spec = "fix",
                                      Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, 1)))
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)
  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list, ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch,
    Use_F_pen = 1, sigmaC_spec = "fix", ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(1, 1, 1)),
    ObsDiscard = sim_data$ObsDiscard, UseDiscard = sim_data$UseDiscard,
    # the operating model carries no discard data, so an estimated dmr mean is
    # unidentified and puts an exactly zero row in the Hessian. That makes the
    # Newton refinement's solve() fail, leaving the fit wherever nlminb stopped
    # and the gradient check platform dependent, so dmr stays fixed here
    sigma_dmr_spec = "fix", dmr_mean_spec = "fix", ln_sigmaD = sim_data$ln_sigmaD))
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
    ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = NULL,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    UseFishLenComps = array(0, dim = dim(sim_data$UseFishAgeComps)),
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = NULL,
    fish_idx_type = "biom", FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none", FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1")
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    SrvIdx_LikeType = SrvIdx_LikeType, SrvIdx_Cov = SrvIdx_Cov,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = NULL,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    UseSrvLenComps = array(0, dim = dim(sim_data$UseSrvAgeComps)),
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = NULL,
    srv_idx_type = "biom", SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none", SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1")
  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list, fish_sel_model = "logist1_Fleet_1",
                                        fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "est_all",
                                        use_fixed_ret_sel = 1)
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list, srv_sel_model = "logist1_Fleet_1",
                                       srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "est_all",
                                       srv_q_type = srv_q_type)
  input_list <- Setup_Mod_Weighting(
    input_list = input_list, Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1,
    Wt_Rec = 1, Wt_F = 1,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)))
  input_list
}

# Fit the base model, run the self-test with the block's own seed, and summarize
# median relative error per reported quantity. Assertions are the caller's job,
# so the testthat files use expect_lt while the dev driver stops on violation.
selftest_run <- function(input_list, what, sim_recruitment = "input", seed = 1) {
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL, silent = TRUE)
  sd_rep <- RTMB::sdreport(fit)
  stopifnot(max(abs(fit$gr(fit$env$last.par.best))) < 1e-3)

  # each block draws its replicates from its own seed, so results do not depend
  # on what ran before it
  set.seed(seed)
  res <- simulation_self_test(
    data = fit$data, parameters = input_list$par, mapping = input_list$map,
    random = NULL, rep = fit$rep, sd_rep = sd_rep,
    n_sims = selftest_cfg$n_sims, what = what, sim_recruitment = sim_recruitment
  )

  summ <- lapply(what, function(w) {
    truth <- as.numeric(fit$rep[[w]])
    est <- res[[which(what == w)]]
    est_mat <- matrix(est, nrow = length(truth), ncol = selftest_cfg$n_sims)
    re <- sweep(est_mat, 1, truth, "-") / matrix(truth, length(truth), selftest_cfg$n_sims)
    c(median_RE = stats::median(re, na.rm = TRUE),
      q10_RE = stats::quantile(re, 0.1, na.rm = TRUE, names = FALSE),
      q90_RE = stats::quantile(re, 0.9, na.rm = TRUE, names = FALSE))
  })
  names(summ) <- what

  list(fit = fit, res = res, summ = summ)
}

# The mvn blocks share this covariance: a survey-wide scaling error with a
# strong common factor and marginal sd at fraction `scale` of the index level.
selftest_mvn_cov <- function(idx_scale, scale = 0.15, lambda = 0.7) {
  n_yrs <- selftest_cfg$n_yrs
  lam <- rep(lambda, n_yrs); d <- rep(scale * idx_scale, n_yrs)
  R <- outer(lam, lam); diag(R) <- 1
  outer(d, d) * R
}
