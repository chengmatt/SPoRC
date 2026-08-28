# Shared machinery for the conditional age-at-length self-tests. The operating
# model carries length bins, a size-age transition built the way SS3 builds one
# (binned normal CDF of length at age with the tails accumulated into the end
# bins), and a length-stratified ageing design: a fixed number of otoliths per
# length bin, which is the sampling scheme CAAL exists to handle. The estimation
# model then fits marginal length compositions plus CAAL and no marginal age
# compositions, so every piece of age information it has comes through the CAAL
# likelihood. Marginal length compositions carry the abundance signal, the CAAL
# rows carry P(age | length).

caal_cfg <- list(
  n_yrs = 30, n_ages = 10, n_lens = 12, n_sims = 50,
  idx_se = 0.1, comp_iss = 300, caal_per_bin = 25,
  len_lower = seq(10, 65, by = 5), # lower edges of the length bins
  f_ramp = c(seq(0.05, 0.5, length.out = 18), seq(0.5, 0.12, length.out = 12))
)

#' Binned normal age-length key, SS3's calc_ALK
#'
#' Columns are ages and sum to one. The lower tail lands in the first bin and
#' the upper tail in the last, so every age has all of its mass somewhere.
#'
#' @keywords internal
caal_alk <- function(len_lower, mean_laa, sd_laa) {
  n_lens <- length(len_lower); n_ages <- length(mean_laa)
  alk <- matrix(0, n_lens, n_ages)
  for(a in 1:n_ages) {
    cdf <- stats::pnorm((len_lower - mean_laa[a]) / sd_laa[a])
    col <- c(cdf[2:n_lens], 1) - cdf
    col[1] <- col[1] + cdf[1]
    alk[,a] <- col
  } # end a loop
  alk
}

# von Bertalanffy mean length at age with a constant CV, per sex
caal_growth <- function(n_ages, n_sexes) {
  linf <- c(60, 52)[1:n_sexes]; k <- c(0.30, 0.35)[1:n_sexes]; t0 <- -0.5; cv <- 0.10
  lapply(1:n_sexes, function(s) {
    mu <- linf[s] * (1 - exp(-k[s] * ((1:n_ages) - t0)))
    list(mean = mu, sd = cv * mu)
  })
}

#' Operating model with length bins and conditional age-at-length sampling
#'
#' @param caal_like "Multinomial" or "Dirichlet-Multinomial" for the CAAL draws.
#' @param caal_type composition type string for the CAAL draws.
#' @param ln_theta DM log overdispersion used when caal_like is DM.
#' @param caal_bins length bins that receive aged fish; others get none.
#' @param n_sexes one or two sexes.
#' @param n_sims replicates to simulate.
#' @keywords internal
caal_make_om <- function(caal_like = "Multinomial", caal_type = "spltRspltS", ln_theta = log(1),
                         caal_bins = NULL, n_sexes = 1, n_sims = 1, seed = 55) {

  n_yrs <- caal_cfg$n_yrs; n_ages <- caal_cfg$n_ages; n_lens <- caal_cfg$n_lens
  if(is.null(caal_bins)) caal_bins <- 1:n_lens

  sim_list <- Setup_Sim_Dim(n_sims = n_sims, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages,
                            n_lens = n_lens, n_sexes = n_sexes, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)

  curve7 <- function(slope, infl, scale = 1, n_fleets = 1) {
    array(rep(scale / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
          dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, n_fleets))
  }

  # a fixed number of otoliths per length bin in the bins that are sampled
  iss_caal <- array(0, dim = c(1, n_yrs, 1, n_lens, n_sexes, 1, n_sims))
  iss_caal[,,,caal_bins,,,] <- caal_cfg$caal_per_bin
  type_mat <- array(c(agg = 0, spltRspltS = 1, spltRjntS = 2)[[caal_type]], dim = c(n_yrs, 1))

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(n_sims, curve7(3, 2)),
    ret_sel_input = replicate(n_sims, curve7(3, 2)),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, n_sims)),
    Fmort_input = array(rep(caal_cfg$f_ramp, each = 1), dim = c(1, n_yrs, 1, 1, n_sims)),
    # marginal ages are drawn but never fit; marginal lengths are fit
    ISS_FishAgeComps = array(caal_cfg$comp_iss, dim = c(1, n_yrs, 1, n_sexes, 1, n_sims)),
    ISS_FishLenComps = array(caal_cfg$comp_iss, dim = c(1, n_yrs, 1, n_sexes, 1, n_sims)),
    FishAgeComps_Type = array(if(n_sexes == 1) 1 else 2, dim = c(n_yrs, 1)),
    FishLenComps_Type = array(if(n_sexes == 1) 1 else 2, dim = c(n_yrs, 1)),
    comp_fish_caal_like = caal_like, ISS_Fish_caal = iss_caal, Fish_caal_Type = type_mat,
    ln_Fish_caal_theta = array(ln_theta, dim = c(1, n_sexes, 1)),
    ln_Fish_caal_theta_agg = ln_theta
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list, srv_sel_input = replicate(n_sims, curve7(1, 3)),
    t_srv = array(0, dim = c(1, 1, 1)), # start of the season, where the fixture's key is read
    ObsSrvIdx_SE = array(caal_cfg$idx_se, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(caal_cfg$comp_iss, dim = c(1, n_yrs, 1, n_sexes, 1, n_sims)),
    ISS_SrvLenComps = array(caal_cfg$comp_iss, dim = c(1, n_yrs, 1, n_sexes, 1, n_sims)),
    SrvAgeComps_Type = array(if(n_sexes == 1) 1 else 2, dim = c(n_yrs, 1)),
    SrvLenComps_Type = array(if(n_sexes == 1) 1 else 2, dim = c(n_yrs, 1)),
    comp_srv_caal_like = caal_like, ISS_Srv_caal = iss_caal, Srv_caal_Type = type_mat,
    ln_Srv_caal_theta = array(ln_theta, dim = c(1, n_sexes, 1)),
    ln_Srv_caal_theta_agg = ln_theta
  )

  # size-age transition from the growth curve, sex specific, constant over time
  growth <- caal_growth(n_ages, n_sexes)
  sat <- array(0, dim = c(1, 1, n_yrs, 1, n_lens, n_ages, n_sexes, n_sims))
  for(s in 1:n_sexes) {
    alk <- caal_alk(caal_cfg$len_lower, growth[[s]]$mean, growth[[s]]$sd)
    for(y in 1:n_yrs) for(i in 1:n_sims) sat[1,1,y,1,,,s,i] <- alk
  } # end s loop

  waa <- 5 / (1 + exp(-3 * ((1:n_ages) - 3))); mat <- 1 / (1 + exp(-3 * ((1:n_ages) - 3)))
  biol6 <- function(val) array(rep(val, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, n_sexes))
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(n_sims, array(0.3, dim = c(1, 1, n_yrs, n_ages, n_sexes))),
    WAA_input = replicate(n_sims, biol6(waa)),
    WAA_fish_input = replicate(n_sims, array(rep(waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 1))),
    WAA_srv_input = replicate(n_sims, array(rep(waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, n_sexes, 1))),
    MatAA_input = replicate(n_sims, biol6(mat)),
    SizeAgeTrans_input = sat
  ))

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, n_sexes, n_sims))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(n_sims, array(5, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(0.3), dim = c(2, 1, 1)),
    sexratio_input = replicate(n_sims, array(1 / n_sexes, dim = c(1, 1, n_yrs, n_sexes))),
    recruitment_opt = "mean_rec", init_age_strc = 1
  )

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}

#' Estimation model that sees lengths and CAAL but no marginal ages
#'
#' @param caal_like likelihood for the CAAL streams.
#' @param caal_type composition type string for the CAAL streams.
#' @param use_caal whether CAAL is fit at all; FALSE gives a lengths only model.
#' @param use_age_comps whether marginal age comps are fit instead of CAAL.
#' @param osa whether to track compositions for one step ahead residuals.
#' @keywords internal
caal_build_input <- function(sim_data, caal_like = "Multinomial", caal_type = "spltRspltS",
                             use_caal = TRUE, use_age_comps = FALSE, osa = FALSE, n_sexes = 1) {

  n_yrs <- caal_cfg$n_yrs; n_ages <- caal_cfg$n_ages; n_lens <- caal_cfg$n_lens
  ct <- if(n_sexes == 1) "spltRspltS" else "spltRjntS"

  input_list <- Setup_Mod_Dim(years = 1:n_yrs, ages = 1:n_ages, lens = caal_cfg$len_lower + 2.5,
                              n_regions = 1, n_sexes = n_sexes, n_fish_fleets = 1, n_srv_fleets = 1,
                              n_pop = 1, natal_region = 1, verbose = FALSE,
                              do_internal_comp_osa = osa)
  input_list <- Setup_Mod_Rec(input_list = input_list, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                              ln_sigmaR = array(log(0.3), c(2, 1, 1)), rec_model = "mean_rec",
                              sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2,
                              ln_global_R0 = log(5))
  input_list <- Setup_Mod_Biologicals(input_list = input_list, WAA = sim_data$WAA,
                                      MatAA = sim_data$MatAA, WAA_fish = sim_data$WAA_fish,
                                      WAA_srv = sim_data$WAA_srv, fit_lengths = 1,
                                      SizeAgeTrans = sim_data$SizeAgeTrans, comp_const_obs = 0,
                                      AgeingError = sim_data$AgeingError, M_spec = "fix",
                                      Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, n_sexes)))
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1,
                                   Fixed_Movement = NA, do_recruits_move = 0)
  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list, ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch,
    Use_F_pen = 1, sigmaC_spec = "fix", ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(1, 1, 1)),
    ObsDiscard = sim_data$ObsDiscard, UseDiscard = sim_data$UseDiscard,
    sigma_dmr_spec = "fix", dmr_mean_spec = "fix", ln_sigmaD = sim_data$ln_sigmaD))

  no_use <- array(0, dim = dim(sim_data$UseFishAgeComps))
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)), fish_idx_type = "biom",
    ObsFishAgeComps = sim_data$ObsFishAgeComps, ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    UseFishAgeComps = if(use_age_comps) sim_data$UseFishAgeComps else no_use,
    FishAgeComps_LikeType = if(use_age_comps) "Multinomial" else "none",
    FishAgeComps_Type = paste0(if(use_age_comps) ct else "none", "_Year_1-terminal_Fleet_1"),
    ObsFishLenComps = sim_data$ObsFishLenComps, ISS_FishLenComps = sim_data$ISS_FishLenComps,
    UseFishLenComps = sim_data$UseFishLenComps,
    FishLenComps_LikeType = "Multinomial", FishLenComps_Type = paste0(ct, "_Year_1-terminal_Fleet_1"),
    ObsFish_caal = if(use_caal) sim_data$ObsFish_caal else NULL,
    UseFish_caal = if(use_caal) sim_data$UseFish_caal else NULL,
    ISS_Fish_caal = if(use_caal) sim_data$ISS_Fish_caal else NULL,
    Fish_caal_LikeType = if(use_caal) caal_like else "none",
    Fish_caal_Type = paste0(if(use_caal) caal_type else "none", "_Year_1-terminal_Fleet_1"))
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx, srv_idx_type = "biom",
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    UseSrvAgeComps = if(use_age_comps) sim_data$UseSrvAgeComps else no_use,
    SrvAgeComps_LikeType = if(use_age_comps) "Multinomial" else "none",
    SrvAgeComps_Type = paste0(if(use_age_comps) ct else "none", "_Year_1-terminal_Fleet_1"),
    ObsSrvLenComps = sim_data$ObsSrvLenComps, ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    UseSrvLenComps = sim_data$UseSrvLenComps,
    SrvLenComps_LikeType = "Multinomial", SrvLenComps_Type = paste0(ct, "_Year_1-terminal_Fleet_1"),
    ObsSrv_caal = if(use_caal) sim_data$ObsSrv_caal else NULL,
    UseSrv_caal = if(use_caal) sim_data$UseSrv_caal else NULL,
    ISS_Srv_caal = if(use_caal) sim_data$ISS_Srv_caal else NULL,
    Srv_caal_LikeType = if(use_caal) caal_like else "none",
    Srv_caal_Type = paste0(if(use_caal) caal_type else "none", "_Year_1-terminal_Fleet_1"))
  input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list, fish_sel_model = "logist1_Fleet_1",
                                        fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "est_all",
                                        use_fixed_ret_sel = 1)
  input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list, srv_sel_model = "logist1_Fleet_1",
                                       srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "est_all",
                                       t_srv = array(0, dim = c(1, 1, 1)))
  input_list <- Setup_Mod_Weighting(
    input_list = input_list, Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
    Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)),
    Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, n_sexes, 1)))
  input_list
}

# Fit, self-test, and summarize median relative error per reported quantity,
# the same shape selftest_run() in helper-selftest_features.R returns.
caal_run <- function(input_list, what = c("SSB", "Rec"), n_sims = caal_cfg$n_sims, seed = 1) {
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL, silent = TRUE)
  sd_rep <- RTMB::sdreport(fit)
  stopifnot(max(abs(fit$gr(fit$env$last.par.best))) < 1e-3)

  set.seed(seed)
  res <- simulation_self_test(
    data = fit$data, parameters = input_list$par, mapping = input_list$map,
    random = NULL, rep = fit$rep, sd_rep = sd_rep,
    n_sims = n_sims, what = what, sim_recruitment = "input"
  )

  summ <- lapply(what, function(w) {
    truth <- as.numeric(fit$rep[[w]])
    est_mat <- matrix(res[[which(what == w)]], nrow = length(truth), ncol = n_sims)
    re <- sweep(est_mat, 1, truth, "-") / matrix(truth, length(truth), n_sims)
    c(median_RE = stats::median(re, na.rm = TRUE),
      q10_RE = stats::quantile(re, 0.1, na.rm = TRUE, names = FALSE),
      q90_RE = stats::quantile(re, 0.9, na.rm = TRUE, names = FALSE),
      n_ok = sum(!is.na(est_mat[1,])))
  })
  names(summ) <- what
  list(fit = fit, res = res, summ = summ)
}
