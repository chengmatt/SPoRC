# Shared routines for the semi-parametric growth self-test. Mean length at age moves by a KNOWN
# year-by-age surface, and the estimating model has to recover it from the parametric curve alone.
#
# The surface is smooth in both directions rather than white noise, which is what the correlated process
# errors are for: a 2D AR(1) or 3D GMRF borrows strength across neighboring ages and years.

spcfg <- list(
  n_yrs = 30,
  n_ages = 12,
  n_sims = 1,
  idx_se = 0.1,
  comp_iss = 400,
  caal_per_bin = 40,
  len_lower = seq(10, 75, by = 5),
  f_ramp = c(seq(0.05, 0.45, length.out = 18), seq(0.45, 0.15, length.out = 12))
)
spcfg$n_lens <- length(spcfg$len_lower)

# the parametric curve underneath, in the Schnute form the growth module reads
sp_pars <- c(L1 = 14, L2 = 68, K = 0.25, CV1 = 0.12, CV2 = 0.07)

#' The true deviation surface: a wave over years scaled by a gradient over ages
#'
#' Older fish depart more than young ones, and the departure swings slowly over
#' the series, so a correlated process error has real structure to find.
#'
#' @param n_yrs,n_ages dimensions
#' @param amp largest departure in mean length at age, on the log scale
#' @keywords internal
sp_true_devs <- function(n_yrs = spcfg$n_yrs, n_ages = spcfg$n_ages, amp = 0.12) {
  wave <- sin(2 * pi * (seq_len(n_yrs) - 1) / (n_yrs / 1.5))
  grad <- seq(0.25, 1, length.out = n_ages)
  amp * outer(wave, grad)
}

#' Binned normal age-length key, as the growth module builds one
#' @keywords internal
sp_alk <- function(len_lower, mu, sd) {
  n_lens <- length(len_lower); n_ages <- length(mu)
  alk <- matrix(0, n_lens, n_ages)
  for(a in 1:n_ages) {
    cdf <- stats::pnorm((len_lower - mu[a]) / sd[a])
    col <- c(cdf[2:n_lens], 1) - cdf
    col[1] <- col[1] + cdf[1]
    alk[,a] <- col
  } # end a loop
  alk
}

#' Mean and spread of length at age under the parametric curve, by year
#'
#' The curve read at integer ages with the Schnute reference ages at the first
#' and last age, times the deviation surface. Uses the package's own
#' get_laa_curve so the operating model and the estimating model agree on the
#' curve and differ only in the deviations.
#'
#' @keywords internal
sp_growth <- function(devs = NULL) {
  n_yrs <- spcfg$n_yrs; n_ages <- spcfg$n_ages; ages <- 1:n_ages
  crv <- get_laa_curve(
    x = ages,
    L0 = spcfg$len_lower[1],
    L1 = sp_pars[["L1"]],
    L2 = sp_pars[["L2"]],
    K = sp_pars[["K"]],
    CV1 = sp_pars[["CV1"]],
    CV2 = sp_pars[["CV2"]],
    A1 = 1,
    A2 = n_ages
  )
  if(is.null(devs)) devs <- matrix(0, n_yrs, n_ages)
  mu <- sweep(exp(devs), 2, crv$L, "*")          # [year, age]
  sd <- sweep(mu, 2, crv$cv, "*")                # the CV at age is untouched
  list(mean = mu, sd = sd, cv = crv$cv)
}

#' Operating model whose size at age moves by a known surface
#'
#' @param seed random seed for the draws
#' @keywords internal
semipar_simulate <- function(seed = 11) {

  n_yrs <- spcfg$n_yrs; n_ages <- spcfg$n_ages; n_lens <- spcfg$n_lens

  sim_list <- Setup_Sim_Dim(
    n_sims = 1,
    n_yrs = n_yrs,
    n_regions = 1,
    n_ages = n_ages,
    n_lens = n_lens,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1
  )
  sim_list <- Setup_Sim_Containers(sim_list)

  curve7 <- function(slope, infl) array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
                                        dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
  iss_caal <- array(spcfg$caal_per_bin, dim = c(1, n_yrs, 1, n_lens, 1, 1, 1))
  type_mat <- array(1, dim = c(n_yrs, 1)) # split by region and sex

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, curve7(3, 2)),
    ret_sel_input = replicate(1, curve7(3, 2)),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(spcfg$f_ramp, dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(spcfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)),
    ISS_FishLenComps = array(spcfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)),
    FishAgeComps_Type = array(1, dim = c(n_yrs, 1)),
    FishLenComps_Type = array(1, dim = c(n_yrs, 1))
  )
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(1, curve7(1.2, 3)),
    t_srv = array(0, dim = c(1, 1, 1)),
    ObsSrvIdx_SE = array(spcfg$idx_se, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(spcfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)),
    ISS_SrvLenComps = array(spcfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)),
    SrvAgeComps_Type = array(1, dim = c(n_yrs, 1)),
    SrvLenComps_Type = array(1, dim = c(n_yrs, 1)),
    comp_srv_caal_like = "Multinomial",
    ISS_Srv_caal = iss_caal,
    Srv_caal_Type = type_mat,
    ln_Srv_caal_theta = array(log(1), dim = c(1, 1, 1)),
    ln_Srv_caal_theta_agg = log(1)
  )

  # the size-age transition the operating model runs on: one key per year, built
  # from the parametric curve times the true deviation surface
  devs <- sp_true_devs()
  g <- sp_growth(devs)
  sat <- array(0, dim = c(1, 1, n_yrs, 1, n_lens, n_ages, 1, 1))
  for(y in 1:n_yrs) sat[1,1,y,1,,,1,1] <- sp_alk(spcfg$len_lower, g$mean[y, ], g$sd[y, ])

  # weight at age follows the same keys, so the biomass moves with growth too
  w_mid <- spcfg$len_lower + 2.5
  waa <- t(sapply(1:n_yrs, function(y) as.vector(t(sat[1,1,y,1,,,1,1]) %*% (1e-5 * w_mid^3))))
  mat <- 1 / (1 + exp(-1.5 * ((1:n_ages) - 4)))
  waa6 <- array(waa, dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(0.25, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, waa6),
    WAA_fish_input = replicate(1, array(waa, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(waa, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    MatAA_input = replicate(1, array(rep(mat, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))),
    SizeAgeTrans_input = sat
  ))
  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(1, array(8, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(0.3), dim = c(2, 1, 1)),
    sexratio_input = replicate(1, array(1, dim = c(1, 1, n_yrs, 1))),
    recruitment_opt = "mean_rec",
    init_age_strc = 1
  )

  set.seed(seed)
  om <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
  obs <- simulation_data_to_SPoRC(sim_env = om, y = n_yrs, sim = 1)
  list(obs = obs, devs = devs, mean_LAA = g$mean, sd_LAA = g$sd, om = om)
}

#' Estimating model with a semi-parametric growth surface
#'
#' Growth is estimated from the parametric curve up, plus a deviation surface
#' under the named process error. Survey lengths and conditional age-at-length
#' are the data that inform it; no marginal age compositions are fit, so the
#' size at age has to come from the length data and the age-at-length rows.
#'
#' @param form one of none, iid, rw, 2dar1, 3dmarg, 3dcond
#' @param obs the observation list from semipar_simulate(); a fresh one is
#'   simulated when absent
#' @keywords internal
semipar_input <- function(form = "2dar1", obs = NULL) {

  if(is.null(obs)) obs <- semipar_simulate()$obs
  n_yrs <- spcfg$n_yrs; n_ages <- spcfg$n_ages

  input_list <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = spcfg$len_lower + 2.5,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(0.3), c(2, 1, 1)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(8)
  )
  input_list <- suppressMessages(Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = NULL,
    MatAA = obs$MatAA,
    fit_lengths = 1,
    SizeAgeTrans = NA,
    comp_const_obs = 0,
    AgeingError = obs$AgeingError,
    M_spec = "fix",
    Fixed_natmort = array(0.25, dim = c(1, 1, n_yrs, n_ages, 1)),
    growth_model = "vb_schnute",
    growth_spec = "est_all",
    ln_growth_pars = array(log(sp_pars), dim = c(1, 1, 1, 5)),
    growth_A1 = 1,
    growth_A2 = n_ages,
    growth_len_lower = spcfg$len_lower,
    growth_plus_group = "curve",
    waa_model = "wt_len",
    wt_len_pars = c(1e-5, 3),
    growth_semipar = form,
    growth_semipar_spec = "fix",
    do_caal = 1
  ))
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )
  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = obs$ObsCatch,
    UseCatch = obs$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    ln_sigmaC = obs$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(1, 1, 1)),
    ObsDiscard = obs$ObsDiscard,
    UseDiscard = obs$UseDiscard,
    sigma_dmr_spec = "fix",
    dmr_mean_spec = "fix",
    ln_sigmaD = obs$ln_sigmaD
  ))

  no_use <- array(0, dim = dim(obs$UseFishAgeComps))
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = obs$ObsFishIdx,
    ObsFishIdx_SE = obs$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(obs$UseFishIdx)),
    fish_idx_type = "biom",
    ObsFishAgeComps = obs$ObsFishAgeComps,
    ISS_FishAgeComps = obs$ISS_FishAgeComps,
    UseFishAgeComps = no_use,
    FishAgeComps_LikeType = "none",
    FishAgeComps_Type = "none_Year_1-terminal_Fleet_1",
    ObsFishLenComps = obs$ObsFishLenComps,
    ISS_FishLenComps = obs$ISS_FishLenComps,
    UseFishLenComps = obs$UseFishLenComps,
    FishLenComps_LikeType = "Multinomial",
    FishLenComps_Type = "spltRspltS_Year_1-terminal_Fleet_1",
    t_fish = array(0, dim = c(1, 1, 1))
  )
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = obs$ObsSrvIdx,
    ObsSrvIdx_SE = obs$ObsSrvIdx_SE,
    UseSrvIdx = obs$UseSrvIdx,
    srv_idx_type = "biom",
    ObsSrvAgeComps = obs$ObsSrvAgeComps,
    ISS_SrvAgeComps = obs$ISS_SrvAgeComps,
    UseSrvAgeComps = array(0, dim = dim(obs$UseSrvAgeComps)),
    SrvAgeComps_LikeType = "none",
    SrvAgeComps_Type = "none_Year_1-terminal_Fleet_1",
    ObsSrvLenComps = obs$ObsSrvLenComps,
    ISS_SrvLenComps = obs$ISS_SrvLenComps,
    UseSrvLenComps = obs$UseSrvLenComps,
    SrvLenComps_LikeType = "Multinomial",
    SrvLenComps_Type = "spltRspltS_Year_1-terminal_Fleet_1",
    ObsSrv_caal = obs$ObsSrv_caal,
    UseSrv_caal = obs$UseSrv_caal,
    ISS_Srv_caal = obs$ISS_Srv_caal,
    Srv_caal_LikeType = "Multinomial",
    Srv_caal_Type = "spltRspltS_Year_1-terminal_Fleet_1"
  )
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "est_all",
    use_fixed_ret_sel = 1
  )
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all",
    t_srv = array(0, dim = c(1, 1, 1))
  )
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
    Wt_FishLenComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
    Wt_SrvLenComps = array(1, dim = c(1, n_yrs, 1, 1, 1))
  )
  input_list
}
