# Operating and estimation model pair with TIME-VARYING weight at age, so the year the per-recruit
# calculation is taken at actually changes S0. Flat biology gives the same phi0 at every reference year.

sr_ref_cfg <- list(n_yrs = 30, n_ages = 12, M = 0.2, sigmaR = 0.3, h = 0.7)

# weight at age doubles linearly across the series
sr_ref_waa <- function(y) {
  seq(0.5, 3, length.out = sr_ref_cfg$n_ages) * (1 + (y - 1) / (sr_ref_cfg$n_yrs - 1))
}
sr_ref_mat <- function() stats::plogis(seq(-3, 3, length.out = sr_ref_cfg$n_ages))

# phi0 straight from its definition, with no SPoRC code involved
sr_ref_phi0 <- function(y) {
  n_ages <- sr_ref_cfg$n_ages; M <- sr_ref_cfg$M
  N <- numeric(n_ages); N[1] <- 1
  for(a in 2:n_ages) N[a] <- N[a - 1] * exp(-M)
  N[n_ages] <- N[n_ages] / (1 - exp(-M))
  sum(N * sr_ref_waa(y) * sr_ref_mat() * 0.5)
}

sr_ref_make_om <- function(SR_ref_yr, R0 = 10, seed = 42) {
  n_yrs <- sr_ref_cfg$n_yrs; n_ages <- sr_ref_cfg$n_ages
  sl <- Setup_Sim_Dim(
    n_sims = 1,
    n_yrs = n_yrs,
    n_regions = 1,
    n_ages = n_ages,
    n_lens = NULL,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1
  )
  sl <- Setup_Sim_Containers(sl)
  curve7 <- function(slope, infl) array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
                                        dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
  sl <- Setup_Sim_Fishing(
    sim_list = sl,
    fish_sel_input = replicate(1, curve7(3, 3)),
    ret_sel_input = replicate(1, curve7(3, 3)),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(rep(seq(0.02, 0.15, length.out = n_yrs), each = 1),
                                              dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(200, dim = c(1, n_yrs, 1, 1, 1, 1))
  )
  sl <- Setup_Sim_Survey(
    sim_list = sl,
    srv_sel_input = replicate(1, curve7(1.5, 4)),
    ObsSrvIdx_SE = array(0.1, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(200, dim = c(1, n_yrs, 1, 1, 1, 1))
  )
  waa6 <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, 1))
  for(y in 1:n_yrs) waa6[1, 1, y, 1, , 1] <- sr_ref_waa(y)
  waa7 <- array(waa6, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))
  mat6 <- array(rep(sr_ref_mat(), each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sl <- Setup_Sim_Biologicals(
    sim_list = sl,
    natmort_input = replicate(1, array(sr_ref_cfg$M, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, waa6),
    WAA_fish_input = replicate(1, waa7),
    WAA_srv_input = replicate(1, waa7),
    MatAA_input = replicate(1, mat6)
  ))
  sl <- Setup_Sim_Tagging(sim_list = sl, use_conv_fish_tagging = 0)
  sl$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  sl <- Setup_Sim_Rec(
    sim_list = sl,
    R0_input = replicate(1, array(R0, dim = c(1, 1, n_yrs))),
    h_input = replicate(1, array(sr_ref_cfg$h, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(sr_ref_cfg$sigmaR), dim = c(2, 1, 1)),
    recruitment_opt = "ricker_rec",
    rec_dd = "global",
    init_age_strc = 1,
    SR_ref_yr = SR_ref_yr
  )
  set.seed(seed)
  suppressMessages(Simulate_Pop_Static(sim_list = sl, output_path = NULL))
}

# Estimation model over the operating model's data, with its own SR_ref_yr so the two can
# be matched or deliberately mismatched.
sr_ref_make_em <- function(sd, SR_ref_yr) suppressWarnings(suppressMessages({
  n_yrs <- sr_ref_cfg$n_yrs; n_ages <- sr_ref_cfg$n_ages
  il <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = NULL,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  il <- Setup_Mod_Rec(
    il,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(sr_ref_cfg$sigmaR), c(2, 1, 1)),
    rec_model = "ricker_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(10),
    rec_lag = 1,
    steepness_h = array(sr_ref_cfg$h, dim = c(1, 1)),
    h_spec = "fix",
    SR_ref_yr = SR_ref_yr
  )
  il <- Setup_Mod_Biologicals(
    il,
    WAA = sd$WAA,
    MatAA = sd$MatAA,
    WAA_fish = sd$WAA_fish,
    WAA_srv = sd$WAA_srv,
    fit_lengths = 0,
    AgeingError = sd$AgeingError,
    M_spec = "fix",
    Fixed_natmort = array(sr_ref_cfg$M, dim = c(1, 1, n_yrs, n_ages, 1))
  )
  il <- Setup_Mod_Tagging(il, use_conv_fish_tagging = 0)
  il <- Setup_Mod_Movement(il, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)
  il <- Setup_Mod_Catch_and_F(
    il,
    ObsCatch = sd$ObsCatch,
    UseCatch = sd$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    ln_sigmaC = sd$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(1, 1, 1))
  )
  il <- Setup_Mod_FishIdx_and_Comps(
    il,
    ObsFishIdx = sd$ObsFishIdx,
    ObsFishIdx_SE = sd$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(sd$UseFishIdx)),
    ObsFishAgeComps = sd$ObsFishAgeComps,
    ObsFishLenComps = NULL,
    UseFishAgeComps = sd$UseFishAgeComps,
    UseFishLenComps = array(0, dim = dim(sd$UseFishAgeComps)),
    ISS_FishAgeComps = sd$ISS_FishAgeComps,
    ISS_FishLenComps = NULL,
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  il <- Setup_Mod_SrvIdx_and_Comps(
    il,
    ObsSrvIdx = sd$ObsSrvIdx,
    ObsSrvIdx_SE = sd$ObsSrvIdx_SE,
    UseSrvIdx = sd$UseSrvIdx,
    SrvIdx_LikeType = "lognormal",
    ObsSrvAgeComps = sd$ObsSrvAgeComps,
    ObsSrvLenComps = NULL,
    UseSrvAgeComps = sd$UseSrvAgeComps,
    UseSrvLenComps = array(0, dim = dim(sd$UseSrvAgeComps)),
    ISS_SrvAgeComps = sd$ISS_SrvAgeComps,
    ISS_SrvLenComps = NULL,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  il <- Setup_Mod_Fishsel_and_Q(
    il,
    fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "est_all",
    use_fixed_ret_sel = 1
  )
  il <- Setup_Mod_Srvsel_and_Q(
    il,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all"
  )
  Setup_Mod_Weighting(
    il,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1))
  )
}))
