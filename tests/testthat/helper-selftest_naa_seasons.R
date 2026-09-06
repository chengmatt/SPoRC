# Two season operating model with a known seasonal process error on the numbers at age, and the
# estimation model over its data.
#
# Seasons are unequal (0.4/0.6) and both are fished and surveyed on purpose: a season with no
# observations has no information about its own state, so a self test would measure the prior.

naaseas_cfg <- list(
  n_yrs = 30,
  n_ages = 8,
  n_seas = 2,
  seasdur = c(0.4, 0.6),
  M = 0.25,
  sigmaR = 0.3,
  idx_se = 0.08,
  comp_iss = 400,
  waa = 5 / (1 + exp(-3 * ((1:8) - 3))),
  mat = 1 / (1 + exp(-3 * ((1:8) - 3))),
  f_ramp = c(seq(0.04, 0.30, length.out = 18), seq(0.30, 0.08, length.out = 12))
)

#' Simulation environment for the innovation draws alone
#'
#' Built directly rather than through the full operating model, because
#' \code{draw_naa_innovations} reads dimensions and process settings and nothing else.
#'
#' @keywords internal
naaseas_env <- function(
  n_pop = 1,
  n_regions = 1,
  n_yrs = 40,
  n_ages = 8,
  n_sexes = 1,
  n_seas = 3,
  ...
) {
  base <- list(
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_ages = n_ages,
    n_sexes = n_sexes,
    n_seas = n_seas,
    naa_re_ages = 1:n_ages,
    naa_re_yrs = 1:n_yrs,
    naa_re_seas = 1:n_seas,
    NAA_re = 1,
    sigmaNAA = 0.4,
    naa_rho = c(age = 0, year = 0, cohort = 0),
    NAA_re_pop = 0,
    NAA_re_region = 0,
    NAA_re_sex = 0,
    NAA_re_season = 0,
    naa_pop_corr = 0,
    naa_region_corr = 0,
    naa_sex_corr = 0,
    naa_season_corr = 0
  )
  list2env(utils::modifyList(base, list(...)))
}

#' Operating model with a seasonal state on the numbers at age
#'
#' @param NAA_re_seasons Passed to \code{Setup_Sim_NAA_state}.
#' @keywords internal
naaseas_make_om <- function(
  NAA_re = "iid",
  sigmaNAA = 0.35,
  NAA_re_seasons = "all",
  NAA_re_season = "iid",
  season_corr = 0,
  seed = 707
) {

  cfg <- naaseas_cfg
  n_yrs <- cfg$n_yrs; n_ages <- cfg$n_ages; n_seas <- cfg$n_seas

  sim_list <- Setup_Sim_Dim(
    n_sims = 1,
    n_yrs = n_yrs,
    n_regions = 1,
    n_ages = n_ages,
    n_lens = NULL,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = n_seas,
    seasdur = cfg$seasdur,
    n_pop = 1,
    natal_region = 1
  )
  sim_list <- Setup_Sim_Containers(sim_list)

  sel <- function(slope, infl) array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))),
                                         each = n_yrs * n_seas),
                                     dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, sel(1.5, 3)),
    ret_sel_input = replicate(1, array(1, dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))),
    dmr_input = array(0, dim = c(1, n_yrs, n_seas, 1, 1)),
    Fmort_input = array(rep(cfg$f_ramp, times = n_seas), dim = c(1, n_yrs, n_seas, 1, 1)),
    ISS_FishAgeComps = array(cfg$comp_iss, dim = c(1, n_yrs, n_seas, 1, 1, 1)),
    ln_sigmaC = array(log(0.01), dim = c(1, n_yrs, n_seas, 1)))

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(1, sel(1, 2.5)),
    ISS_SrvAgeComps = array(cfg$comp_iss, dim = c(1, n_yrs, n_seas, 1, 1, 1)),
    ObsSrvIdx_SE = array(cfg$idx_se, dim = c(1, n_yrs, n_seas, 1)))

  biol <- function(v) array(rep(v, each = n_yrs * n_seas), dim = c(1, 1, n_yrs, n_seas, n_ages, 1))
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, biol(rep(cfg$M, n_ages))),
    WAA_input = replicate(1, biol(cfg$waa)),
    WAA_fish_input = replicate(1, array(rep(cfg$waa, each = n_yrs * n_seas),
                                        dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(rep(cfg$waa, each = n_yrs * n_seas),
                                       dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol(cfg$mat))))

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, n_seas, n_ages, 1, 1))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(1, array(8, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(cfg$sigmaR), dim = c(2, 1, 1)),
    recruitment_opt = "mean_rec",
    init_age_strc = 1
  )

  sim_list <- Setup_Sim_NAA_state(
    sim_list,
    NAA_re = NAA_re,
    sigmaNAA = sigmaNAA,
    NAA_re_seasons = NAA_re_seasons,
    NAA_re_season = NAA_re_season,
    season_corr = season_corr
  )

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}

#' Estimation model over the two season data
#' @keywords internal
naaseas_build_em <- function(sim_data, NAA_re = "none", ...) {

  cfg <- naaseas_cfg
  n_yrs <- dim(sim_data$WAA)[3]; n_ages <- cfg$n_ages; n_seas <- cfg$n_seas

  il <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = NULL,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = n_seas,
    seasdur = cfg$seasdur,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  il <- Setup_Mod_Rec(
    input_list = il,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(cfg$sigmaR), c(2, 1, 1)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(8)
  )
  il <- suppressMessages(do.call(Setup_Mod_Biologicals, utils::modifyList(list(
    input_list = il,
    WAA = sim_data$WAA,
    MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0,
    AgeingError = sim_data$AgeingError,
    M_spec = "fix",
    Fixed_natmort = array(cfg$M, dim = c(1, 1, n_yrs, n_seas, n_ages, 1)),
    NAA_re = NAA_re
  ), list(...))))
  il <- Setup_Mod_Tagging(input_list = il, use_conv_fish_tagging = 0)
  il <- Setup_Mod_Movement(
    input_list = il,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )
  il <- suppressWarnings(Setup_Mod_Catch_and_F(
    input_list = il,
    ObsCatch = sim_data$ObsCatch,
    Catch_Type = sim_data$Catch_Type,
    UseCatch = sim_data$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    ln_sigmaC = array(log(0.01), dim = c(1, n_yrs, n_seas, 1)),
    ln_sigmaF = array(log(1), dim = c(1, n_seas, 1))
  ))
  il <- Setup_Mod_FishIdx_and_Comps(
    input_list = il,
    ObsFishIdx = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
    ObsFishAgeComps = sim_data$ObsFishAgeComps,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ObsFishLenComps = NULL,
    UseFishLenComps = array(0, dim = dim(sim_data$UseFishAgeComps)),
    ISS_FishLenComps = NULL,
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  il <- Setup_Mod_SrvIdx_and_Comps(
    input_list = il,
    ObsSrvIdx = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    SrvIdx_LikeType = "lognormal",
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    ObsSrvLenComps = NULL,
    UseSrvLenComps = array(0, dim = dim(sim_data$UseSrvAgeComps)),
    ISS_SrvLenComps = NULL,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  il <- Setup_Mod_Fishsel_and_Q(
    input_list = il,
    fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "est_all",
    use_fixed_ret_sel = 1
  )
  il <- Setup_Mod_Srvsel_and_Q(
    input_list = il,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all"
  )
  fw <- array(1, dim = c(1, n_yrs, n_seas, 1, 1))
  Setup_Mod_Weighting(
    input_list = il,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_FishAgeComps = fw,
    Wt_SrvAgeComps = fw
  )
}

#' The operating model's observations, in the shape the estimation model reads
#' @keywords internal
naaseas_om_data <- function(om) simulation_data_to_SPoRC(sim_env = om, y = naaseas_cfg$n_yrs, sim = 1)

#' Relative root mean squared error of estimated against true spawning biomass
#' @keywords internal
naaseas_ssb_rmse <- function(fit, om) {
  truth <- as.vector(om$SSB[1, 1, , 1])
  sqrt(mean(((as.vector(fit$rep$SSB[1, 1, ]) - truth) / truth)^2))
}
