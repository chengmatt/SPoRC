# Operating model to estimation model machinery for the state-space numbers at age.
#
# The operating model advances the population with a known process error on the numbers at age,
# generates index and composition data from the result, and the estimation model then has to
# recover the process from those data alone, with the states integrated out as random effects.
#
# This is a different question from the conditional recovery in test-selftest_naa_state.R. There the
# true states were handed over and only the process parameters were estimated. Here nothing about
# the states is known, and the variance has to be separated from observation error, which is the
# part of a state-space fit that actually fails in practice.

naaom_cfg <- list(
  n_yrs = 35, n_ages = 8, idx_se = 0.08, comp_iss = 400, M = 0.25, sigmaR = 0.3,
  waa = 5 / (1 + exp(-3 * ((1:8) - 3))),
  mat = 1 / (1 + exp(-3 * ((1:8) - 3))),
  # a two-way trip so the population is not in a single state through the series
  f_ramp = c(seq(0.05, 0.45, length.out = 20), seq(0.45, 0.12, length.out = 15))
)

#' Operating model carrying a known numbers-at-age process error
#' @keywords internal
naaom_make_om <- function(NAA_re = "2dar1", sigmaNAA = 0.25, rho_age = 0.5, rho_year = 0.4,
                          seed = 808) {
  n_yrs <- naaom_cfg$n_yrs; n_ages <- naaom_cfg$n_ages
  curve <- function(slope, infl, scale = 1)
    array(rep(scale / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
          dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))

  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages,
                            n_lens = NULL, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list, fish_sel_input = replicate(1, curve(3, 2)),
    ret_sel_input = replicate(1, curve(3, 2)),
    dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
    Fmort_input = array(naaom_cfg$f_ramp, dim = c(1, n_yrs, 1, 1, 1)),
    ISS_FishAgeComps = array(naaom_cfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)))
  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list, srv_sel_input = replicate(1, curve(1, 3)),
    ObsSrvIdx_SE = array(naaom_cfg$idx_se, dim = c(1, n_yrs, 1, 1)),
    ISS_SrvAgeComps = array(naaom_cfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)))
  biol <- function(v) array(rep(v, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, array(naaom_cfg$M, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, biol(naaom_cfg$waa)),
    WAA_fish_input = replicate(1, array(rep(naaom_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(rep(naaom_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol(naaom_cfg$mat))))
  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list, R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(naaom_cfg$sigmaR), dim = c(2, 1, 1)),
    recruitment_opt = "mean_rec", init_age_strc = 1)

  # the piece under test
  sim_list <- Setup_Sim_NAA_state(sim_list, NAA_re = NAA_re, sigmaNAA = sigmaNAA,
                                  rho_age = rho_age, rho_year = rho_year)

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}

#' The operating model's observations, in the shape the estimation model reads
#'
#' Strips the replicate margin the simulator carries and peels to the terminal year.
#'
#' @param om Output of \code{naaom_make_om}.
#' @keywords internal
naaom_om_data <- function(om) simulation_data_to_SPoRC(sim_env = om, y = naaom_cfg$n_yrs, sim = 1)

#' Estimation model over that data, with the state on
#' @keywords internal
naaom_build_em <- function(sim_data, NAA_re = "none") {
  # read from the data rather than the fixture so the same builder serves an assessment run inside
  # a closed loop, where the series grows by one year at a time
  n_yrs <- dim(sim_data$WAA)[3]; n_ages <- naaom_cfg$n_ages
  il <- Setup_Mod_Dim(years = 1:n_yrs, ages = 1:n_ages, lens = NULL, n_regions = 1, n_sexes = 1,
                      n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1, natal_region = 1, verbose = FALSE)
  il <- Setup_Mod_Rec(input_list = il, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                      ln_sigmaR = array(log(naaom_cfg$sigmaR), c(2, 1, 1)), rec_model = "mean_rec",
                      sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2,
                      ln_global_R0 = log(5))
  il <- Setup_Mod_Biologicals(
    input_list = il, WAA = sim_data$WAA, MatAA = sim_data$MatAA, WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv, fit_lengths = 0, AgeingError = sim_data$AgeingError,
    M_spec = "fix", Fixed_natmort = array(naaom_cfg$M, dim = c(1, 1, n_yrs, n_ages, 1)),
    NAA_re = NAA_re)
  il <- Setup_Mod_Tagging(input_list = il, use_conv_fish_tagging = 0)
  il <- Setup_Mod_Movement(input_list = il, use_fixed_movement = 1, Fixed_Movement = NA,
                           do_recruits_move = 0)
  suppressWarnings(il <- Setup_Mod_Catch_and_F(
    input_list = il, ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch,
    Use_F_pen = 1, sigmaC_spec = "fix", ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(1, 1, 1)),
    ObsDiscard = sim_data$ObsDiscard, UseDiscard = sim_data$UseDiscard,
    sigma_dmr_spec = "fix", dmr_mean_spec = "fix", ln_sigmaD = sim_data$ln_sigmaD))
  il <- Setup_Mod_FishIdx_and_Comps(
    input_list = il, ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
    ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = NULL,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    UseFishLenComps = array(0, dim = dim(sim_data$UseFishAgeComps)),
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = NULL,
    fish_idx_type = "biom", FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none", FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1")
  il <- Setup_Mod_SrvIdx_and_Comps(
    input_list = il, ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx, SrvIdx_LikeType = "lognormal",
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = NULL,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    UseSrvLenComps = array(0, dim = dim(sim_data$UseSrvAgeComps)),
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = NULL,
    srv_idx_type = "biom", SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none", SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1")
  il <- Setup_Mod_Fishsel_and_Q(input_list = il, fish_sel_model = "logist1_Fleet_1",
                                fish_fixed_sel_pars_spec = "est_all", fish_q_spec = "est_all",
                                use_fixed_ret_sel = 1)
  il <- Setup_Mod_Srvsel_and_Q(input_list = il, srv_sel_model = "logist1_Fleet_1",
                               srv_fixed_sel_pars_spec = "est_all", srv_q_spec = "est_all")
  il <- Setup_Mod_Weighting(input_list = il, Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1,
                            Wt_Rec = 1, Wt_F = 1,
                            Wt_FishAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)),
                            Wt_SrvAgeComps = array(1, dim = c(1, n_yrs, 1, 1, 1)))
  il
}

#' Fixed effects of a fit, read by name
#'
#' Under random effects \code{last.par.best} holds the fixed and random parameters
#' together, and \code{parList} misassigns it: on this model it returns 3.32 for a
#' standard deviation whose value is 0.22, with only a warning about replacement
#' length to say so. Indexing by \code{env$random} to drop the states and then
#' selecting by name is what reads correctly.
#'
#' @param fit Object from \code{fit_model}.
#' @param nm Parameter name, or \code{NULL} for the whole fixed vector.
#' @keywords internal
naaom_fixed <- function(fit, nm = NULL) {
  p <- fit$env$last.par.best
  pf <- if(length(fit$env$random)) p[-fit$env$random] else p
  if(is.null(nm)) pf else unname(pf[names(pf) == nm])
}

#' Relative root mean squared error of estimated against true spawning biomass
#' @keywords internal
naaom_ssb_rmse <- function(fit, om) {
  truth <- as.vector(om$SSB[1, 1, , 1])
  sqrt(mean(((as.vector(fit$rep$SSB[1, 1, ]) - truth) / truth)^2))
}
