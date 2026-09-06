# A small operating and estimation model pair for the block tests. Survey selectivity breaks at a known
# year, and optionally catchability steps at the same year, so a blocked process has a truth to recover.

blocks_sim <- function(n_sims = 1, n_yrs = 30, true_break = 11, seed = 123, q_late = 1) {
  n_ages <- 10
  sl <- Setup_Sim_Dim(
    n_sims = n_sims,
    n_yrs = n_yrs,
    n_regions = 1,
    n_ages = n_ages,
    n_lens = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_seas = 1,
    n_pop = 1
  )
  sl <- Setup_Sim_Containers(sl)
  logi <- function(k, a50) 1 / (1 + exp(-k * ((1:n_ages) - a50)))
  arr7 <- function(v, ny, nf) array(rep(v, each = ny), dim = c(1, 1, ny, 1, n_ages, 1, nf))
  sl <- Setup_Sim_Fishing(sl, fish_sel_input = replicate(n_sims, arr7(logi(3, 5), n_yrs, 1)))
  q_path <- c(rep(1, true_break - 1), rep(q_late, n_yrs - true_break + 1))
  sl <- Setup_Sim_Survey(sl, srv_sel_input = replicate(n_sims, arr7(logi(1, 3), n_yrs, 1)),
                         srv_q_input = array(rep(q_path, times = n_sims), dim = c(1, n_yrs, 1, n_sims)))
  sl$srv_sel[, , true_break:n_yrs, , , , , ] <- replicate(n_sims, arr7(logi(3, 5), length(true_break:n_yrs), 1))
  waa <- 5 * logi(3, 3); mat <- logi(3, 3)
  arr6 <- function(v) array(rep(v, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  sl <- Setup_Sim_Biologicals(
    sl,
    natmort_input = replicate(n_sims, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(n_sims, arr6(waa)),
    WAA_fish_input = replicate(n_sims, arr7(waa, n_yrs, 1)),
    WAA_srv_input = replicate(n_sims, arr7(waa, n_yrs, 1)),
    MatAA_input = replicate(n_sims, arr6(mat))
  )
  sl <- Setup_Sim_Tagging(sl, use_conv_fish_tagging = 0)
  sl$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, n_sims))
  sl <- Setup_Sim_Rec(
    sl,
    R0_input = replicate(n_sims, array(5, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(1), dim = c(2, 1, 1)),
    recruitment_opt = 'mean_rec',
    init_age_strc = 1
  )
  set.seed(seed)
  suppressWarnings(suppressMessages(Simulate_Pop_Static(sim_list = sl, output_path = NULL)))
}

blocks_em <- function(
  sim_obj,
  sim = 1,
  sel_blocks = "none_Fleet_1",
  q_blocks = "none_Fleet_1",
  fish_q_blocks = "none_Fleet_1",
  fish_sel_blocks = "none_Fleet_1"
) suppressWarnings({
  sd <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = sim)
  il <- Setup_Mod_Dim(
    years = 1:sim_obj$n_years,
    ages = 1:sim_obj$n_ages,
    lens = sim_obj$n_lens,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    natal_region = sim_obj$natal_region,
    verbose = FALSE
  )
  il <- Setup_Mod_Rec(
    il,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(1), c(2, 1, 1)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(5)
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
    Fixed_natmort = array(0.3, dim = c(1, 1, sim_obj$n_years, sim_obj$n_ages, 1))
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
    UseFishIdx = sd$UseFishIdx,
    ObsFishAgeComps = sd$ObsFishAgeComps,
    ObsFishLenComps = sd$ObsFishLenComps,
    UseFishAgeComps = sd$UseFishAgeComps,
    UseFishLenComps = sd$UseFishLenComps,
    ISS_FishAgeComps = sd$ISS_FishAgeComps,
    ISS_FishLenComps = sd$ISS_FishLenComps,
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
    ObsSrvAgeComps = sd$ObsSrvAgeComps,
    ObsSrvLenComps = sd$ObsSrvLenComps,
    UseSrvAgeComps = sd$UseSrvAgeComps,
    UseSrvLenComps = sd$UseSrvLenComps,
    ISS_SrvAgeComps = sd$ISS_SrvAgeComps,
    ISS_SrvLenComps = sd$ISS_SrvLenComps,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  il <- Setup_Mod_Fishsel_and_Q(
    il,
    fish_sel_model = "logist2_Fleet_1",
    fish_q_blocks = fish_q_blocks,
    fish_sel_blocks = fish_sel_blocks,
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "est_all"
  )
  il <- Setup_Mod_Srvsel_and_Q(
    il,
    srv_sel_blocks = sel_blocks,
    srv_q_blocks = q_blocks,
    srv_sel_model = "logist2_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all"
  )
  ny <- sim_obj$n_years
  Setup_Mod_Weighting(
    il,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(1, ny, 1, 1, 1)),
    Wt_FishLenComps = array(1, dim = c(1, ny, 1, 1, 1)),
    Wt_SrvAgeComps = array(1, dim = c(1, ny, 1, 1, 1)),
    Wt_SrvLenComps = array(0, dim = c(1, ny, 1, 1, 1))
  )
})
