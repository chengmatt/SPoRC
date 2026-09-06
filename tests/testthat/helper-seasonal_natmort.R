# Small two season OM and the EM over its data, for the seasonal M tests. One pop, region and sex, so
# the only thing under test is what happens within the year.
#
# Seasons are unequal (0.4/0.6) on purpose, so confusing a rate with what it accumulates cannot pass.

seasonal_M_cfg <- list(n_yrs = 25, n_ages = 8, n_seas = 2, seasdur = c(0.4, 0.6), M = 0.3)

#' Two season simulated data set
#'
#' @param M_by_seas Instantaneous rate per year in each season. A single value is
#'   shared across both, which is the constant mortality the estimation model
#'   reproduces under \code{M_seasblk_spec = "constant"}.
#'
#' @keywords internal
seasonal_M_sim <- function(M_by_seas = seasonal_M_cfg$M, seed = 909, sigmaR = 0.3, idx_se = 0.1, iss = 300) {

  cfg <- seasonal_M_cfg
  n_yrs <- cfg$n_yrs; n_ages <- cfg$n_ages; n_seas <- cfg$n_seas
  if(length(M_by_seas) == 1) M_by_seas <- rep(M_by_seas, n_seas)

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

  # selectivity is flat over years and seasons, so it cannot absorb a seasonal
  # difference in mortality
  sel <- function(slope, infl) {
    array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))),
              each = n_yrs * n_seas),
          dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))
  }

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(1, sel(1.5, 3)),
    ret_sel_input = replicate(1, array(1, dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))),
    dmr_input = array(0, dim = c(1, n_yrs, n_seas, 1, 1)),
    Fmort_input = array(rep(seq(0.05, 0.25, length.out = n_yrs), times = n_seas),
                        dim = c(1, n_yrs, n_seas, 1, 1)),
    ISS_FishAgeComps = array(iss, dim = c(1, n_yrs, n_seas, 1, 1, 1)),
    ln_sigmaC = array(log(0.01), dim = c(1, n_yrs, n_seas, 1))
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(1, sel(1, 2.5)),
    ISS_SrvAgeComps = array(iss, dim = c(1, n_yrs, n_seas, 1, 1, 1)),
    ObsSrvIdx_SE = array(idx_se, dim = c(1, n_yrs, n_seas, 1))
  )

  # M has the season dim so the two seasons can differ
  natmort <- array(0, dim = c(1, 1, n_yrs, n_seas, n_ages, 1))
  for(seas in 1:n_seas) natmort[,,,seas,,] <- M_by_seas[seas]

  biol <- function(val) array(rep(val, each = n_yrs * n_seas),
                              dim = c(1, 1, n_yrs, n_seas, n_ages, 1))
  waa <- 5 / (1 + exp(-3 * ((1:n_ages) - 3)))
  mat <- 1 / (1 + exp(-3 * ((1:n_ages) - 3)))

  suppressWarnings(sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = replicate(1, natmort),
    WAA_input = replicate(1, biol(waa)),
    WAA_fish_input = replicate(1, array(rep(waa, each = n_yrs * n_seas),
                                        dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(rep(waa, each = n_yrs * n_seas),
                                       dim = c(1, 1, n_yrs, n_seas, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol(mat))
  ))

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, n_seas, n_ages, 1, 1))
  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = replicate(1, array(8, dim = c(1, 1, n_yrs))),
    ln_sigmaR = array(log(sigmaR), dim = c(2, 1, 1)),
    recruitment_opt = "mean_rec",
    init_age_strc = 1
  )

  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}


#' Estimation model over the two season data
#'
#' @param M_arg Named list passed straight to \code{Setup_Mod_Biologicals}, so a
#'   caller can hold mortality at a five or six dimensional array, or estimate it
#'   over season blocks.
#'
#' @keywords internal
seasonal_M_input <- function(sim_obj, M_arg, sigmaR = 0.3) {

  cfg <- seasonal_M_cfg
  n_yrs <- cfg$n_yrs; n_ages <- cfg$n_ages; n_seas <- cfg$n_seas
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = n_yrs, sim = 1)

  input_list <- Setup_Mod_Dim(
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

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(sigmaR), c(2, 1, 1)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(8)
  )

  input_list <- do.call(Setup_Mod_Biologicals, c(list(
    input_list = input_list,
    WAA = sim_data$WAA,
    MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv,
    AgeingError = sim_data$AgeingError
  ), M_arg))

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )

  input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sim_data$ObsCatch,
    Catch_Type = sim_data$Catch_Type,
    UseCatch = sim_data$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    ln_sigmaC = array(log(0.01), dim = c(1, n_yrs, n_seas, 1))
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = sim_data$UseFishIdx,
    ObsFishAgeComps = sim_data$ObsFishAgeComps,
    UseFishAgeComps = sim_data$UseFishAgeComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishLenComps = sim_data$ISS_FishLenComps,
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx,
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "agg_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = "logist1_Fleet_1",
    fish_fixed_sel_pars_spec = "est_all",
    fish_q_spec = "est_all"
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "est_all",
    srv_q_spec = "est_all"
  )

  fw <- array(1, dim = c(1, n_yrs, n_seas, 1, 1))
  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Discard = 1,
    Wt_D = 1,
    Wt_FishAgeComps = fw,
    Wt_FishLenComps = fw,
    Wt_SrvAgeComps = fw,
    Wt_SrvLenComps = fw
  )

  input_list
}

#' Fixed mortality array at a given per season rate
#' @keywords internal
seasonal_M_fixed <- function(M_by_seas, six_d = TRUE) {
  cfg <- seasonal_M_cfg
  if(length(M_by_seas) == 1) M_by_seas <- rep(M_by_seas, cfg$n_seas)
  if(!six_d) return(array(M_by_seas[1], dim = c(1, 1, cfg$n_yrs, cfg$n_ages, 1)))
  out <- array(0, dim = c(1, 1, cfg$n_yrs, cfg$n_seas, cfg$n_ages, 1))
  for(seas in 1:cfg$n_seas) out[,,,seas,,] <- M_by_seas[seas]
  out
}
