library(SPoRC)
library(testthat)

# The bin restriction has to reach the use flags through the real setup entry
# points, not just through drop_empty_fitted_blocks called directly. Two of the
# three setup paths were wired in the wrong order and the unit tests could not
# see it: the discard data sources reconciled before their bins were parsed, and the
# survey data sources reconciled after their use arrays had already been stored. Both
# were silent no-ops. This drives the restriction end to end instead.

test_that("a restriction reaches the use flags through every setup path", {
  skip_if_not(exists("objective_setup_sim"), "helper-objective_setup.R not loaded")

  sim_obj <- objective_setup_sim(NULL)
  sd <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)
  n_yrs <- sim_obj$n_years; n_ages <- sim_obj$n_ages
  keep <- 1:4
  bad_year <- 5

  # year 5 has mass only outside the bins being fitted, on every data source
  blank_year <- function(arr) {
    arr[1, bad_year, 1, , 1, 1] <- 0
    arr[1, bad_year, 1, 5:n_ages, 1, 1] <- 20
    arr
  }
  sd$ObsFishAgeComps <- blank_year(sd$ObsFishAgeComps)
  sd$ObsSrvAgeComps <- blank_year(sd$ObsSrvAgeComps)
  sd$ObsFishAgeComps_discard[1, , 1, , 1, 1] <- sd$ObsFishAgeComps[1, , 1, , 1, 1]
  sd$UseFishAgeComps_discard[] <- 1
  sd$ISS_FishAgeComps_discard <- sd$ISS_FishAgeComps

  input_list <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = sim_obj$n_lens,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    natal_region = sim_obj$natal_region,
    verbose = FALSE
  )
  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(0.5), c(2, 1, 1)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(5)
  )
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = sd$WAA,
    MatAA = sd$MatAA,
    WAA_fish = sd$WAA_fish,
    WAA_srv = sd$WAA_srv,
    fit_lengths = 0,
    SizeAgeTrans = NA,
    AgeingError = sd$AgeingError,
    M_spec = "fix",
    Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))
  )
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )
  suppressWarnings(input_list <- Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sd$ObsCatch,
    UseCatch = sd$UseCatch,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    ln_sigmaC = sd$ln_sigmaC,
    ln_sigmaF = array(log(1), dim = c(1, 1, 1)),
    ObsDiscard = sd$ObsDiscard,
    UseDiscard = sd$UseDiscard,
    sigma_dmr_spec = "fix",
    dmr_mean_spec = "est_all",
    ln_sigmaD = sd$ln_sigmaD
  ))

  suppressWarnings(input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sd$ObsFishIdx,
    ObsFishIdx_SE = sd$ObsFishIdx_SE,
    UseFishIdx = sd$UseFishIdx,
    ObsFishAgeComps = sd$ObsFishAgeComps,
    UseFishAgeComps = sd$UseFishAgeComps,
    ISS_FishAgeComps = sd$ISS_FishAgeComps,
    ObsFishLenComps = sd$ObsFishLenComps,
    UseFishLenComps = sd$UseFishLenComps,
    ISS_FishLenComps = sd$ISS_FishLenComps,
    fish_idx_type = "biom",
    FishAgeComps_LikeType = "Multinomial",
    FishLenComps_LikeType = "none",
    FishAgeComps_Type = "spltRspltS_Year_1-terminal_Fleet_1",
    FishLenComps_Type = "none_Year_1-terminal_Fleet_1",
    FishAgeComps_bins = list(keep),
    ObsFishAgeComps_discard = sd$ObsFishAgeComps_discard,
    UseFishAgeComps_discard = sd$UseFishAgeComps_discard,
    ISS_FishAgeComps_discard = sd$ISS_FishAgeComps_discard,
    FishAgeComps_discard_LikeType = "Multinomial",
    FishAgeComps_discard_Type = "spltRspltS_Year_1-terminal_Fleet_1",
    FishAgeComps_discard_bins = list(keep)
  ))

  suppressWarnings(input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sd$ObsSrvIdx,
    ObsSrvIdx_SE = sd$ObsSrvIdx_SE,
    UseSrvIdx = sd$UseSrvIdx,
    ObsSrvAgeComps = sd$ObsSrvAgeComps,
    UseSrvAgeComps = sd$UseSrvAgeComps,
    ISS_SrvAgeComps = sd$ISS_SrvAgeComps,
    ObsSrvLenComps = sd$ObsSrvLenComps,
    UseSrvLenComps = sd$UseSrvLenComps,
    ISS_SrvLenComps = sd$ISS_SrvLenComps,
    srv_idx_type = "biom",
    SrvAgeComps_LikeType = "Multinomial",
    SrvLenComps_LikeType = "none",
    SrvAgeComps_Type = "spltRspltS_Year_1-terminal_Fleet_1",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1",
    SrvAgeComps_bins = list(keep)
  ))

  # every restriction took
  expect_equal(which(input_list$data$FishAgeComps_bins[,1] == 1), keep)
  expect_equal(which(input_list$data$FishAgeComps_discard_bins[,1] == 1), keep)
  expect_equal(which(input_list$data$SrvAgeComps_bins[,1] == 1), keep)

  # and every one of the three setup paths cleared the emptied year
  expect_equal(input_list$data$UseFishAgeComps[1, bad_year, 1, 1], 0,
               info = "retained fishery")
  expect_equal(input_list$data$UseFishAgeComps_discard[1, bad_year, 1, 1], 0,
               info = "discard fishery: reconciled before its bins were parsed")
  expect_equal(input_list$data$UseSrvAgeComps[1, bad_year, 1, 1], 0,
               info = "survey: reconciled after its use array was stored")

  # exactly the years holding no mass in the fitted bins are cleared, no more and
  # no fewer. The planted year 5 is one of them; the test setup supplies others of
  # its own, which is the reconciliation doing real work rather than a no-op.
  empty_years <- function(obs) which(apply(obs[1, , 1, keep, 1, 1], 1, sum) == 0)
  expect_equal(which(input_list$data$UseFishAgeComps[1, , 1, 1] == 0),
               empty_years(sd$ObsFishAgeComps))
  expect_equal(which(input_list$data$UseFishAgeComps_discard[1, , 1, 1] == 0),
               empty_years(sd$ObsFishAgeComps_discard))
  expect_equal(which(input_list$data$UseSrvAgeComps[1, , 1, 1] == 0),
               empty_years(sd$ObsSrvAgeComps))
  expect_true(bad_year %in% empty_years(sd$ObsFishAgeComps))
})
