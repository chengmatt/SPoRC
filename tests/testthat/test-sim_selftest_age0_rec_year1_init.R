# Get_Init_NAA already seeds the recruit age class in year one, so recruiting at age 0 must not
# add recruits on top of it. The contamination only shows years later as that cohort matures.

library(SPoRC)
library(testthat)

test_that("Age-0 (rec_lag = 0) recruitment does not double-count Get_Init_NAA's equilibrium seed in year 1", {

  # Get_Init_NAA seeds the recruit age class with an equilibrium value built from R0 whatever
  # the recruitment lag. Recruiting at age 1 overwrites that seed before the season loop.
  #
  # Recruiting at age 0 instead inserts recruitment inside the season loop, so it has to
  # replace that seed rather than add to it, which would double year one's recruitment.
  #
  # Year one's own reported recruitment and spawning biomass would not show it, since the
  # extra age 0 fish do not mature for years, so this reads the numbers at age directly.

  n_pop <- 1
  n_regions <- 1
  n_ages <- 8
  n_sexes <- 1
  n_fish_fleets <- 1
  n_srv_fleets <- 1
  n_seas <- 1
  n_years <- 10
  true_R0 <- 10
  true_h <- 0.7

  inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

  input_list <- Setup_Mod_Dim(
    years = 1:n_years,
    ages = 1:n_ages,
    lens = NULL,
    n_regions = n_regions,
    n_sexes = n_sexes,
    n_fish_fleets = n_fish_fleets,
    n_srv_fleets = n_srv_fleets,
    n_seas = n_seas,
    n_pop = n_pop,
    verbose = FALSE
  )

  input_list <- Setup_Mod_Rec(input_list = input_list, rec_model = "bh_rec", rec_dd = "global",
    rec_lag = 0, spawn_seas = 1, t_spawn = 0, do_rec_bias_ramp = 0, sigmaR_switch = 1,
    ln_sigmaR = array(log(0.4), c(2, n_pop, n_regions)), sigmaR_spec = "fix",
    h_spec = "fix", steepness_h = array(inv_steepness(true_h), dim = c(n_pop, n_regions)),
    init_age_strc = 1, equil_init_age_strc = 0, ln_global_R0 = log(true_R0),
    RecDevs_spec = "fix" # fixed at 0 - fully deterministic, nothing to mask the bug
  )

  MatAA <- array(0, dim = c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes))
  MatAA[,,,,3:n_ages,] <- 1 # ages at index 1-2 immature; required MatAA == 0 at recruit age
  WAA <- array(rep(5 / (1 + exp(-3 * ((1:n_ages) - 3))), each = n_pop * n_regions * n_years * n_seas),
              dim = c(n_pop, n_regions, n_years, n_seas, n_ages, n_sexes))

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = WAA,
    MatAA = MatAA,
    fit_lengths = 0,
    M_spec = "fix",
    Fixed_natmort = array(0.2, dim = c(n_pop, n_regions, n_years, n_ages, n_sexes))
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
    ObsCatch = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
    UseCatch = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
    Use_F_pen = 0,
    sigmaC_spec = "fix"
  ))

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
    ObsFishIdx_SE = array(0.2, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
    UseFishIdx = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
    ObsFishAgeComps = array(0, dim = c(n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets)),
    UseFishAgeComps = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
    ISS_FishAgeComps = array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)),
    ObsFishLenComps = array(0, dim = c(n_regions, n_years, n_seas, 1, n_sexes, n_fish_fleets)),
    UseFishLenComps = array(0, dim = c(n_regions, n_years, n_seas, n_fish_fleets)),
    ISS_FishLenComps = array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)),
    fish_idx_type = c("none"),
    FishAgeComps_LikeType = c("none"),
    FishLenComps_LikeType = c("none"),
    FishAgeComps_Type = c("none_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = array(0, dim = c(n_regions, n_years, n_seas, n_srv_fleets)),
    ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years, n_seas, n_srv_fleets)),
    UseSrvIdx = array(0, dim = c(n_regions, n_years, n_seas, n_srv_fleets)),
    ObsSrvAgeComps = array(0, dim = c(n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets)),
    UseSrvAgeComps = array(0, dim = c(n_regions, n_years, n_seas, n_srv_fleets)),
    ISS_SrvAgeComps = array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets)),
    ObsSrvLenComps = array(0, dim = c(n_regions, n_years, n_seas, 1, n_sexes, n_srv_fleets)),
    UseSrvLenComps = array(0, dim = c(n_regions, n_years, n_seas, n_srv_fleets)),
    ISS_SrvLenComps = array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets)),
    srv_idx_type = c("none"),
    SrvAgeComps_LikeType = c("none"),
    SrvLenComps_LikeType = c("none"),
    SrvAgeComps_Type = c("none_Year_1-terminal_Fleet_1"),
    SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = c("logist2_Fleet_1"),
    fish_fixed_sel_pars_spec = c("fix"),
    fish_q_spec = "fix"
  )
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = c("logist2_Fleet_1"),
    srv_fixed_sel_pars_spec = c("fix"),
    srv_q_spec = c("fix")
  )

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)),
    Wt_FishLenComps = array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets)),
    Wt_SrvLenComps = array(0, dim = c(n_regions, n_years, n_seas, n_sexes, n_srv_fleets))
  )

  model <- suppressWarnings(fit_model(
    input_list$data,
    input_list$par,
    input_list$map,
    random = NULL,
    silent = TRUE,
    do_fit = FALSE
  ))
  rep <- model$rep

  # with the recruitment deviations at zero, year one's recruit-age numbers have to be exactly
  # recruitment times the season share times the sex ratio, and nothing more.
  #
  # this deliberately does not assert that recruitment and spawning biomass are constant across
  # years: the equilibrium age structure is an approximation, and both lags settle the same way
  spawn_seas <- input_list$data$spawn_seas
  expected_age0_naa <- rep$Rec[1,1,1] * rep$rec_seas_prop[1,spawn_seas] * 1 # sexratio = 1 (n_sexes = 1)
  expect_equal(rep$NAA[1,1,1,spawn_seas,1,1], expected_age0_naa, tolerance = 1e-8)
  expect_equal(rep$NAA0[1,1,1,spawn_seas,1,1], expected_age0_naa, tolerance = 1e-8)

})
