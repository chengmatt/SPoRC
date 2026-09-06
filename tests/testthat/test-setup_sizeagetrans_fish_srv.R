# SizeAgeTrans_fish/SizeAgeTrans_srv: fixed per-fleet size-age keys, the growth_model = "none"
# counterpart of derived per-fleet keys and of WAA_fish/WAA_srv overriding the shared WAA.
#
# Before this, every fleet under growth_model = "none" was forced onto the one shared SizeAgeTrans.

library(SPoRC)
library(testthat)

n_yrs <- 4; n_ages <- 6; n_lens <- 8
lens <- seq(10, 45, by = 5)

# two hand-built keys: fleet 1 mass concentrated in the short bins, fleet 2 in
# the long ones, both column-stochastic (age columns sum to one)
key_short <- matrix(0, n_lens, n_ages)
key_long <- matrix(0, n_lens, n_ages)
for(a in 1:n_ages) {
  key_short[, a] <- dnorm(1:n_lens, mean = 2, sd = 1.2)
  key_long[, a] <- dnorm(1:n_lens, mean = 7, sd = 1.2)
  key_short[, a] <- key_short[, a] / sum(key_short[, a])
  key_long[, a] <- key_long[, a] / sum(key_long[, a])
}
shared_key <- (key_short + key_long) / 2 # what both fleets get without an override

mk_sat <- function(key) array(rep(key, n_yrs), dim = c(1, 1, n_yrs, 1, n_lens, n_ages, 1))

build_input <- function(sat_fish = NULL) {

  input_list <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = lens,
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 2,
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
  input_list <- Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = array(1e-4 * (lens[1:n_ages + 1])^3, dim = c(1, 1, n_yrs, 1, n_ages, 1)),
    MatAA = array(rep(c(0, 0.5, 1, 1, 1, 1), n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1)),
    fit_lengths = 1,
    SizeAgeTrans = mk_sat(shared_key),
    SizeAgeTrans_fish = sat_fish,
    AgeingError = NULL,
    M_spec = "fix",
    Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))
  )
  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    use_fixed_movement = 1,
    Fixed_Movement = NA,
    do_recruits_move = 0
  )
  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- suppressWarnings(Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = array(100, dim = c(1, n_yrs, 1, 2)),
    UseCatch = array(1, dim = c(1, n_yrs, 1, 2)),
    Use_F_pen = 0,
    ln_F_mean_spec = "est",
    sigmaC_spec = "fix",
    ln_sigmaC = array(log(0.05), dim = c(1, n_yrs, 1, 2))
  ))
  no_age <- array(0, dim = c(1, n_yrs, 1, 2))
  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = array(NA, dim = c(1, n_yrs, 1, 2)),
    ObsFishIdx_SE = array(NA, dim = c(1, n_yrs, 1, 2)),
    UseFishIdx = array(0, dim = c(1, n_yrs, 1, 2)),
    fish_idx_type = rep("none", 2),
    ObsFishAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, 2)),
    UseFishAgeComps = no_age,
    FishAgeComps_LikeType = rep("none", 2),
    FishAgeComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:2),
    ObsFishLenComps = array(0, dim = c(1, n_yrs, 1, n_lens, 1, 2)),
    UseFishLenComps = array(0, dim = c(1, n_yrs, 1, 2)),
    FishLenComps_LikeType = rep("none", 2),
    FishLenComps_Type = paste0("none_Year_1-terminal_Fleet_", 1:2)
  )
  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = array(NA, dim = c(1, n_yrs, 1, 1)),
    ObsSrvIdx_SE = array(NA, dim = c(1, n_yrs, 1, 1)),
    UseSrvIdx = array(0, dim = c(1, n_yrs, 1, 1)),
    srv_idx_type = "none",
    ObsSrvAgeComps = array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1)),
    UseSrvAgeComps = array(0, dim = c(1, n_yrs, 1, 1)),
    SrvAgeComps_LikeType = "none",
    SrvAgeComps_Type = "none_Year_1-terminal_Fleet_1",
    ObsSrvLenComps = array(0, dim = c(1, n_yrs, 1, n_lens, 1, 1)),
    UseSrvLenComps = array(0, dim = c(1, n_yrs, 1, 1)),
    SrvLenComps_LikeType = "none",
    SrvLenComps_Type = "none_Year_1-terminal_Fleet_1"
  )
  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_selex_type = "length",
    fish_sel_model = paste0("logist1_Fleet_", 1:2),
    fish_fixed_sel_pars_spec = rep("fix", 2),
    fish_q_spec = rep("fix", 2)
  )
  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = "logist1_Fleet_1",
    srv_fixed_sel_pars_spec = "fix",
    srv_q_spec = "fix"
  )
  input_list <- Setup_Mod_Weighting(input_list = input_list, Wt_Catch = 1, Wt_F = 1)

  # a flat logistic (b50 far below every bin) makes selectivity ~1 everywhere,
  # so the predicted catch at length is essentially the key itself. Dims are
  # [region, param, block, sex, fleet]; param 1 is b50, param 2 the slope
  input_list$par$fish_fixed_sel_pars[, 1, , , ] <- -5
  input_list$par$fish_fixed_sel_pars[, 2, , , ] <- 5
  input_list$par$srv_fixed_sel_pars[, 1, , , ] <- -5
  input_list$par$srv_fixed_sel_pars[, 2, , , ] <- 5
  input_list$par$ln_F_mean[] <- log(0.3)

  input_list
}


test_that("SizeAgeTrans_fish requires growth_model = 'none' and the right dimensions", {

  base <- build_input()

  expect_error(
    Setup_Mod_Biologicals(
      input_list = Setup_Mod_Rec(
        Setup_Mod_Dim(
          years = 1:n_yrs,
          ages = 1:n_ages,
          lens = lens,
          n_regions = 1,
          n_sexes = 1,
          n_fish_fleets = 2,
          n_srv_fleets = 1,
          n_pop = 1,
          natal_region = 1,
          verbose = FALSE
        ),
        do_rec_bias_ramp = 0,
        sigmaR_switch = 1,
        ln_sigmaR = array(log(0.3), c(2, 1, 1)),
        rec_model = "mean_rec",
        sigmaR_spec = "fix",
        init_age_strc = 1,
        equil_init_age_strc = 2,
        ln_global_R0 = log(8)
      ),
      WAA = array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1)),
      MatAA = array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1)),
      fit_lengths = 1,
      growth_model = "vb_schnute",
      ln_growth_pars = array(log(c(20, 60, 0.2, 0.1, 0.1)), c(1, 1, 1, 5)),
      growth_A1 = 1,
      growth_A2 = n_ages,
      growth_len_lower = lens,
      SizeAgeTrans = NA,
      SizeAgeTrans_fish = mk_sat(key_short),
      AgeingError = NULL,
      M_spec = "fix",
      Fixed_natmort = array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))
    ),
    "growth_model = 'none'"
  )

  ok_shape <- array(0, dim = c(1, 1, n_yrs, 1, n_lens, n_ages, 1, 2))
  ok_shape[, , , , , , , 1] <- mk_sat(key_short); ok_shape[, , , , , , , 2] <- mk_sat(key_long)
  bad_dim <- ok_shape[, , 1:2, , , , , , drop = FALSE] # wrong year count
  expect_error(build_input(sat_fish = bad_dim), "SizeAgeTrans_fish")
})


test_that("each fleet's catch at length follows its own fixed key", {

  fish_keys <- array(0, dim = c(1, 1, n_yrs, 1, n_lens, n_ages, 1, 2))
  fish_keys[, , , , , , , 1] <- mk_sat(key_short)
  fish_keys[, , , , , , , 2] <- mk_sat(key_long)

  il <- build_input(sat_fish = fish_keys)
  fit <- suppressWarnings(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))

  cal1 <- apply(fit$rep$CAL[1, 1, , 1, , 1, 1], 2, sum)
  cal2 <- apply(fit$rep$CAL[1, 1, , 1, , 1, 2], 2, sum)

  # fleet 1's catch at length peaks in the short bins its key concentrates on,
  # fleet 2's in the long bins, and they are not the same shape
  expect_gt(which.max(cal2), which.max(cal1))
  expect_false(isTRUE(all.equal(cal1 / sum(cal1), cal2 / sum(cal2))))

  # the fixed keys are echoed back in the report, like a supplied WAA_fish would be
  expect_equal(fit$rep$SizeAgeTrans_fish, fish_keys)

  # without an override both fleets share the same SizeAgeTrans, so they no
  # longer differ this way
  il0 <- build_input(sat_fish = NULL)
  fit0 <- suppressWarnings(fit_model(il0$data, il0$par, il0$map, do_optim = FALSE, silent = TRUE))
  cal1_0 <- apply(fit0$rep$CAL[1, 1, , 1, , 1, 1], 2, sum)
  cal2_0 <- apply(fit0$rep$CAL[1, 1, , 1, , 1, 2], 2, sum)
  expect_equal(cal1_0 / sum(cal1_0), cal2_0 / sum(cal2_0), tolerance = 1e-8)
  expect_null(fit0$rep$SizeAgeTrans_fish)
})
