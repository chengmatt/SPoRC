library(SPoRC)
library(testthat)

# Reference point / projection consistency for a Ricker fit: Fmsy and Bmsy come
# from Get_Reference_Points' equilibrium solve, and a long deterministic
# projection at Fmsy must equilibrate at exactly Bmsy. The two go through
# different code paths (the analytic equilibrium against the year-loop
# projection), so agreement pins both. Routines shared through
# helper-selftest_features.R.

test_that("a 500-year projection at the Ricker Fmsy equilibrates at Bmsy", {

  n_yrs <- selftest_cfg$n_yrs; n_ages <- selftest_cfg$n_ages
  sigR <- 0.05
  om <- selftest_make_om(
    recruitment_opt = "ricker_rec",
    sigmaR = sigR,
    idx_se_om = 0.02,
    iss_om = 2000
  )
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = n_yrs, sim = 1)

  input_list <- selftest_build_input(sim_data, rec_model = "ricker_rec", sigmaR = sigR)
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL, silent = TRUE)
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-3)

  data <- fit$data; rep <- fit$rep
  rp <- Get_Reference_Points(
    data = data,
    rep = rep,
    type = "single_region",
    what = "MSY",
    calc_rec_st_yr = 1,
    rec_age = 1
  )

  n_proj_yrs <- 500
  constant_F_HCR <- function(x, frp, brp) frp
  terminal_NAA <- array(rep$NAA[, , n_yrs, , , ], dim = c(1, 1, 1, n_ages, 1))
  terminal_NAA0 <- array(rep$NAA0[, , n_yrs, , , ], dim = c(1, 1, 1, n_ages, 1))
  proj_arr <- function(vals_by_age) {
    out <- array(0, dim = c(1, 1, n_proj_yrs, 1, n_ages, 1))
    for(y in 1:n_proj_yrs) out[1, 1, y, 1, , 1] <- vals_by_age
    out
  }
  WAAp <- proj_arr(selftest_cfg$waa); MatAAp <- proj_arr(selftest_cfg$mat)
  WAA_fishp <- array(WAAp, dim = c(1, 1, n_proj_yrs, 1, n_ages, 1, 1))
  fish_selp <- array(0, dim = c(1, 1, n_proj_yrs, 1, n_ages, 1, 1))
  for(y in 1:n_proj_yrs) fish_selp[1, 1, y, 1, , 1, 1] <- rep$fish_sel[1, 1, n_yrs, 1, , 1, 1]
  natmortp <- array(0.3, dim = c(1, 1, n_proj_yrs, 1, n_ages, 1))
  Movementp <- array(1, dim = c(1, 1, 1, n_proj_yrs, 1, n_ages, 1))

  srr_opt <- list(
    rec_dd = 1,
    rec_lag = 1,
    do_recruits_move = 0,
    R0 = rep$R0,
    h = array(rep$h_trans, dim = c(1, 1)),
    rec_region_prop = rep$rec_region_prop,
    WAA = array(selftest_cfg$waa, dim = c(1, 1, 1, n_ages)),
    MatAA = array(selftest_cfg$mat, dim = c(1, 1, 1, n_ages)),
    SSB = rep$SSB,
    Movement = array(1, dim = c(1, 1, 1, 1, n_ages)),
    sex_ratio_f = array(0.5, dim = c(1, 1)),
    sgl_seas_spawning_movement = NULL,
    stray_rate = array(0, dim = c(1)),
    natmort = array(0.3, dim = c(1, 1, 1, n_ages)), # [pop, region, season, age]
    fish_sel = array(rep$fish_sel[1, 1, n_yrs, 1, , 1, 1], dim = c(1, 1, 1, n_ages, 1)),
    ret_sel = array(1, dim = c(1, 1, 1, n_ages, 1)),
    init_F = array(0, dim = c(1, 1, 1)),
    dmr = array(0, dim = c(1, 1, 1)),
    SR_ref_yr = if(!is.null(data$SR_ref_yr)) data$SR_ref_yr else 1
  )

  # sexratio apportions recruits to sexes, so a single-sex model takes 1 (the
  # 0.5 female spawning fraction is applied inside the SSB calculation,
  # matching the estimation model's convention)
  proj <- Do_Population_Projection(
    n_proj_yrs = n_proj_yrs,
    n_pop = 1,
    n_regions = 1,
    n_ages = n_ages,
    n_sexes = 1,
    sexratio = array(1, dim = c(1, 1, n_proj_yrs, 1)),
    n_fish_fleets = 1,
    do_recruits_move = 0,
    recruitment = rep$Rec,
    terminal_NAA = terminal_NAA,
    terminal_NAA0 = terminal_NAA0,
    terminal_F = array(rep$Fmort[1, n_yrs, 1, 1], dim = c(1, 1, 1)),
    natmort = natmortp,
    WAA = WAAp,
    WAA_fish = WAA_fishp,
    MatAA = MatAAp,
    fish_sel = fish_selp,
    ret_sel = array(1, dim = dim(fish_selp)),
    Movement = Movementp,
    f_ref_pt = array(rp$f_ref_pt, dim = c(1, n_proj_yrs)),
    b_ref_pt = array(rp$b_ref_pt, dim = c(1, 1, n_proj_yrs)),
    HCR_function = constant_F_HCR,
    recruitment_opt = "ricker_rec",
    fmort_opt = "HCR",
    t_spawn = 0,
    srr_opt = srr_opt
  )

  ssb_end <- proj$proj_SSB[1, 1, n_proj_yrs]
  expect_lt(abs((ssb_end - rp$b_ref_pt) / rp$b_ref_pt), 1e-6)
  # and it has actually equilibrated, rather than still settling through Bmsy
  expect_equal(ssb_end, proj$proj_SSB[1, 1, n_proj_yrs - 1], tolerance = 1e-10)

})
