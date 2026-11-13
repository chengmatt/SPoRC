# Growth options no other test checks: seasons, several strata, the plus group choice, cohort growth
# with no deviations, and the parameter, sharing and year maps Setup_Mod_Biologicals builds.

library(SPoRC)
library(testthat)

opt_cfg <- list(n_yrs = 8, n_ages = 10, len_lower = seq(10, 80, by = 5))
opt_pars <- c(L1 = 15, L2 = 70, K = 0.25, CV1 = 0.1, CV2 = 0.08)

# the growth module alone, one stratum and one season unless told otherwise
opt_growth_call <- function(...) {
  args <- list(
    ln_growth_pars = array(log(opt_pars), dim = c(1, 1, 1, 5)),
    growth_A1 = 1,
    growth_A2 = opt_cfg$n_ages,
    growth_L0 = opt_cfg$len_lower[1],
    growth_len_lower = opt_cfg$len_lower,
    growth_cv_type = 0,
    growth_sd_type = 0,
    growth_dist = 0,
    growth_plus_group = 1,
    derive_waa = 1,
    wt_len_pars = array(c(1e-5, 3), dim = c(1, 1, 1, 2)),
    ages = 1:opt_cfg$n_ages,
    seasdur = 1,
    spawn_seas = 1,
    t_spawn = 0,
    n_pop = 1,
    n_regions = 1,
    n_yrs = opt_cfg$n_yrs,
    n_seas = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    t_srv = array(0.5, dim = c(1, 1, 1))
  )
  do.call(Get_Growth, utils::modifyList(args, list(...)))
}

# the setup with two regions and two sexes, stopping after the biologicals
opt_setup_call <- function(...) {
  n_yrs <- opt_cfg$n_yrs
  n_ages <- opt_cfg$n_ages
  input <- Setup_Mod_Dim(
    years = 1:n_yrs,
    ages = 1:n_ages,
    lens = opt_cfg$len_lower + 2.5,
    n_regions = 2,
    n_sexes = 2,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_pop = 1,
    natal_region = 1,
    verbose = FALSE
  )
  input <- Setup_Mod_Rec(
    input_list = input,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    ln_sigmaR = array(log(0.3), c(2, 1, 2)),
    rec_model = "mean_rec",
    sigmaR_spec = "fix",
    init_age_strc = 1,
    equil_init_age_strc = 2,
    ln_global_R0 = log(5)
  )
  args <- list(
    input_list = input,
    WAA = NULL,
    MatAA = array(1, c(1, 2, n_yrs, 1, n_ages, 2)),
    fit_lengths = 1,
    SizeAgeTrans = NA,
    AgeingError = diag(n_ages),
    M_spec = "fix",
    Fixed_natmort = array(0.3, dim = c(1, 2, n_yrs, n_ages, 2)),
    growth_model = "vb_schnute",
    ln_growth_pars = array(log(opt_pars), dim = c(1, 2, 2, 5)),
    growth_A1 = 1,
    growth_A2 = n_ages,
    growth_len_lower = opt_cfg$len_lower,
    waa_model = "wt_len",
    wt_len_pars = c(1e-5, 3)
  )
  do.call(Setup_Mod_Biologicals, utils::modifyList(args, list(...)))
}

# map levels laid out on the parameter's own dims
map_levels <- function(input, par_name) array(as.integer(input$map[[par_name]]), dim = dim(input$par[[par_name]]))

# The Growth Module ----------------------------------------------------------

test_that("splitting the year into seasons leaves size at a given point in the year unchanged", {
  n_ages <- opt_cfg$n_ages

  # two seasons of 0.4 and 0.6 of a year: the second season's mid point and a survey half way
  # through it both fall at 0.7 of the year, and spawning half way through it too
  g2 <- opt_growth_call(
    n_seas = 2,
    seasdur = c(0.4, 0.6),
    spawn_seas = 2,
    t_spawn = 0.5,
    t_srv = array(0.5, dim = c(1, 2, 1))
  )
  # one season read at 0.7 of the year
  g1 <- opt_growth_call(t_spawn = 0.7, t_srv = array(0.7, dim = c(1, 1, 1)))
  expect_equal(g2$SizeAgeTrans_srv[1, 1, , 2, , , 1, 1], g1$SizeAgeTrans_srv[1, 1, , 1, , , 1, 1], tolerance = 1e-12)
  expect_equal(g2$SizeAgeTrans_fish[1, 1, , 2, , , 1, 1], g1$SizeAgeTrans_srv[1, 1, , 1, , , 1, 1], tolerance = 1e-12)
  expect_equal(g2$SizeAgeTrans_spawn, g1$SizeAgeTrans_spawn, tolerance = 1e-12)
  expect_equal(g2$WAA[1, 1, , 2, , 1], g1$WAA[1, 1, , 1, , 1], tolerance = 1e-12)

  # the first season's mid point is 0.2 of the year, plus group included
  g_02 <- opt_growth_call(t_srv = array(0.2, dim = c(1, 1, 1)))
  expect_equal(g2$mean_LAA_fish[1, 1, , 1, , 1, 1], g_02$mean_LAA_srv[1, 1, , 1, , 1, 1], tolerance = 1e-12)
  expect_gt(g2$mean_LAA_fish[1, 1, 1, 2, n_ages, 1, 1], g2$mean_LAA_fish[1, 1, 1, 1, n_ages, 1, 1])
})

test_that("each population, region and sex is built from its own parameters", {
  ln_pars <- array(0, dim = c(1, 2, 2, 5))
  ln_pars[1, 1, 1, ] <- log(opt_pars)
  ln_pars[1, 2, 1, ] <- log(opt_pars * c(1.1, 1, 1, 1, 1)) # region two, longer at A1
  ln_pars[1, 1, 2, ] <- log(opt_pars * c(1, 0.9, 1.2, 1, 1)) # sex two, smaller and faster
  ln_pars[1, 2, 2, ] <- log(opt_pars * c(1, 1, 1, 1.5, 0.8)) # both, a different spread
  wt_len <- array(c(1e-5, 1e-5, 1.2e-5, 1.2e-5, 3, 3, 2.9, 2.9), dim = c(1, 2, 2, 2))

  g <- opt_growth_call(
    ln_growth_pars = ln_pars,
    wt_len_pars = wt_len,
    n_regions = 2,
    n_sexes = 2,
    t_srv = array(0.5, dim = c(2, 1, 1))
  )
  for(r in 1:2) {
    for(s in 1:2) {
      g_rs <- opt_growth_call(
        ln_growth_pars = array(ln_pars[1, r, s, ], dim = c(1, 1, 1, 5)),
        wt_len_pars = array(wt_len[1, r, s, ], dim = c(1, 1, 1, 2))
      )
      expect_equal(g$SizeAgeTrans_srv[1, r, , , , , s, 1], g_rs$SizeAgeTrans_srv[1, 1, , , , , 1, 1], tolerance = 1e-12)
      expect_equal(g$WAA[1, r, , , , s], g_rs$WAA[1, 1, , , , 1], tolerance = 1e-12)
    } # end s loop
  } # end r loop
})

test_that("the plus group reads the curve at the accumulator age or the decay-weighted mixture", {
  n_ages <- opt_cfg$n_ages
  on_curve <- opt_growth_call(growth_plus_group = 0)
  mixture <- opt_growth_call(growth_plus_group = 1)
  crv <- get_laa_curve(1:n_ages, opt_cfg$len_lower[1], opt_pars[["L1"]], opt_pars[["L2"]], opt_pars[["K"]],
                       opt_pars[["CV1"]], opt_pars[["CV2"]], A1 = 1, A2 = n_ages)

  expect_equal(on_curve$L_beg[1, 1, 1, , 1], crv$L, tolerance = 1e-12)
  expect_equal(mixture$L_beg[1, 1, 1, n_ages, 1], plus_group_size(crv$L[n_ages], crv$Linf, n_ages), tolerance = 1e-12)
  # every younger age is the same under both
  expect_equal(mixture$L_beg[1, 1, 1, -n_ages, 1], crv$L[-n_ages], tolerance = 1e-12)
})

test_that("cohort growth with no deviations stays on the curve below the plus group", {
  n_yrs <- opt_cfg$n_yrs
  n_ages <- opt_cfg$n_ages
  zero_devs <- array(0, dim = c(1, 1, n_yrs, 5, 1))
  tv_model <- c(0, 0, 1, 0, 0) # K varies, by zero

  curve_growth <- opt_growth_call()
  cohort_growth <- opt_growth_call(ln_growth_devs = zero_devs, growth_tv_model = tv_model, growth_tv_type = 1, growth_cohort_styr = 1)
  NAA_y <- array(exp(-0.3 * (1:n_ages)), dim = c(1, 1, n_ages, 1))
  for(y in 1:n_yrs) {
    cohort_growth <- Get_Growth_Year(
      growth = cohort_growth,
      y = y,
      NAA_y = NAA_y,
      ln_growth_pars = array(log(opt_pars), dim = c(1, 1, 1, 5)),
      ln_growth_devs = zero_devs,
      growth_tv_model = tv_model,
      growth_tv_link = 0,
      growth_par_bounds = matrix(0, 5, 2),
      growth_A1 = 1,
      growth_A2 = n_ages,
      growth_L0 = opt_cfg$len_lower[1],
      growth_len_lower = opt_cfg$len_lower,
      growth_cv_type = 0,
      growth_sd_type = 0,
      growth_dist = 0,
      growth_plus_group = 1,
      derive_waa = 1,
      wt_len_pars = array(c(1e-5, 3), dim = c(1, 1, 1, 2)),
      ages = 1:n_ages,
      seasdur = 1,
      spawn_seas = 1,
      t_spawn = 0,
      n_pop = 1,
      n_regions = 1,
      n_seas = 1,
      n_sexes = 1,
      t_srv = array(0.5, dim = c(1, 1, 1))
    )
  } # end y loop

  below_plus <- 1:(n_ages - 1)
  expect_equal(cohort_growth$L_beg[1, 1, , below_plus, 1], curve_growth$L_beg[1, 1, , below_plus, 1], tolerance = 1e-12)
  expect_equal(cohort_growth$mean_LAA_srv[1, 1, , 1, below_plus, 1, 1], curve_growth$mean_LAA_srv[1, 1, , 1, below_plus, 1, 1], tolerance = 1e-12)
  expect_equal(cohort_growth$sd_LAA_srv[1, 1, , 1, below_plus, 1, 1], curve_growth$sd_LAA_srv[1, 1, , 1, below_plus, 1, 1], tolerance = 1e-12)
})

# The Setup ------------------------------------------------------------------

test_that("growth_spec shares the growth parameters over the strata it names", {
  est_all <- map_levels(opt_setup_call(), "ln_growth_pars") # [pop, region, sex, parameter]
  expect_equal(length(unique(as.vector(est_all))), 2 * 2 * 5)

  shared_r <- map_levels(opt_setup_call(growth_spec = "est_shared_r"), "ln_growth_pars")
  expect_equal(shared_r[1, 1, , ], shared_r[1, 2, , ])
  expect_equal(length(unique(as.vector(shared_r))), 2 * 5)

  shared_s <- map_levels(opt_setup_call(growth_spec = "est_shared_s"), "ln_growth_pars")
  expect_equal(shared_s[1, , 1, ], shared_s[1, , 2, ])
  expect_equal(length(unique(as.vector(shared_s))), 2 * 5)

  shared_r_s <- map_levels(opt_setup_call(growth_spec = "est_shared_r_s"), "ln_growth_pars")
  for(k in 1:5) expect_equal(length(unique(as.vector(shared_r_s[1, , , k]))), 1)
  expect_equal(length(unique(as.vector(shared_r_s))), 5)

  expect_true(all(is.na(map_levels(opt_setup_call(growth_spec = "fix"), "ln_growth_pars"))))
})

test_that("growth_fix keeps the named parameters at their starting values whatever growth_spec says", {
  fix_K <- map_levels(opt_setup_call(growth_fix = c(FALSE, FALSE, TRUE, FALSE, FALSE)), "ln_growth_pars")
  expect_true(all(is.na(fix_K[, , , 3])))
  expect_false(any(is.na(fix_K[, , , -3])))

  fix_K_shared <- map_levels(opt_setup_call(growth_spec = "est_shared_r_s", growth_fix = c(FALSE, FALSE, TRUE, FALSE, FALSE)), "ln_growth_pars")
  expect_true(all(is.na(fix_K_shared[, , , 3])))
  expect_equal(length(unique(as.vector(fix_K_shared[, , , -3]))), 4)
})

test_that("growth_A2 = 'Linf' reads L2 as the asymptote", {
  input <- opt_setup_call(growth_A2 = "Linf")
  expect_equal(input$data$growth_L2_asymptote, 1)
  expect_equal(input$data$growth_A2, opt_cfg$n_ages) # the accumulator age bounds the spread
  g <- opt_growth_call(growth_L2_asymptote = 1)
  expect_equal(g$Linf[1, 1, 1, 1], opt_pars[["L2"]])
  expect_lt(g$L_beg[1, 1, 1, opt_cfg$n_ages - 1, 1], opt_pars[["L2"]])
})

test_that("time-varying deviations are estimated only in growth_tv_years and shared as growth_tv_spec says", {
  input <- opt_setup_call(
    growth_tv_model = c(K = "iid"),
    growth_tv_years = 3:6,
    growth_tv_spec = "est_shared_r",
    growth_tv_sigma_spec = "est"
  )
  expect_equal(input$data$growth_tv_model, c(0, 0, 1, 0, 0))

  # deviations [pop, region, year, parameter, sex]: K only, years 3 to 6, regions shared, sexes not
  dev_map <- map_levels(input, "ln_growth_devs")
  expect_true(all(is.na(dev_map[, , , -3, ])))
  expect_true(all(is.na(dev_map[, , -(3:6), 3, ])))
  expect_false(any(is.na(dev_map[, , 3:6, 3, ])))
  expect_equal(dev_map[1, 1, , 3, ], dev_map[1, 2, , 3, ])
  expect_equal(length(unique(as.vector(dev_map[!is.na(dev_map)]))), 4 * 2)

  # the process error sd, slot three of the time-varying half, shared the same way
  pe_map <- map_levels(input, "growth_pe_pars") # [pop, region, slot, sex, half]
  expect_false(any(is.na(pe_map[1, , 3, , 1])))
  expect_true(all(is.na(pe_map[, , -3, , ])))
  expect_true(all(is.na(pe_map[, , , , 2])))
  expect_equal(pe_map[1, 1, 3, , 1], pe_map[1, 2, 3, , 1])
  expect_true(all(is.na(map_levels(opt_setup_call(growth_tv_model = c(K = "iid")), "growth_pe_pars")))) # "fix" by default

  # a random walk on the same parameter
  rw <- opt_setup_call(growth_tv_model = c(K = "rw"))
  expect_equal(rw$data$growth_tv_model, c(0, 0, 2, 0, 0))
})
