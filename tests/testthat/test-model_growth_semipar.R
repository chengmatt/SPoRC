# Time-varying and semi-parametric growth. One varies the growth PARAMETERS year by year, off each year's
# curve or cohort by cohort; the other adds a year-by-age surface of DEVIATIONS on mean length at age.
#
# Checks the mechanics of both against hand calculations, then that a deviation surface is recovered.

library(SPoRC)
library(testthat)

# a small test setup: one population, one region, one sex, a von Bertalanffy curve
gcfg <- list(n_yrs = 25, n_ages = 12, len_lower = seq(10, 75, by = 5))
gpars <- c(L1 = 14, L2 = 68, K = 0.25, CV1 = 0.12, CV2 = 0.07)

growth_call <- function(..., n_yrs = gcfg$n_yrs, n_ages = gcfg$n_ages) {
  ages <- 1:n_ages
  args <- list(
    ln_growth_pars = array(log(gpars), dim = c(1, 1, 1, 5)),
    growth_A1 = 1,
    growth_A2 = n_ages,
    growth_L0 = gcfg$len_lower[1],
    growth_len_lower = gcfg$len_lower,
    growth_cv_type = 0,
    growth_sd_type = 0,
    growth_dist = 0,
    growth_plus_group = 0,
    derive_waa = 1,
    wt_len_pars = array(c(1e-5, 3), dim = c(1, 1, 1, 2)),
    ages = ages,
    seasdur = 1,
    spawn_seas = 1,
    t_spawn = 0,
    n_pop = 1,
    n_regions = 1,
    n_yrs = n_yrs,
    n_seas = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    t_fish = array(0.5, dim = c(1, 1, 1)),
    t_srv = array(0.5, dim = c(1, 1, 1))
  )
  do.call(Get_Growth, utils::modifyList(args, list(...)))
}


test_that("a growth parameter's deviations move that parameter and no other", {

  n_yrs <- gcfg$n_yrs
  devs <- array(0, dim = c(1, 1, n_yrs, 5, 1))
  devs[1, 1, , 3, 1] <- seq(-0.3, 0.3, length.out = n_yrs) # K only

  # log link: the parameter is multiplied by exp(deviation)
  g <- growth_call(ln_growth_devs = devs, growth_tv_model = c(0, 0, 1, 0, 0), growth_tv_link = 0)
  expect_equal(g$growth_pars_y[1, 1, , 3, 1], gpars[["K"]] * exp(devs[1, 1, , 3, 1]), tolerance = 1e-12)
  for(k in c(1, 2, 4, 5)) expect_equal(unique(g$growth_pars_y[1, 1, , k, 1]), unname(gpars[k]), tolerance = 1e-12)
  # the von Bertalanffy form leaves the Richards coefficient at one
  expect_equal(unique(g$growth_pars_y[1, 1, , 6, 1]), 1)

  # logit link: the parameter stays inside its bounds however large the deviation,
  # approaching a bound rather than crossing it
  bnds <- matrix(c(1, 20, 40, 90, 0.05, 0.6, 0.01, 0.4, 0.01, 0.4), 5, 2, byrow = TRUE)
  big <- devs; big[1, 1, , 3, 1] <- seq(-40, 40, length.out = n_yrs)
  gl <- growth_call(
    ln_growth_devs = big,
    growth_tv_model = c(0, 0, 1, 0, 0),
    growth_tv_link = 1,
    growth_par_bounds = bnds
  )
  k_y <- gl$growth_pars_y[1, 1, , 3, 1]
  expect_true(all(k_y >= bnds[3, 1] & k_y <= bnds[3, 2]))
  expect_equal(k_y[1], bnds[3, 1], tolerance = 1e-6)              # driven to the lower bound
  expect_equal(k_y[n_yrs], bnds[3, 2], tolerance = 1e-6)          # and to the upper
  expect_true(all(diff(k_y) > 0))                                 # monotone in the deviation
  # a moderate deviation stays strictly inside
  mid <- devs; mid[1, 1, , 3, 1] <- seq(-2, 2, length.out = n_yrs)
  km <- growth_call(
    ln_growth_devs = mid,
    growth_tv_model = c(0, 0, 1, 0, 0),
    growth_tv_link = 1,
    growth_par_bounds = bnds
  )$growth_pars_y[1, 1, , 3, 1]
  expect_true(all(km > bnds[3, 1] & km < bnds[3, 2]))
  # and a zero deviation returns the parameter itself
  expect_equal(km[which.min(abs(mid[1, 1, , 3, 1]))], gpars[["K"]], tolerance = 1e-3)

  # every parameter can vary at once, each on its own series
  all_devs <- array(0, dim = c(1, 1, n_yrs, 5, 1))
  for(k in 1:5) all_devs[1, 1, , k, 1] <- 0.05 * k
  ga <- growth_call(ln_growth_devs = all_devs, growth_tv_model = rep(1, 5))
  for(k in 1:5) expect_equal(unique(ga$growth_pars_y[1, 1, , k, 1]), unname(gpars[k]) * exp(0.05 * k), tolerance = 1e-12)
})


test_that("the Richards coefficient generalizes the von Bertalanffy curve", {

  # rho of one reproduces the von Bertalanffy curve exactly
  vb <- growth_call()
  rich <- growth_call(ln_growth_pars = array(log(c(gpars, rho = 1)), dim = c(1, 1, 1, 6)))
  expect_equal(rich$mean_LAA_spawn[1, 1, 1, 1, , 1], vb$mean_LAA_spawn[1, 1, 1, 1, , 1], tolerance = 1e-10)

  # and a rho other than one bends it while keeping both reference lengths
  r2 <- growth_call(ln_growth_pars = array(log(c(gpars, rho = 1.6)), dim = c(1, 1, 1, 6)))
  L <- r2$mean_LAA_spawn[1, 1, 1, 1, , 1]
  expect_equal(L[1], gpars[["L1"]], tolerance = 1e-8)              # length at A1
  expect_equal(L[gcfg$n_ages], gpars[["L2"]], tolerance = 1e-8)    # length at A2
  expect_true(all(diff(L) > 0))                                    # still monotonic
  expect_false(isTRUE(all.equal(L, vb$mean_LAA_spawn[1, 1, 1, 1, , 1])))
})


test_that("semi-parametric deviations scale mean length at age and leave the CV alone", {

  n_yrs <- gcfg$n_yrs; n_ages <- gcfg$n_ages
  sp <- array(0, dim = c(1, 1, n_yrs, n_ages, 1))
  set.seed(42)
  sp[1, 1, , , 1] <- matrix(rnorm(n_yrs * n_ages, 0, 0.05), n_yrs, n_ages)

  base <- growth_call()
  g <- growth_call(ln_growth_semipar_devs = sp, growth_semipar = 5)

  # the mean at age is the parametric mean times exp(deviation), year by year
  for(y in c(1, 7, n_yrs)) {
    expect_equal(g$mean_LAA_spawn[1, 1, y, 1, , 1],
                 base$mean_LAA_spawn[1, 1, 1, 1, , 1] * exp(sp[1, 1, y, , 1]), tolerance = 1e-12)
    # the spread is that mean times the coefficient of variation the deviated
    # length implies, since the CV is interpolated ON LENGTH here: a deviation
    # that makes a fish longer moves it along the CV ramp as well
    L_y <- g$mean_LAA_spawn[1, 1, y, 1, , 1]
    cv_y <- gpars[["CV1"]] + (L_y - gpars[["L1"]]) * (gpars[["CV2"]] - gpars[["CV1"]]) / (gpars[["L2"]] - gpars[["L1"]])
    cv_y[n_ages] <- gpars[["CV2"]] # the oldest age is at the reference age, which takes CV2 outright
    expect_equal(g$sd_LAA_spawn[1, 1, y, 1, , 1] / L_y, cv_y, tolerance = 1e-6)
  }
  # under a CV that is a function of AGE instead, the deviations move only the
  # mean and leave the spread at age exactly where the parametric curve put it
  ga <- growth_call(ln_growth_semipar_devs = sp, growth_semipar = 5, growth_cv_type = 1)
  ba <- growth_call(growth_cv_type = 1)
  for(y in c(1, 7, n_yrs)) {
    expect_equal(ga$sd_LAA_spawn[1, 1, y, 1, , 1] / ga$mean_LAA_spawn[1, 1, y, 1, , 1],
                 ba$sd_LAA_spawn[1, 1, 1, 1, , 1] / ba$mean_LAA_spawn[1, 1, 1, 1, , 1], tolerance = 1e-12)
  }
  # a zero surface is the parametric model
  z <- growth_call(ln_growth_semipar_devs = array(0, dim = dim(sp)), growth_semipar = 5)
  expect_equal(z$mean_LAA_spawn, base$mean_LAA_spawn, tolerance = 1e-12)
  # and the key moves with the mean, so weight at age does too
  expect_false(isTRUE(all.equal(g$WAA[1, 1, 1, 1, , 1], g$WAA[1, 1, n_yrs, 1, , 1])))
})


test_that("cohort growth has size at age forward and blends the plus group by numbers", {

  n_yrs <- 6; n_ages <- gcfg$n_ages; ages <- 1:n_ages
  devs <- array(0, dim = c(1, 1, n_yrs, 5, 1))
  devs[1, 1, , 3, 1] <- c(0, 0, 0.2, -0.15, 0.1, 0) # K varies from year 3

  g <- growth_call(
    n_yrs = n_yrs,
    ln_growth_devs = devs,
    growth_tv_model = c(0, 0, 1, 0, 0),
    growth_tv_type = 1,
    growth_cohort_styr = 3
  )
  # only the years before the propagation starts are filled here; they all sit on
  # the first year's curve
  expect_equal(g$L_beg[1, 1, 1, , 1], g$L_beg[1, 1, 3, , 1], tolerance = 1e-12)

  # one year of propagation, done by hand: every propagated age grows from the
  # size it reached by the increment this year's parameters imply, and the plus
  # group blends the cohort entering it with the fish already there
  NAA_y <- array(exp(-0.3 * ages), dim = c(1, 1, n_ages, 1))
  g2 <- Get_Growth_Year(
    growth = g,
    y = 3,
    NAA_y = NAA_y,
    ln_growth_pars = array(log(gpars), dim = c(1, 1, 1, 5)),
    ln_growth_devs = devs,
    growth_tv_model = c(0, 0, 1, 0, 0),
    growth_tv_link = 0,
    growth_par_bounds = matrix(0, 5, 2),
    growth_A1 = 1,
    growth_A2 = n_ages,
    growth_L0 = gcfg$len_lower[1],
    growth_len_lower = gcfg$len_lower,
    growth_cv_type = 0,
    growth_sd_type = 0,
    growth_dist = 0,
    growth_plus_group = 0,
    derive_waa = 1,
    wt_len_pars = array(c(1e-5, 3), dim = c(1, 1, 1, 2)),
    ages = ages,
    seasdur = 1,
    spawn_seas = 1,
    t_spawn = 0,
    n_pop = 1,
    n_regions = 1,
    n_seas = 1,
    n_sexes = 1,
    t_fish = array(0.5, dim = c(1, 1, 1)),
    t_srv = array(0.5, dim = c(1, 1, 1))
  )

  # the asymptote is derived from the Schnute pair, so read the year's own
  K3 <- gpars[["K"]] * exp(devs[1, 1, 3, 3, 1]); Linf <- g2$Linf[1, 1, 3, 1]
  L_beg <- g$L_beg[1, 1, 3, , 1]
  grown <- grow_increment(L_beg, 1, K3, Linf, 1)
  nxt <- g2$L_beg[1, 1, 4, , 1]
  # ages past the first propagated one take the previous age's grown size
  for(a in 4:(n_ages - 1)) expect_equal(nxt[a], grown[a - 1], tolerance = 1e-10)
  # the plus group is the numbers-weighted blend
  blend <- ((NAA_y[1, 1, n_ages - 1, 1] + 0.01) * grown[n_ages - 1] +
              (NAA_y[1, 1, n_ages, 1] + 0.01) * grown[n_ages]) /
    (NAA_y[1, 1, n_ages - 1, 1] + NAA_y[1, 1, n_ages, 1] + 0.02)
  expect_equal(nxt[n_ages], blend, tolerance = 1e-10)
  # and it lies between the two sizes it blends
  expect_gte(nxt[n_ages], min(grown[n_ages - 1], grown[n_ages]) - 1e-10)
  expect_lte(nxt[n_ages], max(grown[n_ages - 1], grown[n_ages]) + 1e-10)
})


test_that("every semi-parametric process error form builds and is penalized", {

  source(test_path("helper-selftest_growth_semipar.R"), local = TRUE)
  for(form in c("iid", "rw", "2dar1", "3dmarg", "3dcond")) {
    inp <- semipar_input(form)
    expect_equal(inp$data$growth_semipar, c(iid = 1, rw = 2, `3dmarg` = 3, `3dcond` = 4, `2dar1` = 5)[[form]])
    fit <- fit_model(inp$data, inp$par, inp$map, random = NULL, silent = TRUE, do_optim = FALSE)
    # the deviations are estimated, the surface enters the objective, and the
    # process error is a finite penalty rather than a dropped term
    expect_equal(sum(!is.na(inp$map$ln_growth_semipar_devs)), length(inp$par$ln_growth_semipar_devs))
    expect_true(is.finite(fit$rep$growth_semipar_nLL))
    expect_true(all(is.finite(fit$gr(fit$par))))
    expect_jnLL_decomposes(fit, label = paste("semi-parametric growth,", form))
  }
})


test_that("a simulated deviation surface is recovered by the model that estimates it", {

  source(test_path("helper-selftest_growth_semipar.R"), local = TRUE)
  sim <- semipar_simulate(seed = 11)

  # the estimating model starts from a flat surface and has to find the one the
  # data were simulated under
  inp <- semipar_input("2dar1", obs = sim$obs)
  fit <- fit_model(
    inp$data,
    inp$par,
    inp$map,
    random = NULL,
    silent = TRUE,
    do_optim = TRUE,
    newton_loops = 2
  )
  est <- fit$rep$mean_LAA_srv[1, 1, , 1, , 1, 1]

  # mean length at age comes back close to the truth, and much closer than the
  # parametric curve alone could get
  flat <- fit_model(
    semipar_input("none", obs = sim$obs)$data,
    semipar_input("none", obs = sim$obs)$par,
    semipar_input("none", obs = sim$obs)$map,
    random = NULL,
    silent = TRUE,
    do_optim = TRUE,
    newton_loops = 2
  )
  flat_laa <- flat$rep$mean_LAA_srv[1, 1, , 1, , 1, 1]

  err_semi <- max(abs(est / sim$mean_LAA - 1))
  err_flat <- max(abs(flat_laa / sim$mean_LAA - 1))
  expect_lt(err_semi, 0.10)
  expect_lt(err_semi, err_flat)
})
