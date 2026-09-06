# The model against the equations it is documented by.
#
# c_model_equations.Rmd states what SPoRC computes and nothing checks that it does. A regression test fixes
# the number the code produces, so a model computing something else passes forever.
#
# helper-equation_handwrite.R holds a second implementation written from the vignette. Where the two
# disagree, either the code has diverged from its documentation or the documentation is wrong.

oracle_setup <- function(t_spawn = 0, n_yrs = 10, n_ages = 6) {
  il <- collapse_input(nr = 1, nx = 1, nf = 1)
  if(t_spawn != 0) il$data$t_spawn <- t_spawn
  fit <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  list(il = il, rep = fit$rep, n_yrs = collapse_cfg$n_yrs, n_ages = collapse_cfg$n_ages)
}


test_that("total mortality is natural mortality plus fishing mortality at age", {
  fx <- oracle_setup()
  Z <- oracle_ya(fx$rep$ZAA, fx$n_yrs, fx$n_ages)
  M <- oracle_ya(fx$rep$natmort, fx$n_yrs, fx$n_ages)
  FAA <- oracle_ya(fx$rep$tot_FAA, fx$n_yrs, fx$n_ages)

  expect_equal(Z, M + FAA, tolerance = 1e-12)
})


test_that("numbers at age follow the documented projection and plus group", {
  # N[y+1, a+1] = N[y, a] exp(-Z[y, a]), with the plus group taking survivors of
  # both the last true age and of itself.
  fx <- oracle_setup()
  Z <- oracle_ya(fx$rep$ZAA, fx$n_yrs, fx$n_ages)
  N <- drop(fx$rep$NAA)

  reference <- oracle_project_naa(n1 = N[1, ], rec = oracle_y(fx$rep$Rec, fx$n_yrs), Z = Z)

  expect_equal(N[2:(fx$n_yrs + 1), ], reference[2:(fx$n_yrs + 1), ], tolerance = 1e-10)
})


test_that("the plus group conserves the fish entering it", {
  # Stated separately from the projection because it is the line most often
  # written as a plain assignment rather than an accumulation, and a model that
  # drops the second term still projects every other age correctly.
  fx <- oracle_setup()
  Z <- oracle_ya(fx$rep$ZAA, fx$n_yrs, fx$n_ages)
  N <- drop(fx$rep$NAA)
  A <- fx$n_ages

  for(y in seq_len(fx$n_yrs)) {
    entering <- N[y, A - 1] * exp(-Z[y, A - 1])
    staying <- N[y, A] * exp(-Z[y, A])
    expect_equal(N[y + 1, A], entering + staying, tolerance = 1e-10,
                 label = sprintf("plus group in year %d", y))
  }
})


test_that("catch at age follows Baranov's equation", {
  fx <- oracle_setup()
  Z <- oracle_ya(fx$rep$ZAA, fx$n_yrs, fx$n_ages)
  N <- drop(fx$rep$NAA)[seq_len(fx$n_yrs), ]
  retF <- oracle_ya(fx$rep$ret_FAA, fx$n_yrs, fx$n_ages)
  CAA <- oracle_ya(fx$rep$CAA, fx$n_yrs, fx$n_ages)

  expect_equal(CAA, oracle_baranov(retF, Z, N), tolerance = 1e-10)
})


test_that("spawning biomass weighs the population propagated to spawning", {
  # t_spawn = 0 leaves the population where it is; a non-zero value has it
  # into the season first. Both are checked, because a model that ignores the
  # timing agrees with the reference at zero and only differs away from it.
  for(t_spawn in c(0, 0.35)) {
    fx <- oracle_setup(t_spawn = t_spawn)
    Z <- oracle_ya(fx$rep$ZAA, fx$n_yrs, fx$n_ages)
    N <- drop(fx$rep$NAA)[seq_len(fx$n_yrs), ]
    WAA <- oracle_ya(fx$il$data$WAA, fx$n_yrs, fx$n_ages)
    MatAA <- oracle_ya(fx$il$data$MatAA, fx$n_yrs, fx$n_ages)

    expect_equal(oracle_y(fx$rep$SSB, fx$n_yrs),
                 oracle_ssb(N, Z, WAA, MatAA, t_spawn), tolerance = 1e-10,
                 label = sprintf("SSB at t_spawn = %g", t_spawn))
    expect_equal(oracle_y(fx$rep$Total_Biom, fx$n_yrs),
                 oracle_total_biomass(N, Z, WAA, t_spawn), tolerance = 1e-10,
                 label = sprintf("total biomass at t_spawn = %g", t_spawn))
  }
})


test_that("the spawning timing check is sensitive to the timing", {
  # If SSB ignored t_spawn the test above would still pass at t_spawn = 0 and
  # against a reference that also ignored it. The two values have to differ.
  a <- oracle_y(oracle_setup(t_spawn = 0)$rep$SSB, collapse_cfg$n_yrs)
  b <- oracle_y(oracle_setup(t_spawn = 0.35)$rep$SSB, collapse_cfg$n_yrs)

  expect_gt(max(abs(a - b) / pmax(abs(a), 1e-30)), 0.01)
})


test_that("logistic selectivity matches the documented form", {
  # Sel_b = 1 / (1 + exp(-k (b - b50))), evaluated through the model's own
  # selectivity routine so the comparison is against what the model computes
  # rather than against a second copy of the same call. Parameters are (b50,
  # slope), which the vignette does not say.
  bins <- 1:12
  for(pars in list(c(4, 1.2), c(7, 0.4), c(2, 3))) {
    got <- Get_Selex(
      Selex_Model = 0,
      Bin = bins,
      pars = log(pars),
      TimeVary_Model = 0,
      ln_seldevs = array(0, c(1, 1, 6, 1, 1)),
      Region = 1,
      Year = 1,
      Sex = 1
    )
    expect_equal(as.numeric(got), oracle_selex("logist1", bins, pars), tolerance = 1e-10,
                 label = sprintf("logist1 at b50=%g slope=%g", pars[1], pars[2]))
  }
})


test_that("gamma selectivity squares the peak inside the root", {
  # The vignette prints sqrt(bmax + 4 delta^2); the code and this reference both
  # use sqrt(bmax^2 + 4 delta^2). Pinning it here says which one the model means,
  # so the vignette can be corrected against a test rather than against a reading.
  bins <- 1:12
  for(pars in list(c(5, 2), c(8, 1.5))) {
    got <- Get_Selex(
      Selex_Model = 1,
      Bin = bins,
      pars = log(pars),
      TimeVary_Model = 0,
      ln_seldevs = array(0, c(1, 1, 6, 1, 1)),
      Region = 1,
      Year = 1,
      Sex = 1
    )
    expect_equal(as.numeric(got), oracle_selex("gamma", bins, pars), tolerance = 1e-10,
                 label = sprintf("gamma at bmax=%g delta=%g", pars[1], pars[2]))

    # the unsquared form the vignette prints is a different curve, so the test
    # above is choosing between them rather than passing either way
    unsquared <- {
      p <- 0.5 * (sqrt(pars[1] + 4 * pars[2]^2) - pars[1])
      (bins / pars[1])^(pars[1] / p) * exp((pars[1] - bins) / p)
    }
    expect_gt(max(abs(as.numeric(got) - unsquared)), 1e-6)
  }
})


test_that("the multinomial composition likelihood matches the documented form", {
  # -l = ESS sum_b (O_b + c) [log(O_b + c) - log(E_b + c)]
  #
  # The sweep test setup is used rather than the collapse one because it fits
  # compositions; the collapse test setup switches every data source off.
  n_yrs <- 8; n_ages <- 5
  il <- sweep_input(dims = list(
    n_regions = 1,
    n_sexes = 1,
    n_fish_fleets = 1,
    n_srv_fleets = 1,
    n_yrs = n_yrs,
    n_ages = n_ages
  ))
  # observed proportions that are not flat, so a likelihood ignoring either
  # argument is distinguishable from one reading both
  set.seed(42)
  obs <- array(0, dim = c(1, n_yrs, 1, n_ages, 1, 1))
  for(y in seq_len(n_yrs)) obs[1, y, 1, , 1, 1] <- {
    w <- stats::runif(n_ages, 0.5, 2); w / sum(w)
  }
  il$data$ObsFishAgeComps <- obs
  rep <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep

  # the documented likelihood is stated on the expected proportions, with ageing
  # error a separate step before it. Holding that step at the identity keeps this
  # a test of the likelihood rather than of the two composed.
  n_yrs_ae <- dim(il$data$AgeingError)[1]
  for(y in seq_len(n_yrs_ae)) {
    expect_equal(il$data$AgeingError[y, , ], diag(n_ages), tolerance = 1e-12,
                 label = sprintf("ageing error is the identity in year %d", y))
  }

  reported <- drop(rep$FishAgeComps_nLL)
  # expected proportions come from catch at age, which is what the compositions
  # are predicted from; PredCatchAA is the separate at-age data source
  pred <- drop(rep$CAA)
  iss <- drop(il$data$ISS_FishAgeComps)
  const <- il$data$addtocomp
  const_obs <- isTRUE(il$data$comp_const_obs == 1)

  expect_false(all(reported == 0))

  for(y in seq_len(n_yrs)) {
    o <- obs[1, y, 1, , 1, 1] / sum(obs[1, y, 1, , 1, 1])
    p <- pred[y, ] / sum(pred[y, ])
    expect_equal(as.numeric(reported[y]),
                 oracle_multinomial_nll(o, p, iss[y], const, const_obs),
                 tolerance = 1e-8, label = sprintf("multinomial nLL in year %d", y))
  }
})


test_that("the Dirichlet-multinomial effective sample size matches the documented form", {
  # ESS = 1 / (1 + theta) + ISS theta / (1 + theta)
  for(theta in c(0.1, 1, 7.5)) {
    for(iss in c(25, 100, 400)) {
      expect_equal(oracle_dm_ess(iss, theta),
                   1 / (1 + theta) + iss * theta / (1 + theta), tolerance = 1e-12)
    }
  }
  # the identity worth stating: a large theta recovers the input sample size, and
  # a vanishing one collapses to a single observation
  expect_equal(oracle_dm_ess(300, 1e6), 300, tolerance = 1e-3)
  expect_equal(oracle_dm_ess(300, 1e-9), 1, tolerance = 1e-5)
})
