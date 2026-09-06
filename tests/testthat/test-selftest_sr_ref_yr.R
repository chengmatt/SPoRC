library(SPoRC)
library(testthat)

# Every input to unfished spawning biomass per recruit is taken at SR_ref_yr, in both the
# estimation model and the operating model. R0 is the deliberate exception: it scales phi0
# into S0 and also sits in the curve's numerator, so it is the year's own value.
# Helpers, including a time-varying weight at age that makes the reference year bite, are
# in helper-sr_ref_yr.R.

test_that("phi0 through Get_Det_Recruitment matches its definition at any reference year", {

  # recover phi0 as the SSB at which deterministic recruitment equals R0, then compare with
  # the per-recruit sum written out by hand
  phi0_sporc <- function(yr, R0 = 10) {
    n_ages <- sr_ref_cfg$n_ages
    args <- list(
      recruitment_model = 1,
      rec_dd = 0,
      y = 5,
      rec_lag = 1,
      R0 = R0,
      rec_region_prop = array(1, dim = c(1, 1)),
      rec_seas_prop = array(1, dim = c(1, 1)),
      h = array(sr_ref_cfg$h, dim = c(1, 1)),
      n_pop = 1,
      n_regions = 1,
      n_ages = n_ages,
      n_fish_fleets = 1,
      WAA = array(sr_ref_waa(yr), dim = c(1, 1, 1, n_ages)),
      MatAA = array(sr_ref_mat(), dim = c(1, 1, 1, n_ages)),
      natmort = array(sr_ref_cfg$M, dim = c(1, 1, 1, n_ages)),
      Movement = array(1, dim = c(1, 1, 1, 1, n_ages)),
      sgl_seas_spawning_movement = array(1, dim = c(1, 1, 1, n_ages)),
      stray_rate = array(0, dim = 1),
      do_recruits_move = 0,
      t_spawn = 0,
      init_F = array(0, dim = c(1, 1, 1)),
      dmr = array(0, dim = c(1, 1, 1)),
      fish_sel = array(1, dim = c(1, 1, 1, n_ages, 1)),
      ret_sel = array(1, dim = c(1, 1, 1, n_ages, 1)),
      n_seas = 1,
      spawn_seas = 1,
      natal_region = 1,
      seasdur = 1,
      sexratio_f = array(0.5, dim = c(1, 1))
    )
    f <- function(sv) SPoRC:::Get_Det_Recruitment(
      recruitment_model = args$recruitment_model,
      rec_dd = args$rec_dd,
      y = args$y,
      rec_lag = args$rec_lag,
      R0 = args$R0,
      rec_region_prop = args$rec_region_prop,
      rec_seas_prop = args$rec_seas_prop,
      h = args$h,
      n_pop = args$n_pop,
      n_regions = args$n_regions,
      n_ages = args$n_ages,
      n_fish_fleets = args$n_fish_fleets,
      WAA = args$WAA,
      MatAA = args$MatAA,
      natmort = args$natmort,
      SSB_vals = array(sv, dim = c(1, 1, 10)),
      Movement = args$Movement,
      sgl_seas_spawning_movement = args$sgl_seas_spawning_movement,
      stray_rate = args$stray_rate,
      do_recruits_move = args$do_recruits_move,
      t_spawn = args$t_spawn,
      init_F = args$init_F,
      dmr = args$dmr,
      fish_sel = args$fish_sel,
      ret_sel = args$ret_sel,
      n_seas = args$n_seas,
      spawn_seas = args$spawn_seas,
      natal_region = args$natal_region,
      seasdur = args$seasdur,
      sexratio_f = args$sexratio_f
    )[1, 1] - R0
    stats::uniroot(f, c(1e-8, 1e6))$root / R0
  }

  for(yr in c(1, sr_ref_cfg$n_yrs)) expect_equal(phi0_sporc(yr), sr_ref_phi0(yr), tolerance = 1e-6)
  # the constructed weight at age doubles across the series, so phi0 must too
  expect_equal(sr_ref_phi0(sr_ref_cfg$n_yrs) / sr_ref_phi0(1), 2, tolerance = 1e-8)
})

test_that("S0 is rebuilt from the R0 handed in, so it is not proportional at fixed SSB", {

  # if S0 were kept independent of R0, doubling R0 would exactly double recruitment at a
  # fixed SSB. It does not, because S0 moves with R0.
  n_ages <- sr_ref_cfg$n_ages
  rec_at <- function(R0) SPoRC:::Get_Det_Recruitment(
    recruitment_model = 1,
    rec_dd = 0,
    y = 5,
    rec_lag = 1,
    R0 = R0,
    rec_region_prop = array(1, dim = c(1, 1)),
    rec_seas_prop = array(1, dim = c(1, 1)),
    h = array(sr_ref_cfg$h, dim = c(1, 1)),
    n_pop = 1,
    n_regions = 1,
    n_ages = n_ages,
    n_fish_fleets = 1,
    WAA = array(sr_ref_waa(1), dim = c(1, 1, 1, n_ages)),
    MatAA = array(sr_ref_mat(), dim = c(1, 1, 1, n_ages)),
    natmort = array(sr_ref_cfg$M, dim = c(1, 1, 1, n_ages)),
    SSB_vals = array(20, dim = c(1, 1, 10)),
    Movement = array(1, dim = c(1, 1, 1, 1, n_ages)),
    sgl_seas_spawning_movement = array(1, dim = c(1, 1, 1, n_ages)),
    stray_rate = array(0, dim = 1),
    do_recruits_move = 0,
    t_spawn = 0,
    init_F = array(0, dim = c(1, 1, 1)),
    dmr = array(0, dim = c(1, 1, 1)),
    fish_sel = array(1, dim = c(1, 1, 1, n_ages, 1)),
    ret_sel = array(1, dim = c(1, 1, 1, n_ages, 1)),
    n_seas = 1,
    spawn_seas = 1,
    natal_region = 1,
    seasdur = 1,
    sexratio_f = array(0.5, dim = c(1, 1))
  )[1, 1]
  expect_false(isTRUE(all.equal(rec_at(12) / rec_at(5), 12 / 5, tolerance = 1e-6)))
})

test_that("the reference year changes the operating model", {

  skip_on_cran()
  a <- sr_ref_make_om(SR_ref_yr = 1)
  b <- sr_ref_make_om(SR_ref_yr = sr_ref_cfg$n_yrs)
  expect_false(isTRUE(all.equal(as.vector(a$Rec), as.vector(b$Rec))))
  # doubling phi0 halves depletion, so the terminal recruitment differs substantially
  expect_gt(abs(tail(as.vector(b$Rec), 1) / tail(as.vector(a$Rec), 1) - 1), 0.5)
})

test_that("a matched reference year recovers R0 and a mismatched one biases it", {

  skip_on_cran()
  n <- sr_ref_cfg$n_yrs; true_R0 <- 10
  om <- sr_ref_make_om(SR_ref_yr = n)
  sd <- simulation_data_to_SPoRC(sim_env = om, y = n, sim = 1)

  fit_at <- function(ref) {
    il <- sr_ref_make_em(sd, SR_ref_yr = ref)
    f <- fit_model(il$data, il$par, il$map, silent = TRUE, do_optim = TRUE, n_newton_loops = 1)
    list(R0 = f$rep$R0[1], grad = max(abs(f$gr(f$env$last.par.best))))
  }
  matched <- fit_at(n)
  mismatched <- fit_at(1)

  # R0 is essentially a mean over n_yrs log-recruitments, so a single replicate has a
  # standard error of sigmaR / sqrt(n_yrs), about 5.5% here. Any absolute tolerance has to
  # be set from that rather than picked: at three standard errors a correct model passes
  # this on all but roughly one draw in 400. Measured over ten seeds the mean error is
  # within noise of zero and its sign flips, so there is no systematic bias to catch.
  tol <- 3 * sr_ref_cfg$sigmaR / sqrt(n)
  expect_lt(matched$grad, 1e-3)
  expect_lt(abs(matched$R0 / true_R0 - 1), tol)

  # the load-bearing assertion is PAIRED: both fits see the same simulated data and differ
  # only in SR_ref_yr, so the replicate's own noise cancels and what is left is the effect
  # of the reference year. The wrong year misses by far more than the noise it sits in.
  expect_lt(abs(matched$R0 / true_R0 - 1), abs(mismatched$R0 / true_R0 - 1))
  expect_gt(abs(mismatched$R0 / true_R0 - 1) - abs(matched$R0 / true_R0 - 1), tol)
})
