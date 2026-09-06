library(SPoRC)
library(testthat)
library(Matrix)

# Operating-model / estimation-model parity under move_timing.
#
# A simulation self-test can only recover the truth if the simulator and the estimation
# model implement the same dynamics. The operator-level properties are covered in
# test-transition-operators.R and test-move-timing-selftest.R; what this file guards is
# that the *simulator* honors move_timing everywhere the estimation model does.
#
# Three modules previously defaulted to move_timing = 0 inside the simulator while the
# estimation model ran at the configured timing, which biased the round trip:
#   1. compute_biom_y_sim  (spawning-time biomass)
#   2. Get_Init_NAA        (initial age structure)
#   3. Get_Det_Recruitment (SSB0 / SPR behind Beverton-Holt)

# Shared test setup: a three-region CTMC generator with region-varying mortality, which is
# the only regime where the three timings actually differ.
make_setup <- function(n_regions = 3, n_ages = 8, n_seas = 1, n_sexes = 1, n_yrs = 3) {
  seasdur <- rep(1 / n_seas, n_seas)
  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L
  dat <- expand.grid(
    pop = 1,
    regions = seq_len(n_regions),
    years = seq_len(n_yrs),
    seas = seq_len(n_seas),
    ages = seq_len(n_ages),
    sexes = seq_len(n_sexes)
  )

  mv <- Get_Movement(
    move_type = 1,
    do_recruits_move = 0,
    n_pop = 1,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_proj_yrs_devs = 0,
    n_ages = n_ages,
    n_sexes = n_sexes,
    n_seas = n_seas,
    move_pars = array(0, c(1, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
    use_fixed_movement = 0,
    Fixed_Movement = NULL,
    ctmc_move_dat = dat,
    preference_formula = ~ 1,
    diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.35),
    move_preference_pars = 0,
    area_r = rep(1, n_regions),
    adjacency_mat = adj,
    ctmc_diffusion_bounds = 0,
    seasdur = seasdur,
    ctmc_scale_by_seasdur = 1
  )

  list(n_regions = n_regions, n_ages = n_ages, n_seas = n_seas, n_sexes = n_sexes,
       n_yrs = n_yrs, seasdur = seasdur, adj = adj, ctmc_dat = dat,
       Movement = mv$Movement, Mrate = mv$Mrate,
       # region-varying F, so movement interacts with mortality
       Fr = c(0.02, 0.12, 0.30)[seq_len(n_regions)], M = 0.25)
}

# ---------------------------------------------------------------------------
# 1. Spawning-time biomass
# ---------------------------------------------------------------------------

test_that("simulator and estimation model agree on spawning biomass under every timing", {
  # compute_biom_y_sim() is the simulator's copy of compute_biom_y(). Given identical
  # state they must return identical SSB, Total_Biom and Dynamic_SSB0 -- otherwise the
  # operating model's recruitment is driven off a different SSB than the estimation
  # model reconstructs, and no amount of refitting recovers the truth.
  set.seed(11)
  fx <- make_setup()
  n_regions <- fx$n_regions; n_ages <- fx$n_ages; n_sexes <- fx$n_sexes
  n_pop <- 1; n_seas <- fx$n_seas; n_yrs <- fx$n_yrs
  y <- 2; seas <- 1; t_spawn <- 0.3

  NAA <- array(stats::runif(n_pop * n_regions * (n_yrs + 1) * n_seas * n_ages * n_sexes, 20, 200),
               dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes))
  NAA0 <- NAA * 1.4
  WAA <- array(rep(5 / (1 + exp(-3 * (seq_len(n_ages) - 3))), each = n_pop * n_regions * (n_yrs + 1) * n_seas),
               dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes))
  MatAA <- array(rep(1 / (1 + exp(-3 * (seq_len(n_ages) - 3))), each = n_pop * n_regions * (n_yrs + 1) * n_seas),
                 dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes))
  natmort <- array(fx$M, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes))
  ZAA <- array(0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes))
  for (r in seq_len(n_regions)) ZAA[1, r, , , , ] <- (fx$M + fx$Fr[r]) * fx$seasdur[1]
  sgl <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs + 1, n_ages, n_sexes))
  for (yy in seq_len(n_yrs + 1)) for (a in seq_len(n_ages)) for (s in seq_len(n_sexes)) {
    sgl[1, , , yy, a, s] <- diag(n_regions)
  }
  stray <- array(0, dim = c(n_pop, n_yrs + 1))

  # Pad the movement arrays out to the NAA year dimension
  pad <- function(x) {
    out <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes))
    out[, , , seq_len(n_yrs), , , ] <- x
    out[, , , n_yrs + 1, , , ] <- x[, , , n_yrs, , , ]
    out
  }
  Mov <- pad(fx$Movement); Mra <- pad(fx$Mrate)

  for (tm in 0:2) {
    em <- SPoRC:::compute_biom_y(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, fx$seasdur,
                                 n_seas, n_pop, n_regions, n_ages, n_sexes,
                                 sgl, natal_region = 1, stray,
                                 Mov, Mra, tm, do_recruits_move = 0)

    # Minimal simulation environment: one replicate, everything with a trailing sim dim
    sim_env <- new.env()
    sim_env$NAA <- array(NAA, dim = c(dim(NAA), 1)); sim_env$NAA0 <- array(NAA0, dim = c(dim(NAA0), 1))
    sim_env$WAA <- array(WAA, dim = c(dim(WAA), 1)); sim_env$MatAA <- array(MatAA, dim = c(dim(MatAA), 1))
    sim_env$ZAA <- array(ZAA, dim = c(dim(ZAA), 1)); sim_env$natmort <- array(natmort, dim = c(dim(natmort), 1))
    sim_env$t_spawn <- t_spawn; sim_env$seasdur <- fx$seasdur
    sim_env$n_pop <- n_pop; sim_env$n_regions <- n_regions; sim_env$n_seas <- n_seas
    sim_env$n_ages <- n_ages; sim_env$n_sexes <- n_sexes
    sim_env$natal_region <- 1; sim_env$stray_rate <- array(stray, dim = c(dim(stray), 1))
    sim_env$sgl_seas_spawning_movement <- array(sgl, dim = c(dim(sgl), 1))
    sim_env$Movement <- array(Mov, dim = c(dim(Mov), 1)); sim_env$Mrate <- array(Mra, dim = c(dim(Mra), 1))
    sim_env$move_timing <- tm; sim_env$do_recruits_move <- 0

    om <- SPoRC:::compute_biom_y_sim(y, seas, sim = 1, sim_env = sim_env)

    expect_equal(om$SSB_y, em$SSB_y, tolerance = 1e-10,
                 label = sprintf("SSB parity, move_timing = %d", tm))
    expect_equal(om$Total_Biom_y, em$Total_Biom_y, tolerance = 1e-10,
                 label = sprintf("Total_Biom parity, move_timing = %d", tm))
    expect_equal(om$Dynamic_SSB0_y, em$Dynamic_SSB0_y, tolerance = 1e-10,
                 label = sprintf("Dynamic_SSB0 parity, move_timing = %d", tm))
  }
})

test_that("spawning biomass actually depends on move_timing when mortality varies by region", {
  # Guards the parity test above against being vacuously true: if all three timings gave
  # the same answer, a simulator that ignored move_timing would still pass.
  set.seed(12)
  fx <- make_setup()
  n <- fx$n_regions; n_ages <- fx$n_ages
  Z <- matrix(0, n, n_ages)
  for (r in seq_len(n)) Z[r, ] <- fx$M + fx$Fr[r]
  N <- matrix(stats::runif(n * n_ages, 20, 200), n, n_ages)

  ssb <- sapply(0:2, function(tm) {
    tot <- rep(0, n)
    for (a in seq_len(n_ages)) {
      tot <- tot + spawn_state(N[, a], fx$Movement[1, , , 1, 1, a, 1], Z[, a],
                               fx$Mrate[1, , , 1, 1, a, 1], fx$seasdur[1], 0.3, tm)
    }
    tot
  })

  expect_false(isTRUE(all.equal(ssb[, 1], ssb[, 3], tolerance = 1e-4)))
  expect_false(isTRUE(all.equal(ssb[, 2], ssb[, 3], tolerance = 1e-4)))
})

# ---------------------------------------------------------------------------
# 2. Initial age structure
# ---------------------------------------------------------------------------

test_that("Get_Init_NAA respects move_timing", {
  # The simulator builds its initial age structure through the same Get_Init_NAA the
  # estimation model uses, so it has to pass Mrate and move_timing through. If either
  # were dropped the equilibrium would silently be built at timing 0.
  fx <- make_setup()
  n_regions <- fx$n_regions; n_ages <- fx$n_ages; n_seas <- fx$n_seas

  init_F <- array(0, dim = c(n_regions, n_seas, 1)); init_F[, 1, 1] <- fx$Fr

  call_init <- function(tm) Get_Init_NAA(
    init_age_strc = 2,
    init_iter = n_ages * 5,
    n_pop = 1,
    n_regions = n_regions,
    n_sexes = 1,
    n_ages = n_ages,
    n_seas = n_seas,
    n_fish_fleets = 1,
    seasdur = fx$seasdur,
    natmort = array(fx$M, dim = c(1, n_regions, n_seas, n_ages, 1)),
    init_F = init_F,
    dmr = array(0, dim = c(n_regions, n_seas, 1)),
    fish_sel = array(1, dim = c(1, n_regions, n_seas, n_ages, 1, 1)),
    ret_sel = array(1, dim = c(1, n_regions, n_seas, n_ages, 1, 1)),
    R0_r = array(15, dim = c(1, n_regions)),
    rec_seas_prop = array(1 / n_seas, dim = c(1, n_seas)),
    sexratio = array(1, dim = c(1, n_regions, 1)),
    Movement = array(fx$Movement[, , , 1, , , ], dim = c(1, n_regions, n_regions, n_seas, n_ages, 1)),
    do_recruits_move = 0,
    ln_InitDevs = array(0, dim = c(1, n_regions, n_ages - 1)),
    Mrate = array(fx$Mrate[, , , 1, , , ], dim = c(1, n_regions, n_regions, n_seas, n_ages, 1)),
    move_timing = tm
  )

  naa <- lapply(0:2, call_init)
  for (tm in 1:3) {
    expect_true(all(is.finite(naa[[tm]])))
    expect_true(all(naa[[tm]] >= 0))
  }
  # The three timings must differ, otherwise dropping move_timing would be harmless
  expect_false(isTRUE(all.equal(naa[[1]], naa[[3]], tolerance = 1e-4)))
  expect_false(isTRUE(all.equal(naa[[2]], naa[[3]], tolerance = 1e-4)))
})

# ---------------------------------------------------------------------------
# 3. Simulation environment plumbing
# ---------------------------------------------------------------------------

test_that("Setup_sim_env always binds Mrate and rejects timing 2 without it", {
  # A NULL list element is dropped by list2env, so Mrate would escape to the enclosing
  # frame and error at move_timing = 1. It must resolve to NULL instead, and continuous
  # movement must fail loudly rather than silently reverting.
  base <- list(n_sims = 1, n_yrs = 2, n_regions = 2)

  env0 <- Setup_sim_env(base)
  expect_equal(env0$move_timing, 0)
  expect_true(exists("Mrate", envir = env0, inherits = FALSE))
  expect_null(env0$Mrate)

  env1 <- Setup_sim_env(c(base, list(move_timing = 1)))
  expect_true(exists("Mrate", envir = env1, inherits = FALSE))
  expect_null(env1$Mrate)

  expect_error(Setup_sim_env(c(base, list(move_timing = 2))), "requires Mrate")
})
