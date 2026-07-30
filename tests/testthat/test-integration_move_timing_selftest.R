library(SPoRC)
library(testthat)
library(Matrix)

# Self-consistency checks for continuous movement (move_timing = 2).
#
# The end-to-end shape guards live in test-move-timing-end-to-end.R. What this file
# checks is that the continuous-movement machinery is internally coherent:
#
#   1. reference points computed under move_timing = 2 are an actual fixed point of the
#      timing-2 projection (fishing at F_ref holds SSB at B_ref indefinitely);
#   2. the CTMC diffusion parameter is recoverable from data generated under continuous
#      movement, i.e. the movement fractions the model produces invert back to the
#      generator that produced them;
#   3. SSB is reproduced when the same population is stepped forward two ways.
#
# These are the properties that would break first if any module were still applying
# movement and mortality in the wrong order.

# ---------------------------------------------------------------------------
# 1. Reference points are a fixed point of the projection under every timing
# ---------------------------------------------------------------------------

test_that("equilibrium SSB is a fixed point of the seasonal operator under every timing", {
  # A per-recruit style equilibrium: constant recruitment into region 1, movement and
  # mortality applied via the shared operator, iterated to convergence. The converged
  # state must satisfy N = T(N) + R exactly, for all three timings. Any ordering
  # mismatch between the operator used to iterate and the one used to check would
  # show up immediately here.
  set.seed(2024)
  n <- 3; n_ages <- 12; n_seas <- 2
  seasdur <- c(0.4, 0.6)

  Mv <- array(0, c(n, n, n_seas)); Qv <- array(0, c(n, n, n_seas))
  for (s in seq_len(n_seas)) {
    D <- matrix(stats::runif(n * n, 0.05, 0.4), n, n)
    diag(D) <- 0; diag(D) <- -colSums(D)
    Qv[, , s] <- t(D)
    Mv[, , s] <- t(as.matrix(Matrix::expm(D * seasdur[s])))
  }
  # region-varying mortality, so the three timings genuinely differ
  Z <- outer(c(0.25, 0.45, 0.15), seasdur)   # [region, season]
  R <- c(1000, 0, 0)

  equilibrate <- function(tm, iters = 400) {
    N <- matrix(0, n, n_ages)
    for (i in seq_len(iters)) {
      N[, 1] <- R
      for (s in seq_len(n_seas)) {
        stepped <- N
        for (a in seq_len(n_ages)) {
          stepped[, a] <- advance_seas(N[, a], Mv[, , s], Z[, s], Qv[, , s], seasdur[s], tm)
        }
        if (s < n_seas) {
          N <- stepped
        } else {
          plus <- stepped[, n_ages]
          N[, 2:n_ages] <- stepped[, 1:(n_ages - 1)]
          N[, n_ages] <- N[, n_ages] + plus
        }
      }
    }
    N
  }

  for (tm in 0:2) {
    Neq <- equilibrate(tm)

    # Step the converged state forward one more full year; it must not move
    Nnext <- Neq
    Nnext[, 1] <- R
    for (s in seq_len(n_seas)) {
      stepped <- Nnext
      for (a in seq_len(n_ages)) {
        stepped[, a] <- advance_seas(Nnext[, a], Mv[, , s], Z[, s], Qv[, , s], seasdur[s], tm)
      }
      if (s < n_seas) {
        Nnext <- stepped
      } else {
        plus <- stepped[, n_ages]
        Nnext[, 2:n_ages] <- stepped[, 1:(n_ages - 1)]
        Nnext[, n_ages] <- Nnext[, n_ages] + plus
      }
    }

    expect_equal(Nnext, Neq, tolerance = 1e-8,
                 label = sprintf("equilibrium is a fixed point, move_timing = %d", tm))
    expect_true(all(Neq >= 0))
  }
})

test_that("the analytic plus group agrees with brute-force iteration under every timing", {
  # solve_plus_group's geometric series is what reference points rely on. It must equal
  # simply iterating the plus-group recursion many times, for each timing -- this is the
  # link between the reference point machinery and the projection dynamics.
  set.seed(99)
  n <- 3; n_seas <- 2
  seasdur <- c(0.35, 0.65)
  Mv <- array(0, c(n, n, n_seas)); Qv <- array(0, c(n, n, n_seas))
  for (s in seq_len(n_seas)) {
    D <- matrix(stats::runif(n * n, 0.05, 0.4), n, n)
    diag(D) <- 0; diag(D) <- -colSums(D)
    Qv[, , s] <- t(D); Mv[, , s] <- t(as.matrix(Matrix::expm(D * seasdur[s])))
  }
  M_pen <- c(0.12, 0.2, 0.09); M_plus <- c(0.13, 0.22, 0.1)
  F_pen <- matrix(stats::runif(n * n_seas, 0.02, 0.25), n, n_seas)
  F_plus <- matrix(stats::runif(n * n_seas, 0.02, 0.25), n, n_seas)
  N_penult <- c(50, 20, 35)

  for (tm in 0:2) {
    Ts <- build_plus_group_T(M_penult = M_pen, M_plus = M_plus,
                             F_penult = F_pen, F_plus = F_plus,
                             Mov_penult = Mv, Mov_plus = Mv,
                             n_regions = n, n_seas = n_seas, seasdur = seasdur,
                             Mrate_penult = Qv, Mrate_plus = Qv, move_timing = tm)

    analytic <- as.vector(solve_plus_group(Ts, N_penult, N_penult, n)$fished)

    # Brute force: repeatedly feed penultimate survivors in and let the plus group accrue
    brute <- rep(0, n)
    inflow <- as.vector(Ts$T_penult_fished %*% N_penult)
    for (i in 1:2000) brute <- inflow + as.vector(Ts$T_plus_fished %*% brute)

    expect_equal(analytic, brute, tolerance = 1e-8,
                 label = sprintf("analytic plus group vs iteration, move_timing = %d", tm))
  }
})

test_that("projecting at each timing's own F40% equilibrates at its own B40%", {
  # The end-to-end version of the fixed-point property, through the real pipeline:
  # global_SPR (reference points) -> Do_Population_Projection (forward dynamics).
  # Fishing at F40% forever must drive SSB to exactly the B40% the reference point
  # machinery reports, for every timing. This is the strongest available check that
  # the per-recruit equilibrium, the plus-group solve, spawn_state and the projection
  # dynamics all implement the same seasonal operator -- a mismatch in any one of them
  # shows up here as a projection that settles somewhere other than B40%.
  #
  # Natural mortality varies by region so the three timings genuinely differ; with a
  # region-invariant Z they would coincide and the test would be vacuous.
  n_pop <- 1; n_regions <- 3; n_seas <- 1; n_sexes <- 1; n_ages <- 20; n_flt <- 1
  seasdur <- rep(1 / n_seas, n_seas); ages <- seq_len(n_ages)
  t_spawn <- 0.25; spawn_seas <- 1
  NY <- 120                      # projection years; converged to ~1e-10 well before this
  R_tot <- 100
  rec_prop <- c(0.5, 0.3, 0.2)
  M_r <- c(0.12, 0.20, 0.35)

  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L
  ctmc_dat <- expand.grid(pop = 1, regions = seq_len(n_regions), years = 1,
                          seas = seq_len(n_seas), ages = ages, sexes = 1)
  mv <- Get_Movement(
    move_type = 1, do_recruits_move = 0, n_pop = 1, n_regions = n_regions, n_yrs = 1,
    n_proj_yrs_devs = 0, n_ages = n_ages, n_sexes = 1, n_seas = n_seas,
    move_pars = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, n_ages, 1)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, n_ages, 1)),
    use_fixed_movement = 0, Fixed_Movement = NULL, ctmc_move_dat = ctmc_dat,
    preference_formula = ~ 1, diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.35), move_preference_pars = 0,
    area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
    seasdur = seasdur, ctmc_scale_by_seasdur = 1)
  Mov1 <- mv$Movement[1, , , 1, 1, , 1]; Mra1 <- mv$Mrate[1, , , 1, 1, , 1]

  sel <- 1 / (1 + exp(-1.2 * (ages - 5)))
  waa <- 6 / (1 + exp(-0.5 * (ages - 6)))
  mat <- 1 / (1 + exp(-0.9 * (ages - 6)))

  spr_data <- function(tm) list(
    t_spawn = t_spawn, n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
    n_pop = n_pop, n_ages = n_ages, n_regions = n_regions,
    F_fract_flt = array(1, dim = c(n_regions, n_seas, n_flt)),
    dmr = array(0, dim = c(n_regions, n_seas, n_flt)),
    fish_sel = array(rep(sel, each = n_pop * n_regions * n_seas), dim = c(n_pop, n_regions, n_seas, n_ages, n_flt)),
    ret_sel = array(1, dim = c(n_pop, n_regions, n_seas, n_ages, n_flt)),
    natmort = array(rep(M_r, times = n_ages), dim = c(n_pop, n_regions, n_ages)),
    WAA = array(rep(waa, each = n_pop * n_regions * n_seas), dim = c(n_pop, n_regions, n_seas, n_ages)),
    MatAA = array(rep(mat, each = n_pop * n_regions * n_seas), dim = c(n_pop, n_regions, n_seas, n_ages)),
    Movement = array(Mov1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
    Mrate = array(Mra1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
    move_timing = tm, do_recruits_move = 0,
    sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)),
    rec_seas_prop = array(1 / n_seas, dim = c(n_pop, n_seas)),
    stray_rate = array(0, dim = n_pop),
    sgl_seas_spawning_movement = array(rep(diag(n_regions), n_ages), dim = c(n_pop, n_regions, n_regions, n_ages)),
    natal_region = 1, n_pop_in_region = array(1, dim = n_regions),
    SPR_x = 0.40, rec_region_prop = array(rec_prop, dim = c(n_pop, n_regions)))

  project <- function(tm, Fval) {
    six <- function(v) array(rep(v, each = n_pop * n_regions * NY * n_seas),
                             dim = c(n_pop, n_regions, NY, n_seas, n_ages, n_sexes))
    svn <- function(v) array(rep(v, each = n_pop * n_regions * NY * n_seas),
                             dim = c(n_pop, n_regions, NY, n_seas, n_ages, n_sexes, n_flt))
    Mov <- array(0, c(n_pop, n_regions, n_regions, NY, n_seas, n_ages, n_sexes))
    Mra <- array(0, c(n_pop, n_regions, n_regions, NY, n_seas, n_ages, n_sexes))
    sgl <- array(0, c(n_pop, n_regions, n_regions, NY, n_ages, n_sexes))
    for (y in seq_len(NY)) for (a in seq_len(n_ages)) {
      Mov[1, , , y, 1, a, 1] <- Mov1[, , a]; Mra[1, , , y, 1, a, 1] <- Mra1[, , a]
      sgl[1, , , y, a, 1] <- diag(n_regions)
    }
    tNAA <- array(0, c(n_pop, n_regions, n_seas, n_ages, n_sexes)); tNAA[1, , 1, , 1] <- 10

    Do_Population_Projection(
      n_proj_yrs = NY, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages, n_sexes = n_sexes,
      sexratio = array(1, dim = c(n_pop, n_regions, NY, n_sexes)),
      n_fish_fleets = n_flt, do_recruits_move = 0,
      recruitment = array(rep(R_tot * rec_prop, each = n_pop), dim = c(n_pop, n_regions, NY)),
      terminal_NAA = tNAA, terminal_NAA0 = tNAA,
      terminal_F = array(Fval, dim = c(n_regions, n_seas, n_flt)),
      dmr = array(0, dim = c(n_regions, n_seas, n_flt)),
      natmort = array(rep(M_r, times = NY * n_ages), dim = c(n_pop, n_regions, NY, n_ages, n_sexes)),
      natal_region = 1, WAA = six(waa), MatAA = six(mat), WAA_fish = svn(waa),
      fish_sel = svn(sel), ret_sel = array(1, dim = c(n_pop, n_regions, NY, n_seas, n_ages, n_sexes, n_flt)),
      Movement = Mov, sgl_seas_spawning_movement = sgl, stray_rate = array(0, dim = c(n_pop, NY)),
      f_ref_pt = array(Fval, dim = c(n_regions, NY)), b_ref_pt = NULL, HCR_function = NULL,
      recruitment_opt = "mean_rec", fmort_opt = "Input", t_spawn = t_spawn,
      n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
      rec_seas_prop = array(1 / n_seas, dim = c(n_pop, n_seas)),
      Mrate = Mra, move_timing = tm)
  }

  F40 <- numeric(3); B40 <- numeric(3)
  for (tm in 0:2) {
    rp <- SPoRC:::optim_ref_pts(SPoRC:::global_SPR, spr_data(tm), list(log_F_x = log(0.1)))
    F40[tm + 1] <- rp$rep$F_x
    B40[tm + 1] <- sum(rp$rep$SB) * R_tot

    pj40 <- project(tm, rp$rep$F_x)
    pj0  <- project(tm, 1e-12)

    expect_equal(rp$rep$SPR, 0.40, tolerance = 1e-6,
                 label = sprintf("SPR solve converged, move_timing = %d", tm))
    # fishing at F40% forever holds SSB exactly at B40%
    expect_equal(sum(pj40$proj_SSB[, , NY]), sum(rp$rep$SB) * R_tot, tolerance = 1e-8,
                 label = sprintf("projected SSB == B40%%, move_timing = %d", tm))
    # and the unfished control lands on SSB0
    expect_equal(sum(pj0$proj_SSB[, , NY]), sum(rp$rep$SB0) * R_tot, tolerance = 1e-8,
                 label = sprintf("projected unfished SSB == SSB0, move_timing = %d", tm))
    # depletion is 40% of the same timing's unfished equilibrium
    expect_equal(sum(pj40$proj_SSB[, , NY]) / sum(pj0$proj_SSB[, , NY]), 0.40, tolerance = 1e-8)
    # the spatial split must match too, not just the total
    expect_equal(as.vector(pj40$proj_SSB[1, , NY]), as.vector(rp$rep$SB[1, ]) * R_tot,
                 tolerance = 1e-8,
                 label = sprintf("region-wise SSB == B40%% by region, move_timing = %d", tm))
    # genuinely at equilibrium rather than still drifting
    expect_equal(sum(pj40$proj_SSB[, , NY]), sum(pj40$proj_SSB[, , NY - 1]), tolerance = 1e-8)
  }

  # Guard against vacuity: the timings must actually disagree on the biomass reference
  # points here, otherwise the fixed-point check above could pass on a single shared answer.
  expect_gt(diff(range(B40)) / mean(B40), 1e-3)
})

test_that("MSY and SPR reference points agree on per-recruit biology at every timing", {
  # global_SPR and the Beverton-Holt Fmsy routines build the same unfished per-recruit
  # age structure by different code paths, so their unfished spawning biomass per recruit
  # must agree exactly -- at every move_timing, not just the default.
  #
  # This is the check that exposed the MSY routines honouring move_timing only inside
  # build_plus_group_T while their age loops, spawning propagation and catch equation
  # stayed hard-coded to movement-then-mortality. Before the fix these disagreed by 3.1%
  # at timing 1 and 1.3% at timing 2, while agreeing exactly at timing 0.
  n_regions <- 3; n_ages <- 20; n_seas <- 1; n_flt <- 1
  seasdur <- rep(1 / n_seas, n_seas); ages <- seq_len(n_ages)
  t_spawn <- 0.25; spawn_seas <- 1
  M_r <- c(0.12, 0.20, 0.35)          # region-varying, so the timings genuinely differ
  rec_prop <- c(0.5, 0.3, 0.2)
  Fv <- 0.2

  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L
  ctmc <- expand.grid(pop = 1, regions = seq_len(n_regions), years = 1,
                      seas = seq_len(n_seas), ages = ages, sexes = 1)
  mv <- Get_Movement(
    move_type = 1, do_recruits_move = 0, n_pop = 1, n_regions = n_regions, n_yrs = 1,
    n_proj_yrs_devs = 0, n_ages = n_ages, n_sexes = 1, n_seas = n_seas,
    move_pars = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, n_ages, 1)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, n_ages, 1)),
    use_fixed_movement = 0, Fixed_Movement = NULL, ctmc_move_dat = ctmc,
    preference_formula = ~ 1, diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.35), move_preference_pars = 0,
    area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
    seasdur = seasdur, ctmc_scale_by_seasdur = 1)
  Mov1 <- mv$Movement[1, , , 1, 1, , 1]; Mra1 <- mv$Mrate[1, , , 1, 1, , 1]

  sel <- 1 / (1 + exp(-1.2 * (ages - 5)))
  waa <- 6 / (1 + exp(-0.5 * (ages - 6)))
  mat <- 1 / (1 + exp(-0.9 * (ages - 6)))

  msy_data <- function(tm) list(
    t_spawn = t_spawn, n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
    n_ages = n_ages, n_regions = n_regions,
    F_fract_flt = array(1, dim = c(n_regions, n_seas, n_flt)),
    dmr = array(0, dim = c(n_regions, n_seas, n_flt)),
    fish_sel = array(rep(sel, each = n_regions * n_seas), dim = c(1, n_regions, n_seas, n_ages, n_flt)),
    ret_sel = array(1, dim = c(1, n_regions, n_seas, n_ages, n_flt)),
    Movement = array(Mov1, dim = c(n_regions, n_regions, n_seas, n_ages)),
    Mrate = array(Mra1, dim = c(n_regions, n_regions, n_seas, n_ages)),
    move_timing = tm, do_recruits_move = 0,
    natmort = array(rep(M_r, times = n_ages), dim = c(n_regions, n_ages)),
    WAA = array(rep(waa, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    WAA_fish = array(rep(waa, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    MatAA = array(rep(mat, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    sex_ratio_f = array(0.5, dim = n_regions), rec_seas_prop = array(1 / n_seas, dim = n_seas),
    rec_region_prop = array(rec_prop, dim = n_regions),
    is_discard_fleet = array(0, dim = n_flt),
    SR_type = 1, h = 0.7, R0 = 100, ln_global_R0 = log(100))

  spr_data <- function(tm) list(
    t_spawn = t_spawn, n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
    n_pop = 1, n_ages = n_ages, n_regions = n_regions,
    F_fract_flt = array(1, dim = c(n_regions, n_seas, n_flt)),
    dmr = array(0, dim = c(n_regions, n_seas, n_flt)),
    fish_sel = array(rep(sel, each = n_regions * n_seas), dim = c(1, n_regions, n_seas, n_ages, n_flt)),
    ret_sel = array(1, dim = c(1, n_regions, n_seas, n_ages, n_flt)),
    natmort = array(rep(M_r, times = n_ages), dim = c(1, n_regions, n_ages)),
    WAA = array(rep(waa, each = n_regions * n_seas), dim = c(1, n_regions, n_seas, n_ages)),
    MatAA = array(rep(mat, each = n_regions * n_seas), dim = c(1, n_regions, n_seas, n_ages)),
    Movement = array(Mov1, dim = c(1, n_regions, n_regions, n_seas, n_ages)),
    Mrate = array(Mra1, dim = c(1, n_regions, n_regions, n_seas, n_ages)),
    move_timing = tm, do_recruits_move = 0,
    sex_ratio_f = array(0.5, dim = c(1, n_regions)),
    rec_seas_prop = array(1 / n_seas, dim = c(1, n_seas)),
    stray_rate = array(0, dim = 1),
    sgl_seas_spawning_movement = array(rep(diag(n_regions), n_ages), dim = c(1, n_regions, n_regions, n_ages)),
    natal_region = 1, n_pop_in_region = array(1, dim = n_regions),
    SPR_x = 0.40, rec_region_prop = array(rec_prop, dim = c(1, n_regions)))

  ev <- function(fn, dat, parname, val) {
    p <- list(); p[[parname]] <- log(val)
    RTMB::MakeADFun(SPoRC:::cmb(fn, dat), parameters = p, random = NULL, silent = TRUE)$report(log(val))
  }

  spr_sb0 <- numeric(3)
  for (tm in 0:2) {
    r_spr <- ev(SPoRC:::global_SPR,     spr_data(tm), "log_F_x",  Fv)
    r_msy <- ev(SPoRC:::global_BH_Fmsy, msy_data(tm), "log_Fmsy", Fv)
    spr_sb0[tm + 1] <- sum(r_spr$SB0)

    expect_equal(r_msy$SBPR_0, sum(r_spr$SB0), tolerance = 1e-10,
                 label = sprintf("unfished SBPR, global_BH_Fmsy vs global_SPR, move_timing = %d", tm))
    expect_equal(r_msy$SBPR_F, sum(r_spr$SB), tolerance = 1e-10,
                 label = sprintf("fished SBPR, global_BH_Fmsy vs global_SPR, move_timing = %d", tm))
    expect_equal(r_msy$SPR, r_spr$SPR, tolerance = 1e-10,
                 label = sprintf("SPR agrees, move_timing = %d", tm))
  }

  # Vacuity guard: the per-recruit biology must genuinely depend on move_timing here,
  # otherwise a routine that ignored the flag entirely would still pass the checks above.
  expect_gt(diff(range(spr_sb0)) / mean(spr_sb0), 1e-3)
})

test_that("projecting at Fmsy under Beverton-Holt feedback equilibrates at Bmsy and Req", {
  # The MSY analogue of the F40% fixed-point test above, and the sharper version of it:
  # here recruitment is not pinned, it is regenerated each year from the Beverton-Holt
  # curve, so the per-recruit biology, the catch equation, the plus group and the
  # stock-recruit machinery all have to agree for the loop to settle on Bmsy.
  #
  # Two things this pins that were previously broken:
  #   1. Get_Det_Recruitment's global density-dependence branch skipped the spawning-season
  #      movement on the plus group ("mortality, no movement" = move_timing 1 semantics
  #      hard-coded), so timings 0 and 2 settled slightly off Bmsy while 1 was exact.
  #   2. Do_Population_Projection's rec_lag != 0 Beverton-Holt call did not forward
  #      Mrate/move_timing, so the SSB0 behind the curve was built at timing 0.
  #
  # NOTE: global_BH_Fmsy assumes GLOBAL density dependence, so bh_rec_opt$rec_dd must be 1.
  # Pairing it with rec_dd = 0 (local) compares two different equilibria and looks like a
  # bug when it is not.
  n_pop <- 1; n_regions <- 3; n_seas <- 1; n_sexes <- 1; n_ages <- 20; n_flt <- 1
  seasdur <- rep(1 / n_seas, n_seas); ages <- seq_len(n_ages)
  t_spawn <- 0.25; spawn_seas <- 1
  # The Beverton-Holt feedback loop damps geometrically and converges more slowly than
  # the age structure alone: relative error in equilibrium SSB is 1.3e-9 at 150 years,
  # 1.5e-12 at 200 and ~2e-15 by 250. 280 leaves headroom under the 1e-10 tolerance.
  NY <- 280
  M_r <- c(0.12, 0.20, 0.35); rec_prop <- c(0.5, 0.3, 0.2)
  R0 <- 100; h <- 0.7

  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L
  ctmc <- expand.grid(pop = 1, regions = seq_len(n_regions), years = 1,
                      seas = seq_len(n_seas), ages = ages, sexes = 1)
  mv <- Get_Movement(
    move_type = 1, do_recruits_move = 0, n_pop = 1, n_regions = n_regions, n_yrs = 1,
    n_proj_yrs_devs = 0, n_ages = n_ages, n_sexes = 1, n_seas = n_seas,
    move_pars = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, n_ages, 1)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, n_ages, 1)),
    use_fixed_movement = 0, Fixed_Movement = NULL, ctmc_move_dat = ctmc,
    preference_formula = ~ 1, diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.35), move_preference_pars = 0,
    area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
    seasdur = seasdur, ctmc_scale_by_seasdur = 1)
  Mov1 <- mv$Movement[1, , , 1, 1, , 1]; Mra1 <- mv$Mrate[1, , , 1, 1, , 1]

  sel <- 1 / (1 + exp(-1.2 * (ages - 5)))
  waa <- 6 / (1 + exp(-0.5 * (ages - 6)))
  mat <- 1 / (1 + exp(-0.9 * (ages - 6)))

  msy_data <- function(tm) list(
    t_spawn = t_spawn, n_seas = n_seas, seasdur = seasdur, spawn_seas = spawn_seas,
    n_ages = n_ages, n_regions = n_regions,
    F_fract_flt = array(1, dim = c(n_regions, n_seas, n_flt)),
    dmr = array(0, dim = c(n_regions, n_seas, n_flt)),
    fish_sel = array(rep(sel, each = n_regions * n_seas), dim = c(1, n_regions, n_seas, n_ages, n_flt)),
    ret_sel = array(1, dim = c(1, n_regions, n_seas, n_ages, n_flt)),
    Movement = array(Mov1, dim = c(n_regions, n_regions, n_seas, n_ages)),
    Mrate = array(Mra1, dim = c(n_regions, n_regions, n_seas, n_ages)),
    move_timing = tm, do_recruits_move = 0,
    natmort = array(rep(M_r, times = n_ages), dim = c(n_regions, n_ages)),
    WAA = array(rep(waa, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    WAA_fish = array(rep(waa, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    MatAA = array(rep(mat, each = n_regions * n_seas), dim = c(n_regions, n_seas, n_ages)),
    sex_ratio_f = array(0.5, dim = n_regions), rec_seas_prop = array(1 / n_seas, dim = n_seas),
    rec_region_prop = array(rec_prop, dim = n_regions),
    is_discard_fleet = array(0, dim = n_flt), SR_type = 1, h = h, R0 = R0, ln_global_R0 = log(R0))

  six <- function(v) array(rep(v, each = n_pop * n_regions * NY * n_seas),
                           dim = c(n_pop, n_regions, NY, n_seas, n_ages, n_sexes))
  svn <- function(v) array(rep(v, each = n_pop * n_regions * NY * n_seas),
                           dim = c(n_pop, n_regions, NY, n_seas, n_ages, n_sexes, n_flt))

  project <- function(tm, Fval, rec_mode, Rconst = NULL) {
    Mov <- array(0, c(n_pop, n_regions, n_regions, NY, n_seas, n_ages, n_sexes))
    Mra <- array(0, c(n_pop, n_regions, n_regions, NY, n_seas, n_ages, n_sexes))
    sgl <- array(0, c(n_pop, n_regions, n_regions, NY, n_ages, n_sexes))
    for (y in seq_len(NY)) for (a in seq_len(n_ages)) {
      Mov[1, , , y, 1, a, 1] <- Mov1[, , a]; Mra[1, , , y, 1, a, 1] <- Mra1[, , a]
      sgl[1, , , y, a, 1] <- diag(n_regions)
    }
    tNAA <- array(0, c(n_pop, n_regions, n_seas, n_ages, n_sexes)); tNAA[1, , 1, , 1] <- 10

    bh <- list(R0 = R0, h = array(h, dim = c(n_pop, n_regions)),
               rec_region_prop = array(rec_prop, dim = c(n_pop, n_regions)),
               rec_seas_prop = array(1 / n_seas, dim = c(n_pop, n_seas)),
               SSB = array(1, dim = c(n_pop, n_regions, 1)),
               WAA = array(rep(waa, each = n_pop * n_regions * n_seas), dim = c(n_pop, n_regions, n_seas, n_ages)),
               MatAA = array(rep(mat, each = n_pop * n_regions * n_seas), dim = c(n_pop, n_regions, n_seas, n_ages)),
               natmort = array(rep(M_r, times = n_ages), dim = c(n_pop, n_regions, n_ages)),
               Movement = array(Mov1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
               Mrate = array(Mra1, dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
               sgl_seas_spawning_movement = array(rep(diag(n_regions), n_ages),
                                                  dim = c(n_pop, n_regions, n_regions, n_ages)),
               stray_rate = array(0, dim = n_pop),
               init_F = array(0, dim = c(n_regions, n_seas, n_flt)),
               fish_sel = array(rep(sel, each = n_pop * n_regions * n_seas),
                                dim = c(n_pop, n_regions, n_seas, n_ages, n_flt)),
               ret_sel = array(1, dim = c(n_pop, n_regions, n_seas, n_ages, n_flt)),
               dmr = array(0, dim = c(n_regions, n_seas, n_flt)),
               sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)),
               rec_dd = 1, rec_lag = 1)   # rec_dd = 1: global DD, matching global_BH_Fmsy

    Do_Population_Projection(
      n_proj_yrs = NY, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages, n_sexes = n_sexes,
      sexratio = array(1, dim = c(n_pop, n_regions, NY, n_sexes)), n_fish_fleets = n_flt,
      do_recruits_move = 0,
      recruitment = array(rep(if (is.null(Rconst)) R0 * rec_prop else Rconst * rec_prop, each = n_pop),
                          dim = c(n_pop, n_regions, NY)),
      terminal_NAA = tNAA, terminal_NAA0 = tNAA,
      terminal_F = array(Fval, dim = c(n_regions, n_seas, n_flt)),
      dmr = array(0, dim = c(n_regions, n_seas, n_flt)),
      natmort = array(rep(M_r, times = NY * n_ages), dim = c(n_pop, n_regions, NY, n_ages, n_sexes)),
      natal_region = 1, WAA = six(waa), MatAA = six(mat), WAA_fish = svn(waa),
      fish_sel = svn(sel), ret_sel = array(1, dim = c(n_pop, n_regions, NY, n_seas, n_ages, n_sexes, n_flt)),
      Movement = Mov, sgl_seas_spawning_movement = sgl, stray_rate = array(0, dim = c(n_pop, NY)),
      f_ref_pt = array(Fval, dim = c(n_regions, NY)), b_ref_pt = NULL, HCR_function = NULL,
      recruitment_opt = rec_mode, bh_rec_opt = if (rec_mode == "bh_rec") bh else NULL,
      fmort_opt = "Input", t_spawn = t_spawn, n_seas = n_seas, seasdur = seasdur,
      spawn_seas = spawn_seas, rec_seas_prop = array(1 / n_seas, dim = c(n_pop, n_seas)),
      Mrate = Mra, move_timing = tm)
  }

  Bmsy_all <- numeric(3)
  for (tm in 0:2) {
    rp <- SPoRC:::optim_ref_pts(SPoRC:::global_BH_Fmsy, msy_data(tm), list(log_Fmsy = log(0.1)))
    Fmsy <- rp$rep$Fmsy; Req <- rp$rep$Req
    Bmsy <- sum(apply(rp$rep$SB_age[2, , , drop = FALSE], 2, sum)) * Req
    Bmsy_all[tm + 1] <- Bmsy

    pjA <- project(tm, Fmsy, "mean_rec", Rconst = Req)   # recruitment pinned at Req
    pjB <- project(tm, Fmsy, "bh_rec")                   # full Beverton-Holt feedback

    # (a) with recruitment pinned, SSB must land exactly on Bmsy
    expect_equal(sum(pjA$proj_SSB[, , NY]), Bmsy, tolerance = 1e-10,
                 label = sprintf("pinned-recruitment SSB == Bmsy, move_timing = %d", tm))
    # (b) with the full stock-recruit loop closed, both SSB and recruitment must land there
    expect_equal(sum(pjB$proj_SSB[, , NY]), Bmsy, tolerance = 1e-10,
                 label = sprintf("Beverton-Holt SSB == Bmsy, move_timing = %d", tm))
    expect_equal(sum(pjB$proj_NAA[, , NY, 1, 1, ]), Req, tolerance = 1e-10,
                 label = sprintf("Beverton-Holt recruitment == Req, move_timing = %d", tm))
    # (c) projected yield equals MSY up to the single-sex convention: global_BH_Fmsy
    #     starts its cohort at rec_region_prop * sex_ratio_f (0.5), so its Yield is per
    #     female recruit while the projection puts the whole Req into the population.
    expect_equal(sum(pjA$proj_Catch[, , NY, , ]), 2 * rp$rep$Yield, tolerance = 1e-10,
                 label = sprintf("projected yield == 2 x MSY, move_timing = %d", tm))
    # genuinely at equilibrium rather than still drifting
    expect_equal(sum(pjB$proj_SSB[, , NY]), sum(pjB$proj_SSB[, , NY - 1]), tolerance = 1e-8)
  }

  # Vacuity guard: Bmsy must actually depend on move_timing here
  expect_gt(diff(range(Bmsy_all)) / mean(Bmsy_all), 1e-3)
})

test_that("local_BH_MSY is a fixed point of a two-season projection under every timing", {
  # Widens the Fmsy fixed-point check along the two axes the single-season global test
  # does not reach: LOCAL density dependence (region-specific Fmsy, per-origin equilibrium
  # recruitment solved by Newton-Raphson) and SEASONALITY (two unequal seasons, spawning
  # in the second, so the seasonal operators compose and the "advance into spawning
  # season" branches actually execute).
  #
  # Two harness details that matter and are easy to get wrong:
  #   - bh_rec_opt$rec_dd must be 0 (local) to match local_BH_Fmsy_sglpop.
  #   - terminal_NAA must be seeded across ALL seasons. Projection year 1 IS the terminal
  #     data year (proj_NAA[,,1,,,] <- terminal_NAA, and the real caller passes
  #     rep$NAA[,,n_yrs,,,]), so its later seasons are inputs, not something the
  #     projection recomputes. Seeding only season 1 leaves SSB[,,1] = 0, which with
  #     Beverton-Holt recruitment collapses the whole projection to zero.
  n_pop <- 1; n_regions <- 3; n_sexes <- 1; n_ages <- 20; n_flt <- 1
  NS <- 2; SPAWN <- 2
  seasdur <- c(0.4, 0.6); ages <- seq_len(n_ages); t_spawn <- 0.25
  NY <- 280
  M_r <- c(0.12, 0.20, 0.35); rec_prop <- c(0.5, 0.3, 0.2)
  R0 <- 100; h <- 0.7

  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L
  ctmc <- expand.grid(pop = 1, regions = seq_len(n_regions), years = 1,
                      seas = seq_len(NS), ages = ages, sexes = 1)
  mv <- Get_Movement(
    move_type = 1, do_recruits_move = 0, n_pop = 1, n_regions = n_regions, n_yrs = 1,
    n_proj_yrs_devs = 0, n_ages = n_ages, n_sexes = 1, n_seas = NS,
    move_pars = array(0, c(1, n_regions, n_regions - 1, 1, NS, n_ages, 1)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, 1, NS, n_ages, 1)),
    use_fixed_movement = 0, Fixed_Movement = NULL, ctmc_move_dat = ctmc,
    preference_formula = ~ 1, diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.35), move_preference_pars = 0,
    area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
    seasdur = seasdur, ctmc_scale_by_seasdur = 1)
  Mov1 <- array(mv$Movement[1, , , 1, , , 1], c(n_regions, n_regions, NS, n_ages))
  Mra1 <- array(mv$Mrate[1, , , 1, , , 1],    c(n_regions, n_regions, NS, n_ages))

  sel <- 1 / (1 + exp(-1.2 * (ages - 5)))
  waa <- 6 / (1 + exp(-0.5 * (ages - 6)))
  mat <- 1 / (1 + exp(-0.9 * (ages - 6)))
  r4 <- function(v) array(rep(v, each = n_regions * NS), dim = c(n_regions, NS, n_ages))
  r5 <- function(v) array(rep(v, each = n_regions * NS), dim = c(1, n_regions, NS, n_ages, n_flt))

  msy_l <- function(tm) list(
    t_spawn = t_spawn, n_seas = NS, seasdur = seasdur, spawn_seas = SPAWN,
    n_ages = n_ages, n_regions = n_regions, n_pop = 1,
    F_fract_flt = array(1 / NS, dim = c(n_regions, NS, n_flt)),
    dmr = array(0, dim = c(n_regions, NS, n_flt)),
    fish_sel = r5(sel), ret_sel = array(1, dim = c(1, n_regions, NS, n_ages, n_flt)),
    Movement = Mov1, Mrate = Mra1, move_timing = tm, do_recruits_move = 0,
    natmort = array(rep(M_r, times = n_ages), dim = c(n_regions, n_ages)),
    WAA = r4(waa), WAA_fish = r4(waa), MatAA = r4(mat),
    rec_seas_prop = array(1 / NS, dim = NS), is_discard_fleet = array(0, dim = n_flt),
    sex_ratio_f = array(0.5, dim = n_regions), rec_region_prop = array(rec_prop, dim = n_regions),
    h = array(h, dim = n_regions), R0 = R0, newton_steps = 200, natal_region = 1,
    n_pop_in_region = array(1, dim = n_regions), stray_rate = array(0, dim = 1),
    sgl_seas_spawning_movement = array(rep(diag(n_regions), n_ages),
                                       dim = c(1, n_regions, n_regions, n_ages)))

  six <- function(v) array(rep(v, each = n_pop * n_regions * NY * NS),
                           dim = c(n_pop, n_regions, NY, NS, n_ages, n_sexes))
  svn <- function(v) array(rep(v, each = n_pop * n_regions * NY * NS),
                           dim = c(n_pop, n_regions, NY, NS, n_ages, n_sexes, n_flt))

  project <- function(tm, Fvec) {
    Mov <- array(0, c(n_pop, n_regions, n_regions, NY, NS, n_ages, n_sexes))
    Mra <- array(0, c(n_pop, n_regions, n_regions, NY, NS, n_ages, n_sexes))
    sgl <- array(0, c(n_pop, n_regions, n_regions, NY, n_ages, n_sexes))
    for (y in seq_len(NY)) {
      for (s in seq_len(NS)) for (a in seq_len(n_ages)) {
        Mov[1, , , y, s, a, 1] <- Mov1[, , s, a]; Mra[1, , , y, s, a, 1] <- Mra1[, , s, a]
      }
      for (a in seq_len(n_ages)) sgl[1, , , y, a, 1] <- diag(n_regions)
    }
    tNAA <- array(10, c(n_pop, n_regions, NS, n_ages, n_sexes))  # ALL seasons; see note above

    bh <- list(R0 = R0, h = array(h, dim = c(n_pop, n_regions)),
               rec_region_prop = array(rec_prop, dim = c(n_pop, n_regions)),
               rec_seas_prop = array(1 / NS, dim = c(n_pop, NS)),
               SSB = array(1, dim = c(n_pop, n_regions, 1)),
               WAA = array(r4(waa), c(n_pop, n_regions, NS, n_ages)),
               MatAA = array(r4(mat), c(n_pop, n_regions, NS, n_ages)),
               natmort = array(rep(M_r, times = n_ages), dim = c(n_pop, n_regions, n_ages)),
               Movement = array(Mov1, c(n_pop, n_regions, n_regions, NS, n_ages)),
               Mrate = array(Mra1, c(n_pop, n_regions, n_regions, NS, n_ages)),
               sgl_seas_spawning_movement = array(rep(diag(n_regions), n_ages),
                                                  dim = c(n_pop, n_regions, n_regions, n_ages)),
               stray_rate = array(0, dim = n_pop), init_F = array(0, dim = c(n_regions, NS, n_flt)),
               fish_sel = array(r5(sel), c(n_pop, n_regions, NS, n_ages, n_flt)),
               ret_sel = array(1, dim = c(n_pop, n_regions, NS, n_ages, n_flt)),
               dmr = array(0, dim = c(n_regions, NS, n_flt)),
               sex_ratio_f = array(0.5, dim = c(n_pop, n_regions)),
               rec_dd = 0, rec_lag = 1)   # rec_dd = 0: local DD, matching local_BH_Fmsy_sglpop

    Do_Population_Projection(
      n_proj_yrs = NY, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages, n_sexes = n_sexes,
      sexratio = array(1, dim = c(n_pop, n_regions, NY, n_sexes)), n_fish_fleets = n_flt,
      do_recruits_move = 0,
      recruitment = array(rep(R0 * rec_prop, each = n_pop), dim = c(n_pop, n_regions, NY)),
      terminal_NAA = tNAA, terminal_NAA0 = tNAA,
      terminal_F = array(rep(Fvec / NS, NS * n_flt), dim = c(n_regions, NS, n_flt)),
      dmr = array(0, dim = c(n_regions, NS, n_flt)),
      natmort = array(rep(M_r, times = NY * n_ages), dim = c(n_pop, n_regions, NY, n_ages, n_sexes)),
      natal_region = 1, WAA = six(waa), MatAA = six(mat), WAA_fish = svn(waa),
      fish_sel = svn(sel), ret_sel = array(1, dim = c(n_pop, n_regions, NY, NS, n_ages, n_sexes, n_flt)),
      Movement = Mov, sgl_seas_spawning_movement = sgl, stray_rate = array(0, dim = c(n_pop, NY)),
      f_ref_pt = matrix(Fvec, n_regions, NY), b_ref_pt = NULL, HCR_function = NULL,
      recruitment_opt = "bh_rec", bh_rec_opt = bh, fmort_opt = "Input", t_spawn = t_spawn,
      n_seas = NS, seasdur = seasdur, spawn_seas = SPAWN,
      rec_seas_prop = array(1 / NS, dim = c(n_pop, NS)), Mrate = Mra, move_timing = tm)
  }

  Btot <- numeric(3)
  for (tm in 0:2) {
    rp <- SPoRC:::optim_ref_pts(SPoRC:::local_BH_Fmsy_sglpop, msy_l(tm),
                                list(log_Fmsy = rep(log(0.1), n_regions)))
    Btgt <- sapply(seq_len(n_regions), function(r) sum(rp$rep$SB_fished_mat[, r] * rp$rep$Req_o))
    Btot[tm + 1] <- sum(Btgt)

    pj <- project(tm, rp$rep$Fmsy)
    ssb <- as.vector(pj$proj_SSB[1, , NY])

    # region-by-region, not just the total: local DD makes the spatial split the point
    expect_equal(ssb, Btgt, tolerance = 1e-8,
                 label = sprintf("two-season local Bmsy by region, move_timing = %d", tm))
    expect_equal(sum(ssb), sum(Btgt), tolerance = 1e-8)
    # genuinely settled rather than still drifting
    expect_equal(sum(pj$proj_SSB[1, , NY]), sum(pj$proj_SSB[1, , NY - 1]), tolerance = 1e-8)
    expect_true(all(ssb > 0))   # guards the all-zero collapse mode described above
  }

  # Vacuity guard: the timings must genuinely disagree on Bmsy in this configuration
  expect_gt(diff(range(Btot)) / mean(Btot), 1e-3)
})

test_that("seasonal recruitment is apportioned across regions consistently in MSY and SPR", {
  # global_BH_Fmsy seeds age 1 with rec_region_prop * sex_ratio_f * rec_seas_prop[1] but
  # previously topped up seasons 2..n with rec_seas_prop[seas] * sex_ratio_f only, dropping
  # the regional apportionment. With n_seas > 1 and non-uniform rec_region_prop that put
  # equilibrium SSB ~50% out. global_SPR builds the same quantity correctly, so comparing
  # their unfished per-recruit spawning biomass isolates it.
  #
  # This needs n_seas > 1 AND non-uniform rec_region_prop to bite -- the single-season
  # MSY/SPR agreement test above never enters the seasonal top-up branch.
  n_regions <- 3; n_ages <- 20; NS <- 2; n_flt <- 1
  seasdur <- c(0.4, 0.6); ages <- seq_len(n_ages); t_spawn <- 0.25; SPAWN <- 1
  M_r <- c(0.12, 0.20, 0.35)
  rec_prop <- c(0.5, 0.3, 0.2)            # deliberately non-uniform
  Fv <- 0.2

  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L
  ctmc <- expand.grid(pop = 1, regions = seq_len(n_regions), years = 1,
                      seas = seq_len(NS), ages = ages, sexes = 1)
  mv <- Get_Movement(
    move_type = 1, do_recruits_move = 0, n_pop = 1, n_regions = n_regions, n_yrs = 1,
    n_proj_yrs_devs = 0, n_ages = n_ages, n_sexes = 1, n_seas = NS,
    move_pars = array(0, c(1, n_regions, n_regions - 1, 1, NS, n_ages, 1)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, 1, NS, n_ages, 1)),
    use_fixed_movement = 0, Fixed_Movement = NULL, ctmc_move_dat = ctmc,
    preference_formula = ~ 1, diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.35), move_preference_pars = 0,
    area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
    seasdur = seasdur, ctmc_scale_by_seasdur = 1)
  Mov1 <- array(mv$Movement[1, , , 1, , , 1], c(n_regions, n_regions, NS, n_ages))
  Mra1 <- array(mv$Mrate[1, , , 1, , , 1],    c(n_regions, n_regions, NS, n_ages))

  sel <- 1 / (1 + exp(-1.2 * (ages - 5)))
  waa <- 6 / (1 + exp(-0.5 * (ages - 6)))
  mat <- 1 / (1 + exp(-0.9 * (ages - 6)))
  r4 <- function(v) array(rep(v, each = n_regions * NS), dim = c(n_regions, NS, n_ages))
  r5 <- function(v) array(rep(v, each = n_regions * NS), dim = c(1, n_regions, NS, n_ages, n_flt))

  shared <- function(tm) list(
    t_spawn = t_spawn, n_seas = NS, seasdur = seasdur, spawn_seas = SPAWN,
    n_ages = n_ages, n_regions = n_regions,
    F_fract_flt = array(1 / NS, dim = c(n_regions, NS, n_flt)),
    dmr = array(0, dim = c(n_regions, NS, n_flt)),
    fish_sel = r5(sel), ret_sel = array(1, dim = c(1, n_regions, NS, n_ages, n_flt)),
    Movement = Mov1, Mrate = Mra1, move_timing = tm, do_recruits_move = 0,
    WAA = r4(waa), WAA_fish = r4(waa), MatAA = r4(mat),
    is_discard_fleet = array(0, dim = n_flt))

  msy_data <- function(tm) c(shared(tm), list(
    natmort = array(rep(M_r, times = n_ages), dim = c(n_regions, n_ages)),
    rec_seas_prop = array(1 / NS, dim = NS), sex_ratio_f = array(0.5, dim = n_regions),
    rec_region_prop = array(rec_prop, dim = n_regions),
    SR_type = 1, h = 0.7, R0 = 100, ln_global_R0 = log(100)))

  spr_data <- function(tm) c(shared(tm), list(
    n_pop = 1,
    natmort = array(rep(M_r, times = n_ages), dim = c(1, n_regions, n_ages)),
    WAA = array(r4(waa), c(1, n_regions, NS, n_ages)),
    MatAA = array(r4(mat), c(1, n_regions, NS, n_ages)),
    Movement = array(Mov1, c(1, n_regions, n_regions, NS, n_ages)),
    Mrate = array(Mra1, c(1, n_regions, n_regions, NS, n_ages)),
    sex_ratio_f = array(0.5, dim = c(1, n_regions)),
    rec_seas_prop = array(1 / NS, dim = c(1, NS)), stray_rate = array(0, dim = 1),
    sgl_seas_spawning_movement = array(rep(diag(n_regions), n_ages),
                                       dim = c(1, n_regions, n_regions, n_ages)),
    natal_region = 1, n_pop_in_region = array(1, dim = n_regions), SPR_x = 0.40,
    rec_region_prop = array(rec_prop, dim = c(1, n_regions))))

  ev <- function(fn, dat, parname, val) {
    p <- list(); p[[parname]] <- log(val)
    RTMB::MakeADFun(SPoRC:::cmb(fn, dat), parameters = p, random = NULL, silent = TRUE)$report(log(val))
  }

  for (tm in 0:2) {
    r_spr <- ev(SPoRC:::global_SPR,     spr_data(tm), "log_F_x",  Fv)
    r_msy <- ev(SPoRC:::global_BH_Fmsy, msy_data(tm), "log_Fmsy", Fv)
    expect_equal(r_msy$SBPR_0, sum(r_spr$SB0), tolerance = 1e-10,
                 label = sprintf("two-season unfished SBPR, MSY vs SPR, move_timing = %d", tm))
    expect_equal(r_msy$SBPR_F, sum(r_spr$SB), tolerance = 1e-10,
                 label = sprintf("two-season fished SBPR, MSY vs SPR, move_timing = %d", tm))
  }
})

# ---------------------------------------------------------------------------
# 2. Diffusion parameter recovery
# ---------------------------------------------------------------------------

test_that("CTMC diffusion parameter is recoverable from continuous movement fractions", {
  # Generate movement under a known diffusion parameter, then estimate it back by
  # minimising the discrepancy between observed and predicted movement fractions.
  # This checks that Get_Movement's generator and the seasdur scaling are mutually
  # consistent -- if the two disagreed, the recovered value would be biased.
  n_regions <- 3; n_ages <- 3; n_sexes <- 1; n_yrs <- 2; n_seas <- 2
  seasdur <- c(0.3, 0.7)
  dat <- expand.grid(pop = 1, regions = seq_len(n_regions), years = seq_len(n_yrs),
                     seas = seq_len(n_seas), ages = seq_len(n_ages), sexes = seq_len(n_sexes))
  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L

  get_move <- function(log_theta) {
    Get_Movement(
      move_type = 1, do_recruits_move = 1,
      n_pop = 1, n_regions = n_regions, n_yrs = n_yrs, n_proj_yrs_devs = 0,
      n_ages = n_ages, n_sexes = n_sexes, n_seas = n_seas,
      move_pars = array(0, c(1, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
      move_devs = array(0, c(1, n_regions, n_regions - 1, n_yrs, n_seas, n_ages, n_sexes)),
      use_fixed_movement = 0, Fixed_Movement = NULL,
      ctmc_move_dat = dat, preference_formula = ~ 1, diffusion_formula = ~ 1,
      log_move_diffusion_pars = log_theta, move_preference_pars = 0,
      area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
      seasdur = seasdur, ctmc_scale_by_seasdur = 1
    )
  }

  true_log_theta <- log(0.35)
  obs <- get_move(true_log_theta)$Movement

  # Least squares over the movement fractions
  nll <- function(lt) sum((get_move(lt)$Movement - obs)^2)
  fit <- stats::optimize(nll, interval = c(log(0.02), log(3)), tol = 1e-10)

  expect_equal(fit$minimum, true_log_theta, tolerance = 1e-3)
  expect_lt(fit$objective, 1e-10)
})

test_that("season-duration scaling is identifiable, not absorbed by the diffusion parameter", {
  # With unequal season durations, movement in a short season must differ from a long one.
  # If the seasdur scaling were dropped, no diffusion value could reproduce both.
  n_regions <- 3; n_seas <- 2
  seasdur <- c(0.25, 0.75)
  dat <- expand.grid(pop = 1, regions = seq_len(n_regions), years = 1,
                     seas = seq_len(n_seas), ages = 1:2, sexes = 1)
  adj <- matrix(1L, n_regions, n_regions); diag(adj) <- 0L

  res <- Get_Movement(
    move_type = 1, do_recruits_move = 1,
    n_pop = 1, n_regions = n_regions, n_yrs = 1, n_proj_yrs_devs = 0,
    n_ages = 2, n_sexes = 1, n_seas = n_seas,
    move_pars = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, 2, 1)),
    move_devs = array(0, c(1, n_regions, n_regions - 1, 1, n_seas, 2, 1)),
    use_fixed_movement = 0, Fixed_Movement = NULL,
    ctmc_move_dat = dat, preference_formula = ~ 1, diffusion_formula = ~ 1,
    log_move_diffusion_pars = log(0.4), move_preference_pars = 0,
    area_r = rep(1, n_regions), adjacency_mat = adj, ctmc_diffusion_bounds = 0,
    seasdur = seasdur, ctmc_scale_by_seasdur = 1
  )

  short <- res$Movement[1, , , 1, 1, 2, 1]
  long  <- res$Movement[1, , , 1, 2, 2, 1]
  # Residency should be higher in the short season
  expect_gt(mean(diag(short)), mean(diag(long)))
  # And composing the two must equal one full year of the generator
  Q <- res$Mrate[1, , , 1, 1, 2, 1]
  composed <- t(as.matrix(Matrix::expm(t(Q) * seasdur[1]))) %*% t(as.matrix(Matrix::expm(t(Q) * seasdur[2])))
  expect_equal(unname(composed), unname(t(as.matrix(Matrix::expm(t(Q))))),
               tolerance = 1e-9, ignore_attr = TRUE)
})

# ---------------------------------------------------------------------------
# 3. SSB reproduction
# ---------------------------------------------------------------------------

test_that("continuous movement reproduces SSB when stepped as one season or several", {
  # A population advanced through k sub-seasons under a constant generator and constant
  # annual mortality must land where a single full-year step lands. Spawning biomass is a
  # weighted sum of that state, so this is the SSB-recovery property in its sharpest form.
  set.seed(31)
  n <- 4
  D <- matrix(stats::runif(n * n, 0.05, 0.4), n, n)
  diag(D) <- 0; diag(D) <- -colSums(D)
  Q <- t(D)
  Z_annual <- c(0.3, 0.5, 0.2, 0.4)
  WAA <- c(1.2, 0.9, 1.5, 1.1); MatAA <- c(0.8, 0.6, 0.9, 0.7)
  N0 <- c(500, 300, 150, 250)

  for (k in c(2, 4, 12)) {
    dur <- rep(1 / k, k)
    N <- N0
    for (s in seq_len(k)) {
      Mv <- t(as.matrix(Matrix::expm(t(Q) * dur[s])))
      N <- advance_seas(N, Mv, Z_annual * dur[s], Q, dur[s], 2)
    }
    Mv1 <- t(as.matrix(Matrix::expm(t(Q))))
    one <- advance_seas(N0, Mv1, Z_annual, Q, 1, 2)

    expect_equal(N, one, tolerance = 1e-9, label = sprintf("%d sub-seasons vs one step", k))
    expect_equal(sum(N * WAA * MatAA), sum(one * WAA * MatAA), tolerance = 1e-9)
  }
})

test_that("spawning state is consistent between partial and full propagation", {
  # spawn_state at t_spawn then the remaining fraction must equal a full seasonal step,
  # which is what keeps SSB consistent with the numbers-at-age carried forward.
  set.seed(52)
  n <- 3
  D <- matrix(stats::runif(n * n, 0.05, 0.4), n, n)
  diag(D) <- 0; diag(D) <- -colSums(D)
  Q <- t(D)
  Mv <- t(as.matrix(Matrix::expm(D)))
  Z <- c(0.3, 0.15, 0.45)
  N0 <- c(400, 200, 100)
  t_spawn <- 0.4

  at_spawn <- spawn_state(N0, Mv, Z, Q, 1, t_spawn, move_timing = 2)
  # advancing the remaining (1 - t_spawn) of the season from the spawning state
  rest <- as.vector(as.matrix(Matrix::expm((t(Q) - diag(Z, n)) * (1 - t_spawn))) %*% at_spawn)
  full <- advance_seas(N0, Mv, Z, Q, 1, 2)

  expect_equal(rest, full, tolerance = 1e-9)
})
