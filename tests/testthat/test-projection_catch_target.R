library(SPoRC)
library(testthat)

# Tests for fmort_opt = "Catch" in Do_Population_Projection().
#
# The central test is a round trip. Run the projection at a known input F, take
# the realized catch as the target, and re-run under fmort_opt = "Catch". The
# solver has to recover the original F and reproduce every other projected
# quantity, which pins down the whole seasonal and spatial catch calculation
# rather than just checking that some F produces some catch.

# Builds a small but non-degenerate projection: multi-region, multi-season,
# multi-fleet, two sexes, partial retention on fleet 2 and non-zero discard
# mortality, and mildly mixing movement. Every array is filled by explicit index
# so nothing recycles down the wrong stride.
make_proj_inputs <- function(n_regions, n_seas, n_fish_fleets, n_proj_yrs = 6,
                             n_ages = 12, n_sexes = 2, move_timing = 0) {
  n_pop <- 1
  ages <- 1:n_ages

  natmort <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_ages, n_sexes))
  WAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  MatAA <- array(0, dim = dim(WAA))
  WAA_fish <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  fish_sel <- array(0, dim = dim(WAA_fish))
  ret_sel <- array(1, dim = dim(WAA_fish))

  for (r in 1:n_regions) for (y in 1:n_proj_yrs) for (s in 1:n_sexes) {
    natmort[1, r, y, , s] <- 0.1 + 0.01 * s
    for (seas in 1:n_seas) {
      WAA[1, r, y, seas, , s] <- (1 - exp(-0.2 * ages))^3 * (5 + s)
      MatAA[1, r, y, seas, , s] <- 1 / (1 + exp(-0.8 * (ages - 5)))
      for (f in 1:n_fish_fleets) {
        WAA_fish[1, r, y, seas, , s, f] <- WAA[1, r, y, seas, , s]
        fish_sel[1, r, y, seas, , s, f] <- 1 / (1 + exp(-1.1 * (ages - (3 + f) - 0.3 * r)))
        # fleet 2 discards small fish, so the retention and discard path is live
        if (f == 2) ret_sel[1, r, y, seas, , s, f] <- 1 / (1 + exp(-1.5 * (ages - 7)))
      }
    }
  }

  dmr <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  dmr[, , min(2, n_fish_fleets)] <- 0.4

  Movement <- array(0, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  Mrate <- array(0, dim = dim(Movement))
  if (n_regions > 1) {
    M <- matrix(0.12, n_regions, n_regions); diag(M) <- 0; diag(M) <- 1 - rowSums(M)
    Q <- matrix(0.15, n_regions, n_regions); diag(Q) <- 0; diag(Q) <- -rowSums(Q)
  } else {
    M <- matrix(1, 1, 1); Q <- matrix(0, 1, 1)
  }
  for (y in 1:n_proj_yrs) for (seas in 1:n_seas) for (a in 1:n_ages) for (s in 1:n_sexes) {
    Movement[1, , , y, seas, a, s] <- M
    Mrate[1, , , y, seas, a, s] <- Q
  }

  sgl_move <- array(0, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_ages, n_sexes))
  for (y in 1:n_proj_yrs) for (a in 1:n_ages) for (s in 1:n_sexes) sgl_move[1, , , y, a, s] <- diag(n_regions)

  terminal_NAA <- array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  for (r in 1:n_regions) for (seas in 1:n_seas) for (s in 1:n_sexes) {
    terminal_NAA[1, r, seas, , s] <- 1e6 * exp(-0.25 * (ages - 1)) / n_regions * (1 + 0.1 * r)
  }

  terminal_F <- array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  for (r in 1:n_regions) for (seas in 1:n_seas) for (f in 1:n_fish_fleets) {
    terminal_F[r, seas, f] <- 0.05 * (1 + 0.2 * r) * (1 + 0.1 * seas) / (n_seas * n_fish_fleets) * f
  }

  list(n_proj_yrs = n_proj_yrs, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages,
       n_sexes = n_sexes, n_fish_fleets = n_fish_fleets, n_seas = n_seas,
       seasdur = rep(1 / n_seas, n_seas), spawn_seas = 1, t_spawn = 0.2, natal_region = 1,
       stray_rate = array(0, dim = c(n_pop, n_proj_yrs)),
       sexratio = array(1 / n_sexes, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)),
       recruitment = array(2e5, dim = c(n_pop, n_regions, 20)),
       rec_seas_prop = array(1 / n_seas, dim = c(n_pop, n_seas)),
       terminal_NAA = terminal_NAA, terminal_NAA0 = terminal_NAA * 1.6,
       terminal_F = terminal_F, natmort = natmort, WAA = WAA, WAA_fish = WAA_fish,
       MatAA = MatAA, fish_sel = fish_sel, ret_sel = ret_sel, dmr = dmr,
       Movement = Movement, Mrate = Mrate, move_timing = move_timing,
       sgl_seas_spawning_movement = sgl_move, do_recruits_move = 0,
       recruitment_opt = "mean_rec")
}

PROJ_ARGS <- c("n_proj_yrs", "n_pop", "n_regions", "n_ages", "n_sexes", "sexratio",
               "n_fish_fleets", "do_recruits_move", "recruitment", "terminal_NAA",
               "terminal_NAA0", "terminal_F", "dmr", "natmort", "natal_region", "WAA",
               "WAA_fish", "MatAA", "fish_sel", "ret_sel", "Movement",
               "sgl_seas_spawning_movement", "stray_rate", "recruitment_opt", "t_spawn",
               "n_seas", "seasdur", "spawn_seas", "rec_seas_prop", "Mrate", "move_timing")

run_proj <- function(inp, fmort_opt, ...) {
  over <- list(...)
  base <- inp[setdiff(PROJ_ARGS, names(over))] # anything passed in ... wins
  do.call(Do_Population_Projection, c(base, list(fmort_opt = fmort_opt), over))
}

# regional annual catch, summed over populations, seasons and fleets
annual_catch <- function(out, n_regions, n_proj_yrs) {
  m <- array(0, dim = c(n_regions, n_proj_yrs))
  for (r in 1:n_regions) for (y in 1:n_proj_yrs) m[r, y] <- sum(out$proj_Catch[, r, y, , ])
  m
}

# regional and seasonal catch
seasonal_catch <- function(out, n_regions, n_proj_yrs, n_seas) {
  a <- array(0, dim = c(n_regions, n_proj_yrs, n_seas))
  for (r in 1:n_regions) for (y in 1:n_proj_yrs) for (s in 1:n_seas) {
    a[r, y, s] <- sum(out$proj_Catch[, r, y, s, ])
  }
  a
}

threshold_hcr <- function(x, frp, brp) frp * min(1, x / brp)


test_that("annual catch targets recover the fishing mortality that produced them", {

  # Each configuration runs at a known input F, has its realized catch fed back
  # as a target, and must return the same F and the same projection.
  cases <- list(
    list(lab = "1 region, 1 season, 1 fleet",              nr = 1, ns = 1, nf = 1, mt = 0),
    list(lab = "1 region, 4 seasons, 2 fleets",            nr = 1, ns = 4, nf = 2, mt = 0),
    list(lab = "3 regions, 1 season, 1 fleet",             nr = 3, ns = 1, nf = 1, mt = 0),
    list(lab = "3 regions, 4 seasons, 2 fleets, timing 0", nr = 3, ns = 4, nf = 2, mt = 0),
    list(lab = "3 regions, 4 seasons, 2 fleets, timing 1", nr = 3, ns = 4, nf = 2, mt = 1),
    list(lab = "3 regions, 4 seasons, 2 fleets, timing 2", nr = 3, ns = 4, nf = 2, mt = 2)
  )

  for (cs in cases) {
    inp <- make_proj_inputs(cs$nr, cs$ns, cs$nf, move_timing = cs$mt)
    npy <- inp$n_proj_yrs

    f_ref_pt <- array(0, dim = c(cs$nr, npy))
    for (r in 1:cs$nr) f_ref_pt[r, ] <- seq(0.06, 0.22, length.out = npy) * (1 + 0.25 * r)

    ref <- run_proj(inp, "Input", f_ref_pt = f_ref_pt)

    ci <- array(NA_real_, dim = c(cs$nr, npy))
    ci[, 2:npy] <- annual_catch(ref, cs$nr, npy)[, 2:npy]
    got <- run_proj(inp, "Catch", catch_input = ci)

    yy <- 2:npy
    expect_equal(got$proj_F[, yy, drop = FALSE], ref$proj_F[, yy, drop = FALSE],
                 tolerance = 1e-5, info = cs$lab)
    expect_equal(got$proj_Catch, ref$proj_Catch, tolerance = 1e-5, info = cs$lab)
    expect_equal(got$proj_SSB, ref$proj_SSB, tolerance = 1e-5, info = cs$lab)
    expect_equal(got$proj_NAA, ref$proj_NAA, tolerance = 1e-5, info = cs$lab)
    expect_lt(max(abs(got$proj_catch_resid[, yy, drop = FALSE])), 1e-5)
  }
})


test_that("catch targets are met under Beverton-Holt recruitment, including rec_lag = 0", {

  # rec_lag = 0 is the hard case: the F being solved for changes SSB at
  # spawn_seas, which changes that same year's recruitment, which changes catch.
  for (rec_lag in c(1, 0)) {

    n_regions <- 2; n_seas <- 4; n_fish_fleets <- 1; spawn_seas <- 2
    inp <- make_proj_inputs(n_regions, n_seas, n_fish_fleets)
    npy <- inp$n_proj_yrs; n_pop <- 1; n_ages <- inp$n_ages

    # rec_lag = 0 requires no recruitment entering before the spawning season
    rsp <- array(0, dim = c(n_pop, n_seas))
    rsp[1, spawn_seas:n_seas] <- 1 / length(spawn_seas:n_seas)
    inp$rec_seas_prop <- rsp
    inp$spawn_seas <- spawn_seas

    bh <- list(
      rec_dd = 0, rec_lag = rec_lag, R0 = rep(3e5, n_pop),
      h = array(0.7, dim = c(n_pop, n_regions)),
      rec_region_prop = array(1 / n_regions, dim = c(n_pop, n_regions)),
      rec_seas_prop = rsp,
      SSB = array(4e6, dim = c(n_pop, n_regions, 20)),
      WAA = array(inp$WAA[, , 1, , , 1], dim = c(n_pop, n_regions, n_seas, n_ages)),
      MatAA = array(inp$MatAA[, , 1, , , 1], dim = c(n_pop, n_regions, n_seas, n_ages)),
      natmort = array(inp$natmort[, , 1, , 1], dim = c(n_pop, n_regions, n_ages)),
      Movement = array(inp$Movement[, , , 1, , , 1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
      Mrate = array(inp$Mrate[, , , 1, , , 1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
      sgl_seas_spawning_movement = array(inp$sgl_seas_spawning_movement[, , , 1, , 1],
                                         dim = c(n_pop, n_regions, n_regions, n_ages)),
      stray_rate = rep(0, n_pop), init_F = inp$terminal_F,
      fish_sel = array(inp$fish_sel[, , 1, , , 1, ], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
      ret_sel = array(inp$ret_sel[, , 1, , , 1, ], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
      dmr = inp$dmr, sex_ratio_f = array(0.5, dim = c(n_pop, n_regions))
    )

    f_ref_pt <- array(0, dim = c(n_regions, npy))
    for (r in 1:n_regions) f_ref_pt[r, ] <- seq(0.07, 0.20, length.out = npy) * (1 + 0.2 * r)

    ref <- run_proj(inp, "Input", recruitment_opt = "bh_rec", bh_rec_opt = bh, f_ref_pt = f_ref_pt)
    ci <- array(NA_real_, dim = c(n_regions, npy))
    ci[, 2:npy] <- annual_catch(ref, n_regions, npy)[, 2:npy]
    got <- run_proj(inp, "Catch", recruitment_opt = "bh_rec", bh_rec_opt = bh, catch_input = ci)

    expect_equal(got$proj_F[, 2:npy], ref$proj_F[, 2:npy], tolerance = 1e-5,
                 info = paste("rec_lag =", rec_lag))
    expect_equal(got$proj_SSB, ref$proj_SSB, tolerance = 1e-5, info = paste("rec_lag =", rec_lag))
  }
})


test_that("seasonal catch targets hit a profile the terminal year's seasonal shares cannot", {

  n_regions <- 3; n_seas <- 4; n_fish_fleets <- 2

  for (mt in c(0, 2)) {
    inp <- make_proj_inputs(n_regions, n_seas, n_fish_fleets, move_timing = mt)
    npy <- inp$n_proj_yrs

    ref <- run_proj(inp, "Input", f_ref_pt = array(0.15, dim = c(n_regions, npy)))

    # tilt the seasonal profile hard away from the terminal year shares, which an
    # annual F spread at those shares could not reproduce
    ci <- array(NA_real_, dim = c(n_regions, npy, n_seas))
    ci[, 2:npy, ] <- seasonal_catch(ref, n_regions, npy, n_seas)[, 2:npy, ]
    tilt <- seq(1.8, 0.35, length.out = n_seas)
    for (s in 1:n_seas) ci[, 2:npy, s] <- ci[, 2:npy, s] * tilt[s]

    got <- run_proj(inp, "Catch", catch_input = ci)

    expect_equal(seasonal_catch(got, n_regions, npy, n_seas)[, 2:npy, ],
                 ci[, 2:npy, ], tolerance = 1e-5, info = paste("move_timing", mt))

    # and the solved seasonal profile must actually have left the terminal one
    shr_term <- apply(inp$terminal_F, c(1, 2), sum)
    shr_term <- shr_term / rowSums(shr_term)
    shr_got <- array(got$proj_F_seas[, 3, ], dim = c(n_regions, n_seas))
    shr_got <- shr_got / rowSums(shr_got)
    expect_gt(max(abs(shr_got - shr_term)), 0.05)
  }
})


test_that("seasonal targets matching the terminal shares reproduce the annual solve", {

  inp <- make_proj_inputs(3, 4, 2)
  npy <- inp$n_proj_yrs
  f_ref_pt <- array(0, dim = c(3, npy))
  for (r in 1:3) f_ref_pt[r, ] <- seq(0.06, 0.22, length.out = npy) * (1 + 0.25 * r)

  ref <- run_proj(inp, "Input", f_ref_pt = f_ref_pt)
  ci <- array(NA_real_, dim = c(3, npy, 4))
  ci[, 2:npy, ] <- seasonal_catch(ref, 3, npy, 4)[, 2:npy, ]
  got <- run_proj(inp, "Catch", catch_input = ci)

  expect_equal(got$proj_F[, 2:npy], ref$proj_F[, 2:npy], tolerance = 1e-5)
  expect_equal(got$proj_SSB, ref$proj_SSB, tolerance = 1e-5)
  expect_equal(dim(got$proj_catch_resid), dim(ci)) # resid is shaped like catch_input
})


test_that("years left NA fall back to the harvest control rule", {

  # Splice test: targets for years 2 and 3 taken off a pure fallback run, so a
  # correctly wired mixed run must reproduce that run everywhere.
  for (fallback in c("HCR", "Input", "HCR_global")) {
    for (seasonal in c(FALSE, TRUE)) {

      inp <- make_proj_inputs(3, 4, 2)
      npy <- inp$n_proj_yrs
      extra <- list(f_ref_pt = array(0.13, dim = c(3, npy)),
                    b_ref_pt = array(3e6, dim = c(inp$n_pop, 3, npy)),
                    HCR_function = threshold_hcr)

      ref <- do.call(run_proj, c(list(inp, fallback), extra))

      ci <- if (seasonal) array(NA_real_, dim = c(3, npy, 4)) else array(NA_real_, dim = c(3, npy))
      if (seasonal) {
        ci[, 2:3, ] <- seasonal_catch(ref, 3, npy, 4)[, 2:3, ]
      } else {
        ci[, 2:3] <- annual_catch(ref, 3, npy)[, 2:3]
      }
      got <- do.call(run_proj, c(list(inp, "Catch"), extra,
                                 list(catch_input = ci, catch_fallback_opt = fallback)))

      lab <- paste(fallback, if (seasonal) "seasonal" else "annual")
      expect_equal(got$proj_F[, 2:npy], ref$proj_F[, 2:npy], tolerance = 1e-5, info = lab)
      expect_equal(got$proj_SSB, ref$proj_SSB, tolerance = 1e-5, info = lab)
      expect_equal(got$proj_NAA, ref$proj_NAA, tolerance = 1e-5, info = lab)

      # residuals are recorded only where a target was set
      rs <- array(if (seasonal) got$proj_catch_resid[, , 1] else got$proj_catch_resid,
                  dim = c(3, npy))
      expect_false(any(is.na(rs[, 2:3])))
      expect_true(all(is.na(rs[, c(1, 4:npy)])))
    }
  }
})


test_that("the fallback rule responds to the stock the constrained years produced", {

  inp <- make_proj_inputs(3, 4, 2)
  npy <- inp$n_proj_yrs
  f_ref_pt <- array(0.13, dim = c(3, npy))
  b_ref_pt <- array(3e6, dim = c(inp$n_pop, 3, npy))
  extra <- list(f_ref_pt = f_ref_pt, b_ref_pt = b_ref_pt, HCR_function = threshold_hcr)

  ref <- do.call(run_proj, c(list(inp, "HCR"), extra))
  ci <- array(NA_real_, dim = c(3, npy))
  ci[, 2:3] <- annual_catch(ref, 3, npy)[, 2:3] * 0.5 # halve the constrained years
  got <- do.call(run_proj, c(list(inp, "Catch"), extra,
                             list(catch_input = ci, catch_fallback_opt = "HCR")))

  # constrained years take the halved catch
  expect_equal(annual_catch(got, 3, npy)[, 2:3], ci[, 2:3], tolerance = 1e-5)

  # fishing less leaves more fish, so the later HCR years differ from the pure run
  expect_gt(max(abs(got$proj_F[, 4:npy] - ref$proj_F[, 4:npy])), 1e-6)
  expect_true(all(got$proj_SSB[1, , 4] > ref$proj_SSB[1, , 4]))

  # and those years are still exactly the HCR applied to the mixed run's own SSB
  expect_equal(got$proj_F[, 5],
               sapply(1:3, function(r) threshold_hcr(sum(got$proj_SSB[, r, 4]), 0.13, 3e6)),
               tolerance = 1e-10)
})


test_that("catch_terminal_yr solves projection year 1 against its own target", {

  for (seasonal in c(FALSE, TRUE)) {
    inp <- make_proj_inputs(3, 4, 2)
    npy <- inp$n_proj_yrs
    ref <- run_proj(inp, "Input", f_ref_pt = array(0.15, dim = c(3, npy)))

    ci <- if (seasonal) seasonal_catch(ref, 3, npy, 4) else annual_catch(ref, 3, npy)
    got <- run_proj(inp, "Catch", catch_input = ci, catch_terminal_yr = TRUE)

    # year 1's catch was produced at terminal_F, so the solve must return it
    expect_equal(got$proj_F[, 1], rowSums(inp$terminal_F), tolerance = 1e-6,
                 info = if (seasonal) "seasonal" else "annual")
    expect_equal(got$proj_Catch, ref$proj_Catch, tolerance = 1e-5)
    expect_equal(got$proj_NAA, ref$proj_NAA, tolerance = 1e-5)
  }

  # without it, year 1 is left at terminal_F and its target ignored
  inp <- make_proj_inputs(1, 1, 1)
  npy <- inp$n_proj_yrs
  ci <- array(NA_real_, dim = c(1, npy)); ci[1, ] <- 3e5
  got <- run_proj(inp, "Catch", catch_input = ci)
  expect_equal(got$proj_F[1, 1], sum(inp$terminal_F), tolerance = 1e-12)
  expect_true(is.na(got$proj_catch_resid[1, 1]))
})


test_that("a zero target means no fishing and is distinct from NA", {

  inp <- make_proj_inputs(3, 4, 2)
  npy <- inp$n_proj_yrs

  # region 2 asked for zero catch, regions 1 and 3 for a real catch
  ci <- array(NA_real_, dim = c(3, npy))
  ci[, 2:npy] <- 0
  ci[c(1, 3), 2:npy] <- 4e5
  got <- run_proj(inp, "Catch", catch_input = ci)

  expect_true(all(got$proj_F[2, 2:npy] == 0))
  expect_true(all(annual_catch(got, 3, npy)[2, 2:npy] == 0))
  expect_lt(max(abs(got$proj_catch_resid[c(1, 3), 2:npy])), 1e-5)

  # a year of all zeros is a real instruction, not a fallback
  ci2 <- array(NA_real_, dim = c(3, npy)); ci2[, 2] <- 0
  got2 <- run_proj(inp, "Catch", catch_input = ci2,
                   f_ref_pt = array(0.13, dim = c(3, npy)),
                   b_ref_pt = array(3e6, dim = c(inp$n_pop, 3, npy)),
                   HCR_function = threshold_hcr)
  expect_true(all(got2$proj_F[, 2] == 0))
  expect_true(all(got2$proj_F[, 3] > 0))
})


test_that("an unreachable catch target warns, caps F, and records the shortfall", {

  inp <- make_proj_inputs(3, 4, 2)
  npy <- inp$n_proj_yrs
  ci <- array(NA_real_, dim = c(3, npy)); ci[, 2:npy] <- 1e12

  # every targeted year should warn, so capture them all rather than just the first
  w <- capture_warnings(got <- run_proj(inp, "Catch", catch_input = ci, catch_f_max = 5))
  expect_length(w, npy - 1L)
  expect_true(all(grepl("not reachable", w)))
  expect_true(all(grepl("region\\(s\\) 1, 2, 3", w)))

  expect_equal(max(abs(got$proj_F[, 2:npy] - 5)), 0, tolerance = 1e-9)
  expect_true(all(got$proj_catch_resid[, 2:npy] < -0.9)) # target massively undershot
})


test_that("catch_input is validated before the projection runs", {

  inp <- make_proj_inputs(3, 4, 2)
  npy <- inp$n_proj_yrs
  ok_ci <- array(NA_real_, dim = c(3, npy)); ok_ci[, 2:3] <- 3e5
  extra <- list(f_ref_pt = array(0.13, dim = c(3, npy)),
                b_ref_pt = array(3e6, dim = c(inp$n_pop, 3, npy)),
                HCR_function = threshold_hcr)

  expect_error(run_proj(inp, "Catch"), "requires catch_input")
  expect_error(run_proj(inp, "Catch", catch_input = array(1, dim = c(2, npy))), "dimensioned")
  expect_error(run_proj(inp, "Catch", catch_input = array(-1, dim = c(3, npy))), "negative")
  expect_error(run_proj(inp, "Catch", catch_input = array(NA_real_, dim = c(3, npy))),
               "no projection year carries")
  expect_error(do.call(run_proj, c(list(inp, "Catch"), extra,
                                   list(catch_input = ok_ci, catch_fallback_opt = "bogus"))),
               "fallback options")

  # a year has to be all target or all NA
  part <- array(NA_real_, dim = c(3, npy)); part[1, 2] <- 1e5
  expect_error(do.call(run_proj, c(list(inp, "Catch"), extra, list(catch_input = part))),
               "only partly specified")

  # falling back without the inputs the fallback needs
  expect_error(run_proj(inp, "Catch", catch_input = ok_ci, catch_fallback_opt = "HCR"),
               "needs f_ref_pt")

  # a season the terminal year never fished cannot be given a target
  inp2 <- make_proj_inputs(3, 4, 2)
  inp2$terminal_F[, 2, ] <- 0
  ref2 <- run_proj(inp2, "Input", f_ref_pt = array(0.15, dim = c(3, npy)))
  ci2 <- seasonal_catch(ref2, 3, npy, 4)
  ci2[, , 2] <- 1e4
  expect_error(run_proj(inp2, "Catch", catch_input = ci2), "no fleet split")

  # a bare vector is accepted for a single-region model
  inp3 <- make_proj_inputs(1, 1, 1)
  got <- run_proj(inp3, "Catch", catch_input = c(NA, rep(3e5, npy - 1)))
  expect_lt(max(abs(got$proj_catch_resid[, -1])), 1e-5)
})


test_that("other fmort_opt settings are unaffected by the catch machinery", {

  # A catch projection must not change what the existing options do. These run
  # without any catch arguments at all and simply have to work.
  inp <- make_proj_inputs(3, 4, 2)
  npy <- inp$n_proj_yrs
  extra <- list(f_ref_pt = array(0.12, dim = c(3, npy)),
                b_ref_pt = array(2e6, dim = c(inp$n_pop, 3, npy)),
                HCR_function = threshold_hcr)

  for (opt in c("HCR", "Input", "HCR_global")) {
    got <- do.call(run_proj, c(list(inp, opt), extra))
    expect_true(all(is.na(got$proj_catch_resid))) # no targets, so no residuals
    expect_true(all(got$proj_F[, 2:npy] > 0))
    # proj_F_seas is a view of proj_F, for every year including the trailing one
    expect_equal(max(abs(sapply(1:(npy + 1), function(y)
      rowSums(got$proj_F_seas[, y, ]) - got$proj_F[, y]))), 0, tolerance = 1e-12)
  }

  # Input F is a pure scaling, so doubling it must not change the seasonal shares
  a <- do.call(run_proj, c(list(inp, "Input"), extra))
  b <- do.call(run_proj, c(list(inp, "Input"),
                           utils::modifyList(extra, list(f_ref_pt = array(0.24, dim = c(3, npy))))))
  expect_equal(b$proj_F[, 3] / a$proj_F[, 3], rep(2, 3), tolerance = 1e-12)
})


test_that("projected total biomass continues the estimated series and exceeds SSB", {

  inp <- make_proj_inputs(3, 4, 2)
  npy <- inp$n_proj_yrs
  got <- run_proj(inp, "Input", f_ref_pt = array(0.12, dim = c(3, npy)))

  expect_equal(dim(got$proj_Total_Biom), c(inp$n_pop, 3L, npy))

  # total biomass covers all ages and both sexes, so it must exceed female
  # mature biomass everywhere
  expect_true(all(got$proj_Total_Biom > got$proj_SSB))

  # and it must respond to fishing in the right direction
  low <- run_proj(inp, "Input", f_ref_pt = array(0.02, dim = c(3, npy)))
  expect_true(all(low$proj_Total_Biom[, , npy] > got$proj_Total_Biom[, , npy]))

  # Check the definition directly rather than by a monotonicity heuristic: total
  # biomass is numbers times weight over all ages and sexes at the spawning
  # point, and SSB is the mature-female analog over the same state. A
  # single-region, single-season model is used so that the spawning-point state
  # is exactly proj_NAA discounted by t_spawn.
  inp1 <- make_proj_inputs(1, 1, 1)
  one <- run_proj(inp1, "Input", f_ref_pt = array(0.08, dim = c(1, npy)))

  for (y in 1:npy) {
    disc <- exp(-one$proj_ZAA[1, 1, y, 1, , ] * inp1$t_spawn)
    expect_equal(one$proj_Total_Biom[1, 1, y],
                 sum(one$proj_NAA[1, 1, y, 1, , ] * inp1$WAA[1, 1, y, 1, , ] * disc),
                 tolerance = 1e-10, info = paste("total biomass, year", y))
    expect_equal(one$proj_SSB[1, 1, y],
                 sum(one$proj_NAA[1, 1, y, 1, , 1] * inp1$WAA[1, 1, y, 1, , 1] *
                       inp1$MatAA[1, 1, y, 1, , 1] * disc[, 1]),
                 tolerance = 1e-10, info = paste("SSB, year", y))
  }
})


test_that("returned arrays carry the documented dimensions", {

  n_regions <- 3; n_seas <- 4; n_fish_fleets <- 2
  inp <- make_proj_inputs(n_regions, n_seas, n_fish_fleets)
  npy <- inp$n_proj_yrs; n_pop <- inp$n_pop; n_ages <- inp$n_ages; n_sexes <- inp$n_sexes

  ci <- array(NA_real_, dim = c(n_regions, npy)); ci[, 2:3] <- 3e5
  got <- run_proj(inp, "Catch", catch_input = ci, catch_fallback_opt = "HCR",
                  f_ref_pt = array(0.1, dim = c(n_regions, npy)),
                  b_ref_pt = array(2e6, dim = c(n_pop, n_regions, npy)),
                  HCR_function = threshold_hcr)

  expect_equal(dim(got$proj_F), c(n_regions, npy + 1L))
  expect_equal(dim(got$proj_F_seas), c(n_regions, npy + 1L, n_seas))
  expect_equal(dim(got$proj_ret_FAA), c(n_pop, n_regions, npy + 1L, n_seas, n_ages, n_sexes, n_fish_fleets))
  expect_equal(dim(got$proj_disc_FAA), dim(got$proj_ret_FAA))
  expect_equal(dim(got$proj_Catch), c(n_pop, n_regions, npy, n_seas, n_fish_fleets))
  expect_equal(dim(got$proj_SSB), c(n_pop, n_regions, npy))
  expect_equal(dim(got$proj_eff_SSB), c(n_pop, npy))
  expect_equal(dim(got$proj_Total_Biom), c(n_pop, n_regions, npy))
  expect_equal(dim(got$proj_Dynamic_SSB0), c(n_pop, n_regions, npy))
  expect_equal(dim(got$proj_NAA), c(n_pop, n_regions, npy + 1L, n_seas, n_ages, n_sexes))
  expect_equal(dim(got$proj_NAA0), dim(got$proj_NAA))
  expect_equal(dim(got$proj_ZAA), c(n_pop, n_regions, npy + 1L, n_seas, n_ages, n_sexes))
  expect_equal(dim(got$proj_catch_resid), c(n_regions, npy))

  # documented trailing-slot behavior
  expect_true(any(got$proj_NAA[, , npy + 1, , , ] > 0))   # filled
  expect_true(all(got$proj_ZAA[, , npy + 1, , , ] == 0))  # left at 0
  expect_true(all(got$proj_ret_FAA[, , npy + 1, , , , ] == 0))

  # seasonal targets give a seasonal residual array
  ci3 <- array(NA_real_, dim = c(n_regions, npy, n_seas)); ci3[, 2:3, ] <- 8e4
  got3 <- run_proj(inp, "Catch", catch_input = ci3, catch_fallback_opt = "HCR",
                   f_ref_pt = array(0.1, dim = c(n_regions, npy)),
                   b_ref_pt = array(2e6, dim = c(n_pop, n_regions, npy)),
                   HCR_function = threshold_hcr)
  expect_equal(dim(got3$proj_catch_resid), dim(ci3))
})
