library(SPoRC)
library(testthat)

test_that("Do_Population_Projection handles age-0 (rec_lag = 0) recruitment with spawning after season 1", {

  # Minimal single-population, single-region, 2-season model with
  # spawn_seas = 2 -- recruits can only enter in the spawning season itself,
  # since there's no season after it. Projects forward under F = 0 and checks
  # (a) no recruits ever appear pre-spawn (season 1), and (b) SSB settles to
  # a stable unfished equilibrium, confirming the projection's age-0 timing
  # logic runs correctly and converges.

  n_pop <- 1; n_regions <- 1; n_ages <- 6; n_sexes <- 1
  n_fish_fleets <- 1; n_seas <- 2; spawn_seas <- 2
  n_proj_yrs <- 300

  seasdur <- c(0.5, 0.5)
  t_spawn <- 0
  do_recruits_move <- 0

  M_val <- 0.2
  WAA_vec <- c(0.01, 0.05, 0.15, 0.30, 0.45, 0.55)
  MatAA_vec <- c(0, 0, 1, 1, 1, 1) # age-0 and age-1 immature (required for rec_lag = 0)

  natmort <- array(M_val, dim = c(n_pop, n_regions, n_proj_yrs, n_ages, n_sexes))
  WAA <- array(rep(WAA_vec, each = n_pop * n_regions * n_proj_yrs * n_seas),
              dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  WAA_fish <- array(WAA_vec[1], dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  MatAA <- array(rep(MatAA_vec, each = n_pop * n_regions * n_proj_yrs * n_seas),
                dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  fish_sel <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # no fishing
  ret_sel <- array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  Movement <- array(1, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  sexratio <- array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes))

  # rec_seas_prop must be 0 before spawn_seas -- all recruits enter spawn_seas
  rec_seas_prop <- array(0, dim = c(n_pop, n_seas))
  rec_seas_prop[, spawn_seas] <- 1

  R0 <- 1000
  h <- 0.7

  # Arbitrary starting age structure with age index 1 (age-0) empty, matching
  # the invariant maintained by the rec_lag = 0 population loop. terminal_NAA
  # is season-specific, so both seasons of the terminal year need values.
  terminal_NAA <- array(0, dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  terminal_NAA[1,1,1,2:n_ages,1] <- c(800, 600, 400, 250, 600)
  terminal_NAA[1,1,2,2:n_ages,1] <- c(700, 550, 380, 240, 590)
  terminal_NAA0 <- terminal_NAA

  srr_opt <- list(
    rec_dd = 1,
    rec_lag = 0,
    do_recruits_move = do_recruits_move,
    R0 = R0,
    h = array(h, dim = c(n_pop, n_regions)),
    rec_region_prop = array(1, dim = c(n_pop, n_regions)),
    WAA = array(WAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    MatAA = array(MatAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    SSB = array(0, dim = c(n_pop, n_regions, 1)), # no assessment history; projection starts from year 1
    Movement = array(Movement[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
    sex_ratio_f = array(1, dim = c(n_pop, n_regions)),
    sgl_seas_spawning_movement = NULL,
    stray_rate = array(0, dim = c(n_pop)),
    natmort = array(natmort[,,1,,1], dim = c(n_pop, n_regions, n_ages)),
    fish_sel = array(fish_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    ret_sel = array(ret_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  )

  out <- Do_Population_Projection(n_proj_yrs = n_proj_yrs,
                                  n_pop = n_pop,
                                  n_regions = n_regions,
                                  n_ages = n_ages,
                                  n_sexes = n_sexes,
                                  sexratio = sexratio,
                                  n_fish_fleets = n_fish_fleets,
                                  do_recruits_move = do_recruits_move,
                                  recruitment = array(0, dim = c(n_pop, n_regions, 1)), # unused for bh_rec
                                  terminal_NAA = terminal_NAA,
                                  terminal_NAA0 = terminal_NAA0,
                                  # nonzero just to keep the seasonal F-ratio calc (terminal_F / sum(terminal_F))
                                  # well-defined; fish_sel = 0 below means it contributes no actual mortality
                                  terminal_F = array(0.1, dim = c(n_regions, n_seas, n_fish_fleets)),
                                  dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
                                  natmort = natmort,
                                  WAA = WAA,
                                  WAA_fish = WAA_fish,
                                  MatAA = MatAA,
                                  fish_sel = fish_sel,
                                  ret_sel = ret_sel,
                                  Movement = Movement,
                                  recruitment_opt = "bh_rec",
                                  fmort_opt = "Input",
                                  f_ref_pt = array(0, dim = c(n_regions, n_proj_yrs)), # F = 0 throughout
                                  t_spawn = t_spawn,
                                  n_seas = n_seas,
                                  seasdur = seasdur,
                                  spawn_seas = spawn_seas,
                                  natal_region = 1,
                                  srr_opt = srr_opt
  )

  # No recruits ever appear pre-spawn (season 1, age index 1)
  expect_true(all(out$proj_NAA[1,1,,1,1,1] == 0))
  expect_true(all(out$proj_NAA0[1,1,,1,1,1] == 0))

  # Recruits do appear in spawn_seas from year 2 onward (year 1 carries the
  # terminal state forward with no new recruitment event)
  expect_true(all(out$proj_NAA[1,1,2:n_proj_yrs,spawn_seas,1,1] > 0))

  # Under F = 0, SSB should settle to a stable unfished equilibrium
  terminal_ssb <- out$proj_SSB[1,1,n_proj_yrs]
  penultimate_ssb <- out$proj_SSB[1,1,n_proj_yrs - 1]
  expect_equal(terminal_ssb, penultimate_ssb, tolerance = 1e-6)
  expect_true(terminal_ssb > 0)

})
