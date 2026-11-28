# An unfished model with two populations that home to their natal region to spawn stays at equilibrium in every
# year: spawning biomass does not change and recruitment equals R0 (1e-8). Checked with spawning at the start of
# the year and part way through, natural mortality the same or different by region, and each movement timing.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

# two populations, one born in each region, under Beverton-Holt recruitment by population and no catch
homing_equilibrium_gaps <- function(t_spawn, natmort_by_region, homing, move_timing = 0) {

  n_pop <- 2
  n_regions <- 2
  n_yrs <- 25
  n_ages <- 7

  natmort <- array(0, dim = c(n_pop, n_regions, n_yrs, n_ages, 1))
  for(r in 1:n_regions) natmort[,r,,,] <- natmort_by_region[r] # per year

  # spawning movement [pop, from, to, year, age, sex]: a share of the fish outside the natal region return to spawn
  spawn_move <- array(0, dim = c(n_pop, n_regions, n_regions, n_yrs, n_ages, 1))
  for(p in 1:n_pop) {
    home <- p
    away <- 3 - p
    spawn_move[p,home,home,,,] <- 1
    spawn_move[p,away,home,,,] <- homing
    spawn_move[p,away,away,,,] <- 1 - homing
  } # end p loop

  il <- suppressWarnings(suppressMessages(sweep_input(
    dims = list(n_pop = n_pop, natal_region = c(1, 2), n_regions = n_regions, n_sexes = 1, n_fish_fleets = 1, n_yrs = n_yrs, n_ages = n_ages),
    rec = list(
      rec_model = "bh_rec",
      rec_dd = "local", # one curve per population
      h_spec = "fix",
      init_age_strc = 2, # the equilibrium with movement at every age
      t_spawn = t_spawn,
      sgl_seas_spawning_movement = spawn_move
    ),
    biol = list(Fixed_natmort = natmort),
    move = list(move_timing = move_timing),
    catch = list(UseCatch = array(0, dim = c(n_regions, n_yrs, 1, 1))) # no catch fit, so F is zero
  )))
  rep <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))$rep

  # refuse a run that would make the comparisons below trivial
  if(max(rep$Fmort) != 0) stop("the run is fished")

  # spawning biomass and recruitment by population and year, summed over regions
  SSB <- apply(rep$SSB, c(1, 3), sum)
  Rec <- apply(rep$Rec, c(1, 3), sum)

  c(SSB = max(abs(SSB / SSB[,1] - 1)), # largest change in spawning biomass from year 1
    Rec = max(abs(Rec / as.numeric(rep$R0) - 1))) # largest gap between recruitment and R0
}

# Tests ----------------------------------------------------------------------

test_that("an unfished natal homing model stays at equilibrium when the homing or the spawn timing cannot matter", {

  # spawning at the start of the year
  gaps <- homing_equilibrium_gaps(t_spawn = 0, natmort_by_region = c(0.15, 0.35), homing = 0.6)
  expect_lt(max(gaps), 1e-8)

  # natural mortality the same in both regions
  gaps <- homing_equilibrium_gaps(t_spawn = 0.4, natmort_by_region = c(0.25, 0.25), homing = 0.6)
  expect_lt(max(gaps), 1e-8)

  # no fish return to spawn
  gaps <- homing_equilibrium_gaps(t_spawn = 0.4, natmort_by_region = c(0.15, 0.35), homing = 0)
  expect_lt(max(gaps), 1e-8)

})

test_that("an unfished natal homing model stays at equilibrium with mid-year spawning and mortality that differs by region", {

  # movement before mortality
  gaps <- homing_equilibrium_gaps(t_spawn = 0.4, natmort_by_region = c(0.15, 0.35), homing = 0.6, move_timing = 0)
  expect_lt(gaps[["SSB"]], 1e-8)
  expect_lt(gaps[["Rec"]], 1e-8)

  # mortality before movement
  gaps <- homing_equilibrium_gaps(t_spawn = 0.4, natmort_by_region = c(0.15, 0.35), homing = 0.6, move_timing = 1)
  expect_lt(gaps[["SSB"]], 1e-8)
  expect_lt(gaps[["Rec"]], 1e-8)

})
