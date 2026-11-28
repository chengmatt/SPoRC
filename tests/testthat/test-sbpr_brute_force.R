# Spawning biomass per recruit against a cohort followed for 800 years, with no plus group formula (1e-10).
# Two regions, every movement timing, recruits moving or not, local and global density dependence. Mortality,
# movement, selectivity and retention differ at every age, so the oldest true age and the plus group cannot be swapped.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

sbpr_hand <- list(
  n_regions = 2,
  n_ages = 4, # age 4 is the plus group
  natmort = rbind(c(0.40, 0.30, 0.22, 0.15), # per year [region, age]
                  c(0.35, 0.28, 0.18, 0.10)),
  WAA = rbind(c(0.4, 1.0, 1.7, 2.5), # weight at age (kg) [region, age]
              c(0.5, 1.2, 1.9, 2.2)),
  MatAA = rbind(c(0, 0.4, 0.8, 1), # maturity [region, age]
                c(0.1, 0.5, 0.9, 1)),
  fish_sel = rbind(c(0.2, 0.6, 0.9, 1), # selectivity [region, age]
                   c(0.5, 1, 0.8, 0.4)),
  ret_sel = rbind(c(0.5, 0.8, 0.9, 1), # retention [region, age]
                  c(0.6, 0.7, 0.95, 0.85)),
  dmr = c(0.3, 0.5), # discard mortality by region
  Fmort = c(0.25, 0.12), # fishing mortality by region
  # movement fractions [from, to], one matrix per age
  move = list(matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2, byrow = TRUE),
              matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE),
              matrix(c(0.6, 0.4, 0.35, 0.65), nrow = 2, byrow = TRUE),
              matrix(c(0.85, 0.15, 0.5, 0.5), nrow = 2, byrow = TRUE)),
  # movement rates [from, to] for continuous movement, rows summing to zero
  rate = list(matrix(c(-0.1, 0.1, 0.2, -0.2), nrow = 2, byrow = TRUE),
              matrix(c(-0.5, 0.5, 0.3, -0.3), nrow = 2, byrow = TRUE),
              matrix(c(-0.8, 0.8, 0.6, -0.6), nrow = 2, byrow = TRUE),
              matrix(c(-0.2, 0.2, 0.9, -0.9), nrow = 2, byrow = TRUE)),
  rec_region_prop = c(0.7, 0.3), # share of recruits by region
  sexratio_f = c(0.5, 0.5), # female share of recruits by region
  t_spawn = 0.3
)

sbpr_hand_model <- function(rec_dd, move_timing, do_recruits_move) {

  n_regions <- sbpr_hand$n_regions
  n_ages <- sbpr_hand$n_ages
  Movement <- Mrate <- array(0, dim = c(1, n_regions, n_regions, 1, n_ages))
  for(a in 1:n_ages) {
    Movement[1,,,1,a] <- sbpr_hand$move[[a]]
    Mrate[1,,,1,a] <- sbpr_hand$rate[[a]]
  } # end a loop

  SPoRC:::Get_SBPR_Table(
    rec_dd = rec_dd,
    rec_region_prop = array(sbpr_hand$rec_region_prop, dim = c(1, n_regions)),
    rec_seas_prop = array(1, dim = c(1, 1)),
    n_pop = 1,
    n_regions = n_regions,
    n_ages = n_ages,
    n_fish_fleets = 1,
    WAA = array(sbpr_hand$WAA, dim = c(1, n_regions, 1, n_ages)),
    MatAA = array(sbpr_hand$MatAA, dim = c(1, n_regions, 1, n_ages)),
    natmort = array(sbpr_hand$natmort, dim = c(1, n_regions, 1, n_ages)),
    Movement = Movement,
    sgl_seas_spawning_movement = array(diag(n_regions), dim = c(1, n_regions, n_regions, n_ages)), # no natal homing
    do_recruits_move = do_recruits_move,
    t_spawn = sbpr_hand$t_spawn,
    init_F = array(sbpr_hand$Fmort, dim = c(n_regions, 1, 1)),
    dmr = array(sbpr_hand$dmr, dim = c(n_regions, 1, 1)),
    fish_sel = array(sbpr_hand$fish_sel, dim = c(1, n_regions, 1, n_ages, 1)),
    ret_sel = array(sbpr_hand$ret_sel, dim = c(1, n_regions, 1, n_ages, 1)),
    n_seas = 1,
    spawn_seas = 1,
    seasdur = 1,
    sexratio_f = array(sbpr_hand$sexratio_f, dim = c(1, n_regions)),
    Mrate = Mrate,
    move_timing = move_timing,
    expm_nsub = 0
  )
}

# exp(A) by its power series (A is 2 by 2 and small)
sbpr_hand_expm <- function(A) {
  term <- diag(nrow(A)) # A^k / k!
  expA <- matrix(0, nrow(A), ncol(A))
  for(k in 0:80) {
    expA <- expA + term
    term <- term %*% A / (k + 1)
  } # end k loop
  expA
}

# spawning biomass per recruit by region, from one year class followed until almost none are left
sbpr_hand_calc <- function(move_timing, do_recruits_move, fished) {

  n_regions <- sbpr_hand$n_regions
  n_ages <- sbpr_hand$n_ages
  numbers <- sbpr_hand$sexratio_f * sbpr_hand$rec_region_prop # female recruits by region, for one recruit in total
  sbpr <- rep(0, n_regions)

  for(year in 1:800) {

    a <- min(year, n_ages) # the plus group keeps its own values every year after

    # total mortality by region, with dead discards
    Z <- sbpr_hand$natmort[,a]
    if(fished) Z <- Z + sbpr_hand$Fmort * (sbpr_hand$fish_sel[,a] * sbpr_hand$ret_sel[,a] +
                                             sbpr_hand$fish_sel[,a] * (1 - sbpr_hand$ret_sel[,a]) * sbpr_hand$dmr)

    # recruits stay where they are unless they are set to move
    move_a <- diag(n_regions)
    rate_a <- matrix(0, n_regions, n_regions)
    if(do_recruits_move == 1 || a > 1) {
      move_a <- sbpr_hand$move[[a]]
      rate_a <- sbpr_hand$rate[[a]]
    }

    # numbers at spawning and at the end of the year, under each timing
    if(move_timing == 0) { # move, then die
      at_spawning <- as.vector(numbers %*% move_a) * exp(-sbpr_hand$t_spawn * Z)
      year_end <- as.vector(numbers %*% move_a) * exp(-Z)
    }
    if(move_timing == 1) { # die, then move
      at_spawning <- numbers * exp(-sbpr_hand$t_spawn * Z)
      year_end <- as.vector((numbers * exp(-Z)) %*% move_a)
    }
    if(move_timing == 2) { # both together
      at_spawning <- as.vector(numbers %*% sbpr_hand_expm((rate_a - diag(Z)) * sbpr_hand$t_spawn))
      year_end <- as.vector(numbers %*% sbpr_hand_expm(rate_a - diag(Z)))
    }

    sbpr <- sbpr + at_spawning * sbpr_hand$WAA[,a] * sbpr_hand$MatAA[,a]
    numbers <- year_end

  } # end year loop

  sbpr
}

# Tests ----------------------------------------------------------------------

test_that("spawning biomass per recruit matches a cohort followed to extinction", {

  for(move_timing in 0:2) {
    for(do_recruits_move in 0:1) {

      label <- sprintf("move_timing %d, do_recruits_move %d", move_timing, do_recruits_move)
      unfished <- sbpr_hand_calc(move_timing, do_recruits_move, fished = FALSE)
      fished <- sbpr_hand_calc(move_timing, do_recruits_move, fished = TRUE)

      # local density dependence: by region
      local <- sbpr_hand_model(rec_dd = 0, move_timing, do_recruits_move)
      expect_equal(as.numeric(local$phi0), unfished, tolerance = 1e-10, label = paste("unfished by region,", label))
      expect_equal(as.numeric(local$phiF), fished, tolerance = 1e-10, label = paste("fished by region,", label))

      # global density dependence: summed over regions
      global <- sbpr_hand_model(rec_dd = 1, move_timing, do_recruits_move)
      expect_equal(as.numeric(global$phi0), sum(unfished), tolerance = 1e-10, label = paste("unfished total,", label))
      expect_equal(as.numeric(global$phiF), sum(fished), tolerance = 1e-10, label = paste("fished total,", label))

    } # end do_recruits_move loop
  } # end move_timing loop

})

test_that("the cohort calculation separates the timings, recruit movement and fishing", {

  # if any of these did nothing in this setup, the test above could not tell the cases apart
  base <- sbpr_hand_calc(move_timing = 0, do_recruits_move = 0, fished = FALSE)

  expect_gt(max(abs(sbpr_hand_calc(1, 0, FALSE) / base - 1)), 0.01)
  expect_gt(max(abs(sbpr_hand_calc(2, 0, FALSE) / base - 1)), 0.01)
  expect_gt(max(abs(sbpr_hand_calc(0, 1, FALSE) / base - 1)), 1e-4)
  expect_gt(max(abs(sbpr_hand_calc(0, 0, TRUE) / base - 1)), 0.3)

})
