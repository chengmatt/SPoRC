# Do_Population_Projection at input fishing mortality, against a hand calculation written from the equations (1e-8):
# numbers at age, retained and discard fishing mortality, catch, discards, spawning biomass and dynamic B0. Two fleets that
# discard, mortality that differs by age, and spawning half way through the year.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

proj_hand <- list(
  n_proj_yrs = 5,
  n_ages = 4, # age 4 is the plus group
  n_fleets = 2,
  natmort = c(0.35, 0.25, 0.18, 0.12), # per year, by age
  WAA = c(0.5, 1.1, 1.8, 2.6), # weight at age (kg)
  MatAA = c(0, 0.3, 0.8, 1), # maturity at age
  fish_sel = cbind(c(0.2, 0.6, 1, 1), c(1, 0.8, 0.5, 0.3)), # selectivity [age, fleet]
  ret_sel = cbind(c(0.3, 0.7, 0.9, 1), c(0.5, 0.8, 1, 1)), # retention [age, fleet]
  dmr = c(0.3, 0.6), # discard mortality by fleet
  terminal_F = c(0.12, 0.06), # year 1 fishing mortality by fleet
  f_input = c(0.25, 0.10, 0.30, 0.05, 0.40), # f_input[y] is the total F in year y + 1
  recruitment = c(900, 1100, 1300, 700), # past recruitment, averaged for the projection
  rec_devs = c(1, 1.4, 0.6, 1.1, 0.8), # multiplier on mean recruitment by year
  terminal_NAA = c(1000, 700, 450, 600), # fished numbers at age in year 1
  terminal_NAA0 = c(1000, 900, 750, 1500), # unfished numbers at age in year 1
  t_spawn = 0.5
)
proj_hand$WAA_fish <- cbind(proj_hand$WAA * 1.1, proj_hand$WAA * 0.9) # fleet weights differ from the population's

# the projection's arrays, with every year at the same biology
proj_hand_run <- function() {

  n_proj_yrs <- proj_hand$n_proj_yrs
  n_ages <- proj_hand$n_ages
  n_fleets <- proj_hand$n_fleets
  biol_dim <- c(1, 1, n_proj_yrs, 1, n_ages, 1)

  natmort <- WAA <- MatAA <- array(0, dim = biol_dim)
  WAA_fish <- fish_sel <- ret_sel <- array(0, dim = c(biol_dim, n_fleets))
  for(a in 1:n_ages) {
    natmort[1,1,,1,a,1] <- proj_hand$natmort[a]
    WAA[1,1,,1,a,1] <- proj_hand$WAA[a]
    MatAA[1,1,,1,a,1] <- proj_hand$MatAA[a]
    for(f in 1:n_fleets) {
      WAA_fish[1,1,,1,a,1,f] <- proj_hand$WAA_fish[a,f]
      fish_sel[1,1,,1,a,1,f] <- proj_hand$fish_sel[a,f]
      ret_sel[1,1,,1,a,1,f] <- proj_hand$ret_sel[a,f]
    } # end f loop
  } # end a loop

  Do_Population_Projection(
    n_proj_yrs = n_proj_yrs,
    n_pop = 1,
    n_regions = 1,
    n_ages = n_ages,
    n_sexes = 1,
    n_fish_fleets = n_fleets,
    sexratio = array(1, dim = c(1, 1, n_proj_yrs, 1)),
    do_recruits_move = 0,
    recruitment = array(proj_hand$recruitment, dim = c(1, 1, length(proj_hand$recruitment))),
    terminal_NAA = array(proj_hand$terminal_NAA, dim = c(1, 1, 1, n_ages, 1)),
    terminal_NAA0 = array(proj_hand$terminal_NAA0, dim = c(1, 1, 1, n_ages, 1)),
    terminal_F = array(proj_hand$terminal_F, dim = c(1, 1, n_fleets)),
    dmr = array(proj_hand$dmr, dim = c(1, 1, n_fleets)),
    natmort = natmort,
    WAA = WAA,
    WAA_fish = WAA_fish,
    MatAA = MatAA,
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    Movement = array(1, dim = c(1, 1, 1, n_proj_yrs, 1, n_ages, 1)),
    f_ref_pt = array(proj_hand$f_input, dim = c(1, n_proj_yrs)),
    b_ref_pt = array(1, dim = c(1, 1, n_proj_yrs)),
    recruitment_opt = "mean_rec",
    fmort_opt = "Input",
    rec_devs = array(proj_hand$rec_devs, dim = c(1, 1, n_proj_yrs)),
    t_spawn = proj_hand$t_spawn
  )
}

# the same projection one year, age and fleet at a time
proj_hand_calc <- function() {

  n_proj_yrs <- proj_hand$n_proj_yrs
  n_ages <- proj_hand$n_ages
  n_fleets <- proj_hand$n_fleets

  NAA <- NAA0 <- ZAA <- matrix(0, n_proj_yrs, n_ages)
  ret_FAA <- disc_FAA <- CAA <- DAA <- array(0, dim = c(n_proj_yrs, n_ages, n_fleets))
  Catch <- matrix(0, n_proj_yrs, n_fleets)
  SSB <- Dynamic_SSB0 <- rep(0, n_proj_yrs)
  NAA[1,] <- proj_hand$terminal_NAA
  NAA0[1,] <- proj_hand$terminal_NAA0
  fleet_share <- proj_hand$terminal_F / sum(proj_hand$terminal_F) # the total F is split as the terminal year's was

  for(y in 1:n_proj_yrs) {

    # year 1 fishes at the terminal F, later years at the input for the year before
    F_fleet <- proj_hand$terminal_F
    if(y > 1) F_fleet <- proj_hand$f_input[y - 1] * fleet_share

    # recruits enter both the fished and the unfished population
    if(y > 1) NAA[y,1] <- NAA0[y,1] <- mean(proj_hand$recruitment) * proj_hand$rec_devs[y]

    for(a in 1:n_ages) {

      ZAA[y,a] <- proj_hand$natmort[a]
      for(f in 1:n_fleets) {
        ret_FAA[y,a,f] <- F_fleet[f] * proj_hand$fish_sel[a,f] * proj_hand$ret_sel[a,f] # retained
        disc_FAA[y,a,f] <- F_fleet[f] * proj_hand$fish_sel[a,f] * (1 - proj_hand$ret_sel[a,f]) * proj_hand$dmr[f] # discarded and dead
        ZAA[y,a] <- ZAA[y,a] + ret_FAA[y,a,f] + disc_FAA[y,a,f]
      } # end f loop

      # Baranov catch and dead discards in numbers
      for(f in 1:n_fleets) {
        CAA[y,a,f] <- ret_FAA[y,a,f] / ZAA[y,a] * NAA[y,a] * (1 - exp(-ZAA[y,a]))
        DAA[y,a,f] <- disc_FAA[y,a,f] / ZAA[y,a] * NAA[y,a] * (1 - exp(-ZAA[y,a]))
        Catch[y,f] <- Catch[y,f] + CAA[y,a,f] * proj_hand$WAA_fish[a,f]
      } # end f loop

      # one sex holds both, so half the mature biomass at spawning is female
      SSB[y] <- SSB[y] + 0.5 * NAA[y,a] * exp(-proj_hand$t_spawn * ZAA[y,a]) * proj_hand$WAA[a] * proj_hand$MatAA[a]
      Dynamic_SSB0[y] <- Dynamic_SSB0[y] + 0.5 * NAA0[y,a] * exp(-proj_hand$t_spawn * proj_hand$natmort[a]) * proj_hand$WAA[a] * proj_hand$MatAA[a]

    } # end a loop

    # ageing, with the plus group taking the oldest true age and its own survivors
    if(y < n_proj_yrs) {
      for(a in 1:(n_ages - 1)) {
        NAA[y + 1,a + 1] <- NAA[y,a] * exp(-ZAA[y,a])
        NAA0[y + 1,a + 1] <- NAA0[y,a] * exp(-proj_hand$natmort[a])
      } # end a loop
      NAA[y + 1,n_ages] <- NAA[y + 1,n_ages] + NAA[y,n_ages] * exp(-ZAA[y,n_ages])
      NAA0[y + 1,n_ages] <- NAA0[y + 1,n_ages] + NAA0[y,n_ages] * exp(-proj_hand$natmort[n_ages])
    }

  } # end y loop

  list(NAA = NAA, NAA0 = NAA0, ZAA = ZAA, ret_FAA = ret_FAA, disc_FAA = disc_FAA, CAA = CAA, DAA = DAA,
       Catch = Catch, SSB = SSB, Dynamic_SSB0 = Dynamic_SSB0)
}

# Tests ----------------------------------------------------------------------

test_that("projected mortality and numbers at age match the hand calculation", {

  model <- proj_hand_run()
  hand <- proj_hand_calc()
  yrs <- 1:proj_hand$n_proj_yrs

  expect_equal(model$proj_ret_FAA[1,1,yrs,1,,1,], hand$ret_FAA, tolerance = 1e-8)
  expect_equal(model$proj_disc_FAA[1,1,yrs,1,,1,], hand$disc_FAA, tolerance = 1e-8)
  expect_equal(model$proj_ZAA[1,1,yrs,1,,1], hand$ZAA, tolerance = 1e-8)
  expect_equal(model$proj_NAA[1,1,yrs,1,,1], hand$NAA, tolerance = 1e-8)
  expect_equal(model$proj_NAA0[1,1,yrs,1,,1], hand$NAA0, tolerance = 1e-8)

})

test_that("projected catch, spawning biomass and dynamic B0 match the hand calculation", {

  model <- proj_hand_run()
  hand <- proj_hand_calc()

  expect_equal(model$proj_CAA[1,1,,1,,1,], hand$CAA, tolerance = 1e-8)
  expect_equal(model$proj_Catch[1,1,,1,], hand$Catch, tolerance = 1e-8)
  expect_equal(model$proj_SSB[1,1,], hand$SSB, tolerance = 1e-8)
  expect_equal(model$proj_Dynamic_SSB0[1,1,], hand$Dynamic_SSB0, tolerance = 1e-8)

})

test_that("projected dead discards at age match the hand calculation", {

  model <- proj_hand_run()
  hand <- proj_hand_calc()

  expect_equal(model$proj_DAA[1,1,,1,,1,], hand$DAA, tolerance = 1e-8)
  expect_gt(max(hand$DAA), 10) # the setup does discard

})
