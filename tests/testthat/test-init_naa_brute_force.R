# Initial numbers at age under each init_age_strc against a population run for 200 years at constant recruitment
# (1e-8). Two regions, seasons, sexes and fleets, with discards, under every movement timing. Mortality, movement,
# selectivity and retention differ by region, season, age and sex.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

init_hand <- list(
  n_regions = 2,
  n_ages = 4, # age 4 is the plus group
  n_seas = 2,
  n_sexes = 2,
  n_fleets = 2,
  seasdur = c(0.4, 0.6), # share of the year in each season
  rec_seas_prop = matrix(c(0.7, 0.3), nrow = 1), # share of recruits entering in each season
  R0_r = matrix(c(800, 300), nrow = 1), # recruits by region
  sexratio = array(c(0.45, 0.5, 0.55, 0.5), dim = c(1, 2, 2)) # share of recruits by sex [pop, region, sex]
)

# every input array, with a different value in every cell
init_hand_arrays <- function(dmr_by_fleet) {

  n_regions <- init_hand$n_regions
  n_ages <- init_hand$n_ages
  n_seas <- init_hand$n_seas
  n_sexes <- init_hand$n_sexes
  n_fleets <- init_hand$n_fleets

  natmort <- array(0, dim = c(1, n_regions, n_seas, n_ages, n_sexes)) # per year
  fish_sel <- ret_sel <- array(0, dim = c(1, n_regions, n_seas, n_ages, n_sexes, n_fleets))
  init_F <- dmr <- array(0, dim = c(n_regions, n_seas, n_fleets))
  Movement <- Mrate <- array(0, dim = c(1, n_regions, n_regions, n_seas, n_ages, n_sexes))

  for(r in 1:n_regions) {
    for(seas in 1:n_seas) {

      for(f in 1:n_fleets) {
        init_F[r,seas,f] <- 0.1 * c(1, 0.6)[r] * c(1, 1.5)[seas] * c(1, 0.5)[f] # fishing mortality over the season
        dmr[r,seas,f] <- dmr_by_fleet[f] * c(1, 0.9)[r] # share of discards that die
      } # end f loop

      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {

          natmort[1,r,seas,a,s] <- 0.25 * c(1, 1.2)[r] * c(1, 0.8)[seas] * c(1.4, 1, 0.8, 0.6)[a] * c(1, 1.1)[s]
          for(f in 1:n_fleets) {
            fish_sel[1,r,seas,a,s,f] <- min(1, 0.3 * a * c(1, 0.9)[r] * c(1, 0.8)[f])
            ret_sel[1,r,seas,a,s,f] <- min(1, 0.4 + 0.15 * a + 0.05 * s) # the rest are discarded
          } # end f loop

          # movement fractions [from, to], rows summing to one
          stay <- c(0.9 - 0.1 * a - 0.03 * s - 0.05 * seas, 0.6 + 0.05 * a + 0.02 * seas)
          Movement[1,1,,seas,a,s] <- c(stay[1], 1 - stay[1])
          Movement[1,2,,seas,a,s] <- c(1 - stay[2], stay[2])

          # movement rates [from, to] per year, rows summing to zero
          leave <- c(0.2 * a + 0.1 * seas, 0.5 / a + 0.05 * s)
          Mrate[1,1,,seas,a,s] <- c(-leave[1], leave[1])
          Mrate[1,2,,seas,a,s] <- c(leave[2], -leave[2])

        } # end s loop
      } # end a loop
    } # end seas loop
  } # end r loop

  list(natmort = natmort, fish_sel = fish_sel, ret_sel = ret_sel, init_F = init_F, dmr = dmr, Movement = Movement, Mrate = Mrate)
}

# the same arrays with no movement at the ages given
init_hand_no_move <- function(arrays, ages) {
  for(seas in 1:init_hand$n_seas) {
    for(a in ages) {
      for(s in 1:init_hand$n_sexes) {
        arrays$Movement[1,,,seas,a,s] <- diag(init_hand$n_regions)
        arrays$Mrate[1,,,seas,a,s] <- 0
      } # end s loop
    } # end a loop
  } # end seas loop
  arrays
}

init_hand_model <- function(arrays, init_age_strc, move_timing, do_recruits_move, ln_InitDevs = NULL) {

  if(is.null(ln_InitDevs)) ln_InitDevs <- array(0, dim = c(1, init_hand$n_regions, init_hand$n_ages - 1, init_hand$n_sexes))

  SPoRC:::Get_Init_NAA(
    init_age_strc = init_age_strc,
    init_iter = 150,
    n_regions = init_hand$n_regions,
    n_pop = 1,
    n_sexes = init_hand$n_sexes,
    n_ages = init_hand$n_ages,
    n_seas = init_hand$n_seas,
    n_fish_fleets = init_hand$n_fleets,
    seasdur = init_hand$seasdur,
    rec_seas_prop = init_hand$rec_seas_prop,
    natmort = arrays$natmort,
    init_F = arrays$init_F,
    dmr = arrays$dmr,
    fish_sel = arrays$fish_sel,
    ret_sel = arrays$ret_sel,
    R0_r = init_hand$R0_r,
    sexratio = init_hand$sexratio,
    Movement = arrays$Movement,
    do_recruits_move = do_recruits_move,
    ln_InitDevs = ln_InitDevs,
    Mrate = arrays$Mrate,
    move_timing = move_timing,
    expm_nsub = 0
  )[1,,,]
}

# exp(A) by its power series (A is 2 by 2 and small)
init_hand_expm <- function(A) {
  term <- diag(nrow(A)) # A^k / k!
  expA <- matrix(0, nrow(A), ncol(A))
  for(k in 0:80) {
    expA <- expA + term
    term <- term %*% A / (k + 1)
  } # end k loop
  expA
}

# numbers at the start of the year [region, age, sex] after 200 years at constant recruitment
init_hand_calc <- function(arrays, move_timing, do_recruits_move) {

  n_regions <- init_hand$n_regions
  n_ages <- init_hand$n_ages
  n_seas <- init_hand$n_seas
  start_of_year <- array(0, dim = c(n_regions, n_ages, init_hand$n_sexes))

  for(s in 1:init_hand$n_sexes) {

    # where a region's fish are at the end of each season [from, to], one matrix per season and age
    season_step <- vector("list", n_seas * n_ages)
    for(seas in 1:n_seas) {
      for(a in 1:n_ages) {

        # total mortality over the season by region, each fleet's discards dying at its own rate
        Z <- arrays$natmort[1,,seas,a,s] * init_hand$seasdur[seas]
        for(f in 1:init_hand$n_fleets) {
          retained <- arrays$fish_sel[1,,seas,a,s,f] * arrays$ret_sel[1,,seas,a,s,f]
          dead_discards <- arrays$fish_sel[1,,seas,a,s,f] * (1 - arrays$ret_sel[1,,seas,a,s,f]) * arrays$dmr[,seas,f]
          Z <- Z + arrays$init_F[,seas,f] * (retained + dead_discards)
        } # end f loop

        # recruits stay where they are unless they are set to move
        move_a <- diag(n_regions)
        rate_a <- matrix(0, n_regions, n_regions)
        if(do_recruits_move == 1 || a > 1) {
          move_a <- arrays$Movement[1,,,seas,a,s]
          rate_a <- arrays$Mrate[1,,,seas,a,s]
        }

        if(move_timing == 0) step <- move_a %*% diag(exp(-Z)) # move, then die
        if(move_timing == 1) step <- diag(exp(-Z)) %*% move_a # die, then move
        if(move_timing == 2) step <- init_hand_expm(rate_a * init_hand$seasdur[seas] - diag(Z)) # together
        season_step[[(seas - 1) * n_ages + a]] <- step

      } # end a loop
    } # end seas loop

    numbers <- matrix(0, n_regions, n_ages)

    for(year in 1:200) {
      for(seas in 1:n_seas) {

        # this season's recruits join age 1
        numbers[,1] <- numbers[,1] + init_hand$R0_r[1,] * init_hand$sexratio[1,,s] * init_hand$rec_seas_prop[1,seas]
        if(year == 200 && seas == 1) start_of_year[,,s] <- numbers

        for(a in 1:n_ages) numbers[,a] <- as.vector(numbers[,a] %*% season_step[[(seas - 1) * n_ages + a]])

      } # end seas loop

      # ageing, with the plus group taking the oldest true age and its own survivors
      aged <- matrix(0, n_regions, n_ages)
      aged[,2:n_ages] <- numbers[,1:(n_ages - 1)]
      aged[,n_ages] <- aged[,n_ages] + numbers[,n_ages]
      numbers <- aged

    } # end year loop
  } # end s loop

  start_of_year
}

# every init_age_strc that solves an equilibrium, each against the long run it should reproduce
expect_init_naa_matches <- function(dmr_by_fleet) {

  arrays <- init_hand_arrays(dmr_by_fleet)
  arrays_no_move <- init_hand_no_move(arrays, ages = 1:init_hand$n_ages)
  arrays_top_stay <- init_hand_no_move(arrays, ages = (init_hand$n_ages - 1):init_hand$n_ages)

  for(move_timing in 0:2) {
    for(do_recruits_move in 0:1) {

      label <- sprintf("move_timing %d, do_recruits_move %d", move_timing, do_recruits_move)
      long_run <- init_hand_calc(arrays, move_timing, do_recruits_move)

      # the iterated annual cycle and the matrix series, with movement at every age
      expect_equal(init_hand_model(arrays, 0, move_timing, do_recruits_move), long_run, tolerance = 1e-8, label = paste("init_age_strc 0,", label))
      expect_equal(init_hand_model(arrays, 2, move_timing, do_recruits_move), long_run, tolerance = 1e-8, label = paste("init_age_strc 2,", label))

      # the scalar series, exact when no age moves
      expect_equal(init_hand_model(arrays_no_move, 1, move_timing, do_recruits_move),
                   init_hand_calc(arrays_no_move, move_timing, do_recruits_move), tolerance = 1e-8, label = paste("init_age_strc 1,", label))

      # movement below the plus group and a scalar plus group, exact when the two oldest ages stay
      expect_equal(init_hand_model(arrays_top_stay, 3, move_timing, do_recruits_move),
                   init_hand_calc(arrays_top_stay, move_timing, do_recruits_move), tolerance = 1e-8, label = paste("init_age_strc 3,", label))

    } # end do_recruits_move loop
  } # end move_timing loop
}

# Tests ----------------------------------------------------------------------

test_that("each equilibrium option matches the long run when fleets share one discard mortality", {
  expect_init_naa_matches(dmr_by_fleet = c(0.4, 0.4))
})

test_that("each equilibrium option matches the long run when discard mortality differs by fleet", {
  expect_init_naa_matches(dmr_by_fleet = c(0.3, 0.6))
})

test_that("the scalar series ignores movement it is given", {

  arrays <- init_hand_arrays(dmr_by_fleet = c(0.4, 0.4))
  arrays_no_move <- init_hand_no_move(arrays, ages = 1:init_hand$n_ages)

  expect_equal(init_hand_model(arrays, 1, 0, 0), init_hand_model(arrays_no_move, 1, 0, 0), tolerance = 1e-12)

})

test_that("free initial ages are the deviations split by sex, and deviations multiply an equilibrium", {

  arrays <- init_hand_arrays(dmr_by_fleet = c(0.4, 0.4))
  n_regions <- init_hand$n_regions
  n_ages <- init_hand$n_ages
  n_sexes <- init_hand$n_sexes

  ln_InitDevs <- array(0, dim = c(1, n_regions, n_ages - 1, n_sexes))
  for(r in 1:n_regions) for(a in 2:n_ages) for(s in 1:n_sexes) ln_InitDevs[1,r,a - 1,s] <- 0.2 * r - 0.15 * a + 0.1 * s

  free <- init_hand_model(arrays, 4, 0, 0, ln_InitDevs)
  equilibrium <- init_hand_model(arrays, 2, 0, 0)
  deviated <- init_hand_model(arrays, 2, 0, 0, ln_InitDevs)

  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {

      # age 1 is always the first season's recruits
      recruits <- init_hand$R0_r[1,r] * init_hand$sexratio[1,r,s] * init_hand$rec_seas_prop[1,1]
      expect_equal(free[r,1,s], recruits, tolerance = 1e-12)
      expect_equal(deviated[r,1,s], recruits, tolerance = 1e-12)

      for(a in 2:n_ages) {
        expect_equal(free[r,a,s], init_hand$sexratio[1,r,s] * exp(ln_InitDevs[1,r,a - 1,s]), tolerance = 1e-12)
        expect_equal(deviated[r,a,s], equilibrium[r,a,s] * exp(ln_InitDevs[1,r,a - 1,s]), tolerance = 1e-12)
      } # end a loop

    } # end s loop
  } # end r loop

})
