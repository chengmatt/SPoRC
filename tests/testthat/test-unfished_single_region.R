# Unfished numbers at age and dynamic B0 in a one-region model, against a hand calculation from the reported
# natural mortality (1e-10). Two seasons, spawning part way into the second, and natural mortality that
# differs by age so the oldest true age and the plus group cannot be swapped.

# Setup ----------------------------------------------------------------------

# one region and one fleet, natural mortality falling with age
unfished_one_region_input <- function(n_sexes = 1, fished = TRUE) {

  n_ages <- 6
  n_seas <- 2
  natmort_at_age <- c(0.45, 0.35, 0.28, 0.22, 0.17, 0.12) # per year, a different value at every age

  il <- collapse_input(nr = 1, nx = n_sexes, nf = 1, ns = n_seas)
  for(a in 1:n_ages) il$data$Fixed_natmort[,,,,a,] <- natmort_at_age[a]
  il$data$t_spawn <- 0.4 # spawning 40 pct of the way into the spawning season
  if(!fished) il$par$ln_F_mean[] <- -30 # fishing mortality of effectively zero

  il
}

unfished_one_region_rep <- function(n_sexes = 1, fished = TRUE) {
  il <- unfished_one_region_input(n_sexes = n_sexes, fished = fished)
  fit <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  list(il = il, rep = fit$rep)
}

# Tests ----------------------------------------------------------------------

test_that("unfished numbers at age die at natural mortality alone within a year", {

  out <- unfished_one_region_rep()
  rep <- out$rep
  seasdur <- out$il$data$seasdur
  n_yrs <- length(out$il$data$years)
  n_ages <- length(out$il$data$ages)

  # recruits enter at age 1, so ages 2 and older in season 2 are the survivors of season 1
  for(y in 1:n_yrs) {
    for(a in 2:n_ages) {
      survivors <- rep$NAA0[1,1,y,1,a,1] * exp(-rep$natmort[1,1,y,1,a,1] * seasdur[1])
      expect_equal(rep$NAA0[1,1,y,2,a,1], survivors, tolerance = 1e-10,
                   label = sprintf("unfished numbers, year %d season 2 age %d", y, a))
    } # end a loop
  } # end y loop

})

test_that("unfished numbers at age advance a year and fill the plus group", {

  out <- unfished_one_region_rep()
  rep <- out$rep
  seasdur <- out$il$data$seasdur
  n_yrs <- length(out$il$data$years)
  n_ages <- length(out$il$data$ages)
  n_seas <- out$il$data$n_seas

  for(y in 1:(n_yrs - 1)) {

    # survivors of the last season, before they age
    survivors <- rep(0, n_ages)
    for(a in 1:n_ages) survivors[a] <- rep$NAA0[1,1,y,n_seas,a,1] * exp(-rep$natmort[1,1,y,n_seas,a,1] * seasdur[n_seas])

    # ages 2 to the oldest true age are last year's survivors one age younger
    for(a in 2:(n_ages - 1)) {
      expect_equal(rep$NAA0[1,1,y + 1,1,a,1], survivors[a - 1], tolerance = 1e-10,
                   label = sprintf("unfished numbers, year %d age %d", y + 1, a))
    } # end a loop

    # the plus group takes the oldest true age and its own survivors
    expect_equal(rep$NAA0[1,1,y + 1,1,n_ages,1], survivors[n_ages - 1] + survivors[n_ages], tolerance = 1e-10,
                 label = sprintf("unfished plus group, year %d", y + 1))

  } # end y loop

})

test_that("dynamic B0 weighs the unfished numbers at spawning, halved for one sex", {

  for(n_sexes in 1:2) {

    out <- unfished_one_region_rep(n_sexes = n_sexes)
    rep <- out$rep
    data <- out$il$data
    n_yrs <- length(data$years)
    n_ages <- length(data$ages)
    spawn_seas <- data$spawn_seas
    female_share <- if(n_sexes == 1) 0.5 else 1 # one sex holds both, so half are female

    for(y in 1:n_yrs) {
      dynamic_b0 <- 0
      for(a in 1:n_ages) {
        at_spawning <- rep$NAA0[1,1,y,spawn_seas,a,1] * exp(-rep$natmort[1,1,y,spawn_seas,a,1] * data$t_spawn * data$seasdur[spawn_seas])
        dynamic_b0 <- dynamic_b0 + at_spawning * data$WAA[1,1,y,spawn_seas,a,1] * data$MatAA[1,1,y,spawn_seas,a,1]
      } # end a loop
      expect_equal(rep$Dynamic_SSB0[1,1,y], dynamic_b0 * female_share, tolerance = 1e-10,
                   label = sprintf("dynamic B0, %d sex, year %d", n_sexes, y))
    } # end y loop

  } # end n_sexes loop

})

test_that("with no fishing the unfished and fished populations are the same", {

  rep <- unfished_one_region_rep(fished = FALSE)$rep

  expect_equal(rep$NAA0, rep$NAA, tolerance = 1e-10)
  expect_equal(rep$Dynamic_SSB0, rep$SSB, tolerance = 1e-10)

})

test_that("fishing lowers the fished population below the unfished one", {

  # without this the checks above would pass in a model that never fished
  rep <- unfished_one_region_rep(fished = TRUE)$rep
  n_yrs <- dim(rep$SSB)[3]

  expect_lt(rep$SSB[1,1,n_yrs], 0.95 * rep$Dynamic_SSB0[1,1,n_yrs])

})
