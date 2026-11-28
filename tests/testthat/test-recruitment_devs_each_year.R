# Every recruitment deviation multiplies its own year's recruitment by exp(deviation) and leaves earlier years
# alone (1e-10), the first and the last year included, under recruitment a year after spawning and in the same year.

# Setup ----------------------------------------------------------------------

# recruitment by year from the Beverton-Holt sweep model at a given set of deviations
rec_devs_recruitment <- function(rec_lag, ln_RecDevs) {

  n_yrs <- 10
  n_ages <- 6
  MatAA <- array(0, dim = c(1, 1, n_yrs, 1, n_ages, 1))
  for(a in 1:n_ages) MatAA[,,,,a,] <- 1 / (1 + exp(-(a - 3)))
  if(rec_lag == 0) MatAA[,,,,1,] <- 0 # recruits of the year cannot be mature

  il <- suppressWarnings(suppressMessages(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_yrs = n_yrs, n_ages = n_ages),
    rec = list(rec_model = "bh_rec", h_spec = "fix", rec_lag = rec_lag),
    biol = list(MatAA = MatAA)
  )))
  il$par$ln_RecDevs[1,1,] <- ln_RecDevs
  suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)))$rep$Rec[1,1,]
}

# Tests ----------------------------------------------------------------------

test_that("a recruitment deviation scales its own year's recruitment and no earlier year's", {

  n_yrs <- 10
  dev <- 0.3

  for(rec_lag in c(1, 0)) {

    base <- rec_devs_recruitment(rec_lag, rep(0, n_yrs))

    for(dev_year in c(1, 5, n_yrs - 1, n_yrs)) {

      ln_RecDevs <- rep(0, n_yrs)
      ln_RecDevs[dev_year] <- dev
      change <- rec_devs_recruitment(rec_lag, ln_RecDevs) / base
      label <- sprintf("rec_lag %d, deviation in year %d", rec_lag, dev_year)

      expect_equal(change[dev_year], exp(dev), tolerance = 1e-10, label = paste("that year's recruitment,", label))
      if(dev_year > 1) expect_equal(change[1:(dev_year - 1)], rep(1, dev_year - 1), tolerance = 1e-10, label = paste("earlier years,", label))

    } # end dev_year loop
  } # end rec_lag loop

})
