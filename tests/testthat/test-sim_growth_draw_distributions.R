# Correlated growth deviation surfaces drawn by the operating model, measured over 4000 replicates (10 pct of the
# variance): the 3D GMRF draws have the covariance its precision matrix implies, with the same variance in every
# cell under the marginal form, and the separable AR1 draws have the marginal sd and the two correlations asked for.
# Years drawn after fitted years are centered on what the fitted years imply, under a random walk and both forms.

library(SPoRC)
library(testthat)

# Setup ----------------------------------------------------------------------

growth_draw_n_yrs <- 3
growth_draw_n_bins <- 3
growth_draw_n_sims <- 4000

# every cell of a year by bin surface its own free level, for one population, region and sex
growth_draw_map <- array(seq_len(growth_draw_n_yrs * growth_draw_n_bins), dim = c(1, 1, growth_draw_n_yrs, growth_draw_n_bins, 1))

# replicate surfaces, one column per replicate, cells in year by bin order with the bin moving fastest
growth_draw_surfaces <- function(PE_model, pe_pars, seed) {
  set.seed(seed)
  draws <- matrix(0, growth_draw_n_yrs * growth_draw_n_bins, growth_draw_n_sims)
  for(sim in 1:growth_draw_n_sims) {
    surface <- SPoRC:::draw_growth_pe_surface(PE_model = PE_model, map = growth_draw_map, pe_pars = pe_pars, bins = 1:growth_draw_n_bins)
    draws[,sim] <- as.vector(t(surface[1,1,,,1])) # [year, bin] read with the bin moving fastest
  } # end sim loop
  draws
}

# Tests ----------------------------------------------------------------------

test_that("3D GMRF growth surfaces have the covariance of the precision matrix the penalty uses", {

  pcorr_bin <- 0.3 # partial correlation between neighboring bins in a year
  pcorr_year <- 0.4 # between neighboring years at a bin
  pcorr_cohort <- 0.15 # along the cohort diagonal
  variance <- 0.25
  pe_pars <- array(c(pcorr_bin, pcorr_year, pcorr_cohort, log(variance)), dim = c(1, 1, 4, 1))

  # marginal form (3) and conditional form (4)
  for(PE_model in 3:4) {

    draws <- growth_draw_surfaces(PE_model = PE_model, pe_pars = pe_pars, seed = 20 + PE_model)
    var_type <- PE_model - 3 # the marginal form reads Var_Type 0, the conditional form 1
    precision <- as.matrix(SPoRC:::Get_3d_precision(n_ages = growth_draw_n_bins, n_yrs = growth_draw_n_yrs, pcorr_age = pcorr_bin,
                                                    pcorr_year = pcorr_year, pcorr_cohort = pcorr_cohort, ln_var_value = log(variance), Var_Type = var_type))
    expected_cov <- solve(precision)

    expect_lt(max(abs(stats::cov(t(draws)) - expected_cov)) / variance, 0.10)
    expect_lt(max(abs(rowMeans(draws))), 0.04) # centered on zero

    # the marginal form gives every cell the variance asked for
    if(PE_model == 3) expect_equal(as.numeric(apply(draws, 1, stats::var)), rep(variance, nrow(draws)), tolerance = 0.10)

  } # end PE_model loop

})

test_that("separable AR1 growth surfaces have the marginal sd and the year and bin correlations asked for", {

  rho_bin <- 0.5
  rho_year <- 0.7
  conditional_sd <- 0.2 # sd of a cell given its neighbors
  pe_pars <- array(c(atanh(rho_bin), atanh(rho_year), 0, log(conditional_sd)), dim = c(1, 1, 4, 1)) # correlations stored before their tanh transform

  draws <- growth_draw_surfaces(PE_model = 5, pe_pars = pe_pars, seed = 25)
  marginal_sd <- conditional_sd / sqrt(1 - rho_year^2) / sqrt(1 - rho_bin^2)

  # cell k is year (k - 1) %/% n_bins + 1 and bin (k - 1) %% n_bins + 1
  cell_sd <- apply(draws, 1, stats::sd)
  expect_equal(as.numeric(cell_sd), rep(marginal_sd, nrow(draws)), tolerance = 0.05)

  cell_cor <- stats::cor(t(draws))
  expect_equal(cell_cor[1,2], rho_bin, tolerance = 0.06) # year 1, bins 1 and 2
  expect_equal(cell_cor[1,4], rho_year, tolerance = 0.06) # bin 1, years 1 and 2
  expect_equal(cell_cor[1,5], rho_bin * rho_year, tolerance = 0.1) # one bin and one year apart
  expect_equal(cell_cor[1,7], rho_year^2, tolerance = 0.1) # bin 1, years 1 and 3

})

# Draws Given Fitted Years ---------------------------------------------------

test_that("years drawn after the fitted years are centered on what the fitted years imply", {

  n_yrs <- 3
  n_bins <- 2
  n_cond <- 2 # the first two years keep the fit's deviations and the third is drawn
  n_sims <- 3000
  map <- array(seq_len(n_yrs * n_bins), dim = c(1, 1, n_yrs, n_bins, 1))
  fit_devs <- array(0, dim = c(1, 1, n_yrs, n_bins, 1))
  fit_devs[1,1,1,,1] <- c(0.30, -0.10) # year 1, bins 1 and 2
  fit_devs[1,1,2,,1] <- c(0.50, -0.40) # year 2

  # year 3's draws [bin, replicate], with the fitted years returned unchanged
  drawn_year <- function(PE_model, pe_pars, seed) {
    set.seed(seed)
    draws <- matrix(0, n_bins, n_sims)
    for(sim in 1:n_sims) {
      surface <- SPoRC:::draw_growth_pe_surface(PE_model = PE_model, map = map, pe_pars = pe_pars, bins = 1:n_bins, fit_devs = fit_devs, n_cond = n_cond)
      if(sim == 1) expect_equal(surface[1,1,1:n_cond,,1], fit_devs[1,1,1:n_cond,,1])
      draws[,sim] <- surface[1,1,n_yrs,,1]
    } # end sim loop
    draws
  }

  # random walk: each bin steps on from the fit's last value at its own sd
  sigma <- c(0.2, 0.1)
  draws <- drawn_year(PE_model = 2, pe_pars = array(log(sigma), dim = c(1, 1, n_bins, 1)), seed = 41)
  for(bin in 1:n_bins) {
    expect_lt(abs(mean(draws[bin,]) - fit_devs[1,1,n_cond,bin,1]), 4 * sigma[bin] / sqrt(n_sims))
    expect_equal(stats::sd(draws[bin,]), sigma[bin], tolerance = 0.06)
  } # end bin loop

  # separable AR1: rho_year times last year's deviation, with the spread left after knowing last year
  rho_bin <- 0.5
  rho_year <- 0.7
  conditional_sd <- 0.2
  draws <- drawn_year(PE_model = 5, pe_pars = array(c(atanh(rho_bin), atanh(rho_year), 0, log(conditional_sd)), dim = c(1, 1, 4, 1)), seed = 42)
  marginal_sd <- conditional_sd / sqrt(1 - rho_year^2) / sqrt(1 - rho_bin^2)
  sd_given_last_year <- marginal_sd * sqrt(1 - rho_year^2)
  for(bin in 1:n_bins) {
    expect_lt(abs(mean(draws[bin,]) - rho_year * fit_devs[1,1,n_cond,bin,1]), 4 * sd_given_last_year / sqrt(n_sims))
    expect_equal(stats::sd(draws[bin,]), sd_given_last_year, tolerance = 0.06)
  } # end bin loop
  expect_equal(stats::cor(draws[1,], draws[2,]), rho_bin, tolerance = 0.08)

  # 3D GMRF: the mean and covariance of year 3 given years 1 and 2, worked out from the covariance matrix
  pcorr_bin <- 0.3
  pcorr_year <- 0.4
  pcorr_cohort <- 0.15
  variance <- 0.25
  draws <- drawn_year(PE_model = 4, pe_pars = array(c(pcorr_bin, pcorr_year, pcorr_cohort, log(variance)), dim = c(1, 1, 4, 1)), seed = 43)
  covariance <- solve(as.matrix(SPoRC:::Get_3d_precision(n_ages = n_bins, n_yrs = n_yrs, pcorr_age = pcorr_bin, pcorr_year = pcorr_year,
                                                         pcorr_cohort = pcorr_cohort, ln_var_value = log(variance), Var_Type = 1)))
  fitted_cells <- 1:(n_cond * n_bins) # cells run year by year with the bin moving fastest
  drawn_cells <- (n_cond * n_bins + 1):(n_yrs * n_bins)
  fitted_values <- as.vector(t(fit_devs[1,1,1:n_cond,,1]))
  weights <- covariance[drawn_cells, fitted_cells] %*% solve(covariance[fitted_cells, fitted_cells])
  expected_mean <- as.vector(weights %*% fitted_values)
  expected_cov <- covariance[drawn_cells, drawn_cells] - weights %*% covariance[fitted_cells, drawn_cells]
  for(bin in 1:n_bins) expect_lt(abs(mean(draws[bin,]) - expected_mean[bin]), 4 * sqrt(expected_cov[bin,bin] / n_sims))
  expect_lt(max(abs(stats::cov(t(draws)) - expected_cov)) / variance, 0.08)

  # the fitted years move the mean well away from zero here
  expect_gt(max(abs(expected_mean)), 0.1)

})
