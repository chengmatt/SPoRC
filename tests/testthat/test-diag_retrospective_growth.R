# A retrospective peel has to shorten the growth deviations along with everything
# else. Their penalties re-dimension the deviation and map arrays to the model's
# year count, so an untruncated array is read in the wrong order rather than
# ignored past the peel, and the penalty comes out wrong without a word.

library(SPoRC)
library(testthat)
data("sgl_rg_ebs_pcod_data")

test_that("a peel truncates the growth deviations, their maps and the cohort start year", {

  dat <- sgl_rg_ebs_pcod_data
  n_yrs <- length(dat$years)
  peel <- 3

  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)
  cut <- truncate_yr(peel, input_list$data, input_list$par, input_list$map)

  # the deviations and both map mirrors lose the peeled years
  expect_equal(dim(cut$retro_parameters$ln_growth_devs)[3], n_yrs - peel)
  expect_equal(dim(cut$retro_data$map_ln_growth_devs)[3], n_yrs - peel)
  expect_equal(dim(cut$retro_parameters$ln_growth_semipar_devs)[3], n_yrs - peel)
  expect_equal(dim(cut$retro_data$map_ln_growth_semipar_devs)[3], n_yrs - peel)

  # the mapping factor stays the same length as the parameter it maps
  expect_equal(length(cut$retro_mapping$ln_growth_devs), prod(dim(cut$retro_parameters$ln_growth_devs)))
  expect_equal(length(cut$retro_mapping$ln_growth_semipar_devs), prod(dim(cut$retro_parameters$ln_growth_semipar_devs)))

  # and the years kept are the first ones, not a reshuffle
  expect_equal(cut$retro_parameters$ln_growth_devs[1, 1, , , 1],
               input_list$par$ln_growth_devs[1, 1, 1:(n_yrs - peel), , 1])

  # the cohort propagation cannot start after the end of a peeled series
  expect_lte(cut$retro_data$growth_cohort_styr, n_yrs - peel)
})


test_that("the peeled model penalizes only the deviations it kept", {

  dat <- sgl_rg_ebs_pcod_data
  n_yrs <- length(dat$years)
  peel <- 3

  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)
  full <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)$rep

  cut <- truncate_yr(peel, input_list$data, input_list$par, input_list$map)
  cut_fit <- fit_model(cut$retro_data, cut$retro_parameters, cut$retro_mapping, do_optim = FALSE, silent = TRUE)$rep

  # both growth parameters carry independent deviations with a fixed sigma, so
  # the penalty is a sum of normal densities and the peel drops the last three
  # years of each series
  devs <- input_list$par$ln_growth_devs
  sigma <- exp(input_list$par$growth_pe_pars[1, 1, , 1, 1])
  dropped <- 0
  for(k in which(input_list$data$growth_tv_model > 0)) {
    for(y in (n_yrs - peel + 1):n_yrs) {
      if(!is.na(input_list$data$map_ln_growth_devs[1, 1, y, k, 1])) {
        dropped <- dropped - stats::dnorm(devs[1, 1, y, k, 1], 0, sigma[k], log = TRUE)
      }
    } # end y loop
  } # end k loop

  # the peel takes real deviations out, so the two penalties differ
  expect_false(isTRUE(all.equal(cut_fit$growth_tv_nLL, full$growth_tv_nLL)))
  expect_equal(cut_fit$growth_tv_nLL, full$growth_tv_nLL - dropped, tolerance = 1e-10)
})


test_that("a peel truncates the conditional age-at-length arrays", {

  data("mlt_rg_goa_rex_data")
  dat <- mlt_rg_goa_rex_data
  n_yrs <- length(dat$years)
  peel <- 4

  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  cut <- truncate_yr(peel, input_list$data, input_list$par, input_list$map)

  expect_equal(dim(cut$retro_data$ObsSrv_caal)[2], n_yrs - peel)
  expect_equal(dim(cut$retro_data$UseSrv_caal)[2], n_yrs - peel)
  expect_equal(dim(cut$retro_data$ISS_Srv_caal)[2], n_yrs - peel)
  expect_equal(nrow(cut$retro_data$Srv_caal_Type), n_yrs - peel)
  expect_equal(dim(cut$retro_data$Wt_Srv_caal)[2], n_yrs - peel)

  # the years kept are the first ones
  expect_equal(cut$retro_data$ObsSrv_caal[, , 1, , , , 1],
               input_list$data$ObsSrv_caal[, 1:(n_yrs - peel), 1, , , , 1])
})
