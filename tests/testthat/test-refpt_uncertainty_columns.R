# The interval and CV columns of Get_Reference_Point_Uncertainty against their definitions (1e-10): lognormal
# bounds and CV from the log standard error, and the drawn bounds as quantiles of the draws.

library(SPoRC)
library(testthat)
data("dusky_rtmb_model")

# the fit's standard errors, worked out once and passed to both calls
dusky_sd_rep <- RTMB::sdreport(dusky_rtmb_model)

test_that("the delta method bounds and CV follow from the log standard error", {

  level <- 0.90
  rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model,
                                        SPR_x = 0.4,
                                        type = "single_region",
                                        what = "SPR",
                                        par_subset = "fish_fixed_sel_pars",
                                        sd_rep = dusky_sd_rep,
                                        level = level)

  z <- stats::qnorm(1 - (1 - level) / 2) # 1.645 for a 90 pct interval
  est <- rp$refpts$est
  log_se <- rp$refpts$log_se

  expect_equal(rp$refpts$lwr, est * exp(-z * log_se), tolerance = 1e-10)
  expect_equal(rp$refpts$upr, est * exp(z * log_se), tolerance = 1e-10)
  expect_equal(rp$refpts$cv, sqrt(exp(log_se^2) - 1), tolerance = 1e-10) # CV of a lognormal

  # the fishing mortality reference point has a standard error, so its bounds and CV differ from the estimate
  expect_gt(max(log_se), 1e-6)

})

test_that("the drawn bounds are the quantiles of the draws", {

  level <- 0.80
  rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model,
                                        SPR_x = 0.4,
                                        type = "single_region",
                                        what = "SPR",
                                        method = "mvn",
                                        par_subset = "fish_fixed_sel_pars",
                                        sd_rep = dusky_sd_rep,
                                        n_draw = 12,
                                        seed = 7,
                                        level = level)

  # the draws are kept on the log scale, one column per quantity
  n_quantities <- ncol(rp$draws)
  expect_equal(nrow(rp$mvn), n_quantities)

  for(i in 1:n_quantities) {
    expect_equal(rp$mvn$lwr[i], exp(as.numeric(stats::quantile(rp$draws[,i], probs = 0.10))), tolerance = 1e-10)
    expect_equal(rp$mvn$upr[i], exp(as.numeric(stats::quantile(rp$draws[,i], probs = 0.90))), tolerance = 1e-10)
    expect_equal(rp$mvn$log_sd[i], stats::sd(rp$draws[,i]), tolerance = 1e-10)
  } # end i loop

  # the bounds are apart for the fishing mortality reference point
  f_row <- which(rp$mvn$quantity == "f_ref_pt")
  expect_gt(rp$mvn$upr[f_row], 1.001 * rp$mvn$lwr[f_row])

})

test_that("spawn timing left unset is read from the fitted model", {

  # the same fit, read as spawning 30 pct of the way into the year
  dusky_mid_year <- dusky_rtmb_model
  dusky_mid_year$data$t_spawn <- 0.3

  rp <- Get_Reference_Point_Uncertainty(obj = dusky_mid_year,
                                        SPR_x = 0.4,
                                        type = "single_region",
                                        what = "SPR",
                                        par_subset = "fish_fixed_sel_pars",
                                        sd_rep = dusky_sd_rep)
  mid_year <- Get_Reference_Points(data = dusky_mid_year$data, rep = dusky_mid_year$rep, SPR_x = 0.4, t_spawn = 0.3,
                                   type = "single_region", what = "SPR")
  year_start <- Get_Reference_Points(data = dusky_mid_year$data, rep = dusky_mid_year$rep, SPR_x = 0.4, t_spawn = 0,
                                     type = "single_region", what = "SPR")

  f_est <- rp$refpts$est[rp$refpts$quantity == "f_ref_pt"]
  expect_equal(f_est, as.numeric(mid_year$f_ref_pt), tolerance = 1e-8)
  expect_gt(abs(f_est / as.numeric(year_start$f_ref_pt) - 1), 0.005) # the two timings give different F40

})
