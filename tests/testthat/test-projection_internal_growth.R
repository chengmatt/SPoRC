# Weight at age in the projection read from the fit's projected growth, on the Pacific cod bridge whose growth
# is derived and time-varying: the internal projection matches the post hoc one with the projected years copied
# in (1e-6), the average and propagated options differ, and the projected quantities have a finite gradient.
# The bridge's Hessian is not positive definite, so no standard error is checked here.

library(SPoRC)
library(testthat)
data("sgl_rg_ebs_pcod_data")

n_proj_yrs <- 6
build_ext <- build_ebs_pcod_input
body(build_ext) <- do.call(substitute, list(body(build_ext), list(Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 5)))))
il <- suppressWarnings(seed_ebs_pcod_mle(suppressMessages(build_ext(sgl_rg_ebs_pcod_data)), sgl_rg_ebs_pcod_data))
il$data$growth_tv_type <- 0 # curve mode builds every year up front, which is what projected years need
n_yrs <- length(il$data$years)

# the model at the seeded estimates, with the projection on and nothing re-optimized
build_at_seed <- function(growth_opt, sel_opt) {
  ip <- Setup_Mod_Projection(il, n_proj_yrs = n_proj_yrs, f_opt = "terminal", growth_opt = growth_opt, sel_opt = sel_opt, n_avg_yrs = 1)
  fit_model(ip$data, ip$par, ip$map, do_optim = FALSE, silent = TRUE)
}

test_that("growth_opt = 'propagate' is refused where there is no growth model to propagate", {

  data("sgl_rg_dusky_data", envir = environment())
  dusky <- build_goa_dusky_input(sgl_rg_dusky_data)
  expect_error(Setup_Mod_Projection(dusky, n_proj_yrs = n_proj_yrs, rec_opt = "mean", growth_opt = "propagate"), "waa_model")

})

test_that("projected weight at age and selectivity read the years the fit extended", {

  fit_avg <- build_at_seed("average", "average")
  fit_prj <- build_at_seed("propagate", "propagate")

  # the fit extends growth and selectivity over the projection years
  expect_equal(dim(fit_prj$rep$WAA)[3], n_yrs + 5)
  expect_equal(dim(fit_prj$rep$fish_sel)[3], n_yrs + 5)

  # the post hoc projection with the projected years copied in by hand
  rec_devs <- array(1, dim = c(1, 1, n_proj_yrs))
  rec_devs[1,1,2:n_proj_yrs] <- exp(fit_prj$env$parList()$ln_RecDevs[1,1,n_yrs + 1:(n_proj_yrs - 1)])
  f_term <- sum(fit_prj$rep$Fmort[1,n_yrs,,])
  inputs <- build_projection_inputs(data = fit_prj$data, rep = fit_prj$rep, refpts = list(f_ref_pt = f_term, b_ref_pt = array(1, dim = c(1, 1))), n_proj_yrs = n_proj_yrs)
  inputs$ret_sel <- array(0, dim = dim(inputs$fish_sel)) # the builder leaves retention at its default, the model has its own
  for(y in 1:n_proj_yrs) inputs$ret_sel[,,y,,,,] <- fit_prj$rep$ret_sel[,,n_yrs,,,,]
  inputs_prj <- inputs
  for(y in 1:n_proj_yrs) {
    inputs_prj$WAA[,,y,,,] <- fit_prj$rep$WAA[,,n_yrs + y - 1,,,]
    inputs_prj$WAA_fish[,,y,,,,] <- fit_prj$rep$WAA_fish[,,n_yrs + y - 1,,,,]
    inputs_prj$fish_sel[,,y,,,,] <- fit_prj$rep$fish_sel[,,n_yrs + y - 1,,,,]
    inputs_prj$ret_sel[,,y,,,,] <- fit_prj$rep$ret_sel[,,n_yrs + y - 1,,,,]
  } # end y loop
  ext_avg <- do.call(Do_Population_Projection, c(inputs, list(fmort_opt = "Input", rec_devs = rec_devs)))
  ext_prj <- do.call(Do_Population_Projection, c(inputs_prj, list(fmort_opt = "Input", rec_devs = rec_devs)))

  expect_equal(as.numeric(fit_avg$rep$proj_SSB[1,1,]), as.numeric(ext_avg$proj_SSB[1,1,]), tolerance = 1e-6)
  expect_equal(as.numeric(fit_prj$rep$proj_SSB[1,1,]), as.numeric(ext_prj$proj_SSB[1,1,]), tolerance = 1e-6)

  # the iid growth deviations project at their mean, below the terminal year, so the two options part from year two
  expect_equal(fit_avg$rep$proj_SSB[1,1,1], fit_prj$rep$proj_SSB[1,1,1], tolerance = 1e-8)
  expect_gt(abs(fit_avg$rep$proj_SSB[1,1,3] - fit_prj$rep$proj_SSB[1,1,3]) / fit_avg$rep$proj_SSB[1,1,3], 0.05)

  # projected numbers at age are reported for the projection years
  expect_equal(dim(fit_prj$rep$proj_NAA)[3], n_proj_yrs + 1)
  expect_true(all(is.finite(fit_prj$rep$proj_NAA)))

})

test_that("the projected quantities have a finite gradient in every parameter, growth deviations included", {

  fit_prj <- build_at_seed("propagate", "propagate")

  # with the identity in place of the covariance, sdreport's sd is the length of each quantity's gradient
  sdr <- RTMB::sdreport(fit_prj, hessian.fixed = diag(length(fit_prj$par)))
  grad_len <- sdr$sd[names(sdr$value) == "log_proj_SSB"]
  expect_true(all(is.finite(grad_len)))
  expect_true(all(grad_len > 0))
  expect_true(all(is.finite(sdr$sd[names(sdr$value) == "log_proj_NAA"])))

  # and the projected growth deviations are among the parameters the projection reads
  fit_avg <- build_at_seed("average", "average")
  sdr_avg <- RTMB::sdreport(fit_avg, hessian.fixed = diag(length(fit_avg$par)))
  grad_len_avg <- sdr_avg$sd[names(sdr_avg$value) == "log_proj_SSB"]
  expect_equal(grad_len[1], grad_len_avg[1], tolerance = 1e-8) # year one is the terminal year either way
  expect_false(isTRUE(all.equal(grad_len[-1], grad_len_avg[-1], tolerance = 1e-3)))

})
