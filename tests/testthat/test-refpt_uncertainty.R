library(SPoRC)
library(testthat)
data("dusky_rtmb_model")

# Dusky is single region, single pop, single season, mean recruitment, with M and
# steepness mapped off. Several sensitivities are therefore known analytically.

test_that("Get_Reference_Point_Uncertainty returns intervals with the expected structure", {

  rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model,
                                        SPR_x = 0.4,
                                        type = "single_region",
                                        what = "SPR",
                                        par_subset = "fish_fixed_sel_pars")

  expect_s3_class(rp$refpts, "data.frame")
  expect_true(all(c("quantity", "est", "log_se", "lwr", "upr", "cv") %in% names(rp$refpts)))
  expect_true(all(c("f_ref_pt", "b_ref_pt", "virgin_b_ref_pt") %in% rp$refpts$quantity))

  # log scale keeps bounds positive and bracketing
  expect_true(all(rp$refpts$lwr > 0))
  expect_true(all(rp$refpts$lwr <= rp$refpts$est))
  expect_true(all(rp$refpts$upr >= rp$refpts$est))

  # Only selectivity was perturbed here, and an SPR target pins SBPR(F_x) to a fixed
  # fraction of SBPR(0). Since SBPR(0) is evaluated at F = 0, neither biomass
  # reference point can respond to selectivity at all; only F_x itself moves. The
  # residual is rounding, and comes out as an exact zero on some platforms.
  se <- stats::setNames(rp$refpts$log_se, rp$refpts$quantity)
  expect_equal(se[["virgin_b_ref_pt"]], 0, tolerance = 1e-10)
  expect_equal(se[["b_ref_pt"]], 0, tolerance = 1e-10)
  expect_gt(se[["f_ref_pt"]], 1e-6)

  # point estimate must match Get_Reference_Points
  base <- Get_Reference_Points(data = dusky_rtmb_model$data,
                               rep = dusky_rtmb_model$rep,
                               SPR_x = 0.4,
                               type = "single_region",
                               what = "SPR")
  expect_equal(rp$refpts$est[rp$refpts$quantity == "f_ref_pt"],
               as.numeric(base$f_ref_pt), tolerance = 1e-8)

  expect_equal(dim(rp$d), c(nrow(rp$refpts), length(dusky_rtmb_model$env$last.par.best)))
  expect_equal(dim(rp$log_cov), c(nrow(rp$refpts), nrow(rp$refpts)))

})

test_that("parameters that cannot enter a per-recruit reference point give exactly zero", {

  rp <- Get_Reference_Point_Uncertainty(
    obj = dusky_rtmb_model,
    SPR_x = 0.4,
    type = "single_region",
    what = "SPR",
    par_subset = c("ln_global_R0", "ln_RecDevs",
                                                       "ln_F_devs", "fish_fixed_sel_pars")
  )

  p_names <- colnames(rp$d)
  f_row <- rp$d["f_ref_pt", ]

  # F40 is per-recruit, so recruitment scale cannot move it. With one fleet and one
  # season the fleet F ratio is one whatever the F devs do.
  expect_true(all(f_row[p_names == "ln_global_R0"] == 0))
  expect_true(all(f_row[p_names == "ln_RecDevs"] == 0))
  expect_true(all(f_row[p_names == "ln_F_devs"] == 0))

  # selectivity is all that is left
  expect_true(any(abs(f_row[p_names == "fish_fixed_sel_pars"]) > 1e-6))

  # B40 scales SBPR by mean recruitment, so it does respond to the rec devs
  expect_true(any(abs(rp$d["b_ref_pt", p_names == "ln_RecDevs"]) > 1e-6))

})

test_that("MSY sensitivity to unfished recruitment matches the Beverton-Holt algebra", {

  # Switch onto the Beverton-Holt branch so R0 enters. The yield curve scales with R0
  # without moving its peak, so Bmsy is proportional to it and Fmsy independent of it.
  obj_bh <- dusky_rtmb_model
  obj_bh$data$rec_model <- 1

  rp <- Get_Reference_Point_Uncertainty(obj = obj_bh,
                                        type = "single_region",
                                        what = "MSY",
                                        par_subset = "ln_global_R0")

  j <- which(colnames(rp$d) == "ln_global_R0")

  expect_equal(as.numeric(rp$d["b_ref_pt", j]), 1, tolerance = 1e-4)
  expect_equal(as.numeric(rp$d["f_ref_pt", j]), 0, tolerance = 1e-8)

})

test_that("sensitivities do not depend on the finite difference step size", {

  d_at <- function(step) {
    rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model,
                                          SPR_x = 0.4,
                                          type = "single_region",
                                          what = "SPR",
                                          par_subset = "fish_fixed_sel_pars",
                                          rel_step = step)
    rp$d["f_ref_pt", colnames(rp$d) == "fish_fixed_sel_pars"]
  }

  # Wide plateau between truncation and cancellation error. A failure here means the
  # inner solve has stopped converging tightly.
  expect_equal(d_at(1e-2), d_at(1e-3), tolerance = 1e-5)
  expect_equal(d_at(1e-3), d_at(1e-4), tolerance = 1e-5)

})

test_that("delta method agrees with re-solving the reference point at drawn parameters", {

  status <- function(rep, refpts) {
    ssb <- rep$SSB[1, 1, dim(rep$SSB)[3]]
    c(SSB_terminal = ssb, status = ssb / as.numeric(refpts$b_ref_pt))
  }

  rp <- Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model,
                                        SPR_x = 0.4,
                                        type = "single_region",
                                        what = "SPR",
                                        extra_quantities = status,
                                        method = "both",
                                        n_draw = 150,
                                        seed = 42)

  expect_s3_class(rp$mvn, "data.frame")
  expect_equal(nrow(rp$mvn), nrow(rp$refpts))

  # close but not exact, since only the delta method linearizes
  ratio <- rp$mvn$log_sd / rp$refpts$log_se
  expect_true(all(ratio > 0.75 & ratio < 1.35))

  # status is better determined than either part, since the same parameters move both
  se <- stats::setNames(rp$refpts$log_se, rp$refpts$quantity)
  naive <- sqrt(se[["SSB_terminal"]]^2 + se[["b_ref_pt"]]^2)
  expect_lt(se[["status"]], naive)
  expect_gt(naive / se[["status"]], 1.5)

})

test_that("a report list instead of a fitted object gives a clear error", {

  expect_error(
    Get_Reference_Point_Uncertainty(obj = dusky_rtmb_model$rep,
                                    SPR_x = 0.4,
                                    type = "single_region",
                                    what = "SPR"),
    "fitted model object"
  )

})

