# Reference points read weight at age from the data, which is a zero placeholder
# when the growth module derives it. A model with waa_model = "wt_len" would then
# get an unfished biomass of zero and reference points to match, without a word.

library(SPoRC)
library(testthat)
data("sgl_rg_ebs_pcod_data")

test_that("reference points use the weight at age growth derived", {

  dat <- sgl_rg_ebs_pcod_data
  input_list <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(dat))), dat)
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)

  # the model derives weight at age, so the data list holds a placeholder only
  expect_true(all(input_list$data$WAA == 0))
  expect_false(is.null(fit$rep$WAA))
  expect_true(any(fit$rep$WAA > 0))

  spr <- Get_Reference_Points(data = input_list$data, rep = fit$rep, SPR_x = 0.4,
                              t_spawn = input_list$data$t_spawn, calc_rec_st_yr = 1,
                              rec_age = 1, type = "single_region", what = "SPR",
                              n_avg_yrs = 1)

  # unfished biomass per recruit is built from a real weight at age, not zeros
  expect_true(is.finite(spr$virgin_b_ref_pt[1, 1]))
  expect_gt(spr$virgin_b_ref_pt[1, 1], 0)
  expect_gt(spr$b_ref_pt[1, 1], 0)
  expect_gt(spr$f_ref_pt[1], 0)

  # and the target sits below the unfished level, as a 40 percent SPR target must
  expect_lt(spr$b_ref_pt[1, 1], spr$virgin_b_ref_pt[1, 1])

  # feeding the data placeholder in instead collapses the unfished biomass, which
  # is what the reference points did before they read the report
  rep_zero <- fit$rep
  rep_zero$WAA <- NULL
  # the solve walks into NaN with no weight at age, which is the point
  spr_zero <- suppressWarnings(Get_Reference_Points(data = input_list$data, rep = rep_zero, SPR_x = 0.4,
                                                    t_spawn = input_list$data$t_spawn, calc_rec_st_yr = 1,
                                                    rec_age = 1, type = "single_region", what = "SPR",
                                                    n_avg_yrs = 1))
  expect_equal(spr_zero$virgin_b_ref_pt[1, 1], 0)
})
