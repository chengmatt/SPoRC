# The years mean recruitment is averaged over for SPR biomass reference points: calc_rec_yrs against the
# calc_rec_st_yr / rec_age window, the mean_rec it reports, and refusal of years outside the model.

library(SPoRC)
library(testthat)
data("sgl_rg_sable_rep")
data("sgl_rg_sable_data")
data("mlt_rg_sable_rep")
data("mlt_rg_sable_data")

n_years_sgl <- length(sgl_rg_sable_data$years)

refpt_sgl <- function(...) {
  Get_Reference_Points(data = sgl_rg_sable_data,
                       rep = sgl_rg_sable_rep,
                       SPR_x = 0.4,
                       type = "single_region",
                       what = "SPR",
                       ...)
}


test_that("calc_rec_yrs over the default window gives the same reference points", {
  window <- refpt_sgl(calc_rec_st_yr = 20, rec_age = 2)
  yrs <- refpt_sgl(calc_rec_yrs = 20:(n_years_sgl - 2))

  for(x in c("f_ref_pt", "b_ref_pt", "virgin_b_ref_pt", "pop_b_ref_pt", "virgin_pop_b_ref_pt", "mean_rec")) {
    expect_identical(yrs[[x]], window[[x]], label = x)
  }
})


test_that("mean_rec is the mean over the years asked for, and B40 scales with it", {
  base <- refpt_sgl(calc_rec_st_yr = 20, rec_age = 2)

  # every year including the terminal one, and a set with gaps in it
  for(rec_yrs in list(1:n_years_sgl, c(5, 10, 20:30))) {
    ref <- refpt_sgl(calc_rec_yrs = rec_yrs)

    expect_equal(dim(ref$mean_rec), c(1, 1))
    expect_equal(ref$mean_rec[1, 1], mean(sgl_rg_sable_rep$Rec[1, 1, rec_yrs]), tolerance = 1e-12)

    # F40 and spawning biomass per recruit do not depend on recruitment, so only the scaling moves
    expect_equal(ref$f_ref_pt, base$f_ref_pt, tolerance = 1e-12)
    expect_equal(ref$b_ref_pt / ref$mean_rec[1, 1], base$b_ref_pt / base$mean_rec[1, 1], tolerance = 1e-10)
    expect_equal(ref$virgin_b_ref_pt / ref$mean_rec[1, 1], base$virgin_b_ref_pt / base$mean_rec[1, 1], tolerance = 1e-10)
  }
})


test_that("mean_rec is reported by population and region under every method", {
  n_years_mlt <- length(mlt_rg_sable_data$years)
  rec_yrs <- 20:(n_years_mlt - 2)
  expected <- apply(mlt_rg_sable_rep$Rec[, , rec_yrs, drop = FALSE], c(1, 2), mean)

  # independent SPR solves region by region, global MSY never scales by mean recruitment
  for(what in c("independent_SPR", "global_SPR", "global_MSY")) {
    ref <- Get_Reference_Points(data = mlt_rg_sable_data,
                                rep = mlt_rg_sable_rep,
                                SPR_x = 0.4,
                                type = "multi_region",
                                what = what,
                                calc_rec_yrs = rec_yrs)
    expect_equal(dim(ref$mean_rec), c(mlt_rg_sable_data$n_pop, mlt_rg_sable_data$n_regions), label = what)
    expect_equal(c(ref$mean_rec), c(expected), tolerance = 1e-12, label = what)
  }

  msy <- Get_Reference_Points(data = sgl_rg_sable_data,
                              rep = sgl_rg_sable_rep,
                              type = "single_region",
                              what = "MSY",
                              calc_rec_yrs = 1:n_years_sgl)
  expect_equal(msy$mean_rec[1, 1], mean(sgl_rg_sable_rep$Rec[1, 1, ]), tolerance = 1e-12)
})


test_that("recruitment years outside the model are refused", {
  # calendar years, a year zero, an empty set, and a default window that runs past the terminal year
  expect_error(refpt_sgl(calc_rec_yrs = 1977:2020), "model years")
  expect_error(refpt_sgl(calc_rec_yrs = 0:10), "model years")
  expect_error(refpt_sgl(calc_rec_yrs = integer(0)), "model years")
  expect_error(refpt_sgl(calc_rec_st_yr = 1, rec_age = -1), "model years")
})


test_that("get_key_quants projects from the same recruitment years as calc_rec_yrs", {
  il <- collapse_input(nr = 2, nx = 1)
  kq_data <- il$data
  kq_rep <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep
  n_yrs <- length(kq_data$years)
  kq_rep$Rec[, , n_yrs] <- 3 * kq_rep$Rec[, , n_yrs] # unfitted recruitment is flat, so give the terminal year a strong cohort

  HCR_function <- function(x, frp, brp, alpha = 0.05) {
    stock_status <- x / brp
    if(stock_status >= 1) f <- frp
    if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
    if(stock_status < alpha) f <- 0
    return(f)
  }

  run_kq <- function(rec_opt) {
    ref_opt <- c(list(SPR_x = 0.4,
                      t_spawn = 0,
                      sex_ratio_f = array(0.5, dim = c(kq_data$n_pop, kq_data$n_regions)),
                      calc_rec_st_yr = 1,
                      rec_age = 1,
                      type = "multi_region",
                      what = "global_SPR"),
                 rec_opt)
    get_key_quants(data = list(kq_data),
                   rep = list(kq_rep),
                   reference_points_opt = ref_opt,
                   proj_model_opt = list(n_proj_yrs = 2,
                                         n_avg_yrs = 1,
                                         HCR_function = HCR_function,
                                         recruitment_opt = "mean_rec",
                                         fmort_opt = "HCR_global"),
                   model_names = "mlt_rg")[[1]]
  }

  window <- run_kq(list())
  same_yrs <- run_kq(list(calc_rec_yrs = 1:(n_yrs - 1)))
  expect_identical(same_yrs, window)

  # a different set of years moves both the biomass reference point and the projected catch
  other_yrs <- run_kq(list(calc_rec_yrs = c(1, n_yrs)))
  expect_false(isTRUE(all.equal(other_yrs$B_Ref_Pt, window$B_Ref_Pt)))
  expect_false(isTRUE(all.equal(other_yrs$Catch_Advice, window$Catch_Advice)))
})
