# get_key_quants on a multi-region model. Averaging over years drops the year dim of movement, which
# has to be put back as [pop, from, to, year, seas, age, sex] for the projection to index it, and the
# projection has to run under the fmort_opt asked for.

library(SPoRC)
library(testthat)

# small two-region population, evaluated rather than fitted (helper-collapse.R)
il <- collapse_input(nr = 2, nx = 2)
kq_data <- il$data
kq_rep <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep

kq_ref_pts_opt <- list(SPR_x = 0.4,
                       t_spawn = 0,
                       sex_ratio_f = array(0.5, dim = c(kq_data$n_pop, kq_data$n_regions)),
                       calc_rec_st_yr = 1,
                       rec_age = 1,
                       type = "multi_region",
                       what = "global_SPR")

run_kq <- function(HCR_function, fmort_opt, n_avg_yrs = 1) {
  get_key_quants(data = list(kq_data),
                 rep = list(kq_rep),
                 reference_points_opt = kq_ref_pts_opt,
                 proj_model_opt = list(n_proj_yrs = 2,
                                       n_avg_yrs = n_avg_yrs,
                                       HCR_function = HCR_function,
                                       recruitment_opt = "mean_rec",
                                       fmort_opt = fmort_opt),
                 model_names = "mlt_rg")[[1]]
}

test_that("get_key_quants projects a multi-region model under global SPR", {

  HCR_function <- function(x, frp, brp, alpha = 0.05) {
    stock_status <- x / brp
    if(stock_status >= 1) f <- frp
    if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
    if(stock_status < alpha) f <- 0
    return(f)
  }

  # one and several averaging years, since the year dim is dropped and rebuilt either way
  for(n_avg_yrs in c(1, 3)) {
    kq <- run_kq(HCR_function, "HCR_global", n_avg_yrs)

    expect_equal(nrow(kq), kq_data$n_regions)
    expect_true(all(is.finite(kq$Catch_Advice) & kq$Catch_Advice > 0))
    # year 1 of the projection replays the terminal year
    expect_equal(kq$Terminal_SSB, round(apply(kq_rep$SSB[,,length(kq_data$years), drop = FALSE], 2, sum), 5))
  }
})

test_that("get_key_quants projects under the fmort_opt it is given", {

  # a rule that always halves F, so applying it shows up in the catch whatever the stock status
  half_F <- function(x, frp, brp) frp / 2

  catch_input <- run_kq(half_F, "Input")$Catch_Advice
  catch_hcr <- run_kq(half_F, "HCR")$Catch_Advice
  catch_hcr_global <- run_kq(half_F, "HCR_global")$Catch_Advice

  expect_true(all(catch_hcr < catch_input))
  expect_equal(catch_hcr, catch_hcr_global)
})
