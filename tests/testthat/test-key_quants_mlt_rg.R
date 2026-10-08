# get_key_quants on a multi-region model. Averaging over years drops the year dim of movement, which
# has to be put back as [pop, from, to, year, seas, age, sex] for the projection to index it.

library(SPoRC)
library(testthat)

test_that("get_key_quants projects a multi-region model under global SPR", {

  # small two-region population, evaluated rather than fitted (helper-collapse.R)
  il <- collapse_input(nr = 2, nx = 2)
  data <- il$data
  rep <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$rep
  n_regions <- data$n_regions
  n_yrs <- length(data$years)

  HCR_function <- function(x, frp, brp, alpha = 0.05) {
    stock_status <- x / brp
    if(stock_status >= 1) f <- frp
    if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
    if(stock_status < alpha) f <- 0
    return(f)
  }

  reference_points_opt <- list(SPR_x = 0.4,
                               t_spawn = 0,
                               sex_ratio_f = array(0.5, dim = c(data$n_pop, n_regions)),
                               calc_rec_st_yr = 1,
                               rec_age = 1,
                               type = "multi_region",
                               what = "global_SPR")

  # one and several averaging years, since the year dim is dropped and rebuilt either way
  for(n_avg_yrs in c(1, 3)) {
    proj_model_opt <- list(n_proj_yrs = 2,
                           n_avg_yrs = n_avg_yrs,
                           HCR_function = HCR_function,
                           recruitment_opt = "mean_rec",
                           fmort_opt = "HCR")

    kq <- get_key_quants(data = list(data),
                         rep = list(rep),
                         reference_points_opt = reference_points_opt,
                         proj_model_opt = proj_model_opt,
                         model_names = "mlt_rg")[[1]]

    expect_equal(nrow(kq), n_regions)
    expect_true(all(is.finite(kq$Catch_Advice) & kq$Catch_Advice > 0))
    # year 1 of the projection replays the terminal year
    expect_equal(kq$Terminal_SSB, round(apply(rep$SSB[,,n_yrs, drop = FALSE], 2, sum), 5))
  }
})
