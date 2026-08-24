# Francis weighting for conditional age-at-length. A CAAL row is the age
# composition of the fish aged from one length bin, so a fleet's samples are the
# (year, season, length bin) rows and the weight pools the standardized mean-age
# residuals over all of them.

library(SPoRC)
library(testthat)
data("mlt_rg_goa_rex_data")

test_that("CAAL proportions are the age composition within each length bin", {

  dat <- mlt_rg_goa_rex_data
  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  props <- get_caal_prop(input_list$data, fit$rep)

  # observed ages, not model ages: the ageing error has already been applied
  expect_equal(dim(props$Pred_Srv_caal)[5], dim(input_list$data$ObsSrv_caal)[5])
  expect_null(props$Pred_Fish_caal) # this model has survey CAAL only

  # every row a fleet uses is a composition, so it sums to one within its length bin
  used <- which(input_list$data$UseSrv_caal == 1, arr.ind = TRUE)
  for(i in c(1, round(nrow(used) / 2), nrow(used))) {
    idx <- used[i, ]
    for(s in 1:input_list$data$n_sexes) {
      row <- props$Pred_Srv_caal[idx[1], idx[2], idx[3], idx[4], , s, idx[5]]
      expect_equal(sum(row), 1, tolerance = 1e-10)
      expect_true(all(row >= 0))
    } # end s loop
  } # end i loop

  # and rows nothing was aged from are left alone
  unused <- which(input_list$data$UseSrv_caal == 0, arr.ind = TRUE)[1, ]
  expect_true(all(is.na(props$Pred_Srv_caal[unused[1], unused[2], unused[3], unused[4], , , unused[5]])))
})


test_that("the CAAL weight is one pooled inverse variance per fleet, region and sex", {

  dat <- mlt_rg_goa_rex_data
  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)
  props <- get_caal_prop(input_list$data, fit$rep)

  n_bins <- dim(props$Pred_Srv_caal)[5]
  wts <- input_list$data$Wt_Srv_caal
  wts[] <- NA
  info <- get_francis_weights_caal(
    n_regions = input_list$data$n_regions, n_sexes = input_list$data$n_sexes,
    n_fleets = input_list$data$n_srv_fleets, n_years = length(input_list$data$years),
    n_seas = input_list$data$n_seas, n_lens = dim(input_list$data$UseSrv_caal)[4],
    Use = input_list$data$UseSrv_caal, ISS = input_list$data$ISS_Srv_caal,
    Pred_array = props$Pred_Srv_caal, Obs_array = props$Obs_Srv_caal,
    weights = wts, bins = seq_len(n_bins), comp_type = input_list$data$Srv_caal_Type
  )

  # these compositions are split by region and sex, so the weight is constant
  # over years and length bins within each region, sex and fleet
  w <- info$weights
  expect_true(all(is.finite(w[!is.na(w)])))
  expect_true(all(w[!is.na(w)] > 0))

  # only the fleets that aged fish get a weight
  has_caal <- which(apply(input_list$data$UseSrv_caal, 5, sum) > 0)
  expect_gt(length(has_caal), 0)

  # a cell that aged fish carries one weight across every year and length bin
  for(f in has_caal) {
    any_cell <- FALSE
    for(r in 1:input_list$data$n_regions) {
      for(s in 1:input_list$data$n_sexes) {
        cell <- w[r, , , , s, f]
        if(all(is.na(cell))) next
        expect_equal(length(unique(cell[!is.na(cell)])), 1)
        any_cell <- TRUE
      } # end s loop
    } # end r loop
    expect_true(any_cell)
  } # end f loop
  for(f in setdiff(1:input_list$data$n_srv_fleets, has_caal)) expect_true(all(is.na(w[, , , , , f])))

  # the weight against the definition, pooled by hand over every row used
  f <- has_caal[1]
  cells <- which(!is.na(w[, 1, 1, 1, , f]), arr.ind = TRUE)
  r <- cells[1, 1]
  s <- cells[1, 2]
  res <- c()
  for(y in 1:length(input_list$data$years)) {
    for(l in 1:dim(input_list$data$UseSrv_caal)[4]) {
      if(input_list$data$UseSrv_caal[r, y, 1, l, f] != 1) next
      p_row <- props$Pred_Srv_caal[r, y, 1, l, , s, f]
      o_row <- props$Obs_Srv_caal[r, y, 1, l, , s, f]
      n_row <- input_list$data$ISS_Srv_caal[r, y, 1, l, s, f]
      if(is.na(n_row) || n_row <= 0 || sum(o_row, na.rm = TRUE) == 0) next
      e_bar <- sum(seq_len(n_bins) * p_row)
      o_bar <- sum(seq_len(n_bins) * o_row)
      v <- sum(seq_len(n_bins)^2 * p_row) - e_bar^2
      if(!is.finite(v) || v <= 0) next
      res <- c(res, (o_bar - e_bar) / sqrt(v / n_row))
    } # end l loop
  } # end y loop

  expect_gt(length(res), 50) # the pooling really is over many rows, not a few
  expect_equal(unname(w[r, which(!is.na(w[r, , 1, 1, s, f]))[1], 1, 1, s, f]),
               1 / stats::var(res), tolerance = 1e-10)
})


test_that("do_francis_reweighting returns CAAL weights alongside the rest", {

  dat <- mlt_rg_goa_rex_data
  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)

  n_obs_ages <- dim(input_list$data$ObsSrvAgeComps)[4]
  wts <- do_francis_reweighting(data = input_list$data, rep = fit$rep,
                                age_labels = 1:n_obs_ages, len_labels = input_list$data$lens,
                                year_labels = input_list$data$years)

  expect_true("new_srv_caal_wts" %in% names(wts))
  expect_true("new_fish_caal_wts" %in% names(wts))
  expect_true(any(!is.na(wts$new_srv_caal_wts)))
  expect_gt(nrow(wts$mean_francis_caal), 0)
  expect_true(all(c("Len_Bin", "obs", "pred", "Fleet") %in% names(wts$mean_francis_caal)))

  # a fleet type with no CAAL keeps the weights it had, so applying them is a no-op
  expect_equal(wts$new_fish_caal_wts, input_list$data$Wt_Fish_caal)
})


test_that("get_caal_fits returns one row per aged length bin, with the Francis residual", {

  dat <- mlt_rg_goa_rex_data
  input_list <- suppressWarnings(suppressMessages(build_goa_rex_input(dat)))
  fit <- fit_model(input_list$data, input_list$par, input_list$map, do_optim = FALSE, silent = TRUE)

  fits <- get_caal_fits(input_list$data, fit$rep)
  expect_gt(nrow(fits), 0)
  expect_true(all(c("Region", "Year", "Len_Bin", "Length", "Sex", "Fleet", "Type",
                    "obs", "pred", "sd_pred", "ISS") %in% names(fits)))
  expect_equal(unique(fits$Type), "Survey") # this model ages fish at the survey only
  expect_true(all(fits$Year %in% input_list$data$years))
  expect_true(all(fits$Length %in% input_list$data$lens))

  # a mean age has to sit inside the age bins it is a mean over
  n_bins <- dim(input_list$data$ObsSrv_caal)[5]
  expect_true(all(fits$obs >= 1 & fits$obs <= n_bins))
  expect_true(all(fits$pred >= 1 & fits$pred <= n_bins))

  # longer fish are older, in the data and in the model. Rank correlation rather
  # than Pearson, because mean age flattens at the longest bins as they fill with
  # the plus group
  by_len <- stats::aggregate(cbind(obs, pred) ~ Length, data = fits, FUN = mean)
  expect_gt(stats::cor(by_len$Length, by_len$obs, method = "spearman"), 0.9)
  expect_gt(stats::cor(by_len$Length, by_len$pred, method = "spearman"), 0.9)
  # and the shortest bins hold much younger fish than the longest
  short <- utils::head(by_len, 5)
  long <- utils::tail(by_len, 5)
  expect_gt(mean(long$obs) - mean(short$obs), 2)
  expect_gt(mean(long$pred) - mean(short$pred), 2)

  # residuals are not here: they come from get_osa and plot_resids
  expect_false("resid" %in% names(fits))
})
