# plot() on the results of do_retrospective and get_osa draws what get_retrospective_plot and plot_resids draw.
# The labels that make this work must leave the results reading as the data frame and list they were.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

data <- dusky_rtmb_model$data
parameters <- dusky_rtmb_model$parameters
mapping <- dusky_rtmb_model$mapping

test_that("plot on retrospective results matches get_retrospective_plot", {

  retro <- do_retrospective(n_retro = 2, data = data, parameters = parameters, mapping = mapping, random = NULL,
                            do_par = FALSE, newton_loops = 1)

  # labeled, and still a data frame with the same columns
  expect_s3_class(retro, "sporc_retro")
  expect_s3_class(retro, "data.frame")
  expect_true(all(c("Pop", "Region", "Year", "Type", "peel", "value") %in% names(retro)))

  plots_method <- plot(retro, Rec_Age = 4)
  plots_direct <- get_retrospective_plot(retro_output = retro, Rec_Age = 4)
  expect_equal(length(plots_method), 3)
  for(i in 1:3) {
    expect_s3_class(plots_method[[i]], "ggplot")
    expect_equal(plots_method[[i]]$data, plots_direct[[i]]$data)
  } # end i loop

  # the label does not change the relative differences, checked against the unlabeled data frame
  retro_plain <- as.data.frame(unclass(retro))
  expect_equal(as.data.frame(unclass(get_retrospective_relative_difference(retro))),
               as.data.frame(unclass(get_retrospective_relative_difference(retro_plain))))

})

test_that("with return_models = TRUE the labeled results are in retro_df", {

  retro <- do_retrospective(n_retro = 1, data = data, parameters = parameters, mapping = mapping, random = NULL,
                            do_par = FALSE, newton_loops = 1, return_models = TRUE)

  expect_named(retro, c("retro_df", "retro_models"))
  expect_s3_class(retro$retro_df, "sporc_retro")
  expect_equal(plot(retro$retro_df, Rec_Age = 4)[[1]]$data, get_retrospective_plot(retro$retro_df, 4)[[1]]$data)

})

test_that("plot on OSA residuals matches plot_resids", {

  skip_if_not_installed("compResidual")

  # multinomial ages drawn from a known composition, in the array get_osa reads: region, year, season, bin, sex, fleet
  set.seed(42)
  n_years <- 20
  n_bins <- 10
  true_prop <- rep(1 / n_bins, n_bins)
  obs_mat <- array(data = 0, dim = c(1, n_years, 1, n_bins, 1, 1))
  exp_mat <- array(data = 0, dim = c(1, n_years, 1, n_bins, 1, 1))
  for(y in 1:n_years) {
    obs_mat[1,y,1,,1,1] <- as.vector(rmultinom(1, 200, true_prop))
    exp_mat[1,y,1,,1,1] <- true_prop
  } # end y loop

  osa <- get_osa(obs_mat = obs_mat, exp_mat = exp_mat, N = rep(200, n_years), years = 1:n_years, seas = 1, fleet = 1,
                 bins = 1:n_bins, comp_type = 0, bin_label = "Ages")

  # labeled, and still a list holding the residual data frame
  expect_s3_class(osa, "sporc_osa")
  expect_true(is.list(osa))
  expect_s3_class(osa$res, "data.frame")

  plots_method <- plot(osa)
  plots_direct <- plot_resids(osa_results = osa)
  expect_equal(names(plots_method), names(plots_direct))
  expect_s3_class(plots_method[[1]], "ggplot")
  expect_equal(plots_method[[1]]$data, plots_direct[[1]]$data)
  expect_equal(plots_method[[2]]$data, plots_direct[[2]]$data)

})
