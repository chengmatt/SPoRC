# With random effects, fit_model adds the marginal and the laplace term to obj$rep and get_nLL_plot shows the
# marginal as the total. Without random effects neither is added and the plot shows jnLL as before.

library(SPoRC)
library(testthat)

data("sgl_rg_dusky_data")

test_that("with random effects, obj$rep has the marginal the fit minimized and the plot shows it", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = c("ln_RecDevs", "ln_InitDevs"),
                   newton_loops = 0, silent = TRUE)

  # marginal from the fit, and the laplace term is the rest
  expect_equal(fit$rep$marginal_nLL, as.numeric(fit$optim$objective), tolerance = 1e-6)
  expect_equal(fit$rep$laplace_nLL, fit$rep$marginal_nLL - fit$rep$jnLL)
  expect_true(abs(fit$rep$laplace_nLL) > 0.1) # the marginal and jnLL at the mode differ here

  # plot total is the marginal, with the laplace term next to the components
  plot_df <- get_nLL_plot(list(input_list$data), list(fit$rep), "RE")[[1]]$data
  expect_equal(plot_df$value[plot_df$name == "Marginal nLL"], fit$rep$marginal_nLL)
  expect_equal(plot_df$value[plot_df$name == "Laplace Term"], fit$rep$laplace_nLL)
  expect_false("jnLL" %in% plot_df$name)

})

test_that("without random effects, obj$rep has no marginal and the plot shows jnLL", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL,
                   newton_loops = 0, silent = TRUE)

  expect_null(fit$rep$marginal_nLL)
  expect_null(fit$rep$laplace_nLL)

  plot_df <- get_nLL_plot(list(input_list$data), list(fit$rep), "Fixed")[[1]]$data
  expect_equal(plot_df$value[plot_df$name == "jnLL"], fit$rep$jnLL)
  expect_false(any(c("Marginal nLL", "Laplace Term") %in% plot_df$name))

})

test_that("with random effects but no optimization, obj$rep has no marginal", {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  fit <- fit_model(input_list$data, input_list$par, input_list$map, random = c("ln_RecDevs", "ln_InitDevs"),
                   do_optim = FALSE, silent = TRUE)

  expect_null(fit$rep$marginal_nLL)
  expect_null(fit$rep$laplace_nLL)

})
