# Checks the sporc_fit class fit_model puts on its return: the fit still reads as a list, print shows the fit rather than
# the list, summary matches the sdreport, logLik gives the same AIC as marg_AIC, and plot matches the plotting functions.

library(testthat)
library(SPoRC)
data("dusky_rtmb_model")

data <- dusky_rtmb_model$data
parameters <- dusky_rtmb_model$parameters
mapping <- dusky_rtmb_model$mapping

fit <- suppressMessages(fit_model(data, parameters, mapping, newton_loops = 1, silent = TRUE, do_sdrep = TRUE))

test_that("the fit has the class and still reads as a list", {

  expect_s3_class(fit, "sporc_fit")
  expect_true(is.list(fit))
  expect_true(all(c("par", "fn", "gr", "env", "optim", "rep", "sd_rep", "data", "parameters", "mapping") %in% names(fit)))

  # the class does not get in the way of RTMB reading the fit
  sd_rep <- RTMB::sdreport(fit)
  expect_equal(sd_rep$par.fixed, fit$sd_rep$par.fixed)
  expect_equal(fit$fn(fit$optim$par), fit$optim$objective)

})

test_that("print shows the fit and returns it invisibly", {

  printed <- capture.output(returned <- withVisible(print(fit)))
  expect_false(returned$visible)
  expect_identical(returned$value, fit)

  expect_true(length(printed) < 12) # a short description, not the whole list
  expect_true(any(grepl(paste(min(data$years), "to", max(data$years)), printed)))
  expect_true(any(grepl(paste(length(fit$par), "fixed effects, no random effects"), printed)))
  expect_true(any(grepl("relative convergence", printed)))
  expect_true(any(grepl("jnLL:", printed)))
  max_grad_par <- names(fit$par)[which.max(abs(fit$optim$gradient))] # parameter the largest gradient belongs to
  expect_true(any(grepl(paste0("max gradient:.*\\(", max_grad_par, "\\)"), printed)))
  expect_true(any(grepl("positive definite", printed)))
  expect_true(any(grepl("sanity checks:  passed", printed)))

})

test_that("summary has one row per parameter name and matches the sdreport", {

  par_table <- summary(fit)
  expect_s3_class(par_table, "data.frame")
  expect_named(par_table, c("Par", "Type", "N", "Est_Min", "Est_Max", "SE_Max", "Abs_Grad_Max"))
  expect_equal(par_table$Par, unique(names(fit$sd_rep$par.fixed)))
  expect_equal(sum(par_table$N), length(fit$par))
  expect_true(all(par_table$Type == "fixed"))

  # each row against the sdreport's own values for that parameter
  se_fixed <- sqrt(diag(fit$sd_rep$cov.fixed))
  for(i in seq_len(nrow(par_table))) {
    par_idx <- which(names(fit$sd_rep$par.fixed) == par_table$Par[i])
    expect_equal(par_table$N[i], length(par_idx))
    expect_equal(par_table$Est_Min[i], min(fit$sd_rep$par.fixed[par_idx]))
    expect_equal(par_table$Est_Max[i], max(fit$sd_rep$par.fixed[par_idx]))
    expect_equal(par_table$SE_Max[i], max(se_fixed[par_idx]))
    expect_equal(par_table$Abs_Grad_Max[i], max(abs(fit$sd_rep$gradient.fixed[par_idx])))
  } # end i loop

})

test_that("logLik gives the AIC marg_AIC gives", {

  loglik <- logLik(fit)
  expect_s3_class(loglik, "logLik")
  expect_equal(as.numeric(loglik), -fit$optim$objective)
  expect_equal(attr(loglik, "df"), length(fit$par))
  expect_equal(AIC(fit), marg_AIC(fit$optim))

  # two fits in one call, AIC names its rows by object so the second needs its own name
  fit_copy <- fit
  aic_table <- AIC(fit, fit_copy)
  expect_equal(aic_table$AIC, rep(marg_AIC(fit$optim), 2))
  expect_equal(aic_table$df, rep(length(fit$par), 2))

})

test_that("plot draws what the plotting function draws for this fit alone", {

  # each type against the function it calls, given the fit as a list of one model
  expect_equal(plot(fit, type = "selex")[[1]]$data, get_selex_plot(list(fit$rep), "Model 1")[[1]]$data)
  expect_equal(plot(fit, type = "idx_fits")$data, get_idx_fits_plot(list(fit$data), list(fit$rep), "Model 1")$data)
  expect_equal(plot(fit, type = "data_fitted")$data, get_data_fitted_plot(list(fit$data), "Model 1")$data)
  expect_equal(plot(fit, type = "nLL")[[1]]$data, get_nLL_plot(list(fit$data), list(fit$rep), "Model 1")[[1]]$data)

  # the default is the time series plot
  ts_plots <- plot(fit)
  ts_direct <- get_ts_plot(list(fit$rep), list(fit$sd_rep), "Model 1")
  expect_equal(length(ts_plots), length(ts_direct))
  expect_s3_class(ts_plots[[1]], "ggplot")
  expect_equal(ts_plots[[1]]$data, ts_direct[[1]]$data)

  # every other type runs and returns plots
  expect_true(length(plot(fit, type = "biological")) > 0)
  expect_true(length(plot(fit, type = "catch_fits")) > 0)
  expect_null(plot(fit, type = "at_age_fits")) # dusky fits no at-age data source

  # the legend name and the plotting function's own arguments are passed on
  expect_true(all(plot(fit, type = "selex", model_name = "dusky")[[1]]$data$Model == "dusky"))
  expect_equal(plot(fit, type = "selex", year_indx = 1)[[1]]$data,
               get_selex_plot(list(fit$rep), "Model 1", year_indx = 1)[[1]]$data)

  expect_error(plot(fit, type = "selectivity"), "type must be one of")

})

test_that("a fit without the sdreport prints and summarizes without standard errors", {

  fit_no_sd <- fit_model(data, parameters, mapping, newton_loops = 1, silent = TRUE)
  expect_error(plot(fit_no_sd), "needs the sdreport") # the default type reads it
  expect_s3_class(plot(fit_no_sd, type = "idx_fits"), "ggplot")
  printed <- capture.output(print(fit_no_sd))
  expect_true(any(grepl("Hessian:        not checked", printed)))
  expect_true(any(grepl("optimizer:      converged \\(nlminb code 0", printed))) # convergence does not need the sdreport
  expect_true(any(grepl("max gradient:", printed)))
  expect_false(any(grepl("sanity checks", printed)))

  par_table <- summary(fit_no_sd)
  expect_true(all(is.na(par_table$SE_Max)))
  expect_equal(par_table$Abs_Grad_Max, summary(fit)$Abs_Grad_Max) # same fit, so the same gradients
  expect_equal(par_table$Est_Min, summary(fit)$Est_Min) # and the same estimates

})

test_that("the stored gradient is the gradient at the estimates, and the one the sdreport has", {

  expect_equal(fit$optim$gradient, as.numeric(fit$gr(fit$optim$par)))
  expect_equal(fit$optim$gradient, as.numeric(fit$sd_rep$gradient.fixed), tolerance = 1e-8)

})

test_that("print says so when nlminb stopped before converging", {

  # two iterations is not enough to reach the minimum, so nlminb returns a nonzero code
  fit_stopped <- fit_model(data, parameters, mapping, newton_loops = 0, silent = TRUE, nlminb_control = list(iter.max = 2))
  expect_true(fit_stopped$optim$convergence != 0)

  printed <- capture.output(print(fit_stopped))
  expect_true(any(grepl("optimizer:      not converged \\(nlminb code 1", printed)))
  expect_true(max(abs(fit_stopped$optim$gradient)) > 1) # and the gradient shows it

})

test_that("a fit that was not optimized prints, and summary and logLik stop", {

  fit_no_optim <- fit_model(data, parameters, mapping, do_optim = FALSE, silent = TRUE)
  expect_s3_class(fit_no_optim, "sporc_fit")
  printed <- capture.output(print(fit_no_optim))
  expect_true(any(grepl("optimizer:      not run", printed)))
  expect_true(any(grepl("at the starting values", printed)))

  expect_error(summary(fit_no_optim), "do_optim = FALSE")
  expect_error(logLik(fit_no_optim), "do_optim = FALSE")

})

test_that("with random effects, print names them, summary reports them and logLik is the marginal", {

  data("sgl_rg_dusky_data")
  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  fit_re <- suppressMessages(fit_model(input_list$data, input_list$par, input_list$map, random = c("ln_RecDevs", "ln_InitDevs"),
                                       newton_loops = 0, silent = TRUE, do_sdrep = TRUE))

  printed <- capture.output(print(fit_re))
  expect_true(any(grepl("random effects \\(ln_RecDevs, ln_InitDevs\\)", printed)))
  expect_true(any(grepl("marginal nLL:", printed)))
  expect_false(any(grepl("jnLL:", printed)))

  expect_equal(as.numeric(logLik(fit_re)), -fit_re$rep$marginal_nLL, tolerance = 1e-6)
  expect_equal(attr(logLik(fit_re), "df"), length(fit_re$par)) # fixed effects only

  # random effect rows against RTMB's own summary of the sdreport
  par_table <- summary(fit_re)
  random_rows <- par_table[par_table$Type == "random", ]
  expect_setequal(random_rows$Par, c("ln_RecDevs", "ln_InitDevs"))
  expect_true(all(is.na(random_rows$Abs_Grad_Max)))

  sd_random <- summary(fit_re$sd_rep, "random")
  for(i in seq_len(nrow(random_rows))) {
    par_idx <- which(rownames(sd_random) == random_rows$Par[i])
    expect_equal(random_rows$N[i], length(par_idx))
    expect_equal(random_rows$Est_Max[i], max(sd_random[par_idx, "Estimate"]))
    expect_equal(random_rows$SE_Max[i], max(sd_random[par_idx, "Std. Error"]))
  } # end i loop

})
