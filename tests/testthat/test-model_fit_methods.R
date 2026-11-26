# Checks the sporc_fit class fit_model puts on its return: the fit still reads as a list, print shows the fit rather than
# the list, summary, coef, vcov and confint match the sdreport, logLik gives the same AIC as marg_AIC, as.data.frame
# matches the report, and plot matches the plotting functions.

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

test_that("coef and vcov return the sdreport's fixed effects and their covariance", {

  expect_equal(coef(fit), fit$sd_rep$par.fixed)
  expect_named(coef(fit), names(fit$par)) # names repeat, one per estimated value
  expect_identical(vcov(fit), fit$sd_rep$cov.fixed)
  expect_equal(coef(fit, random = TRUE), coef(fit)) # no random effects in this fit

})

test_that("confint is the estimate plus and minus a normal quantile times its standard error", {

  est_val <- fit$sd_rep$par.fixed
  se_val <- sqrt(diag(fit$sd_rep$cov.fixed))

  # every fixed effect, each row against its own estimate and standard error
  ci_table <- confint(fit)
  expect_equal(dim(ci_table), c(length(est_val), 2))
  expect_equal(colnames(ci_table), c("2.5 %", "97.5 %"))
  expect_equal(rownames(ci_table), names(est_val))
  expect_equal(unname(ci_table[, 1]), unname(est_val - qnorm(0.975) * se_val))
  expect_equal(unname(ci_table[, 2]), unname(est_val + qnorm(0.975) * se_val))

  # a parameter estimated once matches stats::confint.default
  expect_equal(confint(fit, parm = "ln_global_R0"), stats::confint.default(fit, parm = "ln_global_R0"))

  # a name with many values gets one interval per value, where confint.default repeats the first
  init_dev_idx <- which(names(est_val) == "ln_InitDevs")
  expect_true(length(init_dev_idx) > 1)
  expect_equal(confint(fit, parm = "ln_InitDevs"), ci_table[init_dev_idx, ])
  expect_true(length(unique(ci_table[init_dev_idx, 1])) > 1)
  expect_equal(length(unique(stats::confint.default(fit)[init_dev_idx, 1])), 1)

  # positions in coef(fit), and another level
  expect_equal(confint(fit, parm = c(1, 3)), ci_table[c(1, 3), ])
  ci_90 <- confint(fit, parm = 1, level = 0.9)
  expect_equal(colnames(ci_90), c("5 %", "95 %"))
  expect_equal(unname(ci_90[1, ]), unname(est_val[1] + qnorm(c(0.05, 0.95)) * se_val[1]))

  expect_error(confint(fit, parm = "not_a_parameter"), "not estimated as fixed effects")
  expect_error(confint(fit, parm = length(est_val) + 1), "positions in coef")

})

test_that("as.data.frame returns a reported quantity with one row per value and its dims labeled", {

  # against label_dims and melt called directly
  ssb_df <- as.data.frame(fit, what = "SSB")
  expect_named(ssb_df, c("Pop", "Region", "Year", "value"))
  expect_equal(ssb_df, reshape2::melt(label_dims(fit$rep, fit$data)$SSB))
  expect_equal(unique(ssb_df$Year), data$years)

  # a vector gets a column for what it runs over, and a single value is one row
  agg_ssb_df <- as.data.frame(fit, what = "Aggregated_SSB")
  expect_named(agg_ssb_df, c("Year", "value"))
  expect_equal(agg_ssb_df$Year, data$years)
  expect_equal(as.data.frame(fit, what = "jnLL")$value, fit$rep$jnLL)

  # labels are passed to label_dims
  naa_df <- as.data.frame(fit, what = "NAA", labels = list(Sex = "Female"))
  expect_true(all(naa_df$Sex == "Female"))

  # every reported quantity returns its values in the order the report stores them
  for(rep_name in names(fit$rep)) {
    rep_df <- as.data.frame(fit, what = rep_name)
    expect_equal(nrow(rep_df), length(fit$rep[[rep_name]]), info = rep_name)
    expect_equal(rep_df$value, as.vector(fit$rep[[rep_name]]), info = rep_name)
  } # end rep_name loop

  expect_error(as.data.frame(fit), "what must name one quantity")
  expect_error(as.data.frame(fit, what = "not_reported"), "what must name one quantity")

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

  # the estimates and the report do not need the sdreport, the covariance and intervals do
  expect_equal(coef(fit_no_sd), coef(fit))
  expect_equal(as.data.frame(fit_no_sd, what = "SSB"), as.data.frame(fit, what = "SSB"))
  expect_error(vcov(fit_no_sd), "needs the sdreport")
  expect_error(confint(fit_no_sd), "needs the sdreport")

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
  expect_error(coef(fit_no_optim), "do_optim = FALSE")
  expect_error(confint(fit_no_optim), "do_optim = FALSE")
  expect_equal(as.data.frame(fit_no_optim, what = "SSB")$value, as.vector(fit_no_optim$rep$SSB)) # report at the starting values

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

  # coef, vcov and confint are the fixed effects, and coef adds the random effects when asked
  expect_equal(coef(fit_re), fit_re$sd_rep$par.fixed)
  expect_false(any(c("ln_RecDevs", "ln_InitDevs") %in% names(coef(fit_re))))
  expect_equal(dim(vcov(fit_re)), rep(length(fit_re$par), 2))
  expect_equal(nrow(confint(fit_re)), length(fit_re$par))
  expect_error(confint(fit_re, parm = "ln_RecDevs"), "not estimated as fixed effects")

  coef_all <- coef(fit_re, random = TRUE)
  expect_equal(length(coef_all), length(fit_re$par) + nrow(sd_random))
  expect_equal(unname(coef_all[names(coef_all) %in% c("ln_RecDevs", "ln_InitDevs")]), unname(sd_random[, "Estimate"]))
  expect_equal(coef_all[!names(coef_all) %in% c("ln_RecDevs", "ln_InitDevs")], coef(fit_re))

})
