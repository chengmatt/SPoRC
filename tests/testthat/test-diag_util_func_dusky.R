library(SPoRC)
library(testthat)
data("sgl_rg_dusky_data")
data("dusky_rtmb_model")

test_that("Diagnostic and utility functions run on fitted dusky model", {

  data <- dusky_rtmb_model$data
  parameters <- dusky_rtmb_model$parameters
  mapping <- dusky_rtmb_model$mapping

  # Post optimization checks
  test_that("post_optim_sanity_checks runs", {
    optim_checks <- post_optim_sanity_checks(dusky_rtmb_model$sdrep, dusky_rtmb_model$rep)
    expect_type(optim_checks, "logical")
  })

  # Test if parameter estimation info can be pulled out
  test_that("get_par_est_info runs", {
    par_info <- get_par_est_info(parameters, mapping, dusky_rtmb_model$sdrep)
    expect_type(par_info, "list")
  })

  # Test marginal AIC
  test_that("marg_AIC runs", {
    aic_val <- marg_AIC(dusky_rtmb_model$optim)
    expect_type(aic_val, "double")
  })

  # Test if get_optim_param_list can be pulled out
  test_that("get_optim_param_list runs", {
    optim_param_list <- get_optim_param_list(parameters, mapping, dusky_rtmb_model$sdrep, NULL)
    expect_type(optim_param_list, "list")
  })

  # Model fits testing
  test_that("get_idx_fits_plot, theme_sablefish, do_runs_test, get_idx_fits, and get_comp_prop runs", {
    idx_fits <- get_idx_fits_plot(list(data), list(dusky_rtmb_model$rep), 'Dusky Model') + theme_sablefish()
    comp_fits <- get_comp_prop(data, dusky_rtmb_model$rep, age_labels = 1:27, len_labels = data$lens, year_labels = data$years)
    idx_fits_df <- get_idx_fits(data, dusky_rtmb_model$rep, data$years)
    resids <- do_runs_test(idx_fits_df$resid)
    expect_type(idx_fits, "object")
    expect_type(comp_fits, "list")
    expect_type(idx_fits_df, "list")
    expect_type(resids, "list")
  })

  # Plotting basic plots
  test_that("plot_all_basic runs", {
    expect_no_error(
      suppressWarnings(
        plot_all_basic(list(data), list(dusky_rtmb_model$rep),
                       list(dusky_rtmb_model$sdrep),
                       'Dusky', out_path = NULL)
      )
    )
  })

  # Test whether retro function runs
  test_that("do_retrospective, get_retrospective_plot, and truncate_yr runs", {
    retro <- do_retrospective(data, parameters, mapping, do_par = FALSE,
                              random = NULL, n_retro = 1, newton_loops = 1, return_models = TRUE)
    retro_par <- do_retrospective(data, parameters, mapping, do_par = FALSE,
                                  random = NULL, n_retro = 1, newton_loops = 1, n_cores = 1, return_models = TRUE)
    retro_plot <- get_retrospective_plot(retro$retro_df, 4)
    trunc_retro <- SPoRC:::truncate_yr(0, data, parameters, mapping)

    expect_type(retro, "list")
    expect_type(retro_par, "list")
    expect_type(retro_plot, "list")
    expect_type(trunc_retro, "list")
  })

  # Test whether reference points runs
  test_that("Get_Reference_Points runs", {
    rp <- Get_Reference_Points(data = data, rep = dusky_rtmb_model$rep,
                               SPR_x = 0.4, what = 'SPR', type = 'single_region')
    expect_type(rp, "list")
  })

  # Test likelihood profile
  test_that("do_likelihood_profile runs", {
    lp <- do_likelihood_profile(data, parameters, mapping, NULL,
                                what = 'ln_global_R0', min_val = 1, max_val = 3, inc = 1)
    lp_par  <- do_likelihood_profile(data, parameters, mapping, NULL,
                                     what = 'ln_global_R0',
                                     min_val = 1, max_val = 3, inc = 1,
                                     do_par = TRUE, n_cores = 3)
    expect_type(lp, "list")
    expect_type(lp_par, "list")

  })

  # Test jitter analysis
  test_that("do_jitter runs", {
    jit_df <- do_jitter(data, parameters, mapping, NULL, sd = 0.1, n_jitter = 1, do_par = FALSE)
    jit_par_df <- do_jitter(data = data, parameters = parameters, mapping = mapping, random = NULL, sd = 0.1, n_jitter = 1, do_par = TRUE, n_cores = 1)
    expect_type(jit_df, "list")
    expect_type(jit_par_df, "list")
  })

  # Test francis
  test_that("run_francis runs", {
    francis_rewt <- run_francis(data, parameters, mapping, NULL, 1, 0)
    expect_type(francis_rewt, "list")
  })

  # Test self tests
  test_that("simulation_self_test runs", {
    suppressWarnings(
      self_test <-   # Test self tests
        simulation_self_test(
          data, parameters, mapping, NULL,
          dusky_rtmb_model$rep, dusky_rtmb_model$sdrep, 1
        )
    )
    expect_type(self_test, "list")
  })

  # Test condition closed loop simulations
  test_that("condition_closed_loop_simulations runs", {
    condition_cl_sim <- condition_closed_loop_simulations(
      1, 1, data, parameters, mapping,
      dusky_rtmb_model$sdrep, dusky_rtmb_model$rep, NULL
    )
    dusky_sim_env <- Setup_sim_env(condition_cl_sim)
    expect_type(condition_cl_sim, "list")
    expect_type(dusky_sim_env, "environment")

  })



})
