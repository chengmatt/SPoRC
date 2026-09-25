# Checks that n_proj_yrs_devs extends growth like recruitment: the deviation arrays gain the projected
# years, curve-mode Get_Growth runs over them, and the report's weight at age follows the dsem forecast (1e-6).

test_that("the fit reports projected weight at age driven by the dsem's forecast", {

  data("sgl_rg_ebs_pcod_data", envir = environment())
  build_ext <- build_ebs_pcod_input
  body(build_ext) <- do.call(substitute, list(body(build_ext), list(Setup_Mod_Dim = quote(function(...) Setup_Mod_Dim(..., n_proj_yrs_devs = 5)))))
  il <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ext(sgl_rg_ebs_pcod_data))), sgl_rg_ebs_pcod_data)
  il$data$growth_tv_type <- 0 # curve mode builds every year up front, which is what projected years need
  n_yrs <- length(il$data$years); n_proj <- 5
  proj <- n_yrs + 1:n_proj

  # the arrays and their maps reach the projected years, and those years are estimated
  expect_equal(dim(il$par$ln_growth_devs)[3], n_yrs + n_proj)
  expect_equal(dim(il$par$ln_growth_semipar_devs)[3], n_yrs + n_proj)
  expect_equal(dim(il$data$map_ln_growth_devs)[3], n_yrs + n_proj)
  map_L1 <- array(il$map$ln_growth_devs, dim = dim(il$par$ln_growth_devs))[1,1,proj,1,1]
  expect_true(all(!is.na(map_L1)))

  set.seed(4)
  env <- data.frame(year = il$data$years, env = cumsum(rnorm(n_yrs, 0.05, 0.2)))
  target <- "growth_Pop_1_Region_1_Par_1_Sex_1"
  arrows <- hold_arrows(c(paste0("env -> ", target, ", 0, b_env"), "env -> env, 1, rho", "env <-> env, 0, sd_env", paste0(target, " <-> ", target, ", 0, sd_g")),
                        c(sd_g = 0.05, b_env = 0.2, rho = 0.9, sd_env = 0.2))
  d <- suppressMessages(Setup_Mod_DSEM(il, arrows, env, dsem_processes = "growth", dsem_mu_spec = "fix"))
  expect_equal(d$data$dsem_n_grid_yrs, n_yrs + n_proj)

  fit <- fit_model(d$data, d$par, d$map, random = NULL, newton_loops = 0, silent = TRUE)
  pars <- fit$env$parList()

  # the report holds the projected years, and their weight at age is finite
  expect_equal(dim(fit$rep$WAA)[3], n_yrs + n_proj)
  expect_true(all(is.finite(fit$rep$WAA[1,1,proj,,,])))

  # in years with nothing else to read, the projected L1 deviation is the arrow applied to the env forecast
  expect_equal(unname(pars$ln_growth_devs[1,1,proj,1,1]), 0.2 * (pars$dsem_x[proj,1] - pars$dsem_mu[1]), tolerance = 1e-6)

  # and weight at age moves with it: rerun Get_Growth by hand at the fitted deviations and compare
  dd <- fit$data
  by_hand <- Get_Growth(ln_growth_pars = pars$ln_growth_pars, growth_A1 = dd$growth_A1, growth_A2 = dd$growth_A2, growth_L0 = dd$growth_L0,
                        growth_len_lower = dd$growth_len_lower, growth_cv_type = dd$growth_cv_type, growth_sd_type = dd$growth_sd_type, growth_dist = dd$growth_dist,
                        growth_plus_group = dd$growth_plus_group, growth_L2_asymptote = dd$growth_L2_asymptote, derive_waa = dd$derive_waa, wt_len_pars = dd$wt_len_pars,
                        ages = dd$ages, seasdur = dd$seasdur, spawn_seas = dd$spawn_seas, t_spawn = dd$t_spawn, n_pop = 1, n_regions = 1, n_yrs = n_yrs + n_proj, n_seas = dd$n_seas,
                        n_sexes = 1, n_fish_fleets = dd$n_fish_fleets, n_srv_fleets = dd$n_srv_fleets, t_fish = dd$t_fish, t_srv = dd$t_srv,
                        ln_growth_devs = pars$ln_growth_devs, growth_tv_model = dd$growth_tv_model, growth_tv_link = dd$growth_tv_link, growth_par_bounds = dd$growth_par_bounds,
                        growth_tv_type = 0, growth_cohort_styr = dd$growth_cohort_styr, years_eval = NULL,
                        ln_growth_semipar_devs = pars$ln_growth_semipar_devs, growth_semipar = dd$growth_semipar)
  expect_equal(as.numeric(fit$rep$WAA[1,1,proj,,,]), as.numeric(by_hand$WAA[1,1,proj,,,]), tolerance = 1e-10)

})
