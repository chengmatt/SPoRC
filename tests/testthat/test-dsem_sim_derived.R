# Checks a drawn dsem series reaches the population when it is not a recruitment deviation: growth rebuilt
# by Get_Growth, and by Get_Growth_Year from the OM's own numbers (1e-10), movement by Get_Movement, the NAA state kept.

data("sgl_rg_ebs_pcod_data")
data("sgl_rg_dusky_data")

# pcod with the length at age 1.5 deviations linked to a covariate, at the assessment's values, on a curve
# unless asked for its own cohort propagation
derived_growth_om <- function(n_sims = 2, closed_loop_yrs = 3, condition_on_fit = FALSE, seed = 12, cohort = FALSE) {
  il <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(sgl_rg_ebs_pcod_data))), sgl_rg_ebs_pcod_data)
  if(!cohort) il$data$growth_tv_type <- 0
  n_yrs <- length(il$data$years)
  set.seed(seed)
  env <- data.frame(year = il$data$years, env = cumsum(rnorm(n_yrs, 0, 0.3)))
  target <- "growth_Pop_1_Region_1_Par_1_Sex_1"
  arrows <- c(paste0("env -> ", target, ", 0, b_env"), "env -> env, 1, rho", "env <-> env, 0, sd_env", paste0(target, " <-> ", target, ", 0, sd_g"))
  d <- suppressMessages(Setup_Mod_DSEM(il, hold_arrows(arrows, c(sd_g = 0.1, b_env = 0.3, rho = 0.8, sd_env = 0.3)), env, dsem_processes = "growth", dsem_mu_spec = "fix"))
  obj <- fit_model(d$data, d$par, d$map, random = NULL, do_optim = FALSE, silent = TRUE)
  pars <- obj$env$parList()
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = closed_loop_yrs, n_sims = n_sims, data = obj$data, parameters = obj$parameters,
                                                mapping = obj$mapping, sd_rep = list(par.fixed = obj$par, par.random = NULL), rep = obj$rep, random = NULL)
  sim_list <- Setup_Sim_DSEM(sim_list, obj$data, pars, rep = obj$rep, condition_on_fit = condition_on_fit)
  list(obj = obj, pars = pars, sim_list = sim_list, sim_env = Setup_sim_env(sim_list), n_yrs = n_yrs)
}

# what the fit forms after growth, for one year: the fishery weight of what its length selectivity takes, and
# fishery and survey selectivity at age read through the keys
pcod_after_growth <- function(dd, rep, pars, growth, y, n_ages) {
  n_extra <- dim(growth$WAA)[3] - dim(rep$fish_sel_l)[2] # years past the fit read the last fitted selectivity
  fish_sel_l <- extend_years(rep$fish_sel_l, n_extra, 2, "last")
  srv_sel_l <- extend_years(rep$srv_sel_l, n_extra, 2, "last")
  WAA_fish <- growth_selected_waa_year(growth$WAA_fish, growth$SizeAgeTrans_fish, fish_sel_l, dd$wt_len_pars, growth_len_mid(dd$growth_len_lower),
                                       dd$fish_waa_selected, y, 1, 1, dd$n_seas, 1)
  list(WAA_fish = WAA_fish[1,1,y,,,1,1],
       fish_sel = as.numeric(fish_sel_l[1,y,,1,1] %*% growth$SizeAgeTrans_fish[1,1,y,1,,,1,1]),
       srv_sel = as.numeric(srv_sel_l[1,y,,1,1] %*% growth$SizeAgeTrans_srv[1,1,y,1,,,1,1]))
}

# Get_Growth spelled out from the data list, the way the objective calls it, at given deviations
pcod_growth_by_hand <- function(dd, pars, ln_growth_devs, ln_growth_semipar_devs, n_yrs, growth_tv_type = 0) {
  Get_Growth(ln_growth_pars = pars$ln_growth_pars, growth_A1 = dd$growth_A1, growth_A2 = dd$growth_A2, growth_L0 = dd$growth_L0,
             growth_len_lower = dd$growth_len_lower, growth_cv_type = dd$growth_cv_type, growth_sd_type = dd$growth_sd_type, growth_dist = dd$growth_dist,
             growth_plus_group = dd$growth_plus_group, growth_L2_asymptote = dd$growth_L2_asymptote, derive_waa = dd$derive_waa, wt_len_pars = dd$wt_len_pars,
             ages = dd$ages, seasdur = dd$seasdur, spawn_seas = dd$spawn_seas, t_spawn = dd$t_spawn, n_pop = 1, n_regions = 1, n_yrs = n_yrs, n_seas = dd$n_seas,
             n_sexes = 1, n_fish_fleets = dd$n_fish_fleets, n_srv_fleets = dd$n_srv_fleets, t_fish = dd$t_fish, t_srv = dd$t_srv,
             ln_growth_devs = ln_growth_devs, growth_tv_model = dd$growth_tv_model, growth_tv_link = dd$growth_tv_link, growth_par_bounds = dd$growth_par_bounds,
             growth_tv_type = growth_tv_type, growth_cohort_styr = dd$growth_cohort_styr, years_eval = NULL,
             ln_growth_semipar_devs = ln_growth_semipar_devs, growth_semipar = dd$growth_semipar)
}

test_that("a drawn growth series rebuilds weight at age and the keys through Get_Growth", {

  om <- derived_growth_om()
  se <- om$sim_env
  n_sim_yrs <- om$n_yrs + 3

  # the operating model now holds the deviation arrays, sized to its own years, with the linked cells marked
  expect_equal(dim(se$ln_growth_devs), c(dim(om$pars$ln_growth_devs)[1:2], n_sim_yrs, dim(om$pars$ln_growth_devs)[4:5], 2))
  expect_equal(dim(se$ln_growth_semipar_devs)[3], n_sim_yrs)
  drawn <- se$dsem_drawn$ln_growth_devs
  expect_true(all(drawn[1,1,,1,1]))
  expect_false(any(drawn[1,1,,-1,1]))

  # the linked series is the drawn grid column, different in each replicate; the other parameters keep the fit
  expect_equal(se$ln_growth_devs[1,1,,1,1,1], se$dsem_x_sim[,2,1], tolerance = 1e-12)
  expect_gt(max(abs(se$ln_growth_devs[1,1,,1,1,1] - se$ln_growth_devs[1,1,,1,1,2])), 0.05)
  for(sim in 1:2) expect_equal(se$ln_growth_devs[1,1,1:om$n_yrs,-1,1,sim], om$pars$ln_growth_devs[1,1,1:om$n_yrs,-1,1], tolerance = 1e-12)

  # weight at age and the size-age keys are Get_Growth at that replicate's deviations, over every simulated year,
  # and what the fit forms after growth follows: the fishery weighs what it selects, selectivity at age reads the key
  for(sim in 1:2) {
    by_hand <- pcod_growth_by_hand(om$obj$data, om$pars, array(se$ln_growth_devs[,,,,,sim], dim = dim(se$ln_growth_devs)[-6]),
                                   array(se$ln_growth_semipar_devs[,,,,,sim], dim = dim(se$ln_growth_semipar_devs)[-6]), n_sim_yrs)
    expect_equal(as.numeric(se$WAA[,,,,,,sim]), as.numeric(by_hand$WAA), tolerance = 1e-10)
    expect_equal(as.numeric(se$WAA_srv[,,,,,,,sim]), as.numeric(by_hand$WAA_srv), tolerance = 1e-10) # the survey does not weigh by selection
    expect_equal(as.numeric(se$SizeAgeTrans_fish[,,,,,,,,sim]), as.numeric(by_hand$SizeAgeTrans_fish), tolerance = 1e-10)
    expect_equal(as.numeric(se$SizeAgeTrans_srv[,,,,,,,,sim]), as.numeric(by_hand$SizeAgeTrans_srv), tolerance = 1e-10)
    for(y in c(1, om$n_yrs, n_sim_yrs)) {
      after <- pcod_after_growth(om$obj$data, om$obj$rep, om$pars, by_hand, y, length(om$obj$data$ages))
      expect_equal(se$WAA_fish[1,1,y,,,1,1,sim], after$WAA_fish, tolerance = 1e-10)
      expect_equal(se$fish_sel[1,1,y,1,,1,1,sim], after$fish_sel, tolerance = 1e-10)
      expect_equal(se$srv_sel[1,1,y,1,,1,1,sim], after$srv_sel, tolerance = 1e-10)
    } # end y loop
  } # end sim loop

  # so the replicates no longer share the report's weight at age or selectivity at age
  expect_gt(max(abs(se$WAA[1,1,,1,,1,1] - se$WAA[1,1,,1,,1,2])), 1e-3)
  expect_gt(max(abs(se$WAA[1,1,1:om$n_yrs,1,,1,1] - om$obj$rep$WAA[1,1,1:om$n_yrs,1,,1])), 1e-3)
  expect_gt(max(abs(se$fish_sel[1,1,,1,,1,1,1] - se$fish_sel[1,1,,1,,1,1,2])), 1e-3)

})

test_that("conditioned on the fit, fitted years reproduce the report's growth and later years move", {

  om <- derived_growth_om(condition_on_fit = TRUE)
  se <- om$sim_env
  fit_yrs <- 1:om$n_yrs

  for(sim in 1:2) {
    expect_equal(se$ln_growth_devs[1,1,fit_yrs,,1,sim], om$pars$ln_growth_devs[1,1,fit_yrs,,1], tolerance = 1e-12)
    expect_equal(as.numeric(se$WAA[1,1,fit_yrs,,,,sim]), as.numeric(om$obj$rep$WAA[1,1,fit_yrs,,,]), tolerance = 1e-10)
    # the report's fishery weight is the selected one and its selectivities at age come through the fit's keys
    expect_equal(as.numeric(se$WAA_fish[1,1,fit_yrs,,,,,sim]), as.numeric(om$obj$rep$WAA_fish[1,1,fit_yrs,,,,]), tolerance = 1e-10)
    expect_equal(as.numeric(se$fish_sel[1,1,fit_yrs,,,,,sim]), as.numeric(om$obj$rep$fish_sel[1,1,fit_yrs,,,,]), tolerance = 1e-10)
    expect_equal(as.numeric(se$srv_sel[1,1,fit_yrs,,,,,sim]), as.numeric(om$obj$rep$srv_sel[1,1,fit_yrs,,,,]), tolerance = 1e-10)
  } # end sim loop
  proj <- om$n_yrs + 1:3
  expect_gt(max(abs(se$WAA[1,1,proj,1,,1,1] - se$WAA[1,1,proj,1,,1,2])), 1e-4)

})

test_that("cohort growth advances year by year from the operating model's own numbers at age", {

  om <- derived_growth_om(cohort = TRUE, n_sims = 2)
  sim_env <- Setup_sim_env(om$sim_list) # built here, since the annual cycle finds the environment by this name in the frame that built it
  dd <- om$obj$data
  styr <- dd$growth_cohort_styr
  n_sim_yrs <- om$n_yrs + 3
  expect_equal(dd$growth_tv_type, 1) # pcod propagates cohort by cohort from styr
  expect_length(sim_env$growth_state, 2)

  # the years before the propagation are built up front, the rest wait for the population
  built <- seq_len(styr - 1)
  g <- pcod_growth_by_hand(dd, om$pars, array(sim_env$ln_growth_devs[,,,,,1], dim = dim(sim_env$ln_growth_devs)[-6]),
                           array(sim_env$ln_growth_semipar_devs[,,,,,1], dim = dim(sim_env$ln_growth_semipar_devs)[-6]), n_sim_yrs, growth_tv_type = 1)
  expect_equal(as.numeric(sim_env$WAA[1,1,built,,,,1]), as.numeric(g$WAA[1,1,built,,,]), tolerance = 1e-10)
  expect_true(all(sim_env$WAA[1,1,styr:n_sim_yrs,,,,1] == om$sim_list$WAA[1,1,styr:n_sim_yrs,,,,1])) # still the report's copy

  set.seed(16)
  for(y in 1:n_sim_yrs) run_annual_cycle(y, 1, sim_env)

  # the same chain by hand, each year's plus group blended by the numbers this replicate had at the start of it
  for(y in styr:n_sim_yrs) {
    g <- Get_Growth_Year(growth = g, y = y, NAA_y = array(sim_env$NAA[,,y,1,,,1], dim = c(1, 1, length(dd$ages), 1)),
                         ln_growth_pars = om$pars$ln_growth_pars, ln_growth_devs = array(sim_env$ln_growth_devs[,,,,,1], dim = dim(sim_env$ln_growth_devs)[-6]),
                         growth_tv_model = dd$growth_tv_model, growth_tv_link = dd$growth_tv_link, growth_par_bounds = dd$growth_par_bounds,
                         growth_A1 = dd$growth_A1, growth_A2 = dd$growth_A2, growth_L0 = dd$growth_L0, growth_len_lower = dd$growth_len_lower,
                         growth_cv_type = dd$growth_cv_type, growth_sd_type = dd$growth_sd_type, growth_dist = dd$growth_dist, growth_plus_group = dd$growth_plus_group,
                         growth_L2_asymptote = dd$growth_L2_asymptote, derive_waa = dd$derive_waa, wt_len_pars = dd$wt_len_pars, ages = dd$ages, seasdur = dd$seasdur,
                         spawn_seas = dd$spawn_seas, t_spawn = dd$t_spawn, n_pop = 1, n_regions = 1, n_seas = dd$n_seas, n_sexes = 1, t_fish = dd$t_fish, t_srv = dd$t_srv,
                         ln_growth_semipar_devs = array(sim_env$ln_growth_semipar_devs[,,,,,1], dim = dim(sim_env$ln_growth_semipar_devs)[-6]), growth_semipar = dd$growth_semipar)
  } # end y loop
  expect_equal(as.numeric(sim_env$WAA[1,1,,,,,1]), as.numeric(g$WAA[1,1,,,,]), tolerance = 1e-10)
  expect_equal(as.numeric(sim_env$SizeAgeTrans_fish[1,1,,,,,,,1]), as.numeric(g$SizeAgeTrans_fish[1,1,,,,,,]), tolerance = 1e-10)
  for(y in c(styr, n_sim_yrs)) {
    after <- pcod_after_growth(dd, om$obj$rep, om$pars, g, y, length(dd$ages))
    expect_equal(sim_env$WAA_fish[1,1,y,,,1,1,1], after$WAA_fish, tolerance = 1e-10)
    expect_equal(sim_env$fish_sel[1,1,y,1,,1,1,1], after$fish_sel, tolerance = 1e-10)
    expect_equal(sim_env$srv_sel[1,1,y,1,,1,1,1], after$srv_sel, tolerance = 1e-10)
  } # end y loop

  # the plus group moved with the population, so the propagated years are not the report's
  expect_gt(max(abs(sim_env$WAA[1,1,styr:om$n_yrs,1,,1,1] - om$obj$rep$WAA[1,1,styr:om$n_yrs,1,,1])), 1e-4)
  expect_true(all(is.finite(sim_env$SSB[,,,1])))

})

test_that("a growth link under selectivity at length needs the report", {

  il <- seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(sgl_rg_ebs_pcod_data))), sgl_rg_ebs_pcod_data)
  n_yrs <- length(il$data$years)
  target <- "growth_Pop_1_Region_1_Par_1_Sex_1"
  d <- suppressMessages(Setup_Mod_DSEM(il, c(paste0("env -> ", target, ", 0, b_env"), "env <-> env, 0, sd_env", paste0(target, " <-> ", target, ", 0, sd_g")),
                                       data.frame(year = il$data$years, env = rnorm(n_yrs)), dsem_processes = "growth", dsem_mu_spec = "fix"))
  expect_error(Setup_Sim_DSEM(list(n_pop = 1, n_regions = 1, n_yrs = n_yrs, n_sims = 2), d$data, d$par), "needs rep")

})

test_that("a drawn movement series rebuilds the movement matrix through Get_Movement", {

  sweep <- sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA)) # three regions, movement estimated
  n_yrs <- length(sweep$data$years)
  series <- "move_Pop_1_From_1_To_1_Seas_1_Age_3_Sex_1" # recruits stay put in this model, so an older age

  # fixed movement never reads its deviations, so a link there is refused
  fixed <- sweep_input()
  fixed$map$move_devs <- factor(seq_along(fixed$par$move_devs))
  fixed$data$map_move_devs <- array(as.numeric(fixed$map$move_devs), dim = dim(fixed$par$move_devs))
  expect_error(suppressMessages(Setup_Mod_DSEM(fixed, c(paste0("x -> ", series, ", 0, b"), "x <-> x, 0, sx", paste0(series, " <-> ", series, ", 0, sm")),
                                               data.frame(year = fixed$data$years, x = rnorm(n_yrs)), dsem_processes = "move", dsem_mu_spec = "fix")), "movement is fixed")

  # the sweep model maps every movement deviation off, so estimate this one series first
  map_move <- array(as.integer(sweep$map$move_devs), dim = dim(sweep$par$move_devs))
  map_move[1,1,1,,1,3,1] <- seq_len(n_yrs)
  sweep$map$move_devs <- factor(map_move)
  sweep$data$map_move_devs <- array(as.numeric(sweep$map$move_devs), dim = dim(sweep$par$move_devs))
  set.seed(13)
  arrows <- c(paste0("x -> ", series, ", 0, b"), "x -> x, 1, r", "x <-> x, 0, sx", paste0(series, " <-> ", series, ", 0, sm"))
  linked <- suppressMessages(Setup_Mod_DSEM(sweep, hold_arrows(arrows, c(b = 0.5, r = 0.5, sx = 1, sm = 0.3)), data.frame(year = sweep$data$years, x = rnorm(n_yrs)),
                                            dsem_processes = "move", dsem_mu_spec = "fix"))
  obj <- fit_model(linked$data, linked$par, linked$map, random = NULL, do_optim = FALSE, silent = TRUE)
  pars <- obj$env$parList()
  dd <- obj$data

  # the smallest operating model list the draw needs: the dims and the report's movement for every replicate
  n_sims <- 3
  sim_list <- list(n_pop = 1, n_regions = 3, n_yrs = n_yrs, n_seas = 1, n_ages = 7, n_sexes = 2, n_sims = n_sims, Movement = replicate(n_sims, obj$rep$Movement))
  sim_list <- Setup_Sim_DSEM(sim_list, dd, pars)
  se <- Setup_sim_env(sim_list)

  expect_equal(dim(se$move_devs), c(dim(pars$move_devs), n_sims))
  expect_true(all(se$dsem_drawn$move_devs[1,1,1,,1,3,1]))
  expect_equal(sum(se$dsem_drawn$move_devs), n_yrs)
  expect_equal(se$move_devs[1,1,1,,1,3,1,2], se$dsem_x_sim[,2,2], tolerance = 1e-12)

  for(sim in 1:n_sims) {
    by_hand <- Get_Movement(move_type = dd$move_type, do_recruits_move = dd$do_recruits_move, n_pop = 1, n_regions = 3, n_yrs = n_yrs, n_proj_yrs_devs = 0,
                            n_ages = 7, n_sexes = 2, n_seas = 1, move_pars = pars$move_pars, move_devs = array(se$move_devs[,,,,,,,sim], dim = dim(pars$move_devs)),
                            use_fixed_movement = dd$use_fixed_movement, Fixed_Movement = dd$Fixed_Movement, log_move_diffusion_pars = pars$log_move_diffusion_pars,
                            move_preference_pars = pars$move_preference_pars, area_r = dd$area_r, adjacency_mat = dd$adjacency_mat,
                            ctmc_diffusion_bounds = dd$ctmc_diffusion_bounds, ctmc_diffusion_eps = dd$ctmc_diffusion_eps, seasdur = dd$seasdur,
                            ctmc_scale_by_seasdur = dd$ctmc_scale_by_seasdur, expm_nsub = 0)
    expect_equal(as.numeric(se$Movement[,,,,,,,sim]), as.numeric(by_hand$Movement), tolerance = 1e-12)
    for(y in 1:n_yrs) expect_equal(unname(rowSums(se$Movement[1,,,y,1,3,1,sim])), rep(1, 3), tolerance = 1e-12) # still rows of probabilities
    expect_equal(se$Movement[1,,,,1,2,1,sim], obj$rep$Movement[1,,,,1,2,1], tolerance = 1e-12) # an unlinked age keeps the fit's movement
  } # end sim loop
  expect_gt(max(abs(se$Movement[1,1,,,1,3,1,1] - se$Movement[1,1,,,1,3,1,2])), 1e-3) # the linked row moves with the draw

})

test_that("a drawn numbers at age series replaces its cells of the replicate's innovations", {

  build_naa <- build_goa_dusky_input
  body(build_naa) <- do.call(substitute, list(body(build_naa), list(Setup_Mod_Biologicals = quote(function(...) Setup_Mod_Biologicals(..., NAA_re = "iid")))))
  il <- suppressMessages(build_naa(sgl_rg_dusky_data))
  n_yrs <- length(il$data$years)
  series <- "NAA_Pop_1_Region_1_Seas_1_Age_3_Sex_1"
  set.seed(14)
  # a state near the deterministic path, jittered so the fitted innovations are not all zero
  base <- build_goa_dusky_input(sgl_rg_dusky_data)
  base_rep <- fit_model(base$data, base$par, base$map, random = NULL, do_optim = FALSE, silent = TRUE)$rep
  il$par$ln_NAA[] <- log(base_rep$NAA[,,1:n_yrs,,,]) + rnorm(length(il$par$ln_NAA), 0, 0.1)
  arrows <- c(paste0("env -> ", series, ", 0, b"), "env -> env, 1, r", "env <-> env, 0, sd_env", paste0(series, " <-> ", series, ", 0, sn"))
  linked <- suppressMessages(Setup_Mod_DSEM(il, hold_arrows(arrows, c(b = 0.4, r = 0.5, sd_env = 1, sn = 0.2)), data.frame(year = il$data$years, env = rnorm(n_yrs)),
                                            dsem_processes = "NAA", dsem_mu_spec = "fix"))
  obj <- fit_model(linked$data, linked$par, linked$map, random = NULL, do_optim = FALSE, silent = TRUE)
  pars <- obj$env$parList()
  sim_list <- condition_closed_loop_simulations(closed_loop_yrs = 2, n_sims = 2, data = obj$data, parameters = obj$parameters, mapping = obj$mapping,
                                                sd_rep = list(par.fixed = obj$par, par.random = NULL), rep = obj$rep, random = NULL)
  expect_gt(sim_list$NAA_re, 0) # the closed loop conditions the state on from the data

  # the state has to be on, and the report is needed to measure the fitted state against its prediction
  expect_error(Setup_Sim_DSEM(sim_list, obj$data, pars), "needs rep")
  off <- sim_list
  off$NAA_re <- 0
  expect_error(Setup_Sim_DSEM(off, obj$data, pars, rep = obj$rep), "NAA_re is off")

  sim_list <- Setup_Sim_DSEM(sim_list, obj$data, pars, rep = obj$rep, condition_on_fit = TRUE)
  sim_env <- Setup_sim_env(sim_list)
  drawn <- sim_env$dsem_drawn$naa_eta_all
  expect_equal(dim(drawn), c(1, 1, n_yrs + 2, 1, 30, 1))
  expect_true(all(drawn[1,1,2:(n_yrs + 2),1,3,1])) # the linked age, from the first state year through the closed loop
  expect_equal(sum(drawn), n_yrs + 1)

  # the fitted years are the fit's own innovations, the state less its prediction
  cells <- obj$data$dsem_link_cell[[1]]
  eta_fit <- pars$ln_NAA[cells] - log(obj$rep$NAA_pred[cells])
  expect_equal(sim_env$dsem_x_sim[2:n_yrs,2,1], eta_fit, tolerance = 1e-12)
  expect_equal(sim_env$naa_eta_all[1,1,2:n_yrs,1,3,1,1], eta_fit, tolerance = 1e-12)

  # the replicate draws its own innovations everywhere else and keeps the dsem's at the linked cells
  set.seed(15)
  run_annual_cycle(1, 1, sim_env)
  expect_equal(sim_env$naa_eta[drawn], array(sim_env$naa_eta_all[,,,,,,1], dim = dim(drawn))[drawn], tolerance = 1e-12)
  expect_equal(sim_env$naa_eta[1,1,2:n_yrs,1,3,1], eta_fit, tolerance = 1e-12)
  expect_gt(stats::sd(sim_env$naa_eta[1,1,2:n_yrs,1,4,1]), 0) # a neighboring age is the replicate's own draw
  expect_equal(as.numeric(sim_env$naa_eta_all[,,,,,,1]), as.numeric(sim_env$naa_eta)) # and the stored copy is what the population reads

})
