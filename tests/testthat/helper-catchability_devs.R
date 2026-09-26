# An operating and estimation model pair for the catchability recovery tests. Sixty years, one region,
# one survey, flat catchability unless a deviation process or a covariate path is asked for.

q_devs_sim <- function(form = "none", sigma_q = 0.25, rho = 0, n_yrs = 60, seed = 20, q_path = NULL, exact = FALSE) {

  n_ages <- 10
  logi <- function(k, a50) 1 / (1 + exp(-k * ((1:n_ages) - a50)))
  arr7 <- function(v, ny, nf) array(rep(v, each = ny), dim = c(1, 1, ny, 1, n_ages, 1, nf))
  arr6 <- function(v) array(rep(v, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))

  sl <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = 1, n_sexes = 1,
                      n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 1, n_pop = 1)
  sl <- Setup_Sim_Containers(sl)

  # exact: catch and indices measured to a tenth of a percent, compositions from a hundred thousand fish,
  # so an estimate lands on its replicate's own truth and replicates can be compared to each other
  obs_sd <- if(exact) 0.001 else 0.2
  iss <- if(exact) 1e5 else 100
  sl <- Setup_Sim_Fishing(sl, fish_sel_input = replicate(1, arr7(logi(3, 5), n_yrs, 1)),
                          ln_sigmaC = array(log(if(exact) 0.001 else 0.02), dim = c(1, n_yrs, 1, 1)),
                          ObsFishIdx_SE = array(obs_sd, dim = c(1, n_yrs, 1, 1)),
                          ISS_FishAgeComps = array(iss, dim = c(1, n_yrs, 1, 1, 1, 1)))

  if(is.null(q_path)) q_path <- rep(1, n_yrs)
  sl <- Setup_Sim_Survey(sl, srv_sel_input = replicate(1, arr7(logi(1, 3), n_yrs, 1)),
                         srv_q_input = array(q_path, dim = c(1, n_yrs, 1, 1)),
                         ObsSrvIdx_SE = array(obs_sd, dim = c(1, n_yrs, 1, 1)),
                         ISS_SrvAgeComps = array(iss, dim = c(1, n_yrs, 1, 1, 1, 1)))

  waa <- 5 * logi(3, 3)
  sl <- suppressWarnings(Setup_Sim_Biologicals(sl, natmort_input = replicate(1, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
                              WAA_input = replicate(1, arr6(waa)), WAA_fish_input = replicate(1, arr7(waa, n_yrs, 1)),
                              WAA_srv_input = replicate(1, arr7(waa, n_yrs, 1)), MatAA_input = replicate(1, arr6(logi(3, 3)))))
  sl <- Setup_Sim_Tagging(sl, use_conv_fish_tagging = 0)
  sl$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))
  sl <- suppressWarnings(Setup_Sim_Rec(sl, R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
                                       ln_sigmaR = array(log(1), dim = c(2, 1, 1)), recruitment_opt = "mean_rec", init_age_strc = 1))

  if(form != "none") sl <- Setup_Sim_q_devs(sl, srv_q_model = form, sigma_srv_q = sigma_q, srv_q_rho = rho)

  set.seed(seed)
  suppressWarnings(suppressMessages(Simulate_Pop_Static(sim_list = sl, output_path = NULL)))
}

# the blocks estimation model with a catchability deviation process in place of blocks
q_devs_em <- function(sim_obj, form, sigma_start = 0.25) {
  il <- blocks_em(sim_obj)
  il <- suppressMessages(setup_q_devs(il, q_model = form, sigma_q_spec = "est_all", q_rho_spec = "est_all",
                                      q_rw_init_sigma = NA, q_type = "est", prefix = "srv",
                                      fleet_field = "n_srv_fleets", use_field = "UseSrvIdx",
                                      fleet_label = "survey fleet", starting_values = list()))
  il$par$ln_sigma_srv_q[] <- log(sigma_start)
  il
}

# parList() takes no argument on purpose: handed last.par.best it mis-slices a fixed plus random vector
q_devs_fit <- function(il, random = "ln_srv_q_devs") {
  fit <- suppressWarnings(suppressMessages(fit_model(il$data, il$par, il$map, random = random,
                                                     do_optim = TRUE, newton_loops = 1, silent = TRUE)))
  list(fit = fit, pars = fit$env$parList(), grad = max(abs(fit$gr(fit$optim$par))))
}

# one dsem driven survey catchability, fitted once and reused: 40 years, one region, one survey, the
# covariate effect and the process error both carried by the catchability series
q_dsem_cache <- new.env(parent = emptyenv())

q_dsem_fit <- function() {

  if(!is.null(q_dsem_cache$fit)) return(q_dsem_cache$fit)

  n_yrs <- 40
  set.seed(77); env_x <- as.numeric(scale(stats::rnorm(n_yrs)))
  set.seed(31); q_path <- exp(0.5 * env_x + stats::rnorm(n_yrs, 0, 0.15))
  om <- q_devs_sim("none", n_yrs = n_yrs, seed = 20, q_path = q_path)

  il <- q_devs_em(om, "dsem")
  tgt <- dsem_series(il, "srv_q")
  arrows <- c(paste0("env -> ", tgt, ", 0, b_env"), "env <-> env, 0, NA, 1", paste0(tgt, " <-> ", tgt, ", 0, sd_q"))
  d <- suppressMessages(Setup_Mod_DSEM(il, arrows, data.frame(year = 1:n_yrs, env = env_x), dsem_processes = "srv_q",
                                       dsem_family = c(env = "fixed"), dsem_mu_spec = c(env = 0)))

  q_dsem_cache$fit <- suppressWarnings(suppressMessages(fit_model(d$data, d$par, d$map, random = "ln_srv_q_devs",
                                                                  do_optim = TRUE, newton_loops = 1, silent = TRUE)))
  q_dsem_cache$fit
}

q_dsem_sdrep <- function(fit) list(par.fixed = fit$optim$par, par.random = fit$env$last.par.best[fit$env$random])

# the simulation list the self test hands the operating model, taken without paying for the refits
q_selftest_simlist <- function(fit) {
  out <- NULL
  testthat::with_mocked_bindings(
    Simulate_Pop_Static = function(sim_list, ...) { out <<- sim_list; stop("captured") },
    try(suppressWarnings(suppressMessages(
      simulation_self_test(data = fit$data, parameters = fit$env$parList(), mapping = fit$mapping,
                           random = "ln_srv_q_devs", rep = fit$rep, sd_rep = q_dsem_sdrep(fit),
                           n_sims = 1, newton_loops = 0, what = "SSB"))), silent = TRUE),
    .package = "SPoRC")
  out
}

# a bare operating model list holding only what the catchability draw reads
q_cond_sl <- function(n_yrs = 20, n_fleets = 1, devs = NULL, n_cond = NULL) {
  sl <- list(n_regions = 1, n_yrs = n_yrs, n_srv_fleets = n_fleets, n_fish_fleets = 1,
             srv_q = array(0.05, dim = c(1, n_yrs, n_fleets, 1)),
             fish_q = array(0.01, dim = c(1, n_yrs, 1, 1)))
  if(!is.null(devs)) sl$ln_srv_q_devs <- array(devs, dim = c(1, n_yrs, n_fleets, 1))
  if(!is.null(n_cond)) sl$n_cond_yrs <- n_cond
  sl
}

# the same survey catchability driven by a derived series, its sd fixed at zero, so the objective
# works it out from the covariate and nothing is integrated
q_dsem_derived_cache <- new.env(parent = emptyenv())

q_dsem_projected_fit <- function() {

  if(!is.null(q_dsem_derived_cache$fit)) return(q_dsem_derived_cache$fit)

  n_yrs <- 40
  set.seed(77); env_x <- as.numeric(scale(stats::rnorm(n_yrs)))
  om <- q_devs_sim("none", n_yrs = n_yrs, seed = 20, q_path = exp(0.5 * env_x))

  il <- q_devs_em(om, "dsem")
  tgt <- dsem_series(il, "srv_q")
  arrows <- c(paste0("env -> ", tgt, ", 0, b_env"), "env <-> env, 0, NA, 1", paste0(tgt, " <-> ", tgt, ", 0, NA, 0"))
  d <- suppressMessages(Setup_Mod_DSEM(il, arrows, data.frame(year = 1:n_yrs, env = env_x), dsem_processes = "srv_q",
                                       dsem_family = c(env = "fixed"), dsem_mu_spec = c(env = 0)))

  q_dsem_derived_cache$fit <- suppressWarnings(suppressMessages(fit_model(d$data, d$par, d$map, random = NULL,
                                                                          do_optim = TRUE, newton_loops = 1, silent = TRUE)))
  q_dsem_derived_cache$fit
}
