# Catchability deviations: the penalty forms against hand written densities, the two refusals, the
# sharing maps, recovery from simulated data, and the dsem routes. Equivalence checked to 1e-8.

data("sgl_rg_ebs_pcod_data")

pcod_q_input <- function() {
  seed_ebs_pcod_mle(suppressWarnings(suppressMessages(build_ebs_pcod_input(sgl_rg_ebs_pcod_data))), sgl_rg_ebs_pcod_data)
}

add_q_devs <- function(il, q_model, q_type = "est", sigma = 0.2) {
  out <- suppressMessages(setup_q_devs(il, q_model = q_model, sigma_q_spec = "est_all", q_rho_spec = "est_all",
                                       q_rw_init_sigma = NA, q_type = q_type, prefix = "srv",
                                       fleet_field = "n_srv_fleets", use_field = "UseSrvIdx",
                                       fleet_label = "survey fleet", starting_values = list()))
  out$par$ln_sigma_srv_q[] <- log(sigma)
  out
}

test_that("catchability is the block value shifted by its deviation", {

  il <- add_q_devs(pcod_q_input(), "iid")
  set.seed(4)
  il$par$ln_srv_q_devs[] <- stats::rnorm(length(il$par$ln_srv_q_devs), 0, 0.2)

  obj <- fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expected <- exp(obj$parameters$ln_srv_q[1,1,1] + as.numeric(il$par$ln_srv_q_devs[1,,1]))

  expect_equal(as.numeric(obj$rep$srv_q[1,,1]), expected, tolerance = 1e-8)
})

test_that("the iid penalty is the normal density of the deviations", {

  il <- add_q_devs(pcod_q_input(), "iid")
  set.seed(5)
  il$par$ln_srv_q_devs[] <- stats::rnorm(length(il$par$ln_srv_q_devs), 0, 0.2)

  pen <- Get_q_dev_penalty(ln_q_devs = il$par$ln_srv_q_devs, ln_sigma_q = il$par$ln_sigma_srv_q,
                           q_rho = il$par$srv_q_rho, q_model = il$data$srv_q_model,
                           map_ln_q_devs = il$data$map_ln_srv_q_devs, q_rw_init_sigma = NA)

  expect_equal(pen, -sum(stats::dnorm(as.numeric(il$par$ln_srv_q_devs), 0, 0.2, TRUE)), tolerance = 1e-10)
})

test_that("a random walk penalizes each step and an ar1 at zero correlation matches iid", {

  devs <- c(0.1, -0.05, 0.2, 0.0)
  il <- add_q_devs(pcod_q_input(), "rw")
  il$par$ln_srv_q_devs <- array(0, dim = c(1, length(devs), 1))
  il$par$ln_srv_q_devs[1,,1] <- devs
  map_all <- array(seq_along(devs), dim = dim(il$par$ln_srv_q_devs))

  rw <- Get_q_dev_penalty(il$par$ln_srv_q_devs, il$par$ln_sigma_srv_q, il$par$srv_q_rho,
                          q_model = 3, map_ln_q_devs = map_all, q_rw_init_sigma = NA)
  by_hand <- -stats::dnorm(devs[1], 0, 0.2, TRUE) - sum(stats::dnorm(devs[-1], devs[-length(devs)], 0.2, TRUE))
  expect_equal(rw, by_hand, tolerance = 1e-10)

  # an ar1 with no correlation is independent deviations, with year one at the stationary spread
  ar1 <- Get_q_dev_penalty(il$par$ln_srv_q_devs, il$par$ln_sigma_srv_q, il$par$srv_q_rho,
                           q_model = 4, map_ln_q_devs = map_all, q_rw_init_sigma = NA)
  iid <- Get_q_dev_penalty(il$par$ln_srv_q_devs, il$par$ln_sigma_srv_q, il$par$srv_q_rho,
                           q_model = 2, map_ln_q_devs = map_all, q_rw_init_sigma = NA)
  expect_equal(ar1, iid, tolerance = 1e-10)
})

test_that("every form builds a model that evaluates and differentiates", {

  for(form in c("iid", "rw", "ar1")) {
    il <- add_q_devs(pcod_q_input(), form)
    obj <- fit_model(il$data, il$par, il$map, random = "ln_srv_q_devs", do_optim = FALSE, silent = TRUE)
    expect_true(is.finite(obj$fn(obj$par)))
    expect_true(all(is.finite(obj$gr(obj$par))))
    expect_equal(sum(!is.na(as.integer(il$map$ln_srv_q_devs))), length(il$data$years))
  } # end form loop
})

test_that("an analytically solved catchability refuses deviations", {

  for(type in c("arith", "geo")) {
    expect_error(add_q_devs(pcod_q_input(), "rw", q_type = type), "solves catchability from the observations")
  } # end type loop
  expect_silent(add_q_devs(pcod_q_input(), "none", q_type = "geo")) # no deviations, so nothing to absorb
})

test_that("catchability blocks refuse deviations", {

  il <- pcod_q_input()
  il$data$srv_q_blocks[,1:20,1] <- 2L
  expect_error(add_q_devs(il, "rw"), "Both describe time variation in catchability")
  expect_silent(add_q_devs(il, "none"))
})

test_that("a derived dsem series drives catchability as a fixed effect regression", {

  il <- add_q_devs(pcod_q_input(), "dsem")
  set.seed(11)
  yrs <- il$data$years
  env <- data.frame(year = yrs, env = stats::rnorm(length(yrs), 2, 1))
  beta <- 0.4
  tgt <- dsem_series(il, "srv_q")

  arrows <- c(paste0("env -> ", tgt, ", 0, NA, ", beta), "env <-> env, 0, sd_env", paste0(tgt, " <-> ", tgt, ", 0, NA, 0"))
  d <- suppressMessages(Setup_Mod_DSEM(il, arrows, env, dsem_processes = "srv_q",
                                       dsem_family = c(env = "fixed"), dsem_mu_spec = c(env = 0)))

  # nothing is integrated: the deviations come from the arrows and the covariate cells are known
  expect_true(d$data$dsem_model$derived[match(tgt, d$data$dsem_model$variables)])
  expect_equal(sum(!is.na(as.integer(d$map$ln_srv_q_devs))), 0)
  expect_equal(sum(!is.na(as.integer(d$map$dsem_x))), 0)

  obj <- fit_model(d$data, d$par, d$map, random = NULL, do_optim = FALSE, silent = TRUE)

  # the multiplicative covariate effect the catchability formula used to apply
  expected <- exp(obj$parameters$ln_srv_q[1,1,1]) * exp(beta * env$env)
  expect_equal(as.numeric(obj$rep$srv_q[1,,1]), as.numeric(expected), tolerance = 1e-8)
})

test_that("a derived series the objective cannot fill is refused rather than left at zero", {

  il <- pcod_q_input()
  il$data$growth_tv_type <- 0
  env <- data.frame(year = il$data$years, env = stats::rnorm(length(il$data$years), 2, 1))
  tgt <- "growth_Pop_1_Region_1_Par_1_Sex_1"
  arrows <- c(paste0("env -> ", tgt, ", 0, b_env"), "env <-> env, 0, sd_env", paste0(tgt, " <-> ", tgt, ", 0, NA, 0"))

  expect_error(suppressMessages(Setup_Mod_DSEM(il, arrows, env, dsem_processes = "growth",
                                               dsem_family = c(env = "fixed"))),
               "the deviations would stay at zero")
})

test_that("the operating model draws a series per replicate and scales catchability by it", {

  n_yrs <- 30
  sl <- list(n_regions = 1, n_yrs = n_yrs, n_srv_fleets = 1, n_fish_fleets = 1,
             srv_q = array(0.05, dim = c(1, n_yrs, 1, 1)), fish_q = array(0.01, dim = c(1, n_yrs, 1, 1)))
  base_q <- sl$srv_q

  sl_rw <- Setup_Sim_q_devs(sl, srv_q_model = "rw", sigma_srv_q = 0.25)
  env <- Setup_sim_env(sl_rw)
  expect_equal(dim(env$ln_srv_q_devs), dim(base_q))

  set.seed(9)
  draw_sim_q_devs(1, env)
  devs <- env$ln_srv_q_devs[1,,1,1]

  expect_equal(as.numeric(env$srv_q[1,,1,1]), as.numeric(base_q[1,,1,1]) * exp(devs), tolerance = 1e-12)
  expect_equal(stats::sd(diff(devs)), 0.25, tolerance = 0.15) # a walk, so the steps hold the sigma

  # no deviations leaves the supplied catchability alone
  env_none <- Setup_sim_env(Setup_Sim_q_devs(sl, srv_q_model = "none"))
  draw_sim_q_devs(1, env_none)
  expect_equal(env_none$srv_q, base_q)
})

test_that("an ar1 in the operating model starts at its stationary spread", {

  n_yrs <- 400
  sl <- list(n_regions = 1, n_yrs = n_yrs, n_srv_fleets = 1, n_fish_fleets = 1,
             srv_q = array(0.05, dim = c(1, n_yrs, 1, 1)), fish_q = array(0.01, dim = c(1, n_yrs, 1, 1)))
  env <- Setup_sim_env(Setup_Sim_q_devs(sl, srv_q_model = "ar1", sigma_srv_q = 0.3, srv_q_rho = 0.7))

  set.seed(12)
  draw_sim_q_devs(1, env)
  devs <- env$ln_srv_q_devs[1,,1,1]

  expect_equal(stats::sd(devs), 0.3 / sqrt(1 - 0.7^2), tolerance = 0.15)
  expect_equal(stats::cor(devs[-1], devs[-n_yrs]), 0.7, tolerance = 0.15)
})

test_that("a retrospective peel shortens the deviation array and its mirror", {

  il <- add_q_devs(pcod_q_input(), "rw")
  n_dev_yrs <- dim(il$par$ln_srv_q_devs)[2]

  peeled <- truncate_yr(j = 2, data = il$data, parameters = il$par, mapping = il$map)

  expect_equal(dim(peeled$retro_parameters$ln_srv_q_devs)[2], n_dev_yrs - 2)
  expect_equal(dim(peeled$retro_data$map_ln_srv_q_devs)[2], n_dev_yrs - 2)
  expect_equal(length(peeled$retro_mapping$ln_srv_q_devs), prod(dim(peeled$retro_parameters$ln_srv_q_devs)))

  # the years kept are the early ones, unchanged
  expect_equal(peeled$retro_parameters$ln_srv_q_devs[1,,1], il$par$ln_srv_q_devs[1, 1:(n_dev_yrs - 2), 1])
})

test_that("independent catchability deviations are recovered from simulated data", {

  om <- q_devs_sim("iid", sigma_q = 0.25, seed = 20)
  truth <- as.numeric(om$ln_srv_q_devs[1,,1,1])

  out <- q_devs_fit(q_devs_em(om, "iid"))
  est <- as.numeric(out$pars$ln_srv_q_devs[1,,1])

  expect_lt(out$grad, 1e-4)
  expect_equal(exp(as.numeric(out$pars$ln_sigma_srv_q)), 0.25, tolerance = 0.35)
  expect_gt(stats::cor(est, truth), 0.7)
})

test_that("a catchability random walk is recovered from simulated data", {

  om <- q_devs_sim("rw", sigma_q = 0.25, seed = 20)
  truth <- as.numeric(om$ln_srv_q_devs[1,,1,1])

  out <- q_devs_fit(q_devs_em(om, "rw"))
  est <- as.numeric(out$pars$ln_srv_q_devs[1,,1])

  expect_lt(out$grad, 1e-4)
  expect_equal(exp(as.numeric(out$pars$ln_sigma_srv_q)), 0.25, tolerance = 0.35) # the step sd, not the spread of the series
  expect_gt(stats::cor(est, truth), 0.8)
})

test_that("an autoregressive catchability process recovers its correlation as well as its sigma", {

  om <- q_devs_sim("ar1", sigma_q = 0.25, rho = 0.6, seed = 20)
  truth <- as.numeric(om$ln_srv_q_devs[1,,1,1])

  out <- q_devs_fit(q_devs_em(om, "ar1"))
  est <- as.numeric(out$pars$ln_srv_q_devs[1,,1])
  rho_hat <- 2 / (1 + exp(-2 * as.numeric(out$pars$srv_q_rho))) - 1

  expect_lt(out$grad, 1e-4)
  expect_equal(exp(as.numeric(out$pars$ln_sigma_srv_q)), 0.25, tolerance = 0.35)
  expect_equal(rho_hat, 0.6, tolerance = 0.35)
  expect_gt(stats::cor(est, truth), 0.7)
})

test_that("a covariate effect on catchability is recovered as a fixed effect", {

  n_yrs <- 60
  beta_true <- 0.5
  set.seed(77)
  env_x <- as.numeric(scale(stats::rnorm(n_yrs))) # centered, so catchability averages near one
  om <- q_devs_sim("none", n_yrs = n_yrs, seed = 20, q_path = exp(beta_true * env_x))

  il <- q_devs_em(om, "dsem")
  tgt <- dsem_series(il, "srv_q")
  arrows <- c(paste0("env -> ", tgt, ", 0, b_env"), "env <-> env, 0, NA, 1", paste0(tgt, " <-> ", tgt, ", 0, NA, 0"))
  d <- suppressMessages(Setup_Mod_DSEM(il, arrows, data.frame(year = 1:n_yrs, env = env_x), dsem_processes = "srv_q",
                                       dsem_family = c(env = "fixed"), dsem_mu_spec = c(env = 0)))

  # the covariate effect costs one parameter and integrates nothing
  expect_equal(length(d$par$dsem_beta), 1)
  expect_equal(sum(!is.na(as.integer(d$map$ln_srv_q_devs))), 0)
  expect_equal(sum(!is.na(as.integer(d$map$dsem_x))), 0)

  out <- q_devs_fit(d, random = NULL)

  expect_lt(out$grad, 1e-4)
  expect_equal(as.numeric(out$pars$dsem_beta), beta_true, tolerance = 0.2)

  # and the catchability path it implies tracks the one the operating model used
  shape_fit <- as.numeric(out$fit$rep$srv_q[1,,1]) / exp(as.numeric(out$pars$ln_srv_q[1,1,1]))
  expect_gt(stats::cor(shape_fit, exp(beta_true * env_x)), 0.99)
})

test_that("the sigma and correlation sharing strings share what they say", {

  # a two region, two survey fleet shell holding only what the mapping reads
  n_r <- 2; n_f <- 2; n_y <- 5
  il <- list(
    data = list(n_regions = n_r, n_srv_fleets = n_f, n_proj_yrs_devs = 0, years = 1:n_y,
                UseSrvIdx = array(1, dim = c(n_r, n_y, 1, n_f)),
                UseSrvIdx_pop = array(0, dim = c(1, n_r, n_y, 1, n_f))),
    par = list(ln_srv_q_devs = array(0, dim = c(n_r, n_y, n_f)),
               ln_sigma_srv_q = array(0, dim = c(n_r, n_f)),
               srv_q_rho = array(0, dim = c(n_r, n_f))),
    map = list()
  )
  n_levels <- function(x) length(unique(stats::na.omit(as.integer(x))))
  map_of <- function(spec, q_model = c(2, 2)) {
    out <- do_q_devs_mapping(il, q_model = q_model, sigma_q_spec = spec, q_rho_spec = spec,
                             prefix = "srv", fleet_field = "n_srv_fleets", use_field = "UseSrvIdx")
    as.integer(out$map$ln_sigma_srv_q)
  }

  expect_equal(n_levels(map_of("est_all")), n_r * n_f)      # one per region and fleet
  expect_equal(n_levels(map_of("est_shared_r")), n_f)       # shared across regions
  expect_equal(n_levels(map_of("est_shared_f")), n_r)       # shared across fleets
  expect_equal(n_levels(map_of("est_shared_r_f")), 1)
  expect_equal(n_levels(map_of("fix")), 0)

  # a fleet sharing across fleets takes the same parameter as its partner
  shared_f <- matrix(map_of("est_shared_f"), n_r, n_f)
  expect_equal(shared_f[,1], shared_f[,2])

  # a fleet with no deviations reads no sigma, and the partner keeps one
  mixed <- matrix(map_of("est_all", q_model = c(2, 1)), n_r, n_f)
  expect_true(all(is.na(mixed[,2])))
  expect_true(all(!is.na(mixed[,1])))
  expect_equal(sort(stats::na.omit(as.vector(mixed))), 1:n_r) # renumbered with no gaps

  # only an ar1 reads a correlation
  rho_map <- function(q_model) as.integer(do_q_devs_mapping(il, q_model = q_model, sigma_q_spec = "est_all",
                                                            q_rho_spec = "est_all", prefix = "srv",
                                                            fleet_field = "n_srv_fleets", use_field = "UseSrvIdx")$map$srv_q_rho)
  expect_equal(n_levels(rho_map(c(2, 3))), 0)  # iid and a random walk
  expect_equal(n_levels(rho_map(c(4, 4))), n_r * n_f)
  expect_equal(n_levels(rho_map(c(4, 2))), n_r) # one ar1 fleet only

  # a dsem fleet hands its sd to the arrows, so it reads no sigma of its own
  expect_equal(n_levels(map_of("est_all", q_model = c(5, 5))), 0)
})

test_that("a region with no index data reads no catchability sigma or correlation", {

  # two regions and one ar1 survey, with the index in region one only
  n_r <- 2; n_y <- 5
  use <- array(1, dim = c(n_r, n_y, 1, 1))
  use[2,,,] <- 0

  il <- list(
    data = list(n_regions = n_r, n_srv_fleets = 1, n_proj_yrs_devs = 0, years = 1:n_y,
                UseSrvIdx = use, UseSrvIdx_pop = array(0, dim = c(1, n_r, n_y, 1, 1))),
    par = list(ln_srv_q_devs = array(0, dim = c(n_r, n_y, 1)),
               ln_sigma_srv_q = array(0, dim = c(n_r, 1)),
               srv_q_rho = array(0, dim = c(n_r, 1))),
    map = list()
  )

  map_with <- function(spec) do_q_devs_mapping(il, q_model = 4, sigma_q_spec = spec, q_rho_spec = spec,
                                               prefix = "srv", fleet_field = "n_srv_fleets", use_field = "UseSrvIdx")
  out <- map_with("est_all")

  dev_map <- matrix(as.integer(out$map$ln_srv_q_devs), n_r, n_y)
  expect_true(all(!is.na(dev_map[1,])))
  expect_true(all(is.na(dev_map[2,])))

  # the region the penalty skips reads neither, and the region with data keeps both
  expect_equal(as.integer(out$map$ln_sigma_srv_q), c(1L, NA))
  expect_equal(as.integer(out$map$srv_q_rho), c(1L, NA))

  # sharing across regions keeps the level alive on the region that does hold data
  expect_equal(as.integer(map_with("est_shared_r")$map$ln_sigma_srv_q), c(1L, NA))

  # what the blanking prevents: a sigma and a correlation no penalty ever reads
  set.seed(4)
  devs <- array(stats::rnorm(n_r * n_y, 0, 0.2), dim = c(n_r, n_y, 1))
  sig <- array(log(0.25), dim = c(n_r, 1))
  rho <- array(0.3, dim = c(n_r, 1))
  pen <- function(s, p) Get_q_dev_penalty(devs, s, p, q_model = 4, map_ln_q_devs = out$data$map_ln_srv_q_devs)
  sig_r2 <- sig; sig_r2[2,1] <- log(9)
  rho_r2 <- rho; rho_r2[2,1] <- 0.95
  expect_equal(pen(sig_r2, rho), pen(sig, rho))
  expect_equal(pen(sig, rho_r2), pen(sig, rho))
})

test_that("a dsem catchability series can hold a covariate effect and process error together", {

  n_yrs <- 60
  beta_true <- 0.5
  set.seed(77); env_x <- as.numeric(scale(stats::rnorm(n_yrs)))
  set.seed(31); pe <- stats::rnorm(n_yrs, 0, 0.15)
  om <- q_devs_sim("none", n_yrs = n_yrs, seed = 20, q_path = exp(beta_true * env_x + pe))

  il <- q_devs_em(om, "dsem")
  tgt <- dsem_series(il, "srv_q")
  arrows <- c(paste0("env -> ", tgt, ", 0, b_env"), "env <-> env, 0, NA, 1", paste0(tgt, " <-> ", tgt, ", 0, sd_q"))
  d <- suppressMessages(Setup_Mod_DSEM(il, arrows, data.frame(year = 1:n_yrs, env = env_x), dsem_processes = "srv_q",
                                       dsem_family = c(env = "fixed"), dsem_mu_spec = c(env = 0)))

  # with an sd line the series is no longer derived, so the deviations are integrated
  expect_false(any(d$data$dsem_model$derived))
  expect_equal(sum(!is.na(as.integer(d$map$ln_srv_q_devs))), n_yrs)
  expect_equal(length(d$par$ln_dsem_sd), 1)

  # the declaration hands the density to the arrows, so the module's own sigma and penalty go quiet
  expect_equal(sum(!is.na(as.integer(d$map$ln_sigma_srv_q))), 0)

  out <- q_devs_fit(d, random = "ln_srv_q_devs")
  expect_equal(out$fit$rep$srv_q_nLL, 0)
  expect_equal(as.numeric(out$pars$dsem_beta), beta_true, tolerance = 0.2)
})

test_that("a vanishing process error matches the fixed effect regression it becomes", {

  n_yrs <- 60
  set.seed(77); env_x <- as.numeric(scale(stats::rnorm(n_yrs)))
  set.seed(31); om <- q_devs_sim("none", n_yrs = n_yrs, seed = 20,
                                 q_path = exp(0.5 * env_x + stats::rnorm(n_yrs, 0, 0.15)))

  fit_with <- function(sd_line, random) {
    il <- q_devs_em(om, "dsem")
    tgt <- dsem_series(il, "srv_q")
    arrows <- c(paste0("env -> ", tgt, ", 0, b_env"), "env <-> env, 0, NA, 1", paste0(tgt, " <-> ", tgt, sd_line))
    d <- suppressMessages(Setup_Mod_DSEM(il, arrows, data.frame(year = 1:n_yrs, env = env_x), dsem_processes = "srv_q",
                                         dsem_family = c(env = "fixed"), dsem_mu_spec = c(env = 0)))
    o <- suppressWarnings(suppressMessages(fit_model(d$data, d$par, d$map, random = random, do_optim = TRUE,
                                                     newton_loops = 1, silent = TRUE)))
    # as.numeric: a Laplace objective comes back with a logarithm attribute the derived route has no
    list(nll = as.numeric(o$optim$objective), beta = as.numeric(o$env$parList()$dsem_beta))
  }

  zero <- fit_with(", 0, NA, 0", NULL)                # no innovation, so nothing is integrated
  tiny <- fit_with(", 0, NA, 1e-4", "ln_srv_q_devs")  # the same relation as a vanishing random effect

  # the two routes are different code paths onto the same model, so they agree in the limit
  expect_equal(tiny$nll, zero$nll, tolerance = 1e-4)
  expect_equal(tiny$beta, zero$beta, tolerance = 1e-6)
})

test_that("catchability deviations leave a conditioning period alone", {

  n_yrs <- 30
  n_cond <- 20
  sl <- list(n_regions = 1, n_yrs = n_yrs, n_srv_fleets = 1, n_fish_fleets = 1,
             srv_q = array(0.05, dim = c(1, n_yrs, 1, 2)), fish_q = array(0.01, dim = c(1, n_yrs, 1, 2)))
  sl <- Setup_Sim_q_devs(sl, srv_q_model = "rw", sigma_srv_q = 0.3)

  # the catchability handed in for a conditioning period already holds the fit's deviations, so
  # drawing over those years again would square them
  sl$n_cond_yrs <- n_cond
  env <- Setup_sim_env(sl)
  set.seed(4)
  draw_sim_q_devs(1, env)
  draw_sim_q_devs(2, env)

  devs_1 <- env$ln_srv_q_devs[1,,1,1]
  expect_equal(as.numeric(devs_1[seq_len(n_cond)]), rep(0, n_cond))
  expect_equal(as.numeric(env$srv_q[1, seq_len(n_cond), 1, 1]), rep(0.05, n_cond))

  # past them the deviations are drawn, and each replicate draws its own
  later <- (n_cond + 1):n_yrs
  expect_gt(max(abs(devs_1[later])), 0.05)
  expect_false(isTRUE(all.equal(as.numeric(devs_1[later]), as.numeric(env$ln_srv_q_devs[1, later, 1, 2]))))

  # with no conditioning period every year is drawn, as before
  sl$n_cond_yrs <- NULL
  env_free <- Setup_sim_env(sl)
  set.seed(4)
  draw_sim_q_devs(1, env_free)
  expect_gt(max(abs(env_free$ln_srv_q_devs[1, seq_len(n_cond), 1, 1])), 0.05)
})
