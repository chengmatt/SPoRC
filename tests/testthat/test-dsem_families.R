# Checks the covariate families beyond fixed and normal: each observation density against its base R form through
# SPoRC_rtmb (1e-10), the support checks, the maps and link-scale starts, the operating model draws, and dsem's C++.

fam_input <- function(cov, family, link = NULL, fixed_sd = NULL) {
  input_list <- sweep_input()
  cov_data <- data.frame(year = input_list$data$years, env = cov)
  arrows <- c("env -> rec_Pop_1_Region_2, 0, b", "env <-> env, 0, sd_env", "rec_Pop_1_Region_2 <-> rec_Pop_1_Region_2, 0, sd_r")
  Setup_Mod_DSEM(input_list, dsem_arrows = arrows, dsem_data = cov_data, dsem_family = c(env = family),
                 dsem_link = if(is.null(link)) NULL else c(env = link),
                 dsem_fixed_sd = if(is.null(fixed_sd)) NULL else data.frame(year = input_list$data$years, env = fixed_sd))
}

# dsem_obs_nLL reported by the objective at a chosen covariate state, spread and power
fam_obs_nLL <- function(d, x_state, obs_sd = NULL, logit_p = NULL) {
  d$par$dsem_x[,1] <- x_state
  if(!is.null(obs_sd)) d$par$ln_dsem_obs_sd <- log(obs_sd)
  if(!is.null(logit_p)) d$par$logit_dsem_tweedie_p <- logit_p
  data <- sync_dev_map_data(d$data, d$map)
  obj <- RTMB::MakeADFun(cmb(SPoRC_rtmb, data), d$par, map = d$map, silent = TRUE)
  obj$fn(obj$par)
  obj$report()$dsem_obs_nLL
}

n_fam_yrs <- length(sweep_input()$data$years)

test_that("each family's observation likelihood is its base R density, missing years left out", {

  set.seed(31)
  x <- rnorm(n_fam_yrs, 0, 0.5)

  y <- rbinom(n_fam_yrs, 1, 0.4)
  expect_equal(fam_obs_nLL(fam_input(y, "bernoulli"), x), -sum(dbinom(y, 1, plogis(x), log = TRUE)), tolerance = 1e-10)
  expect_equal(fam_obs_nLL(fam_input(y, "binomial"), x), -sum(dbinom(y, 1, plogis(x), log = TRUE)), tolerance = 1e-10)

  y <- rpois(n_fam_yrs, 3); y[c(2, 5)] <- NA
  seen <- !is.na(y)
  expect_equal(fam_obs_nLL(fam_input(y, "poisson"), x + log(3)), -sum(dpois(y[seen], exp(x + log(3))[seen], log = TRUE)), tolerance = 1e-10)

  y <- rgamma(n_fam_yrs, shape = 4, scale = 1)
  expect_equal(fam_obs_nLL(fam_input(y, "gamma"), x + log(4), obs_sd = 0.5), -sum(dgamma(y, shape = 4, scale = exp(x + log(4)) * 0.25, log = TRUE)), tolerance = 1e-10)

  y <- rlnorm(n_fam_yrs, 1, 0.3)
  expect_equal(fam_obs_nLL(fam_input(y, "lognormal"), x + 1, obs_sd = 0.3), -sum(dlnorm(y, x + 1, 0.3, log = TRUE)), tolerance = 1e-10)

  # tweedie: a zero has the closed form mass exp(-mu^(2-p) / (phi (2-p))), the rest RTMB's own density off the tape
  y <- draw_tweedie(rep(2, n_fam_yrs), 1, 1.5); y[1:3] <- 0
  p <- 1 + plogis(0.4); mu <- exp(x + log(2))
  own <- -sum(RTMB::dtweedie(y[y > 0], mu[y > 0], 1.2, p, TRUE)) + sum(mu[y == 0]^(2 - p) / (1.2 * (2 - p)))
  expect_equal(fam_obs_nLL(fam_input(y, "tweedie"), x + log(2), obs_sd = 1.2, logit_p = 0.4), own, tolerance = 1e-8)

  # dsem's names, its other links, and the fixed-sd normal
  y <- rbinom(n_fam_yrs, 1, 0.4)
  expect_equal(fam_obs_nLL(fam_input(y, "bernoulli", link = "cloglog"), x), -sum(dbinom(y, 1, 1 - exp(-exp(x)), log = TRUE)), tolerance = 1e-10)
  y <- rpois(n_fam_yrs, 3)
  expect_equal(fam_obs_nLL(fam_input(y, "poisson", link = "identity"), exp(x) + 2), -sum(dpois(y, exp(x) + 2, log = TRUE)), tolerance = 1e-10)
  y <- rgamma(n_fam_yrs, shape = 4, scale = 1)
  expect_equal(fam_obs_nLL(fam_input(y, "Gamma"), x + log(4), obs_sd = 0.5), fam_obs_nLL(fam_input(y, "gamma"), x + log(4), obs_sd = 0.5), tolerance = 1e-12)
  y <- rnorm(n_fam_yrs, 2, 0.5)
  expect_equal(fam_obs_nLL(fam_input(y, "gaussian"), x + 2, obs_sd = 0.5), -sum(dnorm(y, x + 2, 0.5, log = TRUE)), tolerance = 1e-10)
  known_sd <- seq(0.2, 0.8, length.out = n_fam_yrs)
  expect_equal(fam_obs_nLL(fam_input(y, "gaussian_fixed_sd", fixed_sd = known_sd), x + 2), -sum(dnorm(y, x + 2, known_sd, log = TRUE)), tolerance = 1e-10)
  expect_equal(fam_obs_nLL(fam_input(abs(y), "gaussian_fixed_sd", link = "log", fixed_sd = known_sd), x), -sum(dnorm(abs(y), exp(x), known_sd, log = TRUE)), tolerance = 1e-10)

})

test_that("values outside a family's support and unknown families are refused", {

  y <- rep(c(0, 1), length.out = n_fam_yrs)
  expect_error(fam_input(replace(y, 1, 2), "bernoulli"), "0 or 1")
  expect_error(fam_input(y + 0.5, "poisson"), "non-negative integers")
  expect_error(fam_input(y, "gamma"), "positive")
  expect_error(fam_input(y, "lognormal"), "positive")
  expect_error(fam_input(y - 1, "tweedie"), "non-negative")
  expect_error(fam_input(y, "beta"), "one of")
  expect_error(fam_input(y, "poisson", link = "inverse"), "one of")
  expect_error(fam_input(y, "gaussian_fixed_sd"), "needs dsem_fixed_sd")
  expect_error(fam_input(y, "gaussian_fixed_sd", fixed_sd = replace(rep(0.3, n_fam_yrs), 4, NA)), "every observed year")
  expect_silent(suppressMessages(fam_input(y, "tweedie")))

})

test_that("maps and starting values follow the family", {

  set.seed(32)
  y <- rpois(n_fam_yrs, 5)
  d <- fam_input(y, "poisson")
  expect_true(all(is.na(d$map$ln_dsem_obs_sd)))
  expect_true(all(is.na(d$map$logit_dsem_tweedie_p)))
  expect_equal(d$par$dsem_mu[1], log(mean(y)))
  expect_true(all(d$par$dsem_x[,1] == log(mean(y)))) # link-scale families start every cell at the mean
  expect_true(all(!is.na(matrix(as.integer(d$map$dsem_x), n_fam_yrs)[,1]))) # and every cell is latent

  d <- fam_input(rbinom(n_fam_yrs, 1, 0.3), "bernoulli")
  expect_equal(d$par$dsem_mu[1], qlogis(mean(d$data$dsem_cov_obs[,1])))

  y <- rgamma(n_fam_yrs, 4, scale = 2)
  d <- fam_input(y, "gamma")
  expect_false(is.na(d$map$ln_dsem_obs_sd[1]))
  expect_true(is.na(d$map$logit_dsem_tweedie_p[1]))
  expect_equal(d$par$ln_dsem_obs_sd, log(sd(y) / mean(y) / 2))

  d <- fam_input(y, "tweedie")
  expect_false(is.na(d$map$ln_dsem_obs_sd[1]))
  expect_false(is.na(d$map$logit_dsem_tweedie_p[1]))
  expect_equal(d$par$logit_dsem_tweedie_p, 0)

  y <- rnorm(n_fam_yrs)
  d <- fam_input(y, "normal")
  expect_equal(d$par$dsem_x[,1], y) # a normal covariate starts at its observations
  expect_equal(d$data$dsem_cov_link, 0)
  expect_length(fam_input(y, "fixed")$par$logit_dsem_tweedie_p, 1)

  # links: the family's own by default, the mean through the link at the start, and the fixed sd stored on the grid
  expect_equal(fam_input(rpois(n_fam_yrs, 2), "poisson")$data$dsem_cov_link, 1)
  d <- fam_input(rbinom(n_fam_yrs, 1, 0.3), "bernoulli", link = "cloglog")
  expect_equal(d$data$dsem_cov_link, 3)
  expect_equal(d$par$dsem_mu[1], log(-log(1 - mean(d$data$dsem_cov_obs[,1]))))
  d <- fam_input(y, "gaussian_fixed_sd", fixed_sd = rep(0.3, n_fam_yrs))
  expect_true(all(is.na(d$map$ln_dsem_obs_sd)))
  expect_equal(d$data$dsem_cov_fixed_sd[,1], rep(0.3, n_fam_yrs), ignore_attr = TRUE)
  expect_equal(d$par$dsem_x[,1], y)

})

test_that("the operating model draws each family about the link-scale state", {

  set.seed(33)
  state <- array(log(3), dim = c(1, 1, 40000))
  y <- draw_dsem_cov_obs(state, 2, 2, 1, 1.5)
  expect_true(all(y %in% c(0, 1)))
  expect_equal(mean(y), plogis(log(3)), tolerance = 0.03)
  expect_equal(mean(draw_dsem_cov_obs(state - 2, 2, 3, 1, 1.5)), 1 - exp(-exp(log(3) - 2)), tolerance = 0.03) # cloglog

  y <- draw_dsem_cov_obs(state, 5, 0, 1, 1.5, fixed_sd = 0.2) # the fixed-sd normal, identity link
  expect_equal(mean(y), log(3), tolerance = 0.03)
  expect_equal(sd(y), 0.2, tolerance = 0.05)

  y <- draw_dsem_cov_obs(state, 3, 1, 1, 1.5)
  expect_true(all(y == round(y) & y >= 0))
  expect_equal(mean(y), 3, tolerance = 0.03)

  y <- draw_dsem_cov_obs(state, 4, 1, 0.4, 1.5)
  expect_true(all(y > 0))
  expect_equal(mean(y), 3, tolerance = 0.03)
  expect_equal(sd(y) / mean(y), 0.4, tolerance = 0.05)

  y <- draw_dsem_cov_obs(state, 6, 1, 0.3, 1.5)
  expect_equal(mean(log(y)), log(3), tolerance = 0.03)
  expect_equal(sd(log(y)), 0.3, tolerance = 0.05)

  # tweedie: mean mu, variance phi mu^p, and zeros with probability exp(-mu^(2-p) / (phi (2-p)))
  y <- draw_dsem_cov_obs(state, 7, 1, 1.5, 1.4)
  expect_equal(dim(y), dim(state))
  expect_equal(mean(y), 3, tolerance = 0.03)
  expect_equal(var(y), 1.5 * 3^1.4, tolerance = 0.06)
  expect_equal(mean(y == 0), exp(-3^0.6 / (1.5 * 0.6)), tolerance = 0.05)

})

test_that("an operating model written from scratch takes a family, a spread and a power", {

  sim_list <- list(n_yrs = 20, n_pop = 1, n_regions = 1)
  arrows <- c("env -> rec, 0, b", "env -> env, 1, rho", "env <-> env, 0, sd_env", "rec <-> rec, 0, sd_rec")
  values <- c(b = 0.3, rho = 0.5, sd_env = 0.4, sd_rec = 0.6)

  scratch <- scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = c(env = log(3)), dsem_cov_obs_sd = c(env = NA), mod_var_logscale = FALSE,
                              dsem_cov_family = c(env = "poisson"))
  expect_equal(scratch$data$dsem_cov_family, 3)
  expect_equal(scratch$data$dsem_cov_link, 1)
  expect_equal(scratch$pars$logit_dsem_tweedie_p, 0)

  scratch <- scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = c(env = 0), dsem_cov_obs_sd = c(env = 0.3), mod_var_logscale = FALSE,
                              dsem_cov_family = c(env = "gaussian_fixed_sd"), dsem_cov_link = c(env = "log"))
  expect_equal(scratch$data$dsem_cov_link, 1)
  expect_equal(scratch$data$dsem_cov_fixed_sd[,1], rep(0.3, 20))
  expect_error(scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = NULL, dsem_cov_obs_sd = c(env = 1), mod_var_logscale = FALSE,
                                dsem_cov_link = c(env = "probit")), "one of")

  scratch <- scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = c(env = log(3)), dsem_cov_obs_sd = c(env = 1.2), mod_var_logscale = FALSE,
                              dsem_cov_family = c(env = "tweedie"), dsem_cov_tweedie_p = c(env = 1.7))
  expect_equal(scratch$data$dsem_cov_family, 7)
  expect_equal(scratch$pars$ln_dsem_obs_sd, log(1.2))
  expect_equal(scratch$pars$logit_dsem_tweedie_p, qlogis(0.7))

  # the defaults: an sd means normal error, NA means known
  scratch <- scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = NULL, dsem_cov_obs_sd = c(env = 0.2), mod_var_logscale = FALSE)
  expect_equal(scratch$data$dsem_cov_family, 1)
  expect_error(scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = NULL, dsem_cov_obs_sd = c(env = NA), mod_var_logscale = FALSE,
                                dsem_cov_family = c(env = "gamma")), "dsem_cov_obs_sd is needed")
  expect_error(scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = NULL, dsem_cov_obs_sd = c(env = 1), mod_var_logscale = FALSE,
                                dsem_cov_family = c(env = "tweedie"), dsem_cov_tweedie_p = c(env = 2)), "between 1 and 2")
  expect_error(scratch_dsem_fit(sim_list, arrows, values, "rec", dsem_cov_mu = NULL, dsem_cov_obs_sd = c(env = 1), mod_var_logscale = FALSE,
                                dsem_cov_family = c(env = "beta")), "one of")

})

test_that("the families match dsem's C++ at the same grid", {

  skip_if_not_installed("dsem")
  set.seed(34)
  n_t <- 20
  sem <- c("x -> y, 0, b", "x -> x, 1, rho", "x <-> x, 0, sx", "y <-> y, 0, sy")
  state <- as.numeric(stats::arima.sim(list(ar = 0.5), n_t)) # the latent covariate, on the link scale
  y_dev <- 0.4 * state + rnorm(n_t, 0, 0.7)
  m <- read_dsem_arrows(sem, c("x", "y"))
  grid <- cbind(x = state, y = y_dev)
  gmrf <- as.numeric(get_dsem_nLL(c(0.4, 0.5), log(c(0.9, 0.7)), x_grid = grid, mu_grid = matrix(c(0.1, -0.2), n_t, 2, byrow = TRUE), dsem_model = m, dsem_cells = get_dsem_cells(m, n_t)))

  # dsem's joint density at the same grid: its objective's full parameter vector, random effects included
  theirs <- function(obs, family, lnsigma = numeric(0), grid_used = grid) { # family: a family object with its link, as dsem takes it
    tsdata <- stats::ts(data.frame(x = obs, y = y_dev), start = 1)
    ctl <- dsem::dsem_control(quiet = TRUE, run_model = FALSE, gmrf_parameterization = "full", use_REML = FALSE)
    d <- suppressWarnings(dsem::dsem(sem = sem, tsdata = tsdata, family = list(x = family, y = dsem::fixed()), estimate_mu = c("x", "y"), control = ctl))
    par <- d$obj$env$par # y's fixed column is data there, so x_tj holds the covariate column only
    expect_equal(sum(names(par) == "x_tj"), n_t)
    par[names(par) == "beta_z"] <- c(0.4, 0.5, 0.9, 0.7)
    if(length(lnsigma) > 0) par[names(par) == "lnsigma_z"] <- lnsigma
    par[names(par) == "mu_j"] <- c(0.1, -0.2)
    par[names(par) == "x_tj"] <- grid_used[,"x"]
    as.numeric(d$obj$env$f(par))
  }

  obs <- rpois(n_t, exp(state + 1))
  expect_equal(gmrf + as.numeric(get_dsem_obs_nLL(obs, state, 3, 1, 1, 1.5)), theirs(obs, stats::poisson(link = "log")), tolerance = 1e-8)

  obs <- rbinom(n_t, 1, plogis(state))
  expect_equal(gmrf + as.numeric(get_dsem_obs_nLL(obs, state, 2, 2, 1, 1.5)), theirs(obs, stats::binomial(link = "logit")), tolerance = 1e-8)

  obs <- rgamma(n_t, shape = 4, scale = exp(state) / 4)
  expect_equal(gmrf + as.numeric(get_dsem_obs_nLL(obs, state, 4, 1, 0.5, 1.5)), theirs(obs, stats::Gamma(link = "log"), log(0.5)), tolerance = 1e-8)

  obs <- exp(state + rnorm(n_t, 0, 0.3))
  expect_equal(gmrf + as.numeric(get_dsem_obs_nLL(obs, state, 6, 1, 0.3, 1.5)), theirs(obs, dsem::lognormal(), log(0.3)), tolerance = 1e-8)

  # dsem's tweedie has dispersion exp(exp(lnsigma1)) and power 1 + plogis(exp(lnsigma2)); ours exp(ln_dsem_obs_sd) and 1 + plogis(logit p)
  obs <- draw_tweedie(exp(state), 1, 1.5)
  expect_equal(gmrf + as.numeric(get_dsem_obs_nLL(obs, state, 7, 1, exp(0.5), 1 + plogis(0.8))), theirs(obs, dsem::tweedie(), log(c(0.5, 0.8))), tolerance = 1e-8)

  # the other links, and the fixed-sd normal
  obs <- rbinom(n_t, 1, 1 - exp(-exp(state)))
  expect_equal(gmrf + as.numeric(get_dsem_obs_nLL(obs, state, 2, 3, 1, 1.5)), theirs(obs, stats::binomial(link = "cloglog")), tolerance = 1e-8)
  # the identity link hands the Poisson the cell itself, so the covariate column is shifted positive for this one
  grid_pos <- cbind(x = state + 4, y = y_dev)
  gmrf_pos <- as.numeric(get_dsem_nLL(c(0.4, 0.5), log(c(0.9, 0.7)), x_grid = grid_pos, mu_grid = matrix(c(0.1, -0.2), n_t, 2, byrow = TRUE), dsem_model = m, dsem_cells = get_dsem_cells(m, n_t)))
  obs <- rpois(n_t, state + 4)
  expect_equal(gmrf_pos + as.numeric(get_dsem_obs_nLL(obs, state + 4, 3, 0, 1, 1.5)), theirs(obs, stats::poisson(link = "identity"), grid_used = grid_pos), tolerance = 1e-8)
  known_sd <- seq(0.2, 0.6, length.out = n_t)
  obs <- state + rnorm(n_t, 0, known_sd)
  expect_equal(gmrf + as.numeric(get_dsem_obs_nLL(obs, state, 5, 0, 1, 1.5, fixed_sd = known_sd)), theirs(obs, dsem::gaussian_fixed_sd(sd = known_sd)), tolerance = 1e-8)

})
