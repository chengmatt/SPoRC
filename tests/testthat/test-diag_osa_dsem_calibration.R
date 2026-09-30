# Checks that a dsem covariate's OSA residuals are standard normal under a correct model, one family at a
# time. Data drawn from the model at known parameters, no refit, replicates as the independent unit.

# a small model, since the residual code is the same whatever the assessment around it
cal_input <- sweep_input(dims = list(n_yrs = 15, n_ages = 5, n_regions = 1, n_sexes = 1,
                                     n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 1, n_pop = 1))
cal_yrs <- cal_input$data$years
cal_n <- length(cal_yrs)
cal_rho <- 0.5
cal_sd <- 0.5
cal_known_sd <- 0.25

# rec has to appear but nothing joins it to env, so env alone is a draw from the model's own prior
cal_arrows <- c(paste0("env -> env, 1, NA, ", cal_rho), paste0("env <-> env, 0, NA, ", cal_sd),
                "rec <-> rec, 0, NA, 0.9")

# the spread each family reads, and where its state sits on the link scale
cal_spread <- c(normal = 0.3, gaussian_fixed_sd = 0, lognormal = 0.3, gamma = 0.4,
                bernoulli = 0, poisson = 0, tweedie = 1.2)
cal_offset <- c(normal = 0, gaussian_fixed_sd = 0, lognormal = 0, gamma = log(5),
                bernoulli = 0, poisson = log(8), tweedie = log(3))

# one replicate: the state from the dsem's own prior (year one the innovation alone, then the ar1
# forward), the observations through the operating model's draw, the residuals at the true parameters
cal_replicate <- function(family, seed) {

  set.seed(seed)
  x <- numeric(cal_n)
  x[1] <- rnorm(1, 0, cal_sd)
  for(t in 2:cal_n) x[t] <- cal_rho * x[t - 1] + rnorm(1, 0, cal_sd)
  x <- x + cal_offset[[family]]

  fam_code <- unname(dsem_family_codes()[family])
  obs <- as.numeric(draw_dsem_cov_obs(array(x, dim = c(cal_n, 1, 1)), fam_code, dsem_default_link(fam_code),
                                      cal_spread[[family]], 1.5,
                                      if(family == "gaussian_fixed_sd") rep(cal_known_sd, cal_n) else NULL))

  d <- suppressMessages(Setup_Mod_DSEM(cal_input, cal_arrows, data.frame(year = cal_yrs, env = obs),
                                       dsem_mu_spec = c(env = cal_offset[[family]]), dsem_family = c(env = family),
                                       dsem_fixed_sd = if(family == "gaussian_fixed_sd") data.frame(year = cal_yrs, env = rep(cal_known_sd, cal_n)) else NULL))
  if(cal_spread[[family]] > 0) d$par$ln_dsem_obs_sd[] <- log(cal_spread[[family]])
  d$map$ln_dsem_obs_sd <- factor(rep(NA, length(d$par$ln_dsem_obs_sd)))
  fit <- suppressWarnings(fit_model(d$data, d$par, d$map, random = c("ln_RecDevs", "dsem_x"), do_optim = FALSE, silent = TRUE))

  # a randomisation seed per replicate, or every discrete residual takes the same uniform draw
  osa_fam <- if(family %in% c("bernoulli", "poisson", "tweedie")) family else "continuous"
  r <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE, family = osa_fam, seed = seed))
  r$res$resid

}

test_that("every covariate family's residuals are standard normal under a correct model", {

  skip_on_cran()
  n_reps <- 10

  for(family in names(cal_spread)) {

    rep_mean <- rep_sd <- numeric(0)
    for(rep_i in seq_len(n_reps)) {
      v <- cal_replicate(family, 1000 * match(family, names(cal_spread)) + rep_i)
      expect_true(all(is.finite(v)), info = family)
      rep_mean <- c(rep_mean, mean(v))
      rep_sd <- c(rep_sd, sd(v))
    } # end rep_i loop

    # a tolerance rather than a hypothesis test, loose enough that a correct model never trips it and
    # tight enough to catch the defects this has caught: a 0.3 shift on the tweedie from splining its
    # density, and a 0.18 shift on the bernoulli from one shared randomisation seed
    expect_lt(abs(mean(rep_mean)), 0.25, label = paste(family, "mean of replicate means"))
    expect_lt(abs(mean(rep_sd) - 1), 0.25, label = paste(family, "mean of replicate sds"))

  } # end family loop

})
