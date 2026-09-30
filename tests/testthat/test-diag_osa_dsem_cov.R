# Checks OSA residuals for the covariates a dsem observes with error: that the packed vector and the
# objective read the same years in the same order, the labels, the count and tweedie families,
# and the refusals.

data("sgl_rg_dusky_data")

dsem_random <- c("ln_RecDevs", "dsem_x")

# a normal and a lognormal covariate observed in different years, so a residual put on the wrong
# covariate or the wrong year cannot line up by accident
dsem_cov_setup <- function(family = c(env = "normal", cpi = "lognormal"),
                           cov_data = NULL) {

  input_list <- build_goa_dusky_input(sgl_rg_dusky_data)
  yrs <- input_list$data$years

  if(is.null(cov_data)) {
    set.seed(11)
    cov_data <- data.frame(year = yrs, env = rnorm(length(yrs)), cpi = rlnorm(length(yrs), 0, 0.3))
    cov_data$env[c(3, 8)] <- NA
    cov_data$cpi[seq(2, length(yrs), by = 2)] <- NA # only the odd years, against the normal one's every year but two
  }

  arrows <- c("env -> rec, 0, b_env", "cpi -> rec, 1, b_cpi", "env <-> env, 0, sd_env",
              "cpi <-> cpi, 0, sd_cpi", "rec <-> rec, 0, NA, 0.9")
  suppressMessages(Setup_Mod_DSEM(input_list, arrows, cov_data, dsem_mu_spec = "fix", dsem_family = family))
}

test_that("the packed vector holds the same years, in the same order, the objective reads per covariate", {

  d <- dsem_cov_setup()
  obs_mat <- d$data$dsem_cov_obs
  packed <- pack_dsem_cov_osa(obs_mat, d$data$dsem_cov_family, "continuous")

  # every observed cell of both covariates, and the value at the cell its map names
  expect_equal(nrow(packed$map), sum(!is.na(obs_mat)))
  expect_equal(packed$vec, obs_mat[cbind(packed$map[,"grid_yr"], packed$map[,"cov"])])

  # the objective slices this vector by covariate and takes its mean from which(!is.na(obs)), so the
  # two have to name the same years in the same order or a residual lands on the wrong cell
  for(k in seq_len(ncol(obs_mat))) {
    take <- which(packed$map[,"cov"] == k)
    expect_equal(packed$map[take,"grid_yr"], which(!is.na(obs_mat[,k])), ignore_attr = TRUE)
  }

  # a count covariate goes in its own vector, and a family with no covariate on it is empty
  counts <- data.frame(year = d$data$years, env = stats::rpois(length(d$data$years), 8),
                       cpi = stats::rlnorm(length(d$data$years), 0, 0.3))
  d2 <- dsem_cov_setup(family = c(env = "poisson", cpi = "lognormal"), cov_data = counts)
  expect_equal(unique(pack_dsem_cov_osa(d2$data$dsem_cov_obs, d2$data$dsem_cov_family, "poisson")$map[,"cov"]), 1L)
  expect_equal(unique(pack_dsem_cov_osa(d2$data$dsem_cov_obs, d2$data$dsem_cov_family, "continuous")$map[,"cov"]), 2L)
  expect_null(pack_dsem_cov_osa(d2$data$dsem_cov_obs, d2$data$dsem_cov_family, "bernoulli"))
  expect_null(pack_dsem_cov_osa(d$data$dsem_cov_obs, d$data$dsem_cov_family, "poisson"))
  expect_null(pack_dsem_cov_osa(d$data$dsem_cov_obs, d$data$dsem_cov_family, "tweedie"))
  expect_error(pack_dsem_cov_osa(d$data$dsem_cov_obs, d$data$dsem_cov_family, "discrete"), "must be")

})

test_that("a residual comes back for every observed covariate year, labelled by covariate and family", {

  d <- dsem_cov_setup()
  fit <- suppressWarnings(fit_model(d$data, d$par, d$map, random = dsem_random, newton_loops = 0, silent = TRUE))
  osa <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE)) # the inner optimization probes bad states
  yrs <- fit$data$years

  expect_equal(nrow(osa$res), sum(!is.na(fit$data$dsem_cov_obs))) # one per observed cell, both covariates
  expect_setequal(unique(osa$res$covariate), c("env", "cpi"))
  expect_equal(unique(osa$res$family[osa$res$covariate == "env"]), "normal")
  expect_equal(unique(osa$res$family[osa$res$covariate == "cpi"]), "lognormal")
  expect_equal(unique(osa$res$idx_type), "DsemCov")

  # each covariate's own unobserved years are the ones left out
  expect_equal(sort(setdiff(yrs, osa$res$year[osa$res$covariate == "env"])), yrs[c(3, 8)])
  expect_equal(sort(osa$res$year[osa$res$covariate == "cpi"]), yrs[seq(1, length(yrs), by = 2)])

  expect_true(all(is.finite(osa$res$resid)))
  expect_lt(abs(stats::sd(osa$res$resid) - 1), 0.4) # standard normal under the model, so near one

  p <- plot_resids(osa)
  expect_named(p, c("sdnr_plot", "resid_plot"))

})

test_that("a bernoulli and a poisson covariate in one model each take their own support", {

  set.seed(13)
  yrs <- build_goa_dusky_input(sgl_rg_dusky_data)$data$years
  both <- data.frame(year = yrs, env = stats::rpois(length(yrs), 8), cpi = stats::rbinom(length(yrs), 1, 0.4))

  d <- dsem_cov_setup(family = c(env = "poisson", cpi = "bernoulli"), cov_data = both)
  fit <- suppressWarnings(fit_model(d$data, d$par, d$map, random = dsem_random, newton_loops = 0, silent = TRUE))

  # the two counts cannot share a support, so each is asked for on its own and each comes back labelled
  # with its own covariate. this pairing leaves both states weakly identified, so what is checked here is
  # the routing, and the single covariate tests below are where a residual is read
  pois <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE, family = "poisson"))
  bern <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE, family = "bernoulli"))

  expect_equal(nrow(pois$res), length(yrs))
  expect_equal(nrow(bern$res), length(yrs))
  expect_equal(unique(pois$res$covariate), "env")
  expect_equal(unique(bern$res$covariate), "cpi")
  expect_equal(unique(pois$res$family), "poisson")
  expect_equal(unique(bern$res$family), "bernoulli")

  # and the continuous vector is empty, since neither covariate is on one of its families
  expect_warning(out <- get_osa(model = fit, data = fit$data, dsem = TRUE), "nothing to compute")
  expect_null(out)

})

test_that("a count covariate on its own gives standard normal residuals", {

  set.seed(17)
  yrs <- build_goa_dusky_input(sgl_rg_dusky_data)$data$years
  n <- length(yrs)

  for(fam in c("bernoulli", "poisson")) {

    obs <- if(fam == "bernoulli") stats::rbinom(n, 1, 0.4) else stats::rpois(n, 8)
    one <- data.frame(year = yrs, env = obs, cpi = stats::rlnorm(n, 0, 0.3))
    d <- dsem_cov_setup(family = c(env = fam, cpi = "fixed"), cov_data = one)
    fit <- suppressWarnings(fit_model(d$data, d$par, d$map, random = dsem_random, newton_loops = 0, silent = TRUE))

    r <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE, family = fam))
    expect_equal(nrow(r$res), n)
    expect_true(all(is.finite(r$res$resid)))
    expect_lt(abs(stats::sd(r$res$resid) - 1), 0.4) # standard normal under the model
    expect_named(plot_resids(r), c("sdnr_plot", "resid_plot"))

  } # end fam loop

})

test_that("a method override is refused on the dsem path, whatever the covariate", {

  # every covariate takes oneStepGeneric. the Gaussian methods read the conditional density's curvature as
  # its variance, and the state reaches the observation through the population model, so that is wrong
  # even for a normal covariate on the identity link
  d <- dsem_cov_setup()
  fit <- fit_model(d$data, d$par, d$map, random = dsem_random, do_optim = FALSE, silent = TRUE)
  expect_error(get_osa(model = fit, data = fit$data, dsem = TRUE, osa_method = "oneStepGaussian"), "not read")
  expect_error(get_osa(model = fit, data = fit$data, dsem = TRUE, osa_method = "oneStepGaussianOffMode"), "not read")
  expect_error(get_osa(model = fit, data = fit$data, dsem = TRUE, osa_method = "oneStepGeneric"), "not read")

})

test_that("a tweedie covariate's zeros come back as residuals, through the mixed support", {

  # a tweedie reads its density directly rather than through the spline, one Laplace step per integrator
  # evaluation, which on the dusky model runs to hours. the small model exercises the same path in seconds
  small <- sweep_input(dims = list(n_yrs = 15, n_ages = 5, n_regions = 1, n_sexes = 1,
                                   n_fish_fleets = 1, n_srv_fleets = 1, n_seas = 1, n_pop = 1))
  yrs <- small$data$years
  set.seed(14)
  tw <- data.frame(year = yrs, env = stats::rlnorm(length(yrs), 0, 0.3), cpi = stats::rlnorm(length(yrs), 0, 0.3))
  tw$env[c(4, 9)] <- 0 # the point mass, which is what needs the mixed path

  arrows <- c("env -> rec, 0, b_env", "cpi -> rec, 1, b_cpi", "env <-> env, 0, sd_env",
              "cpi <-> cpi, 0, sd_cpi", "rec <-> rec, 0, NA, 0.9")
  d <- suppressMessages(Setup_Mod_DSEM(small, arrows, tw, dsem_mu_spec = "fix", dsem_family = c(env = "tweedie", cpi = "lognormal")))
  fit <- suppressWarnings(fit_model(d$data, d$par, d$map, random = dsem_random, newton_loops = 0, silent = TRUE))

  osa <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE, family = "tweedie"))
  expect_equal(nrow(osa$res), length(yrs)) # the tweedie covariate alone, its zeros included
  expect_equal(unique(osa$res$family), "tweedie")
  expect_true(all(is.finite(osa$res$resid))) # the zeros are the ones a continuous call would lose

  # the lognormal one is in its own vector, so the tweedie's point mass does not reach it
  cont <- suppressWarnings(get_osa(model = fit, data = fit$data, dsem = TRUE, family = "continuous"))
  expect_equal(unique(cont$res$covariate), "cpi")

  # the mixed path is for the generic method alone
  expect_error(get_osa(model = fit, data = fit$data, dsem = TRUE, family = "tweedie",
                       osa_method = "oneStepGaussian"), "oneStepGeneric")

})

test_that("a fixed covariate and a model without a dsem have no residuals", {

  d <- dsem_cov_setup(family = c(env = "fixed", cpi = "fixed"))
  fit <- fit_model(d$data, d$par, d$map, random = dsem_random, do_optim = FALSE, silent = TRUE)
  expect_warning(out <- get_osa(model = fit, data = fit$data, dsem = TRUE), "nothing to compute")
  expect_null(out)
  expect_warning(out <- get_osa(model = fit, data = fit$data, dsem = TRUE, family = "tweedie"), "nothing to compute")
  expect_null(out)

  no_dsem <- build_goa_dusky_input(sgl_rg_dusky_data)
  plain <- fit_model(no_dsem$data, no_dsem$par, no_dsem$map, do_optim = FALSE, silent = TRUE)
  expect_warning(out <- get_osa(model = plain, data = plain$data, dsem = TRUE), "no dsem")
  expect_null(out)

})
