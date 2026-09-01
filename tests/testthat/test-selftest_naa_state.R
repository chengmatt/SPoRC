library(SPoRC)
library(testthat)

# Operating-model and estimation-model agreement for the state-space numbers at age. The simulator
# draws innovations from a covariance and the penalty is the density of that same covariance, so
# the two have to agree on two things: which margin each correlation runs over, and what sigmaNAA
# means. Both are the kind of disagreement that produces a plausible number rather than an error.

naa_sim_env <- function(n_pop = 1, n_regions = 3, n_yrs = 40, n_ages = 8, n_sexes = 1, ...) {
  base <- list(n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_ages = n_ages,
               n_sexes = n_sexes, naa_re_ages = 1:n_ages, naa_re_yrs = 1:n_yrs,
               NAA_re = 1, sigmaNAA = 0.5, naa_rho = c(age = 0, year = 0, cohort = 0),
               NAA_re_pop = 0, NAA_re_region = 0, NAA_re_sex = 0,
               naa_pop_corr = 0, naa_region_corr = 0, naa_sex_corr = 0)
  list2env(utils::modifyList(base, list(...)))
}

test_that("simulated innovations carry the correlation on the margin they name", {
  # the cross-checks are the point: 1dar1_a has to show nothing over years and 1dar1_y nothing
  # over ages, which is what a transposed margin would fail
  set.seed(3)
  reps <- function(env, n = 1500) replicate(n, SPoRC:::draw_naa_innovations(env), simplify = FALSE)
  cor_age <- function(x) { na <- dim(x)[4]; stats::cor(as.vector(x[1,,,1:(na-1),1]), as.vector(x[1,,,2:na,1])) }
  cor_yr <- function(x) { ny <- dim(x)[3]; stats::cor(as.vector(x[1,,1:(ny-1),,1]), as.vector(x[1,,2:ny,,1])) }
  avg <- function(L, f) mean(vapply(L, f, numeric(1)))

  D <- reps(naa_sim_env(NAA_re = 1, sigmaNAA = 0.5))
  expect_equal(avg(D, function(x) stats::sd(as.vector(x))), 0.5, tolerance = 0.02)

  D <- reps(naa_sim_env(NAA_re = 2, naa_rho = c(age = 0.7, year = 0, cohort = 0)))
  expect_equal(avg(D, cor_age), 0.7, tolerance = 0.03)
  expect_equal(avg(D, cor_yr), 0, tolerance = 0.03)

  D <- reps(naa_sim_env(NAA_re = 3, naa_rho = c(age = 0, year = 0.6, cohort = 0)))
  expect_equal(avg(D, cor_yr), 0.6, tolerance = 0.03)
  expect_equal(avg(D, cor_age), 0, tolerance = 0.03)
})

test_that("a simulated region correlation appears across regions and leaves ages alone", {
  set.seed(3)
  env <- naa_sim_env(NAA_re = 4, naa_rho = c(age = 0.7, year = 0.4, cohort = 0),
                     NAA_re_region = 1, naa_region_corr = c(0.8, 0.2, 0.5))
  D <- replicate(1500, SPoRC:::draw_naa_innovations(env), simplify = FALSE)
  rc <- function(i, j) mean(vapply(D, function(x) stats::cor(as.vector(x[1,i,,,1]), as.vector(x[1,j,,,1])), numeric(1)))
  expect_equal(rc(1, 2), 0.8, tolerance = 0.03)
  expect_equal(rc(1, 3), 0.2, tolerance = 0.03)
  expect_equal(rc(2, 3), 0.5, tolerance = 0.03)
  na <- 8
  expect_equal(mean(vapply(D, function(x) stats::cor(as.vector(x[1,,,1:(na-1),1]), as.vector(x[1,,,2:na,1])), numeric(1))),
               0.7, tolerance = 0.03)
})

test_that("the penalty recovers the parameters the simulator drew from", {
  # Given the true states, maximizing the penalty must return the process parameters. Averaged over
  # replicates rather than checked on one, because what this is guarding against is bias: if the
  # simulator and the penalty disagree about whether sigmaNAA is the marginal or the conditional
  # standard deviation, every estimate is low by exactly sqrt(1 - rho^2) per correlated margin,
  # which on a single draw is indistinguishable from an unlucky realization.
  set.seed(101)
  np <- 1; nr <- 1; ny <- 40; na <- 10; ns <- 1
  n_rep <- 15
  rtinv <- function(x) 2 / (1 + exp(-2 * x)) - 1

  fit_reps <- function(code, sd_true, rho_true) {
    env <- naa_sim_env(n_pop = np, n_regions = nr, n_yrs = ny, n_ages = na, n_sexes = ns,
                       NAA_re = code, sigmaNAA = sd_true,
                       naa_rho = c(age = rho_true[1], year = rho_true[2], cohort = 0))
    pred <- array(exp(5), dim = c(np, nr, ny, na, ns))
    n_rho <- if(code == 1) 0 else if(code %in% c(2, 3)) 1 else 2
    out <- vapply(seq_len(n_rep), function(i) {
      eta <- SPoRC:::draw_naa_innovations(env)
      nll <- function(th) {
        pe <- array(0, dim = c(np, nr, 3, ns))
        if(code == 2) pe[,,1,] <- th[2]
        if(code == 3) pe[,,2,] <- th[2]
        if(code == 4) { pe[,,1,] <- th[2]; pe[,,2,] <- th[3] }
        SPoRC:::Get_NAA_state_penalty(log(pred) + eta, pred, array(exp(th[1]), dim = dim(pred)),
                                      1:na, 1:ny, NAA_re = code, NAA_pe_pars = pe)
      }
      f <- stats::nlminb(c(log(0.2), rep(0, n_rho)), nll)
      c(exp(f$par[1]), if(n_rho) rtinv(f$par[-1]) else numeric(0))
    }, numeric(1 + n_rho))
    if(is.null(dim(out))) out <- matrix(out, nrow = 1)
    rowMeans(out)
  }

  r <- fit_reps(1, 0.45, c(0, 0))
  expect_equal(r[1], 0.45, tolerance = 0.04)

  r <- fit_reps(2, 0.40, c(0.70, 0))
  expect_equal(r[1], 0.40, tolerance = 0.04)
  expect_equal(r[2], 0.70, tolerance = 0.06)

  r <- fit_reps(3, 0.40, c(0, 0.55))
  expect_equal(r[1], 0.40, tolerance = 0.04)
  expect_equal(r[2], 0.55, tolerance = 0.06)

  r <- fit_reps(4, 0.35, c(0.65, 0.40))
  expect_equal(r[1], 0.35, tolerance = 0.04)
  expect_equal(r[2], 0.65, tolerance = 0.06)
  expect_equal(r[3], 0.40, tolerance = 0.08)
})

test_that("the simulator leaves the numbers at age deterministic when the state is off", {
  env <- naa_sim_env(NAA_re = 0)
  expect_equal(env$NAA_re, 0)
  # Setup_Sim_NAA_state defaults to off and validates its tokens
  sl <- Setup_Sim_NAA_state(list(n_ages = 8, n_yrs = 20))
  expect_equal(sl$NAA_re, 0)
  expect_error(Setup_Sim_NAA_state(list(n_ages = 8, n_yrs = 20), NAA_re = "banana"), "Valid options")
  expect_error(Setup_Sim_NAA_state(list(n_ages = 8, n_yrs = 20), NAA_re_region = "banana"), "Valid options")
})
