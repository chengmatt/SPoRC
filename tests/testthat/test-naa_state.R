library(SPoRC)
library(testthat)

# State-space numbers at age. What these guard, in order: that the feature is inert when off,
# that it is a reparameterization rather than a different model when seeded at the deterministic
# solution, that every correlation form reduces to independence when its correlations are zero,
# and that the age and year dims are not transposed. That last one is the only defect here that
# no objective value, gradient or convergence diagnostic can see.

naa_dims <- list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_seas = 1, n_yrs = 22, n_ages = 8)

# a deterministic pass, and the seed it supplies for the state
naa_seed <- local({
  cached <- NULL
  function() {
    if(!is.null(cached)) return(cached)
    il <- sweep_input(dims = naa_dims)
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    rep <- obj$report(obj$par)
    n_yrs <- length(il$data$years)
    cached <<- list(
      il = il,
      obj = obj,
      rep = rep,
      seed = log(array(rep$NAA[,,1:n_yrs,1,,],
                                     dim = c(1, 1, n_yrs, 1, length(il$data$ages), 1)))
    )
    cached
  }
})

naa_on <- function(form = "iid", ...) {
  d <- naa_seed()
  # modifyList rather than c(): a concatenated list keeps duplicate names, and [[ resolves to the
  # first of them, so an override passed through ... would be silently dropped
  base <- list(NAA_re = form, M_spec = "fix", ln_NAA = d$seed)
  sweep_input(dims = naa_dims, biol = utils::modifyList(base, list(...)))
}

test_that("the state is inert when it is off", {
  d <- naa_seed()
  expect_equal(d$rep$NAA_state_nLL, 0)
  expect_true(all(d$rep$NAA_pred == 0))
  expect_equal(d$il$data$n_est_naa_re, 0)
})

test_that("seeded at the deterministic solution the state is a reparameterization", {
  d <- naa_seed()
  il <- naa_on("iid")
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  rep <- obj$report(obj$par)

  # the prediction is captured after the plus group accumulates, so the innovation at the seed is
  # zero to machine precision; capturing it earlier would leave the accumulation applied on top
  eta <- il$par$ln_NAA[,,il$data$naa_re_yrs, il$data$naa_re_seas, il$data$naa_re_ages,, drop = FALSE] -
         log(rep$NAA_pred[,,il$data$naa_re_yrs, il$data$naa_re_seas, il$data$naa_re_ages,, drop = FALSE])
  expect_lt(max(abs(eta)), 1e-10)

  # and the fit to data is untouched: the whole difference is the prior on the states
  expect_equal(obj$fn(obj$par) - rep$NAA_state_nLL, d$obj$fn(d$obj$par), tolerance = 1e-8)

  # the penalty at a zero innovation is the density at its mode, times the number of states
  expect_equal(rep$NAA_state_nLL,
               -sum(stats::dnorm(0, 0, 0.3, log = TRUE)) * il$data$n_est_naa_re,
               tolerance = 1e-8)
})

test_that("every correlation form collapses to independence at zero correlation", {
  # the correlation parameters start at zero, which is zero on the correlation scale too, so the
  # separable and Markov forms must return exactly what the independent form does
  base <- NULL
  for(form in c("iid", "1dar1_a", "2dar1", "3dcond")) {
    il <- naa_on(form)
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    val <- obj$fn(obj$par)
    if(is.null(base)) base <- val else expect_equal(val, base, tolerance = 1e-8, label = form)
  } # end form loop
})

test_that("the age and year dims are not transposed", {
  # Deliberately asymmetric: a test with the same structure on both dims passes whether or not
  # they are swapped, so it would catch nothing. Penalized against an explicit Kronecker covariance
  # built independently of the density code.
  skip_if_not_installed("mvtnorm")
  set.seed(11)
  ny <- 9; na <- 6
  rho_a <- 0.75; rho_y <- 0.25; sd_prs <- 0.4
  rt_inv <- function(r) 0.5 * log((1 + r) / (1 - r))

  pred <- array(exp(matrix(stats::rnorm(ny * na, 5, 0.2), ny, na)), dim = c(1, 1, ny, 1, na, 1))
  eta <- array(stats::rnorm(ny * na, 0, 0.3), dim = c(1, 1, ny, 1, na, 1))
  sig <- array(sd_prs, dim = dim(pred))
  pe <- array(0, dim = c(1, 1, 3, 1))
  pe[1,1,1,1] <- rt_inv(rho_a); pe[1,1,2,1] <- rt_inv(rho_y)

  got <- SPoRC:::Get_NAA_state_penalty(log(pred) + eta, pred, sig, 1:na, 1:ny, 1,
                                       NAA_re = 4, NAA_pe_pars = pe)

  ar1 <- function(n, r) r^abs(outer(1:n, 1:n, "-"))
  scale <- sd_prs / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2)
  v <- as.vector(matrix(eta[1,1,,,,1], ny, na)) # year varies fastest, as column-major order gives

  # Cov(as.vector(M)) is the reverse Kronecker of the dim order, so age is on the left
  right <- -mvtnorm::dmvnorm(v, sigma = scale^2 * kronecker(ar1(na, rho_a), ar1(ny, rho_y)), log = TRUE)
  swapped <- -mvtnorm::dmvnorm(v, sigma = scale^2 * kronecker(ar1(na, rho_y), ar1(ny, rho_a)), log = TRUE)

  expect_equal(got, right, tolerance = 1e-10)
  expect_gt(abs(got - swapped), 0.1) # the swap must be visibly wrong, not a rounding difference
})

test_that("the setup refuses what is not identified", {
  d <- naa_seed()
  # unseeded falls back on an equilibrium decay rather than a constant, so the starting numbers
  # decline with age instead of being flat. The plus group is the exception and has to be: it
  # accumulates every older age, so at equilibrium it sits above the age below it whenever
  # exp(-M) / (1 - exp(-M)) exceeds one.
  unseeded <- suppressMessages(sweep_input(dims = naa_dims, biol = list(NAA_re = "iid", M_spec = "fix")))
  by_age <- unseeded$par$ln_NAA[1, 1, 1, 1, , 1]
  na <- length(by_age)
  expect_true(all(diff(by_age[1:(na - 1)]) < 0))
  M_bar <- mean(exp(as.vector(unseeded$par$ln_M)))
  expect_equal(by_age[na] - by_age[na - 1], -M_bar - log(1 - exp(-M_bar)), tolerance = 1e-10)
  expect_gt(by_age[1], log(1))
  # the penalty covers one rectangular slice, so the active cells have to be contiguous
  expect_error(naa_on("iid", NAA_re_years = c(2:5, 10:15)), "contiguous run of years")
  expect_error(naa_on("iid", NAA_re_ages = 1:8), "first age")
  # only the independent form is defined for a standard deviation that varies by year or age
  expect_error(naa_on("2dar1", NAA_sigma_yearblk_spec = list(1:11, 12:22)),
               "correlated structure")
})


# ---------------------------------------------------------------------------
# Correlation across regions
# ---------------------------------------------------------------------------

naa_rg_dims <- list(n_regions = 3, n_sexes = 2, n_fish_fleets = 1, n_seas = 1, n_yrs = 16, n_ages = 7)

naa_rg_seed <- local({
  cached <- NULL
  function() {
    if(!is.null(cached)) return(cached)
    il <- sweep_input(dims = naa_rg_dims)
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    rep <- obj$report(obj$par)
    ny <- length(il$data$years); na <- length(il$data$ages)
    cached <<- log(array(rep$NAA[,,1:ny,1,,], dim = c(1, 3, ny, 1, na, 2)))
    cached
  }
})

naa_rg_on <- function(...) {
  base <- list(NAA_re = "iid", M_spec = "fix", ln_NAA = naa_rg_seed())
  sweep_input(dims = naa_rg_dims, biol = utils::modifyList(base, list(...)))
}

test_that("a region correlation at zero reduces to independent regions", {
  # whitening by the Cholesky factor of the identity is the identity, and its log determinant is
  # zero, so every age-year structure must return exactly what it returns without the region factor
  set.seed(4)
  nr <- 3; ny <- 7; na <- 5
  pred <- array(exp(stats::rnorm(nr*ny*na, 5, 0.2)), dim = c(1, nr, ny, 1, na, 1))
  eta <- array(stats::rnorm(nr*ny*na, 0, 0.3), dim = c(1, nr, ny, 1, na, 1))
  sig <- array(0.35, dim = dim(pred))
  pe <- array(0.4, dim = c(1, nr, 3, 1))
  zero_rc <- array(0, dim = c(1, nr*(nr-1)/2, 1))

  for(code in c(1, 2, 3, 4, 5)) {
    off <- SPoRC:::Get_NAA_state_penalty(
      log(pred)+eta,
      pred,
      sig,
      1:na,
      1:ny,
      1,
      NAA_re = code,
      NAA_pe_pars = pe,
      NAA_re_region = 0
    )
    on <- SPoRC:::Get_NAA_state_penalty(
      log(pred)+eta,
      pred,
      sig,
      1:na,
      1:ny,
      1,
      NAA_re = code,
      NAA_pe_pars = pe,
      NAA_re_region = 1,
      NAA_region_corr_pars = zero_rc
    )
    expect_equal(on, off, tolerance = 1e-10, label = paste("NAA_re code", code))
  } # end code loop
})

test_that("the region correlation composes as region against the age and year grid", {
  # The whitening route has to reproduce an explicit three-factor Kronecker exactly, and has to be
  # visibly different from the same model with the region factor dropped.
  skip_if_not_installed("mvtnorm")
  set.seed(4)
  nr <- 3; ny <- 7; na <- 5
  sd_prs <- 0.35; rho_a <- 0.7; rho_y <- 0.2
  rt_inv <- function(r) 0.5 * log((1 + r) / (1 - r))

  pred <- array(exp(stats::rnorm(nr*ny*na, 5, 0.2)), dim = c(1, nr, ny, 1, na, 1))
  eta <- array(stats::rnorm(nr*ny*na, 0, 0.3), dim = c(1, nr, ny, 1, na, 1))
  sig <- array(sd_prs, dim = dim(pred))
  pe <- array(0, dim = c(1, nr, 3, 1))
  pe[1,,1,1] <- rt_inv(rho_a); pe[1,,2,1] <- rt_inv(rho_y)
  rc <- array(c(0.6, -0.3, 0.45), dim = c(1, nr*(nr-1)/2, 1))

  got <- SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1,
    NAA_re = 4,
    NAA_pe_pars = pe,
    NAA_re_region = 1,
    NAA_region_corr_pars = rc
  )

  ar1 <- function(n, r) r^abs(outer(1:n, 1:n, "-"))
  C <- SPoRC:::build_us_corr(as.vector(rc), nr)
  scale <- sd_prs / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2)
  v <- as.vector(array(eta[1,,,,,1], dim = c(nr, ny, na))) # region varies fastest

  right <- -mvtnorm::dmvnorm(v, sigma = scale^2 * kronecker(ar1(na, rho_a), kronecker(ar1(ny, rho_y), C)), log = TRUE)
  no_region <- -mvtnorm::dmvnorm(v, sigma = scale^2 * kronecker(ar1(na, rho_a), kronecker(ar1(ny, rho_y), diag(nr))), log = TRUE)

  expect_equal(got, right, tolerance = 1e-10)
  expect_gt(abs(got - no_region), 1)
})

test_that("the region correlation composes with the non-separable three-dimensional field", {
  # the point of whitening rather than forming a Kronecker: the cohort term never has to factor
  set.seed(4)
  nr <- 3; ny <- 7; na <- 5
  pred <- array(exp(stats::rnorm(nr*ny*na, 5, 0.2)), dim = c(1, nr, ny, 1, na, 1))
  eta <- array(stats::rnorm(nr*ny*na, 0, 0.3), dim = c(1, nr, ny, 1, na, 1))
  sig <- array(0.35, dim = dim(pred))
  pe <- array(0.4, dim = c(1, nr, 3, 1))
  rc <- array(c(0.6, -0.3, 0.45), dim = c(1, nr*(nr-1)/2, 1))

  off <- SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1,
    NAA_re = 5,
    NAA_pe_pars = pe,
    NAA_re_region = 0
  )
  on <- SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1,
    NAA_re = 5,
    NAA_pe_pars = pe,
    NAA_re_region = 1,
    NAA_region_corr_pars = rc
  )
  expect_true(is.finite(on))
  expect_false(isTRUE(all.equal(on, off)))
})

test_that("the region correlation sharing specs give the parameter counts they claim", {
  # n_pop is one in this test setup, so sharing over population is a no-op and sharing over sex is not
  n_free <- function(il) length(unique(stats::na.omit(as.integer(as.character(il$map$NAA_region_corr_pars)))))
  expect_equal(n_free(naa_rg_on(NAA_re_region = "us", NAA_re_region_spec = "est_all")), 6)
  expect_equal(n_free(naa_rg_on(NAA_re_region = "us", NAA_re_region_spec = "est_shared_p")), 6)
  expect_equal(n_free(naa_rg_on(NAA_re_region = "us", NAA_re_region_spec = "est_shared_s")), 3)
  expect_equal(n_free(naa_rg_on(NAA_re_region = "us", NAA_re_region_spec = "est_shared_p_s")), 3)
  expect_equal(n_free(naa_rg_on(NAA_re_region = "us", NAA_re_region_spec = "fix")), 0)
  expect_equal(n_free(naa_rg_on(NAA_re_region = "iid")), 0)
})

test_that("a region correlation needs more than one region", {
  d <- naa_seed()
  expect_error(naa_on("iid", NAA_re_region = "us"), "more than one region")
  expect_error(naa_rg_on(NAA_re_region = "banana"), "Valid options")
  expect_error(naa_rg_on(NAA_re_region = "us", NAA_re_region_spec = "banana"), "Valid options")
})

test_that("a multi-region state builds on the tape with both dims live", {
  il <- naa_rg_on(NAA_re = "2dar1", NAA_re_region = "us", NAA_re_region_spec = "est_shared_p_s")
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
  expect_equal(il$data$n_est_naa_re, 3 * (length(il$data$years) - 1) * (length(il$data$ages) - 1) * 2)
})

# ---------------------------------------------------------------------------
# Population and sex dims, the selectivity confound, and tag cohorts
# ---------------------------------------------------------------------------

test_that("all four correlation dims compose as one Kronecker product", {
  # The decisive check on the whitening route: four correlated dims plus a separable structure
  # over the age-year grid, against a five-factor covariance built independently of the density
  # code. Every dim has a different structure, so any permutation of them is visible.
  skip_if_not_installed("mvtnorm")
  set.seed(9)
  np <- 2; nr <- 3; ny <- 6; na <- 4; ns <- 2
  sd0 <- 0.35; rho_a <- 0.6; rho_y <- 0.3
  rt <- function(r) 0.5 * log((1 + r) / (1 - r))
  n <- np * nr * ny * na * ns

  pred <- array(exp(stats::rnorm(n, 5, 0.2)), dim = c(np, nr, ny, 1, na, ns))
  eta <- array(stats::rnorm(n, 0, 0.3), dim = c(np, nr, ny, 1, na, ns))
  sig <- array(sd0, dim = dim(pred))
  pe <- array(0, dim = c(np, nr, 3, ns)); pe[,,1,] <- rt(rho_a); pe[,,2,] <- rt(rho_y)

  rcP <- 0.5; rcS <- -0.4; pair <- c(0.6, -0.3, 0.45)
  # filled explicitly: array() would recycle a length-three vector across pop and sex, giving each
  # a different correlation and quietly invalidating the comparison
  rcR <- array(0, dim = c(np, 3, ns))
  for(p in 1:np) for(s in 1:ns) rcR[p,,s] <- pair

  got <- SPoRC:::Get_NAA_state_penalty(
    log(pred) + eta,
    pred,
    sig,
    1:na,
    1:ny,
    1,
    NAA_re = 4,
    NAA_pe_pars = pe,
    NAA_re_region = 1,
    NAA_region_corr_pars = rcR,
    NAA_re_pop = 1,
    NAA_pop_corr_pars = rcP,
    NAA_re_sex = 1,
    NAA_sex_corr_pars = rcS
  )

  ar1 <- function(k, r) r^abs(outer(1:k, 1:k, "-"))
  Cp <- SPoRC:::build_us_corr(rcP, np)
  Cs <- SPoRC:::build_us_corr(rcS, ns)
  Cr <- SPoRC:::build_us_corr(pair, nr)
  sc <- sd0 / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2)
  K <- function(...) Reduce(kronecker, list(...))

  # dims run pop, region, year, age, sex with pop varying fastest, so the covariance of the
  # vectorized array is the reverse of that order
  ref <- -mvtnorm::dmvnorm(as.vector(eta),
                           sigma = sc^2 * K(Cs, ar1(na, rho_a), ar1(ny, rho_y), Cr, Cp), log = TRUE)
  expect_equal(got, ref, tolerance = 1e-10)
})

test_that("a population or sex correlation at zero reduces to independence", {
  set.seed(9)
  np <- 2; nr <- 2; ny <- 5; na <- 4; ns <- 2
  n <- np * nr * ny * na * ns
  pred <- array(exp(stats::rnorm(n, 5, 0.2)), dim = c(np, nr, ny, 1, na, ns))
  eta <- array(stats::rnorm(n, 0, 0.3), dim = c(np, nr, ny, 1, na, ns))
  sig <- array(0.35, dim = dim(pred))
  pe <- array(0.4, dim = c(np, nr, 3, ns))
  off <- SPoRC:::Get_NAA_state_penalty(log(pred)+eta, pred, sig, 1:na, 1:ny, 1, NAA_re = 4, NAA_pe_pars = pe)
  expect_equal(SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1,
    NAA_re = 4,
    NAA_pe_pars = pe,
    NAA_re_pop = 1,
    NAA_pop_corr_pars = 0
  ), off, tolerance = 1e-10)
  expect_equal(SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1,
    NAA_re = 4,
    NAA_pe_pars = pe,
    NAA_re_sex = 1,
    NAA_sex_corr_pars = 0
  ), off, tolerance = 1e-10)
})

test_that("a correlation dim with one level is refused", {
  expect_error(naa_on("iid", NAA_re_pop = "us"), "more than one pop")
  # the single-region test setup also has one population, so sex is the one that can be exercised there
  expect_error(naa_on("iid", NAA_re_sex = "banana"), "Valid options")
})

test_that("tag cohorts are rescaled by the state and are untouched without it", {
  # A scalar of ones has to leave the tagging model bit-identical, which is what makes the rescale
  # safe to apply unconditionally when the state is on.
  il <- naa_rg_on(NAA_re = "iid")
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  rep <- obj$report(obj$par)
  expect_true(all(rep$NAA_scalar > 0))
  # seeded at the deterministic solution the state changes nothing, so the factor is one throughout
  expect_equal(max(abs(rep$NAA_scalar - 1)), 0, tolerance = 1e-8)
})

test_that("the one-dimensional autoregressions run over the dim they name", {
  # 1dar1 correlates ages with years independent, 1dar1_y the reverse. Checked against explicit
  # Kronecker covariances, and against each other: a form that ran over the wrong dim would
  # still give a finite objective and converge, so only the cross-check catches it.
  skip_if_not_installed("mvtnorm")
  set.seed(21)
  ny <- 9; na <- 6; sd0 <- 0.4; rho_a <- 0.75; rho_y <- 0.25
  rt <- function(r) 0.5 * log((1 + r) / (1 - r))

  pred <- array(exp(stats::rnorm(ny*na, 5, 0.2)), dim = c(1, 1, ny, 1, na, 1))
  eta <- array(stats::rnorm(ny*na, 0, 0.3), dim = c(1, 1, ny, 1, na, 1))
  sig <- array(sd0, dim = dim(pred))
  pe <- array(0, dim = c(1, 1, 3, 1))
  pe[1,1,1,1] <- rt(rho_a); pe[1,1,2,1] <- rt(rho_y)

  P <- function(k) SPoRC:::Get_NAA_state_penalty(log(pred)+eta, pred, sig, 1:na, 1:ny, 1,
                                                 NAA_re = k, NAA_pe_pars = pe)
  ar1 <- function(n, r) r^abs(outer(1:n, 1:n, "-"))
  v <- as.vector(matrix(eta[1,1,,,,1], ny, na)) # year varies fastest
  ll <- function(S) -mvtnorm::dmvnorm(v, sigma = S, log = TRUE)

  # ages independent, years correlated
  expect_equal(P(3), ll((sd0/sqrt(1-rho_y^2))^2 * kronecker(diag(na), ar1(ny, rho_y))), tolerance = 1e-10)
  # years independent, ages correlated
  expect_equal(P(2), ll((sd0/sqrt(1-rho_a^2))^2 * kronecker(ar1(na, rho_a), diag(ny))), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(P(2), P(3))))

  # and both are the independent form when their correlation is zero
  pe0 <- array(0, dim = c(1, 1, 3, 1))
  Q <- function(k) SPoRC:::Get_NAA_state_penalty(log(pred)+eta, pred, sig, 1:na, 1:ny, 1,
                                                 NAA_re = k, NAA_pe_pars = pe0)
  expect_equal(Q(2), Q(1), tolerance = 1e-10)
  expect_equal(Q(3), Q(1), tolerance = 1e-10)
})

test_that("1dar1_y estimates the year correlation slot and leaves the others kept", {
  il <- naa_on("1dar1_y")
  mp <- as.integer(as.character(il$map$NAA_pe_pars))
  live <- array(!is.na(mp), dim = dim(il$par$NAA_pe_pars))
  expect_false(any(live[,,1,])) # age slot kept
  expect_true(all(live[,,2,]))  # year slot estimated
  expect_false(any(live[,,3,])) # cohort slot kept
})

test_that("a retrospective peel truncates the state, its map and its active years", {
  # The penalty slices ln_NAA with naa_re_yrs, so an untruncated index vector reads past the end of
  # the shortened array, and n_est_naa_re gates both the dynamics hook and the penalty. Neither can
  # be reused from the full model.
  il <- naa_on("iid")
  n_yrs <- length(il$data$years)
  peel <- 3
  cut <- SPoRC:::truncate_yr(peel, il$data, il$par, il$map)

  expect_equal(dim(cut$retro_parameters$ln_NAA)[3], n_yrs - peel)
  expect_equal(dim(cut$retro_data$map_ln_NAA)[3], n_yrs - peel)
  expect_equal(dim(cut$retro_data$naa_sigma_blocks)[3], n_yrs - peel)
  expect_true(all(cut$retro_data$naa_re_yrs <= n_yrs - peel))
  expect_equal(length(cut$retro_mapping$ln_NAA), length(cut$retro_parameters$ln_NAA))

  # the state count drops by exactly the peeled years' worth of cells
  per_year <- il$data$n_est_naa_re / length(il$data$naa_re_yrs)
  expect_equal(cut$retro_data$n_est_naa_re, per_year * length(cut$retro_data$naa_re_yrs))

  # and the peeled model actually runs, with a smaller penalty than the full one
  peeled <- fit_model(cut$retro_data, cut$retro_parameters, cut$retro_mapping,
                      do_optim = FALSE, silent = TRUE)
  full <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  rp <- peeled$report(peeled$par); rf <- full$report(full$par)
  expect_true(is.finite(rp$NAA_state_nLL))
  # every state sits at its mode here, so the penalty is the same per-state density summed over
  # however many states there are: it scales with the count rather than simply shrinking
  expect_equal(rp$NAA_state_nLL / rf$NAA_state_nLL,
               cut$retro_data$n_est_naa_re / il$data$n_est_naa_re, tolerance = 1e-6)
})

test_that("a peel that removes every active year switches the state off", {
  late <- tail(naa_seed()$il$data$years, 4)
  il <- naa_on("iid", NAA_re_years = late)
  cut <- SPoRC:::truncate_yr(6, il$data, il$par, il$map)
  expect_equal(length(cut$retro_data$naa_re_yrs), 0)
  expect_equal(cut$retro_data$n_est_naa_re, 0)
  expect_true(all(is.na(as.integer(as.character(cut$retro_mapping$ln_NAA)))))
  obj <- fit_model(cut$retro_data, cut$retro_parameters, cut$retro_mapping, do_optim = FALSE, silent = TRUE)
  expect_equal(obj$report(obj$par)$NAA_state_nLL, 0)
})

test_that("the correlation sharing specs give the parameter counts they claim", {
  # n_pop is one in this test setup, so sharing over population is a no-op; region and sex are not.
  # 2dar1 claims two of the three slots, so est_all is n_regions * n_sexes * 2 = 12.
  n_free <- function(il) length(unique(stats::na.omit(as.integer(as.character(il$map$NAA_pe_pars)))))
  expect_equal(n_free(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_all")), 12)
  expect_equal(n_free(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_shared_p")), 12)
  expect_equal(n_free(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_shared_r")), 4)
  expect_equal(n_free(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_shared_s")), 6)
  expect_equal(n_free(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_shared_r_s")), 2)
  expect_equal(n_free(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_shared_p_r_s")), 2)
  expect_equal(n_free(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "fix")), 0)
  expect_equal(n_free(naa_rg_on(NAA_re = "3dcond", NAA_pe_spec = "est_shared_r_s")), 3) # all three slots
  expect_equal(n_free(naa_rg_on(NAA_re = "iid", NAA_pe_spec = "est_all")), 0)           # no slots at all
  expect_error(naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "banana"), "Valid options")
})

test_that("sharing collapses the correlations without moving which cells are estimated", {
  free <- naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_all")
  shared <- naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_shared_r")

  # est_all must still be the plain column-major sequence the flat seq_len used to give
  mp <- as.integer(as.character(free$map$NAA_pe_pars))
  expect_equal(mp[!is.na(mp)], seq_len(sum(!is.na(mp))))

  # sharing over regions gives every region the same parameter, one per slot and sex
  ms <- array(as.integer(as.character(shared$map$NAA_pe_pars)), dim = dim(shared$par$NAA_pe_pars))
  for(k in 1:2) for(s in 1:2) expect_equal(length(unique(ms[1,,k,s])), 1)
  expect_false(any(duplicated(as.vector(ms[1,1,1:2,]))))  # slots and sexes stay distinct

  # the state itself is untouched: same live cells, same count
  expect_equal(free$data$n_est_naa_re, shared$data$n_est_naa_re)
  expect_equal(is.na(free$map$ln_NAA), is.na(shared$map$ln_NAA))
})

test_that("a shared correlation penalizes the same as repeating one value across regions", {
  # sharing is a parameter-count change, not a model change: the objective at a point where every
  # region already holds the same rho must be identical under est_all and est_shared_r
  il_f <- naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_all")
  il_s <- naa_rg_on(NAA_re = "2dar1", NAA_pe_spec = "est_shared_r")
  pe <- array(0, dim = dim(il_f$par$NAA_pe_pars)); pe[,,1,] <- 0.6; pe[,,2,] <- -0.4
  il_f$par$NAA_pe_pars <- pe; il_s$par$NAA_pe_pars <- pe

  o_f <- fit_model(il_f$data, il_f$par, il_f$map, do_optim = FALSE, silent = TRUE)
  o_s <- fit_model(il_s$data, il_s$par, il_s$map, do_optim = FALSE, silent = TRUE)
  expect_equal(o_f$fn(o_f$par), o_s$fn(o_s$par), tolerance = 1e-12)
  expect_equal(length(o_s$par), length(o_f$par) - 8) # 12 correlations collapse to 4
})
