library(SPoRC)
library(testthat)

# Seasonal state-space numbers at age. Season one is the year boundary, so the whole feature is a
# strict superset of the annual state and the first thing to pin down is that the superset still
# reproduces the subset. After that: that a season index is not silently a year or an age index,
# that the season correlation runs over seasons and nothing else, and that the setup refuses the
# combinations the density cannot represent.

seas_dims <- list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_seas = 3, n_yrs = 11, n_ages = 6)

seas_build <- function(..., dims = seas_dims) {
  suppressMessages(sweep_input(
    dims = dims,
    stop_after = "biol",
    biol = utils::modifyList(list(NAA_re = "iid", M_spec = "fix"), list(...))
  ))
}

annual_count <- function(il) {
  d <- il$data
  d$n_pop * d$n_regions * length(d$naa_re_yrs) * length(d$naa_re_ages) * d$n_sexes
}

test_that("the default is the annual state and every array holds the season dim", {
  il <- seas_build()
  expect_equal(il$data$naa_re_seas, 1L)
  expect_equal(dim(il$par$ln_NAA), c(1, 1, 11, 3, 6, 1))
  expect_equal(dim(il$data$naa_sigma_blocks), c(1, 1, 11, 3, 6, 1))
  expect_equal(il$data$n_est_naa_re, annual_count(il))

  # only season one is estimated; the rest of the array is kept
  m <- array(il$map$ln_NAA, dim = dim(il$par$ln_NAA))
  expect_true(all(!is.na(m[,,il$data$naa_re_yrs, 1, il$data$naa_re_ages,])))
  expect_true(all(is.na(m[,,,2:3,,])))
})

test_that("all seasons multiplies the state count by the season dim", {
  il <- seas_build(NAA_re_seasons = "all")
  expect_equal(il$data$naa_re_seas, 1:3)
  expect_equal(il$data$n_est_naa_re, 3 * annual_count(il))
  m <- array(il$map$ln_NAA, dim = dim(il$par$ln_NAA))
  expect_true(all(!is.na(m[,,il$data$naa_re_yrs,, il$data$naa_re_ages,])))
})

test_that("the active seasons need not be contiguous", {
  # unlike ages and years: the season dim is only ever independent or unstructured, and neither
  # reads adjacency, so a state can sit only in the seasons that have observations
  il <- seas_build(NAA_re_seasons = c(1, 3))
  expect_equal(il$data$naa_re_seas, c(1L, 3L))
  expect_equal(il$data$n_est_naa_re, 2 * annual_count(il))
  m <- array(il$map$ln_NAA, dim = dim(il$par$ln_NAA))
  expect_true(all(is.na(m[,,,2,,])))
  expect_true(all(!is.na(m[,,il$data$naa_re_yrs, c(1, 3), il$data$naa_re_ages,])))
})

test_that("the setup refuses seasons and correlations it cannot represent", {
  expect_error(seas_build(NAA_re_seasons = c(1, 4)), "season indices")
  expect_error(seas_build(NAA_re_seasons = integer(0)), "season indices")
  # a correlation across one season has nothing to describe
  expect_error(seas_build(NAA_re_season = "us"), "more than one active season")
  expect_error(seas_build(NAA_re_seasons = "all", NAA_re_season = "banana"), "Valid options")
  # a correlated dim has one scale across the dim it spans
  expect_error(seas_build(
    NAA_re_seasons = "all",
    NAA_re_season = "us",
    NAA_sigma_seasblk_spec = list(1, 2:3)
  ),
               "one standard deviation across the dim")
  # blocking on its own is fine, because the season dim is whitened outside the age-year density
  expect_silent(suppressMessages(seas_build(NAA_re_seasons = "all",
                                            NAA_sigma_seasblk_spec = list(1, 2:3))))
})

test_that("the season correlation spends parameters over the active seasons only", {
  il <- seas_build(NAA_re_seasons = "all", NAA_re_season = "us")
  expect_equal(dim(il$par$NAA_season_corr_pars), c(1, 3, 1)) # 3 * 2 / 2 pairs
  expect_equal(length(unique(stats::na.omit(as.integer(as.character(il$map$NAA_season_corr_pars))))), 3)

  two <- seas_build(NAA_re_seasons = c(2, 3), NAA_re_season = "us")
  expect_equal(dim(two$par$NAA_season_corr_pars), c(1, 1, 1)) # one pair over two active seasons
  expect_equal(length(unique(stats::na.omit(as.integer(as.character(two$map$NAA_season_corr_pars))))), 1)

  shared <- seas_build(NAA_re_seasons = "all", NAA_re_season = "us", NAA_re_season_spec = "fix")
  expect_true(all(is.na(as.integer(as.character(shared$map$NAA_season_corr_pars)))))
})

# ---- penalty ---------------------------------------------------------------

test_that("an annual state on a seasonal array penalizes what it does with no season dim", {
  # the season dim must be inert when only season one is live, otherwise every existing model
  # changes value the moment the array grows a dimension
  set.seed(21)
  ny <- 7; na <- 5; nk <- 3
  pred1 <- array(exp(stats::rnorm(ny*na, 5, 0.2)), dim = c(1, 1, ny, 1, na, 1))
  eta1 <- array(stats::rnorm(ny*na, 0, 0.3), dim = c(1, 1, ny, 1, na, 1))
  sig1 <- array(0.35, dim = dim(pred1))
  pe <- array(0.4, dim = c(1, 1, 3, 1))

  pred3 <- array(0, dim = c(1, 1, ny, nk, na, 1)); pred3[,,,1,,] <- pred1
  eta3 <- array(0, dim = c(1, 1, ny, nk, na, 1)); eta3[,,,1,,] <- eta1
  # the inactive seasons hold nonsense on purpose: the slice must never reach them
  pred3[,,,2:3,,] <- exp(9); eta3[,,,2:3,,] <- 7
  sig3 <- array(0.35, dim = dim(pred3))

  for(code in c(1, 2, 3, 4, 5, 6)) {
    one <- SPoRC:::Get_NAA_state_penalty(log(pred1)+eta1, pred1, sig1, 1:na, 1:ny, 1,
                                         NAA_re = code, NAA_pe_pars = pe)
    three <- SPoRC:::Get_NAA_state_penalty(log(pred3)+eta3, pred3, sig3, 1:na, 1:ny, 1,
                                           NAA_re = code, NAA_pe_pars = pe)
    expect_equal(three, one, tolerance = 1e-12, label = paste("NAA_re code", code))
  } # end code loop
})

test_that("independent seasons are penalized as independent replicates of the age-year surface", {
  set.seed(22)
  ny <- 6; na <- 4; nk <- 3
  pred <- array(exp(stats::rnorm(ny*na*nk, 5, 0.2)), dim = c(1, 1, ny, nk, na, 1))
  eta <- array(stats::rnorm(ny*na*nk, 0, 0.3), dim = c(1, 1, ny, nk, na, 1))
  sig <- array(0.3, dim = dim(pred))
  pe <- array(0.5, dim = c(1, 1, 3, 1))

  for(code in c(1, 2, 4, 5)) {
    got <- SPoRC:::Get_NAA_state_penalty(log(pred)+eta, pred, sig, 1:na, 1:ny, 1:nk,
                                         NAA_re = code, NAA_pe_pars = pe)
    parts <- vapply(1:nk, function(k) {
      p1 <- array(pred[,,,k,,], dim = c(1, 1, ny, 1, na, 1))
      e1 <- array(eta[,,,k,,], dim = c(1, 1, ny, 1, na, 1))
      SPoRC:::Get_NAA_state_penalty(log(p1)+e1, p1, array(0.3, dim = dim(p1)), 1:na, 1:ny, 1,
                                    NAA_re = code, NAA_pe_pars = pe)
    }, numeric(1))
    expect_equal(got, sum(parts), tolerance = 1e-10, label = paste("NAA_re code", code))
  } # end code loop
})

test_that("a season correlation at zero reduces to independent seasons", {
  set.seed(23)
  ny <- 6; na <- 4; nk <- 3
  pred <- array(exp(stats::rnorm(ny*na*nk, 5, 0.2)), dim = c(1, 1, ny, nk, na, 1))
  eta <- array(stats::rnorm(ny*na*nk, 0, 0.3), dim = c(1, 1, ny, nk, na, 1))
  sig <- array(0.3, dim = dim(pred))
  pe <- array(0.4, dim = c(1, 1, 3, 1))
  zero_kc <- array(0, dim = c(1, nk*(nk-1)/2, 1))

  for(code in c(1, 2, 3, 4, 5)) {
    off <- SPoRC:::Get_NAA_state_penalty(
      log(pred)+eta,
      pred,
      sig,
      1:na,
      1:ny,
      1:nk,
      NAA_re = code,
      NAA_pe_pars = pe,
      NAA_re_season = 0
    )
    on <- SPoRC:::Get_NAA_state_penalty(
      log(pred)+eta,
      pred,
      sig,
      1:na,
      1:ny,
      1:nk,
      NAA_re = code,
      NAA_pe_pars = pe,
      NAA_re_season = 1,
      NAA_season_corr_pars = zero_kc
    )
    expect_equal(on, off, tolerance = 1e-10, label = paste("NAA_re code", code))
  } # end code loop
})

test_that("the season dim is the season dim and not a year or an age", {
  # Penalized against an explicit Kronecker covariance built independently of the density code, with
  # every extent distinct and every correlation different, so a swapped dim lands visibly wrong.
  skip_if_not_installed("mvtnorm")
  set.seed(24)
  nk <- 3; ny <- 7; na <- 5
  sd_prs <- 0.35; rho_a <- 0.7; rho_y <- 0.2
  rt_inv <- function(r) 0.5 * log((1 + r) / (1 - r))

  pred <- array(exp(stats::rnorm(nk*ny*na, 5, 0.2)), dim = c(1, 1, ny, nk, na, 1))
  eta <- array(stats::rnorm(nk*ny*na, 0, 0.3), dim = c(1, 1, ny, nk, na, 1))
  sig <- array(sd_prs, dim = dim(pred))
  pe <- array(0, dim = c(1, 1, 3, 1))
  pe[1,1,1,1] <- rt_inv(rho_a); pe[1,1,2,1] <- rt_inv(rho_y)
  kc <- array(c(0.6, -0.3, 0.45), dim = c(1, nk*(nk-1)/2, 1))

  got <- SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1:nk,
    NAA_re = 4,
    NAA_pe_pars = pe,
    NAA_re_season = 1,
    NAA_season_corr_pars = kc
  )

  ar1 <- function(n, r) r^abs(outer(1:n, 1:n, "-"))
  C <- SPoRC:::build_us_corr(as.vector(kc), nk)
  scale <- sd_prs / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2)
  v <- as.vector(array(eta[1,1,,,,1], dim = c(ny, nk, na))) # year varies fastest, then season

  right <- -mvtnorm::dmvnorm(v, sigma = scale^2 * kronecker(ar1(na, rho_a), kronecker(C, ar1(ny, rho_y))), log = TRUE)
  expect_equal(got, right, tolerance = 1e-10)

  # a transposed season dim would penalize one of these instead
  swap_yr <- -mvtnorm::dmvnorm(v, sigma = scale^2 * kronecker(ar1(na, rho_a), kronecker(ar1(nk, rho_y), diag(ny))), log = TRUE)
  expect_gt(abs(got - swap_yr), 1)
})

test_that("a season correlation composes with the non-separable three-dimensional field", {
  # the point of whitening rather than forming a Kronecker: the cohort term never has to factor
  set.seed(25)
  nk <- 3; ny <- 6; na <- 4
  pred <- array(exp(stats::rnorm(nk*ny*na, 5, 0.2)), dim = c(1, 1, ny, nk, na, 1))
  eta <- array(stats::rnorm(nk*ny*na, 0, 0.3), dim = c(1, 1, ny, nk, na, 1))
  sig <- array(0.35, dim = dim(pred))
  pe <- array(0.4, dim = c(1, 1, 3, 1))
  kc <- array(c(0.6, -0.3, 0.45), dim = c(1, nk*(nk-1)/2, 1))

  off <- SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1:nk,
    NAA_re = 5,
    NAA_pe_pars = pe,
    NAA_re_season = 0
  )
  on <- SPoRC:::Get_NAA_state_penalty(
    log(pred)+eta,
    pred,
    sig,
    1:na,
    1:ny,
    1:nk,
    NAA_re = 5,
    NAA_pe_pars = pe,
    NAA_re_season = 1,
    NAA_season_corr_pars = kc
  )
  expect_true(is.finite(on))
  expect_gt(abs(on - off), 1)
})

test_that("a season-varying standard deviation reaches the season it belongs to", {
  set.seed(26)
  ny <- 6; na <- 4; nk <- 2
  pred <- array(exp(stats::rnorm(ny*na*nk, 5, 0.2)), dim = c(1, 1, ny, nk, na, 1))
  eta <- array(stats::rnorm(ny*na*nk, 0, 0.3), dim = c(1, 1, ny, nk, na, 1))
  pe <- array(0.4, dim = c(1, 1, 3, 1))
  sig <- array(0, dim = dim(pred)); sig[,,,1,,] <- 0.2; sig[,,,2,,] <- 0.5

  # under a correlated age-year form the sigma is read once per season, so swapping the two
  # seasons' standard deviations has to change the answer by the same amount as swapping the data
  got <- SPoRC:::Get_NAA_state_penalty(log(pred)+eta, pred, sig, 1:na, 1:ny, 1:nk,
                                       NAA_re = 4, NAA_pe_pars = pe)
  sig_sw <- sig; sig_sw[,,,1,,] <- 0.5; sig_sw[,,,2,,] <- 0.2
  eta_sw <- eta; eta_sw[,,,1,,] <- eta[,,,2,,]; eta_sw[,,,2,,] <- eta[,,,1,,]
  pred_sw <- pred; pred_sw[,,,1,,] <- pred[,,,2,,]; pred_sw[,,,2,,] <- pred[,,,1,,]
  both <- SPoRC:::Get_NAA_state_penalty(log(pred_sw)+eta_sw, pred_sw, sig_sw, 1:na, 1:ny, 1:nk,
                                        NAA_re = 4, NAA_pe_pars = pe)
  expect_equal(both, got, tolerance = 1e-10)

  one <- SPoRC:::Get_NAA_state_penalty(log(pred)+eta, pred, sig_sw, 1:na, 1:ny, 1:nk,
                                       NAA_re = 4, NAA_pe_pars = pe)
  expect_gt(abs(one - got), 1)
})

# ---- dynamics --------------------------------------------------------------

# a deterministic pass over the three season model, and the seed it supplies for the state
seas_seed <- local({
  cached <- NULL
  function() {
    if(!is.null(cached)) return(cached)
    il <- suppressMessages(sweep_input(dims = seas_dims, biol = list(M_spec = "fix")))
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    rep <- obj$report(obj$par)
    d <- seas_dims
    cached <<- list(
      il = il,
      obj = obj,
      rep = rep,
      seed = log(array(rep$NAA[,,1:d$n_yrs,,,],
                                     dim = c(1, d$n_regions, d$n_yrs, d$n_seas, d$n_ages, d$n_sexes)))
    )
    cached
  }
})

seas_state_on <- function(...) {
  d <- seas_seed()
  suppressMessages(sweep_input(dims = seas_dims, biol = utils::modifyList(
    list(NAA_re = "iid", M_spec = "fix", ln_NAA = d$seed), list(...))))
}

test_that("seeded at the deterministic solution the seasonal state is a reparameterization", {
  d <- seas_seed()
  il <- seas_state_on(NAA_re_seasons = "all")
  obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
  rep <- obj$report(obj$par)

  # the within-year hook fires on the survival step and the year-boundary one after the plus
  # group, so at the seed the innovation is zero at every active season
  eta <- il$par$ln_NAA[,,il$data$naa_re_yrs, il$data$naa_re_seas, il$data$naa_re_ages,, drop = FALSE] -
         log(rep$NAA_pred[,,il$data$naa_re_yrs, il$data$naa_re_seas, il$data$naa_re_ages,, drop = FALSE])
  expect_lt(max(abs(eta)), 1e-10)

  expect_equal(obj$fn(obj$par) - rep$NAA_state_nLL, d$obj$fn(d$obj$par), tolerance = 1e-8)
  expect_equal(rep$NAA_state_nLL,
               -sum(stats::dnorm(0, 0, 0.3, log = TRUE)) * il$data$n_est_naa_re,
               tolerance = 1e-8)
})

test_that("the state writes only into the seasons it is active over", {
  ann <- seas_state_on()
  r_ann <- fit_model(ann$data, ann$par, ann$map, do_optim = FALSE, silent = TRUE)$report()
  # season one holds the factor, the rest stay on the deterministic trajectory
  expect_true(all(r_ann$NAA_scalar[,,,2:3,,] == 1))
  expect_true(all(r_ann$NAA_pred[,,,2:3,,] == 0))
  expect_true(any(r_ann$NAA_pred[,,ann$data$naa_re_yrs,1,,] > 0))

  gap <- seas_state_on(NAA_re_seasons = c(1, 3))
  r_gap <- fit_model(gap$data, gap$par, gap$map, do_optim = FALSE, silent = TRUE)$report()
  expect_true(all(r_gap$NAA_scalar[,,,2,,] == 1))
  expect_true(all(r_gap$NAA_pred[,,,2,,] == 0))
  expect_true(any(r_gap$NAA_pred[,,gap$data$naa_re_yrs,3,,] > 0))
})

test_that("a within-year innovation lands in its own season and cuts the recursion after it", {
  # this is the whole reason the state is a level rather than a deviation: the season after the
  # bump is set by its own state, so it does not move, and only its prediction does. That is the
  # conditional independence the sparsity argument rests on, now holding within a year as well.
  il <- seas_state_on(NAA_re_seasons = "all")
  base <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$report()

  y <- il$data$naa_re_yrs[3]; a <- il$data$naa_re_ages[2]
  bump <- il
  bump$par$ln_NAA[1, 1, y, 2, a, 1] <- bump$par$ln_NAA[1, 1, y, 2, a, 1] + log(1.4)
  moved <- fit_model(bump$data, bump$par, bump$map, do_optim = FALSE, silent = TRUE)$report()

  expect_equal(moved$NAA[1, 1, y, 2, a, 1] / base$NAA[1, 1, y, 2, a, 1], 1.4, tolerance = 1e-8)
  # season one of the same year is upstream of the change, so it is untouched
  expect_equal(moved$NAA[1, 1, y, 1, a, 1], base$NAA[1, 1, y, 1, a, 1], tolerance = 1e-12)
  # season three is centered on its own state, so the numbers hold and the prediction moves
  expect_equal(moved$NAA[1, 1, y, 3, a, 1], base$NAA[1, 1, y, 3, a, 1], tolerance = 1e-12)
  expect_equal(moved$NAA_pred[1, 1, y, 3, a, 1] / base$NAA_pred[1, 1, y, 3, a, 1], 1.4, tolerance = 1e-8)
})

test_that("under the annual state a season one innovation persists through the year", {
  # the complement: with the seasons deterministic there is nothing to cut the recursion, so the
  # same bump propagates all the way to the next year boundary
  il <- seas_state_on()
  base <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)$report()

  y <- il$data$naa_re_yrs[3]; a <- il$data$naa_re_ages[2]
  bump <- il
  bump$par$ln_NAA[1, 1, y, 1, a, 1] <- bump$par$ln_NAA[1, 1, y, 1, a, 1] + log(1.4)
  moved <- fit_model(bump$data, bump$par, bump$map, do_optim = FALSE, silent = TRUE)$report()

  for(k in 1:3) expect_equal(moved$NAA[1, 1, y, k, a, 1] / base$NAA[1, 1, y, k, a, 1], 1.4,
                             tolerance = 1e-8, label = paste("season", k))
})
