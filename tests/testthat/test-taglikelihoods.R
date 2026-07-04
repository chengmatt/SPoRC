library(testthat)
library(RTMB)

make_data <- function(like_type, ln_theta = log(5), n_tags_released = 50) {

  n_conv_tag_cohorts         <- 1
  conv_tag_release_indicator <- matrix(c(1, 1, 1), nrow = 1) # region, year, season
  conv_tag_max_liberty       <- 2
  n_yrs                      <- 2
  n_seas                     <- 1
  conv_tag_mixing_period     <- 1
  n_fish_fleets                <- 1
  use_conv_fish_tagging      <- c(1)
  n_conv_tag_pop_pool        <- 1
  n_regions                  <- 2
  n_conv_tag_age_pool        <- 1
  n_conv_tag_sex_pool        <- 1
  conv_tag_pop_pool          <- list(1)
  conv_tag_age_pool          <- list(1)
  conv_tag_sex_pool          <- list(1)
  addtotag                   <- 0.001

  arr_dim <- c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, 1, n_regions, 1, 1, n_fish_fleets)
  obs  <- array(0, dim = arr_dim)
  pred <- array(0, dim = arr_dim)

  obs[1, 1, 1, 1, 1, 1, 1, 1] <- 3;  obs[1, 1, 1, 1, 2, 1, 1, 1] <- 5
  obs[2, 1, 1, 1, 1, 1, 1, 1] <- 2;  obs[2, 1, 1, 1, 2, 1, 1, 1] <- 4

  pred[1, 1, 1, 1, 1, 1, 1, 1] <- 2.5; pred[1, 1, 1, 1, 2, 1, 1, 1] <- 4.5
  pred[2, 1, 1, 1, 1, 1, 1, 1] <- 1.8; pred[2, 1, 1, 1, 2, 1, 1, 1] <- 3.9

  conv_fish_tag_nLL <- array(0, dim = c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_regions, n_fish_fleets))
  conv_tagged_fish  <- array(n_tags_released, dim = c(n_conv_tag_cohorts, 1, 1, 1))

  list(
    n_conv_tag_cohorts = n_conv_tag_cohorts,
    conv_tag_release_indicator = conv_tag_release_indicator,
    conv_tag_max_liberty = conv_tag_max_liberty,
    n_yrs = n_yrs, n_seas = n_seas, conv_tag_mixing_period = conv_tag_mixing_period,
    n_fish_fleets = n_fish_fleets, use_conv_fish_tagging = use_conv_fish_tagging,
    n_conv_tag_pop_pool = n_conv_tag_pop_pool, n_regions = n_regions,
    n_conv_tag_age_pool = n_conv_tag_age_pool, n_conv_tag_sex_pool = n_conv_tag_sex_pool,
    conv_tag_pop_pool = conv_tag_pop_pool, conv_tag_age_pool = conv_tag_age_pool,
    conv_tag_sex_pool = conv_tag_sex_pool, conv_fish_tag_like = like_type,
    conv_fish_tag_nLL = conv_fish_tag_nLL,
    obs_recap = obs,
    obs_conv_tag_fish_recap = obs, pred_conv_tag_fish_recap = pred,
    addtotag = addtotag, ln_conv_fish_tag_theta = ln_theta,
    conv_tagged_fish = conv_tagged_fish
  )
}

run_fit <- function(d) {
  get_conv_tag_likelihoods(
    n_conv_tag_cohorts          = d$n_conv_tag_cohorts,
    conv_tag_release_indicator  = d$conv_tag_release_indicator,
    conv_tag_max_liberty        = d$conv_tag_max_liberty,
    n_yrs                       = d$n_yrs,
    n_seas                      = d$n_seas,
    conv_tag_mixing_period      = d$conv_tag_mixing_period,
    n_fish_fleets                 = d$n_fish_fleets,
    use_conv_fish_tagging       = d$use_conv_fish_tagging,
    n_conv_tag_pop_pool         = d$n_conv_tag_pop_pool,
    n_regions                   = d$n_regions,
    n_conv_tag_age_pool         = d$n_conv_tag_age_pool,
    n_conv_tag_sex_pool         = d$n_conv_tag_sex_pool,
    conv_tag_pop_pool           = d$conv_tag_pop_pool,
    conv_tag_age_pool           = d$conv_tag_age_pool,
    conv_tag_sex_pool           = d$conv_tag_sex_pool,
    conv_fish_tag_like          = d$conv_fish_tag_like,
    conv_fish_tag_nLL           = d$conv_fish_tag_nLL,
    obs_conv_tag_fish_recap     = d$obs_conv_tag_fish_recap,
    pred_conv_tag_fish_recap    = d$pred_conv_tag_fish_recap,
    addtotag                    = d$addtotag,
    ln_conv_fish_tag_theta      = d$ln_conv_fish_tag_theta,
    conv_tagged_fish            = d$conv_tagged_fish
  )
}

## Type 0: Poisson (count)
test_that("external path: Poisson (like_type = 0) matches manual noint Poisson nLL", {

  d   <- make_data(like_type = 0)
  out <- run_fit(d)

  expected <- array(0, dim = dim(out))
  for (ry in 1:2) for (r in 1:2) {
    o <- d$obs_conv_tag_fish_recap[ry, 1, 1, 1, r, 1, 1, 1] + d$addtotag
    p <- d$pred_conv_tag_fish_recap[ry, 1, 1, 1, r, 1, 1, 1] + d$addtotag
    expected[ry, 1, 1, r, 1] <- -(o * log(p) - p - lgamma(o + 1))
  }

  expect_equal(dim(out), dim(d$conv_fish_tag_nLL))
  expect_equal(as.numeric(out), as.numeric(expected), tolerance = 1e-6)
  expect_true(all(is.finite(as.numeric(out))))
})

## Type 1: Negative binomial (count)
test_that("external path: negative binomial (like_type = 1) matches manual noint NB2 nLL", {

  d   <- make_data(like_type = 1, ln_theta = log(5))
  out <- run_fit(d)

  theta <- exp(d$ln_conv_fish_tag_theta)
  expected <- array(0, dim = dim(out))
  for (ry in 1:2) for (r in 1:2) {
    o  <- d$obs_conv_tag_fish_recap[ry, 1, 1, 1, r, 1, 1, 1] + d$addtotag
    mu <- d$pred_conv_tag_fish_recap[ry, 1, 1, 1, r, 1, 1, 1] + d$addtotag
    ll <- lgamma(o + theta) - lgamma(theta) - lgamma(o + 1) +
      theta * log(theta / (theta + mu)) + o * log(mu / (theta + mu))
    expected[ry, 1, 1, r, 1] <- -ll
  }

  expect_equal(as.numeric(out), as.numeric(expected), tolerance = 1e-6)
  expect_true(all(is.finite(as.numeric(out))))
})

## Type 2: Multinomial, release-conditioned (comp)
test_that("external path: release-conditioned multinomial (like_type = 2) matches manual nLL", {

  d   <- make_data(like_type = 2, n_tags_released = 50)
  out <- run_fit(d)

  n_rel <- 50 + d$addtotag
  expected <- array(0, dim = dim(out))
  for (ry in 1:2) {
    obs_r  <- c(d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag

    obs_p  <- obs_r / n_rel
    pred_p <- pred_r / n_rel
    tmp_obs  <- c(obs_p,  1 - sum(obs_p))
    tmp_pred <- c(pred_p, 1 - sum(pred_p))

    expected[ry, 1, 1, 1, 1] <- -n_rel * sum(tmp_obs * log(tmp_pred))
  }

  expect_equal(as.numeric(out[, , , 1, 1]), as.numeric(expected[, , , 1, 1]), tolerance = 1e-6)
})

## Type 3: Multinomial, recapture-conditioned (comp)
test_that("external path: recapture-conditioned multinomial (like_type = 3) matches manual nLL", {

  d   <- make_data(like_type = 3)
  out <- run_fit(d)

  expected <- array(0, dim = dim(out))
  for (ry in 1:2) {
    obs_r  <- c(d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag

    n_recap <- sum(obs_r)
    tmp_obs_all  <- obs_r / n_recap
    tmp_pred_all <- pred_r / sum(pred_r)

    expected[ry, 1, 1, 1, 1] <- -n_recap * sum(tmp_obs_all * log(tmp_pred_all))
  }

  expect_equal(as.numeric(out[, , , 1, 1]), as.numeric(expected[, , , 1, 1]), tolerance = 1e-6)
})

## Type 4: Dirichlet-multinomial, release-conditioned (comp)
test_that("external path: release-conditioned Dirichlet-multinomial (like_type = 4) matches direct ddirmult() call", {

  skip_if_not(exists("ddirmult", mode = "function"),
              "ddirmult() wrapper not found in package namespace")

  d   <- make_data(like_type = 4, ln_theta = log(8), n_tags_released = 50)
  out <- run_fit(d)

  n_rel <- 50 + d$addtotag
  expected <- array(0, dim = dim(out))
  for (ry in 1:2) {
    obs_r  <- c(d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag

    obs_p  <- obs_r / n_rel
    pred_p <- pred_r / n_rel
    tmp_obs  <- c(obs_p,  1 - sum(obs_p))
    tmp_pred <- c(pred_p, 1 - sum(pred_p))

    expected[ry, 1, 1, 1, 1] <- -1 * ddirmult(
      obs = tmp_obs, pred = tmp_pred, Ntotal = n_rel,
      ln_theta = d$ln_conv_fish_tag_theta, TRUE
    )
  }

  expect_equal(as.numeric(out[, , , 1, 1]), as.numeric(expected[, , , 1, 1]), tolerance = 1e-6)
})

## Type 5: Dirichlet-multinomial, recapture-conditioned (comp)
test_that("external path: recapture-conditioned Dirichlet-multinomial (like_type = 5) matches direct ddirmult() call", {

  skip_if_not(exists("ddirmult", mode = "function"),
              "ddirmult() wrapper not found in package namespace")

  d   <- make_data(like_type = 5, ln_theta = log(8))
  out <- run_fit(d)

  expected <- array(0, dim = dim(out))
  for (ry in 1:2) {
    obs_r  <- c(d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag

    n_recap <- sum(obs_r)
    tmp_obs_all  <- obs_r / n_recap
    tmp_pred_all <- pred_r / sum(pred_r)

    expected[ry, 1, 1, 1, 1] <- -1 * ddirmult(
      obs = tmp_obs_all, pred = tmp_pred_all, Ntotal = n_recap,
      ln_theta = d$ln_conv_fish_tag_theta, TRUE
    )
  }

  expect_equal(as.numeric(out[, , , 1, 1]), as.numeric(expected[, , , 1, 1]), tolerance = 1e-6)
})

## Helper: tag_fam_of()
test_that("tag_fam_of() classifies likelihood codes correctly", {
  expect_equal(tag_fam_of(0), "count")
  expect_equal(tag_fam_of(1), "count")
  expect_equal(tag_fam_of(2), "comp")
  expect_equal(tag_fam_of(3), "comp")
  expect_equal(tag_fam_of(4), "comp")
  expect_equal(tag_fam_of(5), "comp")
  expect_true(is.na(tag_fam_of(6)))
  expect_true(is.na(tag_fam_of(-1)))
})

## Helper: tag_grid() -- event enumeration and mixing-period skip logic
test_that("tag_grid() enumerates recovery events and correctly applies the mixing-period skip", {

  conv_tag_release_indicator <- matrix(c(1, 1, 1), nrow = 1) # region, year, season

  # mixing_period = 1 -> both ry = 1 and ry = 2 events kept
  g_keep_all <- tag_grid(
    n_conv_tag_cohorts = 1, conv_tag_release_indicator = conv_tag_release_indicator,
    conv_tag_max_liberty = 2, n_yrs = 2, n_seas = 1, conv_tag_mixing_period = 1
  )
  expect_equal(nrow(g_keep_all), 2)
  expect_equal(g_keep_all$ry, c(1, 2))
  expect_equal(g_keep_all$rseas, c(1, 1))

  # mixing_period = 2 -> the ry = 1 event (total_seas_at_liberty = 1) is skipped,
  # only ry = 2 (total_seas_at_liberty = 2) survives
  g_skip_first <- tag_grid(
    n_conv_tag_cohorts = 1, conv_tag_release_indicator = conv_tag_release_indicator,
    conv_tag_max_liberty = 2, n_yrs = 2, n_seas = 1, conv_tag_mixing_period = 2
  )
  expect_equal(nrow(g_skip_first), 1)
  expect_equal(g_skip_first$ry, 2)

  # mixing_period so large that no events survive -> empty data.frame
  g_none <- tag_grid(
    n_conv_tag_cohorts = 1, conv_tag_release_indicator = conv_tag_release_indicator,
    conv_tag_max_liberty = 2, n_yrs = 2, n_seas = 1, conv_tag_mixing_period = 99
  )
  expect_equal(nrow(g_none), NULL)
})


do_pack <- function(d) {
  pack_tag_osa(
    family                     = tag_fam_of(d$conv_fish_tag_like),
    like_type                  = d$conv_fish_tag_like,
    obs_recap                  = d$obs_conv_tag_fish_recap,
    pred_recap                 = d$pred_conv_tag_fish_recap,
    tagged_fish                = d$conv_tagged_fish,
    conv_tag_release_indicator = d$conv_tag_release_indicator,
    conv_tag_max_liberty       = d$conv_tag_max_liberty,
    n_conv_tag_cohorts         = d$n_conv_tag_cohorts,
    n_yrs                      = d$n_yrs,
    n_seas                     = d$n_seas,
    n_regions                  = d$n_regions,
    n_fish_fleets                = d$n_fish_fleets,
    n_pop_pool                 = d$n_conv_tag_pop_pool,
    n_age_pool                 = d$n_conv_tag_age_pool,
    n_sex_pool                 = d$n_conv_tag_sex_pool,
    pop_pool                   = d$conv_tag_pop_pool,
    age_pool                   = d$conv_tag_age_pool,
    sex_pool                   = d$conv_tag_sex_pool,
    use_fish_tagging           = d$use_conv_fish_tagging,
    conv_tag_mixing_period     = d$conv_tag_mixing_period,
    addtotag                   = d$addtotag
  )
}

# Runs pack_tag_osa() -> RTMB::OBS() -> eval_tag_osa() inside a minimal AD
# graph and returns the fitted RTMB object plus the packed vector.
run_internal <- function(d) {
  pack   <- do_pack(d)
  family <- tag_fam_of(d$conv_fish_tag_like)

  f <- function(par) {
    tracked <- RTMB::OBS(pack$vec)
    nLL_arr <- eval_tag_osa(
      nLL_arr                    = d$conv_fish_tag_nLL,
      tracked                    = tracked,
      family                     = family,
      like_type                  = d$conv_fish_tag_like,
      pred_recap                 = d$pred_conv_tag_fish_recap,
      tagged_fish                = d$conv_tagged_fish,
      conv_tag_release_indicator = d$conv_tag_release_indicator,
      conv_tag_max_liberty       = d$conv_tag_max_liberty,
      n_conv_tag_cohorts         = d$n_conv_tag_cohorts,
      n_yrs                      = d$n_yrs,
      obs_recap                   = d$obs_recap,
      n_seas                     = d$n_seas,
      n_regions                  = d$n_regions,
      n_fish_fleets                = d$n_fish_fleets,
      n_pop_pool                 = d$n_conv_tag_pop_pool,
      n_age_pool                 = d$n_conv_tag_age_pool,
      n_sex_pool                 = d$n_conv_tag_sex_pool,
      pop_pool                   = d$conv_tag_pop_pool,
      age_pool                   = d$conv_tag_age_pool,
      sex_pool                   = d$conv_tag_sex_pool,
      use_fish_tagging           = d$use_conv_fish_tagging,
      conv_tag_mixing_period     = d$conv_tag_mixing_period,
      addtotag                   = d$addtotag,
      ln_theta                   = d$ln_conv_fish_tag_theta
    )
    RTMB::REPORT(nLL_arr)
    sum(nLL_arr)
  }

  obj <- RTMB::MakeADFun(f, parameters = list(dummy_par = 0), silent = TRUE)
  list(obj = obj, pack = pack, family = family)
}

## Type 0: Poisson (count)
test_that("internal path: pack_tag_osa()/eval_tag_osa() Poisson (like_type = 0) matches discrete dpois nLL", {

  skip_if_not_installed("RTMB")

  d <- make_data(like_type = 0)
  pack <- do_pack(d)

  # 2 events x 2 regions = 4 packed count cells, event order then region order
  expect_equal(pack$vec, round(c(3, 5, 2, 4) + d$addtotag))
  expect_length(pack$lengths, 2)   # placeholder per event for "count" family
  expect_length(pack$grp_end, 0)

  res <- run_internal(d)
  total_nLL <- res$obj$fn(res$obj$par)

  preds  <- c(2.5, 4.5, 1.8, 3.9) + d$addtotag
  obscnt <- round(c(3, 5, 2, 4) + d$addtotag)
  manual_nLL <- -sum(dpois(obscnt, preds, log = TRUE))

  expect_equal(total_nLL, manual_nLL, tolerance = 1e-6)
})

## Type 1: Negative binomial (count)
test_that("internal path: pack_tag_osa()/eval_tag_osa() negative binomial (like_type = 1) matches discrete NB2 nLL", {

  skip_if_not_installed("RTMB")

  d <- make_data(like_type = 1, ln_theta = log(5))
  pack <- do_pack(d)
  expect_equal(pack$vec, round(c(3, 5, 2, 4) + d$addtotag))

  res <- run_internal(d)
  total_nLL <- res$obj$fn(res$obj$par)

  theta  <- exp(d$ln_conv_fish_tag_theta)
  mus    <- c(2.5, 4.5, 1.8, 3.9) + d$addtotag
  obscnt <- round(c(3, 5, 2, 4) + d$addtotag)
  manual_nLL <- -sum(dnbinom(obscnt, size = theta, mu = mus, log = TRUE))

  expect_equal(total_nLL, manual_nLL, tolerance = 1e-6)
})

## Type 2: Multinomial, release-conditioned (comp)
test_that("internal path: pack_tag_osa()/eval_tag_osa() release-conditioned multinomial (like_type = 2) matches discrete dmultinom nLL", {

  skip_if_not_installed("RTMB")

  d <- make_data(like_type = 2, n_tags_released = 50)
  pack <- do_pack(d)

  n_rel <- 50 + d$addtotag

  # 2 events, each with 2 recap cells + 1 tail cell -> group length 3 each
  expect_length(pack$lengths, 2)
  expect_equal(pack$lengths, c(3, 3))
  expect_equal(pack$grp_end, c(3, 6))
  expect_length(pack$vec, 6)

  # manually reconstruct expected packed counts per event
  expected_vec <- numeric(0)
  for (ry in 1:2) {
    obs_r <- c(d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
               d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    prop <- obs_r / n_rel
    tail <- max(1 - sum(prop), 0)
    prop <- c(prop, tail)
    prop <- prop / sum(prop)
    expected_vec <- c(expected_vec, round(prop * n_rel))
  }
  expect_equal(pack$vec, expected_vec)

  res <- run_internal(d)
  total_nLL <- res$obj$fn(res$obj$par)

  # manual: discrete dmultinom on each event's packed group vs. renormalized
  # predicted proportions (with the same non-recap tail)
  manual_nLL <- 0
  for (ry in 1:2) {
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pprop <- pred_r / n_rel
    pprop <- c(pprop, 1 - sum(pprop))
    pprop <- pprop / sum(pprop)

    idx <- if (ry == 1) 1:3 else 4:6
    manual_nLL <- manual_nLL - dmultinom(pack$vec[idx], prob = pprop, log = TRUE)
  }

  expect_equal(total_nLL, manual_nLL, tolerance = 1e-6)
})

## Type 3: Multinomial, recapture-conditioned (comp)
test_that("internal path: pack_tag_osa()/eval_tag_osa() recapture-conditioned multinomial (like_type = 3) matches discrete dmultinom nLL", {

  skip_if_not_installed("RTMB")

  d <- make_data(like_type = 3)
  pack <- do_pack(d)

  # 2 events, 2 recap cells each, no tail -> group length 2 each
  expect_equal(pack$lengths, c(2, 2))
  expect_equal(pack$grp_end, c(2, 4))
  expect_length(pack$vec, 4)

  res <- run_internal(d)
  total_nLL <- res$obj$fn(res$obj$par)

  manual_nLL <- 0
  for (ry in 1:2) {
    obs_r  <- c(d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$obs_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag

    n_recap <- sum(obs_r)
    g_counts <- round((obs_r / n_recap) * n_recap)
    pprop <- pred_r / sum(pred_r)

    idx <- if (ry == 1) 1:2 else 3:4
    expect_equal(pack$vec[idx], g_counts)
    manual_nLL <- manual_nLL - dmultinom(pack$vec[idx], prob = pprop, log = TRUE)
  }

  expect_equal(total_nLL, manual_nLL, tolerance = 1e-6)
})

## Type 4 & 5: Dirichlet-multinomial (comp)
test_that("internal path: pack_tag_osa()/eval_tag_osa() Dirichlet-multinomial (like_type = 4, release-conditioned) matches direct ddirmult() call", {

  skip_if_not_installed("RTMB")
  skip_if_not_installed("RTMBdist")

  d <- make_data(like_type = 4, ln_theta = log(8), n_tags_released = 50)
  pack <- do_pack(d)
  n_rel <- 50 + d$addtotag

  res <- run_internal(d)
  total_nLL <- res$obj$fn(res$obj$par)

  manual_nLL <- 0
  for (ry in 1:2) {
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pprop <- pred_r / n_rel
    pprop <- c(pprop, 1 - sum(pprop))
    pprop <- pprop / sum(pprop)

    idx <- if (ry == 1) 1:3 else 4:6
    tr  <- pack$vec[idx]
    manual_nLL <- manual_nLL - RTMBdist::ddirmult(
      tr, sum(tr), pprop * exp(d$ln_conv_fish_tag_theta) * sum(tr), log = TRUE
    )
  }

  expect_equal(total_nLL, manual_nLL, tolerance = 1e-6)
})

test_that("internal path: pack_tag_osa()/eval_tag_osa() Dirichlet-multinomial (like_type = 5, recapture-conditioned) matches direct ddirmult() call", {

  skip_if_not_installed("RTMB")
  d <- make_data(like_type = 5, ln_theta = log(8))
  pack <- do_pack(d)

  res <- run_internal(d)
  total_nLL <- res$obj$fn(res$obj$par)

  manual_nLL <- 0
  for (ry in 1:2) {
    pred_r <- c(d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 1, 1, 1, 1],
                d$pred_conv_tag_fish_recap[ry, 1, 1, 1, 2, 1, 1, 1]) + d$addtotag
    pprop <- pred_r / sum(pred_r)

    idx <- if (ry == 1) 1:2 else 3:4
    tr  <- pack$vec[idx]

    # N matches fitting/eval_tag_osa: raw recapture total (not sum(tr) of rounded counts)
    n_recap <- sum(d$obs_conv_tag_fish_recap[ry, 1, 1, 1, , 1, 1, 1] + d$addtotag)

    manual_nLL <- manual_nLL - ddirmult_osa(
      tr, pprop * exp(d$ln_conv_fish_tag_theta) * n_recap, log = TRUE
    )
  }

  expect_equal(total_nLL, manual_nLL, tolerance = 1e-6)
})

## Edge case: release-conditioned tail guard (tail < 0 clamped to 0)
test_that("internal path: pack_tag_osa() clamps a negative release-conditioned tail to zero and renormalizes", {

  d <- make_data(like_type = 2, n_tags_released = 5) # deliberately too few releases
  pack <- do_pack(d)

  n_rel <- 5 + d$addtotag
  for (g in 1:2) {
    idx <- if (g == 1) 1:3 else 4:6
    grp <- pack$vec[idx]
    expect_true(all(grp >= 0))
    # counts in each packed group should sum back to (approximately) n_rel
    expect_equal(sum(grp), round(n_rel), tolerance = 1)
  }
})

