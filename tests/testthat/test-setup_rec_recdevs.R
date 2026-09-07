# Process error on the recruitment deviations: independent, random walk, or AR1.
#
# Checks the density against a hand-written walk, the guard rails on the setup arguments, and that a
# data list written before the option existed still gets the independent penalty. Tolerance 1e-10.

library(SPoRC)
library(testthat)

# ── do_RecDevs_rho_mapping ───────────────────────────────────────────────────

test_that("do_RecDevs_rho_mapping only activates RecDevs_rho under RecDevs_model = 'ar1'", {

  make_il <- function(RecDevs_model_code, n_pop = 2, n_regions = 3) {
    messages_list <<- character(0) # collect_message() writes to this global
    list(
      data = list(
        n_pop = n_pop,
        n_regions = n_regions,
        RecDevs_model = RecDevs_model_code
      ),
      par = list(RecDevs_rho = array(0, dim = c(n_pop, n_regions))),
      map = list()
    )
  }

  il <- SPoRC:::do_RecDevs_rho_mapping(make_il(1), "est_all")
  expect_true(all(is.na(il$map$RecDevs_rho)))

  il <- SPoRC:::do_RecDevs_rho_mapping(make_il(2), "est_all")
  expect_true(all(is.na(il$map$RecDevs_rho)))

  il <- SPoRC:::do_RecDevs_rho_mapping(make_il(3), "est_all")
  expect_equal(length(unique(il$map$RecDevs_rho[!is.na(il$map$RecDevs_rho)])), 2 * 3)

  il <- SPoRC:::do_RecDevs_rho_mapping(make_il(3), "est_shared_pop_r")
  expect_equal(length(unique(il$map$RecDevs_rho[!is.na(il$map$RecDevs_rho)])), 1)

  il <- SPoRC:::do_RecDevs_rho_mapping(make_il(3), "est_shared_r")
  expect_equal(length(unique(il$map$RecDevs_rho[!is.na(il$map$RecDevs_rho)])), 2)

  il <- SPoRC:::do_RecDevs_rho_mapping(make_il(3), "fix")
  expect_true(all(is.na(il$map$RecDevs_rho)))
})

# ── get_recdev_pe_nLL ────────────────────────────────────────────────────────

# the walk written out by hand, for one series of deviations
hand_recdev_nLL <- function(devs, is_est, sigma, dev_mu, PE_model, rho = 0, init_sd = 5) {
  out <- rep(0, length(devs))
  for(y in seq_along(devs)) {
    if(is_est[y] == 0) next
    if(PE_model == 1) out[y] <- -dnorm(devs[y], dev_mu[y], sigma[y], log = TRUE)
    if(PE_model == 2) {
      if(y == 1) out[y] <- -dnorm(devs[y], 0, init_sd, log = TRUE)
      else out[y] <- -dnorm(devs[y], devs[y - 1], sigma[y], log = TRUE)
    }
    if(PE_model == 3) {
      if(y == 1) out[y] <- -dnorm(devs[y], 0, sigma[y] / sqrt(1 - rho^2), log = TRUE)
      else out[y] <- -dnorm(devs[y], rho * devs[y - 1], sigma[y], log = TRUE)
    }
  }
  out
}

test_that("get_recdev_pe_nLL matches a hand-written walk for iid, rw and ar1", {

  n_yrs <- 10
  set.seed(11)
  devs <- rnorm(n_yrs, sd = 0.5)

  # years 1, 5 and 6 fixed. the walk passes through them, so year 2 anchors on the fixed
  # year 1 and year 7 on the fixed year 6, and none of the three contributes a density
  is_est <- rep(1, n_yrs)
  is_est[c(1, 5, 6)] <- 0

  # an early and a late regime, so a step reads the sigma of the year it lands on
  sigma <- rep(0.4, n_yrs)
  sigma[1:4] <- 0.9
  dev_mu <- rep(-0.1, n_yrs)

  rho_raw <- 0.4
  rho <- 2 / (1 + exp(-2 * rho_raw)) - 1

  for(PE_model in 1:3) {
    got <- SPoRC:::get_recdev_pe_nLL(devs, is_est, sigma, dev_mu, PE_model, rho = rho)
    want <- hand_recdev_nLL(devs, is_est, sigma, dev_mu, PE_model, rho = rho)
    expect_equal(as.numeric(got), want, tolerance = 1e-10, info = paste("PE_model:", PE_model))
    expect_true(all(got[is_est == 0] == 0), info = paste("fixed years zero, PE_model:", PE_model))
  }
})

test_that("a random walk passes through fixed years rather than stepping over them", {

  n_yrs <- 7
  set.seed(12)
  devs <- rnorm(n_yrs, sd = 0.5)
  is_est <- rep(1, n_yrs)
  is_est[4:6] <- 0 # asserted deviations in the middle of the series
  sigma <- rep(0.35, n_yrs)

  got <- SPoRC:::get_recdev_pe_nLL(devs, is_est, sigma, rep(0, n_yrs), 2)

  # the fixed years contribute nothing themselves
  expect_true(all(got[4:6] == 0))

  # year 7 anchors on the fixed year 6, not on year 3 with an inflated variance
  expect_equal(as.numeric(got[7]), -dnorm(devs[7], devs[6], 0.35, log = TRUE), tolerance = 1e-10)
  expect_equal(as.numeric(got[2]), -dnorm(devs[2], devs[1], 0.35, log = TRUE), tolerance = 1e-10)

  # fixing the terminal deviations removes their terms without pulling the last estimated one
  is_est_term <- c(rep(1, 4), rep(0, 3))
  term <- SPoRC:::get_recdev_pe_nLL(devs, is_est_term, sigma, rep(0, n_yrs), 2)
  expect_true(all(term[5:7] == 0))
  expect_equal(as.numeric(term[4]), -dnorm(devs[4], devs[3], 0.35, log = TRUE), tolerance = 1e-10)
})

test_that("the diffuse start sits on year one and fixing it leaves the level free", {

  n_yrs <- 6
  set.seed(13)
  devs <- rnorm(n_yrs, sd = 0.5)
  sigma <- rep(0.4, n_yrs)
  all_est <- rep(1, n_yrs)

  got <- SPoRC:::get_recdev_pe_nLL(devs, all_est, sigma, rep(0, n_yrs), 2)
  expect_equal(as.numeric(got[1]), -dnorm(devs[1], 0, 5, log = TRUE), tolerance = 1e-10)

  # widening the start weakens that term, and only that term changes
  wide <- SPoRC:::get_recdev_pe_nLL(devs, all_est, sigma, rep(0, n_yrs), 2, init_sd = 500)
  expect_equal(as.numeric(wide[1]), -dnorm(devs[1], 0, 500, log = TRUE), tolerance = 1e-10)
  expect_equal(as.numeric(wide[2:n_yrs]), as.numeric(got[2:n_yrs]), tolerance = 1e-14)

  # NA starts the walk at zero under its own sigma instead
  own <- SPoRC:::get_recdev_pe_nLL(devs, all_est, sigma, rep(0, n_yrs), 2, init_sd = NA)
  expect_equal(as.numeric(own[1]), -dnorm(devs[1], 0, sigma[1], log = TRUE), tolerance = 1e-10)

  # fixing year one removes the only term that reads the level, so every remaining term is a
  # difference and shifting the whole series does not change the penalty. this is SAM's flat prior
  is_est <- c(0, rep(1, n_yrs - 1))
  base <- sum(SPoRC:::get_recdev_pe_nLL(devs, is_est, sigma, rep(0, n_yrs), 2))
  shifted <- sum(SPoRC:::get_recdev_pe_nLL(devs + 2.5, is_est, sigma, rep(0, n_yrs), 2))
  expect_equal(as.numeric(base), as.numeric(shifted), tolerance = 1e-10)
})

# ── get_recruitment_penalty ──────────────────────────────────────────────────

# the arguments the recruitment penalty needs that this test is not varying
rec_pen_args <- function(ln_RecDevs, map_ln_RecDevs, ln_sigmaR, sigmaR_switch, ...) {
  n_pop <- dim(ln_RecDevs)[1]
  n_regions <- dim(ln_RecDevs)[2]
  n_ages <- 5
  list(
    n_pop = n_pop,
    n_regions = n_regions,
    n_ages = n_ages,
    n_est_rec_devs = dim(ln_RecDevs)[3],
    rec_dd = 999,
    natal_region = seq_len(n_pop),
    rec_region_prop_spec = 0,
    rec_region_prop = array(1 / n_regions, dim = c(n_pop, n_regions)),
    equil_init_age_strc = 0, # no initial age penalty, so only Rec_nLL is under test
    ln_InitDevs = array(0, dim = c(n_pop, n_regions, n_ages - 1, 1)),
    init_age_devs_shared = NULL,
    ln_sigmaR = ln_sigmaR,
    bias_ramp = rep(0, dim(ln_RecDevs)[3]),
    sigmaR_switch = sigmaR_switch,
    ln_RecDevs = ln_RecDevs,
    sigmaR2_early = array(exp(ln_sigmaR[1,,])^2, dim = c(n_pop, n_regions)),
    sigmaR2_late = array(exp(ln_sigmaR[2,,])^2, dim = c(n_pop, n_regions)),
    do_rec_bias_ramp = 0,
    map_ln_RecDevs = map_ln_RecDevs,
    ...
  )
}

test_that("get_recruitment_penalty walks the deviations under rw and ar1", {

  n_pop <- 1; n_regions <- 2; n_yrs <- 9
  set.seed(14)
  ln_RecDevs <- array(rnorm(n_pop * n_regions * n_yrs, sd = 0.6), dim = c(n_pop, n_regions, n_yrs))

  # region 2 has its first three deviations fixed, so its walk contributes from year 4 on,
  # anchored on the asserted year 3 value
  map_ln_RecDevs <- array(seq_len(length(ln_RecDevs)), dim = dim(ln_RecDevs))
  map_ln_RecDevs[1, 2, 1:3] <- NA

  sigmaR_switch <- 5 # the early sigma runs to year 4, the late sigma from year 5
  ln_sigmaR <- array(log(c(0.8, 0.35)), dim = c(2, n_pop, n_regions))

  rho_raw <- 0.3
  rho <- 2 / (1 + exp(-2 * rho_raw)) - 1

  for(PE_model in 2:3) {

    got <- do.call(SPoRC:::get_recruitment_penalty, rec_pen_args(
      ln_RecDevs, map_ln_RecDevs, ln_sigmaR, sigmaR_switch,
      RecDevs_model = PE_model,
      RecDevs_rho = array(rho_raw, dim = c(n_pop, n_regions))
    ))$Rec_nLL

    for(r in 1:n_regions) {
      sigma_yr <- rep(exp(ln_sigmaR[2, 1, r]), n_yrs)
      sigma_yr[1:(sigmaR_switch - 1)] <- exp(ln_sigmaR[1, 1, r])
      want <- hand_recdev_nLL(
        devs = ln_RecDevs[1, r, ],
        is_est = as.numeric(!is.na(map_ln_RecDevs[1, r, ])),
        sigma = sigma_yr,
        dev_mu = rep(0, n_yrs),
        PE_model = PE_model,
        rho = rho
      )
      expect_equal(as.numeric(got[1, r, ]), want, tolerance = 1e-10,
                   info = paste("PE_model:", PE_model, "region:", r))
    }
  }
})

test_that("get_recruitment_penalty leaves the independent penalty unchanged", {

  n_pop <- 1; n_regions <- 1; n_yrs <- 8
  set.seed(15)
  ln_RecDevs <- array(rnorm(n_yrs, sd = 0.6), dim = c(n_pop, n_regions, n_yrs))
  map_ln_RecDevs <- array(seq_len(n_yrs), dim = dim(ln_RecDevs))
  ln_sigmaR <- array(log(c(0.7, 0.4)), dim = c(2, n_pop, n_regions))

  args <- rec_pen_args(ln_RecDevs, map_ln_RecDevs, ln_sigmaR, sigmaR_switch = 4)

  # not supplying RecDevs_model at all is the path a data list written before the
  # option existed takes, and it must give the independent penalty
  default <- do.call(SPoRC:::get_recruitment_penalty, args)$Rec_nLL
  explicit <- do.call(SPoRC:::get_recruitment_penalty, c(args, list(RecDevs_model = 1)))$Rec_nLL
  expect_equal(default, explicit, tolerance = 1e-14)

  sigma_yr <- rep(exp(ln_sigmaR[2, 1, 1]), n_yrs)
  sigma_yr[1:3] <- exp(ln_sigmaR[1, 1, 1])
  want <- -dnorm(ln_RecDevs[1, 1, ], 0, sigma_yr, log = TRUE)
  expect_equal(as.numeric(default[1, 1, ]), want, tolerance = 1e-10)
})

test_that("dont_pen_recdev_first frees the level of a walk without breaking the series", {

  n_pop <- 1; n_regions <- 1; n_yrs <- 8
  set.seed(16)
  ln_RecDevs <- array(rnorm(n_yrs, sd = 0.6), dim = c(n_pop, n_regions, n_yrs))
  ln_sigmaR <- array(log(c(0.5, 0.5)), dim = c(2, n_pop, n_regions))

  # do_RecDevs_mapping drops the leading years from the data mirror alone, so they stay
  # estimated and only lose their own penalty term
  map_all <- array(seq_len(n_yrs), dim = dim(ln_RecDevs))
  map_first_out <- map_all
  map_first_out[1, 1, 1] <- NA

  walk_args <- function(map) rec_pen_args(ln_RecDevs, map, ln_sigmaR, sigmaR_switch = 1, RecDevs_model = 2)

  with_first <- do.call(SPoRC:::get_recruitment_penalty, walk_args(map_all))$Rec_nLL
  without <- do.call(SPoRC:::get_recruitment_penalty, walk_args(map_first_out))$Rec_nLL

  # only the first year's term goes; every later step is unchanged, so year two still
  # anchors on the first year rather than restarting the walk
  expect_equal(as.numeric(without[1, 1, 1]), 0)
  expect_equal(as.numeric(without[1, 1, 2:n_yrs]), as.numeric(with_first[1, 1, 2:n_yrs]), tolerance = 1e-14)
  expect_equal(as.numeric(with_first[1, 1, 1]), -dnorm(ln_RecDevs[1, 1, 1], 0, 5, log = TRUE), tolerance = 1e-10)

  # with the first year out, the penalty is a sum of differences and the level is free
  shifted <- array(ln_RecDevs + 3, dim = dim(ln_RecDevs))
  args_shift <- rec_pen_args(shifted, map_first_out, ln_sigmaR, sigmaR_switch = 1, RecDevs_model = 2)
  moved <- do.call(SPoRC:::get_recruitment_penalty, args_shift)$Rec_nLL
  expect_equal(sum(as.numeric(moved)), sum(as.numeric(without)), tolerance = 1e-10)
})

test_that("Setup_Mod_Rec leaves a dont_pen_recdev_first year estimated but unpenalized", {

  il <- suppressMessages(suppressWarnings(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
    rec = list(RecDevs_model = "rw", sigmaR_spec = "est_all", dont_pen_recdev_first = 1,
               ln_global_R0_spec = "fix"), # the only identified way to ask for this combination
    stop_after = "rec"
  )))

  # the parameter map keeps the deviation, the data mirror the penalty reads drops it
  expect_false(is.na(il$map$ln_RecDevs[1]))
  expect_true(is.na(il$data$map_ln_RecDevs[1, 1, 1]))
  expect_false(any(is.na(il$data$map_ln_RecDevs[1, 1, -1])))
})

# ── Setup_Mod_Rec guard rails ────────────────────────────────────────────────

test_that("Setup_Mod_Rec validates RecDevs_model and its combinations", {

  mk <- function(...) suppressMessages(suppressWarnings(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
    rec = list(...),
    stop_after = "rec"
  )))

  expect_error(mk(RecDevs_model = "randomwalk"), "RecDevs_model incorrectly specified")
  expect_error(mk(RecDevs_model = "rw", do_rec_bias_ramp = 1, bias_year = c(2, 4, 6, 8)), "bias ramp")
  expect_error(mk(RecDevs_model = "rw", RecDevs_pen_center = "own_mean"), "own_mean")

  expect_warning(
    suppressMessages(sweep_input(
      dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
      rec = list(RecDevs_model = "rw", sigmaR_spec = "fix"),
      stop_after = "rec"
    )),
    "sigmaR_spec"
  )

  expect_warning(
    suppressMessages(sweep_input(
      dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
      rec = list(RecDevs_model = "ar1", sigmaR_spec = "est_all", RecDevs_rho_spec = "fix"),
      stop_after = "rec"
    )),
    "RecDevs_rho_spec"
  )

  il <- mk(RecDevs_model = "rw", sigmaR_spec = "est_all")
  expect_equal(il$data$RecDevs_model, 2)
  expect_true(all(is.na(il$map$RecDevs_rho)))
  expect_equal(dim(il$par$RecDevs_rho), c(1, 1))

  il <- mk(RecDevs_model = "ar1", sigmaR_spec = "est_all", RecDevs_rho_spec = "est_all")
  expect_equal(il$data$RecDevs_model, 3)
  expect_false(any(is.na(il$map$RecDevs_rho)))
})

# ── equil_init_age_strc = "stoch_all_no_pen" ─────────────────────────────────

test_that("stoch_all_no_pen estimates every initial deviation and penalizes none", {

  mk <- function(...) suppressMessages(suppressWarnings(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
    rec = list(init_age_strc = "free", ...),
    stop_after = "rec"
  )))

  n_dev <- 5 # n_ages - 1

  # "equil" means no estimation and no penalty, which is one statement about an equilibrium
  # age structure and the wrong one about a free age structure
  eq <- mk(equil_init_age_strc = "equil")
  expect_true(all(is.na(eq$map$ln_InitDevs)))

  # "stoch_all_no_pen" estimates the same cells "stoch_all" does
  all_pen <- mk(equil_init_age_strc = "stoch_all")
  no_pen <- mk(equil_init_age_strc = "stoch_all_no_pen")
  expect_equal(no_pen$data$equil_init_age_strc, 4)
  expect_equal(length(unique(no_pen$map$ln_InitDevs[!is.na(no_pen$map$ln_InitDevs)])), n_dev)
  expect_equal(as.numeric(no_pen$map$ln_InitDevs), as.numeric(all_pen$map$ln_InitDevs))

  # the sharing specs reach it the same way
  shared <- mk(equil_init_age_strc = "stoch_all_no_pen", InitDevs_spec = "est_shared_pop_r")
  expect_equal(length(unique(shared$map$ln_InitDevs[!is.na(shared$map$ln_InitDevs)])), n_dev)
})

test_that("the initial age penalty is zero under stoch_all_no_pen and not under stoch_all", {

  mk_rep <- function(setting) {
    il <- suppressMessages(suppressWarnings(sweep_input(
      dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
      rec = list(init_age_strc = "free", equil_init_age_strc = setting, sigmaR_spec = "est_all")
    )))
    il$par$ln_InitDevs[] <- 0.4 # something the penalty would notice
    obj <- fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE)
    obj$report(obj$par)
  }

  expect_equal(sum(mk_rep("stoch_all_no_pen")$Init_Rec_nLL), 0)
  expect_true(sum(mk_rep("stoch_all")$Init_Rec_nLL) != 0)
})

# ── ln_global_R0_spec ────────────────────────────────────────────────────────

test_that("ln_global_R0_spec maps the recruitment level and refuses the flat combination", {

  mk <- function(...) suppressMessages(suppressWarnings(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
    rec = list(...),
    stop_after = "rec"
  )))

  expect_error(mk(ln_global_R0_spec = "nope"), "must be est or fix")

  # the default estimates it, "fix" maps it off
  expect_false(all(is.na(mk()$map$ln_global_R0)))
  expect_true(all(is.na(mk(ln_global_R0_spec = "fix")$map$ln_global_R0)))

  # a walk with the first year out of the penalty leaves the level and the deviations exactly
  # unidentified, so the combination is refused rather than fitted to a singular hessian
  expect_error(
    mk(RecDevs_model = "rw", dont_pen_recdev_first = 1, sigmaR_spec = "est_all"),
    "mutually unidentified"
  )

  # fixing the level is the way out, and is how SAM writes recruitment
  expect_true(all(is.na(
    mk(RecDevs_model = "rw", dont_pen_recdev_first = 1, sigmaR_spec = "est_all",
       ln_global_R0_spec = "fix")$map$ln_global_R0)))

  # with the first year still penalized the level is readable, weakly, so it warns rather than stops
  expect_warning(
    suppressMessages(sweep_input(
      dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
      rec = list(RecDevs_model = "rw", sigmaR_spec = "est_all"),
      stop_after = "rec"
    )),
    "weakly identified"
  )

  # an ar1 reverts toward zero, so the deviations do have a level and the level is estimable
  expect_false(all(is.na(
    mk(RecDevs_model = "ar1", dont_pen_recdev_first = 1, sigmaR_spec = "est_all",
       RecDevs_rho_spec = "est_all")$map$ln_global_R0)))
})

test_that("ln_global_R0 still reaches the dots as a starting value", {

  # the spec sits after the dots on purpose. before them, R would partially match a supplied
  # ln_global_R0 to the longer formal and take the starting value as the spec
  il <- suppressMessages(suppressWarnings(sweep_input(
    dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 10, n_ages = 6),
    rec = list(ln_global_R0 = log(1234)),
    stop_after = "rec"
  )))

  expect_equal(as.numeric(il$par$ln_global_R0), log(1234), tolerance = 1e-12)

  fm <- names(formals(SPoRC::Setup_Mod_Rec))
  expect_gt(which(fm == "ln_global_R0_spec"), which(fm == "..."))
})

# ── End to end ───────────────────────────────────────────────────────────────

test_that("a fitted model reports the random walk recruitment penalty", {

  mk_obj <- function(...) {
    il <- suppressMessages(suppressWarnings(sweep_input(
      dims = list(n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_yrs = 12, n_ages = 6),
      rec = list(...)
    )))
    list(il = il, obj = fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  }

  walk <- mk_obj(RecDevs_model = "rw", sigmaR_spec = "est_all")
  rep_walk <- walk$obj$report(walk$obj$par)

  n_yrs <- length(walk$il$data$years)
  sigma <- exp(walk$il$par$ln_sigmaR[2, 1, 1])
  devs <- walk$il$par$ln_RecDevs[1, 1, ]
  is_est <- as.numeric(!is.na(walk$il$data$map_ln_RecDevs[1, 1, ]))

  want <- hand_recdev_nLL(devs, is_est, rep(sigma, n_yrs), rep(0, n_yrs), 2)
  expect_equal(as.numeric(rep_walk$Rec_nLL[1, 1, ]), want, tolerance = 1e-10)

  # the objective still evaluates and differentiates
  expect_true(is.finite(walk$obj$fn(walk$obj$par)))
  expect_true(all(is.finite(walk$obj$gr(walk$obj$par))))

  # and the independent model is not the same objective
  iid <- mk_obj(RecDevs_model = "iid", sigmaR_spec = "est_all")
  expect_false(isTRUE(all.equal(as.numeric(walk$obj$fn(walk$obj$par)),
                                as.numeric(iid$obj$fn(iid$obj$par)))))
})
