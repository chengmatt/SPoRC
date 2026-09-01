# Closed-loop conditioning against the assessment it is conditioned on.
#
# condition_closed_loop_simulations builds the operating model a management
# strategy evaluation is run through: it takes a fitted assessment, extends every
# input forward over the closed-loop years, and replicates the whole thing across
# simulations. Everything downstream of it inherits whatever it produces.
#
# The one test it had asserted that it returns a list. That is worth having, and
# it is the whole of what was checked, over about twelve hundred lines. A silent
# error here does not crash: it produces management advice computed from a stock
# that is not the one that was assessed.
#
# What can be checked without knowing the right answer is that conditioning
# preserves what it was given. The historical period of the operating model must
# be the assessment, the projected years must extend it by the stated rule, and
# the replicates must agree wherever conditioning is deterministic.

data("dusky_rtmb_model")

closed_loop_cfg <- list(closed_loop_yrs = 5, n_sims = 3)

conditioned_loop <- local({
  cached <- NULL
  function() {
    if(!is.null(cached)) return(cached)
    m <- dusky_rtmb_model
    cached <<- condition_closed_loop_simulations(
      closed_loop_yrs = closed_loop_cfg$closed_loop_yrs, n_sims = closed_loop_cfg$n_sims,
      m$data, m$parameters, m$mapping,
      sd_rep = m$sdrep, rep = m$rep, random = NULL)
    cached
  }
})

#' One simulation's slice of a conditioned array over chosen years
#'
#' Conditioned arrays carry year on the third margin and simulation on the last,
#' whatever sits between, so both are addressed by position rather than by
#' spelling out each array's shape.
#'
#' @keywords internal
cl_slice <- function(arr, yrs, sim) {
  idx <- rep(list(bquote()), length(dim(arr)))
  idx[[3]] <- yrs
  idx[[length(dim(arr))]] <- sim
  do.call(`[`, c(list(arr), idx, list(drop = TRUE)))
}

n_hist <- function() length(dusky_rtmb_model$data$years)


test_that("conditioning extends the time series by the closed-loop years", {
  cl <- conditioned_loop()
  want <- n_hist() + closed_loop_cfg$closed_loop_yrs

  # every conditioned array carries the extended length on its year margin
  for(nm in c("WAA", "MatAA", "natmort", "R0", "h", "sexratio")) {
    expect_equal(dim(cl[[nm]])[3], want, label = sprintf("year extent of %s", nm))
  }
  # movement carries a second region margin, so its year sits one place later
  expect_equal(dim(cl$Movement)[4], want)
  # and the simulation count on its last
  for(nm in c("WAA", "MatAA", "natmort")) {
    expect_equal(dim(cl[[nm]])[length(dim(cl[[nm]]))], closed_loop_cfg$n_sims,
                 label = sprintf("simulation extent of %s", nm))
  }
})


test_that("the historical period of the operating model is the assessment", {
  # The years the assessment covers are not projected, they are carried over. If
  # conditioning perturbs them, the operating model is a different stock from the
  # one that was fitted and every result computed from it is about that other
  # stock.
  cl <- conditioned_loop()
  m <- dusky_rtmb_model
  hist <- seq_len(n_hist())

  for(nm in c("WAA", "MatAA", "WAA_fish", "WAA_srv")) {
    expect_equal(as.vector(cl_slice(cl[[nm]], hist, 1)), as.vector(m$data[[nm]]),
                 tolerance = 1e-12, label = sprintf("%s over the historical period", nm))
  }
  # natural mortality is estimated rather than supplied, so the assessment's own
  # reported values are what conditioning has to reproduce
  expect_equal(as.vector(cl_slice(cl$natmort, hist, 1)), as.vector(m$rep$natmort),
               tolerance = 1e-12)
})


test_that("projected years carry the terminal year forward", {
  # The stated rule for extending an input past the assessment is to hold the
  # last year. A projection that instead repeated the first year, or averaged, or
  # left zeros, changes what the management procedure is tested against.
  cl <- conditioned_loop()
  NH <- n_hist(); CL <- closed_loop_cfg$closed_loop_yrs

  for(nm in c("WAA", "MatAA", "natmort", "WAA_fish")) {
    expect_equal(as.vector(cl_slice(cl[[nm]], (NH + 1):(NH + CL), 1)),
                 as.vector(cl_slice(cl[[nm]], rep(NH, CL), 1)), tolerance = 1e-12,
                 label = sprintf("%s in the projected years", nm))
  }
})


test_that("replicates share the conditioning they are not meant to differ in", {
  # Simulations differ in the data they draw, not in the biology they are handed.
  # A replicate whose weight at age or mortality differs from another's is being
  # given a different operating model, and results across simulations are then
  # not comparable.
  cl <- conditioned_loop()
  hist <- seq_len(n_hist())

  for(s in 2:closed_loop_cfg$n_sims) {
    for(nm in c("WAA", "MatAA", "natmort")) {
      expect_equal(as.vector(cl_slice(cl[[nm]], hist, s)),
                   as.vector(cl_slice(cl[[nm]], hist, 1)), tolerance = 1e-12,
                   label = sprintf("%s, simulation %d against simulation 1", nm, s))
    }
  }
})


test_that("these checks are reading arrays that actually vary", {
  # Every comparison above is an equality. If the arrays were constant the tests
  # would hold against almost any conditioning, so at least one has to carry real
  # structure across ages and years.
  cl <- conditioned_loop()
  waa <- cl_slice(cl$WAA, seq_len(n_hist()), 1)

  expect_gt(diff(range(waa)) / max(abs(waa)), 0.5)
  expect_gt(length(unique(as.vector(round(waa, 8)))), 5)
})


# The checks above run on the packaged dusky model, which carries one population,
# one region, one sex, and no at-age observations. Those are the paths the single
# existing test covered too, so the spatial, sexed and at-age conditioning was
# reached by nothing at all.

cl_dims <- list(
  list(nr = 1, nx = 1, label = "1 region, 1 sex"),
  list(nr = 2, nx = 1, label = "2 regions, 1 sex"),
  list(nr = 1, nx = 2, label = "1 region, 2 sexes"),
  list(nr = 2, nx = 2, label = "2 regions, 2 sexes")
)

#' Condition a closed loop on a freshly fitted model of chosen dimensions
#'
#' @keywords internal
condition_at_dims <- function(nr, nx, closed_loop_yrs = 3, n_sims = 2, np = 1) {
  m <- fitted_small_model(nr = nr, nx = nx, np = np)
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs, n_sims = n_sims,
    m$il$data, m$il$par, m$il$map, sd_rep = m$sdrep, rep = m$fit$rep, random = NULL)
  list(cl = cl, il = m$il, fit = m$fit, n_hist = length(m$il$data$years),
       closed_loop_yrs = closed_loop_yrs)
}


test_that("conditioning works across regions and sexes", {
  # The one model that ships fitted is single-region and single-sex, so every
  # spatial and sexed path through conditioning was previously unreached.
  for(d in cl_dims) {
    out <- expect_no_error(condition_at_dims(d$nr, d$nx))
    want <- out$n_hist + out$closed_loop_yrs

    expect_equal(dim(out$cl$WAA)[3], want, label = sprintf("year extent, %s", d$label))
    expect_equal(dim(out$cl$WAA)[2], d$nr, label = sprintf("region extent, %s", d$label))
    expect_equal(dim(out$cl$WAA)[6], d$nx, label = sprintf("sex extent, %s", d$label))
  }
})


test_that("the historical period survives conditioning at every dimension", {
  # The same claim the dusky checks make, asked of the spatial and sexed paths.
  for(d in cl_dims) {
    out <- condition_at_dims(d$nr, d$nx)
    hist <- seq_len(out$n_hist)

    for(nm in c("WAA", "MatAA")) {
      expect_equal(as.vector(cl_slice(out$cl[[nm]], hist, 1)),
                   as.vector(out$il$data[[nm]]), tolerance = 1e-12,
                   label = sprintf("%s over history, %s", nm, d$label))
    }
  }
})


test_that("the at-age observation streams are extended over the closed-loop years", {
  # Every other year-dimensioned input went through extend_years and these did
  # not, so conditioning failed outright for any model carrying at-age data. It
  # went unnoticed because the only model it was ever run on has none: the
  # failure is a shape error, not a wrong number, so it could not have been
  # hiding in a result.
  #
  # These arrays are region by year by season by age by sex by fleet, with no
  # simulation margin, so the year is the second.
  out <- condition_at_dims(2, 2)
  want <- out$n_hist + out$closed_loop_yrs

  for(nm in c("UseCatchAA", "UseDiscardAA", "UseSrvIdxAA",
              "ObsCatchAA_SE", "ObsDiscardAA_SE", "ObsSrvIdxAA_SE")) {
    arr <- out$cl[[nm]]
    expect_false(is.null(arr), label = sprintf("%s is present after conditioning", nm))
    if(is.null(arr)) next
    expect_equal(dim(arr)[2], want, label = sprintf("year extent of %s", nm))
  }
})


# The checks above run with the at-age observation flags all at zero, which is
# what every fixture in the package happens to carry. At zero the rule for
# extending them past the assessment cannot be observed at all: holding the last
# year and filling with zeros produce the same array. These switch them on.

aa_dims <- list(
  list(nr = 1, nx = 1, label = "1 region, 1 sex"),
  list(nr = 2, nx = 1, label = "2 regions, 1 sex"),
  list(nr = 2, nx = 2, label = "2 regions, 2 sexes")
)

#' Condition and run a closed loop on a model carrying at-age catch
#'
#' @keywords internal
run_at_age_loop <- function(nr, nx, closed_loop_yrs = 3, n_sims = 2, seed = 11) {
  m <- fitted_at_age_model(nr = nr, nx = nx)
  set.seed(seed)
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = closed_loop_yrs, n_sims = n_sims,
    m$il$data, m$il$par, m$il$map, sd_rep = m$sdrep, rep = m$fit$rep, random = NULL)
  sim_env <- Setup_sim_env(cl)
  for(s in seq_len(sim_env$n_sims)) {
    for(y in seq_len(sim_env$n_yrs)) run_annual_cycle(y, s, sim_env)
  }
  list(cl = cl, sim_env = sim_env, m = m, closed_loop_yrs = closed_loop_yrs)
}


test_that("at-age observation flags are carried forward, not zeroed", {
  # This is the rule the extension was written to follow, and it is only visible
  # on a model that observes at age: the projected years keep observing whatever
  # the terminal year observed.
  for(d in aa_dims) {
    out <- run_at_age_loop(d$nr, d$nx)
    NY <- out$m$n_yrs; CL <- out$closed_loop_yrs
    u <- out$cl$UseCatchAA

    per_year <- d$nr * out$m$n_ages * d$nx
    expect_equal(sum(u[, seq_len(NY), , , , , drop = FALSE]), per_year * NY,
                 label = sprintf("historical at-age flags, %s", d$label))
    expect_equal(sum(u[, (NY + 1):(NY + CL), , , , , drop = FALSE]), per_year * CL,
                 label = sprintf("projected at-age flags, %s", d$label))
  }
})


test_that("the operating model draws at-age observations across regions and sexes", {
  # Conditioning producing the right shapes is one claim; the simulation actually
  # generating observations into them is another, and it is the one that says the
  # at-age path runs at these dimensions rather than merely being allocated.
  for(d in aa_dims) {
    out <- run_at_age_loop(d$nr, d$nx)
    NY <- out$m$n_yrs
    drawn <- out$sim_env$ObsCatchAA[, seq_len(NY), , , , , 1, drop = FALSE]

    expect_equal(sum(is.finite(drawn) & drawn != 0), d$nr * NY * out$m$n_ages * d$nx,
                 label = sprintf("at-age observations drawn over the burn-in, %s", d$label))
  }
})


test_that("the projection years are unfished until a management procedure sets them", {
  # run_annual_cycle alone carries the population forward; the fishing mortality
  # for the years after feedback_start_yr comes from the management procedure the
  # caller supplies in the loop around it. Recording that here keeps the empty
  # projected observations above from being read as a missing draw.
  #
  # feedback_start_yr is itself still fished at the conditioned rate: it is the
  # last year carried over, and the procedure sets the year after it onwards.
  out <- run_at_age_loop(1, 1)
  fmort <- out$sim_env$Fmort[1, , 1, 1, 1]
  start <- out$sim_env$feedback_start_yr

  expect_true(all(fmort[seq_len(start)] > 0),
              label = "conditioned years, up to and including feedback_start_yr, are fished")
  expect_true(all(fmort[(start + 1):length(fmort)] == 0),
              label = "years after feedback_start_yr are unfished without a management procedure")
})


test_that("conditioning works for more than one population", {
  # Natal homing had no execution anywhere: every packaged dataset is one
  # population, and the only worked example has eval = FALSE on the chunk that
  # fits. Conditioning a closed loop on a multi-population model is therefore
  # reached here for the first time.
  for(cfg in list(list(np = 2, nr = 1), list(np = 2, nr = 2))) {
    out <- expect_no_error(condition_at_dims(cfg$nr, 1, np = cfg$np))
    want <- out$n_hist + out$closed_loop_yrs

    expect_equal(dim(out$cl$WAA)[1], cfg$np,
                 label = sprintf("population extent, %d populations", cfg$np))
    expect_equal(dim(out$cl$WAA)[3], want,
                 label = sprintf("year extent, %d populations", cfg$np))
  }
})


test_that("a multi-population closed loop runs its annual cycle", {
  # Conditioning producing the right shapes is one claim; the population dynamics
  # stepping forward through every year and replicate is the one that says the
  # loop is usable.
  m <- fitted_small_model(nr = 2, nx = 1, np = 2)
  cl <- condition_closed_loop_simulations(
    closed_loop_yrs = 3, n_sims = 2, m$il$data, m$il$par, m$il$map,
    sd_rep = m$sdrep, rep = m$fit$rep, random = NULL)
  sim_env <- Setup_sim_env(cl)

  for(s in seq_len(sim_env$n_sims)) {
    for(y in seq_len(sim_env$n_yrs)) {
      expect_no_error(run_annual_cycle(y, s, sim_env),
                      message = sprintf("annual cycle, simulation %d year %d", s, y))
    }
  }
})


# State-space numbers at age through the closed loop.
#
# The state is generated by the simulator, not by the estimation model, and the closed loop does
# not use the same driver the open loop does: it calls Setup_sim_env once and then drives
# run_annual_cycle itself. Anything allocated or drawn inside Simulate_Pop_Static is therefore
# invisible to it. These check the two halves of that: that the drivers agree, and that a loop
# conditioned on a fitted state-space model projects with the process error it estimated rather
# than silently with none.

naa_cl_fit <- local({
  cached <- NULL
  function() {
    if(!is.null(cached)) return(cached)
    om <- naaom_make_om(NAA_re = "2dar1", sigmaNAA = 0.30, rho_age = 0.5, rho_year = 0.4, seed = 808)
    il <- naaom_build_em(naaom_om_data(om), NAA_re = "2dar1")
    fit <- suppressWarnings(fit_model(il$data, il$par, il$map, random = "ln_NAA", silent = TRUE))
    cached <<- list(il = il, fit = fit, sdr = suppressWarnings(RTMB::sdreport(fit)))
    cached
  }
})

naa_cl_condition <- function(proj_yrs = 6, n_sims = 2) {
  d <- naa_cl_fit()
  condition_closed_loop_simulations(closed_loop_yrs = proj_yrs, n_sims = n_sims,
                                    data = d$il$data, parameters = d$il$par, mapping = d$il$map,
                                    sd_rep = d$sdr, rep = d$fit$rep, random = "ln_NAA")
}

test_that("both simulation drivers produce the same operating model", {
  # Simulate_Pop_Static runs its own loop; the closed loop drives run_annual_cycle directly. The
  # innovations are drawn once per replicate in the annual cycle so that both paths reach them,
  # and at the same point in the random number stream so the two remain comparable.
  b <- body(naaom_make_om); b[[length(b)]] <- quote(sim_list)
  make_list <- naaom_make_om; body(make_list) <- b
  sl <- make_list(NAA_re = "2dar1", sigmaNAA = 0.25)

  set.seed(1); static <- Simulate_Pop_Static(sim_list = sl)
  set.seed(1)
  sim_env <- Setup_sim_env(sl)
  for(y in 1:sim_env$n_yrs) run_annual_cycle(y = y, sim = 1, sim_env = sim_env)

  expect_equal(sim_env$naa_eta_all, static$naa_eta)
  expect_equal(sim_env$SSB, static$SSB)
  expect_gt(stats::sd(static$naa_eta), 0)
})

test_that("conditioning carries the fitted state across and extends it over the projection", {
  d <- naa_cl_fit()
  cl <- naa_cl_condition(proj_yrs = 6)
  n_hist <- length(d$il$data$years)

  expect_equal(cl$NAA_re, d$il$data$NAA_re)
  expect_equal(cl$n_yrs, n_hist + 6)

  # the blocked standard deviations are gathered through their blocks the way the objective
  # gathers them, then held at the terminal year over the projection
  expect_equal(dim(cl$sigmaNAA)[3], cl$n_yrs)
  fitted_sd <- exp(naaom_fixed(d$fit, "ln_sigmaNAA"))
  expect_equal(as.vector(cl$sigmaNAA[1, 1, n_hist, 1, 1]), as.vector(fitted_sd), tolerance = 1e-8)
  expect_equal(cl$sigmaNAA[1, 1, cl$n_yrs, 1, 1], cl$sigmaNAA[1, 1, n_hist, 1, 1])

  # the correlations come back on the natural scale, not the unconstrained one
  expect_true(all(abs(cl$naa_rho) < 1))
  expect_gt(cl$naa_rho[["age"]], 0)
  expect_gt(cl$naa_rho[["year"]], 0)

  # the active years run on through the projection: the operating model keeps generating process
  # error in the future rather than stopping where the data stop
  expect_equal(max(cl$naa_re_yrs), cl$n_yrs)
  expect_equal(min(cl$naa_re_yrs), min(d$il$data$naa_re_yrs))
  expect_equal(cl$naa_re_ages, d$il$data$naa_re_ages)
})

test_that("the conditioned operating model generates process error over the projection years", {
  d <- naa_cl_fit()
  cl <- naa_cl_condition(proj_yrs = 6)
  n_hist <- length(d$il$data$years)

  set.seed(2)
  sim_env <- Setup_sim_env(cl)
  for(y in 1:sim_env$n_yrs) run_annual_cycle(y = y, sim = 1, sim_env = sim_env)

  eta <- sim_env$naa_eta_all[1, 1, , , 1, 1]
  ages <- cl$naa_re_ages
  proj <- eta[(n_hist + 1):cl$n_yrs, ages]
  expect_true(all(is.finite(proj)))
  expect_gt(stats::sd(proj), 0)
  # and at the same magnitude as the historical period rather than tapering off
  expect_equal(stats::sd(proj), stats::sd(eta[2:n_hist, ages]), tolerance = 0.4)
  # year one and age one stay deterministic, in the projection as in the past
  expect_equal(max(abs(eta[, 1])), 0)
  expect_equal(max(abs(eta[1, ])), 0)
})

test_that("a fit without the state conditions to an operating model without it", {
  # The failure this guards is silent: before the state was carried across, every closed loop
  # projected deterministically whatever the assessment had estimated.
  cl <- conditioned_loop()
  expect_equal(cl$NAA_re, 0)
  expect_equal(cl$sigmaNAA, 0)
  expect_equal(cl$naa_re_yrs, integer(0))

  sim_env <- Setup_sim_env(cl)
  expect_equal(dim(sim_env$naa_eta_all)[3], cl$n_yrs)
  run_annual_cycle(y = 1, sim = 1, sim_env = sim_env)
  expect_equal(max(abs(sim_env$naa_eta_all)), 0)
})
