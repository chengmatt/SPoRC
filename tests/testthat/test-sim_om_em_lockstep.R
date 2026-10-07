# The operating model a self test builds reproduces the estimation model at the fit, across option configurations and
# case-study models: every true state and fitted observation is the prediction (1e-8), compositions are drawn from the
# expectation the model builds, and perfect data from a jittered start return a model whose design identifies it.

library(SPoRC)
library(testthat)

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# configurations spanning the observation sources, their population and season forms, and the population processes
lockstep_configs <- list(
  sweep = function() list(il = quiet(sweep_input())),
  populations = function() list(il = quiet(pop_sources_input(n_pop = 3, n_regions = 3))),
  populations_seasons = function() list(il = quiet(pop_sources_input(n_seas = 2))),
  mvn_index = function() list(il = quiet(pop_sources_input(like = "mvn"))),
  normal_index = function() list(il = quiet(pop_sources_input(like = "normal"))),
  year_totals = function() list(il = quiet(year_total_input())),
  at_age_by_population = function() list(il = quiet(pop_seas_at_age_input())$input),
  lengths_discards = function() list(il = quiet(objective_setup_input(n_lens = 10))),
  age_at_length = function() list(il = quiet(caal_build_input(simulation_data_to_SPoRC(caal_make_om(), caal_cfg$n_yrs, 1)))),
  tags_movement = function() list(il = quiet(build_em(quiet(build_om(move_timing = 0)), move_timing = 0))),
  numbers_at_age_state = function() list(il = quiet(naaom_build_em(naaom_om_data(naaom_make_om(NAA_re = "iid")), NAA_re = "iid")), random = "ln_NAA"),
  goa_dusky = function() list(il = quiet(build_goa_dusky_input(sgl_rg_dusky_data))),
  bsai_northern = function() list(il = quiet(seed_bsai_nork_mle(build_bsai_nork_input(sgl_rg_bsai_nork_data), sgl_rg_bsai_nork_data))),
  goa_rex = function() list(il = quiet(seed_goa_rex_mle(build_goa_rex_input(mlt_rg_goa_rex_data), mlt_rg_goa_rex_data))),
  ebs_pcod = function() list(il = quiet(seed_ebs_pcod_mle(build_ebs_pcod_input(sgl_rg_ebs_pcod_data), sgl_rg_ebs_pcod_data)))
)

test_that("the operating model's true state and every fitted observation are the estimation model's prediction", {

  for(name in names(lockstep_configs)) {
    cfg <- lockstep_configs[[name]]()
    diffs <- lockstep_diffs(lockstep_om(cfg$il, random = cfg$random))
    diffs <- diffs[!is.na(diffs)]
    expect_true(length(diffs) > 8, info = name) # the state and at least one data source were compared
    expect_lt(max(diffs), 1e-8, label = paste(name, names(diffs)[which.max(diffs)]))
  } # end name loop
})

test_that("each multinomial composition is drawn from the expectation the estimation model builds", {

  # (bins - 1) / 2 is the expected nLL of an exact draw at a sample size of 1e6; another expectation gives thousands
  for(name in c("sweep", "populations", "year_totals", "lengths_discards", "goa_dusky")) {
    cfg <- lockstep_configs[[name]]()
    comps <- lockstep_comp_nll(lockstep_om(cfg$il, random = cfg$random))
    expect_true(nrow(comps) > 0, info = name)
    for(i in seq_len(nrow(comps))) expect_lt(comps$nll_per_comp[i], 3 * comps$reference[i] + 2, label = paste(name, comps$source[i]))
  } # end name loop
})

test_that("perfect data return a model whose design identifies it, from a jittered start", {

  for(name in c("populations", "age_at_length", "tags_movement", "goa_dusky", "bsai_northern")) {

    cfg <- lockstep_configs[[name]]()
    il <- cfg$il
    obj <- quiet(fit_model(il$data, il$par, il$map, random = cfg$random, do_optim = FALSE, silent = TRUE))
    full <- obj$env$last.par.best

    # start every freely estimated fixed effect off the truth, a shared level moving together
    set.seed(42)
    start <- il$par
    for(par_name in names(start)) {
      levels <- if(is.null(il$map[[par_name]])) seq_along(start[[par_name]]) else as.integer(il$map[[par_name]])
      if(all(is.na(levels))) next
      shift <- stats::rnorm(max(levels, na.rm = TRUE), 0, 0.1)
      free <- !is.na(levels)
      start[[par_name]][free] <- start[[par_name]][free] + shift[levels[free]]
    } # end par_name loop

    res <- quiet(simulation_self_test(data = obj$data, parameters = start, mapping = il$map, random = cfg$random,
                                      rep = obj$report(full), sd_rep = list(par.fixed = full, par.random = NULL),
                                      n_sims = 1, newton_loops = 3, perfect_data = TRUE, what = c("SSB", "Rec", "tot_FAA"),
                                      what_par = "ln_sigmaC"))
    # the catch was drawn at the perfect sd, which is the truth and the value the refit uses
    expect_equal(as.vector(res$truth$ln_sigmaC), rep(log(1e-3), length(res$truth$ln_sigmaC)), info = name)
    expect_equal(as.vector(res$ln_sigmaC), as.vector(res$truth$ln_sigmaC), info = name)
    # the last years' recruits are seen by no data yet, so recruitment is checked on its median alone
    for(quantity in c("SSB", "Rec", "tot_FAA")) {
      err <- abs(as.vector(res[[quantity]]) / as.vector(res$truth[[quantity]]) - 1)
      err <- err[is.finite(err)]
      expect_lt(stats::median(err), 0.01, label = paste(name, quantity, "median"))
      if(quantity != "Rec") expect_lt(max(err), 0.05, label = paste(name, quantity, "max"))
    } # end quantity loop

  } # end name loop
})

test_that("numbers at age drawn fresh are drawn at the fitted sd, and each replicate is compared with its draws", {

  # the operating model the self test builds, with the innovations left to it to draw as joint does
  cfg <- lockstep_configs$numbers_at_age_state()
  run <- lockstep_om(cfg$il, random = cfg$random, perfect_data = FALSE, n_sims = 40)
  run$sim_list$naa_eta_input <- NULL
  run$om <- suppressWarnings(suppressMessages(Simulate_Pop_Static(run$sim_list)))

  # every active cell drawn at the fit's process sd, rather than kept at the fit's own innovations
  sigma <- unique(signif(run$sim_list$sigmaNAA[run$sim_list$sigmaNAA > 0], 8)) # one process sd in this model
  expect_length(sigma, 1)
  eta <- run$om$naa_eta[,,run$sim_list$naa_re_yrs,1,run$sim_list$naa_re_ages,,]
  expect_equal(stats::sd(as.vector(eta)), sigma, tolerance = 0.05)

  # the truth each replicate is compared with is what it ran on, population by region by year
  fit_ssb <- run$rep$SSB
  for(sim in 1:3) {
    truth <- om_truth(run$om$SSB, fit_ssb, sim, 40)
    expect_equal(as.vector(truth), as.vector(run$om$SSB[,,,sim]))
    expect_false(isTRUE(all.equal(as.vector(truth), as.vector(fit_ssb))))
  } # end sim loop

  # a fit's array a year longer takes the operating model's values in the years they share
  n_yrs <- dim(fit_ssb)[3]
  longer <- array(c(fit_ssb, 99), dim = c(1, 1, n_yrs + 1))
  truth <- om_truth(run$om$SSB, longer, 2, 40)
  expect_equal(as.vector(truth), c(run$om$SSB[1,1,,2], 99))

  # deviations the fit stops before the last year take the operating model's leading years
  shorter <- array(0, dim = c(1, 1, n_yrs - 2))
  truth <- om_truth(run$om$SSB, shorter, 3, 40)
  expect_equal(as.vector(truth), run$om$SSB[1,1,seq_len(n_yrs - 2),3])
})

test_that("recruitment drawn from the penalty is what the estimation model computes from those deviations", {

  # each penalty form, the bias ramp, the own-mean center, a walk with an unpenalized first year, and region sharing
  rec_forms <- list(
    ramp = list(sigmaR_spec = "est_all", do_rec_bias_ramp = 1, bias_year = c(2, 5, 9, 12)),
    own_mean = list(sigmaR_spec = "est_all", RecDevs_pen_center = "own_mean"),
    walk_unpenalized_first = list(RecDevs_model = "rw", sigmaR_spec = "est_all", dont_pen_recdev_first = 1, ln_global_R0_spec = "fix"),
    ar1 = list(RecDevs_model = "ar1", sigmaR_spec = "est_all", RecDevs_rho_spec = "est_all"),
    shared_regions = list(sigmaR_spec = "est_all", RecDevs_spec = "est_shared_r")
  )

  set.seed(1)
  for(name in names(rec_forms)) {

    il <- quiet(sweep_input(rec = rec_forms[[name]]))
    il$par$ln_sigmaR[] <- log(0.6)
    if(!is.null(il$par$RecDevs_rho)) il$par$RecDevs_rho[] <- 0.5
    il$par$ln_RecDevs[] <- stats::rnorm(length(il$par$ln_RecDevs), 0, 0.3) # a fitted series to draw on from
    run <- lockstep_om(il, draw_rec = TRUE)
    om_devs <- array(run$om$ln_RecDevs[,,seq_len(dim(il$par$ln_RecDevs)[3]),1], dim = dim(il$par$ln_RecDevs))

    # every penalized deviation is redrawn, one draw per shared level, and an unpenalized one keeps the fit's value
    pen_map <- array(il$data$map_ln_RecDevs, dim = dim(om_devs))
    redrawn <- abs(om_devs - il$par$ln_RecDevs) > 1e-12
    expect_true(all(redrawn[!is.na(pen_map)]), info = name)
    expect_false(any(redrawn[is.na(pen_map)]), info = name)
    for(y in seq_len(dim(om_devs)[3])) {
      penalized <- !is.na(pen_map[,,y])
      expect_equal(length(unique(round(om_devs[,,y][penalized], 12))), length(unique(pen_map[,,y][penalized])), info = paste(name, y))
    } # end y loop

    # the estimation model at the drawn deviations, the initial ages drawn with them, reproduces the operating model
    il$par$ln_RecDevs[] <- om_devs
    il$par$ln_InitDevs[] <- array(run$om$ln_InitDevs[,,,,1], dim = dim(il$par$ln_InitDevs))
    obj <- quiet(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
    run$rep <- obj$report(obj$env$last.par.best)
    diffs <- lockstep_diffs(run)
    expect_true(sum(!is.na(diffs)) > 8, info = name)
    expect_lt(max(diffs, na.rm = TRUE), 1e-8, label = paste(name, names(diffs)[which.max(diffs)]))

  } # end name loop
})

test_that("drawn recruitment deviations follow the penalty the estimation model puts on them", {

  # standardized by that penalty, each form's draws are independent with mean zero and sd one; the independent
  # form carries a bias ramp, an early and late sigma, and a likelihood weight that narrows the draw
  rec_forms <- list(
    iid = list(rec = list(sigmaR_spec = "est_all", do_rec_bias_ramp = 1, bias_year = c(2, 5, 9, 12), sigmaR_switch = 6), wt = list(Wt_Rec = 0.5)),
    walk = list(rec = list(RecDevs_model = "rw", sigmaR_spec = "est_all", RecDevs_rw_init_sigma = NA), wt = list()),
    ar1 = list(rec = list(RecDevs_model = "ar1", sigmaR_spec = "est_all", RecDevs_rho_spec = "est_all"), wt = list())
  )

  set.seed(2)
  n_draws <- 600
  for(name in names(rec_forms)) {

    il <- quiet(sweep_input(dims = list(n_regions = 1), rec = rec_forms[[name]]$rec, wt = rec_forms[[name]]$wt))
    il$par$ln_sigmaR[1,,] <- log(0.4)
    il$par$ln_sigmaR[2,,] <- log(0.7)
    if(name == "ar1") il$par$RecDevs_rho[] <- 0.6 # every parameter list has the correlation, read only under an ar1
    obj <- quiet(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
    rep <- obj$report(obj$env$last.par.best)

    n_yrs <- dim(il$par$ln_RecDevs)[3]
    sigma <- ifelse(seq_len(n_yrs) < obj$data$sigmaR_switch, 0.4, 0.7)
    rho <- if(name == "ar1") rho_trans(il$par$RecDevs_rho[1]) else 0
    center <- -sigma^2 / (2 * (1 - rho^2)) * rep$bias_ramp # minus half the variance, stationary under an ar1

    z <- matrix(NA, n_draws, n_yrs)
    for(k in seq_len(n_draws)) {
      devs <- rec_devs_past_fit(obj$data, il$par, rep, 0, n_yrs)[1,1,]
      if(name == "iid") z[k,] <- (devs - center) / (sigma / sqrt(obj$data$Wt_Rec))
      if(name == "walk") z[k,] <- (devs - c(0, devs[-n_yrs])) / sigma
      if(name == "ar1") z[k,] <- c((devs[1] - center[1]) * sqrt(1 - rho^2) / sigma[1],
                                   (devs[-1] - center[-1] - rho * (devs[-n_yrs] - center[-n_yrs])) / sigma[-1])
    } # end k loop

    expect_lt(abs(mean(z)), 0.05, label = paste(name, "mean"))
    expect_lt(abs(stats::sd(as.vector(z)) - 1), 0.04, label = paste(name, "sd"))
    expect_lt(abs(stats::cor(as.vector(z[,-1]), as.vector(z[,-n_yrs]))), 0.05, label = paste(name, "lag one correlation"))

  } # end name loop
})

test_that("initial ages drawn from the penalty are what the estimation model computes from them", {

  # the default curve shared by sex, the bias ramp, a curve per sex tied by its penalty, and an ar1
  init_forms <- list(
    shared_sexes = list(sigmaR_spec = "est_all"),
    ramp = list(sigmaR_spec = "est_all", do_rec_bias_ramp = 1, bias_year = c(2, 5, 9, 12)),
    sex_tie = list(sigmaR_spec = "est_all", InitDevs_sex_spec = "est_all", Use_init_sex_pen = 1),
    ar1 = list(sigmaR_spec = "est_all", RecDevs_model = "ar1", RecDevs_rho_spec = "est_all")
  )

  set.seed(1)
  for(name in names(init_forms)) {

    il <- quiet(sweep_input(rec = init_forms[[name]]))
    il$par$ln_sigmaR[] <- log(0.6)
    if(name == "ar1") il$par$RecDevs_rho[] <- 0.5
    il$par$ln_InitDevs[] <- stats::rnorm(length(il$par$ln_InitDevs), 0, 0.3) # a fitted curve to draw away from
    run <- lockstep_om(il, draw_rec = TRUE)
    dev_dim <- dim(il$par$ln_InitDevs)
    om_devs <- array(run$om$ln_InitDevs[,,,,1], dim = dev_dim)

    # every penalized cell is redrawn, a fixed one keeps the fit's value, and a shared level takes one draw
    pen_map <- array(run$data$map_ln_InitDevs, dim = dev_dim)
    pen_use <- array(run$data$init_devs_pen_use, dim = dev_dim)
    redrawn <- abs(om_devs - il$par$ln_InitDevs) > 1e-12
    expect_true(all(redrawn[!is.na(pen_map) & pen_use == 1]), info = name)
    expect_false(any(redrawn[is.na(pen_map)]), info = name)
    expect_equal(length(unique(round(om_devs[!is.na(pen_map)], 12))), length(unique(pen_map[!is.na(pen_map)])), info = name)

    # the estimation model at the drawn initial ages and recruitment reproduces the operating model
    il$par$ln_InitDevs[] <- om_devs
    il$par$ln_RecDevs[] <- array(run$om$ln_RecDevs[,,seq_len(dim(il$par$ln_RecDevs)[3]),1], dim = dim(il$par$ln_RecDevs))
    obj <- quiet(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
    run$rep <- obj$report(obj$env$last.par.best)
    diffs <- lockstep_diffs(run)
    expect_lt(max(diffs, na.rm = TRUE), 1e-8, label = paste(name, names(diffs)[which.max(diffs)]))

  } # end name loop
})

test_that("drawn initial ages follow the penalty the estimation model puts on them", {

  # standardized by that penalty the draws have mean zero and sd one, under the bias ramp and, for the second sex,
  # given the first sex's draw through the tie between them
  il_ramp <- quiet(sweep_input(dims = list(n_regions = 1), rec = list(sigmaR_spec = "est_all", do_rec_bias_ramp = 1, bias_year = c(2, 5, 9, 12)),
                               wt = list(Wt_Init_Rec = 0.5)))
  il_tie <- quiet(sweep_input(dims = list(n_regions = 1), rec = list(sigmaR_spec = "est_all", InitDevs_sex_spec = "est_all", Use_init_sex_pen = 1)))

  set.seed(3)
  for(name in c("ramp", "sex_tie")) {

    il <- if(name == "ramp") il_ramp else il_tie
    il$par$ln_sigmaR[1,,] <- log(0.5)
    obj <- quiet(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
    rep <- obj$report(obj$env$last.par.best)
    pen_use <- array(obj$data$init_devs_pen_use, dim = dim(il$par$ln_InitDevs))
    ages <- which(!is.na(array(obj$data$map_ln_InitDevs, dim = dim(il$par$ln_InitDevs))[1,1,,1]) & pen_use[1,1,,1] == 1)
    wt <- array(obj$data$Wt_Init_Rec, dim = dim(il$par$ln_InitDevs))[1,1,ages,1]
    center <- -0.5^2 / 2 * rep$init_bias_ramp[ages]

    z <- NULL
    for(k in 1:500) {
      devs <- init_devs_past_fit(obj$data, il$par, rep)
      first <- (devs[1,1,ages,1] - center) / (0.5 / sqrt(wt))
      if(name == "sex_tie") {
        tie_sd <- exp(obj$data$ln_sigma_init_sex)
        precision <- 1 / 0.5^2 + 1 / tie_sd^2
        cond_mean <- (center / 0.5^2 + devs[1,1,ages,1] / tie_sd^2) / precision
        z <- c(z, (devs[1,1,ages,2] - cond_mean) * sqrt(precision))
      } else z <- c(z, first)
    } # end k loop

    expect_lt(abs(mean(z)), 0.05, label = paste(name, "mean"))
    expect_lt(abs(stats::sd(z) - 1), 0.04, label = paste(name, "sd"))

  } # end name loop
})

test_that("a catch the fit has missing stays missing in the refit, so its fleet keeps fishing", {

  # fleet two fits compositions with no catch recorded, so the model estimates its F rather than closing it
  il <- quiet(sweep_input(dims = list(n_regions = 1, n_fish_fleets = 2, n_sexes = 1)))
  il$data$ObsCatch[,,,2] <- NA
  il$data$UseCatch[,,,2] <- 0
  obj <- quiet(fit_model(il$data, il$par, il$map, do_optim = FALSE, silent = TRUE))
  expect_true(all(obj$rep$Fmort[,,,2] > 0))

  res <- quiet(simulation_self_test(data = obj$data, parameters = il$par, mapping = il$map, random = NULL, rep = obj$rep,
                                    sd_rep = list(par.fixed = obj$par, par.random = NULL), n_sims = 1, newton_loops = 0,
                                    what = c("Fmort", "FishAgeComps_nLL")))
  expect_true(all(res$Fmort[,,,2,1] > 0))
  expect_true(all(is.finite(res$FishAgeComps_nLL)))
})
