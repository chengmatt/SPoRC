library(SPoRC)
library(testthat)

# The growth module against the SS3 conventions it reproduces, then against the
# data-driven model it replaces: with parameters chosen to regenerate the size-age
# transition a fixture was built from, the objective and gradient have to match
# the model that read that transition as data.

L0 <- 9; L1 <- 20; L2 <- 55; K <- 0.25; CV1 <- 0.12; CV2 <- 0.08; A1 <- 2; A2 <- 20

# The CAAL fixture's full input list, with the growth module's fields and parameter
# taken from a biologicals setup that switched growth on. Everything the fixture
# set up downstream of the biologicals (comps, CAAL, selectivity, weights) is kept.
merge_growth <- function(rest, growth_input) {
  keep <- grep("^growth_|^derive_waa$|^wt_len_pars$|^SizeAgeTrans$", names(growth_input$data), value = TRUE)
  for(nm in keep) rest$data[[nm]] <- growth_input$data[[nm]]
  rest$par$ln_growth_pars <- growth_input$par$ln_growth_pars
  rest$map$ln_growth_pars <- growth_input$map$ln_growth_pars
  rest
}


test_that("the curve is linear below A1, continuous at A1, and hits L2 at A2", {
  x <- seq(0, 25, by = 0.5)
  crv <- get_laa_curve(x, L0, L1, L2, K, CV1, CV2, A1, A2)
  # linear phase from L0 at age zero to L1 at A1
  expect_equal(crv$L[x == 0], L0)
  expect_equal(crv$L[x == 1], L0 + 0.5 * (L1 - L0))
  expect_equal(crv$L[x == A1], L1)
  # von Bertalanffy beyond A1, hitting L2 exactly at A2
  expect_equal(crv$L[x == A2], L2)
  linf <- L1 + (L2 - L1) / (1 - exp(-K * (A2 - A1)))
  expect_equal(crv$Linf, linf)
  expect_equal(crv$L[x == 10], linf + (L1 - linf) * exp(-K * (10 - A1)))
  # monotone throughout
  expect_true(all(diff(crv$L) > 0))
})


test_that("L2_asymptote reads L2 as the asymptote and A1 = 0 puts age zero at L1", {
  x <- 0:30
  crv <- get_laa_curve(x, L0, L1, L2, K, CV1, CV2, A1, A2 = max(x), L2_asymptote = 1)
  expect_equal(crv$Linf, L2)
  expect_lt(crv$L[31], L2)
  crv0 <- get_laa_curve(x, L0, L1, L2, K, CV1, CV2, A1 = 0, A2)
  expect_equal(crv0$L[1], L1)
  expect_true(all(is.finite(crv0$L)))
})


test_that("the CV is CV1 below A1, CV2 from A2, and interpolates on length between", {
  x <- seq(0, 25, by = 0.5)
  crv <- get_laa_curve(x, L0, L1, L2, K, CV1, CV2, A1, A2)
  expect_equal(crv$cv[x < A1], rep(CV1, sum(x < A1)))
  expect_equal(crv$cv[x >= A2], rep(CV2, sum(x >= A2)))
  mid <- x >= A1 & x < A2
  expect_equal(crv$cv[mid], CV1 + (crv$L[mid] - L1) * (CV2 - CV1) / (L2 - L1))
  expect_equal(crv$sd, crv$cv * crv$L)
  # on age instead
  crv_a <- get_laa_curve(x, L0, L1, L2, K, CV1, CV2, A1, A2, cv_type = 1)
  expect_equal(crv_a$cv[mid], CV1 + (x[mid] - A1) * (CV2 - CV1) / (A2 - A1))
  # as standard deviations
  crv_s <- get_laa_curve(x, L0, L1, L2, K, CV1, CV2, A1, A2, sd_type = 1)
  expect_equal(crv_s$sd, crv_s$cv)
})


test_that("the plus group adjustment is the decay-weighted mean SS3 uses", {
  L_acc <- 50; linf <- 60; n_acc <- 20
  a <- 0:n_acc; w <- exp(-0.2 * a)
  expect_equal(plus_group_size(L_acc, linf, n_acc), sum(w * (L_acc + a / n_acc * (linf - L_acc))) / sum(w))
  # it sits between the curve at the accumulator age and the asymptote
  expect_gt(plus_group_size(L_acc, linf, n_acc), L_acc)
  expect_lt(plus_group_size(L_acc, linf, n_acc), linf)
})


test_that("the age-length key is column stochastic with the tails in the end bins", {
  lower <- seq(9, 65, by = 2)
  x <- 0:20
  crv <- get_laa_curve(x, L0, L1, L2, K, CV1, CV2, A1, A2)
  alk <- get_alk(lower, crv$L, crv$sd)
  expect_equal(dim(alk), c(length(lower), length(x)))
  expect_equal(colSums(alk), rep(1, length(x)))
  expect_true(all(alk >= 0))
  # a tiny fish at age zero: almost all of its mass is below the first edge and lands in bin 1
  expect_gt(alk[1, 1], 0.5)
  # direct check on one column against the CDF construction
  a <- 10
  cdf <- pnorm((lower - crv$L[a]) / crv$sd[a])
  expect_equal(alk[, a], c(cdf[2:length(lower)], 1) - cdf + c(cdf[1], rep(0, length(lower) - 1)))
  # lognormal variant is also column stochastic
  alk_ln <- get_alk(lower, crv$L, crv$cv, dist = 1)
  expect_equal(colSums(alk_ln), rep(1, length(x)))
})


test_that("with parameters that regenerate the fixture's key, the objective matches the data-driven model", {
  source(test_path("helper-selftest_caal.R"), local = TRUE)
  om <- caal_make_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)
  n_ages <- caal_cfg$n_ages

  # data-driven model: the fixture's key built from VB(linf 60, k 0.3, t0 -0.5), CV 0.10,
  # evaluated at integer ages at the start of the year
  inp_data <- caal_build_input(sim_data)
  fit_data <- fit_model(inp_data$data, inp_data$par, inp_data$map, random = NULL, silent = TRUE, do_optim = FALSE)

  # the same curve in Schnute form with A1 the first age and A2 the last, so the
  # linear phase never applies and L2 is the curve's own value at A2
  g <- caal_growth(n_ages, 1)[[1]]
  pars <- c(L1 = g$mean[1], L2 = g$mean[n_ages], K = 0.3, CV1 = 0.10, CV2 = 0.10)

  build_growth <- function(spec = "fix") {
    input <- Setup_Mod_Dim(years = 1:caal_cfg$n_yrs, ages = 1:n_ages, lens = caal_cfg$len_lower + 2.5,
                           n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                           n_pop = 1, natal_region = 1, verbose = FALSE)
    input <- Setup_Mod_Rec(input_list = input, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                           ln_sigmaR = array(log(0.3), c(2, 1, 1)), rec_model = "mean_rec",
                           sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2, ln_global_R0 = log(5))
    input <- Setup_Mod_Biologicals(input_list = input, WAA = sim_data$WAA, MatAA = sim_data$MatAA,
                                   WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv,
                                   fit_lengths = 1, SizeAgeTrans = NA, comp_const_obs = 0,
                                   AgeingError = sim_data$AgeingError, M_spec = "fix",
                                   Fixed_natmort = array(0.3, dim = c(1, 1, caal_cfg$n_yrs, n_ages, 1)),
                                   growth_model = "vb_schnute", growth_spec = spec,
                                   ln_growth_pars = array(log(pars), dim = c(1, 1, 1, 5)),
                                   growth_A1 = 1, growth_A2 = n_ages, growth_len_lower = caal_cfg$len_lower,
                                   growth_plus_group = "curve")
    # the rest of the fixture, with only the growth fields taken from this input
    merge_growth(caal_build_input(sim_data), input)
  }

  inp_g <- build_growth("fix")
  expect_equal(inp_g$data$growth_model, 1)
  expect_true(all(is.na(inp_g$map$ln_growth_pars)))
  fit_g <- fit_model(inp_g$data, inp_g$par, inp_g$map, random = NULL, silent = TRUE, do_optim = FALSE)

  # the derived key regenerates the fixture's key
  # one key per fleet, both read at the season start here, so both are the fixture's
  expect_equal(as.vector(fit_g$rep$SizeAgeTrans_fish[,,,,,,,1]), as.vector(sim_data$SizeAgeTrans), tolerance = 1e-10)
  expect_equal(as.vector(fit_g$rep$SizeAgeTrans_srv[,,,,,,,1]), as.vector(sim_data$SizeAgeTrans), tolerance = 1e-10)
  expect_equal(fit_g$rep$mean_LAA_srv[1, 1, 1, 1, , 1, 1], g$mean, tolerance = 1e-10)
  # and so the objective and gradient are those of the data-driven model
  expect_equal(fit_g$fn(fit_g$par), fit_data$fn(fit_data$par), tolerance = 1e-10)
  expect_equal(as.numeric(fit_g$gr(fit_g$par)), as.numeric(fit_data$gr(fit_data$par)), tolerance = 1e-8)

  # estimating growth: the gradient is finite and the map has one entry per parameter
  inp_e <- build_growth("est_all")
  expect_equal(sum(!is.na(inp_e$map$ln_growth_pars)), 5)
  fit_e <- fit_model(inp_e$data, inp_e$par, inp_e$map, random = NULL, silent = TRUE, do_optim = FALSE)
  expect_true(all(is.finite(fit_e$gr(fit_e$par))))
  expect_true("ln_growth_pars" %in% names(fit_e$par))
})


test_that("growth parameters are recovered from lengths plus conditional age-at-length", {
  source(test_path("helper-selftest_caal.R"), local = TRUE)
  om <- caal_make_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)
  n_ages <- caal_cfg$n_ages
  g <- caal_growth(n_ages, 1)[[1]]
  truth <- c(L1 = g$mean[1], L2 = g$mean[n_ages], K = 0.3, CV1 = 0.10, CV2 = 0.10)

  input <- Setup_Mod_Dim(years = 1:caal_cfg$n_yrs, ages = 1:n_ages, lens = caal_cfg$len_lower + 2.5,
                         n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                         n_pop = 1, natal_region = 1, verbose = FALSE)
  input <- Setup_Mod_Rec(input_list = input, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                         ln_sigmaR = array(log(0.3), c(2, 1, 1)), rec_model = "mean_rec",
                         sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2, ln_global_R0 = log(5))
  # start away from the truth
  input <- Setup_Mod_Biologicals(input_list = input, WAA = sim_data$WAA, MatAA = sim_data$MatAA,
                                 WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv,
                                 fit_lengths = 1, SizeAgeTrans = NA, comp_const_obs = 0,
                                 AgeingError = sim_data$AgeingError, M_spec = "fix",
                                 Fixed_natmort = array(0.3, dim = c(1, 1, caal_cfg$n_yrs, n_ages, 1)),
                                 growth_model = "vb_schnute",
                                 ln_growth_pars = array(log(truth * c(1.2, 0.9, 0.7, 1.5, 0.6)), dim = c(1, 1, 1, 5)),
                                 growth_spec = "est_all", growth_A1 = 1, growth_A2 = n_ages,
                                 growth_len_lower = caal_cfg$len_lower, growth_plus_group = "curve")
  input <- merge_growth(caal_build_input(sim_data), input)

  fit <- fit_model(input$data, input$par, input$map, random = NULL, silent = TRUE)
  expect_lt(max(abs(fit$gr(fit$env$last.par.best))), 1e-2)
  est <- exp(fit$env$last.par.best[names(fit$env$last.par.best) == "ln_growth_pars"])
  names(est) <- names(truth)
  # one replicate, 30 years of lengths and CAAL: the mean lengths and K come back
  # within a few percent, the CVs within ten
  expect_lt(max(abs(est[c("L1", "L2", "K")] / truth[c("L1", "L2", "K")] - 1)), 0.05)
  expect_lt(max(abs(est[c("CV1", "CV2")] / truth[c("CV1", "CV2")] - 1)), 0.10)
})


test_that("each fleet's key and weight are read at that fleet's own timing", {
  ages <- 0:6; lower <- seq(10, 60, by = 5)
  pars <- array(log(c(15, 55, 0.3, 0.15, 0.08)), dim = c(1, 1, 1, 5))
  wl <- array(c(1e-5, 3), dim = c(1, 1, 1, 2))
  run <- function(t_fish, t_srv) {
    Get_Growth(ln_growth_pars = pars, growth_A1 = 1, growth_A2 = 6, growth_L0 = 10, growth_len_lower = lower,
               growth_cv_type = 0, growth_sd_type = 0, growth_dist = 0, growth_plus_group = 0,
               derive_waa = 1, wt_len_pars = wl, ages = ages, seasdur = 1, spawn_seas = 1, t_spawn = 0,
               n_pop = 1, n_regions = 1, n_yrs = 2, n_seas = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 2,
               t_fish = array(t_fish, dim = c(1, 1, 1)), t_srv = array(t_srv, dim = c(1, 1, 2)))
  }
  g <- run(0.5, c(0.25, 0.5))

  # the second survey is taken when the fishery is, the first is not
  expect_equal(g$SizeAgeTrans_srv[,,,,,,,2], g$SizeAgeTrans_fish[,,,,,,,1])
  expect_equal(g$WAA_srv[1,1,1,1,,1,2], g$WAA_fish[1,1,1,1,,1,1])
  expect_false(isTRUE(all.equal(g$WAA_srv[1,1,1,1,,1,1], g$WAA_fish[1,1,1,1,,1,1])))

  # a fleet's weight is its own key times weight at the bin midpoints
  w_key <- as.vector(t(g$SizeAgeTrans_srv[1,1,1,1,,,1,1]) %*% (1e-5 * (lower + 2.5)^3))
  expect_equal(g$WAA_srv[1,1,1,1,,1,1], w_key, tolerance = 1e-12)

  # the mean length a fleet's key is built on is the curve at its real age
  crv <- get_laa_curve(ages + 0.25, L0 = 10, L1 = 15, L2 = 55, K = 0.3, CV1 = 0.15, CV2 = 0.08, A1 = 1, A2 = 6)
  expect_equal(g$mean_LAA_srv[1,1,1,1,,1,1], crv$L, tolerance = 1e-12)
  expect_equal(g$mean_LAA_spawn[1,1,1,1,,1], get_laa_curve(ages, L0 = 10, L1 = 15, L2 = 55, K = 0.3, CV1 = 0.15, CV2 = 0.08, A1 = 1, A2 = 6)$L, tolerance = 1e-12)
})


test_that("a model estimating growth runs through the simulation self-test", {
  # The operating model the self-test builds takes its size-age transition from
  # the fitted model's report when growth is estimated, since the data list holds
  # none. Three replicates are enough to show the loop runs and that the refits
  # land near the truth they were simulated from.
  source(test_path("helper-selftest_caal.R"), local = TRUE)
  om <- caal_make_om()
  sim_data <- simulation_data_to_SPoRC(sim_env = om, y = caal_cfg$n_yrs, sim = 1)
  n_ages <- caal_cfg$n_ages
  g <- caal_growth(n_ages, 1)[[1]]
  truth <- c(L1 = g$mean[1], L2 = g$mean[n_ages], K = 0.3, CV1 = 0.10, CV2 = 0.10)

  input <- Setup_Mod_Dim(years = 1:caal_cfg$n_yrs, ages = 1:n_ages, lens = caal_cfg$len_lower + 2.5,
                         n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1,
                         n_pop = 1, natal_region = 1, verbose = FALSE)
  input <- Setup_Mod_Rec(input_list = input, do_rec_bias_ramp = 0, sigmaR_switch = 1,
                         ln_sigmaR = array(log(0.3), c(2, 1, 1)), rec_model = "mean_rec",
                         sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = 2, ln_global_R0 = log(5))
  input <- Setup_Mod_Biologicals(input_list = input, WAA = sim_data$WAA, MatAA = sim_data$MatAA,
                                 WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv,
                                 fit_lengths = 1, SizeAgeTrans = NA, comp_const_obs = 0,
                                 AgeingError = sim_data$AgeingError, M_spec = "fix",
                                 Fixed_natmort = array(0.3, dim = c(1, 1, caal_cfg$n_yrs, n_ages, 1)),
                                 growth_model = "vb_schnute",
                                 ln_growth_pars = array(log(truth), dim = c(1, 1, 1, 5)),
                                 growth_spec = "est_all", growth_A1 = 1, growth_A2 = n_ages,
                                 growth_len_lower = caal_cfg$len_lower, growth_plus_group = "curve")
  input <- merge_growth(caal_build_input(sim_data), input)

  fit <- fit_model(input$data, input$par, input$map, random = NULL, silent = TRUE)
  expect_false(is.null(fit$rep$SizeAgeTrans_srv))

  sd_rep <- RTMB::sdreport(fit)
  st <- simulation_self_test(data = input$data, parameters = input$par, mapping = input$map,
                             random = NULL, rep = fit$rep, sd_rep = sd_rep, n_sims = 3,
                             what = c("SSB", "mean_LAA_srv"))

  # every replicate refit: the replicates are stacked in the last dimension and a
  # failed one would be NA
  expect_false(anyNA(st$SSB))

  # spawning biomass comes back around the truth it was simulated from
  ssb_true <- as.vector(fit$rep$SSB)[1:caal_cfg$n_yrs]
  rel <- apply(st$SSB[1, 1, , ], 2, function(x) median(abs(x / ssb_true - 1)))
  expect_lt(median(rel), 0.10)

  # and so does mean length at age at the composition time
  laa_true <- fit$rep$mean_LAA_srv[1, 1, 1, 1, , 1, 1]
  laa_err <- apply(st$mean_LAA_srv[1, 1, 1, 1, , 1, 1, ], 2, function(x) max(abs(x / laa_true - 1)))
  expect_lt(median(laa_err), 0.03)
})


test_that("the selection-weighted weight at age is the key averaged by selectivity", {

  # three length bins and two ages, the first age mostly short fish and the
  # second mostly long, so the two ages overlap in the bins they occupy
  key <- cbind(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6)) # [n_lens x n_ages], columns sum to one
  w_len <- c(0.5, 1.0, 2.0)
  pop <- as.vector(t(key) %*% w_len) # the population mean weight at age

  # flat selectivity takes a random draw from the age, so the selected weight is
  # the population weight at age whatever height the flat value sits at
  for(val in c(0.05, 1, 7)) expect_equal(get_selected_waa(key, rep(val, 3), w_len), pop, tolerance = 1e-12)

  # knife edge on one bin: every fish taken weighs that bin's weight, at every age
  for(l in 1:3) expect_equal(get_selected_waa(key, replace(rep(0, 3), l, 1), w_len), rep(w_len[l], 2), tolerance = 1e-12)

  # a gear that favors the long fish takes heavier fish than the age average and
  # one that favors the short takes lighter, both inside the bins the age occupies
  big <- get_selected_waa(key, c(0.1, 0.5, 1), w_len)
  small <- get_selected_waa(key, c(1, 0.5, 0.1), w_len)
  expect_true(all(big > pop))
  expect_true(all(small < pop))
  expect_true(all(big <= max(w_len)))
  expect_true(all(small >= min(w_len)))

  # against the definition written out age by age
  sel <- c(0.1, 0.5, 1)
  hand <- sapply(1:2, function(a) sum(key[, a] * sel * w_len) / sum(key[, a] * sel))
  expect_equal(get_selected_waa(key, sel, w_len), hand, tolerance = 1e-12)

  # an age nothing is selected from returns zero rather than dividing by zero
  expect_equal(get_selected_waa(key, rep(0, 3), w_len), c(0, 0))
})


test_that("only the flagged fleets take the selection-weighted weight, in the year asked for", {

  n_lens <- 3
  n_ages <- 2
  n_fish <- 2
  n_yrs <- 2
  key <- cbind(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))
  len_mid <- c(10, 20, 30)
  wl <- c(1e-5, 3)
  sel <- c(0.1, 0.5, 1)

  # -1 stands in for whatever the population weight at age was, so anything the
  # call leaves alone is visible
  WAA_fish <- array(-1, dim = c(1, 1, n_yrs, 1, n_ages, 1, n_fish))
  SizeAgeTrans_fish <- array(0, dim = c(1, 1, n_yrs, 1, n_lens, n_ages, 1, n_fish))
  fish_sel_l <- array(0, dim = c(1, n_yrs, n_lens, 1, n_fish))
  for(y in 1:n_yrs) {
    for(f in 1:n_fish) {
      SizeAgeTrans_fish[1, 1, y, 1, , , 1, f] <- key
      fish_sel_l[1, y, , 1, f] <- sel
    } # end f loop
  } # end y loop
  wt_len_pars <- array(wl, dim = c(1, 1, 1, 2))

  out <- growth_selected_waa_year(WAA_fish, SizeAgeTrans_fish, fish_sel_l, wt_len_pars, len_mid,
                                  waa_selected = c(1, 0), y = 2, n_pop = 1, n_regions = 1,
                                  n_seas = 1, n_sexes = 1)

  w_mid <- wl[1] * len_mid^wl[2]
  expect_equal(out[1, 1, 2, 1, , 1, 1], get_selected_waa(key, sel, w_mid), tolerance = 1e-12)
  # the unflagged fleet keeps the population weight, as does the year not asked for
  expect_equal(out[1, 1, 2, 1, , 1, 2], rep(-1, n_ages))
  expect_equal(out[1, 1, 1, 1, , 1, 1], rep(-1, n_ages))
  # and this gear takes the long fish, so the weight it lands is the heavier one
  expect_true(all(out[1, 1, 2, 1, , 1, 1] > as.vector(t(key) %*% w_mid)))
})
