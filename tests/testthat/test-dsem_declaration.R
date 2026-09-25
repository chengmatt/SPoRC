# Checks the "dsem" declaration at the process modules: it holds the sigma nothing reads any more (or refuses an
# explicit estimate), Setup_Mod_DSEM takes the declared processes and insists every series has an arrow, and a
# declaration with no dsem cannot reach the fit.

data("sgl_rg_dusky_data")

dusky_with <- function(rec = list(), biol = list()) {
  build <- build_goa_dusky_input
  subs <- list()
  if(length(rec) > 0) subs$Setup_Mod_Rec <- bquote(function(...) Setup_Mod_Rec(..., ..(rec)), splice = TRUE)
  if(length(biol) > 0) subs$Setup_Mod_Biologicals <- bquote(function(...) Setup_Mod_Biologicals(..., ..(biol)), splice = TRUE)
  body(build) <- do.call(substitute, list(body(build), subs))
  suppressMessages(build(sgl_rg_dusky_data))
}

test_that("RecDevs_model = 'dsem' reads sigmaR off the arrows, refuses an explicit estimate, and needs a dsem to fit", {

  # the sweep model estimates nothing about sigmaR by default and has the ramp off, so the declaration takes
  il <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"), dims = list(n_regions = 1)))
  expect_equal(il$data$dsem_declared, "rec")
  expect_equal(il$data$RecDevs_model, 1) # the iid code
  expect_true(all(is.na(il$map$ln_sigmaR)))

  # an explicit request to estimate sigmaR contradicts the handover, the arrows describe every year, and a random
  # effect takes the full correction or none, so a ramp is refused while a ramp at zero (no correction) is not
  expect_error(suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem", sigmaR_spec = "est_all"), dims = list(n_regions = 1))), "cannot be estimated")
  expect_error(suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem", dont_est_recdev_last = 1), dims = list(n_regions = 1))), "dont_est_recdev_last")
  expect_error(suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem", do_rec_bias_ramp = 1, bias_year = c(2, 4, 10, 12)), dims = list(n_regions = 1))), "bias ramp")
  none <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem", do_rec_bias_ramp = 1, bias_year = rep(999, 4)), dims = list(n_regions = 1)))
  expect_equal(none$data$dsem_declared, "rec")

  # declared but never set up: the deviations would have no density
  expect_error(fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE), "no dsem was set up")

  # the declaration is the default process list, and the arrows then need every recruitment series
  n_yrs <- length(il$data$years)
  env <- data.frame(year = il$data$years, env = rnorm(n_yrs))
  d <- suppressMessages(Setup_Mod_DSEM(il, c("env -> rec, 0, b", "env <-> env, 0, s", "rec <-> rec, 0, sr"), env, dsem_mu_spec = "fix"))
  expect_equal(d$data$dsem_link_par, "ln_RecDevs")
  expect_error(Setup_Mod_DSEM(il, c("env -> env, 1, r", "env <-> env, 0, s"), env, dsem_mu_spec = "fix"), "No deviation series")
  obj <- fit_model(d$data, d$par, d$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj$fn(obj$par)))

  # sigmaR is the arrows' sd line: ln_sigmaR's own value is not read, and the sd line moves the initial age deviations' penalty
  expect_equal(d$data$dsem_link_sd_arrow, which(d$data$dsem_model$arrows$type == "sd" & d$data$dsem_model$arrows$to == "rec"))
  other <- d
  other$par$ln_sigmaR[] <- log(99)
  expect_equal(fit_model(other$data, other$par, other$map, random = NULL, do_optim = FALSE, silent = TRUE)$fn(obj$par), obj$fn(obj$par))
  wider <- d
  wider$par$ln_dsem_sd[match("sr", d$data$dsem_model$ln_sd_names)] <- log(3)
  wider$par$ln_InitDevs[] <- 0.5
  d$par$ln_InitDevs[] <- 0.5
  init_pen <- function(il) fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)$report()$Init_Rec_nLL
  expect_false(isTRUE(all.equal(init_pen(wider), init_pen(d)))) # the sd line reaches the initial age deviations' penalty

})

test_that("a recruitment link refuses terminal years left out of the deviations", {

  # without the declaration the terminal years can be dropped, and a link then meets grid rows with no deviation behind them
  il <- suppressMessages(sweep_input(rec = list(dont_est_recdev_last = 2), dims = list(n_regions = 1)))
  n_yrs <- length(il$data$years)
  expect_equal(dim(il$par$ln_RecDevs)[3], n_yrs - 2)
  env <- data.frame(year = il$data$years, env = rnorm(n_yrs))
  arrows <- c("env -> rec, 0, b", "env <-> env, 0, s", "rec <-> rec, 0, sr")
  expect_error(Setup_Mod_DSEM(il, arrows, env, dsem_processes = "rec", dsem_mu_spec = "fix"), "dont_est_recdev_last")

})

test_that("a declared process with a series left out of the arrows is refused", {

  sweep <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"))) # three regions, three recruitment series
  expect_equal(sweep$data$dsem_declared, "rec")
  n_yrs <- length(sweep$data$years)
  env <- data.frame(year = sweep$data$years, env = rnorm(n_yrs))
  one <- c("env -> rec_Pop_1_Region_1, 0, b", "env <-> env, 0, s", "rec_Pop_1_Region_1 <-> rec_Pop_1_Region_1, 0, sr")
  expect_error(Setup_Mod_DSEM(sweep, one, env, dsem_mu_spec = "fix"), "rec_Pop_1_Region_2")

  all_three <- c(one, "rec_Pop_1_Region_2 <-> rec_Pop_1_Region_2, 0, sr", "rec_Pop_1_Region_3 <-> rec_Pop_1_Region_3, 0, sr")
  d <- suppressMessages(Setup_Mod_DSEM(sweep, all_three, env, dsem_mu_spec = "fix"))
  expect_equal(sort(d$data$dsem_link_col), 2:4)

})

test_that("linking every cell without the declaration is refused while the sigma is estimated", {

  # the sweep model estimates sigmaR by default; every recruitment series linked leaves it read by nothing but the initial deviations
  sweep <- suppressMessages(sweep_input(rec = list(sigmaR_spec = "est_all")))
  expect_true(any(!is.na(sweep$map$ln_sigmaR)))
  n_yrs <- length(sweep$data$years)
  env <- data.frame(year = sweep$data$years, env = rnorm(n_yrs))
  all_three <- c("env <-> env, 0, s", paste0("rec_Pop_1_Region_", 1:3, " <-> rec_Pop_1_Region_", 1:3, ", 0, sr"))
  expect_error(Setup_Mod_DSEM(sweep, all_three, env, dsem_mu_spec = "fix"), "RecDevs_model = 'dsem'")

  # a partial link keeps the penalty reading sigmaR on the other regions, so it may stay estimated
  d <- suppressMessages(Setup_Mod_DSEM(sweep, c("env <-> env, 0, s", "rec_Pop_1_Region_1 <-> rec_Pop_1_Region_1, 0, sr"), env, dsem_mu_spec = "fix"))
  expect_equal(d$data$dsem_link_par, "ln_RecDevs")

})

test_that("NAA_re = 'dsem' keeps the state, holds sigmaNAA and defaults the process list", {

  il <- dusky_with(biol = list(NAA_re = "dsem", NAA_re_ages = 5:9))
  expect_equal(il$data$dsem_declared, "NAA")
  expect_equal(il$data$NAA_re, 1) # the iid code
  expect_true(all(is.na(il$map$ln_sigmaNAA)))
  expect_equal(sum(!is.na(il$map$ln_NAA)), 5 * (length(il$data$years) - 1))
  expect_error(dusky_with(biol = list(NAA_re = "dsem", NAA_re_ages = 5:9, NAA_sigma_spec = "est")), "cannot be estimated")

  naa <- function(a) sprintf("NAA_Pop_1_Region_1_Seas_1_Age_%d_Sex_1", a)
  four <- paste0(naa(2:5), " <-> ", naa(2:5), ", 0, sd_naa")
  expect_error(Setup_Mod_DSEM(il, four, dsem_data = NULL), naa(6))
  five <- paste0(naa(2:6), " <-> ", naa(2:6), ", 0, sd_naa")
  d <- suppressMessages(Setup_Mod_DSEM(il, hold_arrows(five, c(sd_naa = 0.3)), dsem_data = NULL))
  expect_equal(unique(d$data$dsem_link_par), "ln_NAA")

})

test_that("growth declares parameter by parameter, holds that parameter's sd, and owes only its series", {

  data("sgl_rg_ebs_pcod_data", envir = environment())
  base <- suppressWarnings(suppressMessages(build_ebs_pcod_input(sgl_rg_ebs_pcod_data))) # L1 and K vary as iid
  build <- build_ebs_pcod_input
  body(build) <- parse(text = gsub('L1 = "iid"', 'L1 = "dsem"', paste(deparse(body(build))), fixed = TRUE))[[1]]
  il <- suppressWarnings(suppressMessages(build(sgl_rg_ebs_pcod_data)))
  expect_equal(il$data$dsem_declared, "growth")
  expect_equal(il$data$growth_tv_dsem, c(1, 0, 0, 0, 0, 0)) # L1 declared, K (the third parameter) not
  expect_equal(il$data$growth_tv_model, base$data$growth_tv_model) # dsem reuses the iid code, so the codes do not move
  # L1's process error slot is off now; K's is whatever the base build made it
  pe_map <- array(as.integer(il$map$growth_pe_pars), dim = dim(il$par$growth_pe_pars))
  base_pe <- array(as.integer(base$map$growth_pe_pars), dim = dim(base$par$growth_pe_pars))
  expect_true(is.na(pe_map[1,1,1,1,1]))
  expect_equal(is.na(pe_map[1,1,3,1,1]), is.na(base_pe[1,1,3,1,1]))

  n_yrs <- length(il$data$years)
  env <- data.frame(year = il$data$years, env = rnorm(n_yrs))
  L1 <- "growth_Pop_1_Region_1_Par_1_Sex_1"
  K <- "growth_Pop_1_Region_1_Par_3_Sex_1"
  # K's series is not owed: only the declared parameter has to be in the arrows
  d <- suppressMessages(Setup_Mod_DSEM(il, c(paste0("env -> ", L1, ", 0, b"), "env <-> env, 0, s", paste0(L1, " <-> ", L1, ", 0, sg")), env, dsem_mu_spec = "fix"))
  expect_equal(d$data$dsem_link_par, "ln_growth_devs")
  expect_error(Setup_Mod_DSEM(il, c(paste0("env -> ", K, ", 0, b"), "env <-> env, 0, s", paste0(K, " <-> ", K, ", 0, sg")), env, dsem_mu_spec = "fix"), L1)

})

test_that("movement declares through the option's prefix and holds its process error", {

  # a dim the form leaves out is shared, and a shared deviation cannot be linked, so the form names every dim
  expect_error(suppressMessages(sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "dsem_y_a"))), "dsem_y_a_s")
  sweep <- suppressMessages(sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "dsem")))
  expect_equal(sweep$data$dsem_declared, "move")
  expect_equal(sweep$data$move_dsem, 1)
  expect_equal(sweep$data$cont_vary_movement, "iid_y_a_s") # "dsem" named year, age and sex itself
  spelled <- suppressMessages(sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "dsem_y_a_s")))
  expect_equal(spelled$data$cont_vary_movement, "iid_y_a_s")
  any_order <- suppressMessages(sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "dsem_s_a_y")))
  expect_equal(any_order$data$cont_vary_movement, "iid_y_a_s") # the dims are read in the order p, y, seas, a, s however they are written

  # the form is stored as written and the penalty parses its dims back out, so any combination of dims runs
  seasonal <- suppressMessages(sweep_input(dims = list(n_seas = 2, n_sexes = 1), move = list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "iid_y_seas_a", Movement_cont_pe_pars_spec = "est_shared")))
  expect_equal(seasonal$data$cont_vary_movement, "iid_y_seas_a")
  obj_seasonal <- fit_model(seasonal$data, seasonal$par, seasonal$map, random = NULL, do_optim = FALSE, silent = TRUE)
  expect_true(is.finite(obj_seasonal$fn(obj_seasonal$par)))
  expect_true(all(is.na(sweep$map$move_pe_pars)))
  expect_error(suppressMessages(sweep_input(move = list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "dsem_y_a_s", Movement_cont_pe_pars_spec = "est_all"))), "cannot be estimated")

  # every series with an estimated cell is owed: here one per origin, destination and age
  n_yrs <- length(sweep$data$years)
  offered <- get_dsem_link(sweep, "move", "", n_yrs)
  est <- vapply(seq_along(offered$offered), function(i) any(!is.na(as.integer(sweep$map$move_devs)[offered$offered_cell[[i]]])), logical(1))
  owed <- offered$offered[est]
  expect_gt(length(owed), 1)
  env <- data.frame(year = sweep$data$years, env = rnorm(n_yrs))
  one <- c("env <-> env, 0, s", paste0(owed[1], " <-> ", owed[1], ", 0, sm"))
  expect_error(Setup_Mod_DSEM(sweep, one, env, dsem_mu_spec = "fix"), owed[2])
  all_owed <- c("env <-> env, 0, s", paste0(owed, " <-> ", owed, ", 0, sm"))
  d <- suppressMessages(Setup_Mod_DSEM(sweep, all_owed, env, dsem_mu_spec = "fix"))
  expect_equal(unique(d$data$dsem_link_par), "move_devs")
  expect_length(d$data$dsem_link_par, length(owed))

})

test_that("a year varying preference term is refused once a movement series is linked", {

  # the CTMC deviations are the year to year part of preference, so a preference covariate that
  # varies over years fits that part twice, and the formula's version of it has no penalty
  n_regions <- 3; n_yrs <- 13; n_ages <- 7; n_sexes <- 2
  A <- matrix(1, n_regions, n_regions); diag(A) <- 0
  dat <- expand.grid(pop = 1, regions = 1:n_regions, years = 1:n_yrs, seas = 1, ages = 1:n_ages, sexes = 1:n_sexes)
  set.seed(4)
  dat$temp <- rnorm(n_yrs)[dat$years]

  ctmc <- function(pref) suppressMessages(sweep_input(move = list(
    use_fixed_movement = 0, Fixed_Movement = NA, move_type = 1, adjacency_mat = A,
    area_r = rep(1, n_regions), ctmc_move_dat = dat, diffusion_formula = ~1,
    preference_formula = pref, ctmc_diffusion_bounds = "upwind", cont_vary_movement = "dsem")))
  sd_lines <- function(il) {
    s <- dsem_series(il, "move")
    paste0(s, " <-> ", s, ", 0, sd_move")
  }

  # region factors and an age spline stay the same from year to year, so they sit beside the dsem
  flat <- ctmc(~0 + factor(regions))
  expect_equal(unique(suppressMessages(Setup_Mod_DSEM(flat, sd_lines(flat), NULL))$data$dsem_link_par), "move_devs")
  aged <- ctmc(~0 + factor(regions):splines2::bSpline(ages, df = 4, intercept = TRUE))
  expect_no_error(suppressMessages(Setup_Mod_DSEM(aged, sd_lines(aged), NULL)))

  # a covariate read by year, and a year factor, are the two ways of writing it twice
  covar <- ctmc(~0 + factor(regions):temp)
  expect_error(Setup_Mod_DSEM(covar, sd_lines(covar), NULL), "vary over years")
  expect_error(Setup_Mod_DSEM(covar, sd_lines(covar), NULL), "temp")
  yrs <- ctmc(~0 + factor(regions) + factor(years))
  expect_error(Setup_Mod_DSEM(yrs, sd_lines(yrs), NULL), "factor\\(years\\)")

  # the check reads the preference design alone, so it is silent without a movement series
  expect_length(SPoRC:::get_yr_varying_pref_terms(flat), 0)
  expect_gt(length(SPoRC:::get_yr_varying_pref_terms(covar)), 0)

})

test_that("the semi-parametric surface declares like the rest", {

  data("sgl_rg_ebs_pcod_data", envir = environment())
  build <- build_ebs_pcod_input
  body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Biologicals = quote(function(...) Setup_Mod_Biologicals(..., growth_semipar = "dsem")))))
  il <- suppressWarnings(suppressMessages(build(sgl_rg_ebs_pcod_data)))
  expect_true("growth_semipar" %in% il$data$dsem_declared)
  expect_equal(il$data$growth_semipar, 1) # the iid code
  expect_equal(il$data$growth_semipar_dsem, 1)
  expect_true(all(is.na(array(as.integer(il$map$growth_pe_pars), dim = dim(il$par$growth_pe_pars))[,,,,2]))) # the surface's slots are off

  n_yrs <- length(il$data$years)
  env <- data.frame(year = il$data$years, env = rnorm(n_yrs))
  a3 <- "growth_semipar_Pop_1_Region_1_Age_3_Sex_1"
  expect_error(Setup_Mod_DSEM(il, c("env <-> env, 0, s", paste0(a3, " <-> ", a3, ", 0, sg")), env, dsem_mu_spec = "fix"), "growth_semipar_Pop_1_Region_1_Age_1_Sex_1")

})

test_that("a linked series refuses cells that are mapped off", {

  # the series covers every year, so a mapped-off cell would enter it as a fixed value rather than drop out
  il <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem"), dims = list(n_regions = 1)))
  n_yrs <- length(il$data$years)
  env <- data.frame(year = il$data$years, env = rnorm(n_yrs))
  arrows <- c("env -> rec, 0, b", "env <-> env, 0, s", "rec <-> rec, 0, sr")
  expect_no_error(suppressMessages(Setup_Mod_DSEM(il, arrows, env, dsem_mu_spec = "fix")))

  lev <- if(is.null(il$map$ln_RecDevs)) seq_along(il$par$ln_RecDevs) else as.integer(il$map$ln_RecDevs)
  lev[1:3] <- NA
  il$map$ln_RecDevs <- factor(lev)
  expect_error(suppressMessages(Setup_Mod_DSEM(il, arrows, env, dsem_mu_spec = "fix")), "mapped off")
})
