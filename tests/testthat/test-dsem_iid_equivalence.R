# A dsem written as bare sd lines (no covariate, no path) is the iid penalty under another name, so the joint
# negative log likelihood has to come out the same (1e-8). Every process a dsem can link, both routes.

data("sgl_rg_dusky_data")
data("sgl_rg_ebs_pcod_data")

# dusky with the arguments of the rec or biol stage overridden, so one build can declare "dsem" or add a state.
# an override replaces an argument the builder passes itself, such as the bias ramp
dusky_built_with <- function(rec = list(), biol = list()) {
  build <- build_goa_dusky_input
  with_overrides <- function(stage, overrides) function(...) { given <- list(...); given[names(overrides)] <- NULL; do.call(stage, c(given, overrides)) }
  subs <- list()
  if(length(rec) > 0) subs$Setup_Mod_Rec <- with_overrides(Setup_Mod_Rec, rec)
  if(length(biol) > 0) subs$Setup_Mod_Biologicals <- with_overrides(Setup_Mod_Biologicals, biol)
  body(build) <- do.call(substitute, list(body(build), subs))
  suppressMessages(build(sgl_rg_dusky_data))
}

# objective, gradient and report at the starting values, through the same fit_model route a user takes
value_of <- function(il) {
  obj <- fit_model(il$data, il$par, il$map, random = NULL, do_optim = FALSE, silent = TRUE)
  list(fn = as.numeric(obj$fn(obj$par)), gr = as.numeric(obj$gr(obj$par)), par = obj$par, rep = obj$rep)
}

# every reported likelihood summed to one number, so the parts can be compared and not only the total
component_sums <- function(rep) {
  nm <- grep("_nLL$", names(rep), value = TRUE)
  stats::setNames(vapply(nm, function(n) sum(as.numeric(rep[[n]])), 0), nm)
}

# the total agrees, the linked penalty is off, the dsem holds what it held, and nothing else moved
expect_same_jnll <- function(via, plain, moved, tolerance = 1e-8) {
  expect_equal(via$fn, plain$fn, tolerance = tolerance)
  cp <- component_sums(plain$rep)
  cv <- component_sums(via$rep)
  expect_equal(unname(cv[moved]), 0)
  expect_equal(unname(cv["dsem_nLL"]), unname(cp[moved]), tolerance = tolerance)
  others <- setdiff(names(cp), c(moved, "dsem_nLL", "dsem_obs_nLL"))
  expect_equal(cv[others], cp[others], tolerance = tolerance)
}

# sd lines for a set of series, each fixed at one value: "name <-> name, 0, NA, value"
sd_lines <- function(series, sd) sprintf("%s <-> %s, 0, NA, %.17g", series, series, sd)

# pcod with the biological stage's arguments overridden, which is where growth is set up
pcod_built_with <- function(biol) {
  build <- build_ebs_pcod_input
  with_overrides <- function(stage, overrides) function(...) { given <- list(...); given[names(overrides)] <- NULL; do.call(stage, c(given, overrides)) }
  body(build) <- do.call(substitute, list(body(build), list(Setup_Mod_Biologicals = with_overrides(Setup_Mod_Biologicals, biol))))
  suppressWarnings(suppressMessages(build(sgl_rg_ebs_pcod_data)))
}

# three regions with catchability deviations on every fleet, their sd fixed so an sd line can stand in for it
q_built_with <- function(fish, srv) {
  suppressMessages(sweep_input(
    fishsel = list(fish_q_model = rep(fish, sweep_dims$n_fish_fleets), sigma_fish_q_spec = "fix"),
    srvsel = list(srv_q_model = rep(srv, sweep_dims$n_srv_fleets), sigma_srv_q_spec = "fix")
  ))
}

# the same deviations in every build, and sds set apart from each other by region and fleet
q_devs_at <- function(il) {
  set.seed(21)
  il$par$ln_fish_q_devs[] <- rnorm(length(il$par$ln_fish_q_devs), 0, 0.15)
  il$par$ln_srv_q_devs[] <- rnorm(length(il$par$ln_srv_q_devs), 0, 0.15)
  il$par$ln_sigma_fish_q[] <- log(seq(0.06, 0.2, length.out = length(il$par$ln_sigma_fish_q)))
  il$par$ln_sigma_srv_q[] <- log(seq(0.08, 0.15, length.out = length(il$par$ln_sigma_srv_q)))
  il
}

test_that("dusky: an sd-only link on recruitment is the iid penalty, declared or not", {

  plain <- dusky_built_with()
  set.seed(1)
  plain$par$ln_RecDevs[] <- rnorm(length(plain$par$ln_RecDevs), 0, 0.6) # off zero, so the penalty is not trivially zero
  sigmaR <- exp(plain$par$ln_sigmaR[2,1,1])
  base <- value_of(plain)

  # not declared: the rec module keeps iid, and the link switches its penalty off cell by cell
  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines("rec", sigmaR), dsem_data = NULL, dsem_processes = "rec"))
  expect_length(linked$par$dsem_beta, 0) # a fixed sd line leaves no dsem parameter behind
  expect_length(linked$par$ln_dsem_sd, 0)
  expect_true(all(is.na(linked$map$dsem_x))) # every grid cell is a linked deviation
  via <- value_of(linked)
  expect_same_jnll(via, base, moved = "Rec_nLL")
  expect_equal(via$gr, base$gr, tolerance = 1e-6) # same parameters in the same order, so the gradients match too

  # declared: ln_sigmaR is read by nothing, the arrows' sd line stands in for it
  declared <- dusky_built_with(rec = list(RecDevs_model = "dsem"))
  declared$par$ln_RecDevs <- plain$par$ln_RecDevs
  expect_error(fit_model(declared$data, declared$par, declared$map, do_optim = FALSE, silent = TRUE), "no dsem was set up")
  linked_decl <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines("rec", sigmaR), dsem_data = NULL))
  via_decl <- value_of(linked_decl)
  expect_same_jnll(via_decl, base, moved = "Rec_nLL")
  expect_equal(via_decl$gr, base$gr, tolerance = 1e-6)

  # the sd line really is what is read: moving ln_sigmaR changes nothing, moving the line does
  other <- linked_decl
  other$par$ln_sigmaR[] <- log(99)
  expect_equal(value_of(other)$fn, via_decl$fn)
  wider <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines("rec", 2 * sigmaR), dsem_data = NULL))
  expect_false(isTRUE(all.equal(value_of(wider)$fn, via_decl$fn)))

})

test_that("dusky: estimating the sd through the arrows or through ln_sigmaR is the same objective", {

  # plain iid with sigmaR free. the array is [2, pop, region]: row 2 is read by the recruitment penalty and
  # row 1 by the initial age deviations. the arrows' sd line stands in for both, so tie both rows to one parameter
  plain <- dusky_built_with()
  plain$map$ln_sigmaR <- factor(c(1, 1))
  set.seed(3)
  plain$par$ln_RecDevs[] <- rnorm(length(plain$par$ln_RecDevs), 0, 0.6)
  sigmaR <- exp(plain$par$ln_sigmaR[2,1,1])

  # declared, with a named sd line started at the same value, so ln_dsem_sd is the one free sd
  declared <- dusky_built_with(rec = list(RecDevs_model = "dsem"))
  declared$par$ln_RecDevs <- plain$par$ln_RecDevs
  linked <- suppressMessages(Setup_Mod_DSEM(declared, sprintf("rec <-> rec, 0, sd_rec, %.17g", sigmaR), dsem_data = NULL))
  expect_equal(linked$data$dsem_model$ln_sd_names, "sd_rec")
  expect_equal(exp(linked$par$ln_dsem_sd), sigmaR)

  base <- value_of(plain)
  via <- value_of(linked)
  expect_equal(via$fn, base$fn, tolerance = 1e-8)

  # the free sd sits at a different position in each parameter vector, so take it out and compare the rest by position
  expect_equal(via$gr[names(via$par) == "ln_dsem_sd"], base$gr[names(base$par) == "ln_sigmaR"], tolerance = 1e-6)
  expect_equal(via$gr[names(via$par) != "ln_dsem_sd"], base$gr[names(base$par) != "ln_sigmaR"], tolerance = 1e-6)

  # marginal fits, the deviations integrated out, land on the same objective and the same sd
  skip_on_cran()
  fit_plain <- fit_model(plain$data, plain$par, plain$map, random = "ln_RecDevs", silent = TRUE)
  fit_via <- fit_model(linked$data, linked$par, linked$map, random = "ln_RecDevs", silent = TRUE)
  expect_equal(fit_via$optim$objective, fit_plain$optim$objective, tolerance = 1e-6)
  expect_equal(exp(fit_via$env$parList()$ln_dsem_sd), exp(fit_plain$env$parList()$ln_sigmaR[2,1,1]), tolerance = 1e-4)

})

test_that("dusky: iid arrows on a numbers at age state are the native iid penalty, declared or not", {

  # a state on ages 5 to 9 (indices 2 to 6 of the model ages), sigma fixed at the module's start of 0.3
  biol <- list(NAA_re = "iid", NAA_re_ages = 5:9, NAA_sigma_spec = "fix")
  plain <- dusky_built_with(biol = biol)
  n_yrs <- length(plain$data$years)
  age_idx <- match(5:9, plain$data$ages)
  sigmaNAA <- exp(plain$par$ln_sigmaNAA[1])
  series <- dsem_series(plain, "NAA") # the state's ages, by array index
  expect_equal(series, sprintf("NAA_Pop_1_Region_1_Seas_1_Age_%d_Sex_1", age_idx))

  # put the state off its prediction, otherwise every penalty is zero
  pred <- value_of(plain)$rep$NAA[,,1:n_yrs,,,]
  set.seed(4)
  plain$par$ln_NAA[] <- log(pred) + rnorm(length(plain$par$ln_NAA), 0, 0.2)
  base <- value_of(plain)

  # the state starts in year two, so the grid's year one rows are fixed at the mean and add a known
  # constant, one dnorm(0, 0, sigma) per series. take it off and the density is the native one
  year_one <- length(series) * (log(sigmaNAA) + 0.5 * log(2 * pi))
  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, sigmaNAA), dsem_data = NULL, dsem_processes = "NAA"))
  via <- value_of(linked)
  expect_equal(via$fn - year_one, base$fn, tolerance = 1e-8)
  expect_equal(sum(via$rep$NAA_state_nLL), 0, tolerance = 1e-12)
  expect_equal(via$rep$dsem_nLL - year_one, sum(base$rep$NAA_state_nLL), tolerance = 1e-8)

  declared <- dusky_built_with(biol = utils::modifyList(biol, list(NAA_re = "dsem")))
  declared$par$ln_NAA <- plain$par$ln_NAA
  linked_decl <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series, sigmaNAA), dsem_data = NULL))
  expect_equal(value_of(linked_decl)$fn - year_one, base$fn, tolerance = 1e-8)

})

test_that("three regions: sd lines on every recruitment series are the iid penalty, declared or not", {

  plain <- suppressMessages(sweep_input())
  set.seed(5)
  plain$par$ln_RecDevs[] <- rnorm(length(plain$par$ln_RecDevs), 0, 0.5)
  sigmaR <- exp(plain$par$ln_sigmaR[2,1,1])
  expect_true(all(plain$par$ln_sigmaR[2,1,] == log(sigmaR))) # one value serves every region here
  series <- dsem_series(plain, "rec")
  expect_equal(series, paste0("rec_Pop_1_Region_", 1:3))
  base <- value_of(plain)

  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, sigmaR), dsem_data = NULL, dsem_processes = "rec"))
  expect_equal(linked$data$dsem_var_names, series)
  expect_same_jnll(value_of(linked), base, moved = "Rec_nLL")

  # declared, and a declared process refuses a series left out
  declared <- suppressMessages(sweep_input(rec = list(RecDevs_model = "dsem")))
  declared$par$ln_RecDevs <- plain$par$ln_RecDevs
  expect_error(suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series[1:2], sigmaR), dsem_data = NULL)), "rec_Pop_1_Region_3")
  linked_decl <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series, sigmaR), dsem_data = NULL))
  expect_same_jnll(value_of(linked_decl), base, moved = "Rec_nLL")

  # one region linked, two left alone: the penalty stays on the other two and the total still agrees
  one <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series[2], sigmaR), dsem_data = NULL, dsem_processes = "rec"))
  via_one <- value_of(one)
  expect_equal(via_one$fn, base$fn, tolerance = 1e-8)
  expect_equal(sum(via_one$rep$Rec_nLL[1,2,]), 0)
  expect_equal(sum(via_one$rep$Rec_nLL[1,c(1, 3),]), sum(base$rep$Rec_nLL[1,c(1, 3),]), tolerance = 1e-8)

})

test_that("three regions: sd lines on every movement series are the movement penalty, declared or not", {

  # movement estimated with iid deviations by year, age and sex, their sd fixed at its start
  move <- list(use_fixed_movement = 0, Fixed_Movement = NA, cont_vary_movement = "iid_y_a_s", Movement_cont_pe_pars_spec = "fix")
  plain <- suppressMessages(sweep_input(move = move))
  set.seed(6)
  plain$par$move_devs[] <- rnorm(length(plain$par$move_devs), 0, 0.3)
  sd_move <- exp(plain$par$move_pe_pars[1]) # one value for every from region, age and sex
  expect_true(all(plain$par$move_pe_pars == plain$par$move_pe_pars[1]))
  series <- dsem_series(plain, "move") # every series with a cell the map estimates
  dims <- dim(plain$par$move_devs) # [pop, from, to, year, seas, age, sex], the last destination being the reference
  expect_equal(length(series), dims[2] * dims[3] * (dims[6] - 1) * dims[7]) # every free cell but age one, since recruits do not move
  base <- value_of(plain)
  expect_gt(sum(base$rep$Movement_nLL), 0)

  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, sd_move), dsem_data = NULL, dsem_processes = "move"))
  expect_equal(linked$data$dsem_link_par, rep("move_devs", length(series)))
  expect_same_jnll(value_of(linked), base, moved = "Movement_nLL")

  # declared: cont_vary_movement = "dsem" names every dim itself (year, age and sex here)
  declared <- suppressMessages(sweep_input(move = utils::modifyList(move, list(cont_vary_movement = "dsem", Movement_cont_pe_pars_spec = NULL))))
  declared$par$move_devs <- plain$par$move_devs
  expect_equal(declared$data$dsem_declared, "move")
  expect_equal(dsem_series(declared, "move"), series)
  linked_decl <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series, sd_move), dsem_data = NULL))
  expect_same_jnll(value_of(linked_decl), base, moved = "Movement_nLL")

})

test_that("dusky: the random walk and AR1 penalties are the arrows under the right variance form", {

  # both links go through dsem_processes = "rec" with the rec module left on its own model, so the native penalty
  # comes off cell by cell and everything else, the initial age deviations included, still reads ln_sigmaR
  walk_devs <- function(il, seed) {
    set.seed(seed)
    il$par$ln_RecDevs[] <- cumsum(rnorm(length(il$par$ln_RecDevs), 0, 0.3))
    il
  }

  # dusky holds its ramp at zero through the switch, which the walk and the ar1 refuse, so the switch goes off for them:
  # neither penalty takes a correction, and a link under a rec module that is not iid takes none either
  # random walk: year one at sigma (RecDevs_rw_init_sigma = NA), every later year about the one before at sigma.
  # the arrows' walk under the conditional form is that density, since a walk's innovation sd is the sd line
  plain_rw <- walk_devs(suppressWarnings(dusky_built_with(rec = list(RecDevs_model = "rw", RecDevs_rw_init_sigma = NA, do_rec_bias_ramp = 0))), 7)
  sigmaR <- exp(plain_rw$par$ln_sigmaR[2,1,1])
  base_rw <- value_of(plain_rw)
  expect_gt(sum(base_rw$rep$Rec_nLL), 0)
  walk <- suppressMessages(Setup_Mod_DSEM(plain_rw, c("rec -> rec, 1, NA, 1", sd_lines("rec", sigmaR)), dsem_data = NULL, dsem_processes = "rec"))
  expect_same_jnll(value_of(walk), base_rw, moved = "Rec_nLL")

  # the native default draws year one wide (sd 5), which the arrows have no line for, so that setting differs by year one alone
  wide <- walk_devs(suppressWarnings(dusky_built_with(rec = list(RecDevs_model = "rw", do_rec_bias_ramp = 0))), 7)
  d1 <- wide$par$ln_RecDevs[1,1,1]
  expect_equal(value_of(wide)$fn - base_rw$fn, dnorm(d1, 0, sigmaR, log = TRUE) - dnorm(d1, 0, 5, log = TRUE), tolerance = 1e-8)

  # AR1: year one at the stationary sd sigma / sqrt(1 - rho^2), every later year about rho times the one before at sigma.
  # under the diagonal form an sd line at the stationary sd gives exactly that: year one at the line, later years at
  # line * sqrt(1 - rho^2) = sigma. under the conditional form year one sits at sigma instead
  rho <- 0.6
  plain_ar1 <- walk_devs(suppressWarnings(dusky_built_with(rec = list(RecDevs_model = "ar1", do_rec_bias_ramp = 0))), 8)
  plain_ar1$par$RecDevs_rho[] <- atanh(rho) # the penalty reads 2 / (1 + exp(-2 x)) - 1, which is tanh
  base_ar1 <- value_of(plain_ar1)
  stationary_sd <- sigmaR / sqrt(1 - rho^2)
  self_path <- sprintf("rec -> rec, 1, NA, %.17g", rho)
  ar1_diag <- suppressMessages(Setup_Mod_DSEM(plain_ar1, c(self_path, sd_lines("rec", stationary_sd)), dsem_data = NULL, dsem_processes = "rec", dsem_variance = "diagonal"))
  expect_same_jnll(value_of(ar1_diag), base_ar1, moved = "Rec_nLL")

  ar1_cond <- suppressMessages(Setup_Mod_DSEM(plain_ar1, c(self_path, sd_lines("rec", sigmaR)), dsem_data = NULL, dsem_processes = "rec"))
  d1 <- plain_ar1$par$ln_RecDevs[1,1,1]
  expect_equal(value_of(ar1_cond)$fn - base_ar1$fn, dnorm(d1, 0, stationary_sd, log = TRUE) - dnorm(d1, 0, sigmaR, log = TRUE), tolerance = 1e-8)

  # and the marginal form is the diagonal form here, there being no covariance line
  ar1_marg <- suppressMessages(Setup_Mod_DSEM(plain_ar1, c(self_path, sd_lines("rec", stationary_sd)), dsem_data = NULL, dsem_processes = "rec", dsem_variance = "marginal"))
  expect_equal(value_of(ar1_marg)$fn, base_ar1$fn, tolerance = 1e-8)

})

test_that("pcod: sd lines on the time-varying growth parameters are the iid growth penalty, declared or not", {

  # pcod varies L1 and K as iid over part of the series. every year is estimated here, so the dsem grid and
  # the native penalty read the same cells, and each parameter keeps its own fixed sd
  every_yr <- list(L1 = sgl_rg_ebs_pcod_data$years, K = sgl_rg_ebs_pcod_data$years)
  plain <- pcod_built_with(list(growth_tv_years = every_yr))
  n_yrs <- length(plain$data$years)
  set.seed(11)
  plain$par$ln_growth_devs[1,1,,1,1] <- rnorm(n_yrs, 0, 0.2)
  plain$par$ln_growth_devs[1,1,,3,1] <- rnorm(n_yrs, 0, 0.05)

  series <- dsem_series(plain, "growth")
  expect_equal(series, paste0("growth_Pop_1_Region_1_Par_", c(1, 3), "_Sex_1")) # L1 and K, by array index
  sds <- exp(plain$par$growth_pe_pars[1,1,c(1, 3),1,1]) # the first data source of the shared process error array
  base <- value_of(plain)
  expect_gt(abs(sum(base$rep$growth_tv_nLL)), 1) # the penalty has something in it, so the comparison is not empty

  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, sds), dsem_data = NULL, dsem_processes = "growth"))
  via <- value_of(linked)
  expect_same_jnll(via, base, moved = "growth_tv_nLL")
  expect_equal(via$gr, base$gr, tolerance = 1e-6)

  # the two parameters have different sds, so reading each series under the other's line has to move the total
  swapped <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, rev(sds)), dsem_data = NULL, dsem_processes = "growth"))
  expect_false(isTRUE(all.equal(value_of(swapped)$fn, base$fn)))

  # declared: the growth codes stay on iid, so blanking the map mirror is the only thing taking the penalty off
  declared <- pcod_built_with(list(growth_tv_years = every_yr, growth_tv_model = c(L1 = "dsem", K = "dsem")))
  declared$par$ln_growth_devs <- plain$par$ln_growth_devs
  expect_equal(declared$data$dsem_declared, "growth")
  expect_equal(declared$data$growth_tv_model, plain$data$growth_tv_model)
  expect_equal(declared$data$growth_tv_dsem, c(1, 0, 1, 0, 0, 0))
  linked_decl <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series, sds), dsem_data = NULL))
  expect_same_jnll(value_of(linked_decl), base, moved = "growth_tv_nLL")

})

test_that("pcod: sd lines on the semi-parametric surface are its iid penalty, one series per age", {

  # the surface is one deviation per year and age, its sd one value per age. the sds are set apart from each
  # other here so the arrows have to reach the right age
  plain <- pcod_built_with(list(growth_semipar = "iid"))
  n_ages <- length(plain$data$ages)
  plain$par$growth_pe_pars[1,1,,1,2] <- log(seq(0.03, 0.08, length.out = n_ages)) # second data source, fixed
  set.seed(12)
  plain$par$ln_growth_semipar_devs[] <- rnorm(length(plain$par$ln_growth_semipar_devs), 0, 0.03)

  series <- dsem_series(plain, "growth_semipar")
  expect_equal(series, sprintf("growth_semipar_Pop_1_Region_1_Age_%d_Sex_1", seq_len(n_ages)))
  sds <- exp(plain$par$growth_pe_pars[1,1,,1,2])
  base <- value_of(plain)
  expect_gt(abs(sum(base$rep$growth_semipar_nLL)), 1)

  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, sds), dsem_data = NULL, dsem_processes = "growth_semipar"))
  via <- value_of(linked)
  expect_same_jnll(via, base, moved = "growth_semipar_nLL")
  expect_equal(via$gr, base$gr, tolerance = 1e-6)

  # the sds rise with age, so reading them backwards has to move the total
  reversed <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, rev(sds)), dsem_data = NULL, dsem_processes = "growth_semipar"))
  expect_false(isTRUE(all.equal(value_of(reversed)$fn, base$fn)))

})

test_that("three regions: sd lines on the fishery catchability series are its deviation penalty, declared or not", {

  plain <- q_devs_at(q_built_with("iid", "iid"))
  series <- dsem_series(plain, "fish_q")
  expect_equal(series, paste0("fish_q_Region_", 1:3, "_Fleet_", rep(seq_len(sweep_dims$n_fish_fleets), each = 3)))
  sds <- exp(as.vector(plain$par$ln_sigma_fish_q)) # the sigma array is [region, fleet], the series order region first
  base <- value_of(plain)
  expect_gt(abs(sum(base$rep$fish_q_nLL)), 1)
  expect_equal(plain$data$Use_fish_q_prior, 0) # so the component is the deviation penalty alone

  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, sds), dsem_data = NULL, dsem_processes = "fish_q"))
  via <- value_of(linked)
  expect_same_jnll(via, base, moved = "fish_q_nLL") # the survey's own penalty sits in the untouched set
  expect_equal(via$gr, base$gr, tolerance = 1e-6)

  reversed <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, rev(sds)), dsem_data = NULL, dsem_processes = "fish_q"))
  expect_false(isTRUE(all.equal(value_of(reversed)$fn, base$fn)))

  # declared: fish_q_model = "dsem" fixes the sd itself, and every fleet owes its series
  declared <- q_devs_at(q_built_with("dsem", "iid"))
  expect_equal(declared$data$dsem_declared, "fish_q")
  expect_error(suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series[-1], sds[-1]), dsem_data = NULL)), series[1], fixed = TRUE)
  linked_decl <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series, sds), dsem_data = NULL))
  expect_same_jnll(value_of(linked_decl), base, moved = "fish_q_nLL")

})

test_that("three regions: sd lines on the survey catchability series are its deviation penalty, declared or not", {

  plain <- q_devs_at(q_built_with("iid", "iid"))
  series <- dsem_series(plain, "srv_q")
  expect_equal(series, paste0("srv_q_Region_", 1:3, "_Fleet_1"))
  sds <- exp(as.vector(plain$par$ln_sigma_srv_q))
  base <- value_of(plain)
  expect_gt(abs(sum(base$rep$srv_q_nLL)), 1)
  expect_equal(plain$data$Use_srv_q_prior, 0)

  linked <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, sds), dsem_data = NULL, dsem_processes = "srv_q"))
  via <- value_of(linked)
  expect_same_jnll(via, base, moved = "srv_q_nLL")
  expect_equal(via$gr, base$gr, tolerance = 1e-6)

  reversed <- suppressMessages(Setup_Mod_DSEM(plain, sd_lines(series, rev(sds)), dsem_data = NULL, dsem_processes = "srv_q"))
  expect_false(isTRUE(all.equal(value_of(reversed)$fn, base$fn)))

  declared <- q_devs_at(q_built_with("iid", "dsem"))
  expect_equal(declared$data$dsem_declared, "srv_q")
  linked_decl <- suppressMessages(Setup_Mod_DSEM(declared, sd_lines(series, sds), dsem_data = NULL))
  expect_same_jnll(value_of(linked_decl), base, moved = "srv_q_nLL")

})
