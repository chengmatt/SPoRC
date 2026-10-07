# Selectivity, F and discard mortality deviations in the operating model: drawn from the process the fit penalizes, F only
# when integrated out, and the estimation model at a replicate's draws reproduces its selectivity and F exactly.

# the state-space test model with time-varying fishery selectivity and the F deviation sd at sd_F
fleet_devs_input <- function(cont_tv = "iid_Fleet_1", sel_sd = c(0.3, 0.2), sd_F = 0.3) {
  set.seed(11)
  il <- naaom_build_em(naaom_om_data(naaom_make_om(NAA_re = "none")), NAA_re = "none")
  il <- Setup_Mod_Fishsel_and_Q(input_list = il, fish_sel_model = "logist1_Fleet_1", fish_fixed_sel_pars_spec = "est_all",
                                fish_q_spec = "est_all", use_fixed_ret_sel = 1, cont_tv_fish_sel = cont_tv,
                                fishsel_pe_pars_spec = "fix", fish_sel_devs_spec = "est_all")
  il$par$fishsel_pe_pars[1,1:2,1,1] <- log(sel_sd)
  il$par$ln_sigmaF[] <- log(sd_F)
  il$map$ln_fish_q <- factor(rep(NA, length(il$par$ln_fish_q))) # no fishery index
  il
}

# the operating model list the self test builds from a model at its starting values, taken before any replicate runs,
# with the years past n_cond_yrs left to the operating model to draw, as the closed loop draws its projection years
fleet_devs_sim_list <- function(il, random, n_sims, n_cond_yrs = 0, F_devs_draw = "random") {
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, random = random, do_optim = FALSE, silent = TRUE))
  captured <- new.env()
  testthat::with_mocked_bindings(
    Simulate_Pop_Static = function(sim_list, ...) { captured$sim_list <- sim_list; stop("captured") },
    try(suppressWarnings(suppressMessages(
      simulation_self_test(data = obj$data, parameters = obj$parameters, mapping = obj$mapping, random = random, rep = obj$rep,
                           sd_rep = list(par.fixed = obj$par, par.random = obj$env$par[obj$env$random]),
                           n_sims = n_sims, newton_loops = 0, what = "SSB"))), silent = TRUE),
    .package = "SPoRC")
  sim_list <- captured$sim_list
  sim_list$n_cond_yrs <- n_cond_yrs
  if(F_devs_draw == "all") sim_list <- Setup_Sim_Fleet_Devs(sim_list, obj$data, obj$env$parList(), random = random, F_devs_draw = "all")
  list(obj = obj, sim_list = sim_list)
}

test_that("selectivity and integrated F deviations are drawn at the fitted sd, and the estimation model at a draw agrees", {

  il <- fleet_devs_input()
  built <- fleet_devs_sim_list(il, random = c("ln_fishsel_devs", "ln_F_devs"), n_sims = 400)
  expect_equal(unname(built$sim_list$fleet_devs_on[c("F", "dmr", "fish", "ret", "srv")]), c(TRUE, FALSE, TRUE, FALSE, FALSE))

  set.seed(3)
  env <- suppressMessages(Setup_sim_env(built$sim_list))
  expect_equal(sd(env$ln_fishsel_devs[1,,1,1,1,]), 0.3, tolerance = 0.05) # 14000 draws per parameter
  expect_equal(sd(env$ln_fishsel_devs[1,,2,1,1,]), 0.2, tolerance = 0.05)
  expect_equal(sd(env$ln_F_devs[1,,1,1,]), 0.3, tolerance = 0.05)
  expect_lt(abs(mean(env$ln_F_devs[1,,1,1,])), 0.02)

  # F is exp(ln_F_mean + ln_F_devs), and the estimation model at one replicate's draws gives that replicate's selectivity and F
  pl <- built$obj$env$parList(par = built$obj$env$par)
  expect_equal(env$Fmort[1,,1,1,7], exp(pl$ln_F_mean[1,1,1] + env$ln_F_devs[1,,1,1,7]), tolerance = 1e-12)
  at_draw <- pl
  at_draw$ln_fishsel_devs[] <- env$ln_fishsel_devs[,,,,,7]
  at_draw$ln_F_devs[] <- env$ln_F_devs[,,,,7]
  rep_7 <- fit_model(il$data, at_draw, il$map, random = built$obj$random, do_optim = FALSE, silent = TRUE)$rep
  n_yrs <- length(il$data$years)
  expect_equal(as.vector(rep_7$fish_sel[,,1:n_yrs,,,,,drop = FALSE]), as.vector(env$fish_sel[,,1:n_yrs,,,,,7]), tolerance = 1e-12)
  expect_equal(as.vector(rep_7$Fmort[,1:n_yrs,,,drop = FALSE]), as.vector(env$Fmort[,1:n_yrs,,,7]), tolerance = 1e-12)

})

test_that("a weighted F penalty is drawn at its sd over the root of the weight, and a zero weight draws nothing", {

  # the fit penalizes at Wt_F times the density at sd_F, which is the density at sd_F over the root of Wt_F
  built <- fleet_devs_sim_list(fleet_devs_input(sd_F = 0.3), random = c("ln_fishsel_devs", "ln_F_devs"), n_sims = 400)
  built$sim_list$fleet_dev_data$Wt_F <- 4
  set.seed(6)
  env <- suppressMessages(Setup_sim_env(built$sim_list))
  expect_equal(sd(env$ln_F_devs[1,,1,1,]), 0.15, tolerance = 0.05)

  # a penalty at weight zero is no density, so F keeps the fit's deviations
  il <- fleet_devs_input()
  il$data$Wt_F <- 0
  built <- fleet_devs_sim_list(il, random = c("ln_fishsel_devs", "ln_F_devs"), n_sims = 2)
  expect_false(built$sim_list$fleet_devs_on[["F"]])

})

test_that("F deviations kept as fixed effects keep the fit's F", {
  built <- fleet_devs_sim_list(fleet_devs_input(), random = "ln_fishsel_devs", n_sims = 3)
  expect_false(built$sim_list$fleet_devs_on[["F"]])
  env <- suppressMessages(Setup_sim_env(built$sim_list))
  for(sim in 1:3) expect_equal(env$Fmort[,,,,sim], built$obj$rep$Fmort[,seq_along(built$obj$data$years),,], tolerance = 1e-12)
  expect_null(env$ln_F_devs)
})

test_that("F deviations kept as fixed effects are drawn from their penalty when asked for", {
  built <- fleet_devs_sim_list(fleet_devs_input(sd_F = 0.25), random = "ln_fishsel_devs", n_sims = 400, F_devs_draw = "all")
  expect_true(built$sim_list$fleet_devs_on[["F"]])
  set.seed(5)
  env <- suppressMessages(Setup_sim_env(built$sim_list))
  expect_equal(sd(env$ln_F_devs[1,,1,1,]), 0.25, tolerance = 0.05)
})

test_that("a selectivity walk keeps the fit's conditioned years and steps on from the last of them", {

  il <- fleet_devs_input(cont_tv = "rw_Fleet_1")
  il$par$ln_fishsel_devs[1,,1,1,1] <- seq(-0.5, 0.5, length.out = dim(il$par$ln_fishsel_devs)[2]) # the fit's walk, recognizable
  built <- fleet_devs_sim_list(il, random = "ln_fishsel_devs", n_sims = 2000, n_cond_yrs = 10)
  set.seed(4)
  env <- suppressMessages(Setup_sim_env(built$sim_list))
  fit_walk <- il$par$ln_fishsel_devs[1,,1,1,1]
  for(sim in 1:5) expect_equal(env$ln_fishsel_devs[1,1:10,1,1,1,sim], fit_walk[1:10])
  step_11 <- env$ln_fishsel_devs[1,11,1,1,1,] - fit_walk[10]
  expect_lt(abs(mean(step_11)), 0.03)
  expect_equal(sd(step_11), 0.3, tolerance = 0.06)
})

test_that("an F walk and an ar1 bridge a closed year at the elapsed years", {
  is_est <- c(TRUE, TRUE, FALSE, TRUE) # year three closed
  set.seed(8)
  rw <- replicate(20000, draw_F_dev_series(c(1.5, 0, 0, 0), is_est, PE_model = 2, sigma = 0.2, rho = 0, n_cond = 0))
  expect_true(all(rw[1,] == 1.5)) # the diffuse first year keeps its value, a level the data set
  expect_equal(mean(rw[2,]), 1.5, tolerance = 0.01) # and the walk steps on from it
  expect_equal(sd(rw[2,]), 0.2, tolerance = 0.03)
  expect_equal(sd(rw[4,] - rw[2,]), 0.2 * sqrt(2), tolerance = 0.03)
  expect_true(all(rw[3,] == 0)) # a closed year is not drawn
  ar <- replicate(20000, draw_F_dev_series(rep(0, 4), is_est, PE_model = 3, sigma = 0.2, rho = 0.6, n_cond = 0))
  stationary <- 0.2 / sqrt(1 - 0.6^2)
  expect_equal(sd(ar[1,]), stationary, tolerance = 0.03)
  expect_equal(sd(ar[4,]), stationary, tolerance = 0.03) # a stationary series stays stationary across the gap
  expect_equal(cor(ar[2,], ar[4,]), 0.6^2, tolerance = 0.05)
})

test_that("a level shared over bins takes one drawn value", {
  map <- array(NA, dim = c(1, 5, 3, 1, 1))
  map[1,,1,1,1] <- 1:5
  map[1,,2,1,1] <- 1:5 # bin two shares bin one's series
  map[1,,3,1,1] <- 6:10
  set.seed(2)
  devs <- draw_sel_dev_surface(PE_codes = 1, map = map, pe_pars = array(log(0.3), dim = c(1, 3, 1, 1)),
                               fit_devs = array(0, dim = dim(map)), bins = 1:3, n_cond = 0)
  expect_equal(devs[1,,1,1,1], devs[1,,2,1,1])
  expect_false(isTRUE(all.equal(devs[1,,1,1,1], devs[1,,3,1,1])))
})

test_that("a level shared over regions, sexes and fleets takes one drawn value at the sd its first cell reads", {
  il <- suppressMessages(sweep_input(dims = list(n_regions = 2, n_sexes = 2, n_fish_fleets = 2, n_srv_fleets = 1),
    fishsel = list(cont_tv_fish_sel = c("iid_Fleet_1", "iid_Fleet_2"), fishsel_pe_pars_spec = c("est_all", "est_all"),
                   fish_sel_devs_spec = c("est_shared_r_s", "est_shared_f_1"))))
  map <- il$data$map_ln_fishsel_devs
  expect_identical(map[,,,,1], map[,,,,2]) # fleet two holds fleet one's levels
  pe_pars <- il$par$fishsel_pe_pars
  pe_pars[,1,,] <- log(0.3) # first logistic parameter
  pe_pars[,2,,] <- log(0.2) # second
  set.seed(6)
  draws <- replicate(200, draw_sel_dev_surface(PE_codes = il$data$cont_tv_fish_sel, map = map, pe_pars = pe_pars,
                                               fit_devs = array(0, dim = dim(map)), bins = il$data$fishsel_devs_min_shared_bins,
                                               n_cond = 0, pe_wt = il$data$fishsel_pe_wt))
  levels <- sort(unique(as.vector(map)))
  spread <- sapply(levels, function(level) max(apply(draws, 6, function(d) diff(range(d[which(map == level)])))))
  expect_equal(max(spread), 0)
  expect_equal(sd(as.vector(draws[1,,1,1,1,])), 0.3, tolerance = 0.05)
  expect_equal(sd(as.vector(draws[1,,2,1,1,])), 0.2, tolerance = 0.05)
})

test_that("selectivity at length is drawn at length and read through the size-age key, as the estimation model reads it", {
  il <- suppressWarnings(suppressMessages(objective_setup_input(n_lens = 12,
    fishsel = list(fish_selex_type = "length", cont_tv_fish_sel = "iid_Fleet_1", fishsel_pe_pars_spec = "fix", fish_sel_devs_spec = "est_all"))))
  il$par$fishsel_pe_pars[1,1:2,,1] <- log(0.25)
  built <- fleet_devs_sim_list(il, random = "ln_fishsel_devs", n_sims = 200)
  expect_true(built$sim_list$fleet_devs_on[["fish"]])
  set.seed(3)
  env <- suppressMessages(Setup_sim_env(built$sim_list))
  expect_equal(sd(env$ln_fishsel_devs[1,,1,1,1,]), 0.25, tolerance = 0.05)

  # the estimation model at replicate 5's deviations gives its selectivity at age, which differs from the fit's
  at_draw <- built$obj$env$parList(par = built$obj$env$par)
  at_draw$ln_fishsel_devs[] <- env$ln_fishsel_devs[,,,,,5]
  rep_5 <- fit_model(il$data, at_draw, il$map, random = "ln_fishsel_devs", do_optim = FALSE, silent = TRUE)$rep
  n_yrs <- length(il$data$years)
  expect_equal(as.vector(rep_5$fish_sel[,,1:n_yrs,,,,,drop = FALSE]), as.vector(env$fish_sel[,,1:n_yrs,,,,,5]), tolerance = 1e-12)
  expect_gt(max(abs(as.vector(built$obj$rep$fish_sel[,,1:n_yrs,,,,,drop = FALSE]) - as.vector(env$fish_sel[,,1:n_yrs,,,,,5]))), 0.01)
})

test_that("a penalty weighted by w is drawn at the sd divided by the square root of w", {
  map <- array(1:60, dim = c(1, 30, 2, 1, 1))
  pe_pars <- array(log(0.4), dim = c(1, 2, 1, 1))
  set.seed(7)
  draws <- replicate(300, draw_sel_dev_surface(PE_codes = 1, map = map, pe_pars = pe_pars, fit_devs = array(0, dim = dim(map)),
                                               bins = 1:2, n_cond = 0, pe_wt = 4))
  expect_equal(sd(as.vector(draws)), 0.4 / 2, tolerance = 0.03)
  expect_equal(pe_pars_at_weight(array(log(0.4), dim = c(1, 1, 4, 1)), 3, 4)[1,1,4,1], log(0.4) - log(4)) # a 3D GMRF's log variance
})

test_that("a closed loop keeps the fit's selectivity and draws its deviations on through the projection", {

  il <- fleet_devs_input() # iid deviations at sds 0.3 and 0.2 on the two logistic parameters
  obj <- suppressWarnings(fit_model(il$data, il$par, il$map, random = "ln_fishsel_devs", do_optim = FALSE, silent = TRUE))
  n_fit <- length(il$data$years)
  n_proj <- 6
  cl <- suppressWarnings(suppressMessages(condition_closed_loop_simulations(
    closed_loop_yrs = n_proj,
    n_sims = 300,
    obj$data,
    obj$parameters,
    obj$mapping,
    sd_rep = list(par.fixed = obj$par, par.random = obj$env$par[obj$env$random]),
    rep = obj$rep,
    random = "ln_fishsel_devs"
  )))
  expect_equal(unname(cl$fleet_devs_on[c("F", "fish")]), c(FALSE, TRUE)) # F is the control rule's past the fit

  set.seed(21)
  env <- suppressMessages(Setup_sim_env(cl))

  # the fitted years are the fit's, every replicate
  for(sim in c(1, 300)) expect_equal(as.vector(env$fish_sel[,,1:n_fit,,,,1,sim]), as.vector(obj$rep$fish_sel[,,1:n_fit,,,,1]), tolerance = 1e-12)

  # the projection years are drawn at the fitted sds, independent of the last fitted year under iid
  proj_devs <- env$ln_fishsel_devs_proj # region, projection year, slot, sex, fleet, replicate
  expect_equal(dim(proj_devs)[2], n_proj)
  expect_equal(sd(proj_devs[1,,1,1,1,]), 0.3, tolerance = 0.06)
  expect_equal(sd(proj_devs[1,,2,1,1,]), 0.2, tolerance = 0.06)

  # and the selectivity they give varies by replicate where it was held at the last fitted year before
  proj_yrs <- n_fit + seq_len(n_proj)
  sel_age3 <- env$fish_sel[1,1,proj_yrs,1,3,1,1,]
  expect_gt(sd(as.vector(sel_age3)), 0.01)
  expect_true(all(is.finite(env$fish_sel[,,proj_yrs,,,,,])))
})
