# The operating model takes the fit's bias_correct_pe, bias_correct_oe and sigmaR_switch under both the self test and
# the closed loop, and draws discards and lognormal at-age data at their mean under bias_correct_oe = 1.

data("dusky_rtmb_model")

# the dusky fit with every switch off its default
switched_dusky <- function() {
  m <- dusky_rtmb_model
  m$data$bias_correct_pe <- 2 # "all"
  m$data$bias_correct_oe <- 1
  m$data$sigmaR_switch <- 5
  m
}

# the operating model list the self test builds, taken before any replicate runs
self_test_sim_list <- function(m) {
  captured <- new.env()
  testthat::with_mocked_bindings(
    Simulate_Pop_Static = function(sim_list, ...) { captured$sim_list <- sim_list; stop("captured") },
    try(suppressWarnings(suppressMessages(
      simulation_self_test(data = m$data, parameters = m$parameters, mapping = m$mapping, random = NULL, rep = m$rep,
                           sd_rep = m$sdrep, n_sims = 1, newton_loops = 0, what = "SSB"))), silent = TRUE),
    .package = "SPoRC")
  captured$sim_list
}

closed_loop_sim_list <- function(m, ...) {
  suppressWarnings(suppressMessages(condition_closed_loop_simulations(closed_loop_yrs = 2, n_sims = 1, m$data, m$parameters, m$mapping,
                                                                      sd_rep = m$sdrep, rep = m$rep, random = NULL, ...)))
}

test_that("the self test runs the operating model at the fit's centering and sigmaR switch", {
  sim_list <- self_test_sim_list(switched_dusky())
  expect_equal(sim_list$bias_correct_pe, 2)
  expect_equal(sim_list$bias_correct_oe, 1)
  expect_equal(sim_list$sigmaR_switch, 5)
})

test_that("a data list from before the switches existed gives the operating model's defaults", {
  m <- dusky_rtmb_model
  m$data$bias_correct_pe <- m$data$bias_correct_oe <- m$data$sigmaR_switch <- NULL
  sim_list <- self_test_sim_list(m)
  expect_equal(sim_list$bias_correct_pe, 1)
  expect_equal(sim_list$bias_correct_oe, 0)
  expect_equal(sim_list$sigmaR_switch, 1)
})

test_that("the closed loop takes the fit's switches unless given its own", {
  fit_switches <- closed_loop_sim_list(switched_dusky())
  expect_equal(fit_switches$bias_correct_pe, 2)
  expect_equal(fit_switches$bias_correct_oe, 1)
  expect_equal(fit_switches$sigmaR_switch, 5)

  own_switches <- closed_loop_sim_list(switched_dusky(), bias_correct_pe = "none", bias_correct_oe = 0, sigmaR_switch = 1)
  expect_equal(own_switches$bias_correct_pe, 0)
  expect_equal(own_switches$bias_correct_oe, 0)
  expect_equal(own_switches$sigmaR_switch, 1)
})

test_that("a lognormal at-age draw sits at its mean under bias_correct_oe = 1 and its median otherwise", {
  numbers <- array(c(100, 50, 25, 10), dim = c(1, 1, 4, 1))
  use <- array(1, dim = c(1, 4, 1))
  se <- array(0, dim = c(1, 4, 1))
  sigma <- 0.5
  log_ratio <- function(oe) {
    as.vector(replicate(4000, {
      drawn <- SPoRC:::sim_at_age_cell(numbers, numbers, use, se, array(log(sigma), dim = c(4, 1)),
                                       1, 0, 0, FALSE, 1, bias_correct_oe = oe)
      log(drawn$obs / drawn$true)
    }))
  }
  set.seed(5)
  expect_lt(abs(mean(log_ratio(1)) + sigma^2 / 2), 0.02) # se of the mean is 0.004
  expect_lt(abs(mean(log_ratio(0))), 0.02)
})

test_that("discards are drawn at their mean under bias_correct_oe = 1 and their median otherwise", {

  n_yrs <- 3
  n_ages <- 5
  n_sims <- 1000
  sigma_d <- 0.5

  discard_om <- function(oe) {
    sim_list <- Setup_Sim_Dim(n_sims = n_sims, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = NULL, n_sexes = 1,
                              n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1, bias_correct_oe = oe)
    sim_list <- Setup_Sim_Containers(sim_list)
    curve <- function(slope, infl) replicate(n_sims, array(rep(1 / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
                                                           dim = c(1, 1, n_yrs, 1, n_ages, 1, 1)))
    sim_list <- Setup_Sim_Fishing(
      sim_list = sim_list,
      fish_sel_input = curve(3, 1),
      ret_sel_input = curve(3, 3), # the youngest ages caught are mostly discarded
      dmr_input = array(0.5, dim = c(1, n_yrs, 1, 1, n_sims)),
      Fmort_input = array(0.2, dim = c(1, n_yrs, 1, 1, n_sims)),
      ln_sigmaD = array(log(sigma_d), dim = c(1, n_yrs, 1, 1)),
      ln_sigmaD_pop = array(log(sigma_d), dim = c(1, 1, n_yrs, 1, 1)),
      ISS_FishAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 1, n_sims))
    )
    sim_list <- Setup_Sim_Survey(sim_list = sim_list, srv_sel_input = curve(1, 2), ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)),
                                 ISS_SrvAgeComps = array(50, dim = c(1, n_yrs, 1, 1, 1, n_sims)))
    biol <- replicate(n_sims, array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1)))
    sim_list <- suppressWarnings(Setup_Sim_Biologicals(
      sim_list = sim_list,
      natmort_input = replicate(n_sims, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
      WAA_input = biol,
      WAA_fish_input = replicate(n_sims, array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
      WAA_srv_input = replicate(n_sims, array(1, dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
      MatAA_input = biol
    ))
    sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)
    sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, n_sims))
    sim_list <- suppressMessages(Setup_Sim_Rec(sim_list = sim_list, R0_input = replicate(n_sims, array(5, dim = c(1, 1, n_yrs))),
                                               ln_sigmaR = array(log(0.3), dim = c(2, 1, 1)), recruitment_opt = "mean_rec", init_age_strc = 1,
                                               ln_InitDevs_input = array(0, dim = c(1, 1, n_ages - 1, 1, n_sims))))
    suppressMessages(Simulate_Pop_Static(sim_list = sim_list, output_path = NULL))
  }

  set.seed(9)
  for(oe in c(1, 0)) {
    om <- discard_om(oe)
    expect_true(all(om$TrueDiscard > 0))
    centered <- if(oe == 1) -sigma_d^2 / 2 else 0
    expect_lt(abs(mean(log(om$ObsDiscard / om$TrueDiscard)) - centered), 0.04) # se of the mean is 0.009
    expect_lt(abs(mean(log(om$ObsDiscard_pop / om$TrueDiscard_pop)) - centered), 0.04)
  } # end oe loop

})
