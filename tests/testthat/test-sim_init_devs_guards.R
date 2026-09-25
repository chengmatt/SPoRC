# Checks the operating model says when it draws initial age deviations, centers them and recruitment on the bias
# correction switch, and that the estimation setup notes rinit against a deviation on every initial age.

init_devs_sim <- function(rec_bias_correct = 1, ln_InitDevs_input = NULL, n_sims = 200, seed = 5) {
  sim_list <- Setup_Sim_Dim(n_sims = n_sims, n_yrs = 3, n_regions = 1, n_ages = 6, n_lens = NULL, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  yearly <- function(v, dims) replicate(n_sims, array(rep(v, each = 3), dim = dims))
  sim_list <- Setup_Sim_Fishing(sim_list, fish_sel_input = yearly(rep(1, 6), c(1, 1, 3, 1, 6, 1, 1)))
  sim_list <- Setup_Sim_Survey(sim_list, srv_sel_input = yearly(rep(1, 6), c(1, 1, 3, 1, 6, 1, 1)))
  sim_list <- suppressWarnings(Setup_Sim_Biologicals(sim_list, natmort_input = replicate(n_sims, array(0.3, dim = c(1, 1, 3, 6, 1))),
                                                     WAA_input = yearly(1:6, c(1, 1, 3, 1, 6, 1)), WAA_fish_input = yearly(1:6, c(1, 1, 3, 1, 6, 1, 1)),
                                                     WAA_srv_input = yearly(1:6, c(1, 1, 3, 1, 6, 1, 1)), MatAA_input = yearly(rep(1, 6), c(1, 1, 3, 1, 6, 1))))
  sim_list <- Setup_Sim_Tagging(sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, 3, 1, 6, 1, n_sims))
  sim_list <- Setup_Sim_Rec(sim_list, R0_input = replicate(n_sims, array(5, dim = c(1, 1, 3))), rinit_input = array(5, dim = c(1, 1, n_sims)), use_rinit = 1,
                            ln_sigmaR = array(log(1), dim = c(2, 1, 1)), recruitment_opt = "mean_rec", init_age_strc = 1,
                            ln_InitDevs_input = ln_InitDevs_input, rec_bias_correct = rec_bias_correct)
  set.seed(seed)
  Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
}

test_that("the setup says when the initial deviations are drawn, and zeros keep the population in equilibrium", {

  expect_message(om <- init_devs_sim(), "draws its own initial age deviations")
  expect_gt(sd(om$ln_InitDevs), 0.5) # drawn at sigma 1
  expect_silent(om0 <- suppressWarnings(init_devs_sim(ln_InitDevs_input = array(0, dim = c(1, 1, 5, 1, 200)))))
  expect_true(all(om0$ln_InitDevs == 0))

})

test_that("the bias correction switch centers the initial deviations and recruitment", {

  om1 <- suppressMessages(init_devs_sim(rec_bias_correct = 1))
  om0 <- suppressMessages(init_devs_sim(rec_bias_correct = 0))
  # 200 replicates by 5 ages of N(center, 1): the mean sits at the center within 0.1
  expect_equal(mean(om1$ln_InitDevs), -0.5, tolerance = 0.1)
  expect_equal(mean(om0$ln_InitDevs), 0, tolerance = 0.1)
  # recruitment is R0 exp(dev) under 0 and R0 exp(dev - 1/2) under 1
  expect_equal(om0$Rec[1,1,2,], 5 * exp(om0$ln_RecDevs[1,1,2,]), tolerance = 1e-10)
  expect_equal(om1$Rec[1,1,2,], 5 * exp(om1$ln_RecDevs[1,1,2,] - 0.5), tolerance = 1e-10)
  expect_error(suppressMessages(init_devs_sim(rec_bias_correct = 2)), "0 or 1")

})

test_that("the estimation setup notes rinit against a deviation on every initial age", {

  il <- Setup_Mod_Dim(years = 1:10, ages = 1:6, lens = NULL, n_regions = 1, n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1, natal_region = 1, verbose = FALSE)
  rec <- function(equil) Setup_Mod_Rec(il, do_rec_bias_ramp = 0, sigmaR_switch = 1, ln_sigmaR = array(log(1), c(2, 1, 1)), rec_model = "mean_rec", use_rinit = 1,
                                       sigmaR_spec = "fix", init_age_strc = 1, equil_init_age_strc = equil, ln_global_R0 = log(5), ln_rinit = log(2))
  rec(2) # every initial age has a deviation, and rinit is estimated beside them
  expect_true(any(grepl("separated only by the initial deviation penalty", messages_list)))
  rec(0) # no deviations at all, nothing to note
  expect_false(any(grepl("separated only by the initial deviation penalty", messages_list)))

})
