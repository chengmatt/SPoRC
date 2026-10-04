# The annual cycle writes into the environment it is given, whatever the caller names it, so two
# operating models built from one seed and run side by side stay identical.

cycle_env_sim_list <- function() {
  n_yrs <- naaom_cfg$n_yrs
  n_ages <- naaom_cfg$n_ages
  curve <- function(slope, infl, scale = 1)
    array(rep(scale / (1 + exp(-slope * ((1:n_ages) - infl))), each = n_yrs),
          dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))

  sl <- Setup_Sim_Dim(n_sims = 1, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = NULL,
                      n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sl <- Setup_Sim_Containers(sl)
  sl <- Setup_Sim_Fishing(sl,
                          fish_sel_input = replicate(1, curve(3, 2)),
                          ret_sel_input = replicate(1, curve(3, 2)),
                          dmr_input = array(0, dim = c(1, n_yrs, 1, 1, 1)),
                          Fmort_input = array(naaom_cfg$f_ramp, dim = c(1, n_yrs, 1, 1, 1)),
                          ISS_FishAgeComps = array(naaom_cfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)))
  sl <- Setup_Sim_Survey(sl,
                         srv_sel_input = replicate(1, curve(1, 3)),
                         ObsSrvIdx_SE = array(naaom_cfg$idx_se, dim = c(1, n_yrs, 1, 1)),
                         ISS_SrvAgeComps = array(naaom_cfg$comp_iss, dim = c(1, n_yrs, 1, 1, 1, 1)))
  biol <- function(v) array(rep(v, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1))
  suppressWarnings(sl <- Setup_Sim_Biologicals(sl,
    natmort_input = replicate(1, array(naaom_cfg$M, dim = c(1, 1, n_yrs, n_ages, 1))),
    WAA_input = replicate(1, biol(naaom_cfg$waa)),
    WAA_fish_input = replicate(1, array(rep(naaom_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    WAA_srv_input = replicate(1, array(rep(naaom_cfg$waa, each = n_yrs), dim = c(1, 1, n_yrs, 1, n_ages, 1, 1))),
    MatAA_input = replicate(1, biol(naaom_cfg$mat))))
  sl <- Setup_Sim_Tagging(sl, use_conv_fish_tagging = 0)
  sl$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, 1))

  Setup_Sim_Rec(sl,
                R0_input = replicate(1, array(5, dim = c(1, 1, n_yrs))),
                ln_sigmaR = array(log(naaom_cfg$sigmaR), dim = c(2, 1, 1)),
                recruitment_opt = "mean_rec",
                init_age_strc = 1)
}

test_that("the annual cycle runs on an environment not named sim_env", {

  sim_list <- cycle_env_sim_list()
  n_yrs <- sim_list$n_yrs

  set.seed(11)
  static <- suppressMessages(Simulate_Pop_Static(sim_list = sim_list))

  # the same operating model by hand, under another name. This used to stop with
  # "object 'sim_env' not found", or write into whatever sim_env the caller happened to have
  set.seed(11)
  om <- suppressMessages(Setup_sim_env(sim_list))
  for(y in 1:n_yrs) run_annual_cycle(y, 1, om)

  expect_true(all(om$NAA[1, 1, 1:n_yrs, 1, , 1, 1] > 0)) # every year was filled in
  expect_identical(om$NAA, static$NAA)
  expect_identical(om$TrueCatch, static$TrueCatch)

})

test_that("two operating models from one seed run side by side stay identical", {

  sim_list <- cycle_env_sim_list()
  n_yrs <- sim_list$n_yrs

  # one of them named sim_env, as in a closed loop script that keeps an earlier run around
  set.seed(12)
  sim_env <- suppressMessages(Setup_sim_env(sim_list))
  set.seed(12)
  env_old <- suppressMessages(Setup_sim_env(sim_list))

  # alternate years between the two, reseeding so both see the same draws in a year. The second
  # used to write into the first, then stop with "NA in probability vector" on its own empty catch
  for(y in 1:n_yrs) {
    set.seed(100 + y)
    run_annual_cycle(y, 1, sim_env)
    set.seed(100 + y)
    run_annual_cycle(y, 1, env_old)
  } # end y loop

  expect_true(all(env_old$NAA[1, 1, 1:n_yrs, 1, , 1, 1] > 0)) # so the comparisons below are not of two empty arrays
  expect_true(all(env_old$TrueCatch[1, 1:n_yrs, 1, 1, 1] > 0))
  expect_identical(env_old$TrueCatch, sim_env$TrueCatch)
  expect_identical(env_old$NAA, sim_env$NAA)

})

test_that("a copy of an operating model runs without touching the original", {

  sim_list <- cycle_env_sim_list()
  fork_yr <- 20

  set.seed(13)
  sim_env <- suppressMessages(Setup_sim_env(sim_list))
  for(y in 1:(fork_yr - 1)) run_annual_cycle(y, 1, sim_env)

  # copy at the fork and double the fishing in the copy's fork year
  env_copy <- list2env(as.list(sim_env, all.names = TRUE), envir = new.env(parent = parent.env(sim_env)))
  env_copy$Fmort[1, fork_yr, 1, 1, 1] <- 2 * sim_env$Fmort[1, fork_yr, 1, 1, 1]
  NAA_at_fork <- sim_env$NAA
  catch_at_fork <- sim_env$TrueCatch

  # the copy's run leaves the original as it was. The copy's sim_env used to point back at the original
  set.seed(200)
  run_annual_cycle(fork_yr, 1, env_copy)
  expect_identical(sim_env$NAA, NAA_at_fork)
  expect_identical(sim_env$TrueCatch, catch_at_fork)

  # then the original runs the same year on the same draws, and the two differ only through the fishing
  set.seed(200)
  run_annual_cycle(fork_yr, 1, sim_env)
  expect_identical(env_copy$NAA[1, 1, 1:fork_yr, 1, , 1, 1], sim_env$NAA[1, 1, 1:fork_yr, 1, , 1, 1]) # same start of year numbers
  expect_gt(env_copy$TrueCatch[1, fork_yr, 1, 1, 1], sim_env$TrueCatch[1, fork_yr, 1, 1, 1])
  expect_lt(sum(env_copy$NAA[1, 1, fork_yr + 1, 1, , 1, 1]), sum(sim_env$NAA[1, 1, fork_yr + 1, 1, , 1, 1]))

})
