# A dsem written from scratch on catchability: the drawn series lands in the fleet's deviations,
# scales the catchability the index is formed at, and is refused where the operating model has nothing to read it.

# a one region, one fleet operating model with a catchability series and nothing else linked
scratch_q_om <- function(process = "srv_q", n_sims = 2, n_yrs = 20, n_ages = 6, seed = 321) {

  logistic <- function(slope, infl) 1 / (1 + exp(-slope * ((1:n_ages) - infl)))
  yearly <- function(v, dims) replicate(n_sims, array(rep(v, each = n_yrs), dim = dims))

  sim_list <- Setup_Sim_Dim(n_sims = n_sims, n_yrs = n_yrs, n_regions = 1, n_ages = n_ages, n_lens = NULL,
                            n_sexes = 1, n_fish_fleets = 1, n_srv_fleets = 1, n_pop = 1)
  sim_list <- Setup_Sim_Containers(sim_list)
  sim_list <- Setup_Sim_Fishing(sim_list, fish_sel_input = yearly(logistic(3, 3), c(1, 1, n_yrs, 1, n_ages, 1, 1)),
                                ln_sigmaC = array(log(0.02), dim = c(1, n_yrs, 1, 1)),
                                ObsFishIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)),
                                ISS_FishAgeComps = array(200, dim = c(1, n_yrs, 1, 1, 1, n_sims)))
  sim_list <- Setup_Sim_Survey(sim_list, srv_sel_input = yearly(logistic(1, 2), c(1, 1, n_yrs, 1, n_ages, 1, 1)),
                               ObsSrvIdx_SE = array(0.2, dim = c(1, n_yrs, 1, 1)),
                               ISS_SrvAgeComps = array(200, dim = c(1, n_yrs, 1, 1, 1, n_sims)))
  sim_list <- suppressWarnings(Setup_Sim_Biologicals(sim_list,
                                                     natmort_input = replicate(n_sims, array(0.3, dim = c(1, 1, n_yrs, n_ages, 1))),
                                                     WAA_input = yearly(5 * logistic(3, 2), c(1, 1, n_yrs, 1, n_ages, 1)),
                                                     WAA_fish_input = yearly(5 * logistic(3, 2), c(1, 1, n_yrs, 1, n_ages, 1, 1)),
                                                     WAA_srv_input = yearly(5 * logistic(3, 2), c(1, 1, n_yrs, 1, n_ages, 1, 1)),
                                                     MatAA_input = yearly(logistic(3, 2), c(1, 1, n_yrs, 1, n_ages, 1))))
  sim_list <- Setup_Sim_Tagging(sim_list, use_conv_fish_tagging = 0)
  sim_list$Movement <- array(1, dim = c(1, 1, 1, n_yrs, 1, n_ages, 1, n_sims))
  sim_list <- Setup_Sim_Rec(sim_list, R0_input = replicate(n_sims, array(5, dim = c(1, 1, n_yrs))),
                            rinit_input = array(2, dim = c(1, 1, n_sims)), use_rinit = 1,
                            ln_sigmaR = array(log(0.6), dim = c(2, 1, 1)), recruitment_opt = "mean_rec", init_age_strc = 1,
                            ln_InitDevs_input = array(0, dim = c(1, 1, n_ages - 1, 1, n_sims)))

  list(sim_list = sim_list,
       arrows = c("env -> env, 1, rho_env", "env <-> env, 0, sd_env",
                  paste0("env -> ", process, ", 0, b_env"), paste0(process, " <-> ", process, ", 0, sd_q")),
       values = c(rho_env = 0.6, sd_env = 1, b_env = 0.3, sd_q = 0.15),
       seed = seed)

}

test_that("a scratch catchability series reaches the fleet's catchability", {

  for(process in c("srv_q", "fish_q")) {

    om <- scratch_q_om(process)
    prefix <- sub("_q$", "", process)
    sim_list <- do.call(Setup_Sim_q_devs, stats::setNames(list(om$sim_list, "dsem"), c("sim_list", paste0(prefix, "_q_model"))))
    sim_list <- Setup_Sim_DSEM(sim_list, dsem_arrows = om$arrows, dsem_values = om$values, dsem_processes = process)

    set.seed(om$seed)
    sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
    devs <- sim_obj[[paste0("ln_", process, "_devs")]]
    q <- sim_obj[[process]]

    # the drawn column is the fleet's deviations, and catchability is the mean of one scaled by them
    expect_equal(as.numeric(devs[1,,1,1]), as.numeric(sim_obj$dsem_x_sim[,2,1]), tolerance = 1e-12)
    expect_equal(as.numeric(log(q[1,,1,1])), as.numeric(sim_obj$dsem_x_sim[,2,1]), tolerance = 1e-12)
    expect_gt(stats::sd(devs[1,,1,1]), 0.1) # the series varies rather than sitting at its mean
    expect_gt(max(abs(q[1,,1,1] - q[1,,1,2])), 0) # and each replicate draws its own

  } # end process loop

})

test_that("a scratch catchability series is refused where nothing would read it", {

  om <- scratch_q_om("srv_q")

  # the fleet was never given a catchability model, so draw_sim_q_devs would pass over the drawn cells
  expect_error(Setup_Sim_DSEM(om$sim_list, dsem_arrows = om$arrows, dsem_values = om$values, dsem_processes = "srv_q"),
               "holds no srv_q_model")

  # growth still needs a fit, since the operating model builds no deviations of its own for it
  expect_error(Setup_Sim_DSEM(om$sim_list, dsem_arrows = om$arrows, dsem_values = om$values, dsem_processes = "growth"),
               "growth and movement need a fit")

})
