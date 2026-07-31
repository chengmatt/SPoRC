library(SPoRC)
library(testthat)

test_that("get_osa(model = ..., pop = TRUE) + plot_resids() work for a multi-sex, joint-sex, population-specific model", {

  set.seed(123)
  sim_list <- Setup_Sim_Dim(
    n_sims        = 1,
    n_yrs         = 10,
    n_regions     = 2,
    n_ages        = 8,
    n_lens        = NULL,
    n_sexes       = 2,  
    n_fish_fleets = 1,
    n_srv_fleets  = 1,
    n_seas        = 2,
    n_pop         = 3,
    natal_region  = c(1, 1, 2)   # pops 1 & 2 natal to region 1; pop 3 natal to region 2
  )

  sim_list <- Setup_Sim_Containers(sim_list)

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(
      n = sim_list$n_sims,
      {
        arr <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                 sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))
        for (r in 1:sim_list$n_regions)
          for (y in 1:sim_list$n_yrs)
            for (s in 1:sim_list$n_sexes) {
              for(p in 1:sim_list$n_pop) {
                for(seas in 1:sim_list$n_seas) {
                  arr[p,r, y,seas, , s, 1] <-  1 / (1 + exp(-1.5 * (1:sim_list$n_ages - 3)))
                }
              }
            }
        arr
      }
    ),
    Fmort_input = {
      n = sim_list$n_yrs * sim_list$n_seas * sim_list$n_sims * sim_list$n_fish_fleets
      t = seq(0, 2*pi, length.out = n)
      arr <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                               sim_list$n_fish_fleets, sim_list$n_sims))
      arr[1,,,,] <- 0.15 * exp(sin(t) + rnorm(n, 0, 0.1))
      arr[2,,,,] <- 0.05 * exp(-sin(t) + rnorm(n, 0, 0.1))
      arr
    },
    ISS_FishAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                          sim_list$n_fish_fleets, sim_list$n_sims)),
    ISS_FishAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                      sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                      sim_list$n_fish_fleets, sim_list$n_sims)),
    ln_sigmaC = array(log(0.01), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
    ln_sigmaC_pop = array(log(0.01), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(
      n = sim_list$n_sims,
      {
        arr <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                 sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))
        for (r in 1:sim_list$n_regions)
          for (y in 1:sim_list$n_yrs)
            for (s in 1:sim_list$n_sexes) {
              for(p in 1:sim_list$n_pop) {
                for(seas in 1:sim_list$n_seas) {
                  arr[p,r, y,seas, , s, 1] <-  1 / (1 + exp(-1 * (1:sim_list$n_ages - 2.5)))
                }
              }
            }
        arr
      }
    ),
    ISS_SrvAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                         sim_list$n_srv_fleets, sim_list$n_sims)),
    ISS_SrvAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                     sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                     sim_list$n_srv_fleets, sim_list$n_sims)),
    ObsSrvIdx_SE = array(0.15, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets)),
    ObsSrvIdx_pop_SE = array(0.15, dim = c(sim_list$n_pop, sim_list$n_regions,
                                           sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))
  )

  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,
    natmort_input = array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions,
                                       sim_list$n_yrs, sim_list$n_ages,
                                       sim_list$n_sexes, sim_list$n_sims)),
    WAA_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
            each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
            times = sim_list$n_sexes),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)
      )
    ),
    WAA_fish_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
            each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
            times = sim_list$n_sexes * sim_list$n_fish_fleets),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets)
      )
    ),
    WAA_srv_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
            each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
            times = sim_list$n_sexes * sim_list$n_srv_fleets),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets)
      )
    ),
    MatAA_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(1 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
            each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
            times = sim_list$n_sexes),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)
      )
    )
  )

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)

  sim_list$Movement <- array(
    0,
    dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
            sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
            sim_list$n_sexes, sim_list$n_sims)
  )

  stay_prob <- c(0.7, 0.3, 0.7)
  disperse_prob <- (1 - stay_prob) / (sim_list$n_regions - 1)
  non_natal_rate <- 0.15

  for (p in seq_len(sim_list$n_pop)) {
    nr <- sim_list$natal_region[p]
    for (r_from in seq_len(sim_list$n_regions)) {
      for (r_to in seq_len(sim_list$n_regions)) {
        prob <- if (r_to == r_from) stay_prob[p] else disperse_prob[p]
        sim_list$Movement[p, r_from, r_to, , 1, , , ] <- prob
      }
      for (r_to in seq_len(sim_list$n_regions)) {
        prob <- if (r_to == nr) 1 - non_natal_rate * (sim_list$n_regions - 1) else non_natal_rate
        sim_list$Movement[p, r_from, r_to, , 2, , , ] <- prob
      }
    }
  }

  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = {
      R0_arry <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
      R0_arry[1, 1, , ] <- 7
      R0_arry[2, 1, , ] <- 7
      R0_arry[3, 2, , ] <- 7
      R0_arry
    },
    ln_sigmaR = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
    init_age_strc = "matrix",
    recruitment_opt = 'bh_rec',
    rec_dd = 'local',
    spawn_seas = 2,
    t_spawn = 0.5,
    h_input = {
      h_arry <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
      h_arry[1, 1, , ] <- 0.75
      h_arry[2, 1, , ] <- 0.75
      h_arry[3, 2, , ] <- 0.75
      h_arry
    }
  )

  sim_obj <- suppressWarnings(Simulate_Pop_Static(sim_list = sim_list, output_path = NULL))
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)

  # Sanity: pop-specific joint-sex comps should sum to ~ISS (one joint draw
  # across both sexes), confirming the simulation side actually exercises
  # the ct==2 code path the two bugs above touched.
  expect_equal(
    sum(sim_data$ObsFishAgeComps_pop[1,1,5,1,,,1]),
    sim_data$ISS_FishAgeComps_pop[1,1,5,1,1,1]
  )

  input_list <- Setup_Mod_Dim(
    years = 1:sim_obj$n_years,
    ages = 1:sim_obj$n_ages,
    lens = sim_obj$n_lens,
    n_regions = sim_obj$n_regions,
    n_sexes = sim_obj$n_sexes,
    n_fish_fleets = sim_obj$n_fish_fleets,
    n_srv_fleets = sim_obj$n_srv_fleets,
    n_seas = sim_obj$n_seas,
    n_pop = sim_obj$n_pop,
    seasdur = sim_obj$seasdur,
    natal_region = c(1,1,2),
    verbose = FALSE,
    do_internal_comp_osa = TRUE
  )

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0,
    sigmaR_switch = 1,
    init_age_strc = "matrix",
    equil_init_age_strc = "stoch_all",
    spawn_seas = sim_obj$spawn_seas,
    t_spawn = sim_obj$t_spawn,
    rec_model = "bh_rec",
    sigmaR_spec = "fix",
    rec_dd = 'local',
    InitDevs_spec = "est_shared_r",
    RecDevs_spec = "est_shared_r",
    sexratio_spec = "fix",
    rec_region_prop_spec = 'no_dispersal',
    h_spec = 'fix',
    steepness_h = {
      h_arry <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions))
      h_arry[1,1] <- qlogis((0.75 - 0.2) / 0.8)
      h_arry[2,1] <- qlogis((0.75 - 0.2) / 0.8)
      h_arry[3,2] <- qlogis((0.75 - 0.2) / 0.8)
      h_arry
    },
    ln_sigmaR = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
    ln_global_R0 = array(c(log(7), log(5), log(10)), dim = input_list$data$n_pop)
  )

  input_list <- suppressWarnings(Setup_Mod_Biologicals(
    input_list = input_list,
    WAA = sim_data$WAA,
    MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0,
    AgeingError = sim_data$AgeingError,
    M_spec = "est_ln_M"
  ))

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)

  input_list <- Setup_Mod_Movement(
    input_list = input_list,
    do_recruits_move = 0,
    use_fixed_movement = 0,
    Movement_popblk_spec = list(1,2,3),
    Movement_seasblk_spec = list(1,2)
  )

  input_list <- suppressWarnings(Setup_Mod_Catch_and_F(
    input_list = input_list,
    ObsCatch = sim_data$ObsCatch,
    UseCatch = array(0, dim = dim(sim_data$UseCatch)),
    ObsCatch_pop = sim_data$ObsCatch_pop,
    UseCatch_pop = sim_data$UseCatch_pop,
    Use_F_pen = 1,
    sigmaC_spec = "fix",
    sigmaC_pop_spec = 'fix',
    ln_sigmaC = sim_data$ln_sigmaC,
    ln_sigmaC_pop = sim_data$ln_sigmaC_pop,
    ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions,
                                      input_list$data$n_seas,
                                      input_list$data$n_fish_fleets))
  ))

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list,
    ObsFishIdx = sim_data$ObsFishIdx,
    ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
    ObsFishAgeComps = sim_data$ObsFishAgeComps,
    ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = array(0, dim = dim(sim_data$UseFishAgeComps)),
    UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
    ISS_FishLenComps = sim_data$ISS_FishLenComps,
    ObsFishAgeComps_pop = sim_data$ObsFishAgeComps_pop,
    UseFishAgeComps_pop = sim_data$UseFishAgeComps_pop,
    ISS_FishAgeComps_pop = sim_data$ISS_FishAgeComps_pop,
    FishAgeComps_pop_LikeType = c("Multinomial"),
    FishAgeComps_pop_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    fish_idx_type = 'none',
    FishAgeComps_LikeType = c("Multinomial"),
    FishLenComps_LikeType = c("none"),
    FishAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list,
    ObsSrvIdx = sim_data$ObsSrvIdx,
    ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = array(0, dim = dim(sim_data$UseSrvIdx)),
    ObsSrvIdx_pop = sim_data$ObsSrvIdx_pop,
    ObsSrvIdx_pop_SE = sim_data$ObsSrvIdx_pop_SE,
    UseSrvIdx_pop = array(1, dim = dim(sim_data$UseSrvIdx_pop)),
    ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
    ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = array(0, dim = dim(sim_data$UseSrvAgeComps)),
    UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
    ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    ObsSrvAgeComps_pop = sim_data$ObsSrvAgeComps_pop,
    UseSrvAgeComps_pop = sim_data$UseSrvAgeComps_pop,
    ISS_SrvAgeComps_pop = sim_data$ISS_SrvAgeComps_pop,
    SrvAgeComps_pop_LikeType = c("Multinomial"),
    SrvAgeComps_pop_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    srv_idx_type = c("biom"),
    SrvAgeComps_LikeType = c("Multinomial"),
    SrvLenComps_LikeType = c("none"),
    SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
    SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list,
    fish_sel_model = c("logist1_Fleet_1"),
    fish_fixed_sel_pars_spec = c("est_shared_r"),
    fish_q_spec = c("fix")
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list,
    srv_sel_model = c("logist1_Fleet_1"),
    srv_fixed_sel_pars_spec = c("est_shared_r"),
    srv_q_spec = c("est_shared_r")
  )

  input_list <- Setup_Mod_Weighting(
    input_list = input_list,
    Wt_Catch = 1,
    Wt_Catch_pop = 1,
    Wt_SrvIdx_pop = 1,
    Wt_FishIdx = 1,
    Wt_SrvIdx = 1,
    Wt_Rec = 1,
    Wt_F = 1,
    Wt_FishAgeComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                       input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                      input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_FishAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                           input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                          input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
  )

  model <- fit_model(
    input_list$data, input_list$par, input_list$map,
    random = NULL, silent = TRUE, do_optim = TRUE, newton_loops = 5
  )
  expect_true(RTMB::sdreport(model)$pdHess)
  expect_jnLL_decomposes(model)

  osa_pop_fish <- get_osa(model = model, data = input_list$data, comp_source = "FishAge", pop = TRUE,
                          family = "discrete", bins = input_list$data$ages, bin_label = "Age")
  expect_equal(unique(osa_pop_fish$res$comp_type), "SpltR_JntS")
  expect_setequal(unique(osa_pop_fish$res$region), 1:2)
  sdnr_fish <- sd(osa_pop_fish$res$resid)
  expect_true(is.finite(sdnr_fish))
  expect_gt(sdnr_fish, 0.7)
  expect_lt(sdnr_fish, 1.3)

  osa_pop_srv <- get_osa(model = model, data = input_list$data, comp_source = "SrvAge", pop = TRUE,
                         family = "discrete", bins = input_list$data$ages, bin_label = "Age")
  expect_equal(unique(osa_pop_srv$res$comp_type), "SpltR_JntS")
  sdnr_srv <- sd(osa_pop_srv$res$resid)
  expect_true(is.finite(sdnr_srv))
  expect_gt(sdnr_srv, 0.7)
  expect_lt(sdnr_srv, 1.3)

  p1 <- plot_resids(osa_pop_fish)
  expect_s3_class(p1[[1]], "ggplot")
  expect_s3_class(p1[[2]], "ggplot")

  p2 <- plot_resids(osa_pop_srv)
  expect_s3_class(p2[[1]], "ggplot")
  expect_s3_class(p2[[2]], "ggplot")

  # index-type sources (pop = TRUE): regression check for a dplyr data-mask
  # shadowing bug where `map`'s own "pop" column shadowed the `pop`
  # function argument inside `if(pop)`, breaking every pop = TRUE call
  for(idx_src in c("Catch", "SrvIdx")) {
    osa_idx_pop <- get_osa(model = model, data = input_list$data, index_source = idx_src, pop = TRUE)
    expect_false(is.null(osa_idx_pop))
    res_idx <- osa_idx_pop$res
    expect_true(all(c("fleet","region","year","season","pop","resid","idx_type") %in% names(res_idx)))
    expect_equal(unique(res_idx$idx_type), idx_src)
    expect_setequal(unique(res_idx$pop), 1:3)
    expect_true(all(is.finite(res_idx$resid)))

    p_idx <- plot_resids(osa_idx_pop)
    expect_s3_class(p_idx[[1]], "ggplot")
    expect_s3_class(p_idx[[2]], "ggplot")
    # force lazy facet evaluation (the bug above only manifested at render time)
    ggplot2::ggplot_build(p_idx[[1]])
    ggplot2::ggplot_build(p_idx[[2]])
  }
})
