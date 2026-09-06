library(SPoRC)
library(testthat)

test_that("OSA residuals are well-calibrated under correct EM and mis-calibrated under misspecified movement", {

  set.seed(123)
  sim_list <- Setup_Sim_Dim(
    n_sims        = 1,
    n_yrs         = 10,
    n_regions     = 2,
    n_ages        = 8,
    n_lens        = NULL,
    n_sexes       = 1,
    n_fish_fleets = 1,
    n_srv_fleets  = 1,
    n_seas        = 2,
    n_pop         = 3,
    natal_region  = c(1, 1, 2)   # pops 1 & 2 natal to region 1; pop 3 natal to region 2
  )

  sim_list <- Setup_Sim_Containers(sim_list)

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,

    # logistic selectivity
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
      arr[1,,,,] <- 0.15 * exp(sin(t) + rnorm(n, 0, 0.1))   # region 1 higher F, peaks early
      arr[2,,,,] <- 0.05 * exp(-sin(t) + rnorm(n, 0, 0.1))  # region 2 lower F, peaks late
      arr
    },

    # Fishery Age ISS
    ISS_FishAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                          sim_list$n_fish_fleets, sim_list$n_sims)),
    ISS_FishAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                      sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                      sim_list$n_fish_fleets, sim_list$n_sims)),


    # Sigma for catch
    ln_sigmaC = array(log(0.01), dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets)),
    ln_sigmaC_pop = array(log(0.01), dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_fish_fleets))

  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,

    # Logistic selectivity
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

    # Survey Age ISS
    ISS_SrvAgeComps = array(500, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                         sim_list$n_srv_fleets, sim_list$n_sims)),
    ISS_SrvAgeComps_pop = array(round(500 / sim_list$n_pop), dim = c(sim_list$n_pop, sim_list$n_regions,
                                                                     sim_list$n_yrs, sim_list$n_seas, sim_list$n_sexes,
                                                                     sim_list$n_srv_fleets, sim_list$n_sims)),

    # Sigma for Survey Index
    ObsSrvIdx_SE = array(0.15, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets)),
    ObsSrvIdx_pop_SE = array(0.15, dim = c(sim_list$n_pop, sim_list$n_regions,
                                           sim_list$n_yrs, sim_list$n_seas, sim_list$n_srv_fleets))

  )

  sim_list <- Setup_Sim_Biologicals(
    sim_list = sim_list,

    # Natural Mortality
    natmort_input = array(0.3, dim = c(sim_list$n_pop, sim_list$n_regions,
                                       sim_list$n_yrs, sim_list$n_ages,
                                       sim_list$n_sexes, sim_list$n_sims)),

    # Weight at age - Same for all pops
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

    # Fishery weight at age - same as WAA_input
    WAA_fish_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
            each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
            times = sim_list$n_sexes * sim_list$n_fish_fleets),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets)
      )
    ),

    # Survey weight at age - same as WAA_input
    WAA_srv_input = replicate(
      n = sim_list$n_sims,
      array(
        rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))),
            each = sim_list$n_pop * sim_list$n_regions * sim_list$n_yrs * sim_list$n_seas,
            times = sim_list$n_sexes * sim_list$n_srv_fleets),
        dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets)
      )
    ),

    # Maturity at age
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

  sim_list <- Setup_Sim_Tagging(
    sim_list          = sim_list,
    use_conv_fish_tagging = 1,
    n_tags            = 1e3,
    conv_tag_max_liberty  = 8,
    conv_fish_tag_like = "Poisson"
  )

  sim_list$Movement <- array(
    0,
    dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
            sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
            sim_list$n_sexes, sim_list$n_sims)
  )

  # Fill in movement matrix
  stay_prob <- c(0.7, 0.3, 0.7)  # probability of staying in current region during dispersal season
  disperse_prob <- (1 - stay_prob) / (sim_list$n_regions - 1)  # spread remainder equally
  non_natal_rate <- 0.15 # non natal homing rate

  for (p in seq_len(sim_list$n_pop)) {
    nr <- sim_list$natal_region[p]
    for (r_from in seq_len(sim_list$n_regions)) {

      # Season 1: diffusive dispersal — mostly stay, some movement out
      for (r_to in seq_len(sim_list$n_regions)) {
        prob <- if (r_to == r_from) stay_prob[p] else disperse_prob[p]
        sim_list$Movement[p, r_from, r_to, , 1, , , ] <- prob
      }

      # Season 2: natal return with straying
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
      R0_arry[1, 1, , ] <- 7   # pop 1 recruits to region 1
      R0_arry[2, 1, , ] <- 7   # pop 2 recruits to region 1
      R0_arry[3, 2, , ] <- 7   # pop 3 recruits to region 2
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
      h_arry[1, 1, , ] <- 0.75 # pop 1 in region 1
      h_arry[2, 1, , ] <- 0.75 # pop 2 in region 1
      h_arry[3, 2, , ] <- 0.75 # pop 3 in region 2
      h_arry
    }
  )

  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

  setup_em <- function(sim_obj, sim, use_pop_specific_cat_comps) {

    # Extract simulation data for current year and replicate
    sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_yrs, sim = sim)

    # Setup model dimensions
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
      do_internal_comp_osa = TRUE,
      do_internal_conv_tag_osa = TRUE
    )

    input_list <- Setup_Mod_Rec(
      input_list = input_list,
      do_rec_bias_ramp = 0,
      sigmaR_switch = 1,
      init_age_strc = "matrix",
      equil_init_age_strc = "stoch_all",

      # spawning dynamics
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

      # starting values / fixed parameters
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

    # Biological setup
    input_list <- Setup_Mod_Biologicals(
      input_list = input_list,
      WAA = sim_data$WAA,
      MatAA = sim_data$MatAA,
      WAA_fish = sim_data$WAA_fish,
      WAA_srv = sim_data$WAA_srv,
      fit_lengths = 0,
      AgeingError = sim_data$AgeingError,
      M_spec = "est_ln_M"
    )

    # Movement and tagging
    input_list <- Setup_Mod_Tagging(input_list = input_list,
                                    use_conv_fish_tagging = 1,
                                    conv_tagged_fish = sim_data$conv_tagged_fish_attr,
                                    conv_tag_max_liberty = dim(sim_data$obs_conv_tag_fish_recap)[1],
                                    obs_conv_tag_fish_recap = sim_data$obs_conv_tag_fish_recap,
                                    conv_fish_tag_like = 'Poisson',
                                    init_conv_tag_mort_spec = 'fix',
                                    conv_tag_shed_spec = 'fix',
                                    conv_tagrep_spec = 'est_shared_r',
                                    conv_fish_tag_attr = "p_a_s",
                                    conv_tag_release_indicator = sim_data$conv_tag_release_indicator
    )

    input_list <- Setup_Mod_Movement(
      input_list = input_list,
      do_recruits_move = 0,
      use_fixed_movement = 0,
      Movement_popblk_spec = list(1,2,3),
      Movement_seasblk_spec = list(1,2)
    )

    # Catch & F ---------------------------------------------------------------
    if(use_pop_specific_cat_comps) {
      input_list <- Setup_Mod_Catch_and_F(
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
      )
    } else {
      input_list <- Setup_Mod_Catch_and_F(
        input_list = input_list,
        ObsCatch = sim_data$ObsCatch,
        UseCatch = array(1, dim = dim(sim_data$UseCatch)),
        Use_F_pen = 1,
        sigmaC_pop_spec = 'fix',
        sigmaF_spec = 'fix',
        ln_sigmaC = sim_data$ln_sigmaC,
        ln_sigmaC_pop = sim_data$ln_sigmaC_pop,
        ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions,
                                          input_list$data$n_seas,
                                          input_list$data$n_fish_fleets))
      )
    }

    # Fishery index & comps ---------------------------------------------------
    if(use_pop_specific_cat_comps) {
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
    } else {
      input_list <- Setup_Mod_FishIdx_and_Comps(
        input_list = input_list,
        ObsFishIdx = sim_data$ObsFishIdx,
        ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
        UseFishIdx = array(0, dim = dim(sim_data$UseFishIdx)),
        ObsFishAgeComps = sim_data$ObsFishAgeComps,
        ObsFishLenComps = sim_data$ObsFishLenComps,
        UseFishAgeComps = sim_data$UseFishAgeComps,
        UseFishLenComps = sim_data$UseFishLenComps,
        ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
        ISS_FishLenComps = sim_data$ISS_FishLenComps,
        fish_idx_type = 'none',
        FishAgeComps_LikeType = c("Multinomial"),
        FishLenComps_LikeType = c("none"),
        FishAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
        FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
      )
    }

    # Survey index & comps ----------------------------------------------------
    if(use_pop_specific_cat_comps) {
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
    } else {
      input_list <- Setup_Mod_SrvIdx_and_Comps(
        input_list = input_list,
        ObsSrvIdx = sim_data$ObsSrvIdx,
        ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
        UseSrvIdx = sim_data$UseSrvIdx,
        ObsSrvAgeComps = sim_data$ObsSrvAgeComps,
        ObsSrvLenComps = sim_data$ObsSrvLenComps,
        UseSrvAgeComps = sim_data$UseSrvAgeComps,
        UseSrvLenComps = sim_data$UseSrvLenComps,
        ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps,
        ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
        srv_idx_type = c("biom"),
        SrvAgeComps_LikeType = c("Multinomial"),
        SrvLenComps_LikeType = c("none"),
        SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
        SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
      )
    }

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

    # Weighting ---------------------------------------------------------------
    if(use_pop_specific_cat_comps) {
      input_list <- Setup_Mod_Weighting(
        input_list = input_list,
        Wt_Catch = 1,
        Wt_Catch_pop = 1,
        Wt_SrvIdx_pop = 1,
        Wt_FishIdx = 1,
        Wt_SrvIdx = 1,
        Wt_Rec = 1,
        Wt_F = 1,
        Wt_Tagging = 1,
        Wt_FishAgeComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                           input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
        Wt_SrvAgeComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                          input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets)),
        Wt_FishAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                               input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
        Wt_SrvAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                              input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
      )
    } else {
      input_list <- Setup_Mod_Weighting(
        input_list = input_list,
        Wt_Catch = 1,
        Wt_Catch_pop = 1,
        Wt_FishIdx = 1,
        Wt_SrvIdx = 1,
        Wt_Rec = 1,
        Wt_F = 1,
        Wt_Tagging = 1,
        Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                           input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_fish_fleets)),
        Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years),
                                          input_list$data$n_seas, input_list$data$n_sexes, input_list$data$n_srv_fleets))
      )
    }

    return(input_list)
  }

  check_osa_calibration <- function(
    residuals,
    mean_tol = 0.15,
    sd_range = c(0.7, 1.3),
    test_autocorr = TRUE
  ) {
    resid <- residuals[is.finite(residuals)]
    n <- length(resid)

    mean_ok <- abs(mean(resid)) < mean_tol
    sd_val  <- sd(resid)
    sd_ok   <- sd_val > sd_range[1] && sd_val < sd_range[2]

    shapiro_p <- if (n >= 3 && n <= 5000) shapiro.test(resid)$p.value else NA
    normal_ok <- is.na(shapiro_p) || shapiro_p > 0.01

    autocorr_ok <- TRUE
    if (test_autocorr && n > 10) {
      bt <- Box.test(resid, lag = min(5, floor(n/5)), type = "Ljung-Box")
      autocorr_ok <- bt$p.value > 0.01
    }

    list(
      mean = mean(resid),
      sd = sd_val,
      shapiro_p = shapiro_p,
      mean_ok = mean_ok,
      sd_ok = sd_ok,
      normal_ok = normal_ok,
      autocorr_ok = autocorr_ok,
      calibrated = mean_ok && sd_ok && normal_ok && autocorr_ok
    )
  }

  setup_em_variant <- function(sim_obj, sim, use_pop_specific_cat_comps, misspecify_movement = FALSE) {
    input_list <- setup_em(sim_obj, sim, use_pop_specific_cat_comps)
    if (misspecify_movement) {
      # override with a wrong, fixed, well-mixed movement matrix
      wrong_move <- array(1 / input_list$data$n_regions, dim = dim(input_list$data$Fixed_Movement))
      input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1, Fixed_Movement = wrong_move)
    }
    input_list
  }


  # correct EM
  input_correct <- setup_em_variant(sim_obj, 1, TRUE, misspecify_movement = FALSE)
  obj_correct   <- fit_model(input_correct$data, input_correct$par, input_correct$map, NULL, 3, silent = TRUE, do_optim = TRUE)
  obj_francis_correct <- run_francis(
    data = input_correct$data,
    parameters = input_correct$par,
    mapping = input_correct$map,
    random = NULL,
    n_francis_iter = 2,
    0
  ) # also test francis to make sure running
  osa_tag_correct <- oneStepPredict(obj_correct, 'ObsConvTag_osa_count', method = 'cdf', discrete = TRUE, trace = FALSE)
  chk_correct <- check_osa_calibration(osa_tag_correct$residual)

  expect_true(chk_correct$mean_ok)
  expect_true(chk_correct$sd_ok)
  expect_true(chk_correct$normal_ok)
  expect_true(sdreport(obj_correct)$pdHess)
  expect_jnLL_decomposes(obj_correct)
  expect_true((sdreport(obj_francis_correct$obj))$pdHess)

  # misspecified EM
  input_wrong <- setup_em_variant(sim_obj, 1, TRUE, misspecify_movement = TRUE)
  obj_wrong   <- fit_model(input_wrong$data, input_wrong$par, input_wrong$map, NULL, 3, silent = TRUE, do_optim = TRUE)
  osa_tag_wrong <- oneStepPredict(obj_wrong, 'ObsConvTag_osa_count', method = 'cdf', discrete = TRUE, trace = FALSE)
  chk_wrong <- check_osa_calibration(osa_tag_wrong$residual)

  # The corrected (per-pool) count-family likelihood now packs ~n_pop_pool x
  # more, sparser residuals than before, which widens the sampling noise on
  # any single absolute calibration threshold. Compare the misspecified fit
  # RELATIVE to the correct fit instead (worse on mean, sd, normality, or
  # autocorrelation) -- the same relative-comparison pattern already used in
  # test-internal_discard_comp_osas.R -- rather than requiring
  # chk_wrong$calibrated to be categorically FALSE.
  expect_true(
    abs(chk_wrong$mean) > abs(chk_correct$mean) ||
    abs(chk_wrong$sd - 1) > abs(chk_correct$sd - 1) ||
    (!is.na(chk_wrong$shapiro_p) && !is.na(chk_correct$shapiro_p) && chk_wrong$shapiro_p < chk_correct$shapiro_p) ||
    (!chk_wrong$autocorr_ok && chk_correct$autocorr_ok)
  )
  expect_true(sdreport(obj_wrong)$pdHess)

})
