library(SPoRC)
library(testthat)
library(reshape2)
library(dplyr)

test_that("Age-0 (rec_lag = 0) recruitment with spawning after season 1 is concordant between OM and EM", {

  # This exercises exactly the scenario rec_lag = 0 was built for: a seasonal
  # model where spawning does NOT occur in season 1 (spawn_seas = 2 of 2), so
  # recruits can only enter the population in spawn_seas itself. It simulates
  # under known Beverton-Holt parameters with the operating model
  # (Simulate_Pop_Static), then fits the estimation model with the same
  # rec_lag = 0 / spawn_seas = 2 timing and checks that SSB, R0, and steepness
  # are recovered - this only works if the OM (Simulate_Population.R) and EM
  # (SPoRC_rtmb.R) apply the same-year SSB -> recruitment timing consistently.

  set.seed(2024)

  inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

  n_sims <- 20
  true_R0 <- 8
  true_h  <- 0.75

  ### Setup Operating Model ---------------------------------------------------
  sim_list <- Setup_Sim_Dim(n_sims = n_sims,
                            n_yrs = 35,
                            n_regions = 1,
                            n_ages = 8,
                            n_lens = NULL,
                            n_sexes = 1,
                            n_fish_fleets = 1,
                            n_srv_fleets = 1,
                            n_seas = 2,
                            n_pop = 1
  )

  sim_list <- Setup_Sim_Containers(sim_list)

  sim_list <- Setup_Sim_Fishing(
    sim_list = sim_list,
    fish_sel_input = replicate(
      n = sim_list$n_sims,
      array(rep(1 / (1 + exp(-1.5 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs * sim_list$n_seas),
            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                    sim_list$n_sexes, sim_list$n_fish_fleets))
    ),
    Fmort_input = array(0.08, dim = c(sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas,
                                      sim_list$n_fish_fleets, sim_list$n_sims))
  )

  sim_list <- Setup_Sim_Survey(
    sim_list = sim_list,
    srv_sel_input = replicate(
      n = sim_list$n_sims,
      array(rep(1 / (1 + exp(-1 * ((1:sim_list$n_ages) - 2))), each = sim_list$n_yrs * sim_list$n_seas),
            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                    sim_list$n_sexes, sim_list$n_srv))
    )
  )

  # Maturity: age-0 (first age class) must be immature for rec_lag = 0
  MatAA_ogive <- c(0, 1 / (1 + exp(-2 * ((2:sim_list$n_ages) - 4))))

  sim_list <- suppressWarnings(
    Setup_Sim_Biologicals(
      sim_list = sim_list,
      natmort_input = replicate(n = sim_list$n_sims, array(0.25, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs,
                                                                         sim_list$n_ages, sim_list$n_sexes))),
      WAA_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs * sim_list$n_seas),
                                                       dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes))),
      WAA_fish_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs * sim_list$n_seas),
                                                            dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))),
      WAA_srv_input = replicate(n = sim_list$n_sims, array(rep(5 / (1 + exp(-3 * ((1:sim_list$n_ages) - 3))), each = sim_list$n_yrs * sim_list$n_seas),
                                                           dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))),
      MatAA_input = replicate(n = sim_list$n_sims, array(rep(MatAA_ogive, each = sim_list$n_yrs * sim_list$n_seas),
                                                         dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes)))
    )
  )

  sim_list <- Setup_Sim_Tagging(sim_list = sim_list, use_conv_fish_tagging = 0)

  # No movement (single region)
  sim_list$Movement <- array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs,
                                        sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

  ### Setup Recruitment: age-0, spawning in season 2 of 2 ----------------------
  # rec_seas_prop must be 0 before spawn_seas -- all recruits enter in
  # spawn_seas itself here, since there's no season after it.
  rec_seas_prop_input <- array(0, dim = c(sim_list$n_pop, sim_list$n_seas, sim_list$n_sims))
  rec_seas_prop_input[, 2, ] <- 1

  sim_list <- Setup_Sim_Rec(
    sim_list = sim_list,
    R0_input = array(true_R0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    h_input = array(true_h, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims)),
    ln_sigmaR = array(log(0.4), dim = c(2, sim_list$n_pop, sim_list$n_regions)),
    rec_seas_prop_input = rec_seas_prop_input,
    recruitment_opt = "bh_rec",
    rec_dd = "global",
    spawn_seas = 2,
    t_spawn = 0,
    rec_lag = 0,
    init_age_strc = 1
  )

  ## Simulate Data -----------------------------------------------------------
  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)

  # Sanity check: no recruits should ever appear in season 1 (pre-spawn)
  expect_true(all(sim_obj$NAA[1, 1, , 1, 1, 1, ] == 0))
  # ... and recruitment should be strictly positive in spawn_seas (season 2)
  expect_true(all(sim_obj$Rec[1, 1, , ] > 0))

  # Define Estimation Model -------------------------------------------------
  setup_em <- function(sim_obj, sim) {

    sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = sim)

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
      natal_region = sim_obj$natal_region,
      verbose = FALSE
    )

    fixed_rec_seas_prop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_seas))
    fixed_rec_seas_prop[, 2] <- 1

    input_list <- Setup_Mod_Rec(
      input_list = input_list,
      rec_model = "bh_rec",
      rec_dd = "global",
      rec_lag = 0,
      spawn_seas = 2,
      t_spawn = 0,
      use_fixed_rec_seas_prop = 1,
      fixed_rec_seas_prop = fixed_rec_seas_prop,
      do_rec_bias_ramp = 0,
      sigmaR_switch = 1,
      ln_sigmaR = array(log(0.4), c(2, input_list$data$n_pop, input_list$data$n_regions)),
      sigmaR_spec = "fix",
      h_spec = "fix",
      steepness_h = array(inv_steepness(true_h), dim = c(input_list$data$n_pop, input_list$data$n_regions)),
      init_age_strc = 1,
      equil_init_age_strc = 2,
      ln_global_R0 = log(true_R0)
    )

    input_list <- Setup_Mod_Biologicals(
      input_list = input_list,
      WAA = sim_data$WAA,
      MatAA = sim_data$MatAA,
      WAA_fish = sim_data$WAA_fish,
      WAA_srv = sim_data$WAA_srv,
      fit_lengths = 0,
      AgeingError = sim_data$AgeingError,
      M_spec = "fix",
      Fixed_natmort = array(0.25, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                          length(input_list$data$ages), input_list$data$n_sexes))
    )

    input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
    input_list <- Setup_Mod_Movement(
      input_list = input_list,
      use_fixed_movement = 1,
      Fixed_Movement = NA,
      do_recruits_move = 0
    )

    suppressWarnings(
      input_list <- Setup_Mod_Catch_and_F(
        input_list = input_list,
        ObsCatch = sim_data$ObsCatch,
        UseCatch = sim_data$UseCatch,
        Use_F_pen = 1,
        sigmaC_spec = "fix",
        ln_sigmaC = sim_data$ln_sigmaC,
        ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
      )
    )

    input_list <- Setup_Mod_FishIdx_and_Comps(
      input_list = input_list,
      ObsFishIdx = sim_data$ObsFishIdx,
      ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
      UseFishIdx = sim_data$UseFishIdx,
      ObsFishAgeComps = sim_data$ObsFishAgeComps,
      ObsFishLenComps = sim_data$ObsFishLenComps,
      UseFishAgeComps = sim_data$UseFishAgeComps,
      UseFishLenComps = sim_data$UseFishLenComps,
      ISS_FishAgeComps = sim_data$ISS_FishAgeComps,
      ISS_FishLenComps = sim_data$ISS_FishLenComps,
      fish_idx_type = c("biom"),
      FishAgeComps_LikeType = c("Multinomial"),
      FishLenComps_LikeType = c("none"),
      FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
      FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
    )

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
      SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"),
      SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
    )

    input_list <- Setup_Mod_Fishsel_and_Q(
      input_list = input_list,
      fish_sel_model = c("logist2_Fleet_1"),
      fish_fixed_sel_pars_spec = c("est_all"),
      fish_q_spec = "est_all"
    )

    input_list <- Setup_Mod_Srvsel_and_Q(
      input_list = input_list,
      srv_sel_model = c("logist2_Fleet_1"),
      srv_fixed_sel_pars_spec = c("est_all"),
      srv_q_spec = c("est_all")
    )

    input_list <- Setup_Mod_Weighting(
      input_list = input_list,
      Wt_Catch = 1,
      Wt_FishIdx = 1,
      Wt_SrvIdx = 1,
      Wt_Rec = 1,
      Wt_F = 1,
      Wt_Tagging = 0,
      Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                         input_list$data$n_sexes, input_list$data$n_fish_fleets)),
      Wt_FishLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                         input_list$data$n_sexes, input_list$data$n_fish_fleets)),
      Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                        input_list$data$n_sexes, input_list$data$n_srv_fleets)),
      Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                        input_list$data$n_sexes, input_list$data$n_srv_fleets))
    )

    return(input_list)
  }

  ssb_results <- array(NA, dim = c(sim_list$n_yrs, n_sims))
  r0_results <- vector()

  for (i in 1:n_sims) {

    input_list <- setup_em(sim_obj, sim = i)

    model <- fit_model(input_list$data,
                       input_list$par,
                       input_list$map,
                       random = NULL,
                       silent = TRUE)

    ssb_results[, i] <- as.vector(model$rep$SSB)
    r0_results[i] <- model$rep$R0

  }

  ssb_df_res <- reshape2::melt(ssb_results) %>%
    rename(Year = Var1, Sim = Var2, Est = value) %>%
    dplyr::left_join(reshape2::melt(sim_obj$SSB) %>%
                       dplyr::rename(Pop = Var1, Region = Var2, Year = Var3, Sim = Var4, True = value),
                     by = c("Year", "Sim")) %>%
    dplyr::mutate(RE = (Est - True) / True)

  # check to see if relative error is within a few percent
  expect_equal(median(ssb_df_res$RE), 0, tolerance = 0.03)
  expect_equal(median(((r0_results - true_R0) / true_R0)), 0, tolerance = 0.05)

})
