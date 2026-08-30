library(SPoRC)
library(testthat)

test_that("Do_Population_Projection: age-0 (rec_lag = 0) and classic (rec_lag = 1) recruitment converge to the same fished equilibrium", {

  # This checks concordance directly rather than against Get_Reference_Points:
  # while investigating this test, comparing against Get_Reference_Points'
  # BH_MSY reference points exposed a separate, PRE-EXISTING bug (confirmed
  # reproducible with the classic rec_lag = 1 path, completely untouched by
  # the age-0 recruitment work) where its per-recruit calculation doesn't
  # match actual population dynamics when rec_seas_prop concentrates
  # recruitment in a season other than the first (spawn_task filed
  # separately - not something to paper over here).
  #
  # So instead: since Get_Det_Recruitment's SBPR/S0 math is provably
  # rec_lag-invariant (rec_lag only selects which year's SSB feeds the BH
  # curve, never the per-recruit ratios themselves), projecting the SAME
  # fitted population under the SAME constant F with rec_lag = 0 vs
  # rec_lag = 1 should converge to the IDENTICAL equilibrium SSB. That's a
  # direct, decisive test of whether Do_Population_Projection's age-0 code
  # path is concordant with the pre-existing classic path, without depending
  # on the separately-broken reference point calculation.

  set.seed(4041)

  inv_steepness <- function(s) qlogis((s - 0.2) / 0.8)

  true_R0 <- 8
  true_h  <- 0.75

  ### Simulate one operating-model replicate, then fit an age-0 (rec_lag = 0) EM
  sim_list <- Setup_Sim_Dim(n_sims = 1, n_yrs = 35, n_regions = 1, n_ages = 8,
                            n_lens = NULL, n_sexes = 1, n_fish_fleets = 1,
                            n_srv_fleets = 1, n_seas = 2, n_pop = 1)

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

  sim_list$Movement <- array(1, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions, sim_list$n_yrs,
                                        sim_list$n_seas, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims))

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

  sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL)
  sim_data <- simulation_data_to_SPoRC(sim_env = sim_obj, y = sim_obj$n_years, sim = 1)

  input_list <- Setup_Mod_Dim(
    years = 1:sim_obj$n_years, ages = 1:sim_obj$n_ages, lens = sim_obj$n_lens,
    n_regions = sim_obj$n_regions, n_sexes = sim_obj$n_sexes, n_fish_fleets = sim_obj$n_fish_fleets,
    n_srv_fleets = sim_obj$n_srv_fleets, n_seas = sim_obj$n_seas, n_pop = sim_obj$n_pop,
    natal_region = sim_obj$natal_region, verbose = FALSE
  )

  fixed_rec_seas_prop <- array(0, dim = c(input_list$data$n_pop, input_list$data$n_seas))
  fixed_rec_seas_prop[, 2] <- 1

  input_list <- Setup_Mod_Rec(
    input_list = input_list, rec_model = "bh_rec", rec_dd = "global", rec_lag = 0,
    spawn_seas = 2, t_spawn = 0, use_fixed_rec_seas_prop = 1, fixed_rec_seas_prop = fixed_rec_seas_prop,
    do_rec_bias_ramp = 0, sigmaR_switch = 1,
    ln_sigmaR = array(log(0.4), c(2, input_list$data$n_pop, input_list$data$n_regions)),
    sigmaR_spec = "fix", h_spec = "fix",
    steepness_h = array(inv_steepness(true_h), dim = c(input_list$data$n_pop, input_list$data$n_regions)),
    init_age_strc = 1, equil_init_age_strc = 2, ln_global_R0 = log(true_R0)
  )

  input_list <- Setup_Mod_Biologicals(
    input_list = input_list, WAA = sim_data$WAA, MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish, WAA_srv = sim_data$WAA_srv, fit_lengths = 0,
    AgeingError = sim_data$AgeingError, M_spec = "fix",
    Fixed_natmort = array(0.25, dim = c(input_list$data$n_pop, input_list$data$n_regions, length(input_list$data$years),
                                        length(input_list$data$ages), input_list$data$n_sexes))
  )

  input_list <- Setup_Mod_Tagging(input_list = input_list, use_conv_fish_tagging = 0)
  input_list <- Setup_Mod_Movement(input_list = input_list, use_fixed_movement = 1, Fixed_Movement = NA, do_recruits_move = 0)

  suppressWarnings(
    input_list <- Setup_Mod_Catch_and_F(
      input_list = input_list, ObsCatch = sim_data$ObsCatch, UseCatch = sim_data$UseCatch,
      Use_F_pen = 1, sigmaC_spec = "fix", ln_sigmaC = sim_data$ln_sigmaC,
      ln_sigmaF = array(log(1), dim = c(input_list$data$n_regions, input_list$data$n_seas, input_list$data$n_fish_fleets))
    )
  )

  input_list <- Setup_Mod_FishIdx_and_Comps(
    input_list = input_list, ObsFishIdx = sim_data$ObsFishIdx, ObsFishIdx_SE = sim_data$ObsFishIdx_SE,
    UseFishIdx = sim_data$UseFishIdx, ObsFishAgeComps = sim_data$ObsFishAgeComps, ObsFishLenComps = sim_data$ObsFishLenComps,
    UseFishAgeComps = sim_data$UseFishAgeComps, UseFishLenComps = sim_data$UseFishLenComps,
    ISS_FishAgeComps = sim_data$ISS_FishAgeComps, ISS_FishLenComps = sim_data$ISS_FishLenComps,
    fish_idx_type = c("biom"), FishAgeComps_LikeType = c("Multinomial"), FishLenComps_LikeType = c("none"),
    FishAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"), FishLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_SrvIdx_and_Comps(
    input_list = input_list, ObsSrvIdx = sim_data$ObsSrvIdx, ObsSrvIdx_SE = sim_data$ObsSrvIdx_SE,
    UseSrvIdx = sim_data$UseSrvIdx, ObsSrvAgeComps = sim_data$ObsSrvAgeComps, ObsSrvLenComps = sim_data$ObsSrvLenComps,
    UseSrvAgeComps = sim_data$UseSrvAgeComps, UseSrvLenComps = sim_data$UseSrvLenComps,
    ISS_SrvAgeComps = sim_data$ISS_SrvAgeComps, ISS_SrvLenComps = sim_data$ISS_SrvLenComps,
    srv_idx_type = c("biom"), SrvAgeComps_LikeType = c("Multinomial"), SrvLenComps_LikeType = c("none"),
    SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1"), SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1")
  )

  input_list <- Setup_Mod_Fishsel_and_Q(
    input_list = input_list, fish_sel_model = c("logist2_Fleet_1"),
    fish_fixed_sel_pars_spec = c("est_all"), fish_q_spec = "est_all"
  )

  input_list <- Setup_Mod_Srvsel_and_Q(
    input_list = input_list, srv_sel_model = c("logist2_Fleet_1"),
    srv_fixed_sel_pars_spec = c("est_all"), srv_q_spec = c("est_all")
  )

  input_list <- Setup_Mod_Weighting(
    input_list = input_list, Wt_Catch = 1, Wt_FishIdx = 1, Wt_SrvIdx = 1, Wt_Rec = 1, Wt_F = 1, Wt_Tagging = 0,
    Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_FishLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                       input_list$data$n_sexes, input_list$data$n_fish_fleets)),
    Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets)),
    Wt_SrvLenComps = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
                                      input_list$data$n_sexes, input_list$data$n_srv_fleets))
  )

  model <- fit_model(input_list$data, input_list$par, input_list$map, random = NULL, silent = TRUE)
  expect_jnLL_decomposes(model)

  data <- input_list$data
  rep <- model$rep

  ### Set up shared projection inputs ------------------------------------------
  n_proj_yrs <- 500
  n_regions <- data$n_regions; n_ages <- length(data$ages); n_sexes <- data$n_sexes
  n_fish_fleets <- data$n_fish_fleets; n_seas <- data$n_seas; n_pop <- data$n_pop
  t_spawn <- data$t_spawn
  do_recruits_move <- 0
  n_yrs_hist <- length(data$years)
  const_F <- 0.15 # arbitrary constant fishing mortality, applied identically in both projections

  terminal_NAA <- array(rep$NAA[,,n_yrs_hist,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))
  terminal_NAA0 <- array(rep$NAA0[,,n_yrs_hist,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes))

  WAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  for(y in 1:n_proj_yrs) WAA[,,y,,,] <- data$WAA[,,n_yrs_hist,,,]
  WAA_fish <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  for(y in 1:n_proj_yrs) WAA_fish[,,y,,,,] <- data$WAA[,,n_yrs_hist,,,]

  MatAA <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))
  for(y in 1:n_proj_yrs) MatAA[,,y,,,] <- data$MatAA[,,n_yrs_hist,,,]

  fish_sel <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  for(y in 1:n_proj_yrs) fish_sel[,,y,,,,] <- rep$fish_sel[,,n_yrs_hist,,,,]

  ret_sel <- array(0, dim = c(n_pop, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes, n_fish_fleets))
  for(y in 1:n_proj_yrs) ret_sel[,,y,,,,] <- rep$ret_sel[,,n_yrs_hist,,,,]

  Movement <- array(1, dim = c(n_pop, n_regions, n_regions, n_proj_yrs, n_seas, n_ages, n_sexes))

  terminal_F <- array(const_F / n_seas, dim = c(n_regions, n_seas, n_fish_fleets)) # just defines the seasonal F ratio
  terminal_dmr <- array(rep$dmr[,n_yrs_hist,,], dim = c(n_regions, n_seas, n_fish_fleets))

  natmort_slice <- rep$natmort[,, n_yrs_hist, , ]
  natmort <- array(rep(natmort_slice, each = n_proj_yrs), dim = c(n_pop, n_regions, n_proj_yrs, n_ages, n_sexes))

  sexratio <- array(1, dim = c(n_pop, n_regions, n_proj_yrs, n_sexes)) # single-sex model (n_sexes = 1)

  srr_opt_base <- list(
    rec_dd = 1,
    do_recruits_move = do_recruits_move,
    R0 = rep$R0,
    h = array(rep$h_trans, dim = c(n_pop, n_regions)),
    rec_region_prop = rep$rec_region_prop,
    WAA = array(data$WAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    MatAA = array(data$MatAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
    SSB = rep$SSB,
    Movement = array(Movement[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
    sex_ratio_f = array(1, dim = c(n_pop, n_regions)), # single-sex model (n_sexes = 1)
    sgl_seas_spawning_movement = NULL,
    stray_rate = array(0, dim = c(n_pop)),
    natmort = array(natmort_slice, dim = c(n_pop, n_regions, n_ages)),
    fish_sel = array(rep$fish_sel[,,n_yrs_hist,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    ret_sel = array(rep$ret_sel[,,n_yrs_hist,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)),
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)),
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  )

  run_proj <- function(rec_lag) {
    srr_opt <- srr_opt_base
    srr_opt$rec_lag <- rec_lag
    Do_Population_Projection(
      n_proj_yrs = n_proj_yrs, n_pop = n_pop, n_regions = n_regions, n_ages = n_ages, n_sexes = n_sexes,
      sexratio = sexratio, n_fish_fleets = n_fish_fleets, do_recruits_move = do_recruits_move,
      rec_seas_prop = array(rep$rec_seas_prop, dim = c(n_pop, n_seas)),
      recruitment = array(rep$Rec, dim = c(n_pop, n_regions, n_yrs_hist)),
      terminal_NAA = terminal_NAA, terminal_NAA0 = terminal_NAA0, terminal_F = terminal_F, dmr = terminal_dmr,
      natmort = natmort, WAA = WAA, WAA_fish = WAA_fish, MatAA = MatAA, fish_sel = fish_sel, ret_sel = ret_sel,
      Movement = Movement,
      f_ref_pt = array(const_F, dim = c(n_regions, n_proj_yrs)), # hold F constant - no HCR
      recruitment_opt = "bh_rec", fmort_opt = "Input",
      t_spawn = t_spawn, n_seas = n_seas, seasdur = data$seasdur, spawn_seas = data$spawn_seas,
      natal_region = data$natal_region, srr_opt = srr_opt
    )
  }

  out_age0 <- run_proj(rec_lag = 0)  # the new path
  out_lag1 <- run_proj(rec_lag = 1)  # the pre-existing, unmodified path

  # Under identical constant F, both should settle at their own stable
  # equilibrium...
  expect_equal(as.numeric(out_age0$proj_SSB[1,1,n_proj_yrs]), as.numeric(out_age0$proj_SSB[1,1,n_proj_yrs - 1]), tolerance = 1e-6)
  expect_equal(as.numeric(out_lag1$proj_SSB[1,1,n_proj_yrs]), as.numeric(out_lag1$proj_SSB[1,1,n_proj_yrs - 1]), tolerance = 1e-6)

  # ...and, since Get_Det_Recruitment's per-recruit math never depends on
  # rec_lag, that equilibrium should be THE SAME regardless of whether
  # recruitment lags a year or is age-0.
  expect_equal(as.numeric(out_age0$proj_SSB[1,1,n_proj_yrs]), as.numeric(out_lag1$proj_SSB[1,1,n_proj_yrs]), tolerance = 1e-6)

  # No recruits ever appear pre-spawn (season 1, age index 1) in the age-0 path
  expect_true(all(out_age0$proj_NAA[1,1,,1,1,1] == 0))

})
