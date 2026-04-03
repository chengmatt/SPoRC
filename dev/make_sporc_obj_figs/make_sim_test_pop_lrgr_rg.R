# Purpose: To demonstrate a natal homing model with 3 populations but 2 regions.
# Creator: Matthew LH. Cheng
# Date: 3/23/26

library(SPoRC)
library(ggplot2)
library(here)
devtools::load_all(here("R"))
set.seed(555)

### Setup Operating Model ---------------------------------------------------
sim_list <- Setup_Sim_Dim(n_sims = 1, # number of simulations
                          n_yrs = 30, # number of years
                          n_regions = 2,  # number of regions
                          n_ages = 10, # number of ages
                          n_lens = NULL, # number of lengths
                          n_sexes = 1, # number of sexes
                          n_fish_fleets = 1, # number of fishery fleets
                          n_srv_fleets = 1,  # number of survey fleets
                          n_seas = 2,  # number of seasons
                          n_pop = 3, # number of populaitons
                          natal_region = c(1, 1, 2) # natal regions
)

### Setup Simulation Containers ---------------------------------------------
sim_list <- Setup_Sim_Containers(sim_list)

### Setup Fishing Processes -------------------------------------------------
sim_list <- Setup_Sim_Fishing(
  sim_list = sim_list,

  # logistic selectivity
  fish_sel_input = replicate(
    n = sim_list$n_sims,
    {
      arr <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs,
                               sim_list$n_ages, sim_list$n_sexes, sim_list$n_fish_fleets))
      for (r in 1:sim_list$n_regions)
        for (y in 1:sim_list$n_yrs)
          for (s in 1:sim_list$n_sexes) {
            arr[r, y, , s, 1] <-  1 / (1 + exp(-1.5 * (1:sim_list$n_ages - 4)))
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

### Setup Survey Process ----------------------------------------------------
sim_list <- Setup_Sim_Survey(
  sim_list = sim_list,

  # Logistic selectivity
  srv_sel_input = replicate(
    n = sim_list$n_sims,
    {
      arr <- array(NA, dim = c(sim_list$n_regions, sim_list$n_yrs,
                               sim_list$n_ages, sim_list$n_sexes, sim_list$n_srv_fleets))
      for (r in 1:sim_list$n_regions)
        for (y in 1:sim_list$n_yrs)
          for (s in 1:sim_list$n_sexes) {
            arr[r, y, , s, 1] <-  1 / (1 + exp(-1 * (1:sim_list$n_ages - 2.5)))
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

### Setup Biological Dynamics -----------------------------------------------
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

### Setup Tagging and Movement -----------------------------------------------------------
sim_list <- Setup_Sim_Tagging(
  sim_list = sim_list,
  use_conv_fish_tagging = 1,
  n_tags = 1e3,
  conv_tag_max_liberty = 10,
  conv_fish_tag_like = "Poisson"
)

# Movement
sim_list$Movement <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
                                      sim_list$n_yrs, sim_list$n_seas, sim_list$n_ages,
                                      sim_list$n_sexes, sim_list$n_sims))

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

# Single season movement rates
sim_list$sgl_seas_spawning_movement <- array(
  NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_regions,
          sim_list$n_yrs, sim_list$n_ages, sim_list$n_sexes, sim_list$n_sims)
)

### Setup Recruitment Processes ---------------------------------------------
sim_list <- Setup_Sim_Rec(
  sim_list = sim_list,
  R0_input = {
    R0_arry <- array(0, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
    R0_arry[1, 1, , ] <- 7   # pop 1 recruits to region 1
    R0_arry[2, 1, , ] <- 7   # pop 2 recruits to region 1
    R0_arry[3, 2, , ] <- 7  # pop 3 recruits to region 2
    R0_arry
  },
  ln_sigmaR = array(log(0.5), dim = c(2, sim_list$n_pop, sim_list$n_regions)), # sigma R
  init_age_strc = "matrix", # matrix geometic series to initialize population
  recruitment_opt = 'bh_rec', # BH recruitment
  rec_dd = 'local', # localized DD
  spawn_seas = 2, # spawning season
  t_spawn = 0.5, # spawn timing
  h_input = {
    h_arry <- array(NA, dim = c(sim_list$n_pop, sim_list$n_regions, sim_list$n_yrs, sim_list$n_sims))
    h_arry[1, 1, , ] <- 0.75 # pop 1 in region 1
    h_arry[2, 1, , ] <- 0.75 # pop 2 in region 1
    h_arry[3, 2, , ] <- 0.75  # pop 3 in region 2
    h_arry
  }
)


## Simulate Data -----------------------------------------------------------
sim_obj <- Simulate_Pop_Static(sim_list = sim_list, output_path = NULL) # get simulated datasets


# Define Estimation Model -------------------------------------------------
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
    verbose = FALSE
  )

  input_list <- Setup_Mod_Rec(
    input_list = input_list,
    do_rec_bias_ramp = 0, # not doing bias ramp
    sigmaR_switch = 1, # when to switch from early to late sigmaR (switch in first year)
    init_age_strc = "matrix", # scalar geometric series to derive initial age structure
    equil_init_age_strc = "stoch_all", # estimating all intial age deviations

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
    # Data inputs
    WAA = sim_data$WAA,
    MatAA = sim_data$MatAA,
    WAA_fish = sim_data$WAA_fish,
    WAA_srv = sim_data$WAA_srv,
    fit_lengths = 0, # not fitting lengths
    AgeingError = sim_data$AgeingError,
    M_spec = "est_ln_M" # esimate constant M
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
    Movement_popblk_spec = list(1,2,3), # movement population blocks
    Movement_seasblk_spec = list(1,2) # movement seas blocks
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
      pop_FishAgeComps_LikeType = c("Multinomial"),
      pop_FishAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
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
      pop_SrvAgeComps_LikeType = c("Multinomial"),
      pop_SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1"),
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



# Run EM w/ single sim ------------------------------------------------------------------

# setup EM (w/o pop-specific data)
input_list <- setup_em(sim_obj, sim = 1, use_pop_specific_cat_comps = FALSE)

# extract data, pars, and mapping
data <- input_list$data
pars <- input_list$par
map <- input_list$map

# make AD model function
non_pop_obj <- fit_model(data, pars, map, NULL, 3, silent =F, do_optim = T)
non_pop_obj$sd_rep <- sdreport(non_pop_obj)

# setup EM (w/ pop-specific data)
input_list <- setup_em(sim_obj, sim = 1, use_pop_specific_cat_comps = TRUE)

# extract data, pars, and mapping
data <- input_list$data
pars <- input_list$par
map <- input_list$map

# make AD model function
pop_obj <- fit_model(data, pars, map, NULL, 3, silent =F, do_optim = T)
pop_obj$sd_rep <- sdreport(pop_obj)

# Comparison w/ single sim --------------------------------------------------------------

# SSB Comparison
ssb_comp <- rbind(reshape2::melt(non_pop_obj$rep$SSB) %>% mutate(type = 'no_pop_data', se = non_pop_obj$sd_rep$sd[names( non_pop_obj$sd_rep$value) == 'log_SSB'] ),
                  reshape2::melt(pop_obj$rep$SSB) %>% mutate(type = 'pop_data', se = pop_obj$sd_rep$sd[names( pop_obj$sd_rep$value) == 'log_SSB']),
                  reshape2::melt(sim_obj$SSB[,,,1]) %>% mutate(type = 'truth',  se = NA)) %>%
  rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  mutate(Pop = paste("Pop", Pop), Region = paste("Region", Region))

png(here("vignettes", "figures", "r_natal_home_ssb_plot.png"), width = 1000, height = 1000)
ggplot(ssb_comp, aes(x = Year, y = value, color = type, fill = type)) +
  geom_line(aes(group = interaction(type, Pop, Region))) +
  geom_ribbon(
    data = ssb_comp,
    aes(
      ymin = exp(log(value) - 1.96 * se),
      ymax = exp(log(value) + 1.96 * se),
      group = interaction(type, Pop, Region)
    ),
    alpha = 0.2,
    color = NA
  ) +
  ggh4x::facet_grid2(Pop ~ Region, scales = 'free', independent = 'all') +
  theme_sablefish() +
  labs(y = 'SSB')
dev.off()

# Recruitment Comparison
rec_comp <- rbind(reshape2::melt(non_pop_obj$rep$Rec) %>% mutate(type = 'no_pop_data', se = non_pop_obj$sd_rep$sd[names( non_pop_obj$sd_rep$value) == 'log_Rec'] ),
                  reshape2::melt(pop_obj$rep$Rec) %>% mutate(type = 'pop_data', se = pop_obj$sd_rep$sd[names( pop_obj$sd_rep$value) == 'log_Rec']),
                  reshape2::melt(sim_obj$Rec[,,,1]) %>% mutate(type = 'truth',  se = NA)) %>%
  rename(Pop = Var1, Region = Var2, Year = Var3) %>%
  mutate(Pop = paste("Pop", Pop), Region = paste("Region", Region))

png(here("vignettes", "figures", "r_natal_home_rec_plot.png"), width = 1000, height = 1000)
ggplot(rec_comp, aes(x = Year, y = value, color = type, fill = type)) +
  geom_line(aes(group = interaction(type, Pop, Region))) +
  geom_ribbon(
    data = rec_comp,
    aes(
      ymin = exp(log(value) - 1.96 * se),
      ymax = exp(log(value) + 1.96 * se),
      group = interaction(type, Pop, Region)
    ),
    alpha = 0.2,
    color = NA
  ) +
  ggh4x::facet_grid2(Pop ~ Region, scales = 'free', independent = 'all') +
  theme_sablefish() +
  labs(y = 'Recruitment')
dev.off()



# Run EMs w/ multiple sims (parrallel) ------------------------------------

# Helper function to run sims
run_sims <- function(sims = 50,
                     sim_list,
                     n_cores = parallel::detectCores() - 2,
                     use_pop_specific_cat_comps = FALSE) {

  sim_list_1 <- sim_list
  sim_list_1$n_sims <- 1

  # set up cores and parralellizaiton framework
  future::plan(future::multisession, workers = n_cores)
  options(future.globals.maxSize = 5e3 * 1024^2)
  on.exit(future::plan(future::sequential), add = TRUE)

  results <- future.apply::future_lapply(seq_len(sims), function(i) {

    devtools::load_all(here::here("R")) # load in via dev tools

    tryCatch({

      # simulate single sim
      sim_obj    <- Simulate_Pop_Static(sim_list = sim_list_1, output_path = NULL)

      # get EM data
      input_list <- setup_em(sim_obj, sim = 1,
                             use_pop_specific_cat_comps = use_pop_specific_cat_comps)

      # fit model
      obj <- fit_model(input_list$data, input_list$par, input_list$map,
                       NULL, 3, silent = TRUE, do_optim = TRUE)

      # get sdreport
      obj$sd_rep <- RTMB::sdreport(obj)

      # output object
      list(
        em_Rec    = obj$rep$Rec,
        em_SSB    = obj$rep$SSB,
        om_Rec    = sim_obj$Rec[,,, 1],
        om_SSB    = sim_obj$SSB[,,, 1],
        max_grad  = max(abs(obj$gr(obj$opt$par)))
      )

    }, error = function(e) {

      # return NA if faile
      message(sprintf("sim %d failed: %s", i, conditionMessage(e)))
      list(em_Rec = NA, em_SSB = NA, om_Rec = NA, om_SSB = NA,
           converged = FALSE, max_grad = NA)
    })

  }, future.seed = TRUE)

  # Combine results
  collate <- function(field) {
    vals <- lapply(results, `[[`, field)
    if (any(sapply(vals, \(x) identical(x, NA)))) vals
    else simplify2array(vals)  # last dim = n_sims
  }

  list(
    em_Rec    = collate("em_Rec"),
    em_SSB    = collate("em_SSB"),
    om_Rec    = collate("om_Rec"),
    om_SSB    = collate("om_SSB"),
    max_grad  = sapply(results, `[[`, "max_grad")
  )

}

# Helper to tidy one array with labels
tidy_bias <- function(em_arr, om_arr, quantity, model_type) {
  melt(em_arr) %>%
    left_join(
      melt(om_arr) %>% rename(true = value),
      by = c("Var1", "Var2", "Var3", "Var4")
    ) %>%
    mutate(
      rel_bias  = (value - true) / true,
      quantity   = quantity,
      model_type = model_type
    )
}

# Run Sims ----------------------------------------------------------------
non_pop_obj <- run_sims(50,sim_list, n_cores = 16, use_pop_specific_cat_comps = FALSE) # no population-specific data
pop_obj <- run_sims(50,sim_list, n_cores = 16, use_pop_specific_cat_comps = TRUE) # population-specific data

# Plot Sim Results --------------------------------------------------------
df_all <- bind_rows(
  tidy_bias(non_pop_obj$em_Rec, non_pop_obj$om_Rec, "Recruitment", "Non-Population"),
  tidy_bias(pop_obj$em_Rec,     pop_obj$om_Rec,     "Recruitment", "Population"),
  tidy_bias(non_pop_obj$em_SSB, non_pop_obj$om_SSB, "SSB",         "Non-Population"),
  tidy_bias(pop_obj$em_SSB,     pop_obj$om_SSB,     "SSB",         "Population")
) %>%
  rename(pop = Var1, region = Var2, year = Var3, sim = Var4) %>%
  mutate(pop = paste("Population", pop),
         region = paste("Region", region)
         )

# summarize
df_summary <- df_all %>%
  group_by(pop, region, year, quantity, model_type) %>%
  summarise(
    med   = median(rel_bias, na.rm = TRUE),
    lo    = quantile(rel_bias, 0.025, na.rm = TRUE),
    hi    = quantile(rel_bias, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# SSB plot
png(here("vignettes", "figures", "r_natal_home_ssb_plot_multisim.png"), width = 1000, height = 1000)
ggplot(df_summary %>% filter(quantity == 'SSB'), aes(x = year, colour = model_type, fill = model_type)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(aes(y = med), linewidth = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, colour = "black") +
  scale_colour_manual(values = c("#E07B39", "#3A86C8")) +
  scale_fill_manual(  values = c("#E07B39", "#3A86C8")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(pop ~ region, scales = "free_y") +
  labs(x = "Year", y = "Relative error (SSB)", color = 'Data Scenario', fill = 'Data Scenario')+
  theme_bw(base_size = 15) +
  theme(legend.position = 'top')
dev.off()

# Rec plot
png(here("vignettes", "figures", "r_natal_home_rec_plot_multisim.png"), width = 1000, height = 1000)
ggplot(df_summary %>% filter(quantity == 'Recruitment'), aes(x = year, colour = model_type, fill = model_type)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(aes(y = med), linewidth = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, colour = "black") +
  scale_colour_manual(values = c("#E07B39", "#3A86C8")) +
  scale_fill_manual(  values = c("#E07B39", "#3A86C8")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(pop ~ region, scales = "free_y") +
  labs(x = "Year", y = "Relative error (Recruitment)", color = 'Data Scenario', fill = 'Data Scenario')+
  theme_bw(base_size = 15) +
  theme(legend.position = 'top')
dev.off()




