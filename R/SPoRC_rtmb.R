# version 1 - 12/22/24 (M.LH Cheng)
# Bridge model 23.5 from ADMB to RTMB
# Changed code to be more modular, accommodating any number of fishery and survey fleets
# Rectified errors in fitting to length composition data (normalize proportions at length after conversion from age-length matrix)
# Changed survey composition data to be calculated using survival midyear
# Added options for continuous time-varying selectivity
# Added options for TMB / R-like likelihoods (e.g., dnorm) instead of custom likelihoods
# Added options for dirichlet multinomial likelihood

# version 2 - 12/23/24 (M.LH Cheng)
# Incorporated options to fit age and length composition data as sex-aggregated, split by sex (no sex ratio
# information), and jointly by sex (implicit sex ratio information)
# Added in option for Dirichlet Multinomial likelihood

# version 3 - (M.LH Cheng)
# Coded in spatial dimensions
# Parameters (mean recruitment, recruitment devs, initial age devs,
# selectivity, composition likelihood parameters, catchability,
# mean fishing mortality, fishing mortality deviates) can be estimated spatially
# Incorporated options to allow for estimation of movement parameters across years, ages, and sexes
# Tag integrated model incorporated using a Brownie Tag Attrition Model
# Tag Reporting Rates, Tag Shedding, and Tag Induced Mortality are parameters that can be estimated
# Beta priors for tag reporting rates, dirichlet priors for movement rates
# Incorporated iid, random walk, 2d and 3d correaltions for fishery and survey selectivity
# Added in options for Logistic Normal likelihood

# version 4 - (M.LH Cheng)
# Added in capabilities for length-based selectivity processes
# Removed unncessary constants added to likelihoods
# Refactored natural mortality module
# Removed ADMB likelihoods and Sablefish specific calculations
# Added equilibrium plus group calculations to initial abundance when movement occurs
# Refactored Initial Numbers at Age module

# version 5 - (M.LH Cheng)
# Refactored movement priors and movement setup
# Added in CTMC movement


#' Generalized RTMB model
#'
#' @param pars Parameter List
#' @param data Data List
#' @import RTMB
#' @keywords internal
SPoRC_rtmb = function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data, warn = FALSE) # load in starting values and data

  # Model Set Up (Containers) -----------------------------------------------
  n_ages = length(ages) # number of ages
  n_yrs = length(years) # number of years
  n_lens = length(lens) # number of lengths

  # Recruitment stuff
  n_est_rec_devs = dim(ln_RecDevs)[2] # number of recruitment deviates estimated
  Rec = array(0, dim = c(n_regions, n_yrs)) # Recruitment
  R0 = rep(0, n_regions) # R0 or mean recruitment
  sexratio = array(0, dim = c(n_regions, n_yrs, n_sexes)) # recruitment sex ratio

  # Population Dynamics
  NAA = array(data = 0, dim = c(n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Numbers at age
  NAA_bef = array(data = 0, dim = c(n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Numbers at age before movement
  NAA_aft = array(data = 0, dim = c(n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Numbers at age after movement
  NAA0 = array(data = 0, dim = c(n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Unfished Numbers at age
  ZAA = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes)) # Total mortality at age
  natmort = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # natural mortality at age
  Total_Biom = array(0, dim = c(n_regions, n_yrs)) # Total biomass
  SSB = array(0, dim = c(n_regions, n_yrs)) # Spawning stock biomass
  Dynamic_SSB0 = array(0, dim = c(n_regions, n_yrs)) # Dynamic Unfished Spawning stock biomass
  Aggregated_SSB = array(0, dim = c(n_yrs)) # Aggregated Spawning stock biomass
  Dynamic_Aggregated_SSB0 = array(0, dim = c(n_yrs)) # Dynamic Unfished Aggregated Spawning stock biomass

  # Movement Stuff
  Movement = array(data = 0, dim = c(n_regions, n_regions, n_yrs + n_proj_yrs_devs, n_seas, n_ages, n_sexes)) # movement "matrix"
  n_move_age_tag_pool = length(move_age_tag_pool) # number of ages to pool for tagging data
  n_move_sex_tag_pool = length(move_sex_tag_pool) # number of sexes to pool for tagging data

  # Tagging Stuff
  Tags_Avail = array(data = 0, dim = c(max_tag_liberty + 1, n_seas, n_tag_cohorts, n_regions, n_ages, n_sexes)) # Tags availiable for recapture
  Tag_Reporting = array(data = 0, dim = c(n_regions, n_yrs)) # Tag reporting rate
  Pred_Tag_Recap = array(data = 0, dim = c(max_tag_liberty, n_seas, n_tag_cohorts, n_regions, n_ages, n_sexes)) # predicted recaptures

  # Fishery Processes
  init_F = init_F_prop * exp(ln_F_mean[1]) # initial F for age structure
  Fmort = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Fishing mortality scalar
  FAA = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Fishing mortality at age
  CAA = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Catch at age
  CAL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets)) # Catch at length
  PredCatch = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted catch in weight
  PredFishIdx = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted fishery index
  fish_sel = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_ages, n_sexes, n_fish_fleets)) # Fishery selectivity
  fish_sel_l = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_lens, n_sexes, n_fish_fleets)) # Fishery selectivity (lengths)
  fish_q = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery catchability

  # Survey Processes
  SrvIAA = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)) # Survey index at age
  SrvIAL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets)) # Survey index at length
  PredSrvIdx = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets)) # Predicted survey index
  srv_sel = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_ages, n_sexes, n_srv_fleets)) # Survey selectivity ages
  srv_sel_l = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_lens, n_sexes, n_srv_fleets)) # Survey selectivity lengths
  srv_q = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Survey catchability

  # Likelihoods
  Catch_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Fishery Catch Likelihoods
  FishIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Fishery Index Likelihoods
  FishAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Fishery Age Comps Likelihoods
  FishLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Fishery Length Comps Likelihoods
  SrvIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets)) # Survey Index Likelihoods
  SrvAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)) # Survey Age Comps Likelihoods
  SrvLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods
  Tag_nLL = array(data = 0, dim = c(max_tag_liberty, n_seas, n_tag_cohorts, n_regions, n_ages, n_sexes)) # Tagging Likelihoods

  # Penalties and Priors
  Fmort_nLL = array(0, dim = dim(ln_F_devs)) # Fishing Mortality Deviation penalty
  Rec_nLL = array(0, dim = dim(ln_RecDevs)) # Recruitment penalty
  Init_Rec_nLL = array(0, dim = dim(ln_InitDevs)) # Initial Recruitment penalty
  sel_nLL = 0 # Penalty for selectivity deviations
  fish_q_nLL = 0 # Prior/penalty for fishery q
  srv_q_nLL = 0 # Prior/penalty for survey q
  M_nLL = 0 # Penalty/Prior for natural mortality
  h_nLL = 0 # Prior for steepness
  Movement_nLL = 0 # Penalty for movement rates
  TagRep_nLL = 0 # penalty for tag reporting rate
  Rec_prop_nLL = # penalty / prior for recruitment proportions
  jnLL = 0 # Joint negative log likelihood

  # Model Process Equations -------------------------------------------------
  ## Movement Parameters (Set up) --------------------------------------------
  out_move = Get_Movement(
    move_type = move_type, # movement type (unstructured markov, or continuous time markov chain)
    do_recruits_move = do_recruits_move,

    # Dimensions
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_proj_yrs_devs = n_proj_yrs_devs,
    n_ages = n_ages,
    n_sexes = n_sexes,
    n_seas = n_seas,

    # If move_type == 0
    move_pars = move_pars, # movement parameters for unstructred markov
    move_devs = move_devs, # logit movement deviations
    use_fixed_movement = use_fixed_movement, # indicator for fixed movement
    Fixed_Movement = Fixed_Movement, # fixed movement matrix

    # If move_type == 1
    ctmc_move_dat = ctmc_move_dat,
    preference_formula = preference_formula,
    diffusion_formula = diffusion_formula,
    log_move_diffusion_pars = log_move_diffusion_pars,
    move_preference_pars = move_preference_pars,
    area_r = area_r,
    adjacency_mat = adjacency_mat,
    ctmc_diffusion_bounds = ctmc_diffusion_bounds
  )

  # output movement stuff into model
  Mrate = out_move$Mrate
  Movement = out_move$Movement
  Movement_nLL = Movement_nLL + out_move$move_pen

  ## Natural Mortality Parameters (Set up) -----------------------------------
  if(use_fixed_natmort == 0) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            idx = M_blocks[r,y,a,s] # Get indexing from blocking structure on base natural mortality
            natmort[r,y,a,s] = exp(ln_M[idx]) # input into natural mortality array
          } # end s loop
        } # end a loop
      } # end y loop
    } # end r loop
  } # if using fixed natural mortality

  if(use_fixed_natmort == 1) natmort = Fixed_natmort # Using fixed natural mortality

  ## Fishery Selectivity -----------------------------------------------------
  if(Selex_Type == 0) selex_bins = ages # if age-based selectivity
  if(Selex_Type == 1) selex_bins = lens # if length-based selectivity

  for(r in 1:n_regions) {
    for(y in 1:(n_yrs + n_proj_yrs_devs)) {
      for(f in 1:n_fish_fleets) {
        for(s in 1:n_sexes) {

          # Extract variables
          if(y <= n_yrs) { # non-projection years
            fish_sel_blk_idx = fish_sel_blocks[r,y,f] # selectivity block indices
            tmp_fish_sel_model = fish_sel_model[r,y,f] # fishery selectivity model
            if(Selex_Type == 1) tmp_sizeage = SizeAgeTrans[r,y,seas,,,s] # size age transition matrix to use
          } else {
            fish_sel_blk_idx = fish_sel_blocks[r,n_yrs,f] # selectivity block indices
            tmp_fish_sel_model = fish_sel_model[r,n_yrs,f] # fishery selectivity model
            if(Selex_Type == 1) tmp_sizeage = SizeAgeTrans[r,n_yrs,seas,,,s] # size age transition matrix to use
          }

          # Extract out fixed-effect selectivity parameters for a given block
          tmp_fish_sel_vec = ln_fish_fixed_sel_pars[r,,fish_sel_blk_idx,s,f]

          # Compute selectivity functional form
          tmp_sel = Get_Selex(Selex_Model = tmp_fish_sel_model, # selectivity model
                              TimeVary_Model = cont_tv_fish_sel[r,f], # time varying model
                              ln_Pars = tmp_fish_sel_vec, # fixed effect selectivity parameters
                              ln_seldevs = ln_fishsel_devs[,,,,f, drop = FALSE], # Selectivity deviations
                              Region = r, # region index
                              Year = y, # year index
                              Bin = selex_bins, # bin vector
                              Sex = s # sex index
                              )

          # Calculate selectivity
          if(Selex_Type == 0) fish_sel[r,y,,s,f] = tmp_sel # age-based selectivity
          if(Selex_Type == 1) {
            fish_sel_l[r,y,,s,f] = tmp_sel # input into length-based fishery selectivity
            fish_sel[r,y,,s,f] = tmp_sel %*% tmp_sizeage # length-based selectivity (dot product of size age transition)
          }

        } # end s loop
      } # end f loop
    } # end y loop
  } # end r loop


  ## Survey Selectivity ------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:(n_yrs + n_proj_yrs_devs)) {
      for(sf in 1:n_srv_fleets) {
        for(s in 1:n_sexes) {

          # Extract variables
          if(y <= n_yrs) { # non-projection years
            srv_sel_blk_idx = srv_sel_blocks[r,y,sf] # selectivity block indices
            tmp_srv_sel_model = srv_sel_model[r,y,sf] # survey selectivity model
            if(Selex_Type == 1) tmp_sizeage = SizeAgeTrans[r,y,seas,,,s] # size age transition matrix to use
          } else {
            srv_sel_blk_idx = srv_sel_blocks[r,n_yrs,sf] # selectivity block indices
            tmp_srv_sel_model = srv_sel_model[r,n_yrs,sf] # survey selectivity model
            if(Selex_Type == 1) tmp_sizeage = SizeAgeTrans[r,n_yrs,seas,,,s] # size age transition matrix to use
          }

          # Extract out fixed-effect selectivity parameters for a given block
          tmp_srv_sel_vec = ln_srv_fixed_sel_pars[r,,srv_sel_blk_idx,s,sf]

          # Compute selectivity functional form
          tmp_sel = Get_Selex(Selex_Model = tmp_srv_sel_model, # selectivity model
                              TimeVary_Model = cont_tv_srv_sel[r,sf], # time varying model
                              ln_Pars = tmp_srv_sel_vec, # fixed effect selectivity parameters
                              ln_seldevs = ln_srvsel_devs[,,,,sf, drop = FALSE], # Selectivity deviations
                              Region = r, # region index
                              Year = y, # year index
                              Bin = selex_bins, # bin vector
                              Sex = s # sex index
          )

          # Calculate selectivity
          if(Selex_Type == 0) srv_sel[r,y,,s,sf] = tmp_sel # age-based selectivity
          if(Selex_Type == 1) {
            srv_sel_l[r,y,,s,sf] = tmp_sel # input into length-based survey selectivity
            srv_sel[r,y,,s,sf] = tmp_sel %*% tmp_sizeage # length-based selectivity (dot product of size age transition)
          }

        } # end s loop
      } # end sf loop
    } # end y loop
  } # end r loop

  ## Mortality ---------------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {

      for(seas in 1:n_seas) {
        # Fishing Mortality at Age
        for(f in 1:n_fish_fleets) {
          if(UseCatch[r,y,seas,f] == 0) {
            Fmort[r,y,seas,f] = 0 # Set F to zero when no catch data
          } else {
            Fmort[r,y,seas,f] = exp(ln_F_mean[r,seas,f] + ln_F_devs[r,y,seas,f]) # If catch is region-specific
          }
          FAA[r,y,seas,,,f] = Fmort[r,y,seas,f] * fish_sel[r,y,,,f] # Fishing mortality at age
        } # f loop

        for(a in 1:n_ages) for(s in 1:n_sexes) ZAA[r,y,seas,a,s] = sum(FAA[r,y,seas,a,s,]) + (natmort[r,y,a,s] * seasdur[seas]) # Total Mortality at age
      } # end seas loop

    } # end y loop
  } # end r loop


  ## Recruitment Transformations and Bias Ramp (Methot and Taylor) -------------------------------
  ### Parameter Transformations -----------------------------------------------
  # Mean or virgin recruitment proportions
  if(n_regions > 1) {
    Rec_trans_prop = c(0, Rec_prop) # set up vector for transformation
    Rec_trans_prop = exp(Rec_trans_prop) / sum(exp(Rec_trans_prop)) # do multinomial logit to get recruitment proportions
  } else Rec_trans_prop = 1

  # Global recruitment scalar
  R0 = exp(ln_global_R0) # exponentiate

  # Steepness
  h_trans = 0.2 + (1 - 0.2) * RTMB::plogis(steepness_h) # bound steepness between 0.2 and 1

  # Recruitment SD
  sigmaR2_early = exp(ln_sigmaR[1])^2 # recruitment variability for early period
  sigmaR2_late = exp(ln_sigmaR[2])^2 # recruitment variability for late period

  # Recruitment sex-ratio
  if(n_sexes == 2) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        sexratio_blk_idx = sexratio_blocks[r,y] # extract sex ratio block
        sexratio_f = RTMB::plogis(sexratio_pars[r,sexratio_blk_idx]) # get female recruitment sex-ratio
        sexratio[r,y,] = c(sexratio_f, 1 - sexratio_f) # input total sex ratio
      } # end y loop
    } # end r loop
  } else sexratio[] = 1 # set recruitment sex ratio at 1

  ### Bias ramp ---------------------------------------------------------------
  if (do_rec_bias_ramp == 0) {
    bias_ramp = rep(1, n_est_rec_devs) # don't do bias ramp, set values to 1
  } else if (do_rec_bias_ramp == 1) {

    bias_ramp = rep(0, n_est_rec_devs) # set up bias ramp values

    # setup bias ramp year ranges
    years = 1:n_est_rec_devs # years for indexing
    range1 = which(years >= bias_year[1] & years < bias_year[2])  # ascending limb
    range2 = which(years >= bias_year[2] & years < bias_year[3])  # full bias correction
    range3 = which(years >= bias_year[3] & years < bias_year[4])  # descending limb

    # Apply bias ramp to the different ramp year ranges
    if (length(range1) > 0) bias_ramp[range1] = (years[range1] - bias_year[1]) / (bias_year[2] - bias_year[1]) # ascending limb
    if (length(range2) > 0) bias_ramp[range2] = 1 # full bias correction
    if (length(range3) > 0) bias_ramp[range3] = 1 - ((years[range3] - bias_year[3]) / (bias_year[4] - bias_year[3])) # descending limb

    bias_ramp = bias_ramp * max_bias_ramp_fct # scale bias ramp by a factor

  } # end if doing bias ramp

  ## Initial Age Structure ---------------------------------------------------

  # Get initial fished NAA
  Init_Fished_NAA = Get_Init_NAA(
    init_age_strc = init_age_strc, # initial age structure
    init_iter = n_ages * 5, # if init_age_strc == 0, number of iterations to run
    n_regions = n_regions, # regions
    n_sexes = n_sexes, # sexes
    n_ages = n_ages, # ages
    n_seas = n_seas, # seasons
    seasdur = seasdur, # seasonal duration
    natmort = array(natmort[,1,,], dim = c(n_regions, n_ages, n_sexes)), # natural mortality in first year
    init_F = init_F, # initial F applied
    fish_sel = array(fish_sel[,1,,,], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)), # fishery selectivity in first year
    R0_r = R0 * Rec_trans_prop, # regional mean or virgin recruitment
    sexratio = array(sexratio[,1,], dim = c(n_regions, n_sexes)), # sex ratio in first year
    Movement = array(Movement[,,1,,,], dim = c(n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
    do_recruits_move = do_recruits_move, # whether recruits move
    ln_InitDevs = ln_InitDevs # initial deviations
  )

  # Get initial unfished NAA
  Init_Unfished_NAA = Get_Init_NAA(
    init_age_strc = init_age_strc, # initial age structure
    init_iter = n_ages * 5, # if init_age_strc == 0, number of iterations to run
    n_regions = n_regions, # regions
    n_sexes = n_sexes, # sexes
    n_ages = n_ages, # ages
    n_seas = n_seas, # seasons
    seasdur = seasdur, # seasonal duration
    natmort = array(natmort[,1,,], dim = c(n_regions, n_ages, n_sexes)), # natural mortality in first year
    init_F = rep(0, n_seas), # initial F applied (0 for unfished)
    fish_sel = array(fish_sel[,1,,,], dim = c(n_regions, n_ages, n_sexes, n_fish_fleets)), # fishery selectivity in first year
    R0_r = R0 * Rec_trans_prop, # regional mean or virgin recruitment
    sexratio = array(sexratio[,1,], dim = c(n_regions, n_sexes)), # sex ratio in first year
    Movement = array(Movement[,,1,,,], dim = c(n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
    do_recruits_move = do_recruits_move, # whether recruits move
    ln_InitDevs = ln_InitDevs # initial deviations
  )

  # Input into model arrays (first year and season)
  NAA[,1,1,,] = Init_Fished_NAA
  NAA0[,1,1,,] = Init_Unfished_NAA

  ## Population Projection ---------------------------------------------------
  for(y in 1:n_yrs) {

    ### Annual Recruitment ------------------------------------------------------
    # Get Deterministic Recruitment
    tmp_Det_Rec = Get_Det_Recruitment(recruitment_model = rec_model,
                                      recruitment_dd = rec_dd,
                                      R0 = R0,
                                      Rec_Prop = Rec_trans_prop,
                                      h = h_trans,
                                      n_ages = n_ages,
                                      n_regions = n_regions,
                                      # Note: Using first year and female quantities to compute unfished SSB0
                                      sexratio_f = if(n_sexes == 1) rep(0.5, n_regions) else sexratio[,1,1],
                                      WAA = array(WAA[,1,,1], dim = c(n_regions, n_ages)),
                                      MatAA = array(MatAA[,1,,1], dim = c(n_regions, n_ages)),
                                      natmort = array(natmort[,1,,1], dim = c(n_regions, n_ages)),
                                      Movement = array(Movement[,,1,,1], dim = c(n_regions, n_regions, n_ages)),
                                      do_recruits_move = do_recruits_move,
                                      t_spawn = t_spawn,
                                      SSB_vals = SSB,
                                      y = y,
                                      n_seas = n_seas,
                                      spawn_seas = spawn_seas,
                                      seasdur = seasdur,
                                      rec_lag = rec_lag,
                                      init_F = init_F, # initF for dominant fleet
                                      fish_sel = array(fish_sel[,1,,1,1], dim = c(n_regions, n_ages)) # uses dominant fleet
    )

    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        if(y < sigmaR_switch) NAA[r,y,1,1,s] = tmp_Det_Rec[r] * exp(ln_RecDevs[r,y] - (sigmaR2_early/2 * bias_ramp[y])) * sexratio[r,y,s] # early period recruitment
        if(y >= sigmaR_switch && y <= n_est_rec_devs) NAA[r,y,1,1,s] = tmp_Det_Rec[r] * exp(ln_RecDevs[r,y] - (sigmaR2_late/2 * bias_ramp[y])) * sexratio[r,y,s] # late period recruitment
        # Dealing with terminal year recruitment
        if(y > n_est_rec_devs) NAA[r,y,1,1,s] = tmp_Det_Rec[r] * sexratio[s] # mean recruitment in terminal year (not estimate last year rec dev)
      } # end s loop
      Rec[r,y] = sum(NAA[r,y,1,1,]) # get annual recruitment container here
      NAA0[r,y,1,1,] = NAA[r,y,1,1,] # populate unfished NAA
    } # end r loop

    ### Movement ----------------------------------------------------------------
    for(seas in 1:n_seas) {
      if(n_regions > 1) {

        # Record values prior to movement
        NAA_bef = NAA[,y,seas,,]

        # Recruits don't move
        if(do_recruits_move == 0) {
          # Apply movement after ageing processes - start movement at age 2
          for(a in 2:n_ages) {
            for(s in 1:n_sexes) {
              NAA[,y,seas,a,s] = t(NAA[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # Fished
              NAA0[,y,seas,a,s] = t(NAA0[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # Unfished
            } # end s loop
          } # end a loop
        } # end if recruits don't move

        # Recruits move here
        if(do_recruits_move == 1) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              NAA[,y,seas,a,s] = t(NAA[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # Fished
              NAA0[,y,seas,a,s] = t(NAA0[,y,seas,a,s]) %*% Movement[,,y,seas,a,s] # Unfished
            } # end s loop
          } # end a loop
        } # end if

        # Record values after movement
        NAA_aft = NAA[,y,seas,,]

      } # only compute if spatial

      ### Mortality and Ageing ------------------------------------------------------
      if(seas < n_seas) {
        # within year / seasonal mortality
        NAA[,y,seas+1,1:n_ages,] = NAA[,y,seas,1:n_ages,] * exp(-ZAA[,y,seas,1:n_ages,])
        NAA0[,y,seas+1,1:n_ages,] = NAA0[,y,seas,1:n_ages,] * exp(-(natmort[,y,1:n_ages,] * seasdur[seas]))
      } else {
        # age advancement and enter into first season of next year
        # Fished
        NAA[,y+1,1,2:n_ages,] = NAA[,y,n_seas,1:(n_ages-1),] * exp(-ZAA[,y,seas,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group
        NAA[,y+1,1,n_ages,] = NAA[,y+1,1,n_ages,] + NAA[,y,seas,n_ages,] * exp(-ZAA[,y,seas,n_ages,]) # Acuumulate plus group
        # Unfished
        NAA0[,y+1,1,2:n_ages,] = NAA0[,y,n_seas,1:(n_ages-1),] * exp(-natmort[,y,1:(n_ages-1),] * seasdur[n_seas]) # Exponential mortality for individuals not in plus group
        NAA0[,y+1,1,n_ages,] = NAA0[,y+1,1,n_ages,] + NAA0[,y,n_seas,n_ages,] * exp(-natmort[,y,n_ages,] * seasdur[n_seas]) # Acuumulate plus group
      }

      ### Compute Biomass Quantities ----------------------------------------------
      if(seas == spawn_seas) {
        Total_Biom[,y] = apply(NAA[,y,spawn_seas,,,drop = FALSE] * WAA[,y,spawn_seas,,,drop = FALSE] * exp(-ZAA[,y,spawn_seas,,,drop = FALSE] * t_spawn), 1, sum) # Total biomass
        SSB[,y] = apply(NAA[,y,spawn_seas,,1,drop = FALSE] * WAA[,y,spawn_seas,,1,drop = FALSE] * MatAA[,y,spawn_seas,,1,drop = FALSE] * exp(-ZAA[,y,spawn_seas,,1,drop = FALSE] * t_spawn), 1, sum)

        # Get dynamic B0
        SSB0_array <- NAA0[,y,spawn_seas,,1,drop = FALSE] *  WAA[,y,spawn_seas,,1,drop = FALSE] * MatAA[,y,spawn_seas,,1,drop = FALSE]
        mort_spawn <- exp(-natmort[,y,,1,drop = FALSE] * t_spawn * seasdur[seas])
        mort_spawn <- array(mort_spawn, dim = dim(SSB0_array) ) # coerce array
        Dynamic_SSB0[,y] <- apply(SSB0_array * mort_spawn, 1, sum) # Dynamic B0

        # If single sex model, multiply SSB calculations by 0.5
        if(n_sexes == 1) {
          SSB[,y] = SSB[,y] * 0.5
          Dynamic_SSB0[,y] = Dynamic_SSB0[,y] * 0.5
        }
      }

    } # end seas loop
  } # end y loop

  # Get aggregated SSB values
  Aggregated_SSB = colSums(SSB)
  Dynamic_Aggregated_SSB0 = colSums(Dynamic_SSB0)

  ## Fishery Observation Model -----------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        fish_q_blk_idx = fish_q_blocks[r,y,f] # get time-block catchability index
        fish_q[r,y,f] = exp(ln_fish_q[r,fish_q_blk_idx,f]) # Input into fishery catchability container

        for(seas in 1:n_seas) {
          CAA[r,y,seas,,,f] = FAA[r,y,seas,,,f] / ZAA[r,y,seas,,] * NAA[r,y,seas,,] * (1 - exp(-ZAA[r,y,seas,,])) # Catch at age (Baranov's)

          if(fit_lengths == 1) {
            for(s in 1:n_sexes) {
              CAL[r,y,seas,,s,f] = SizeAgeTrans[r,y,seas,,,s] %*% CAA[r,y,seas,,s,f] # Catch at length
            } # end s loop
          } # fitting lengths

          # Get catch
          if(catch_units[f] == 0) PredCatch[r,y,seas,f] = sum(CAA[r,y,seas,,,f]) # abundance
          if(catch_units[f] == 1) PredCatch[r,y,seas,f] = sum(CAA[r,y,seas,,,f] * WAA_fish[r,y,seas,,,f]) # biomass

          # Get fishery index
          if(fish_idx_type[f] == 0) PredFishIdx[r,y,seas,f] = fish_q[r,y,f] * sum(NAA[r,y,seas,,] * fish_sel[r,y,,,f]) # abundance
          if(fish_idx_type[f] == 1) PredFishIdx[r,y,seas,f] = fish_q[r,y,f] * sum(NAA[r,y,seas,,] * fish_sel[r,y,,,f] * WAA_fish[r,y,seas,,,f]) # biomass
        } # end seas loop

      } # end f loop
    } # end y loop
  } # end r loop

  ## Survey Observation Model ------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {

        srv_q_blk_idx = srv_q_blocks[r,y,sf] # get time-block catchability index
        srv_q[r,y,sf] = exp(ln_srv_q[r,srv_q_blk_idx,sf]) # Input into survey catchability container
        if(do_srv_q_cov == 1) srv_q[r,y,sf] = srv_q[r,y,sf] * exp(sum(srv_q_cov[r,y,sf,] * srv_q_coeff[r,sf,])) # adding covariate effects

        for(seas in 1:n_seas) {
          SrvIAA[r,y,seas,,,sf] = NAA[r,y,seas,,] * srv_sel[r,y,,,sf] * exp(-t_srv[r,seas,sf] * ZAA[r,y,seas,,]) # Survey index at age

          if(fit_lengths == 1) {
            for(s in 1:n_sexes) {
              SrvIAL[r,y,seas,,s,sf] = SizeAgeTrans[r,y,seas,,,s] %*% SrvIAA[r,y,seas,,s,sf] # Survey index at length
            } # end s loop
          } # fitting lengths

          # Get predicted survey index
          if(srv_idx_type[sf] == 0) PredSrvIdx[r,y,seas,sf] = srv_q[r,y,sf] * sum(SrvIAA[r,y,seas,,,sf]) # abundance
          if(srv_idx_type[sf] == 1) PredSrvIdx[r,y,seas,sf] = srv_q[r,y,sf] * sum(SrvIAA[r,y,seas,,,sf] * WAA_srv[r,y,seas,,,sf]) # biomass

        } # end seas loop

      } # end sf loop
    } # end y loop
  } # end r loop


  ## Tagging Observation Model -----------------------------------------------
  if(UseTagging == 1) {

    # Set up tag reporting rates
    for(r in 1:n_regions) {
      TagRep_blk_idx = Tag_Reporting_blocks[r,]  # Get all blocks for this region
      Tag_Reporting[r,] = RTMB::plogis(Tag_Reporting_Pars[r,TagRep_blk_idx])  # inverse logit transform
    } # end r loop

    # Transform tagging parameters here
    Tag_Shed = exp(ln_Tag_Shed) # Transform shedding rates
    Init_Tag_Mort = exp(ln_Init_Tag_Mort) # Transform initial mortality

    # Extract out indexing for tag cohorts
    tr_vec = tag_release_indicator[,1] # tag release region vector
    ty_vec = tag_release_indicator[,2] # tag release year vector
    tseas_vec = tag_release_indicator[,3] # tag release season vector

    for(rseas in 1:n_seas) {
      for(tc in 1:n_tag_cohorts) {

        tr = tr_vec[tc] # extract tag release region
        ty = ty_vec[tc] # extract tag release year
        tseas = tseas_vec[tc] # extract tag release season

        for(ry in 1:min(max_tag_liberty, n_yrs - ty + 1)) {

          y = ty + ry - 1 # Get index for actual year in the model (instead of tag year)

          # Get total mortality
          tmp_Z = Get_Tagging_Mortality(tag_selex = tag_selex, tag_natmort = tag_natmort,
                                        Fmort = Fmort, natmort = natmort, fish_sel = fish_sel,
                                        Tag_Shed = Tag_Shed, n_regions = n_regions, n_ages = n_ages,
                                        seas = rseas, seasdur = seasdur,
                                        n_sexes = n_sexes, n_fish_fleets = n_fish_fleets, y = y, what = "Z")

          # Get fishing mortality
          tmp_F = Get_Tagging_Mortality(tag_selex = tag_selex, tag_natmort = tag_natmort,
                                        Fmort = Fmort, natmort = natmort, fish_sel = fish_sel,
                                        Tag_Shed = Tag_Shed, n_regions = n_regions, n_ages = n_ages,
                                        seas = rseas, seasdur = seasdur,
                                        n_sexes = n_sexes, n_fish_fleets = n_fish_fleets, y = y, what = "F")

          # Handle tagging dynamics
          # Discount with tagging time (t_tagging) if it doesn't happen at the start of the season / year
          if(ry == 1 && rseas == tseas) {
            if(t_tagging != 1) tmp_Z = tmp_Z * t_tagging
            # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
            Tags_Avail[1, rseas, tc, tr, , ] = Tagged_Fish[tc, , ] * exp(-exp(ln_Init_Tag_Mort))
          }

          # Move tagged fish around (skip only in first release year + tagging season when tagging occurs mid-season)
          if(t_tagging == 1 || ry != 1 || rseas != tseas) {
            # Movement of tag cohorts
            if(do_recruits_move == 0) {
              for(a in 2:n_ages) for(s in 1:n_sexes) {
                Tags_Avail[ry, rseas, tc, , a, s] =
                  t(Tags_Avail[ry, rseas, tc, , a, s]) %*%
                  Movement[, , y, rseas, a, s]
              } # end s loop
            } else { # if recruits move
              for(a in 1:n_ages) for(s in 1:n_sexes) {
                Tags_Avail[ry, rseas, tc, , a, s] =
                  t(Tags_Avail[ry, rseas, tc, , a, s]) %*%
                  Movement[, , y, rseas, a, s]
              } # end s loop
            } # end else
          } # end if

          # Mortality
          if(rseas < n_seas) {
            # within year mortality / seasonal mortality
            Tags_Avail[ry,rseas+1,tc,,1:n_ages,] = Tags_Avail[ry,rseas,tc,,1:n_ages,] * exp(-tmp_Z[,1,1:n_ages,]) # not in plus group
          } else {
            # Mortality and ageing of tagged fish
            Tags_Avail[ry+1,1,tc,,2:n_ages,] = Tags_Avail[ry,n_seas,tc,,1:(n_ages-1),] * exp(-tmp_Z[,1,1:(n_ages-1),]) # not in plus group
            Tags_Avail[ry+1,1,tc,,n_ages,] = Tags_Avail[ry+1,1,tc,,n_ages,] + (Tags_Avail[ry,n_seas,tc,,n_ages,] * exp(-tmp_Z[,1,n_ages,])) # accumulate plus group
          }

          # Get predicted recaptures
          Pred_Tag_Recap[ry,rseas,tc,,,] = Tag_Reporting[,y] * (tmp_F[,1,,] / tmp_Z[,1,,]) * Tags_Avail[ry,rseas,tc,,,] * (1 - exp(-tmp_Z[,1,,]))

        } # end ry loop
      } # end tc loop
    } # end rseas loop

  } # end if for using tagging data


  # Likelihood Equations -------------------------------------------------------------
  ## Fishery Likelihoods -----------------------------------------------------
  ### Fishery Catches ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(r in 1:n_regions) {

        for(seas in 1:n_seas) {
          if(UseCatch[r,y,seas,f] == 1) {
            Catch_nLL[r,y,seas,f] = UseCatch[r,y,seas,f] -1 * RTMB::dnorm(log(ObsCatch[r,y,seas,f]),
                                                                          log(PredCatch[r,y,seas,f]),
                                                                          exp(ln_sigmaC[r,y,seas,f]), TRUE)
          } # if no 0s for fishery catches

        } # end seas loop
      } # end r loop
    } # end f loop
  } # end y loop


  ### Fishery Indices ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(r in 1:n_regions) {
      for(f in 1:n_fish_fleets) {
        for(seas in 1:n_seas) {

          if(UseFishIdx[r,y,seas,f] == 1) {
            FishIdx_nLL[r,y,seas,f] = -1 * RTMB::dnorm(log(ObsFishIdx[r,y,seas,f] + addtofishidx),
                                                       log(PredFishIdx[r,y,seas,f] + addtofishidx),
                                                       ObsFishIdx_SE[r,y,seas,f], TRUE)
          } # if we have fishery indices

        } # end seas loop
      } # end f loop
    } # end r loop
  } # end y loop

  ### Fishery Compositions ------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      for(seas in 1:n_seas) {
        # Fishery Age Compositions
        if(sum(UseFishAgeComps[,y,seas,f]) >= 1) {
          FishAgeComps_nLL[,y,seas,,f] = Get_Comp_Likelihoods(

            # Expected and Observed values
            Exp = CAA[,y,seas,,,f],
            Obs = ObsFishAgeComps[,y,seas,,,f],

            # Input sample size and multinomial weight
            ISS = ISS_FishAgeComps[,y,seas,,f],
            Wt_Mltnml = Wt_FishAgeComps[,y,seas,,f],
            # Composition and Likelihood Type
            Comp_Type = FishAgeComps_Type[y,f],
            Likelihood_Type = FishAgeComps_LikeType[f],

            # overdispersion pars, Number of sexes, regions, age or length comps, and ageing error
            ln_theta = ln_FishAge_theta[,,f],
            ln_theta_agg = ln_FishAge_theta_agg[f],
            LN_corr_pars = FishAge_corr_pars[,,f,],
            LN_corr_pars_agg = FishAge_corr_pars_agg[f],
            n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0,
            AgeingError = AgeingError[y,,],
            use = UseFishAgeComps[,y,seas,f],
            n_model_bins = n_ages,
            n_obs_bins = dim(ObsFishAgeComps)[4],
            addtocomp = addtocomp
          )

        } # if we have fishery age comps

        # Fishery Length Compositions
        if(sum(UseFishLenComps[,y,seas,f]) >= 1 && fit_lengths == 1) {
          FishLenComps_nLL[,y,seas,,f] = Get_Comp_Likelihoods(

            # Expected and Observed values
            Exp = CAL[,y,seas,,,f],
            Obs = ObsFishLenComps[,y,seas,,,f],

            # Input sample size and multinomial weight
            ISS = ISS_FishLenComps[,y,seas,,f],
            Wt_Mltnml = Wt_FishLenComps[,y,seas,,f],

            # Composition and Likelihood Type
            Comp_Type = FishLenComps_Type[y,f],
            Likelihood_Type = FishLenComps_LikeType[f],

            # overdispersion, Number of sexes, regions age or length comps, and ageing error
            ln_theta = ln_FishLen_theta[,,f],
            ln_theta_agg = ln_FishLen_theta_agg[f],
            LN_corr_pars = FishLen_corr_pars[,,f,],
            LN_corr_pars_agg = FishLen_corr_pars_agg[f],
            n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1,
            AgeingError = NA,
            use = UseFishLenComps[,y,seas,f],
            n_model_bins = n_lens,
            n_obs_bins = dim(ObsFishLenComps)[4],
            addtocomp = addtocomp
          )

        } # if we have fishery length comps
      } # end seas loop

    } # end f loop
  } # end y loop


  ## Survey Likelihoods ------------------------------------------------------
  ### Survey Indices ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(r in 1:n_regions) {
      for(sf in 1:n_srv_fleets) {
        for(seas in 1:n_seas) {

          if(UseSrvIdx[r,y,seas,sf] == 1) {
            SrvIdx_nLL[r,y,seas,sf] = -1 * RTMB::dnorm(log(ObsSrvIdx[r,y,seas,sf] + addtosrvidx),
                                                       log(PredSrvIdx[r,y,seas,sf] + addtosrvidx),
                                                       ObsSrvIdx_SE[r,y,seas,sf], TRUE)
          } # if we have survey indices

        } # end seas loop
      } # end sf loop
    } # end r loop
  } # end y loop

  ### Survey Compositions ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(sf in 1:n_srv_fleets) {
      for(seas in 1:n_seas) {

        # Survey Age Compositions
        if(sum(UseSrvAgeComps[,y,seas,sf]) >= 1) {
          SrvAgeComps_nLL[,y,seas,,sf] = Get_Comp_Likelihoods(

            # Expected and Observed values
            Exp = SrvIAA[,y,seas,,,sf],
            Obs = ObsSrvAgeComps[,y,seas,,,sf],

            # Input sample size and multinomial weight
            ISS = ISS_SrvAgeComps[,y,seas,,sf],
            Wt_Mltnml = Wt_SrvAgeComps[,y,seas,,sf],

            # Composition and Likelihood Type
            Comp_Type = SrvAgeComps_Type[y,sf],
            Likelihood_Type = SrvAgeComps_LikeType[sf],

            # overdispersion, Number of sexes, regions, age or length comps, and ageing error
            ln_theta = ln_SrvAge_theta[,,sf],
            ln_theta_agg = ln_SrvAge_theta_agg[sf],
            LN_corr_pars = SrvAge_corr_pars[,,sf,],
            LN_corr_pars_agg = SrvAge_corr_pars_agg[sf],
            n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0,
            AgeingError = AgeingError[y,,],
            use = UseSrvAgeComps[,y,seas,sf],
            n_model_bins = n_ages,
            n_obs_bins = dim(ObsSrvAgeComps)[4],
            addtocomp = addtocomp
          )

        } # if we have survey age comps

        # Survey Length Compositions
        if(sum(UseSrvLenComps[,y,seas,sf]) >= 1 && fit_lengths == 1) {
          SrvLenComps_nLL[,y,seas,,sf] = Get_Comp_Likelihoods(

            # Expected and Observed values
            Exp = SrvIAL[,y,seas,,,sf],
            Obs = ObsSrvLenComps[,y,seas,,,sf],

            # Input sample size and multinomial weight
            ISS = ISS_SrvLenComps[,y,seas,,sf],
            Wt_Mltnml = Wt_SrvLenComps[,y,seas,,sf],

            # Composition and Likelihood Type
            Comp_Type = SrvLenComps_Type[y,sf],
            Likelihood_Type = SrvLenComps_LikeType[sf],

            # overdispersion, Number of sexes, regions, age or length comps, and ageing error
            ln_theta = ln_SrvLen_theta[,,sf],
            ln_theta_agg = ln_SrvLen_theta_agg[sf],
            LN_corr_pars = SrvLen_corr_pars[,,sf,],
            LN_corr_pars_agg = SrvLen_corr_pars_agg[sf],
            n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1,
            AgeingError = NA,
            use = UseSrvLenComps[,y,seas,sf],
            n_model_bins = n_lens,
            n_obs_bins = dim(ObsSrvLenComps)[4],
            addtocomp = addtocomp
          )

        } # if we have survey length comps

      } # end seas loop
    } # end sf loop
  } # end y loop


  ## Tag Likelihoods ---------------------------------------------------------
  if(UseTagging == 1) {
    for(tc in 1:n_tag_cohorts) {

      # set up tagging cohort indexing
      tr = tr_vec[tc] # extract tag release region
      ty = ty_vec[tc] # extract tag release year
      tseas = tseas_vec[tc] # extract tag release season

      for(ry in mixing_period:min(max_tag_liberty, n_yrs - ty + 1)) { # loop through recapture years
        for(rseas in 1:n_seas) { # loop through recapture seasons

          for(r in 1:n_regions) {
            for(a in 1:n_move_age_tag_pool) {
              for(s in 1:n_move_sex_tag_pool) {

                move_age_pool_idx = move_age_tag_pool[[a]] # extract movement age pool indices
                move_sex_pool_idx = move_sex_tag_pool[[s]] # extract movement sex pool indices

                # Poisson likelihood
                if(Tag_LikeType == 0) {
                  Tag_nLL[ry,rseas,tc,r,1,1] = Tag_nLL[ry,rseas,tc,r,1,1] + -dpois_noint(sum(Obs_Tag_Recap[ry,rseas,tc,r,move_age_pool_idx,move_sex_pool_idx] + addtotag),
                                                                             sum(Pred_Tag_Recap[ry,rseas,tc,r,move_age_pool_idx,move_sex_pool_idx] + addtotag),
                                                                             give_log = TRUE)
                } # end if poisson likelihood

                # Negative binomial likelihood
                if(Tag_LikeType == 1) {
                  log_mu = log(sum(Pred_Tag_Recap[ry,rseas,tc,r,move_age_pool_idx,move_sex_pool_idx] + addtotag)) # log mu
                  log_var_minus_mu = 2 * log_mu - ln_tag_theta # log var minus mu
                  Tag_nLL[ry,rseas,tc,r,1,1] = Tag_nLL[ry,rseas,tc,r,1,1] + -dnbinom_robust_noint(x = sum(Obs_Tag_Recap[ry,rseas,tc,r,move_age_pool_idx,move_sex_pool_idx] + addtotag),
                                                                                                  log_mu = log_mu, log_var_minus_mu = log_var_minus_mu, give_log = TRUE)
                } # end if for negative binomial likelihood

              } # end s loop
            } # end a loop
          } # end r loop

          # Release Conditioned for Multinomial or Dirichlet-Multinomial
          if(Tag_LikeType %in% c(2, 4)) {

            # Temporary vectors for recaptured individuals
            tmp_pred_c_all = vector()
            tmp_obs_c_all = vector()

            # number of tags released for a given tag cohort
            tmp_n_tags_released = sum(Tagged_Fish[tc,,] + addtotag)

            # Loop through age and sex pooling and combine vectors into the correct format
            for(a in 1:n_move_age_tag_pool) {
              for(s in 1:n_move_sex_tag_pool) {
                move_age_pool_idx = move_age_tag_pool[[a]] # extract movement age pool indices
                move_sex_pool_idx = move_sex_tag_pool[[s]] # extract movement sex pool indices

                # Pool observed and expected if any pooling
                for (r in 1:n_regions) {
                  pred_val = sum(Pred_Tag_Recap[ry, rseas, tc, r, move_age_pool_idx, move_sex_pool_idx] + addtotag) # sum across age and sex groups
                  obs_val  = sum(Obs_Tag_Recap[ry, rseas, tc, r, move_age_pool_idx, move_sex_pool_idx] + addtotag) # sum across age and sex groups
                  tmp_pred_c_all = c(tmp_pred_c_all, pred_val) # combine predicted recaptures for a given age sex pooled group
                  tmp_obs_c_all  = c(tmp_obs_c_all,  obs_val) # combine observed recaptures for a given age sex pooled group
                } # end r loop

              } # end a loop
            } # end s loop

            # Normalize observed and predicted recaptures
            tmp_pred_c_all = tmp_pred_c_all / tmp_n_tags_released
            tmp_obs_c_all = tmp_obs_c_all / tmp_n_tags_released

            # Add in observed and predicted non-recaptures
            tmp_pred = c(tmp_pred_c_all, 1 - sum(tmp_pred_c_all))
            tmp_obs = c(tmp_obs_c_all, 1 - sum(tmp_obs_c_all))

            if(Tag_LikeType == 2) Tag_nLL[ry,rseas,tc,1,1,1] = -tmp_n_tags_released * sum((tmp_obs) * log(tmp_pred)) # multinomial
            if(Tag_LikeType == 4) Tag_nLL[ry,rseas,tc,1,1,1] =  -1 * ddirmult(obs = tmp_obs, pred = tmp_pred, Ntotal = tmp_n_tags_released, ln_theta = ln_tag_theta, TRUE) # Dirichlet Multinomial

          } # end if release conditioned

          # Recapture Conditioned (Multinomial or Dirichlet-Multinomial)
          if(Tag_LikeType %in% c(3,5)) {
            # Temporary vectors for recaptured individuals
            tmp_pred_all = vector()
            tmp_obs_all = vector()

            # number of recaptures
            tmp_n_tags_recap = sum(Obs_Tag_Recap[ry,rseas,tc,,,] + addtotag)

            # Loop through age and sex pooling and combine vectors into the correct format
            for(a in 1:n_move_age_tag_pool) {
              for(s in 1:n_move_sex_tag_pool) {
                move_age_pool_idx = move_age_tag_pool[[a]] # extract movement age pool indices
                move_sex_pool_idx = move_sex_tag_pool[[s]] # extract movement sex pool indices

                for (r in 1:n_regions) {
                  pred_val = sum(Pred_Tag_Recap[ry, rseas, tc, r, move_age_pool_idx, move_sex_pool_idx] + addtotag) # sum across age and sex groups
                  obs_val  = sum(Obs_Tag_Recap[ry, rseas, tc, r, move_age_pool_idx, move_sex_pool_idx] + addtotag) # sum across age and sex groups
                  tmp_pred_all = c(tmp_pred_all, pred_val) # combine predicted recaptures for a given age sex pooled group
                  tmp_obs_all  = c(tmp_obs_all,  obs_val) # combine observed recaptures for a given age sex pooled group
                } # end r loop

              } # end a loop
            } # end s loop

            # Normalize observed and predicted recaptures
            tmp_pred_all = tmp_pred_all / sum(tmp_pred_all)
            tmp_obs_all = tmp_obs_all / tmp_n_tags_recap

            if(Tag_LikeType == 3) Tag_nLL[ry,rseas,tc,1,1,1] = -1 * tmp_n_tags_recap * sum(((tmp_obs_all) * log(tmp_pred_all))) # Multinomial
            if(Tag_LikeType == 5) Tag_nLL[ry,rseas,tc,1,1,1] =  -1 * ddirmult(obs = tmp_obs_all, pred = tmp_pred_all, Ntotal = tmp_n_tags_recap, ln_theta = ln_tag_theta, TRUE) # Dirichlet Multinomial

          } # end if recapture conditioned

        } # end if rseas loop
      } # end ry loop

    } # end tc loop
  } # if we are using tagging data

  ## Priors and Penalties ----------------------------------------------------
  ### Fishing Mortality (Penalty) ---------------------------------------------
  if(Use_F_pen == 1) {
    for(f in 1:n_fish_fleets) {
      for(y in 1:n_yrs) {
        for(r in 1:n_regions) {
          for(seas in 1:n_seas) {

            if(UseCatch[r,y,seas,f] == 1) {
              Fmort_nLL[r,y,seas,f] = -RTMB::dnorm(ln_F_devs[r,y,seas,f], 0, exp(ln_sigmaF[r,seas,f]), TRUE)
            } # end if have catch

          } # end seas loop
        } # end r loop
      } # y loop
    } # f loop
  } #  if using fishing mortality penalty

  ### Selectivity (Penalty) ---------------------------------------------------
  for(r in 1:n_regions) {

    # Fishery Selectivity Deviations
    for(f in 1:n_fish_fleets) {

      if(cont_tv_fish_sel[r,f] > 0) {

        if(Selex_Type == 0) tmp_sel_vals = fish_sel[r,,,,f, drop = FALSE] # age-based selectivity
        if(Selex_Type == 1) tmp_sel_vals = fish_sel_l[r,,,,f, drop = FALSE] # length-based selectivity

        sel_nLL = sel_nLL + - Get_sel_PE_loglik(PE_model = cont_tv_fish_sel[r,f], # process error model
                                                PE_pars = fishsel_pe_pars[r,,,f, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_fishsel_devs[r,,,,f, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_fishsel_devs[r,,,,f, drop = FALSE],
                                                sel_vals = tmp_sel_vals,
                                                do_sel_pen = cont_tv_fish_sel_penalty
                                                )
      } # end if

      # Mean Standardizing to help with interpretability
      if(Selex_Type == 0) if(cont_tv_fish_sel[r,f] %in% 3:5) for(s in 1:n_sexes) fish_sel[r,,,s,f] = exp(log(fish_sel[r,,,s,f]) - log(mean(fish_sel[r,,,s,f]))) # age-based selectivity
      if(Selex_Type == 1) if(cont_tv_fish_sel[r,f] %in% 3:5) for(s in 1:n_sexes) fish_sel_l[r,,,s,f] = exp(log(fish_sel_l[r,,,s,f]) - log(mean(fish_sel_l[r,,,s,f]))) # length-based selectivity

    } # end f loop

    # Survey Selectivity Deviations
    for(sf in 1:n_srv_fleets) {

      if(cont_tv_srv_sel[r,sf] > 0) {

        if(Selex_Type == 0) tmp_sel_vals = srv_sel[r,,,,sf, drop = FALSE] # age-based selectivity
        if(Selex_Type == 1) tmp_sel_vals = srv_sel_l[r,,,,sf, drop = FALSE] # length-based selectivity

        sel_nLL = sel_nLL + - Get_sel_PE_loglik(PE_model = cont_tv_srv_sel[r,sf], # process error model
                                                PE_pars = srvsel_pe_pars[r,,,sf, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_srvsel_devs[r,,,,sf, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_srvsel_devs[r,,,,sf, drop = FALSE],
                                                sel_vals = tmp_sel_vals,
                                                do_sel_pen = cont_tv_srv_sel_penalty
                                                )
      } # end if

      # Mean Standardizing to help with interpretability
      if(Selex_Type == 0) if(cont_tv_srv_sel[r,sf] %in% 3:5) for(s in 1:n_sexes) srv_sel[r,,,s,sf] = exp(log(srv_sel[r,,,s,sf]) - log(mean(srv_sel[r,,,s,sf])))
      if(Selex_Type == 1) if(cont_tv_srv_sel[r,sf] %in% 3:5) for(s in 1:n_sexes) srv_sel_l[r,,,s,sf] = exp(log(srv_sel_l[r,,,s,sf]) - log(mean(srv_sel_l[r,,,s,sf])))

    } # end sf loop
  } # end r loop


  ### Selectivity (Prior) -----------------------------------------------------
  # Fishery selectivity parameters
  if(Use_fish_selex_prior == 1) {
    for(i in 1:nrow(fish_selex_prior)) {
      # Extract indices
      r = fish_selex_prior$region[i]
      p = fish_selex_prior$par[i]
      b = fish_selex_prior$block[i]
      s = fish_selex_prior$sex[i]
      f = fish_selex_prior$fleet[i]
      # Compute penalty / prior here
      sel_nLL = sel_nLL - RTMB::dnorm(ln_fish_fixed_sel_pars[r,p,b,s,f], log(fish_selex_prior$mu[i]), fish_selex_prior$sd[i], TRUE)
    } # end i loop
  } # end if using selex priors

  # Survey selectivity parameters
  if(Use_srv_selex_prior == 1) {
    for(i in 1:nrow(srv_selex_prior)) {
      # Extract indices
      r = srv_selex_prior$region[i]
      p = srv_selex_prior$par[i]
      b = srv_selex_prior$block[i]
      s = srv_selex_prior$sex[i]
      sf = srv_selex_prior$fleet[i]
      # Compute penalty / prior here
      sel_nLL = sel_nLL - RTMB::dnorm(ln_srv_fixed_sel_pars[r,p,b,s,sf], log(srv_selex_prior$mu[i]), srv_selex_prior$sd[i], TRUE)
    } # end i loop
  } # end if using selex priors

  ### Recruitment (Penalty) ----------------------------------------------------
  for(r in 1:n_regions) {
    # if initial age structure is stochastic
    if(equil_init_age_strc %in% c(1,2)) Init_Rec_nLL[r,] = -RTMB::dnorm(ln_InitDevs[r,], 0, exp(ln_sigmaR[1]), TRUE) # initial age structure penalty
    if(sigmaR_switch > 1) for(y in 1:(sigmaR_switch-1)) {
      Rec_nLL[r,y] = -RTMB::dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR[1]), TRUE)
      if(do_rec_bias_ramp == 1) Rec_nLL[r,y] = Rec_nLL[r,y] - (1 - 0.5 * bias_ramp[y]) * ln_sigmaR[1] # early period with bias ramp correction
    } # first y loop
    # Note that this penalizes the terminal year rec devs, which is estimated in this case
    for(y in sigmaR_switch:n_est_rec_devs) {
      Rec_nLL[r,y] = -RTMB::dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR[2]), TRUE)
      if(do_rec_bias_ramp == 1) Rec_nLL[r,y] = Rec_nLL[r,y] - (1 - 0.5 * bias_ramp[y]) * ln_sigmaR[2] # late period with bias ramp correction
    } # end second y loop
  } # end r loop


  ### Fishery Catchability (Prior) -----------------------------------------------
  if(Use_fish_q_prior == 1) {
    for(i in 1:nrow(fish_q_prior)) {
      # Extract indices
      r = fish_q_prior$region[i]
      b = fish_q_prior$block[i]
      f = fish_q_prior$fleet[i]
      # Compute penalty / prior here
      fish_q_nLL = fish_q_nLL - sum(RTMB::dnorm(ln_fish_q[r,b,f], log(fish_q_prior$mu[i]), fish_q_prior$sd[i], TRUE))
    } # end i loop
  } # end if using survey catchability prior

  ### Survey Catchability (Prior) -----------------------------------------------
  if(Use_srv_q_prior == 1) {
    for(i in 1:nrow(srv_q_prior)) {
      # Extract indices
      r = srv_q_prior$region[i]
      b = srv_q_prior$block[i]
      sf = srv_q_prior$fleet[i]
      # Compute penalty / prior here
      srv_q_nLL = srv_q_nLL - sum(RTMB::dnorm(ln_srv_q[r,b,sf], log(srv_q_prior$mu[i]), srv_q_prior$sd[i], TRUE))
    } # end i loop
  } # end if using survey catchability prior

  ## Natural Mortality (Prior) -----------------------------------------------
  if(Use_M_prior == 1) {
    for(i in 1:nrow(M_prior)) {
      # Extract indices
      r = M_prior$regionblk[i]
      b = M_prior$yearblk[i]
      a = M_prior$ageblk[i]
      s = M_prior$sexblk[i]
      # Compute prior
      M_nLL = M_nLL + -RTMB::dnorm(ln_M[r,b,a,s], log(M_prior$mu[i]), M_prior$sd[i], TRUE) # TMB likelihood
    }
  } # end if using natural mortality prior

  ### Steepness (Prior) -----------------------------------------------
  if(Use_h_prior == 1) {
    for(i in 1:nrow(h_prior)) {
      # Extract indicies
      r = h_prior$region[i]
      tmp_h_beta_pars = get_beta_scaled_pars(low = 0.2, high = 1, mu = h_prior$mu[i], sigma = h_prior$sd[i]) # get alpha and beta parameters
      tmp_h_trans = (h_trans[r] - 0.2) / (1 - 0.2) # transform random variable
      h_nLL = h_nLL - RTMB::dbeta(x = tmp_h_trans, shape1 = tmp_h_beta_pars[1], shape2 = tmp_h_beta_pars[2], log = TRUE) # penalize
    } # end i
  } # end if using steepness prior


  ### Movement Rates (Penalty) ------------------------------------------------
  if(cont_vary_movement > 0) {
    Movement_nLL = Movement_nLL + - Get_move_PE_loglik(PE_model = cont_vary_movement,
                                                       PE_pars = move_pe_pars,
                                                       move_devs = move_devs,
                                                       map_move_devs = map_move_devs,
                                                       do_recruits_move = do_recruits_move,
                                                       adjacency_collapsed = adjacency_collapsed,
                                                       move_type = move_type
    )
  }

  ### Movement Rates (Prior) ------------------------------------------------
  if(Use_Movement_Prior == 1) {
    for(i in 1:nrow(Movement_prior)) {
      region_from = Movement_prior$region_from[i] # region from
      y = Movement_prior$year[i] # year
      seas = Movement_prior$seas[i] # seas
      a = Movement_prior$age[i] # age
      s = Movement_prior$sex[i] # sex
      alpha = Movement_prior$alpha[[i]] # get prior values
      Movement_nLL = Movement_nLL - ddirichlet(x = Movement[region_from,,y,seas,a,s], alpha = alpha, log = TRUE) # dirichlet prior
    } # end i loop
  }

  ### Recruitment Proportions (Prior) -----------------------------------------
  if(Use_Rec_prop_Prior == 1) {
    Rec_prop_nLL = -ddirichlet(x = Rec_trans_prop, alpha = Rec_prop_prior, log = TRUE) # dirichlet prior
  }

  ### Tag Reporting Rate (Prior) --------------------------------------------
  if(Use_TagRep_Prior == 1) {
    for(i in 1:nrow(TagRep_Prior)) {
      # Extract indices
      r = TagRep_Prior$region[i]
      b = TagRep_Prior$block[i]
      TagRep_val = RTMB::plogis(Tag_Reporting_Pars[r,b]) # extract tag reporting rate value
      if(TagRep_Prior$type[i] == 0) {
        TagRep_nLL = TagRep_nLL - dbeta_symmetric(p_val = TagRep_val, p_ub = 1, p_lb = 0, p_prsd = TagRep_Prior$sd[i], log = TRUE) # penalize
      } # end if symmetric beta

      if(TagRep_Prior$type[i] == 1) {
        a = TagRep_Prior$mu[i] / (TagRep_Prior$sd[i] * TagRep_Prior$sd[i]) # alpha parameter
        b = (1 - TagRep_Prior$mu[i]) / (TagRep_Prior$sd[i] * TagRep_Prior$sd[i]) # beta parameter
        TagRep_nLL = TagRep_nLL -RTMB::dbeta(x = TagRep_val, shape1 = a, shape2 = b, log = TRUE) # penalize
      } # end if for full beta
    } # end i loop
  } # if use tag reporting prior

  # Apply likelihood weights here and compute joint negative log likelihood
  jnLL = sum(Wt_Catch * Catch_nLL) + # Catch likelihoods
         sum(Wt_FishIdx * FishIdx_nLL) + # Fishery Index likelihood
         sum(Wt_SrvIdx * SrvIdx_nLL) + # Survey Index likelihood
         sum(FishAgeComps_nLL) + # Fishery Age likelihood
         sum(FishLenComps_nLL) + # Fishery Length likelihood
         sum(SrvAgeComps_nLL) + # Survey Age likelihood
         sum(SrvLenComps_nLL) + # Survey Length likelihood
         (Wt_Tagging * sum(Tag_nLL)) + # Tagging likelihood
         (Wt_F * sum(Fmort_nLL)) + # Fishery Mortality Penalty
         (Wt_Rec * sum(Rec_nLL)) + # Recruitment Penalty
         (Wt_Rec * sum(Init_Rec_nLL)) + #  Initial Age Penalty
         sel_nLL + #  selectivity penalty
         M_nLL + # Natural Mortality Prior
         h_nLL + # Steepness Prior
         Movement_nLL + # movement Prior
         TagRep_nLL + # tag reporting rate Prior
         fish_q_nLL + # fishery q prior
         srv_q_nLL +  # survey q prior
         Rec_prop_nLL # recruitment proportion prior



  # Report Section ----------------------------------------------------------
  # Biological Processes
  RTMB::REPORT(R0)
  RTMB::REPORT(sexratio)
  RTMB::REPORT(Rec_trans_prop)
  RTMB::REPORT(h_trans)
  RTMB::REPORT(NAA)
  RTMB::REPORT(NAA0)
  RTMB::REPORT(NAA_bef)
  RTMB::REPORT(NAA_aft)
  RTMB::REPORT(ZAA)
  RTMB::REPORT(natmort)
  RTMB::REPORT(bias_ramp)
  RTMB::REPORT(Movement)
  RTMB::REPORT(Mrate)

  # Fishery Processes
  RTMB::REPORT(init_F)
  RTMB::REPORT(ln_sigmaC)
  RTMB::REPORT(Fmort)
  RTMB::REPORT(FAA)
  RTMB::REPORT(CAA)
  RTMB::REPORT(CAL)
  RTMB::REPORT(PredCatch)
  RTMB::REPORT(PredFishIdx)
  RTMB::REPORT(fish_sel)
  RTMB::REPORT(fish_q)

  # Survey Processes
  RTMB::REPORT(PredSrvIdx)
  RTMB::REPORT(srv_sel)
  RTMB::REPORT(srv_q)
  RTMB::REPORT(SrvIAA)
  RTMB::REPORT(SrvIAL)

  # Report length-based selectivity
  if(Selex_Type == 1) {
    RTMB::REPORT(fish_sel_l)
    RTMB::REPORT(srv_sel_l)
  }

  # Tagging Processes
  if(UseTagging == 1) {
    RTMB::REPORT(Pred_Tag_Recap)
    RTMB::REPORT(Tags_Avail)
    RTMB::REPORT(Tag_Reporting)
  }

  # Parameter Deviations
  RTMB::REPORT(ln_RecDevs)
  RTMB::REPORT(move_devs)
  RTMB::REPORT(ln_fishsel_devs)
  RTMB::REPORT(ln_srvsel_devs)

  # Likelihoods
  RTMB::REPORT(Catch_nLL)
  RTMB::REPORT(FishIdx_nLL)
  RTMB::REPORT(SrvIdx_nLL)
  RTMB::REPORT(FishAgeComps_nLL)
  RTMB::REPORT(SrvAgeComps_nLL)
  RTMB::REPORT(FishLenComps_nLL)
  RTMB::REPORT(SrvLenComps_nLL)
  RTMB::REPORT(M_nLL)
  RTMB::REPORT(Fmort_nLL)
  RTMB::REPORT(Rec_nLL)
  RTMB::REPORT(Init_Rec_nLL)
  RTMB::REPORT(Tag_nLL)
  RTMB::REPORT(h_nLL)
  RTMB::REPORT(fish_q_nLL)
  RTMB::REPORT(sel_nLL)
  RTMB::REPORT(srv_q_nLL)
  RTMB::REPORT(Movement_nLL)
  RTMB::REPORT(TagRep_nLL)
  RTMB::REPORT(Rec_prop_nLL)
  RTMB::REPORT(jnLL)

  # Report for derived quantities
  RTMB::REPORT(Total_Biom)
  RTMB::REPORT(SSB)
  RTMB::REPORT(Dynamic_SSB0)
  RTMB::REPORT(Aggregated_SSB)
  RTMB::REPORT(Dynamic_Aggregated_SSB0)
  RTMB::REPORT(Rec)

  # Report these in log space because can't be < 0
  RTMB::ADREPORT(log(Total_Biom))
  RTMB::ADREPORT(log(SSB))
  RTMB::ADREPORT(log(Dynamic_SSB0))
  RTMB::ADREPORT(log(Rec))
  RTMB::ADREPORT(log(Aggregated_SSB))
  RTMB::ADREPORT(log(Dynamic_Aggregated_SSB0))

  return(jnLL)
} # end function
