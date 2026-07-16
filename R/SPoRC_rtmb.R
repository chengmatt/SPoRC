# version 1 - (M.LH Cheng)
# Bridge model 23.5 from ADMB to RTMB
# Changed code to be more modular, accommodating any number of fishery and survey fleets
# Rectified errors in fitting to length composition data (normalize proportions at length after conversion from age-length matrix)
# Changed survey composition data to be calculated using survival midyear
# Added options for continuous time-varying selectivity
# Added options for TMB / R-like likelihoods (e.g., dnorm) instead of custom likelihoods
# Added options for dirichlet multinomial likelihood

# version 2 - (M.LH Cheng)
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
# Re-factored Initial Numbers at Age module

# version 5 - (M.LH Cheng & J.T Thorson)
# Refactored movement priors and movement setup
# Added in CTMC movement

# verison 6 - (M.LH Cheng)
# Re-coded tagging module to include fleet-specific dynamics
# Expanded model to include populations and season partitions (natal homing
# and seasonality dynamics)

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
  n_est_rec_devs = dim(ln_RecDevs)[3] # number of recruitment deviates estimated
  Rec = array(0, dim = c(n_pop, n_regions, n_yrs)) # Recruitment
  R0 = array(0, dim = c(n_pop, n_regions)) # R0 or mean recruitment
  sexratio = array(0, dim = c(n_pop, n_regions, n_yrs, n_sexes)) # recruitment sex ratio
  rec_region_prop = array(0, dim = c(n_pop, n_regions)) # recruitment regional apportionment
  rec_seas_prop = array(0, dim = c(n_pop, n_seas)) # recruitment seasonal apportionment
  stray_rate = array(0, dim = c(n_pop, n_yrs)) # stray rate

  # Population Dynamics
  NAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Numbers at age
  NAA_bef = array(data = 0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Numbers at age before movement
  NAA_aft = array(data = 0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Numbers at age after movement
  NAA0 = array(data = 0, dim = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes)) # Unfished Numbers at age
  ZAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes)) # Total mortality at age
  natmort = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_ages, n_sexes)) # natural mortality at age
  Total_Biom = array(0, dim = c(n_pop, n_regions, n_yrs)) # Total biomass
  SSB = array(0, dim = c(n_pop, n_regions, n_yrs)) # Spawning stock biomass
  eff_SSB = array(0, dim = c(n_pop, n_yrs)) # effective SSB for a given population
  Dynamic_SSB0 = array(0, dim = c(n_pop, n_regions, n_yrs)) # Dynamic Unfished Spawning stock biomass
  Aggregated_SSB = array(0, dim = c(n_pop, n_yrs)) # Aggregated Spawning stock biomass
  Dynamic_Aggregated_SSB0 = array(0, dim = c(n_pop, n_yrs)) # Dynamic Unfished Aggregated Spawning stock biomass

  # Movement Stuff
  Movement = array(data = 0, dim = c(n_pop, n_regions, n_regions, n_yrs + n_proj_yrs_devs, n_seas, n_ages, n_sexes)) # movement "matrix"
  n_conv_tag_pop_pool = length(conv_tag_pop_pool) # number of populations to pool for tagging data
  n_conv_tag_age_pool = length(conv_tag_age_pool) # number of ages to pool for tagging data
  n_conv_tag_sex_pool = length(conv_tag_sex_pool) # number of sexes to pool for tagging data

  # Tagging Stuff
  conv_tag_fish_avail = array(data = 0, dim = c(conv_tag_max_liberty + 1, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes)) # Tags availiable for recapture
  conv_tag_fish_reporting = array(data = 0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Tag reporting rate
  pred_conv_tag_fish_recap = array(data = 0, dim = c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets)) # predicted recaptures

  # Fishery Processes
  Fmort = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Fishing mortality scalar
  dmr = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Discard mortality rate
  tot_FAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Total Fishing mortality at age
  ret_FAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Retained Fishing mortality at age
  disc_FAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Discarded Fishing mortality at age
  CAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Retained Catch at age
  DAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Discarded Catch at age
  CAL = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets)) # Retained Catch at length
  DAL = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets)) # Discarded Catch at length
  PredCatch = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted retained catch (can be abundance or biomass)
  PredDiscard = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted discarded catch (can be abundance, biomass, or abdunance or biomass fraction of retained catch)
  PredFishIdx = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted fishery index
  fish_sel = array(data = 0, dim = c(n_pop, n_regions, n_yrs + n_proj_yrs_devs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Total Fishery selectivity
  fish_sel_l = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_lens, n_sexes, n_fish_fleets)) # Retained Fishery selectivity (lengths)
  ret_sel = array(data = 0, dim = c(n_pop, n_regions, n_yrs + n_proj_yrs_devs, n_seas, n_ages, n_sexes, n_fish_fleets)) # Fishery selectivity
  ret_sel_l = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_lens, n_sexes, n_fish_fleets)) # Fishery selectivity (lengths)
  fish_q = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery catchability

  # Survey Processes
  SrvIAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)) # Survey index at age
  SrvIAL = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets)) # Survey index at length
  PredSrvIdx = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_srv_fleets)) # Predicted survey index
  srv_sel = array(data = 0, dim = c(n_pop, n_regions, n_yrs + n_proj_yrs_devs, n_seas, n_ages, n_sexes, n_srv_fleets)) # Survey selectivity ages
  srv_sel_l = array(data = 0, dim = c(n_regions, n_yrs + n_proj_yrs_devs, n_lens, n_sexes, n_srv_fleets)) # Survey selectivity lengths
  srv_q = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Survey catchability

  # Likelihoods (Not population-specific)
  Catch_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Retained Fishery Catch Likelihoods
  Discard_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Discarded Fishery Likelihoods
  FishIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_fish_fleets)) # Fishery Index Likelihoods
  FishAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Retained Fishery Age Comps Likelihoods
  FishLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Retained Fishery Length Comps Likelihoods
  FishAgeComps_discard_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Discarded Fishery Age Comps Likelihoods
  FishLenComps_discard_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Discarded Fishery Length Comps Likelihoods
  SrvIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_srv_fleets)) # Survey Index Likelihoods
  SrvAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)) # Survey Age Comps Likelihoods
  SrvLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods
  conv_fish_tag_nLL = array(data = 0, dim = c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_regions, n_fish_fleets)) # Tagging Likelihoods

  # Likelihoods (population-specific)
  Catch_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Pop-specific Catch Likelihoods
  Discard_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Pop-specific Discarded Fishery Likelihoods
  FishIdx_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Pop-specific Fishery Index Likelihoods
  FishAgeComps_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Pop-specific Retained Fishery Age Comps Likelihoods
  FishLenComps_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Pop-specific Retained Fishery Length Comps Likelihoods
  FishAgeComps_discard_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Pop-specific Discarded Fishery Age Comps Likelihoods
  FishLenComps_discard_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)) # Pop-specific Discarded Fishery Length Comps Likelihoods
  SrvIdx_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_srv_fleets)) # Pop-specific Survey Index Likelihoods
  SrvAgeComps_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)) # Pop-specific Survey Age Comps Likelihoods
  SrvLenComps_pop_nLL = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)) # Pop-specific Survey Length Comps Likelihoods

  # Penalties and Priors
  Fmort_nLL = array(0, dim = dim(ln_F_devs)) # Fishing Mortality Deviation penalty
  dmr_nLL = array(0, dim = dim(logit_dmr_devs)) # Discard Mortality Deviation penalty
  Rec_nLL = array(0, dim = dim(ln_RecDevs)) # Recruitment penalty
  Init_Rec_nLL = array(0, dim = dim(ln_InitDevs)) # Initial Recruitment penalty
  sel_nLL = 0 # Penalty for selectivity deviations
  fish_q_nLL = 0 # Prior/penalty for fishery q
  srv_q_nLL = 0 # Prior/penalty for survey q
  M_nLL = 0 # Penalty/Prior for natural mortality
  h_nLL = 0 # Prior for steepness
  Movement_nLL = 0 # Penalty for movement rates
  TagRep_nLL = 0 # penalty for tag reporting rate
  rec_prop_nLL = # penalty / prior for recruitment proportions
  jnLL = 0 # Joint negative log likelihood

  # Model Process Equations -------------------------------------------------
  ## Movement Parameters (Set up) --------------------------------------------
  out_move = Get_Movement(
    move_type = move_type, # movement type (unstructured markov, or continuous time markov chain)
    do_recruits_move = do_recruits_move,

    # Dimensions
    n_pop = n_pop,
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
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(y in 1:n_yrs) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              idx = M_blocks[p,r,y,a,s] # Get indexing from blocking structure on base natural mortality
              natmort[p,r,y,a,s] = exp(ln_M[idx]) # input into natural mortality array
            } # end s loop
          } # end a loop
        } # end y loop
      } # end r loop
    } # end p loop
  } # if not using fixed natural mortality

  if(use_fixed_natmort == 1) natmort = Fixed_natmort # Using fixed natural mortality

  ## Total Fishery Selectivity -----------------------------------------------------
  if(srv_selex_type == 0) srv_selex_bins = ages # if age-based selectivity
  if(srv_selex_type == 1) srv_selex_bins = lens # if length-based selectivity
  if(fish_selex_type == 0) fish_selex_bins = ages # if age-based selectivity
  if(fish_selex_type == 1) fish_selex_bins = lens # if length-based selectivity
  if(ret_selex_type == 0) ret_selex_bins = ages # if age-based selectivity
  if(ret_selex_type == 1) ret_selex_bins = lens # if length-based selectivity

  for(r in 1:n_regions) {
    for(y in 1:(n_yrs + n_proj_yrs_devs)) {
      for(f in 1:n_fish_fleets) {

        # estimating fishery selex
        if(use_fixed_fish_sel[f] == 0) {
          for(s in 1:n_sexes) {

            # Extract variables
            if(y <= n_yrs) { # non-projection years
              fish_sel_blk_idx = fish_sel_blocks[r,y,f] # selectivity block indices
              tmp_fish_sel_model = fish_sel_model[r,y,f] # fishery selectivity model
              tmp_n_bin_nodes_bicubic = fish_sel_bicubic_binnodes[r,y,f] # true bin-node count for this block (Selex_Model == 8 only)
              tmp_n_yr_nodes_bicubic = fish_sel_bicubic_yrnodes[r,y,f] # true year-node count for this block (Selex_Model == 8 only)
            } else {
              fish_sel_blk_idx = fish_sel_blocks[r,n_yrs,f] # selectivity block indices
              tmp_fish_sel_model = fish_sel_model[r,n_yrs,f] # fishery selectivity model
              tmp_n_bin_nodes_bicubic = fish_sel_bicubic_binnodes[r,n_yrs,f] # true bin-node count for this block (Selex_Model == 8 only)
              tmp_n_yr_nodes_bicubic = fish_sel_bicubic_yrnodes[r,n_yrs,f] # true year-node count for this block (Selex_Model == 8 only)
            }

            # Extract out fixed-effect selectivity parameters for a given block
            tmp_fish_sel_vec = fish_fixed_sel_pars[r,,fish_sel_blk_idx,s,f]

            # Extract bicubic spline interpolation weight matrices for this block.
            tmp_Wbin_bicubic = array(fish_sel_bicubic_Wbin[r,,,fish_sel_blk_idx,f], dim = dim(fish_sel_bicubic_Wbin)[c(2,3)])
            tmp_Wyr_bicubic = array(fish_sel_bicubic_Wyr[r,,,fish_sel_blk_idx,f], dim = dim(fish_sel_bicubic_Wyr)[c(2,3)])

            # Compute selectivity functional form
            tmp_sel = Get_Selex(Selex_Model = tmp_fish_sel_model, # selectivity model
                                TimeVary_Model = cont_tv_fish_sel[r,f], # time varying model
                                pars = tmp_fish_sel_vec, # fixed effect selectivity parameters
                                ln_seldevs = ln_fishsel_devs[,,,,f, drop = FALSE], # Selectivity deviations
                                Region = r, # region index
                                Year = y, # year index
                                Bin = fish_selex_bins, # bin vector
                                Sex = s, # sex index
                                Wbin_bicubic = tmp_Wbin_bicubic, # bicubic spline bin-node weight matrix (Selex_Model == 8 only)
                                Wyr_bicubic = tmp_Wyr_bicubic, # bicubic spline year-node weight matrix (Selex_Model == 8 only)
                                n_bin_nodes_bicubic = tmp_n_bin_nodes_bicubic, # true bin-node count for this block (Selex_Model == 8 only)
                                n_yr_nodes_bicubic = tmp_n_yr_nodes_bicubic # true year-node count for this block (Selex_Model == 8 only)
            )

            # Compute selectivity
            for(p in 1:n_pop) {
              for(seas in 1:n_seas) {
                if(fish_selex_type == 0) fish_sel[p,r,y,seas,,s,f] = tmp_sel # age-based selectivity
              } # end seas loop
            } # end p loop
            if(fish_selex_type == 1) fish_sel_l[r,y,,s,f] = tmp_sel # input into length-based fishery selectivity

          } # end s loop
        } else {

          # Input fixed selectivity
          for(p in 1:n_pop) {
            for(seas in 1:n_seas) {
              if(fish_selex_type == 0) fish_sel[p,r,y,seas,,,f] = fish_sel_input[p,r,y,seas,,,f] # age-based selectivity
            } # end seas loop
          } # end p loop
          if(fish_selex_type == 1) fish_sel_l[r,y,,,f] = fish_sel_input[r,y,seas,,,f] # input into length-based fishery selectivity
        } # end if else for whether or not using fixed selex inputs

      } # end f loop
    } # end y loop
  } # end r loop

  for(r in 1:n_regions) {
    for(y in 1:(n_yrs + n_proj_yrs_devs)) {
      for(f in 1:n_fish_fleets) {

        if(use_fixed_ret_sel[f] == 0) {
          for(s in 1:n_sexes) {

            # Extract variables
            if(y <= n_yrs) { # non-projection years
              ret_sel_blk_idx = ret_sel_blocks[r,y,f] # selectivity block indices
              tmp_ret_sel_model = ret_sel_model[r,y,f] # fishery selectivity model
            } else {
              ret_sel_blk_idx = ret_sel_blocks[r,n_yrs,f] # selectivity block indices
              tmp_ret_sel_model = ret_sel_model[r,n_yrs,f] # fishery selectivity model
            }

            # Extract out fixed-effect selectivity parameters for a given block
            tmp_ret_sel_vec = ret_fixed_sel_pars[r,,ret_sel_blk_idx,s,f]

            # Extract bicubic spline interpolation weight matrices for this block
            tmp_Wbin_bicubic = array(ret_sel_bicubic_Wbin[r,,,ret_sel_blk_idx,f], dim = dim(ret_sel_bicubic_Wbin)[c(2,3)])
            tmp_Wyr_bicubic = array(ret_sel_bicubic_Wyr[r,,,ret_sel_blk_idx,f], dim = dim(ret_sel_bicubic_Wyr)[c(2,3)])

            # Compute selectivity functional form
            tmp_sel = Get_Selex(Selex_Model = tmp_ret_sel_model, # selectivity model
                                TimeVary_Model = cont_tv_ret_sel[r,f], # time varying model
                                pars = tmp_ret_sel_vec, # fixed effect selectivity parameters
                                ln_seldevs = ln_retsel_devs[,,,,f, drop = FALSE], # Selectivity deviations
                                Region = r, # region index
                                Year = y, # year index
                                Bin = ret_selex_bins, # bin vector
                                Sex = s, # sex index
                                Wbin_bicubic = tmp_Wbin_bicubic, # bicubic spline bin-node weight matrix (Selex_Model == 8 only)
                                Wyr_bicubic = tmp_Wyr_bicubic # bicubic spline year-node weight matrix (Selex_Model == 8 only)
            )

            # Compute selectivity
            for(p in 1:n_pop) {
              for(seas in 1:n_seas) {
                if(ret_selex_type == 0) ret_sel[p,r,y,seas,,s,f] = tmp_sel # age-based selectivity
              } # end seas loop
            } # end p loop
            if(ret_selex_type == 1) ret_sel_l[r,y,,s,f] = tmp_sel # input into length-based fishery selectivity

          } # end s loop
        } else {

          # Input fixed retention selectivity
          for(p in 1:n_pop) {
            for(seas in 1:n_seas) {
              if(ret_selex_type == 0) ret_sel[p,r,y,seas,,,f] = ret_sel_input[p,r,y,seas,,,f] # age-based selectivity
            } # end seas loop
          } # end p loop
          if(ret_selex_type == 1) ret_sel_l[r,y,,,f] = ret_sel_input[r,y,seas,,,f] # input into length-based fishery selectivity
        } # end if else for fixed retention selectivity

      } # end f loop
    } # end y loop
  } # end r loop

  # Mean standardization for selectivity
  for(r in 1:n_regions) {
    for(f in 1:n_fish_fleets) {

      # Mean Standardizing to help with interpretability (fishery total)
      if(fish_selex_type == 0) if(cont_tv_fish_sel[r,f] %in% 3:5 || any(fish_sel_model[r,,f] == 5)) for(s in 1:n_sexes) {
        tmp_mean = log(mean(fish_sel[1,r,,1,,s,f])) # indexing 1 for pop and season because those dims not estimated (depends on growth transition)
        for(p in 1:n_pop) for(seas in 1:n_seas)
          fish_sel[p,r,,seas,,s,f] = exp(log(fish_sel[p,r,,seas,,s,f]) - tmp_mean)
      }
      if(fish_selex_type == 1) if(cont_tv_fish_sel[r,f] %in% 3:5 || any(fish_sel_model[r,,f] == 5)) for(s in 1:n_sexes) fish_sel_l[r,,,s,f] = exp(log(fish_sel_l[r,,,s,f]) - log(mean(fish_sel_l[r,,,s,f]))) # length-based selectivity


      # Mean Standardizing to help with interpretability (fishery retention)
      if(ret_selex_type == 0) if(cont_tv_ret_sel[r,f] %in% 3:5 || any(ret_sel_model[r,,f] == 5)) for(s in 1:n_sexes) {
        tmp_mean = log(mean(ret_sel[1,r,,1,,s,f])) # indexing 1 for pop and season because those dims not estimated (depends on growth transition)
        for(p in 1:n_pop) for(seas in 1:n_seas)
          ret_sel[p,r,,seas,,s,f] = exp(log(ret_sel[p,r,,seas,,s,f]) - tmp_mean)
      }
      if(ret_selex_type == 1) if(cont_tv_ret_sel[r,f] %in% 3:5 || any(ret_sel_model[r,,f] == 5)) for(s in 1:n_sexes) ret_sel_l[r,,,s,f] = exp(log(ret_sel_l[r,,,s,f]) - log(mean(ret_sel_l[r,,,s,f]))) # length-based selectivity

    } # end f loop
  } # end r loop

  ## Survey Selectivity ------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:(n_yrs + n_proj_yrs_devs)) {
      for(sf in 1:n_srv_fleets) {

        if(use_fixed_srv_sel[sf] == 0) {
          for(s in 1:n_sexes) {

            # Extract variables
            if(y <= n_yrs) { # non-projection years
              srv_sel_blk_idx = srv_sel_blocks[r,y,sf] # selectivity block indices
              tmp_srv_sel_model = srv_sel_model[r,y,sf] # survey selectivity model
              tmp_n_bin_nodes_bicubic = srv_sel_bicubic_binnodes[r,y,sf] # true bin-node count for this block (Selex_Model == 8 only)
              tmp_n_yr_nodes_bicubic = srv_sel_bicubic_yrnodes[r,y,sf] # true year-node count for this block (Selex_Model == 8 only)
            } else {
              srv_sel_blk_idx = srv_sel_blocks[r,n_yrs,sf] # selectivity block indices
              tmp_srv_sel_model = srv_sel_model[r,n_yrs,sf] # survey selectivity model
              tmp_n_bin_nodes_bicubic = srv_sel_bicubic_binnodes[r,n_yrs,sf] # true bin-node count for this block (Selex_Model == 8 only)
              tmp_n_yr_nodes_bicubic = srv_sel_bicubic_yrnodes[r,n_yrs,sf] # true year-node count for this block (Selex_Model == 8 only)
            }

            # Extract out fixed-effect selectivity parameters for a given block
            tmp_srv_sel_vec = srv_fixed_sel_pars[r,,srv_sel_blk_idx,s,sf]

            # Extract bicubic spline interpolation weight matrices for this block.
            tmp_Wbin_bicubic = array(srv_sel_bicubic_Wbin[r,,,srv_sel_blk_idx,sf], dim = dim(srv_sel_bicubic_Wbin)[c(2,3)])
            tmp_Wyr_bicubic = array(srv_sel_bicubic_Wyr[r,,,srv_sel_blk_idx,sf], dim = dim(srv_sel_bicubic_Wyr)[c(2,3)])

            # Compute selectivity functional form
            tmp_sel = Get_Selex(Selex_Model = tmp_srv_sel_model, # selectivity model
                                TimeVary_Model = cont_tv_srv_sel[r,sf], # time varying model
                                pars = tmp_srv_sel_vec, # fixed effect selectivity parameters
                                ln_seldevs = ln_srvsel_devs[,,,,sf, drop = FALSE], # Selectivity deviations
                                Region = r, # region index
                                Year = y, # year index
                                Bin = srv_selex_bins, # bin vector
                                Sex = s, # sex index
                                Wbin_bicubic = tmp_Wbin_bicubic, # bicubic spline bin-node weight matrix (Selex_Model == 8 only)
                                Wyr_bicubic = tmp_Wyr_bicubic, # bicubic spline year-node weight matrix (Selex_Model == 8 only)
                                n_bin_nodes_bicubic = tmp_n_bin_nodes_bicubic, # true bin-node count for this block (Selex_Model == 8 only)
                                n_yr_nodes_bicubic = tmp_n_yr_nodes_bicubic # true year-node count for this block (Selex_Model == 8 only)
            )

            # Calculate selectivity
            for(p in 1:n_pop) {
              for(seas in 1:n_seas) {
                if(srv_selex_type == 0) srv_sel[p,r,y,seas,,s,sf] = tmp_sel # age-based selectivity
              } # end seas loop
            } # end p loop
            if(srv_selex_type == 1) srv_sel_l[r,y,,s,sf] = tmp_sel # input into length-based fishery selectivity

          } # end s loop
        } else {
          # Input fixed survey selectivity
          for(p in 1:n_pop) {
            for(seas in 1:n_seas) {
              if(srv_selex_type == 0) srv_sel[p,r,y,seas,,,sf] = srv_sel_input[p,r,y,seas,,,sf] # age-based selectivity
            } # end seas loop
          } # end p loop
          if(srv_selex_type == 1) srv_sel_l[r,y,,,sf] = srv_sel_input[r,y,seas,,,sf] # input into length-based fishery selectivity
        } # end if for whether to estiamte or fix survey selex

      } # end sf loop
    } # end y loop
  } # end r loop

  # Mean standardization for survey selectivity
  for(r in 1:n_regions) {
    for(sf in 1:n_srv_fleets) {

      # Mean Standardizing to help with interpretability (survey)
      if(srv_selex_type == 0) if(cont_tv_srv_sel[r,sf] %in% 3:5 || any(srv_sel_model[r,,sf] == 5)) for(s in 1:n_sexes) {
        tmp_mean = log(mean(srv_sel[1,r,,1,,s,sf])) # indexing 1 for pop and season because those dims not estimated (depends on growth transition)
        for(p in 1:n_pop) for(seas in 1:n_seas)
          srv_sel[p,r,,seas,,s,sf] = exp(log(srv_sel[p,r,,seas,,s,sf]) - tmp_mean)
      }
      if(srv_selex_type == 1) if(cont_tv_srv_sel[r,sf] %in% 3:5 || any(srv_sel_model[r,,sf] == 5)) for(s in 1:n_sexes) srv_sel_l[r,,,s,sf] = exp(log(srv_sel_l[r,,,s,sf]) - log(mean(srv_sel_l[r,,,s,sf])))

    } # end sf
  } # end r

  ## Mortality ---------------------------------------------------------------
  missing_catch = is.na(ObsCatch) # TRUE = aggregate catch observation is missing (not a true recorded zero)

  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(seas in 1:n_seas) {
        for(f in 1:n_fish_fleets) {

          # A cell is a true closure only when no catch is fit
          is_closed = (UseCatch[r,y,seas,f] == 0) && all(UseCatch_pop[,r,y,seas,f] == 0) && !missing_catch[r,y,seas,f]

          if(is_closed) {
            Fmort[r,y,seas,f] = 0
            dmr[r,y,seas,f] = 0
          } else {
            Fmort[r,y,seas,f] = exp(ln_F_mean[r,seas,f] + ln_F_devs[r,y,seas,f])
            dmr[r,y,seas,f] = RTMB::plogis(logit_dmr_mean[r,seas,f] + logit_dmr_devs[r,y,seas,f])
          }

          # get fishing mortality at age
          for(p in 1:n_pop) {
            if(fish_selex_type == 1) for(s in 1:n_sexes) fish_sel[p,r,y,seas,,s,f] = fish_sel_l[r,y,,s,f] %*% SizeAgeTrans[p,r,y,seas,,,s]
            if(ret_selex_type == 1) for(s in 1:n_sexes) ret_sel[p,r,y,seas,,s,f] = ret_sel_l[r,y,,s,f] %*% SizeAgeTrans[p,r,y,seas,,,s]
            ret_FAA[p,r,y,seas,,,f] = Fmort[r,y,seas,f] * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f] # Retained fishing mortality at age
            disc_FAA[p,r,y,seas,,,f] = Fmort[r,y,seas,f] * fish_sel[p,r,y,seas,,,f] * (1 - ret_sel[p,r,y,seas,,,f]) * dmr[r,y,seas,f] # Discarded fishing mortality at age
            tot_FAA[p,r,y,seas,,,f] = ret_FAA[p,r,y,seas,,,f] +  disc_FAA[p,r,y,seas,,,f]# Total fishing mortality at age

          }

        } # f loop

        # get total mortality
        for(p1 in 1:n_pop) for(a in 1:n_ages) for(s in 1:n_sexes)
          ZAA[p1,r,y,seas,a,s] = sum(ret_FAA[p1,r,y,seas,a,s,]) + sum(disc_FAA[p1,r,y,seas,a,s,]) + (natmort[p1,r,y,a,s] * seasdur[seas])
      } # end seas loop
    } # end y loop
  } # end r loop


  ## Recruitment Transformations and Bias Ramp (Methot and Taylor) -------------------------------
  ### Parameter Transformations -----------------------------------------------
  # Mean or virgin recruitment area proportions
  if(n_regions > 1) {
    if(n_pop == 1) {  # if spatial model, with recruitment dispersal
      tmp_rec_region_prop = c(0, rec_region_prop_pars[1,]) # set up vector for transformation
      rec_region_prop[1,] = exp(tmp_rec_region_prop) / sum(exp(tmp_rec_region_prop)) # do multinomial logit to get recruitment regional proportions
    } else {
      for(p in 1:n_pop) {
        if(rec_region_prop_spec == 0) { # Recruitment dispersal with natal homing
          tmp_rec_region_prop = c(0, rec_region_prop_pars[p,]) # set up vector for transformation
          rec_region_prop[p,] = exp(tmp_rec_region_prop) / sum(exp(tmp_rec_region_prop)) # do multinomial logit to get recruitment regional proportions
        }
        # No recruitment dispersal
        if(rec_region_prop_spec == 1) rec_region_prop[p,natal_region[p]] = 1
      } # end p loop
    }
  } else rec_region_prop[] = 1 # non-spatial model

  # Mean or virgin recrutment seasonal proportions
  if(use_fixed_rec_seas_prop == 1) { # use input fixed proportions
    rec_seas_prop[] = fixed_rec_seas_prop
  } else if(n_seas > 1) {
    for(p in 1:n_pop) {
      if(rec_lag == 0 && spawn_seas > 1) {
        # Age-0 (rec_lag = 0) recruitment seasonal proportion
        n_allowed = n_seas - spawn_seas + 1
        tmp_rec_seas_prop_pars = c(0, rec_seas_prop_pars[p, 1:(n_allowed - 1)])
        allowed_prop = exp(tmp_rec_seas_prop_pars) / sum(exp(tmp_rec_seas_prop_pars))
        rec_seas_prop[p,] = 0
        rec_seas_prop[p, spawn_seas:n_seas] = allowed_prop
      } else {
        tmp_rec_seas_prop_pars = c(0, rec_seas_prop_pars[p,]) # set up vector for transformation
        rec_seas_prop[p,] = exp(tmp_rec_seas_prop_pars) / sum(exp(tmp_rec_seas_prop_pars)) # do multinomial logit to get recruitment area proportions
      }
    } # end p loop
  } else rec_seas_prop[] = 1 # non-seasonal model

  # Global recruitment
  R0_r = array(0, dim = c(n_pop, n_regions)) # container
  R0 = exp(ln_global_R0) # exponentiate
  for(p in 1:n_pop) R0_r[p,] = R0[p] * rec_region_prop[p,]

  # Global rinit
  rinit_r = array(0, dim = c(n_pop, n_regions)) # container
  rinit = exp(ln_rinit) # exponentiate
  for(p in 1:n_pop) rinit_r[p,] = rinit[p] * rec_region_prop[p,]

  # Steepness
  h_trans = array(0, dim = c(n_pop, n_regions))
  for(p in 1:n_pop) for(r in 1:n_regions) h_trans[p,r] = 0.2 + (1 - 0.2) * RTMB::plogis(steepness_h[p,r]) # bound steepness between 0.2 and 1

  # Recruitment SD
  sigmaR2_early = array(exp(ln_sigmaR[1,,])^2, dim = c(n_pop, n_regions)) # recruitment variability for early period
  sigmaR2_late = array(exp(ln_sigmaR[2,,])^2, dim = c(n_pop, n_regions)) # recruitment variability for late period

  # Recruitment sex-ratio
  if(n_sexes == 2) {
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        for(y in 1:n_yrs) {
          sexratio_blk_idx = sexratio_blocks[p,r,y] # extract sex ratio block
          sexratio_f = RTMB::plogis(sexratio_pars[p,r,sexratio_blk_idx]) # get female recruitment sex-ratio
          sexratio[p,r,y,] = c(sexratio_f, 1 - sexratio_f) # input total sex ratio
        } # end y loop
      } # end r loop
    } # end p loop
  } else sexratio[] = 1 # set recruitment sex ratio at 1

  # Stray rates
  if(use_fixed_stray_rate == 0) {
    for(p in 1:n_pop) {
      for(y in 1:n_yrs) {
        idx              <- stray_rate_blocks[p, y] # get idx
        stray_rate[p, y] <- RTMB::plogis(stray_rate_pars[p, idx])  # p dimension explicit
      } # end y loop
    } # end p loop
  } # if not using fixed rates

  if(use_fixed_stray_rate == 1) stray_rate = fixed_stray_rate # Using fixed stray rates

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
  # Get initial F for age structure
  catch_flag_base = array(UseCatch[,1,,], dim = c(n_regions, n_seas, n_fish_fleets))
  catch_flag_pop = apply(UseCatch_pop[,,1,,,drop = FALSE], c(2,4,5), max)
  catch_flag = pmax(catch_flag_base, catch_flag_pop)
  init_F = init_F_prop * exp(ln_F_mean) * catch_flag

  # Get initial fished NAA
  Init_Fished_NAA = Get_Init_NAA(
    init_age_strc = init_age_strc, # initial age structure
    init_iter = n_ages * 5, # if init_age_strc == 0, number of iterations to run
    n_pop = n_pop, # populations
    n_regions = n_regions, # regions
    n_sexes = n_sexes, # sexes
    n_ages = n_ages, # ages
    n_seas = n_seas, # seasons
    n_fish_fleets = n_fish_fleets, # fleets
    seasdur = seasdur, # seasonal duration
    rec_seas_prop = rec_seas_prop,
    natmort = array(natmort[,,1,,], dim = c(n_pop, n_regions, n_ages, n_sexes)), # natural mortality in first year
    init_F = init_F, # initial F applied
    fish_sel = array(fish_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # total fishery selectivity in first year
    R0_r = if(use_rinit == 0) R0_r else rinit_r, # regional mean or virgin recruitment
    sexratio = array(sexratio[,,1,], dim = c(n_pop, n_regions, n_sexes)), # sex ratio in first year
    Movement = array(Movement[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
    do_recruits_move = do_recruits_move, # whether recruits move
    ln_InitDevs = ln_InitDevs, # initial deviations
    ret_sel = array(ret_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # retained fishery selectivity in first year
    dmr = array(dmr[,1,,], dim = c(n_regions, n_seas, n_fish_fleets))
  )

  # Get initial unfished NAA
  Init_Unfished_NAA = Get_Init_NAA(
    init_age_strc = init_age_strc, # initial age structure
    init_iter = n_ages * 5, # if init_age_strc == 0, number of iterations to run
    n_pop = n_pop, # populations
    n_regions = n_regions, # regions
    n_sexes = n_sexes, # sexes
    n_ages = n_ages, # ages
    n_fish_fleets = n_fish_fleets, # fleets
    n_seas = n_seas, # seasons
    seasdur = seasdur, # seasonal duration
    rec_seas_prop = rec_seas_prop,
    natmort = array(natmort[,,1,,], dim = c(n_pop, n_regions, n_ages, n_sexes)), # natural mortality in first year
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)), # initial F applied (0 for unfished)
    fish_sel = array(fish_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # total fishery selectivity in first year
    R0_r = if(use_rinit == 0) R0_r else rinit_r, # regional mean or virgin recruitment
    sexratio = array(sexratio[,,1,], dim = c(n_pop, n_regions, n_sexes)), # sex ratio in first year
    Movement = array(Movement[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
    do_recruits_move = do_recruits_move, # whether recruits move
    ln_InitDevs = ln_InitDevs, # initial deviations
    ret_sel = array(ret_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # retained fishery selectivity in first year
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)) # unfished
  )

  # Input into model arrays (first year and season) - and add lognormal mean adjustment
  NAA[,,1,1,,] = Init_Fished_NAA
  NAA0[,,1,1,,] = Init_Unfished_NAA

  ## Population Projection ---------------------------------------------------
  for(y in 1:n_yrs) {

    ### Annual Recruitment (rec_lag != 0 only) -----------------------------------
    if(rec_lag != 0) {

      # Get Deterministic Recruitment
      tmp_Det_Rec = Get_Det_Recruitment(recruitment_model = rec_model,
                                        rec_dd = rec_dd,
                                        R0 = R0,
                                        rec_region_prop = rec_region_prop,
                                        rec_seas_prop = rec_seas_prop,
                                        h = h_trans,
                                        n_pop = n_pop,
                                        n_ages = n_ages,
                                        n_regions = n_regions,
                                        # Note: Using first year and female quantities to compute unfished SSB0
                                        sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,1,1], dim = c(n_pop, n_regions)),
                                        WAA = array(WAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                        MatAA = array(MatAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                        natmort = array(natmort[,,1,,1], dim = c(n_pop, n_regions, n_ages)),
                                        Movement = array(Movement[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                        stray_rate = array(stray_rate[,1], dim = c(n_pop)),
                                        sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,1,,1], dim = c(n_pop, n_regions, n_regions, n_ages)),
                                        do_recruits_move = do_recruits_move,
                                        natal_region = natal_region,
                                        t_spawn = t_spawn,
                                        SSB_vals = SSB,
                                        y = y,
                                        n_seas = n_seas,
                                        spawn_seas = spawn_seas,
                                        seasdur = seasdur,
                                        rec_lag = rec_lag,
                                        n_fish_fleets = n_fish_fleets,
                                        init_F = init_F, # initF
                                        fish_sel = array(fish_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # total fishery selectivity
                                        ret_sel = array(ret_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # retained fishery selectivity in first year
                                        dmr = array(dmr[,1,,], dim = c(n_regions, n_seas, n_fish_fleets))
      )

      for(p in 1:n_pop) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            if(y <= n_est_rec_devs) tmp_total_rec = tmp_Det_Rec[p,r] * exp(ln_RecDevs[p,r,y])
            if(y > n_est_rec_devs) tmp_total_rec = tmp_Det_Rec[p,r]
            # season 1 fraction
            NAA[p,r,y,1,1,s] = tmp_total_rec * rec_seas_prop[p,1] * sexratio[p,r,y,s]
          }
          Rec[p,r,y] = tmp_total_rec  # store total before seasonal split
          NAA0[p,r,y,1,1,] = NAA[p,r,y,1,1,]
        } # end r loop
      } # end p loop
    } # end if rec_lag != 0

    for(seas in 1:n_seas) {

      # Insert seasonal recruits
      if(if(rec_lag != 0) seas > 1 else seas > spawn_seas) {
        for(p in 1:n_pop) {
          for(r in 1:n_regions) {
            for(s in 1:n_sexes) {
              NAA[p,r,y,seas,1,s]  = NAA[p,r,y,seas,1,s]  + Rec[p,r,y] * rec_seas_prop[p,seas] * sexratio[p,r,y,s]
              NAA0[p,r,y,seas,1,s] = NAA0[p,r,y,seas,1,s] + Rec[p,r,y] * rec_seas_prop[p,seas] * sexratio[p,r,y,s]
            } # end s loop
          } # end r loop
        } # end p loop
      }

      ### Movement ----------------------------------------------------------------
      # Record values prior to movement
      NAA_bef[,,y,seas,,] = NAA[,,y,seas,,]

      if(n_regions > 1) {
        for(p in 1:n_pop) {
          # Recruits don't move
          if(do_recruits_move == 0) {
            # Apply movement after ageing processes - start movement at age 2
            for(a in 2:n_ages) {
              for(s in 1:n_sexes) {
                NAA[p,,y,seas,a,s] = t(NAA[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # Fished
              } # end s loop
            } # end a loop
          } # end if recruits don't move

          # Recruits move here
          if(do_recruits_move == 1) {
            for(a in 1:n_ages) {
              for(s in 1:n_sexes) {
                NAA[p,,y,seas,a,s] = t(NAA[p,,y,seas,a,s]) %*% Movement[p,,,y,seas,a,s] # Fished
              } # end s loop
            } # end a loop
          } # end if
        } # end p loop

        # Record values after movement
        NAA_aft[,,y,seas,,] = NAA[,,y,seas,,]

      } # only compute if spatial

      ### Compute Biomass Quantities + Recruitment (rec_lag == 0 only) ------------
      if(rec_lag == 0 && seas == spawn_seas) {

        # SSB from survivors only
        biom = compute_biom_y(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, seasdur,
                              n_seas, n_pop, n_regions, n_ages, n_sexes,
                              sgl_seas_spawning_movement, natal_region, stray_rate)
        SSB[,, y] = biom$SSB_y

        tmp_Det_Rec = Get_Det_Recruitment(recruitment_model = rec_model,
                                          rec_dd = rec_dd,
                                          R0 = R0,
                                          rec_region_prop = rec_region_prop,
                                          rec_seas_prop = rec_seas_prop,
                                          h = h_trans,
                                          n_pop = n_pop,
                                          n_ages = n_ages,
                                          n_regions = n_regions,
                                          sexratio_f = if(n_sexes == 1) array(0.5, dim = c(n_pop, n_regions)) else array(sexratio[,,1,1], dim = c(n_pop, n_regions)),
                                          WAA = array(WAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                          MatAA = array(MatAA[,,1,,,1], dim = c(n_pop, n_regions, n_seas, n_ages)),
                                          natmort = array(natmort[,,1,,1], dim = c(n_pop, n_regions, n_ages)),
                                          Movement = array(Movement[,,,1,,,1], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages)),
                                          stray_rate = array(stray_rate[,1], dim = c(n_pop)),
                                          sgl_seas_spawning_movement = array(sgl_seas_spawning_movement[,,,1,,1], dim = c(n_pop, n_regions, n_regions, n_ages)),
                                          do_recruits_move = do_recruits_move,
                                          natal_region = natal_region,
                                          t_spawn = t_spawn,
                                          SSB_vals = SSB,
                                          y = y,
                                          n_seas = n_seas,
                                          spawn_seas = spawn_seas,
                                          seasdur = seasdur,
                                          rec_lag = rec_lag,
                                          n_fish_fleets = n_fish_fleets,
                                          init_F = init_F, # initF
                                          fish_sel = array(fish_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # total fishery selectivity
                                          ret_sel = array(ret_sel[,,1,,,1,], dim = c(n_pop, n_regions, n_seas, n_ages, n_fish_fleets)), # retained fishery selectivity in first year
                                          dmr = array(dmr[,1,,], dim = c(n_regions, n_seas, n_fish_fleets))
        )

        for(p in 1:n_pop) {
          for(r in 1:n_regions) {
            for(s in 1:n_sexes) {
              if(y <= n_est_rec_devs) tmp_total_rec = tmp_Det_Rec[p,r] * exp(ln_RecDevs[p,r,y])
              if(y > n_est_rec_devs) tmp_total_rec = tmp_Det_Rec[p,r]
              NAA[p,r,y,spawn_seas,1,s]  = tmp_total_rec * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,s]
              NAA0[p,r,y,spawn_seas,1,s] = tmp_total_rec * rec_seas_prop[p,spawn_seas] * sexratio[p,r,y,s]
            }
            Rec[p,r,y] = tmp_total_rec # store total before seasonal split
          } # end r loop
        } # end p loop

        # Recruits just inserted above missed this season's movement step. Move recruits if allowed.
        if(do_recruits_move == 1 && n_regions > 1) {
          for(p in 1:n_pop) {
            for(s in 1:n_sexes) {
              NAA[p,,y,seas,1,s] = t(NAA[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
              NAA0[p,,y,seas,1,s] = t(NAA0[p,,y,seas,1,s]) %*% Movement[p,,,y,seas,1,s]
            } # end s loop
          } # end p loop
          NAA_aft[,,y,seas,1,] = NAA[,,y,seas,1,]
        }

        # Recompute biomass quantities now that this year's recruits are included
        biom = compute_biom_y(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, seasdur,
                              n_seas, n_pop, n_regions, n_ages, n_sexes,
                              sgl_seas_spawning_movement, natal_region, stray_rate)
        Total_Biom[,, y] = biom$Total_Biom_y
        SSB[,, y] = biom$SSB_y
        Dynamic_SSB0[,,y] = biom$Dynamic_SSB0_y
        eff_SSB[,y] = biom$eff_SSB_y

      } # end if rec_lag == 0 && seas == spawn_seas

      ### Mortality and Ageing ------------------------------------------------------
      if(seas < n_seas) {
        # within year / seasonal mortality
        NAA[,,y,seas+1,1:n_ages,] = NAA[,,y,seas,1:n_ages,] * exp(-ZAA[,,y,seas,1:n_ages,])
        NAA0[,,y,seas+1,1:n_ages,] = NAA0[,,y,seas,1:n_ages,] * exp(-(natmort[,,y,1:n_ages,] * seasdur[seas]))
      } else {
        # age advancement and enter into first season of next year
        # Fished
        NAA[,,y+1,1,2:n_ages,] = NAA[,,y,n_seas,1:(n_ages-1),] * exp(-ZAA[,,y,seas,1:(n_ages-1),]) # Exponential mortality for individuals not in plus group
        NAA[,,y+1,1,n_ages,] = NAA[,,y+1,1,n_ages,] + NAA[,,y,seas,n_ages,] * exp(-ZAA[,,y,seas,n_ages,]) # Acuumulate plus group
        # Unfished
        NAA0[,,y+1,1,2:n_ages,] = NAA0[,,y,n_seas,1:(n_ages-1),] * exp(-natmort[,,y,1:(n_ages-1),] * seasdur[n_seas]) # Exponential mortality for individuals not in plus group
        NAA0[,,y+1,1,n_ages,] = NAA0[,,y+1,1,n_ages,] + NAA0[,,y,n_seas,n_ages,] * exp(-natmort[,,y,n_ages,] * seasdur[n_seas]) # Acuumulate plus group
      }

      ### Compute Biomass Quantities (rec_lag != 0: unchanged original timing) ----
      if(rec_lag != 0 && seas == spawn_seas) {
        spawn_biom = compute_biom_y(y, seas, NAA, NAA0, WAA, MatAA, ZAA, natmort, t_spawn, seasdur,
                                    n_seas, n_pop, n_regions, n_ages, n_sexes,
                                    sgl_seas_spawning_movement, natal_region, stray_rate)
        Total_Biom[,, y] = spawn_biom$Total_Biom_y
        SSB[,, y] = spawn_biom$SSB_y
        Dynamic_SSB0[,,y] = spawn_biom$Dynamic_SSB0_y
        eff_SSB[,y] = spawn_biom$eff_SSB_y
      }

    } # end seas loop
  } # end y loop

  # Get aggregated SSB values
  Aggregated_SSB = apply(SSB, 3, sum)
  Dynamic_Aggregated_SSB0 = apply(Dynamic_SSB0, 3, sum)

  ## Fishery Observation Model -----------------------------------------------
  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(f in 1:n_fish_fleets) {

          fish_q_blk_idx = fish_q_blocks[r,y,f] # get time-block catchability index
          fish_q[r,y,f] = exp(ln_fish_q[r,fish_q_blk_idx,f]) # Input into fishery catchability container

          for(seas in 1:n_seas) {

            # Retained Catch at Age
            CAA[p,r,y,seas,,,f] = ret_FAA[p,r,y,seas,,,f] / ZAA[p,r,y,seas,,] * NAA[p,r,y,seas,,] * (1 - exp(-ZAA[p,r,y,seas,,]))

            # Dead Discarded Catch at Age
            DAA[p,r,y,seas,,,f] = disc_FAA[p,r,y,seas,,,f] / ZAA[p,r,y,seas,,] * NAA[p,r,y,seas,,] * (1 - exp(-ZAA[p,r,y,seas,,]))

            if(fit_lengths == 1) {
              for(s in 1:n_sexes) {
                CAL[p,r,y,seas,,s,f] = SizeAgeTrans[p,r,y,seas,,,s] %*% CAA[p,r,y,seas,,s,f] # Retained Catch at length
                DAL[p,r,y,seas,,s,f] = SizeAgeTrans[p,r,y,seas,,,s] %*% DAA[p,r,y,seas,,s,f] # Discarded Catch at length
              } # end s loop
            } # fitting lengths

            # Get catch
            if(catch_units[f] == 0) PredCatch[p,r,y,seas,f] = sum(CAA[p,r,y,seas,,,f]) # abundance
            if(catch_units[f] == 1) PredCatch[p,r,y,seas,f] = sum(CAA[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f]) # biomass

            # Get discards
            if(discard_units[f] == 0) PredDiscard[p,r,y,seas,f] = sum(DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f])  # total discard abundance
            if(discard_units[f] == 1) PredDiscard[p,r,y,seas,f] = sum(DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f] * WAA_fish[p,r,y,seas,,,f])  # total discard biomass
            if(discard_units[f] == 2) {
              total_catch = CAA[p,r,y,seas,,,f] + DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f]
              PredDiscard[p,r,y,seas,f] = 1 - sum(CAA[p,r,y,seas,,,f]) / sum(total_catch)
            } # abundance fraction
            if(discard_units[f] == 3) {
              total_catch = CAA[p,r,y,seas,,,f] + DAA[p,r,y,seas,,,f] / dmr[r,y,seas,f]
              PredDiscard[p,r,y,seas,f] = 1 - sum(CAA[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f]) / sum(total_catch * WAA_fish[p,r,y,seas,,,f])
            } # biomass fraction

            # Get fishery index
            if(fish_idx_type[f] == 0) PredFishIdx[p,r,y,seas,f] = fish_q[r,y,f] * sum(NAA[p,r,y,seas,,] * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f]) # retained abundance
            if(fish_idx_type[f] == 1) PredFishIdx[p,r,y,seas,f] = fish_q[r,y,f] * sum(NAA[p,r,y,seas,,] * fish_sel[p,r,y,seas,,,f] * ret_sel[p,r,y,seas,,,f] * WAA_fish[p,r,y,seas,,,f]) # retained biomass
          } # end seas loop

        } # end f loop
      } # end y loop
    } # end r loop
  } # end p loop

  ## Survey Observation Model ------------------------------------------------
  for(p in 1:n_pop) {
    for(r in 1:n_regions) {
      for(y in 1:n_yrs) {
        for(sf in 1:n_srv_fleets) {

          srv_q_blk_idx = srv_q_blocks[r,y,sf] # get time-block catchability index
          srv_q[r,y,sf] = exp(ln_srv_q[r,srv_q_blk_idx,sf]) # Input into survey catchability container
          if(do_srv_q_cov == 1) srv_q[r,y,sf] = srv_q[r,y,sf] * exp(sum(srv_q_cov[r,y,sf,] * srv_q_coeff[r,sf,])) # adding covariate effects

          for(seas in 1:n_seas) {

            # Convert length-selex to age-selex
            if(srv_selex_type == 1) for(s in 1:n_sexes) srv_sel[p,r,y,seas,,s,sf] = srv_sel_l[r,y,,s,sf] %*% SizeAgeTrans[p,r,y,seas,,,s]

            SrvIAA[p,r,y,seas,,,sf] = NAA[p,r,y,seas,,] * srv_sel[p,r,y,seas,,,sf] * exp(-t_srv[r,seas,sf] * ZAA[p,r,y,seas,,]) # Survey index at age

            if(fit_lengths == 1) {
              for(s in 1:n_sexes) {
                SrvIAL[p,r,y,seas,,s,sf] = SizeAgeTrans[p,r,y,seas,,,s] %*% SrvIAA[p,r,y,seas,,s,sf] # Survey index at length
              } # end s loop
            } # fitting lengths

            if(srv_idx_type[sf] == 0) PredSrvIdx[p,r,y,seas,sf] = srv_q[r,y,sf] * sum(SrvIAA[p,r,y,seas,,,sf]) # abundance
            if(srv_idx_type[sf] == 1) PredSrvIdx[p,r,y,seas,sf] = srv_q[r,y,sf] * sum(SrvIAA[p,r,y,seas,,,sf] * WAA_srv[p,r,y,seas,,,sf]) # biomass

          } # end seas loop

        } # end sf loop
      } # end y loop
    } # end r loop
  } # end p loop


  ## Tagging Observation Model -----------------------------------------------
  if(any(use_conv_fish_tagging == 1)) {

    # Set up tag reporting rates
    for(f in 1:n_fish_fleets) {
      for(r in 1:n_regions) {
        conv_tagrep_blk_idx = conv_tag_fish_reporting_blocks[r,,f]  # Get all blocks for this region
        conv_tag_fish_reporting[r,,f] = RTMB::plogis(conv_tag_fish_reporting_pars[r,conv_tagrep_blk_idx,f])  # inverse logit transform
      } # end r loop
    } # end f loop

    for(tc in 1:n_conv_tag_cohorts) {

      tr = conv_tag_release_indicator[tc,1] # extract tag release region
      ty = conv_tag_release_indicator[tc,2] # extract tag release year
      tseas = conv_tag_release_indicator[tc,3] # extract tag release season

      for(ry in 1:min(conv_tag_max_liberty, n_yrs - ty + 1)) {   # years
        y = ty + ry - 1 # get real year
        for(rseas in 1:n_seas) { # seasons

          # get fishing mortality
          tmp_FAA = array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
          tmp_ret_FAA = array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
          tmp_disc_DAA = array(0, dim = c(n_pop, n_regions, 1, n_ages, n_sexes, n_fish_fleets))
          for(p in 1:n_pop) for(f in 1:n_fish_fleets) {
            if(use_conv_fish_tagging[f] == 1) {
              tmp_ret_FAA[p,,1,,,f] = Fmort[, y, rseas, f] * fish_sel[p,,y,rseas,,,f] * ret_sel[p,,y,rseas,,,f]  # Retained fishing mortality
              tmp_disc_DAA[p,,1,,,f] = Fmort[, y, rseas, f] * fish_sel[p,,y,rseas,,,f] * (1 - ret_sel[p,,y,rseas,,,f]) * dmr[,y,rseas,f] # Dead discard fishing mortality
              tmp_FAA[p,,1,,,f] = tmp_ret_FAA[p,,1,,,f] + tmp_disc_DAA[p,,1,,,f] # Total fishing mortality
            } # end if
          } # end p loop

          # get total mortality
          tmp_natmort = array(natmort[,,y,,], dim = c(n_pop, n_regions, 1, n_ages, n_sexes))
          tmp_ZAA = (tmp_natmort * seasdur[rseas]) + apply(tmp_FAA, 1:5, sum) + (exp(ln_conv_tag_shed[tc]) * seasdur[rseas])

          # Discount with tagging time (conv_tag_t_tagging) if it doesn't happen at the start of the season / year
          if(ry == 1 && rseas == tseas) {

            if(conv_tag_t_tagging[tc] != 1) tmp_ZAA = tmp_ZAA * conv_tag_t_tagging[tc]

            # apportion tagged fish out to appropriate dimensions if necessary
            tmp_tagged_fish = release_conv_tag_attr(array(conv_tagged_fish[tc, , , ], dim = c(n_pop, n_ages, n_sexes)),
                                                    conv_fish_tag_attr,
                                                    conv_tag_release_platform[tc,],
                                                    srv_sel, fish_sel, NAA_bef,
                                                    ty, tseas, tr, n_pop,
                                                    n_ages, n_sexes)

            # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
            conv_tag_fish_avail[1, rseas, tc, , tr, , ] = array(tmp_tagged_fish * exp(-exp(ln_init_conv_tag_mort[tc])), dim = c(n_pop, n_ages, n_sexes))
          }

          # get temporary survival value
          tmp_SAA = exp(-tmp_ZAA)

          # Move tagged fish around (skip only in first release year + tagging season when tagging occurs mid-season)
          if(conv_tag_t_tagging[tc] == 1 || ry != 1 || rseas != tseas) {
            for(p in 1:n_pop) {
              # Movement of tag cohorts
              if(do_recruits_move == 0) {
                for(a in 2:n_ages) for(s in 1:n_sexes) {
                  conv_tag_fish_avail[ry, rseas, tc, p, , a, s] =
                    t(conv_tag_fish_avail[ry, rseas, tc, p, , a, s]) %*%
                    Movement[p, , , y, rseas, a, s]
                }
              } else { # if recruits move
                for(a in 1:n_ages) for(s in 1:n_sexes) {
                  conv_tag_fish_avail[ry, rseas, tc, p, , a, s] =
                    t(conv_tag_fish_avail[ry, rseas, tc, p, , a, s]) %*%
                    Movement[p, , , y, rseas, a, s]
                } # end s loop
              } # end else
            } # end p loop
          } # end if

          # Apply mortality and ageing to tagged fish
          if(rseas < n_seas) {

            # Season mortality within a given year, advance to next season same year/age
            conv_tag_fish_avail[ry, rseas + 1, tc, , , , ] =
              conv_tag_fish_avail[ry, rseas, tc, , , , ] *
              tmp_SAA[,,1,,]

          } else {

            # End of year mortality and age advancement (end of season)
            conv_tag_fish_avail[ry + 1, 1, tc, , , 2:n_ages, ] =
              conv_tag_fish_avail[ry, n_seas, tc, , , 1:(n_ages-1), ] *
              tmp_SAA[,,1,1:(n_ages - 1),]

            # Accumulate plus group
            conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, ] =
              conv_tag_fish_avail[ry + 1, 1, tc, , , n_ages, ] +
              conv_tag_fish_avail[ry, n_seas, tc, , , n_ages, ] *
              tmp_SAA[,,1,n_ages,]
          }

          # # Apply Baranov's to get predicted recaptures
          for(f in 1:n_fish_fleets) {
            for(p in 1:n_pop) {
              pred_conv_tag_fish_recap[ry,rseas,tc,p,,,,f] = conv_tag_fish_reporting[,y,f] *
                (tmp_ret_FAA[p,,1,,,f] / tmp_ZAA[p,,1,,]) *
                conv_tag_fish_avail[ry,rseas,tc,p,,,] *
                (1 - tmp_SAA[p,,1,,])
            } # end p loop
          } # end f loop


        }
      }
    }

  } # end if for using tagging data


  # Likelihood Equations -------------------------------------------------------------
  ## Fishery Likelihoods -----------------------------------------------------
  ### Retained Fishery Catches (Regional) ---------------------------------------------------------
  if(any(UseCatch == 1)) { # setup OSA

    valid_idx = which(UseCatch == 1)
    ObsCatch_map = arrayInd(valid_idx, dim(UseCatch))
    ObsCatch = log(ObsCatch[valid_idx])
    ObsCatch = RTMB::OBS(ObsCatch)

    # compute nLL
    for(i in seq_along(ObsCatch)) {
      r    = ObsCatch_map[i, 1]
      y    = ObsCatch_map[i, 2]
      seas = ObsCatch_map[i, 3]
      f    = ObsCatch_map[i, 4]

      Catch_nLL[r,y,seas,f] = -1 * RTMB::dnorm(ObsCatch[i],
                                               log(sum(PredCatch[,r,y,seas,f])),
                                               exp(ln_sigmaC[r,y,seas,f]), TRUE)
    }
  }

  ### Retained Fishery Catches (Population-Specific) -----------------------------------------------
  if(any(UseCatch_pop == 1)) { # setup OSA

    valid_idx_cp = which(UseCatch_pop == 1)
    ObsCatch_pop_map = arrayInd(valid_idx_cp, dim(UseCatch_pop))
    ObsCatch_pop = log(ObsCatch_pop[valid_idx_cp])
    ObsCatch_pop = RTMB::OBS(ObsCatch_pop)

    # compute nLL
    for(i in seq_along(ObsCatch_pop)) {
      p    = ObsCatch_pop_map[i, 1]
      r    = ObsCatch_pop_map[i, 2]
      y    = ObsCatch_pop_map[i, 3]
      seas = ObsCatch_pop_map[i, 4]
      f    = ObsCatch_pop_map[i, 5]

      Catch_pop_nLL[p,r,y,seas,f] = -1 * RTMB::dnorm(ObsCatch_pop[i],
                                                     log(PredCatch[p,r,y,seas,f]),
                                                     exp(ln_sigmaC_pop[p,r,y,seas,f]), TRUE)
    }
  }


  ### Discarded Fishery Discards (Regional) --------------------------------------------------------
  if(any(UseDiscard == 1)) { # setup OSA

    valid_idx_dr = which(UseDiscard == 1)
    ObsDiscard_map = arrayInd(valid_idx_dr, dim(UseDiscard))
    ObsDiscard = log(ObsDiscard[valid_idx_dr])
    ObsDiscard = RTMB::OBS(ObsDiscard)

    # compute nLL
    for(i in seq_along(ObsDiscard)) {
      r    = ObsDiscard_map[i, 1]
      y    = ObsDiscard_map[i, 2]
      seas = ObsDiscard_map[i, 3]
      f    = ObsDiscard_map[i, 4]

      Discard_nLL[r,y,seas,f] = -1 * RTMB::dnorm(ObsDiscard[i],
                                                 log(sum(PredDiscard[,r,y,seas,f])),
                                                 exp(ln_sigmaD[r,y,seas,f]), TRUE)
    }
  }


  ### Discarded Fishery Discards (Population-Specific) ----------------------------------------------
  if(any(UseDiscard_pop == 1)) { # setup OSA

    valid_idx_dp = which(UseDiscard_pop == 1)
    ObsDiscard_pop_map = arrayInd(valid_idx_dp, dim(UseDiscard_pop))
    ObsDiscard_pop = log(ObsDiscard_pop[valid_idx_dp])
    ObsDiscard_pop = RTMB::OBS(ObsDiscard_pop)

    # compute nLL
    for(i in seq_along(ObsDiscard_pop)) {
      p    = ObsDiscard_pop_map[i, 1]
      r    = ObsDiscard_pop_map[i, 2]
      y    = ObsDiscard_pop_map[i, 3]
      seas = ObsDiscard_pop_map[i, 4]
      f    = ObsDiscard_pop_map[i, 5]

      Discard_pop_nLL[p,r,y,seas,f] = -1 * RTMB::dnorm(ObsDiscard_pop[i],
                                                       log(PredDiscard[p,r,y,seas,f]),
                                                       exp(ln_sigmaD_pop[p,r,y,seas,f]), TRUE)
    }
  }


  ### Retained Fishery Indices (Regional) -----------------------------------------------------------
  if(any(UseFishIdx == 1)) { # setup OSA

    valid_idx_ir = which(UseFishIdx == 1)
    ObsFishIdx_map = arrayInd(valid_idx_ir, dim(UseFishIdx))
    ObsFishIdx = log(ObsFishIdx[valid_idx_ir] + addtofishidx)
    ObsFishIdx = RTMB::OBS(ObsFishIdx)

    # compute nLL
    for(i in seq_along(ObsFishIdx)) {
      r    = ObsFishIdx_map[i, 1]
      y    = ObsFishIdx_map[i, 2]
      seas = ObsFishIdx_map[i, 3]
      f    = ObsFishIdx_map[i, 4]

      FishIdx_nLL[r,y,seas,f] = -1 * RTMB::dnorm(ObsFishIdx[i],
                                                 log(sum(PredFishIdx[,r,y,seas,f] + addtofishidx)),
                                                 ObsFishIdx_SE[r,y,seas,f], TRUE)
    }
  }


  ### Retained Fishery Indices (Population-Specific) ------------------------------------------------
  if(any(UseFishIdx_pop == 1)) { # setup OSA

    valid_idx_ip = which(UseFishIdx_pop == 1)
    ObsFishIdx_pop_map = arrayInd(valid_idx_ip, dim(UseFishIdx_pop))
    ObsFishIdx_pop = log(ObsFishIdx_pop[valid_idx_ip] + addtofishidx)
    ObsFishIdx_pop = RTMB::OBS(ObsFishIdx_pop)

    # compute nLL
    for(i in seq_along(ObsFishIdx_pop)) {
      p    = ObsFishIdx_pop_map[i, 1]
      r    = ObsFishIdx_pop_map[i, 2]
      y    = ObsFishIdx_pop_map[i, 3]
      seas = ObsFishIdx_pop_map[i, 4]
      f    = ObsFishIdx_pop_map[i, 5]

      FishIdx_pop_nLL[p,r,y,seas,f] = -1 * RTMB::dnorm(ObsFishIdx_pop[i],
                                                       log(PredFishIdx[p,r,y,seas,f] + addtofishidx),
                                                       ObsFishIdx_pop_SE[p,r,y,seas,f], TRUE)
    }
  }

  ### Retained Fishery Compositions (Region-Specific) ------------------------------------------------
  if(do_internal_comp_osa == FALSE) {

    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        for(seas in 1:n_seas) {
          # Fishery Age Compositions
          if(sum(UseFishAgeComps[,y,seas,f]) >= 1) {
            FishAgeComps_nLL[,y,seas,,f] = Get_Comp_Likelihoods(
              Exp = apply(CAA[,,y,seas,,,f, drop = FALSE], 2:7, sum),
              Obs = ObsFishAgeComps[,y,seas,,,f], ISS = ISS_FishAgeComps[,y,seas,,f], Wt_Mltnml = Wt_FishAgeComps[,y,seas,,f],
              Comp_Type = FishAgeComps_Type[y,f], Likelihood_Type = FishAgeComps_LikeType[f], ln_theta = ln_FishAge_theta[,,f],
              ln_theta_agg = ln_FishAge_theta_agg[f], LN_corr_pars = FishAge_corr_pars[,,f,], LN_corr_pars_agg = FishAge_corr_pars_agg[f],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError[y,,], use = UseFishAgeComps[,y,seas,f],
              n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps)[4], addtocomp = addtocomp
            )
          } # if we have fishery age comps

          # Fishery Length Compositions
          if(sum(UseFishLenComps[,y,seas,f]) >= 1 && fit_lengths == 1) {
            FishLenComps_nLL[,y,seas,,f] = Get_Comp_Likelihoods(
              Exp = apply(CAL[,,y,seas,,,f, drop = FALSE], 2:7, sum), Obs = ObsFishLenComps[,y,seas,,,f],
              ISS = ISS_FishLenComps[,y,seas,,f], Wt_Mltnml = Wt_FishLenComps[,y,seas,,f], Comp_Type = FishLenComps_Type[y,f],
              Likelihood_Type = FishLenComps_LikeType[f], ln_theta = ln_FishLen_theta[,,f], ln_theta_agg = ln_FishLen_theta_agg[f],
              LN_corr_pars = FishLen_corr_pars[,,f,], LN_corr_pars_agg = FishLen_corr_pars_agg[f],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = NA,
              use = UseFishLenComps[,y,seas,f], n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps)[4], addtocomp = addtocomp
            )
          } # if we have fishery length comps
        } # end seas loop

      } # end f loop
    } # end y loop

  } else {

    # Using OSA Compositions
    # Fishery Age Compositions
    if(any(UseFishAgeComps == 1)) {

      # Discrete Likelihoods
      ObsFishAgeComps_osa_discrete = NULL
      ObsFishAgeComps_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishAgeComps, ISSArr = ISS_FishAgeComps, WtArr = Wt_FishAgeComps,
        UseArr = UseFishAgeComps, TypeMat = FishAgeComps_Type, LikeTypeVec = FishAgeComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete"
      )

      if(!is.null(ObsFishAgeComps_osa_discrete)) {
        ObsFishAgeComps_osa_discrete = RTMB::OBS(ObsFishAgeComps_osa_discrete)
        FishAgeComps_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_nLL, tracked = ObsFishAgeComps_osa_discrete,
          ExpArrFn = function(p, y, seas, f) apply(CAA[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishAgeComps, TypeMat = FishAgeComps_Type, LikeTypeVec = FishAgeComps_LikeType,
          ISSArr = ISS_FishAgeComps, lnThetaArr = ln_FishAge_theta, lnThetaAggVec = ln_FishAge_theta_agg,
          LNcorrArr = FishAge_corr_pars, LNcorrAggVec = FishAge_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps)[4], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE
        )
      }

      # Continuous
      ObsFishAgeComps_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishAgeComps, ISSArr = ISS_FishAgeComps, WtArr = Wt_FishAgeComps,
        UseArr = UseFishAgeComps, TypeMat = FishAgeComps_Type, LikeTypeVec = FishAgeComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous"
      )

      if(!is.null(ObsFishAgeComps_osa_continuous)) {
        ObsFishAgeComps_osa_continuous = RTMB::OBS(ObsFishAgeComps_osa_continuous)
        FishAgeComps_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_nLL, tracked = ObsFishAgeComps_osa_continuous,
          ExpArrFn = function(p, y, seas, f) apply(CAA[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishAgeComps, TypeMat = FishAgeComps_Type, LikeTypeVec = FishAgeComps_LikeType,
          ISSArr = ISS_FishAgeComps, lnThetaArr = ln_FishAge_theta, lnThetaAggVec = ln_FishAge_theta_agg,
          LNcorrArr = FishAge_corr_pars, LNcorrAggVec = FishAge_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps)[4], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishAgeComps_osa_discrete)
        )
      }
    }


    # Fishery Length Compositions
    if(fit_lengths == 1 && any(UseFishLenComps == 1)) {

      # Discrete
      ObsFishLenComps_osa_discrete = NULL
      ObsFishLenComps_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishLenComps, ISSArr = ISS_FishLenComps, WtArr = Wt_FishLenComps,
        UseArr = UseFishLenComps, TypeMat = FishLenComps_Type, LikeTypeVec = FishLenComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete"
      )
      if(!is.null(ObsFishLenComps_osa_discrete)) {
        ObsFishLenComps_osa_discrete = RTMB::OBS(ObsFishLenComps_osa_discrete)
        FishLenComps_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_nLL, tracked = ObsFishLenComps_osa_discrete,
          ExpArrFn = function(p, y, seas, f) apply(CAL[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishLenComps, TypeMat = FishLenComps_Type, LikeTypeVec = FishLenComps_LikeType,
          ISSArr = ISS_FishLenComps, lnThetaArr = ln_FishLen_theta, lnThetaAggVec = ln_FishLen_theta_agg,
          LNcorrArr = FishLen_corr_pars, LNcorrAggVec = FishLen_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps)[4], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE
        )
      }

      # Continuous
      ObsFishLenComps_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishLenComps, ISSArr = ISS_FishLenComps, WtArr = Wt_FishLenComps,
        UseArr = UseFishLenComps, TypeMat = FishLenComps_Type, LikeTypeVec = FishLenComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous"
      )

      if(!is.null(ObsFishLenComps_osa_continuous)) {
        ObsFishLenComps_osa_continuous = RTMB::OBS(ObsFishLenComps_osa_continuous)
        FishLenComps_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_nLL, tracked = ObsFishLenComps_osa_continuous,
          ExpArrFn = function(p, y, seas, f) apply(CAL[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishLenComps, TypeMat = FishLenComps_Type, LikeTypeVec = FishLenComps_LikeType,
          ISSArr = ISS_FishLenComps, lnThetaArr = ln_FishLen_theta, lnThetaAggVec = ln_FishLen_theta_agg,
          LNcorrArr = FishLen_corr_pars, LNcorrAggVec = FishLen_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps)[4], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishLenComps_osa_discrete)
        )
      }
    }

  }


  ### Retained Fishery Compositions (Population-Specific) ------------------------------------------------
  if(do_internal_comp_osa == FALSE) {

    # Fishery Age Compositions
    for(p in 1:n_pop) {
      for(y in 1:n_yrs) {
        for(f in 1:n_fish_fleets) {
          for(seas in 1:n_seas) {
            if(sum(UseFishAgeComps_pop[p,,y,seas,f]) >= 1) {
              FishAgeComps_pop_nLL[p,,y,seas,,f] = Get_Comp_Likelihoods(
                Exp = CAA[p,,y,seas,,,f], Obs = ObsFishAgeComps_pop[p,,y,seas,,,f], ISS = ISS_FishAgeComps_pop[p,,y,seas,,f],
                Wt_Mltnml = Wt_FishAgeComps_pop[p,,y,seas,,f], Comp_Type = FishAgeComps_pop_Type[y,f],
                Likelihood_Type = FishAgeComps_pop_LikeType[f], ln_theta = ln_FishAge_pop_theta[p,,,f],
                ln_theta_agg = ln_FishAge_pop_theta_agg[p,f], LN_corr_pars = FishAge_pop_corr_pars[p,,,f,],
                LN_corr_pars_agg = FishAge_pop_corr_pars_agg[p,f], n_regions = n_regions, n_sexes = n_sexes,
                age_or_len = 0, AgeingError = AgeingError[y,,], use = UseFishAgeComps_pop[p,,y,seas,f],
                n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_pop)[5], addtocomp = addtocomp
              )
            }

            # Fishery Length Compositions
            if(sum(UseFishLenComps_pop[p,,y,seas,f]) >= 1 && fit_lengths == 1) {
              FishLenComps_pop_nLL[p,,y,seas,,f] = Get_Comp_Likelihoods(
                Exp = CAL[p,,y,seas,,,f], Obs = ObsFishLenComps_pop[p,,y,seas,,,f], ISS = ISS_FishLenComps_pop[p,,y,seas,,f],
                Wt_Mltnml = Wt_FishLenComps_pop[p,,y,seas,,f], Comp_Type = FishLenComps_pop_Type[y,f],
                Likelihood_Type = FishLenComps_pop_LikeType[f], ln_theta = ln_FishLen_pop_theta[p,,,f],
                ln_theta_agg = ln_FishLen_pop_theta_agg[p,f], LN_corr_pars = FishLen_pop_corr_pars[p,,,f,],
                LN_corr_pars_agg = FishLen_pop_corr_pars_agg[p,f], n_regions = n_regions,
                n_sexes = n_sexes, age_or_len = 1, AgeingError = NA, use = UseFishLenComps_pop[p,,y,seas,f],
                n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_pop)[5], addtocomp = addtocomp
              )
            }
          }
        }
      }
    }
  } else {

    # Using OSA Compositions
    # Fishery Age Compositions
    if(any(UseFishAgeComps_pop == 1)) {

      # Discrete
      ObsFishAgeComps_pop_osa_discrete = NULL
      ObsFishAgeComps_pop_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishAgeComps_pop, ISSArr = ISS_FishAgeComps_pop, WtArr = Wt_FishAgeComps_pop,
        UseArr = UseFishAgeComps_pop, TypeMat = FishAgeComps_pop_Type, LikeTypeVec = FishAgeComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishAgeComps_pop_osa_discrete)) {
        ObsFishAgeComps_pop_osa_discrete = RTMB::OBS(ObsFishAgeComps_pop_osa_discrete)
        FishAgeComps_pop_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_pop_nLL, tracked = ObsFishAgeComps_pop_osa_discrete,
          ExpArrFn = function(p, y, seas, f) { e = CAA[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_ages, n_sexes); e },
          UseArr = UseFishAgeComps_pop, TypeMat = FishAgeComps_pop_Type, LikeTypeVec = FishAgeComps_pop_LikeType,
          ISSArr = ISS_FishAgeComps_pop, lnThetaArr = ln_FishAge_pop_theta, lnThetaAggVec = ln_FishAge_pop_theta_agg,
          LNcorrArr = FishAge_pop_corr_pars, LNcorrAggVec = FishAge_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_pop)[5], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE, pop = TRUE, n_pop = n_pop
        )
      }

      # Continuous
      ObsFishAgeComps_pop_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishAgeComps_pop, ISSArr = ISS_FishAgeComps_pop, WtArr = Wt_FishAgeComps_pop,
        UseArr = UseFishAgeComps_pop, TypeMat = FishAgeComps_pop_Type, LikeTypeVec = FishAgeComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishAgeComps_pop_osa_continuous)) {
        ObsFishAgeComps_pop_osa_continuous = RTMB::OBS(ObsFishAgeComps_pop_osa_continuous)
        FishAgeComps_pop_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_pop_nLL, tracked = ObsFishAgeComps_pop_osa_continuous,
          ExpArrFn = function(p, y, seas, f) { e = CAA[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_ages, n_sexes); e },
          UseArr = UseFishAgeComps_pop, TypeMat = FishAgeComps_pop_Type, LikeTypeVec = FishAgeComps_pop_LikeType,
          ISSArr = ISS_FishAgeComps_pop, lnThetaArr = ln_FishAge_pop_theta, lnThetaAggVec = ln_FishAge_pop_theta_agg,
          LNcorrArr = FishAge_pop_corr_pars, LNcorrAggVec = FishAge_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_pop)[5], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishAgeComps_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }

    # Fishery Length Compositions
    if(fit_lengths == 1 && any(UseFishLenComps_pop == 1)) {

      # Discrete
      ObsFishLenComps_pop_osa_discrete = NULL
      ObsFishLenComps_pop_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishLenComps_pop, ISSArr = ISS_FishLenComps_pop, WtArr = Wt_FishLenComps_pop,
        UseArr = UseFishLenComps_pop, TypeMat = FishLenComps_pop_Type, LikeTypeVec = FishLenComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishLenComps_pop_osa_discrete)) {
        ObsFishLenComps_pop_osa_discrete = RTMB::OBS(ObsFishLenComps_pop_osa_discrete)
        FishLenComps_pop_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_pop_nLL, tracked = ObsFishLenComps_pop_osa_discrete,
          ExpArrFn = function(p, y, seas, f) { e = CAL[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_lens, n_sexes); e },
          UseArr = UseFishLenComps_pop, TypeMat = FishLenComps_pop_Type, LikeTypeVec = FishLenComps_pop_LikeType,
          ISSArr = ISS_FishLenComps_pop, lnThetaArr = ln_FishLen_pop_theta, lnThetaAggVec = ln_FishLen_pop_theta_agg,
          LNcorrArr = FishLen_pop_corr_pars, LNcorrAggVec = FishLen_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_pop)[5], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE, pop = TRUE, n_pop = n_pop
        )
      }

      # Continuous
      ObsFishLenComps_pop_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishLenComps_pop, ISSArr = ISS_FishLenComps_pop, WtArr = Wt_FishLenComps_pop,
        UseArr = UseFishLenComps_pop, TypeMat = FishLenComps_pop_Type, LikeTypeVec = FishLenComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishLenComps_pop_osa_continuous)) {
        ObsFishLenComps_pop_osa_continuous = RTMB::OBS(ObsFishLenComps_pop_osa_continuous)
        FishLenComps_pop_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_pop_nLL, tracked = ObsFishLenComps_pop_osa_continuous,
          ExpArrFn = function(p, y, seas, f) { e = CAL[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_lens, n_sexes); e },
          UseArr = UseFishLenComps_pop, TypeMat = FishLenComps_pop_Type, LikeTypeVec = FishLenComps_pop_LikeType,
          ISSArr = ISS_FishLenComps_pop, lnThetaArr = ln_FishLen_pop_theta, lnThetaAggVec = ln_FishLen_pop_theta_agg,
          LNcorrArr = FishLen_pop_corr_pars, LNcorrAggVec = FishLen_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_pop)[5], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishLenComps_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }

  }

  ### Discarded Fishery Compositions (Region-Specific) ------------------------------------------------
  if(do_internal_comp_osa == FALSE) {

    # Discarded Fishery Age Compositions
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        for(seas in 1:n_seas) {
          if(sum(UseFishAgeComps_discard[,y,seas,f]) >= 1) {
            FishAgeComps_discard_nLL[,y,seas,,f] = Get_Comp_Likelihoods(
              Exp = apply(DAA[,,y,seas,,,f, drop = FALSE], 2:7, sum), Obs = ObsFishAgeComps_discard[,y,seas,,,f],
              ISS = ISS_FishAgeComps_discard[,y,seas,,f], Wt_Mltnml = Wt_FishAgeComps_discard[,y,seas,,f],
              Comp_Type = FishAgeComps_discard_Type[y,f], Likelihood_Type = FishAgeComps_discard_LikeType[f],
              ln_theta = ln_FishAge_discard_theta[,,f], ln_theta_agg = ln_FishAge_discard_theta_agg[f],
              LN_corr_pars = FishAge_discard_corr_pars[,,f,], LN_corr_pars_agg = FishAge_discard_corr_pars_agg[f],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError[y,,],
              use = UseFishAgeComps_discard[,y,seas,f], n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_discard)[4],
              addtocomp = addtocomp
            )
          }

          # Discarded Fishery Length Compositions
          if(sum(UseFishLenComps_discard[,y,seas,f]) >= 1 && fit_lengths == 1) {
            FishLenComps_discard_nLL[,y,seas,,f] = Get_Comp_Likelihoods(
              Exp = apply(DAL[,,y,seas,,,f, drop = FALSE], 2:7, sum), Obs = ObsFishLenComps_discard[,y,seas,,,f],
              ISS = ISS_FishLenComps_discard[,y,seas,,f], Wt_Mltnml = Wt_FishLenComps_discard[,y,seas,,f],
              Comp_Type = FishLenComps_discard_Type[y,f], Likelihood_Type = FishLenComps_discard_LikeType[f],
              ln_theta = ln_FishLen_discard_theta[,,f], ln_theta_agg = ln_FishLen_discard_theta_agg[f],
              LN_corr_pars = FishLen_discard_corr_pars[,,f,], LN_corr_pars_agg = FishLen_discard_corr_pars_agg[f],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = NA,
              use = UseFishLenComps_discard[,y,seas,f], n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_discard)[4],
              addtocomp = addtocomp
            )
          }
        }
      }
    }
  } else {

    # Doing OSA Compositions
    # Discarded Fishery Age Compositions
    if(any(UseFishAgeComps_discard == 1)) {

      # Discrete
      ObsFishAgeComps_discard_osa_discrete = NULL
      ObsFishAgeComps_discard_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishAgeComps_discard, ISSArr = ISS_FishAgeComps_discard, WtArr = Wt_FishAgeComps_discard,
        UseArr = UseFishAgeComps_discard, TypeMat = FishAgeComps_discard_Type, LikeTypeVec = FishAgeComps_discard_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete"
      )
      if(!is.null(ObsFishAgeComps_discard_osa_discrete)) {
        ObsFishAgeComps_discard_osa_discrete = RTMB::OBS(ObsFishAgeComps_discard_osa_discrete)
        FishAgeComps_discard_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_discard_nLL, tracked = ObsFishAgeComps_discard_osa_discrete,
          ExpArrFn = function(p, y, seas, f) apply(DAA[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishAgeComps_discard, TypeMat = FishAgeComps_discard_Type, LikeTypeVec = FishAgeComps_discard_LikeType,
          ISSArr = ISS_FishAgeComps_discard, lnThetaArr = ln_FishAge_discard_theta, lnThetaAggVec = ln_FishAge_discard_theta_agg,
          LNcorrArr = FishAge_discard_corr_pars, LNcorrAggVec = FishAge_discard_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_discard)[4], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE
        )
      }

      # Continuous
      ObsFishAgeComps_discard_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishAgeComps_discard, ISSArr = ISS_FishAgeComps_discard, WtArr = Wt_FishAgeComps_discard,
        UseArr = UseFishAgeComps_discard, TypeMat = FishAgeComps_discard_Type, LikeTypeVec = FishAgeComps_discard_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous"
      )
      if(!is.null(ObsFishAgeComps_discard_osa_continuous)) {
        ObsFishAgeComps_discard_osa_continuous = RTMB::OBS(ObsFishAgeComps_discard_osa_continuous)
        FishAgeComps_discard_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_discard_nLL, tracked = ObsFishAgeComps_discard_osa_continuous,
          ExpArrFn = function(p, y, seas, f) apply(DAA[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishAgeComps_discard, TypeMat = FishAgeComps_discard_Type, LikeTypeVec = FishAgeComps_discard_LikeType,
          ISSArr = ISS_FishAgeComps_discard, lnThetaArr = ln_FishAge_discard_theta, lnThetaAggVec = ln_FishAge_discard_theta_agg,
          LNcorrArr = FishAge_discard_corr_pars, LNcorrAggVec = FishAge_discard_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_discard)[4], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishAgeComps_discard_osa_discrete)
        )
      }
    }

    # Fishery Length Compositions discards
    if(fit_lengths == 1 && any(UseFishLenComps_discard == 1)) {

      # Discrete
      ObsFishLenComps_discard_osa_discrete = NULL
      ObsFishLenComps_discard_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishLenComps_discard, ISSArr = ISS_FishLenComps_discard, WtArr = Wt_FishLenComps_discard,
        UseArr = UseFishLenComps_discard, TypeMat = FishLenComps_discard_Type, LikeTypeVec = FishLenComps_discard_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete"
      )
      if(!is.null(ObsFishLenComps_discard_osa_discrete)) {
        ObsFishLenComps_discard_osa_discrete = RTMB::OBS(ObsFishLenComps_discard_osa_discrete)
        FishLenComps_discard_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_discard_nLL, tracked = ObsFishLenComps_discard_osa_discrete,
          ExpArrFn = function(p, y, seas, f) apply(DAL[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishLenComps_discard, TypeMat = FishLenComps_discard_Type, LikeTypeVec = FishLenComps_discard_LikeType,
          ISSArr = ISS_FishLenComps_discard, lnThetaArr = ln_FishLen_discard_theta, lnThetaAggVec = ln_FishLen_discard_theta_agg,
          LNcorrArr = FishLen_discard_corr_pars, LNcorrAggVec = FishLen_discard_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_discard)[4], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE
        )
      }

      # Continuous
      ObsFishLenComps_discard_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishLenComps_discard, ISSArr = ISS_FishLenComps_discard, WtArr = Wt_FishLenComps_discard,
        UseArr = UseFishLenComps_discard, TypeMat = FishLenComps_discard_Type, LikeTypeVec = FishLenComps_discard_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous"
      )
      if(!is.null(ObsFishLenComps_discard_osa_continuous)) {
        ObsFishLenComps_discard_osa_continuous = RTMB::OBS(ObsFishLenComps_discard_osa_continuous)
        FishLenComps_discard_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_discard_nLL, tracked = ObsFishLenComps_discard_osa_continuous,
          ExpArrFn = function(p, y, seas, f) apply(DAL[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseFishLenComps_discard, TypeMat = FishLenComps_discard_Type, LikeTypeVec = FishLenComps_discard_LikeType,
          ISSArr = ISS_FishLenComps_discard, lnThetaArr = ln_FishLen_discard_theta, lnThetaAggVec = ln_FishLen_discard_theta_agg,
          LNcorrArr = FishLen_discard_corr_pars, LNcorrAggVec = FishLen_discard_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_discard)[4], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishLenComps_discard_osa_discrete)
        )
      }
    }
  }


  ### Discarded Fishery Compositions (Population-Specific) ------------------------------------------------
  if(do_internal_comp_osa == FALSE) {
    for(p in 1:n_pop) {
      for(y in 1:n_yrs) {
        for(f in 1:n_fish_fleets) {
          for(seas in 1:n_seas) {

            # Discarded Fishery Age Compositions
            if(sum(UseFishAgeComps_discard_pop[p,,y,seas,f]) >= 1) {
              FishAgeComps_discard_pop_nLL[p,,y,seas,,f] = Get_Comp_Likelihoods(
                Exp = DAA[p,,y,seas,,,f], Obs = ObsFishAgeComps_discard_pop[p,,y,seas,,,f], ISS = ISS_FishAgeComps_discard_pop[p,,y,seas,,f],
                Wt_Mltnml = Wt_FishAgeComps_discard_pop[p,,y,seas,,f], Comp_Type = FishAgeComps_discard_pop_Type[y,f],
                Likelihood_Type = FishAgeComps_discard_pop_LikeType[f], ln_theta = ln_FishAge_discard_pop_theta[p,,,f],
                ln_theta_agg = ln_FishAge_discard_pop_theta_agg[p,f], LN_corr_pars = FishAge_discard_pop_corr_pars[p,,,f,],
                LN_corr_pars_agg = FishAge_discard_pop_corr_pars_agg[p,f], n_regions = n_regions, n_sexes = n_sexes,
                age_or_len = 0, AgeingError = AgeingError[y,,], use = UseFishAgeComps_discard_pop[p,,y,seas,f],
                n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_discard_pop)[5], addtocomp = addtocomp
              )
            }

            # Discarded Fishery Length Compositions
            if(sum(UseFishLenComps_discard_pop[p,,y,seas,f]) >= 1 && fit_lengths == 1) {
              FishLenComps_discard_pop_nLL[p,,y,seas,,f] = Get_Comp_Likelihoods(
                Exp = DAL[p,,y,seas,,,f], Obs = ObsFishLenComps_discard_pop[p,,y,seas,,,f], ISS = ISS_FishLenComps_discard_pop[p,,y,seas,,f],
                Wt_Mltnml = Wt_FishLenComps_discard_pop[p,,y,seas,,f], Comp_Type = FishLenComps_discard_pop_Type[y,f],
                Likelihood_Type = FishLenComps_discard_pop_LikeType[f], ln_theta = ln_FishLen_discard_pop_theta[p,,,f],
                ln_theta_agg = ln_FishLen_discard_pop_theta_agg[p,f], LN_corr_pars = FishLen_discard_pop_corr_pars[p,,,f,],
                LN_corr_pars_agg = FishLen_discard_pop_corr_pars_agg[p,f], n_regions = n_regions, n_sexes = n_sexes,
                age_or_len = 1, AgeingError = NA, use = UseFishLenComps_discard_pop[p,,y,seas,f], n_model_bins = n_lens,
                n_obs_bins = dim(ObsFishLenComps_discard_pop)[5], addtocomp = addtocomp
              )
            }
          }
        }
      }
    }
  } else {

    # Doing OSA Compositions
    # Fishery Age Discard Compositions
    if(any(UseFishAgeComps_discard_pop == 1)) {

      # Discrete
      ObsFishAgeComps_discard_pop_osa_discrete = NULL
      ObsFishAgeComps_discard_pop_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishAgeComps_discard_pop, ISSArr = ISS_FishAgeComps_discard_pop, WtArr = Wt_FishAgeComps_discard_pop,
        UseArr = UseFishAgeComps_discard_pop, TypeMat = FishAgeComps_discard_pop_Type, LikeTypeVec = FishAgeComps_discard_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishAgeComps_discard_pop_osa_discrete)) {
        ObsFishAgeComps_discard_pop_osa_discrete = RTMB::OBS(ObsFishAgeComps_discard_pop_osa_discrete)
        FishAgeComps_discard_pop_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_discard_pop_nLL, tracked = ObsFishAgeComps_discard_pop_osa_discrete,
          ExpArrFn = function(p, y, seas, f) { e = DAA[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_ages, n_sexes); e },
          UseArr = UseFishAgeComps_discard_pop, TypeMat = FishAgeComps_discard_pop_Type, LikeTypeVec = FishAgeComps_discard_pop_LikeType,
          ISSArr = ISS_FishAgeComps_discard_pop, lnThetaArr = ln_FishAge_discard_pop_theta, lnThetaAggVec = ln_FishAge_discard_pop_theta_agg,
          LNcorrArr = FishAge_discard_pop_corr_pars, LNcorrAggVec = FishAge_discard_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_discard_pop)[5], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE, pop = TRUE, n_pop = n_pop
        )
      }

      # Continuous
      ObsFishAgeComps_discard_pop_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishAgeComps_discard_pop, ISSArr = ISS_FishAgeComps_discard_pop, WtArr = Wt_FishAgeComps_discard_pop,
        UseArr = UseFishAgeComps_discard_pop, TypeMat = FishAgeComps_discard_pop_Type, LikeTypeVec = FishAgeComps_discard_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishAgeComps_discard_pop_osa_continuous)) {
        ObsFishAgeComps_discard_pop_osa_continuous = RTMB::OBS(ObsFishAgeComps_discard_pop_osa_continuous)
        FishAgeComps_discard_pop_nLL = eval_comp_osa(
          nLL_arr = FishAgeComps_discard_pop_nLL, tracked = ObsFishAgeComps_discard_pop_osa_continuous,
          ExpArrFn = function(p, y, seas, f) { e = DAA[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_ages, n_sexes); e },
          UseArr = UseFishAgeComps_discard_pop, TypeMat = FishAgeComps_discard_pop_Type, LikeTypeVec = FishAgeComps_discard_pop_LikeType,
          ISSArr = ISS_FishAgeComps_discard_pop, lnThetaArr = ln_FishAge_discard_pop_theta, lnThetaAggVec = ln_FishAge_discard_pop_theta_agg,
          LNcorrArr = FishAge_discard_pop_corr_pars, LNcorrAggVec = FishAge_discard_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps_discard_pop)[5], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishAgeComps_discard_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }

    # Fishery Length Discarded Compositions
    if(fit_lengths == 1 && any(UseFishLenComps_discard_pop == 1)) {

      # Discrete
      ObsFishLenComps_discard_pop_osa_discrete = NULL
      ObsFishLenComps_discard_pop_osa_discrete = pack_comp_osa(
        ObsArr = ObsFishLenComps_discard_pop, ISSArr = ISS_FishLenComps_discard_pop, WtArr = Wt_FishLenComps_discard_pop,
        UseArr = UseFishLenComps_discard_pop, TypeMat = FishLenComps_discard_pop_Type, LikeTypeVec = FishLenComps_discard_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishLenComps_discard_pop_osa_discrete)) {
        ObsFishLenComps_discard_pop_osa_discrete = RTMB::OBS(ObsFishLenComps_discard_pop_osa_discrete)
        FishLenComps_discard_pop_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_discard_pop_nLL, tracked = ObsFishLenComps_discard_pop_osa_discrete,
          ExpArrFn = function(p, y, seas, f) { e = DAL[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_lens, n_sexes); e },
          UseArr = UseFishLenComps_discard_pop, TypeMat = FishLenComps_discard_pop_Type, LikeTypeVec = FishLenComps_discard_pop_LikeType,
          ISSArr = ISS_FishLenComps_discard_pop, lnThetaArr = ln_FishLen_discard_pop_theta, lnThetaAggVec = ln_FishLen_discard_pop_theta_agg,
          LNcorrArr = FishLen_discard_pop_corr_pars, LNcorrAggVec = FishLen_discard_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_discard_pop)[5], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE, pop = TRUE, n_pop = n_pop
        )
      }

      # Continuous
      ObsFishLenComps_discard_pop_osa_continuous = pack_comp_osa(
        ObsArr = ObsFishLenComps_discard_pop, ISSArr = ISS_FishLenComps_discard_pop, WtArr = Wt_FishLenComps_discard_pop,
        UseArr = UseFishLenComps_discard_pop, TypeMat = FishLenComps_discard_pop_Type, LikeTypeVec = FishLenComps_discard_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsFishLenComps_discard_pop_osa_continuous)) {
        ObsFishLenComps_discard_pop_osa_continuous = RTMB::OBS(ObsFishLenComps_discard_pop_osa_continuous)
        FishLenComps_discard_pop_nLL = eval_comp_osa(
          nLL_arr = FishLenComps_discard_pop_nLL, tracked = ObsFishLenComps_discard_pop_osa_continuous,
          ExpArrFn = function(p, y, seas, f) { e = DAL[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_lens, n_sexes); e },
          UseArr = UseFishLenComps_discard_pop, TypeMat = FishLenComps_discard_pop_Type, LikeTypeVec = FishLenComps_discard_pop_LikeType,
          ISSArr = ISS_FishLenComps_discard_pop, lnThetaArr = ln_FishLen_discard_pop_theta, lnThetaAggVec = ln_FishLen_discard_pop_theta_agg,
          LNcorrArr = FishLen_discard_pop_corr_pars, LNcorrAggVec = FishLen_discard_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsFishLenComps_discard_pop)[5], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishLenComps_discard_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }
  }

  ## Survey Likelihoods ------------------------------------------------------
  ### Survey Indices (Regional) ---------------------------------------------------------
  if(any(UseSrvIdx == 1)) { # setup OSA
    valid_idx_sr = which(UseSrvIdx == 1)
    ObsSrvIdx_map = arrayInd(valid_idx_sr, dim(UseSrvIdx))
    ObsSrvIdx = log(ObsSrvIdx[valid_idx_sr] + addtosrvidx)
    ObsSrvIdx = RTMB::OBS(ObsSrvIdx)

    # compute nLL
    for(i in seq_along(ObsSrvIdx)) {
      r    = ObsSrvIdx_map[i, 1]
      y    = ObsSrvIdx_map[i, 2]
      seas = ObsSrvIdx_map[i, 3]
      sf   = ObsSrvIdx_map[i, 4]

      SrvIdx_nLL[r,y,seas,sf] = -1 * RTMB::dnorm(ObsSrvIdx[i],
                                                 log(sum(PredSrvIdx[,r,y,seas,sf] + addtosrvidx)),
                                                 ObsSrvIdx_SE[r,y,seas,sf], TRUE)
    }
  }

  ### Survey Indices (Population-Specific) ---------------------------------------------------------
  if(any(UseSrvIdx_pop == 1)) { # setup OSA
    valid_idx_sp = which(UseSrvIdx_pop == 1)
    ObsSrvIdx_pop_map = arrayInd(valid_idx_sp, dim(UseSrvIdx_pop))
    ObsSrvIdx_pop = log(ObsSrvIdx_pop[valid_idx_sp] + addtosrvidx)
    ObsSrvIdx_pop = RTMB::OBS(ObsSrvIdx_pop)

    # compute nLL
    for(i in seq_along(ObsSrvIdx_pop)) {
      p    = ObsSrvIdx_pop_map[i, 1]
      r    = ObsSrvIdx_pop_map[i, 2]
      y    = ObsSrvIdx_pop_map[i, 3]
      seas   = ObsSrvIdx_pop_map[i, 4]
      sf = ObsSrvIdx_pop_map[i, 5]

      SrvIdx_pop_nLL[p,r,y,seas,sf] = -1 * RTMB::dnorm(ObsSrvIdx_pop[i],
                                                       log(PredSrvIdx[p,r,y,seas,sf] + addtosrvidx),
                                                       ObsSrvIdx_pop_SE[p,r,y,seas,sf], TRUE)
    }
  }

  ### Survey Compositions (Region-Specific) ---------------------------------------------------------
  if(do_internal_comp_osa == FALSE) {

    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {
        for(seas in 1:n_seas) {
          if(sum(UseSrvAgeComps[,y,seas,sf]) >= 1) {

            # Survey Age Compositions
            SrvAgeComps_nLL[,y,seas,,sf] = Get_Comp_Likelihoods(
              Exp = apply(SrvIAA[,,y,seas,,,sf, drop = FALSE], 2:7, sum), Obs = ObsSrvAgeComps[,y,seas,,,sf],
              ISS = ISS_SrvAgeComps[,y,seas,,sf], Wt_Mltnml = Wt_SrvAgeComps[,y,seas,,sf],
              Comp_Type = SrvAgeComps_Type[y,sf], Likelihood_Type = SrvAgeComps_LikeType[sf],
              ln_theta = ln_SrvAge_theta[,,sf], ln_theta_agg = ln_SrvAge_theta_agg[sf],
              LN_corr_pars = SrvAge_corr_pars[,,sf,], LN_corr_pars_agg = SrvAge_corr_pars_agg[sf],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError[y,,],
              use = UseSrvAgeComps[,y,seas,sf], n_model_bins = n_ages, n_obs_bins = dim(ObsSrvAgeComps)[4], addtocomp = addtocomp
            )
          }

          # Survey Length Compositions
          if(sum(UseSrvLenComps[,y,seas,sf]) >= 1 && fit_lengths == 1) {
            SrvLenComps_nLL[,y,seas,,sf] = Get_Comp_Likelihoods(
              Exp = apply(SrvIAL[,,y,seas,,,sf, drop = FALSE], 2:7, sum), Obs = ObsSrvLenComps[,y,seas,,,sf],
              ISS = ISS_SrvLenComps[,y,seas,,sf], Wt_Mltnml = Wt_SrvLenComps[,y,seas,,sf], Comp_Type = SrvLenComps_Type[y,sf],
              Likelihood_Type = SrvLenComps_LikeType[sf], ln_theta = ln_SrvLen_theta[,,sf], ln_theta_agg = ln_SrvLen_theta_agg[sf],
              LN_corr_pars = SrvLen_corr_pars[,,sf,], LN_corr_pars_agg = SrvLen_corr_pars_agg[sf], n_regions = n_regions,
              n_sexes = n_sexes, age_or_len = 1, AgeingError = NA, use = UseSrvLenComps[,y,seas,sf],
              n_model_bins = n_lens, n_obs_bins = dim(ObsSrvLenComps)[4], addtocomp = addtocomp
            )
          }
        }
      }
    }
  } else {

    # Doing OSA Compositions
    # Survey Age Compositions
    if(any(UseSrvAgeComps == 1)) {

      # Discrete
      ObsSrvAgeComps_osa_discrete = NULL
      ObsSrvAgeComps_osa_discrete = pack_comp_osa(
        ObsArr = ObsSrvAgeComps, ISSArr = ISS_SrvAgeComps, WtArr = Wt_SrvAgeComps,
        UseArr = UseSrvAgeComps, TypeMat = SrvAgeComps_Type, LikeTypeVec = SrvAgeComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete"
      )
      if(!is.null(ObsSrvAgeComps_osa_discrete)) {
        ObsSrvAgeComps_osa_discrete = RTMB::OBS(ObsSrvAgeComps_osa_discrete)
        SrvAgeComps_nLL = eval_comp_osa(
          nLL_arr = SrvAgeComps_nLL, tracked = ObsSrvAgeComps_osa_discrete,
          ExpArrFn = function(p, y, seas, f) apply(SrvIAA[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseSrvAgeComps, TypeMat = SrvAgeComps_Type, LikeTypeVec = SrvAgeComps_LikeType,
          ISSArr = ISS_SrvAgeComps, lnThetaArr = ln_SrvAge_theta, lnThetaAggVec = ln_SrvAge_theta_agg,
          LNcorrArr = SrvAge_corr_pars, LNcorrAggVec = SrvAge_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsSrvAgeComps)[4], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE
        )
      }

      # Continuous
      ObsSrvAgeComps_osa_continuous = pack_comp_osa(
        ObsArr = ObsSrvAgeComps, ISSArr = ISS_SrvAgeComps, WtArr = Wt_SrvAgeComps,
        UseArr = UseSrvAgeComps, TypeMat = SrvAgeComps_Type, LikeTypeVec = SrvAgeComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous"
      )
      if(!is.null(ObsSrvAgeComps_osa_continuous)) {
        ObsSrvAgeComps_osa_continuous = RTMB::OBS(ObsSrvAgeComps_osa_continuous)
        SrvAgeComps_nLL = eval_comp_osa(
          nLL_arr = SrvAgeComps_nLL, tracked = ObsSrvAgeComps_osa_continuous,
          ExpArrFn = function(p, y, seas, f) apply(SrvIAA[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseSrvAgeComps, TypeMat = SrvAgeComps_Type, LikeTypeVec = SrvAgeComps_LikeType,
          ISSArr = ISS_SrvAgeComps, lnThetaArr = ln_SrvAge_theta, lnThetaAggVec = ln_SrvAge_theta_agg,
          LNcorrArr = SrvAge_corr_pars, LNcorrAggVec = SrvAge_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsSrvAgeComps)[4], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsSrvAgeComps_osa_discrete)
        )
      }
    }

    # Survey Length Compositions
    if(fit_lengths == 1 && any(UseSrvLenComps == 1)) {

      # Discrete
      ObsSrvLenComps_osa_discrete = NULL
      ObsSrvLenComps_osa_discrete = pack_comp_osa(
        ObsArr = ObsSrvLenComps, ISSArr = ISS_SrvLenComps, WtArr = Wt_SrvLenComps,
        UseArr = UseSrvLenComps, TypeMat = SrvLenComps_Type, LikeTypeVec = SrvLenComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete"
      )
      if(!is.null(ObsSrvLenComps_osa_discrete)) {
        ObsSrvLenComps_osa_discrete = RTMB::OBS(ObsSrvLenComps_osa_discrete)
        SrvLenComps_nLL = eval_comp_osa(
          nLL_arr = SrvLenComps_nLL, tracked = ObsSrvLenComps_osa_discrete,
          ExpArrFn = function(p, y, seas, f) apply(SrvIAL[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseSrvLenComps, TypeMat = SrvLenComps_Type, LikeTypeVec = SrvLenComps_LikeType,
          ISSArr = ISS_SrvLenComps, lnThetaArr = ln_SrvLen_theta, lnThetaAggVec = ln_SrvLen_theta_agg,
          LNcorrArr = SrvLen_corr_pars, LNcorrAggVec = SrvLen_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsSrvLenComps)[4], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE
        )
      }

      # Continuous
      ObsSrvLenComps_osa_continuous = pack_comp_osa(
        ObsArr = ObsSrvLenComps, ISSArr = ISS_SrvLenComps, WtArr = Wt_SrvLenComps,
        UseArr = UseSrvLenComps, TypeMat = SrvLenComps_Type, LikeTypeVec = SrvLenComps_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous"
      )
      if(!is.null(ObsSrvLenComps_osa_continuous)) {
        ObsSrvLenComps_osa_continuous = RTMB::OBS(ObsSrvLenComps_osa_continuous)
        SrvLenComps_nLL = eval_comp_osa(
          nLL_arr = SrvLenComps_nLL, tracked = ObsSrvLenComps_osa_continuous,
          ExpArrFn = function(p, y, seas, f) apply(SrvIAL[,,y,seas,,,f, drop=FALSE], 2:7, sum),
          UseArr = UseSrvLenComps, TypeMat = SrvLenComps_Type, LikeTypeVec = SrvLenComps_LikeType,
          ISSArr = ISS_SrvLenComps, lnThetaArr = ln_SrvLen_theta, lnThetaAggVec = ln_SrvLen_theta_agg,
          LNcorrArr = SrvLen_corr_pars, LNcorrAggVec = SrvLen_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsSrvLenComps)[4], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsSrvLenComps_osa_discrete)
        )
      }
    }
  }

  ### Survey Compositions (Population-Specific) -------------------------------
  if(do_internal_comp_osa == FALSE) {
    for(p in 1:n_pop) {
      for(y in 1:n_yrs) {
        for(sf in 1:n_srv_fleets) {
          for(seas in 1:n_seas) {

            # Survey Age Compositions
            if(sum(UseSrvAgeComps_pop[p,,y,seas,sf]) >= 1) {
              SrvAgeComps_pop_nLL[p,,y,seas,,sf] = Get_Comp_Likelihoods(
                Exp = SrvIAA[p,,y,seas,,,sf], Obs = ObsSrvAgeComps_pop[p,,y,seas,,,sf],
                ISS = ISS_SrvAgeComps_pop[p,,y,seas,,sf], Wt_Mltnml = Wt_SrvAgeComps_pop[p,,y,seas,,sf],
                Comp_Type = SrvAgeComps_pop_Type[y,sf], Likelihood_Type = SrvAgeComps_pop_LikeType[sf],
                ln_theta = ln_SrvAge_pop_theta[p,,,sf], ln_theta_agg = ln_SrvAge_pop_theta_agg[p,sf],
                LN_corr_pars = SrvAge_pop_corr_pars[p,,,sf,], LN_corr_pars_agg = SrvAge_pop_corr_pars_agg[p,sf],
                n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError[y,,],
                use = UseSrvAgeComps_pop[p,,y,seas,sf], n_model_bins = n_ages, n_obs_bins = dim(ObsSrvAgeComps_pop)[5],
                addtocomp = addtocomp
              )
            }

            # Survey Length Compositions
            if(sum(UseSrvLenComps_pop[p,,y,seas,sf]) >= 1 && fit_lengths == 1) {
              SrvLenComps_pop_nLL[p,,y,seas,,sf] = Get_Comp_Likelihoods(
                Exp = SrvIAL[p,,y,seas,,,sf], Obs = ObsSrvLenComps_pop[p,,y,seas,,,sf],
                ISS = ISS_SrvLenComps_pop[p,,y,seas,,sf], Wt_Mltnml = Wt_SrvLenComps_pop[p,,y,seas,,sf],
                Comp_Type = SrvLenComps_pop_Type[y,sf], Likelihood_Type = SrvLenComps_pop_LikeType[sf],
                ln_theta = ln_SrvLen_pop_theta[p,,,sf], ln_theta_agg = ln_SrvLen_pop_theta_agg[p,sf],
                LN_corr_pars = SrvLen_pop_corr_pars[p,,,sf,], LN_corr_pars_agg = SrvLen_pop_corr_pars_agg[p,sf],
                n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = NA,
                use = UseSrvLenComps_pop[p,,y,seas,sf], n_model_bins = n_lens,
                n_obs_bins = dim(ObsSrvLenComps_pop)[5], addtocomp = addtocomp
              )
            }
          }
        }
      }
    }
  } else {

    # Doing OSA Compositions
    # Survey Age Compositions
    if(any(UseSrvAgeComps_pop == 1)) {

      # Discrete
      ObsSrvAgeComps_pop_osa_discrete = NULL
      ObsSrvAgeComps_pop_osa_discrete = pack_comp_osa(
        ObsArr = ObsSrvAgeComps_pop, ISSArr = ISS_SrvAgeComps_pop, WtArr = Wt_SrvAgeComps_pop,
        UseArr = UseSrvAgeComps_pop, TypeMat = SrvAgeComps_pop_Type, LikeTypeVec = SrvAgeComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsSrvAgeComps_pop_osa_discrete)) {
        ObsSrvAgeComps_pop_osa_discrete = RTMB::OBS(ObsSrvAgeComps_pop_osa_discrete)
        SrvAgeComps_pop_nLL = eval_comp_osa(
          nLL_arr = SrvAgeComps_pop_nLL, tracked = ObsSrvAgeComps_pop_osa_discrete,
          ExpArrFn = function(p, y, seas, f) { e = SrvIAA[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_ages, n_sexes); e },
          UseArr = UseSrvAgeComps_pop, TypeMat = SrvAgeComps_pop_Type, LikeTypeVec = SrvAgeComps_pop_LikeType,
          ISSArr = ISS_SrvAgeComps_pop, lnThetaArr = ln_SrvAge_pop_theta, lnThetaAggVec = ln_SrvAge_pop_theta_agg,
          LNcorrArr = SrvAge_pop_corr_pars, LNcorrAggVec = SrvAge_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsSrvAgeComps_pop)[5], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE, pop = TRUE, n_pop = n_pop
        )
      }

      # Continuous
      ObsSrvAgeComps_pop_osa_continuous = pack_comp_osa(
        ObsArr = ObsSrvAgeComps_pop, ISSArr = ISS_SrvAgeComps_pop, WtArr = Wt_SrvAgeComps_pop,
        UseArr = UseSrvAgeComps_pop, TypeMat = SrvAgeComps_pop_Type, LikeTypeVec = SrvAgeComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsSrvAgeComps_pop_osa_continuous)) {
        ObsSrvAgeComps_pop_osa_continuous = RTMB::OBS(ObsSrvAgeComps_pop_osa_continuous)
        SrvAgeComps_pop_nLL = eval_comp_osa(
          nLL_arr = SrvAgeComps_pop_nLL, tracked = ObsSrvAgeComps_pop_osa_continuous,
          ExpArrFn = function(p, y, seas, f) { e = SrvIAA[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_ages, n_sexes); e },
          UseArr = UseSrvAgeComps_pop, TypeMat = SrvAgeComps_pop_Type, LikeTypeVec = SrvAgeComps_pop_LikeType,
          ISSArr = ISS_SrvAgeComps_pop, lnThetaArr = ln_SrvAge_pop_theta, lnThetaAggVec = ln_SrvAge_pop_theta_agg,
          LNcorrArr = SrvAge_pop_corr_pars, LNcorrAggVec = SrvAge_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsSrvAgeComps_pop)[5], age_or_len = 0,
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsSrvAgeComps_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }

    # Survey Lengths
    if(fit_lengths == 1 && any(UseSrvLenComps_pop == 1)) {

      # Discrete
      ObsSrvLenComps_pop_osa_discrete = NULL
      ObsSrvLenComps_pop_osa_discrete = pack_comp_osa(
        ObsArr = ObsSrvLenComps_pop, ISSArr = ISS_SrvLenComps_pop, WtArr = Wt_SrvLenComps_pop,
        UseArr = UseSrvLenComps_pop, TypeMat = SrvLenComps_pop_Type, LikeTypeVec = SrvLenComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "discrete", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsSrvLenComps_pop_osa_discrete)) {
        ObsSrvLenComps_pop_osa_discrete = RTMB::OBS(ObsSrvLenComps_pop_osa_discrete)
        SrvLenComps_pop_nLL = eval_comp_osa(
          nLL_arr = SrvLenComps_pop_nLL, tracked = ObsSrvLenComps_pop_osa_discrete,
          ExpArrFn = function(p, y, seas, f) { e = SrvIAL[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_lens, n_sexes); e },
          UseArr = UseSrvLenComps_pop, TypeMat = SrvLenComps_pop_Type, LikeTypeVec = SrvLenComps_pop_LikeType,
          ISSArr = ISS_SrvLenComps_pop, lnThetaArr = ln_SrvLen_pop_theta, lnThetaAggVec = ln_SrvLen_pop_theta_agg,
          LNcorrArr = SrvLen_pop_corr_pars, LNcorrAggVec = SrvLen_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsSrvLenComps_pop)[5], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "discrete", zero_init = TRUE, pop = TRUE, n_pop = n_pop
        )
      }

      # Continuous
      ObsSrvLenComps_pop_osa_continuous = pack_comp_osa(
        ObsArr = ObsSrvLenComps_pop, ISSArr = ISS_SrvLenComps_pop, WtArr = Wt_SrvLenComps_pop,
        UseArr = UseSrvLenComps_pop, TypeMat = SrvLenComps_pop_Type, LikeTypeVec = SrvLenComps_pop_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
        addtocomp = addtocomp, family = "continuous", pop = TRUE, n_pop = n_pop
      )
      if(!is.null(ObsSrvLenComps_pop_osa_continuous)) {
        ObsSrvLenComps_pop_osa_continuous = RTMB::OBS(ObsSrvLenComps_pop_osa_continuous)
        SrvLenComps_pop_nLL = eval_comp_osa(
          nLL_arr = SrvLenComps_pop_nLL, tracked = ObsSrvLenComps_pop_osa_continuous,
          ExpArrFn = function(p, y, seas, f) { e = SrvIAL[p,,y,seas,,,f, drop=FALSE]; dim(e) = c(n_regions, n_lens, n_sexes); e },
          UseArr = UseSrvLenComps_pop, TypeMat = SrvLenComps_pop_Type, LikeTypeVec = SrvLenComps_pop_LikeType,
          ISSArr = ISS_SrvLenComps_pop, lnThetaArr = ln_SrvLen_pop_theta, lnThetaAggVec = ln_SrvLen_pop_theta_agg,
          LNcorrArr = SrvLen_pop_corr_pars, LNcorrAggVec = SrvLen_pop_corr_pars_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_lens, n_obs_bins = dim(ObsSrvLenComps_pop)[5], age_or_len = 1,
          AgeingErrorFn = NULL, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsSrvLenComps_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }
  }

  ## Tag Likelihoods ---------------------------------------------------------
  if(do_internal_conv_tag_osa == FALSE) {
    if(any(use_conv_fish_tagging == 1)) {

      conv_fish_tag_nLL <- get_conv_tag_likelihoods(
        n_conv_tag_cohorts          = n_conv_tag_cohorts,
        conv_tag_release_indicator  = conv_tag_release_indicator,
        conv_tag_max_liberty        = conv_tag_max_liberty,
        n_yrs                       = n_yrs,
        n_seas                      = n_seas,
        conv_tag_mixing_period      = conv_tag_mixing_period,
        n_fish_fleets               = n_fish_fleets,
        use_conv_fish_tagging       = use_conv_fish_tagging,
        n_conv_tag_pop_pool         = n_conv_tag_pop_pool,
        n_regions                   = n_regions,
        n_conv_tag_age_pool         = n_conv_tag_age_pool,
        n_conv_tag_sex_pool         = n_conv_tag_sex_pool,
        conv_tag_pop_pool           = conv_tag_pop_pool,
        conv_tag_age_pool           = conv_tag_age_pool,
        conv_tag_sex_pool           = conv_tag_sex_pool,
        conv_fish_tag_like          = conv_fish_tag_like,
        conv_fish_tag_nLL           = conv_fish_tag_nLL,
        obs_conv_tag_fish_recap     = obs_conv_tag_fish_recap,
        pred_conv_tag_fish_recap    = pred_conv_tag_fish_recap,
        addtotag                    = addtotag,
        ln_conv_fish_tag_theta      = ln_conv_fish_tag_theta,
        conv_tagged_fish            = conv_tagged_fish
      )

    } # if we are using tagging data
  } else {

    if(any(use_conv_fish_tagging == 1)) {

      tag_pack = pack_tag_osa(
        family                     = tag_fam_of(conv_fish_tag_like),
        like_type                  = conv_fish_tag_like,
        obs_recap                  = obs_conv_tag_fish_recap,
        pred_recap                 = pred_conv_tag_fish_recap,
        tagged_fish                = conv_tagged_fish,
        conv_tag_release_indicator = conv_tag_release_indicator,
        conv_tag_max_liberty       = conv_tag_max_liberty,
        n_conv_tag_cohorts         = n_conv_tag_cohorts,
        n_yrs                      = n_yrs,
        n_seas                     = n_seas,
        n_regions                  = n_regions,
        n_fish_fleets              = n_fish_fleets,
        n_pop_pool                 = n_conv_tag_pop_pool,
        n_age_pool                 = n_conv_tag_age_pool,
        n_sex_pool                 = n_conv_tag_sex_pool,
        pop_pool                   = conv_tag_pop_pool,
        age_pool                   = conv_tag_age_pool,
        sex_pool                   = conv_tag_sex_pool,
        use_fish_tagging           = use_conv_fish_tagging,
        conv_tag_mixing_period     = conv_tag_mixing_period,
        addtotag                   = addtotag
      )

      if(!is.null(tag_pack$vec)) {

        # Name the tracked object by family (used later in oneStepPredict)
        if(tag_fam_of(conv_fish_tag_like) == "count") {
          ObsConvTag_osa_count = tag_pack$vec
          ObsConvTag_osa_count = RTMB::OBS(ObsConvTag_osa_count)
        }
        if(tag_fam_of(conv_fish_tag_like) == "comp") {
          ObsConvTag_osa_comp = tag_pack$vec
          ObsConvTag_osa_comp = RTMB::OBS(ObsConvTag_osa_comp)
        }

        # Get tagging OSAs
        conv_fish_tag_nLL = eval_tag_osa(
          nLL_arr                    = conv_fish_tag_nLL,
          tracked                    = switch(tag_fam_of(conv_fish_tag_like), count = ObsConvTag_osa_count, comp  = ObsConvTag_osa_comp),
          family                     = tag_fam_of(conv_fish_tag_like),
          like_type                  = conv_fish_tag_like,
          pred_recap                 = pred_conv_tag_fish_recap,
          tagged_fish                = conv_tagged_fish,
          obs_recap                  = obs_conv_tag_fish_recap,
          conv_tag_release_indicator = conv_tag_release_indicator,
          conv_tag_max_liberty       = conv_tag_max_liberty,
          n_conv_tag_cohorts         = n_conv_tag_cohorts,
          n_yrs                      = n_yrs,
          n_seas                     = n_seas,
          n_regions                  = n_regions,
          n_fish_fleets              = n_fish_fleets,
          n_pop_pool                 = n_conv_tag_pop_pool,
          n_age_pool                 = n_conv_tag_age_pool,
          n_sex_pool                 = n_conv_tag_sex_pool,
          pop_pool                   = conv_tag_pop_pool,
          age_pool                   = conv_tag_age_pool,
          sex_pool                   = conv_tag_sex_pool,
          use_fish_tagging           = use_conv_fish_tagging,
          conv_tag_mixing_period     = conv_tag_mixing_period,
          addtotag                   = addtotag,
          ln_theta                   = ln_conv_fish_tag_theta,
          zero_init                  = TRUE
        )

      }
    }

  }

  ## Priors and Penalties ----------------------------------------------------
  ### Fishing Mortality (Penalty) ---------------------------------------------
  if(Use_F_pen == 1) {
    Fmort_nLL = Get_Fdev_PE_loglik(PE_model = Fdev_model, ln_sigmaF = ln_sigmaF, Fdev_rho = Fdev_rho,
                                   ln_F_devs = ln_F_devs, UseCatch = UseCatch, UseCatch_pop = UseCatch_pop,
                                   missing_catch = missing_catch)
  } #  if using fishing mortality penalty

  ### Discard Mortality Rate (Penalty) ---------------------------------------------
  if(Use_dmr_pen == 1) {
    for(f in 1:n_fish_fleets) {
      for(y in 1:n_yrs) {
        for(r in 1:n_regions) {
          for(seas in 1:n_seas) {

            if(UseDiscard[r,y,seas,f] == 1 || any(UseDiscard_pop[,r,y,seas,f] == 1)) {
              dmr_nLL[r,y,seas,f] = -RTMB::dnorm(logit_dmr_devs[r,y,seas,f], 0, exp(ln_sigma_dmr[r,seas,f]), TRUE)
            } # end if have catch

          } # end seas loop
        } # end r loop
      } # y loop
    } # f loop
  } #  if using discard mortality rate penalty


  ### Selectivity (Penalty) ---------------------------------------------------
  for(r in 1:n_regions) {

    for(f in 1:n_fish_fleets) {

      # Total Fishery Selectivity Deviations
      if(cont_tv_fish_sel[r,f] > 0) {

        sel_nLL = sel_nLL + - Get_sel_PE_loglik(PE_model = cont_tv_fish_sel[r,f], # process error model
                                                PE_pars = fishsel_pe_pars[r,,,f, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_fishsel_devs[r,,,,f, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_fishsel_devs[r,,,,f, drop = FALSE],
                                                min_sel_devs_shared_bins = fishsel_devs_min_shared_bins

        )
      } # end if

      # Retained Fishery Selectivity Deviations
      if(cont_tv_ret_sel[r,f] > 0) {

        sel_nLL = sel_nLL + - Get_sel_PE_loglik(PE_model = cont_tv_ret_sel[r,f], # process error model
                                                PE_pars = retsel_pe_pars[r,,,f, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_retsel_devs[r,,,,f, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_retsel_devs[r,,,,f, drop = FALSE],
                                                min_sel_devs_shared_bins = retsel_devs_min_shared_bins

        )
      } # end if
    } # end f loop

    # Survey Selectivity Deviations
    for(sf in 1:n_srv_fleets) {

      if(cont_tv_srv_sel[r,sf] > 0) {

        sel_nLL = sel_nLL + - Get_sel_PE_loglik(PE_model = cont_tv_srv_sel[r,sf], # process error model
                                                PE_pars = srvsel_pe_pars[r,,,sf, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_srvsel_devs[r,,,,sf, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_srvsel_devs[r,,,,sf, drop = FALSE],
                                                min_sel_devs_shared_bins = srvsel_devs_min_shared_bins

        )
      } # end if

    } # end sf loop
  } # end r loop


  ### Selectivity Smoothness (Penalty) --------------------------------------------------
  smooth_pen_terms = c("smooth_bin_curve", "smooth_bin_diff", "smooth_yr_diff", "smooth_yr_curve", "smooth_dome", "smooth_mean_center")
  for(r in 1:n_regions) {

    for(f in 1:n_fish_fleets) {

      # If Bicubic spline
      bicubic_yrs = which(fish_sel_model[r,,f] == 8)
      has_bicubic = length(bicubic_yrs) > 0
      has_nonzero_pen = any(sapply(smooth_pen_terms, function(nm) safe_extract(fish_sel_pen_wts, nm)) != 0)

      if(has_bicubic || has_nonzero_pen) {
        if(has_bicubic) {
          block_yrs = min(bicubic_yrs):max(bicubic_yrs)
          # Restrict to the actual fit range and bins from the penalty
          selstyr_this = unique(fish_sel_bicubic_selstyr[r, block_yrs, f])
          y_range = if(selstyr_this == 0) block_yrs else block_yrs[block_yrs >= which(data$years == selstyr_this)]
          nselbins_this = unique(fish_sel_bicubic_nselbins[r, block_yrs, f])
          n_fit_bins = if(nselbins_this == 0) (if(fish_selex_type == 0) n_ages else dim(fish_sel_l)[3]) else nselbins_this
        } else {
          # non-bicubic fleet: no sub-range restriction, use the fleet's whole modeled history
          y_range = 1:n_yrs
          n_fit_bins = if(fish_selex_type == 0) n_ages else dim(fish_sel_l)[3]
        }

        # get sel values
        if(fish_selex_type == 0) tmp_sel_vals = array(fish_sel[1,r,y_range,1,1:n_fit_bins,,f, drop = FALSE], dim = c(1, length(y_range), n_fit_bins, n_sexes, 1))
        if(fish_selex_type == 1) tmp_sel_vals = array(fish_sel_l[r,y_range,1:n_fit_bins,,f, drop = FALSE], dim = c(1, length(y_range), n_fit_bins, n_sexes, 1))

        sel_nLL = sel_nLL - Get_Selex_Smoothness_Penalty(tmp_sel_vals,
                                                         wt_bin_curve = safe_extract(fish_sel_pen_wts, "smooth_bin_curve"),
                                                         wt_bin_diff = safe_extract(fish_sel_pen_wts, "smooth_bin_diff"),
                                                         wt_yr_diff = safe_extract(fish_sel_pen_wts, "smooth_yr_diff"),
                                                         wt_yr_curve = safe_extract(fish_sel_pen_wts, "smooth_yr_curve"),
                                                         wt_dome = safe_extract(fish_sel_pen_wts, "smooth_dome"),
                                                         wt_mean_center = safe_extract(fish_sel_pen_wts, "smooth_mean_center"),
                                                         normalize = TRUE)
      } # end if
    } # end f loop

    for(f in 1:n_fish_fleets) {

      # If Bicubic spline (retention)
      bicubic_yrs = which(ret_sel_model[r,,f] == 8)
      has_bicubic = length(bicubic_yrs) > 0
      has_nonzero_pen = any(sapply(smooth_pen_terms, function(nm) safe_extract(ret_sel_pen_wts, nm)) != 0)

      if(has_bicubic || has_nonzero_pen) {
        if(has_bicubic) {
          block_yrs = min(bicubic_yrs):max(bicubic_yrs)
          # Restrict to the actual fit range and bins from the penalty
          selstyr_this = unique(ret_sel_bicubic_selstyr[r, block_yrs, f])
          y_range = if(selstyr_this == 0) block_yrs else block_yrs[block_yrs >= which(data$years == selstyr_this)]
          nselbins_this = unique(ret_sel_bicubic_nselbins[r, block_yrs, f])
          n_fit_bins = if(nselbins_this == 0) (if(ret_selex_type == 0) n_ages else dim(ret_sel_l)[3]) else nselbins_this
        } else {
          # non-bicubic fleet: no sub-range restriction, use the fleet's whole modeled history
          y_range = 1:n_yrs
          n_fit_bins = if(ret_selex_type == 0) n_ages else dim(ret_sel_l)[3]
        }

        # get sel values
        if(ret_selex_type == 0) tmp_sel_vals = array(ret_sel[1,r,y_range,1,1:n_fit_bins,,f, drop = FALSE], dim = c(1, length(y_range), n_fit_bins, n_sexes, 1))
        if(ret_selex_type == 1) tmp_sel_vals = array(ret_sel_l[r,y_range,1:n_fit_bins,,f, drop = FALSE], dim = c(1, length(y_range), n_fit_bins, n_sexes, 1))

        sel_nLL = sel_nLL - Get_Selex_Smoothness_Penalty(tmp_sel_vals,
                                                         wt_bin_curve = safe_extract(ret_sel_pen_wts, "smooth_bin_curve"),
                                                         wt_bin_diff = safe_extract(ret_sel_pen_wts, "smooth_bin_diff"),
                                                         wt_yr_diff = safe_extract(ret_sel_pen_wts, "smooth_yr_diff"),
                                                         wt_yr_curve = safe_extract(ret_sel_pen_wts, "smooth_yr_curve"),
                                                         wt_dome = safe_extract(ret_sel_pen_wts, "smooth_dome"),
                                                         wt_mean_center = safe_extract(ret_sel_pen_wts, "smooth_mean_center"),
                                                         normalize = TRUE)
      } # end if
    } # end f loop

    for(sf in 1:n_srv_fleets) {

      # if bicubic
      bicubic_yrs = which(srv_sel_model[r,,sf] == 8)
      has_bicubic = length(bicubic_yrs) > 0
      has_nonzero_pen = any(sapply(smooth_pen_terms, function(nm) safe_extract(srv_sel_pen_wts, nm)) != 0)
      if(has_bicubic || has_nonzero_pen) {
        if(has_bicubic) {
          block_yrs = min(bicubic_yrs):max(bicubic_yrs)
          selstyr_this = unique(srv_sel_bicubic_selstyr[r, block_yrs, sf])
          y_range = if(selstyr_this == 0) block_yrs else block_yrs[block_yrs >= which(data$years == selstyr_this)]
          nselbins_this = unique(srv_sel_bicubic_nselbins[r, block_yrs, sf])
          n_fit_bins = if(nselbins_this == 0) (if(srv_selex_type == 0) n_ages else dim(srv_sel_l)[3]) else nselbins_this
        } else {
          # non-bicubic fleet: no sub-range restriction, use the fleet's whole modeled history
          y_range = 1:n_yrs
          n_fit_bins = if(srv_selex_type == 0) n_ages else dim(srv_sel_l)[3]
        }

        # get sel values
        if(srv_selex_type == 0) tmp_sel_vals = array(srv_sel[1,r,y_range,1,1:n_fit_bins,,sf, drop = FALSE], dim = c(1, length(y_range), n_fit_bins, n_sexes, 1))
        if(srv_selex_type == 1) tmp_sel_vals = array(srv_sel_l[r,y_range,1:n_fit_bins,,sf, drop = FALSE], dim = c(1, length(y_range), n_fit_bins, n_sexes, 1))

        sel_nLL = sel_nLL - Get_Selex_Smoothness_Penalty(tmp_sel_vals,
                                                         wt_bin_curve = safe_extract(srv_sel_pen_wts, "smooth_bin_curve"),
                                                         wt_bin_diff = safe_extract(srv_sel_pen_wts, "smooth_bin_diff"),
                                                         wt_yr_diff = safe_extract(srv_sel_pen_wts, "smooth_yr_diff"),
                                                         wt_yr_curve = safe_extract(srv_sel_pen_wts, "smooth_yr_curve"),
                                                         wt_dome = safe_extract(srv_sel_pen_wts, "smooth_dome"),
                                                         wt_mean_center = safe_extract(srv_sel_pen_wts, "smooth_mean_center"),
                                                         normalize = TRUE)
      } # end if
    } # end sf loop

  } # end r loop


  ### Selectivity (Prior) -----------------------------------------------------
  # Total Fishery selectivity parameters
  if(Use_fish_selex_prior == 1) {
    for(i in 1:nrow(fish_selex_prior)) {
      # Extract indices
      r = fish_selex_prior$region[i]
      p = fish_selex_prior$par[i]
      b = fish_selex_prior$block[i]
      s = fish_selex_prior$sex[i]
      f = fish_selex_prior$fleet[i]
      # Compute penalty / prior here
      sel_nLL = sel_nLL - RTMB::dnorm(fish_fixed_sel_pars[r,p,b,s,f], log(fish_selex_prior$mu[i]), fish_selex_prior$sd[i], TRUE)
    } # end i loop
  } # end if using selex priors

  # Retained Fishery selectivity parameters
  if(Use_ret_selex_prior == 1) {
    for(i in 1:nrow(ret_selex_prior)) {
      # Extract indices
      r = ret_selex_prior$region[i]
      p = ret_selex_prior$par[i]
      b = ret_selex_prior$block[i]
      s = ret_selex_prior$sex[i]
      f = ret_selex_prior$fleet[i]
      # Compute penalty / prior here
      sel_nLL = sel_nLL - RTMB::dnorm(ret_fixed_sel_pars[r,p,b,s,f], log(ret_selex_prior$mu[i]), ret_selex_prior$sd[i], TRUE)
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
      sel_nLL = sel_nLL - RTMB::dnorm(srv_fixed_sel_pars[r,p,b,s,sf], log(srv_selex_prior$mu[i]), srv_selex_prior$sd[i], TRUE)
    } # end i loop
  } # end if using selex priors

  ### Recruitment (Penalty) ----------------------------------------------------
  for(p in 1:n_pop) {
    for(r in 1:n_regions) {

      # get sigma index
      sigma_idx = ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])

      # Skip penalty if no dispersal and p = r has no recruits
      if(rec_region_prop_spec == 1 && as.numeric(rec_region_prop[p,r]) == 0) next

      # Initial age deviations (if equil_init_age_strc == 0; don't penalize at all)
      if(equil_init_age_strc == 1) Init_Rec_nLL[p,r,1:(n_ages - 2)] = -RTMB::dnorm(ln_InitDevs[p,r,1:(n_ages - 2)], -exp(ln_sigmaR[1,p,sigma_idx])^2/2 * bias_ramp[1], exp(ln_sigmaR[1,p,sigma_idx]), TRUE) # only penalize non plus group
      if(equil_init_age_strc == 2) Init_Rec_nLL[p,r,] = -RTMB::dnorm(ln_InitDevs[p,r,], -exp(ln_sigmaR[1,p,sigma_idx])^2/2 * bias_ramp[1], exp(ln_sigmaR[1,p,sigma_idx]), TRUE) # penalize all
      if(equil_init_age_strc == 3) Init_Rec_nLL[p,r,unique(init_age_devs_shared)] = -RTMB::dnorm(ln_InitDevs[p,r,unique(init_age_devs_shared)], -exp(ln_sigmaR[1,p,sigma_idx])^2/2 * bias_ramp[1], exp(ln_sigmaR[1,p,sigma_idx]), TRUE) # penalize all

      # Early recruitment deviations
      if(sigmaR_switch > 1) {
        Rec_nLL[p,r,1:(sigmaR_switch-1)] = -RTMB::dnorm(ln_RecDevs[p,r,1:(sigmaR_switch-1)], -sigmaR2_early[p,sigma_idx]/2 * bias_ramp[1:(sigmaR_switch-1)], exp(ln_sigmaR[1,p,sigma_idx]), TRUE)
        if(do_rec_bias_ramp == 1 && any(bias_ramp != 0)) Rec_nLL[p,r,1:(sigmaR_switch-1)] = Rec_nLL[p,r,1:(sigmaR_switch-1)] - (1 - 0.5 * bias_ramp[1:(sigmaR_switch-1)]) * ln_sigmaR[1,p,sigma_idx] # adjust w/ bias correction
      }

      # Late recruitment deviations
      Rec_nLL[p,r,sigmaR_switch:n_est_rec_devs] = -RTMB::dnorm(ln_RecDevs[p,r,sigmaR_switch:n_est_rec_devs], -sigmaR2_late[p,sigma_idx]/2 * bias_ramp[sigmaR_switch:n_est_rec_devs], exp(ln_sigmaR[2,p,sigma_idx]), TRUE)
      if(do_rec_bias_ramp == 1 && any(bias_ramp != 0)) Rec_nLL[p,r,sigmaR_switch:n_est_rec_devs] = Rec_nLL[p,r,sigmaR_switch:n_est_rec_devs] - (1 - 0.5 * bias_ramp[sigmaR_switch:n_est_rec_devs]) * ln_sigmaR[2,p,sigma_idx] # adjust w/ bias correction

    } # end r loop
  } # end p loop

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
      p = M_prior$popblk[i]
      r = M_prior$regionblk[i]
      b = M_prior$yearblk[i]
      a = M_prior$ageblk[i]
      s = M_prior$sexblk[i]
      idx = M_blocks[p,r,y,a,s]
      # Compute prior
      M_nLL = M_nLL + -RTMB::dnorm(ln_M[idx], log(M_prior$mu[i]), M_prior$sd[i], TRUE) # TMB likelihood
    }
  } # end if using natural mortality prior

  ### Steepness (Prior) -----------------------------------------------
  if(Use_h_prior == 1) {
    for(i in 1:nrow(h_prior)) {
      # Extract indicies
      p = h_prior$pop[i]
      r = h_prior$region[i]
      tmp_h_beta_pars = get_beta_scaled_pars(low = 0.2, high = 1, mu = h_prior$mu[i], sigma = h_prior$sd[i]) # get alpha and beta parameters
      tmp_h_trans = (h_trans[p,r] - 0.2) / (1 - 0.2) # transform random variable
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
      p = Movement_prior$pop[i] # population
      region_from = Movement_prior$region_from[i] # region from
      y = Movement_prior$year[i] # year
      seas = Movement_prior$seas[i] # seas
      a = Movement_prior$age[i] # age
      s = Movement_prior$sex[i] # sex
      alpha = Movement_prior$alpha[[i]] # get prior values
      Movement_nLL = Movement_nLL - ddirichlet(x = Movement[p, region_from,,y,seas,a,s], alpha = alpha, log = TRUE) # dirichlet prior
    } # end i loop
  }

  ### Recruitment R0 and Proportions (Prior) -----------------------------------------
  if(use_rec_region_prop_prior == 1) { # recruitment regional apportionment
    for(i in 1:nrow(rec_region_prop_prior)) {
      p = rec_region_prop_prior$pop[i] # population
      alpha = rec_region_prop_prior$alpha[[i]] # get concentration values
      rec_prop_nLL = -ddirichlet(x = rec_region_prop[p,], alpha = alpha, log = TRUE) # dirichlet prior
    }
  }

  if(use_rec_seas_prop_prior == 1 && use_fixed_rec_seas_prop == 0) { # recruitment seasonal apportionment
    for(i in 1:nrow(rec_seas_prop_prior)) { # recruitment seasonal apportionment
      p = rec_seas_prop_prior$pop[i] # population
      alpha = rec_seas_prop_prior$alpha[[i]] # get concentration values
      if(rec_lag == 0 && spawn_seas > 1) {
        # seasons before spawn_seas are structurally zero (age-0 recruits
        # can't predate the spawning event that produced them); not evaluating then ...
        rec_prop_nLL = rec_prop_nLL -ddirichlet(x = rec_seas_prop[p, spawn_seas:n_seas], alpha = alpha, log = TRUE) # dirichlet prior
      } else {
        rec_prop_nLL = rec_prop_nLL -ddirichlet(x = rec_seas_prop[p,], alpha = alpha, log = TRUE) # dirichlet prior
      }
    }
  }

  if(use_r0_prior == 1) { # recruitment R0
    for(i in 1:nrow(r0_prior)) {
      p = r0_prior$pop[i] # population
      rec_nLL = rec_nLL - RTMB::dnorm(ln_global_R0[p], log(r0_prior$mu[i]), r0_prior$sd[i], TRUE) # normal prior
    }
  }

  ### Stray Rates (Prior) -----------------------------------------------------
  if(use_stray_rate_prior == 1) {
    for(i in 1:nrow(stray_rate_prior)) {
      # extract indices
      p = stray_rate_prior$pop[i]
      b = stray_rate_prior$block[i]
      # extract beta pars
      mu = stray_rate_prior$mu[i]
      sd = stray_rate_prior$sd[i]
      # derive beta pars
      concentration = mu * (1 - mu) / sd^2 - 1
      alpha = mu * concentration
      beta  = (1 - mu) * concentration
      # extract values
      stray_rate_val = 1e-4 + (1 - 2*1e-4) * RTMB::plogis(stray_rate_pars[p,b])
      rec_prop_nLL = rec_prop_nLL -RTMB::dbeta(x = stray_rate_val, shape1 = alpha, shape2 = beta, log = TRUE) # penalize
    }
  }

  ### Tag Reporting Rate (Prior) --------------------------------------------
  if(use_conv_tag_fishrep_prior == 1) {
    for(i in 1:nrow(conv_tag_fishrep_prior)) {

      # Extract indices
      r = conv_tag_fishrep_prior$region[i]
      b = conv_tag_fishrep_prior$block[i]
      f = conv_tag_fishrep_prior$fleet[i]

      conv_tag_fishrep_val = RTMB::plogis(conv_tag_fish_reporting_pars[r,b,f]) # extract tag reporting rate value
      if(conv_tag_fishrep_prior$type[i] == 0) {
        TagRep_nLL = TagRep_nLL - dbeta_symmetric(p_val = conv_tag_fishrep_val, p_ub = 1, p_lb = 0, p_prsd = conv_tag_fishrep_prior$sd[i], log = TRUE) # penalize
      } # end if symmetric beta

      if(conv_tag_fishrep_prior$type[i] == 1) {
        # extract pars
        mu = conv_tag_fishrep_prior$mu[i]
        sd = conv_tag_fishrep_prior$sd[i]
        # derive beta pars
        concentration = mu * (1 - mu) / sd^2 - 1
        alpha = mu * concentration
        beta  = (1 - mu) * concentration
        TagRep_nLL = TagRep_nLL -RTMB::dbeta(x = conv_tag_fishrep_val, shape1 = alpha, shape2 = beta, log = TRUE) # penalize
      } # end if for full beta

    } # end i loop
  } # if use tag reporting prior

  # Sum up nLL
  jnLL = sum(Wt_Catch * Catch_nLL) +             # Aggregated catch likelihoods
    sum(Wt_Catch_pop * Catch_pop_nLL) +      # Pop-specific catch likelihoods
    sum(Wt_Discard * Discard_nLL) +           # Aggregated discard likelihoods
    sum(Wt_Discard_pop * Discard_pop_nLL) +   # Pop-specific discard likelihoods
    sum(Wt_FishIdx * FishIdx_nLL) +           # Aggregated fishery index likelihoods
    sum(Wt_FishIdx_pop * FishIdx_pop_nLL) +   # Pop-specific fishery index likelihoods
    sum(Wt_SrvIdx * SrvIdx_nLL) +             # Aggregated survey index likelihoods
    sum(Wt_SrvIdx_pop * SrvIdx_pop_nLL) +     # Pop-specific survey index likelihoods
    sum(FishAgeComps_nLL) +                   # Aggregated fishery age likelihoods
    sum(FishAgeComps_pop_nLL) +               # Pop-specific fishery age likelihoods
    sum(FishLenComps_nLL) +                   # Aggregated fishery length likelihoods
    sum(FishLenComps_pop_nLL) +               # Pop-specific fishery length likelihoods
    sum(FishAgeComps_discard_nLL) +            # Aggregated discard age likelihoods
    sum(FishAgeComps_discard_pop_nLL) +        # Pop-specific discard age likelihoods
    sum(FishLenComps_discard_nLL) +            # Aggregated discard length likelihoods
    sum(FishLenComps_discard_pop_nLL) +        # Pop-specific discard length likelihoods
    sum(SrvAgeComps_nLL) +                    # Aggregated survey age likelihoods
    sum(SrvAgeComps_pop_nLL) +                # Pop-specific survey age likelihoods
    sum(SrvLenComps_nLL) +                    # Aggregated survey length likelihoods
    sum(SrvLenComps_pop_nLL) +                # Pop-specific survey length likelihoods
    (Wt_Tagging * sum(conv_fish_tag_nLL)) +   # Tagging likelihood
    (Wt_F * sum(Fmort_nLL)) +                 # Fishing mortality penalty
    (Wt_D * sum(dmr_nLL)) +                   # Discard mortality rate penalty
    (Wt_Rec * sum(Rec_nLL)) +                 # Recruitment penalty
    (Wt_Rec * sum(Init_Rec_nLL)) +            # Initial age penalty
    sel_nLL +                                  # Selectivity penalty
    M_nLL +                                    # Natural mortality prior
    h_nLL +                                    # Steepness prior
    Movement_nLL +                             # Movement prior
    TagRep_nLL +                               # Tag reporting rate prior
    fish_q_nLL +                               # Fishery q prior
    srv_q_nLL +                                # Survey q prior
    rec_prop_nLL                               # Recruitment proportion prior

  # Report Section ----------------------------------------------------------
  # Biological Processes
  RTMB::REPORT(R0)
  RTMB::REPORT(rinit)
  RTMB::REPORT(sexratio)
  RTMB::REPORT(rec_region_prop)
  RTMB::REPORT(rec_seas_prop)
  RTMB::REPORT(stray_rate)
  RTMB::REPORT(h_trans)
  RTMB::REPORT(NAA)
  RTMB::REPORT(NAA0)
  RTMB::REPORT(NAA_bef)
  RTMB::REPORT(NAA_aft)
  RTMB::REPORT(ZAA)
  RTMB::REPORT(natmort)
  RTMB::REPORT(bias_ramp)
  RTMB::REPORT(Movement)
  RTMB::REPORT(sgl_seas_spawning_movement)
  RTMB::REPORT(Mrate)

  # Fishery Processes
  RTMB::REPORT(init_F)
  RTMB::REPORT(ln_sigmaC)
  RTMB::REPORT(ln_sigmaC_pop)
  RTMB::REPORT(Fmort)
  RTMB::REPORT(dmr)
  RTMB::REPORT(tot_FAA)
  RTMB::REPORT(ret_FAA)
  RTMB::REPORT(disc_FAA)
  RTMB::REPORT(CAA)
  RTMB::REPORT(DAA)
  RTMB::REPORT(CAL)
  RTMB::REPORT(DAL)
  RTMB::REPORT(PredCatch)
  RTMB::REPORT(PredDiscard)
  RTMB::REPORT(PredFishIdx)
  RTMB::REPORT(fish_sel)
  RTMB::REPORT(ret_sel)
  RTMB::REPORT(fish_q)

  # Survey Processes
  RTMB::REPORT(PredSrvIdx)
  RTMB::REPORT(srv_sel)
  RTMB::REPORT(srv_q)
  RTMB::REPORT(SrvIAA)
  RTMB::REPORT(SrvIAL)

  # Report length-based selectivity
  if(fish_selex_type == 1) RTMB::REPORT(fish_sel_l)
  if(ret_selex_type == 1) RTMB::REPORT(ret_sel_l)
  if(srv_selex_type == 1) RTMB::REPORT(srv_sel_l)

  # Tagging Processes
  if(any(use_conv_fish_tagging == 1)) {
    RTMB::REPORT(pred_conv_tag_fish_recap)
    RTMB::REPORT(conv_tag_fish_avail)
    RTMB::REPORT(conv_tag_fish_reporting)
  }

  # Parameter Deviations
  RTMB::REPORT(ln_RecDevs)
  RTMB::REPORT(move_devs)
  RTMB::REPORT(ln_fishsel_devs)
  RTMB::REPORT(ln_srvsel_devs)

  # Aggregated Likelihoods
  RTMB::REPORT(Catch_nLL)
  RTMB::REPORT(Discard_nLL)
  RTMB::REPORT(FishIdx_nLL)
  RTMB::REPORT(SrvIdx_nLL)
  RTMB::REPORT(FishAgeComps_nLL)
  RTMB::REPORT(FishAgeComps_discard_nLL)
  RTMB::REPORT(SrvAgeComps_nLL)
  RTMB::REPORT(FishLenComps_nLL)
  RTMB::REPORT(FishLenComps_discard_nLL)
  RTMB::REPORT(SrvLenComps_nLL)

  # Population-specific Likelihoods
  RTMB::REPORT(Catch_pop_nLL)
  RTMB::REPORT(Discard_pop_nLL)
  RTMB::REPORT(FishIdx_pop_nLL)
  RTMB::REPORT(SrvIdx_pop_nLL)
  RTMB::REPORT(FishAgeComps_pop_nLL)
  RTMB::REPORT(FishAgeComps_discard_pop_nLL)
  RTMB::REPORT(SrvAgeComps_pop_nLL)
  RTMB::REPORT(FishLenComps_pop_nLL)
  RTMB::REPORT(FishLenComps_discard_pop_nLL)
  RTMB::REPORT(SrvLenComps_pop_nLL)

  # Penalties and priors
  RTMB::REPORT(M_nLL)
  RTMB::REPORT(Fmort_nLL)
  RTMB::REPORT(dmr_nLL)
  RTMB::REPORT(Rec_nLL)
  RTMB::REPORT(Init_Rec_nLL)
  RTMB::REPORT(conv_fish_tag_nLL)
  RTMB::REPORT(h_nLL)
  RTMB::REPORT(fish_q_nLL)
  RTMB::REPORT(sel_nLL)
  RTMB::REPORT(srv_q_nLL)
  RTMB::REPORT(Movement_nLL)
  RTMB::REPORT(TagRep_nLL)
  RTMB::REPORT(rec_prop_nLL)
  RTMB::REPORT(jnLL)

  # Report for derived quantities
  RTMB::REPORT(Total_Biom)
  RTMB::REPORT(SSB)
  RTMB::REPORT(eff_SSB)
  RTMB::REPORT(Dynamic_SSB0)
  RTMB::REPORT(Aggregated_SSB)
  RTMB::REPORT(Dynamic_Aggregated_SSB0)
  RTMB::REPORT(Rec)

  # Report these in log space and add constant because can't be < 0 or == 0
  log_Total_Biom = log(Total_Biom + 1e-5)
  log_SSB = log(SSB + 1e-5)
  log_eff_SSB = log(eff_SSB + 1e-5)
  log_Dynamic_SSB0 = log(Dynamic_SSB0 + 1e-5)
  log_Rec = log(Rec + 1e-5)
  log_Aggregated_SSB = log(Aggregated_SSB + 1e-5)
  log_Dynamic_Aggregated_SSB0 = log(Dynamic_Aggregated_SSB0 + 1e-5)

  RTMB::ADREPORT(log_Total_Biom)
  RTMB::ADREPORT(log_SSB)
  RTMB::ADREPORT(log_eff_SSB)
  RTMB::ADREPORT(log_Dynamic_SSB0)
  RTMB::ADREPORT(log_Rec)
  RTMB::ADREPORT(log_Aggregated_SSB)
  RTMB::ADREPORT(log_Dynamic_Aggregated_SSB0)

  return(jnLL)
} # end function
