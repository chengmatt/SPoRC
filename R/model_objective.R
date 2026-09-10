# Stage 2 of 3: objective function
#
# The one RTMB objective function, evaluated by fit_model. It runs end to end and is regulated by
# the switches set during setup; the array arithmetic lives in the model_*.R modules it calls.

# Development history of the objective function, kept as provenance for the
# assessment. Current behavior is documented at each section below, not here.

# Stages, in the order they run:
#   1 Model Set Up (Containers)   containers, and the dim shapes named once
#   2 Parameter Transformations   movement, M, growth, selectivity, recruitment, mortality at age
#   3 Population Dynamics         initial age structure, cohort growth, projection, numbers at age state
#   4 Observation Models          fishery, survey and tagging predictions
#   5 Likelihood Equations        catch, discards, indices, compositions, age at length, tags
#   6 Priors and Penalties        deviation penalties, selectivity, recruitment and growth, priors
#   7 Joint Negative Log Likelihood  every term above, weighted and summed
#   8 Report Section                 REPORT and ADREPORT

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
# Natural mortality split out into its own module
# Removed ADMB likelihoods and Sablefish specific calculations
# Added equilibrium plus group calculations to initial abundance when movement occurs
# Initial numbers at age split out into its own module

# version 5 - (M.LH Cheng & J.T Thorson)
# Reworked movement priors and movement setup
# Added in CTMC movement

# verison 6 - (M.LH Cheng)
# Re-coded tagging module to include fleet-specific dynamics
# Expanded model to include populations and season partitions (natal homing
# and seasonality dynamics)
# Incorporated functionality to allow movement to be a continuous process, in addition to movement after mortality

#' Generalized RTMB spatial age-structured model
#'
#' @param pars Parameter List
#' @param data Data List
#' @import RTMB
#' @keywords internal
SPoRC_rtmb = function(pars, data) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data, warn = FALSE) # load in starting values and data

  maintain_backwards_compatibility() # defaults for input lists built by older SPoRC versions

  # Model Set Up (Containers) -----------------------------------------------
  # a length composition recorded on coarser bins than the model tracks is mapped through LenBinMap.
  # that is the operation ageing error does on ages, so it is read by year and fleet to match
  LenBinMap_fn = if(is.null(LenBinMap)) function(y, f) NA else function(y, f) LenBinMap # NA leaves the bins alone

  n_ages = length(ages) # number of ages
  n_yrs = length(years) # number of years
  n_lens = length(lens) # number of lengths

  # container shapes, named once so a dim ordering is written down in one place.
  dim_naa = c(n_pop, n_regions, n_yrs + 1, n_seas, n_ages, n_sexes) # numbers at age, one extra year
  dim_mort_age = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes) # mortality rates at age
  dim_pop_reg_yr = c(n_pop, n_regions, n_yrs) # a series by population and region
  dim_pop_yr = c(n_pop, n_yrs) # the same, aggregated over regions
  dim_pop_reg = c(n_pop, n_regions) # one value per population and region

  # fishery shapes
  dim_pop_fish = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets) # one value per fleet
  dim_fish = c(n_regions, n_yrs, n_seas, n_fish_fleets)
  dim_pop_fish_age = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets) # at age
  dim_fish_age = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_fish_fleets)
  dim_pop_fish_len = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets) # at length
  dim_pop_fish_caal = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_fish_fleets) # at length and age
  dim_fish_caal = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets) # one value per length bin
  dim_pop_fish_comp = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets) # one value per sex
  dim_fish_comp = c(n_regions, n_yrs, n_seas, n_sexes, n_fish_fleets)
  dim_fish_yr = c(n_regions, n_yrs, n_fish_fleets) # annual, no season

  # survey shapes, the same layout with survey fleets
  dim_pop_srv = c(n_pop, n_regions, n_yrs, n_seas, n_srv_fleets)
  dim_srv = c(n_regions, n_yrs, n_seas, n_srv_fleets)
  dim_pop_srv_age = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)
  dim_srv_age = c(n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)
  dim_pop_srv_len = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets)
  dim_pop_srv_caal = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_srv_fleets)
  dim_srv_caal = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets)
  dim_pop_srv_comp = c(n_pop, n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)
  dim_srv_comp = c(n_regions, n_yrs, n_seas, n_sexes, n_srv_fleets)
  dim_srv_yr = c(n_regions, n_yrs, n_srv_fleets)

  # Recruitment stuff
  n_est_rec_devs = dim(ln_RecDevs)[3] # number of recruitment deviates estimated
  Rec = array(0, dim = dim_pop_reg_yr) # Recruitment
  R0 = array(0, dim = dim_pop_reg) # R0 or mean recruitment
  sexratio = array(0, dim = c(n_pop, n_regions, n_yrs, n_sexes)) # recruitment sex ratio
  rec_region_prop = array(0, dim = dim_pop_reg) # recruitment regional apportionment
  rec_seas_prop = array(0, dim = c(n_pop, n_seas)) # recruitment seasonal apportionment
  stray_rate = array(0, dim = dim_pop_yr) # stray rate

  # Population Dynamics
  NAA = array(data = 0, dim = dim_naa) # Numbers at age
  NAA_bef = array(data = 0, dim = dim_naa) # Numbers at age before movement
  NAA_aft = array(data = 0, dim = dim_naa) # Numbers at age after movement
  NAA0 = array(data = 0, dim = dim_naa) # Unfished Numbers at age
  ZAA = array(data = 0, dim = dim_mort_age) # Total mortality at age
  natmort = array(data = 0, dim = dim_mort_age) # natural mortality at age, a rate per year in each season
  Total_Biom = array(0, dim = dim_pop_reg_yr) # Total biomass
  SSB = array(0, dim = dim_pop_reg_yr) # Spawning stock biomass
  eff_SSB = array(0, dim = dim_pop_yr) # effective SSB for a given population
  Dynamic_SSB0 = array(0, dim = dim_pop_reg_yr) # Dynamic Unfished Spawning stock biomass
  Aggregated_SSB = array(0, dim = dim_pop_yr) # Aggregated Spawning stock biomass
  Dynamic_Aggregated_SSB0 = array(0, dim = dim_pop_yr) # Dynamic Unfished Aggregated Spawning stock biomass

  # Movement Stuff
  Movement = array(data = 0, dim = c(n_pop, n_regions, n_regions, n_yrs + n_proj_yrs_devs, n_seas, n_ages, n_sexes)) # movement "matrix"
  n_conv_tag_pop_pool = length(conv_tag_pop_pool) # number of populations to pool for tagging data
  n_conv_tag_age_pool = length(conv_tag_age_pool) # number of ages to pool for tagging data
  n_conv_tag_sex_pool = length(conv_tag_sex_pool) # number of sexes to pool for tagging data

  # Tagging Stuff
  conv_tag_fish_avail = array(data = 0, dim = c(conv_tag_max_liberty + 1, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes)) # Tags availiable for recapture
  conv_tag_fish_reporting = array(data = 0, dim = dim_fish_yr) # Tag reporting rate
  pred_conv_tag_fish_recap = array(data = 0, dim = c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets)) # predicted recaptures

  # Fishery Processes
  Fmort = array(0, dim = dim_fish) # Fishing mortality scalar
  dmr = array(0, dim = dim_fish) # Discard mortality rate
  tot_FAA = array(data = 0, dim = dim_pop_fish_age) # Total Fishing mortality at age
  ret_FAA = array(data = 0, dim = dim_pop_fish_age) # Retained Fishing mortality at age
  disc_FAA = array(data = 0, dim = dim_pop_fish_age) # Discarded Fishing mortality at age
  CAA = array(data = 0, dim = dim_pop_fish_age) # Retained Catch at age
  DAA = array(data = 0, dim = dim_pop_fish_age) # Discarded Catch at age
  CAL = array(data = 0, dim = dim_pop_fish_len) # Retained Catch at length
  DAL = array(data = 0, dim = dim_pop_fish_len) # Discarded Catch at length
  Fish_caal = if(do_caal == 1) array(data = 0, dim = dim_pop_fish_caal) else NULL # Retained catch at length and age
  Fish_caal_discard = if(do_caal == 1) array(data = 0, dim = dim_pop_fish_caal) else NULL # Discarded catch at length and age
  PredCatch = array(0, dim = dim_pop_fish) # Predicted retained catch (can be abundance or biomass)
  PredDiscard = array(0, dim = dim_pop_fish) # Predicted discarded catch (can be abundance, biomass, or abdunance or biomass fraction of retained catch)
  PredFishIdx = array(0, dim = dim_pop_fish) # Predicted fishery index
  fish_q = array(0, dim = dim_fish_yr) # Fishery catchability

  # Survey Processes
  SrvIAA = array(data = 0, dim = dim_pop_srv_age) # Survey index at age
  SrvIAL = array(data = 0, dim = dim_pop_srv_len) # Survey index at length
  Srv_caal = if(do_caal == 1) array(data = 0, dim = dim_pop_srv_caal) else NULL # Survey index at length and age
  PredSrvIdx = array(0, dim = dim_pop_srv) # Predicted survey index
  srv_q = array(0, dim = dim_srv_yr) # Survey catchability

  # Likelihoods (Not population-specific)
  Catch_nLL = array(0, dim = dim_fish) # Retained Fishery Catch Likelihoods
  Discard_nLL = array(0, dim = dim_fish) # Discarded Fishery Likelihoods
  FishIdx_nLL = array(0, dim = dim_fish) # Fishery Index Likelihoods
  FishAgeComps_nLL = array(data = 0, dim = dim_fish_comp) # Retained Fishery Age Comps Likelihoods
  FishLenComps_nLL = array(data = 0, dim = dim_fish_comp) # Retained Fishery Length Comps Likelihoods
  FishAgeComps_discard_nLL = array(data = 0, dim = dim_fish_comp) # Discarded Fishery Age Comps Likelihoods
  FishLenComps_discard_nLL = array(data = 0, dim = dim_fish_comp) # Discarded Fishery Length Comps Likelihoods
  CatchAA_nLL = array(0, dim = dim_fish_age) # Catch at age likelihoods
  CatchAA_pop_nLL = array(0, dim = dim_pop_fish_age) # Population-specific catch at age
  DiscardAA_nLL = array(0, dim = dim_fish_age) # Discard at age likelihoods
  DiscardAA_pop_nLL = array(0, dim = dim_pop_fish_age) # Population-specific discard at age
  SrvIdxAA_nLL = array(0, dim = dim_srv_age) # Survey index at age likelihoods
  SrvIdxAA_pop_nLL = array(0, dim = dim_pop_srv_age) # Population-specific survey index at age
  SrvIdx_nLL = array(0, dim = dim_srv) # Survey Index Likelihoods
  SrvAgeComps_nLL = array(data = 0, dim = dim_srv_comp) # Survey Age Comps Likelihoods
  SrvLenComps_nLL = array(data = 0, dim = dim_srv_comp) # Survey Length Comps Likelihoods
  conv_fish_tag_nLL = array(data = 0, dim = c(conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_regions, n_fish_fleets)) # Tagging Likelihoods

  # Likelihoods (population-specific)
  Catch_pop_nLL = array(0, dim = dim_pop_fish) # Pop-specific Catch Likelihoods
  Discard_pop_nLL = array(0, dim = dim_pop_fish) # Pop-specific Discarded Fishery Likelihoods
  FishIdx_pop_nLL = array(0, dim = dim_pop_fish) # Pop-specific Fishery Index Likelihoods
  FishAgeComps_pop_nLL = array(0, dim = dim_pop_fish_comp) # Pop-specific Retained Fishery Age Comps Likelihoods
  FishLenComps_pop_nLL = array(0, dim = dim_pop_fish_comp) # Pop-specific Retained Fishery Length Comps Likelihoods
  FishAgeComps_discard_pop_nLL = array(0, dim = dim_pop_fish_comp) # Pop-specific Discarded Fishery Age Comps Likelihoods
  FishLenComps_discard_pop_nLL = array(0, dim = dim_pop_fish_comp) # Pop-specific Discarded Fishery Length Comps Likelihoods
  SrvIdx_pop_nLL = array(0, dim = dim_pop_srv) # Pop-specific Survey Index Likelihoods
  SrvAgeComps_pop_nLL = array(0, dim = dim_pop_srv_comp) # Pop-specific Survey Age Comps Likelihoods
  SrvLenComps_pop_nLL = array(0, dim = dim_pop_srv_comp) # Pop-specific Survey Length Comps Likelihoods
  Fish_caal_nLL = array(0, dim = dim_fish_caal) # Fishery Conditional Age-at-Length Likelihoods
  Srv_caal_nLL = array(0, dim = dim_srv_caal) # Survey Conditional Age-at-Length Likelihoods

  # Conditional age-at-length is on only where use flags say so
  do_fish_caal = !is.null(UseFish_caal) && any(UseFish_caal == 1)
  do_srv_caal = !is.null(UseSrv_caal) && any(UseSrv_caal == 1)
  if((do_fish_caal || do_srv_caal) && do_caal != 1) stop("Conditional age-at-length data are being fit, but do_caal is 0, so the joint arrays at length and age were never built. Set do_caal = 1 in Setup_Mod_Biologicals.")

  # Penalties and Priors
  Fmort_nLL = array(0, dim = dim(ln_F_devs)) # Fishing Mortality Deviation penalty
  dmr_nLL = array(0, dim = dim(logit_dmr_devs)) # Discard Mortality Deviation penalty
  Rec_nLL = array(0, dim = dim(ln_RecDevs)) # Recruitment penalty
  Init_Rec_nLL = array(0, dim = dim(ln_InitDevs)) # Initial Recruitment penalty
  Init_Sex_nLL = array(0, dim = dim(ln_InitDevs)) # Initial age deviations, tie between sexes
  sel_nLL = 0 # Penalty for selectivity deviations
  fish_q_nLL = 0 # Prior/penalty for fishery q
  srv_q_nLL = 0 # Prior/penalty for survey q
  M_nLL = 0 # Penalty/Prior for natural mortality
  h_nLL = 0 # Prior for steepness
  Movement_nLL = 0 # Penalty for movement rates
  TagRep_nLL = 0 # penalty for tag reporting rate
  rec_prop_nLL = 0 # penalty / prior for recruitment proportions
  R0_nLL = 0 # prior for global R0
  growth_tv_nLL = 0 # penalty for time-varying growth deviations
  growth_semipar_nLL = 0 # process error for the semi-parametric deviations on mean length at age
  rinit_nLL = 0 # penalty on the initial recruitment offset from R0
  jnLL = 0 # Joint negative log likelihood

  # Parameter Transformations -----------------------------------------------
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

    # If move_type == 0 (unstructured markov)
    move_pars = move_pars, # movement parameters for unstructred markov
    move_devs = move_devs, # logit movement deviations
    use_fixed_movement = use_fixed_movement, # indicator for fixed movement
    Fixed_Movement = Fixed_Movement, # fixed movement matrix

    # If move_type == 1 (CTMC movement)
    ctmc_move_dat = ctmc_move_dat,
    preference_formula = preference_formula,
    diffusion_formula = diffusion_formula,
    log_move_diffusion_pars = log_move_diffusion_pars,
    move_preference_pars = move_preference_pars,
    area_r = area_r,
    adjacency_mat = adjacency_mat,
    ctmc_diffusion_bounds = ctmc_diffusion_bounds,
    ctmc_diffusion_eps = ctmc_diffusion_eps,
    seasdur = seasdur,
    ctmc_scale_by_seasdur = ctmc_scale_by_seasdur,
    expm_nsub = move_expm_nsub
  )

  # output movement stuff into model
  Mrate = out_move$Mrate
  Movement = out_move$Movement
  Movement_nLL = Movement_nLL + out_move$move_pen

  ## Natural Mortality Parameters (Set up) -----------------------------------
  if(use_fixed_natmort == 0) natmort = array(exp(ln_M[as.vector(M_blocks)]), dim = dim(M_blocks))
  if(use_fixed_natmort == 1) natmort = Fixed_natmort # Using fixed natural mortality

  ## Growth ------------------------------------------------------------------
  use_cohort_growth = growth_model != 0 && growth_tv_type == 1 # whether using growth cohort mode or not
  precomputed_mortality_yrs = if(use_cohort_growth) seq_len(max(1, growth_cohort_styr - 1)) else 1:n_yrs

  if(growth_model != 0) {

    # growth settings are the same for every call, whether the whole series is built
    # up front or one year at a time from inside the population loop
    growth_args = list(
      ln_growth_pars = ln_growth_pars,
      growth_A1 = growth_A1,
      growth_A2 = growth_A2,
      growth_L0 = growth_L0,
      growth_len_lower = growth_len_lower,
      growth_cv_type = growth_cv_type,
      growth_sd_type = growth_sd_type,
      growth_dist = growth_dist,
      growth_plus_group = growth_plus_group,
      growth_L2_asymptote = growth_L2_asymptote,
      derive_waa = derive_waa,
      wt_len_pars = wt_len_pars,
      ages = ages,
      seasdur = seasdur,
      spawn_seas = spawn_seas,
      t_spawn = t_spawn,
      n_pop = n_pop,
      n_regions = n_regions,
      n_seas = n_seas,
      n_sexes = n_sexes,
      t_fish = t_fish,
      t_srv = t_srv,
      ln_growth_devs = ln_growth_devs,
      growth_tv_model = growth_tv_model,
      growth_tv_link = growth_tv_link,
      growth_par_bounds = growth_par_bounds,
      ln_growth_semipar_devs = ln_growth_semipar_devs,
      growth_semipar = growth_semipar
    )

    # get growth function
    tmp_growth = do.call(Get_Growth, c(growth_args, list(
      n_yrs = n_yrs,
      n_fish_fleets = n_fish_fleets,
      n_srv_fleets = n_srv_fleets,
      growth_tv_type = growth_tv_type,
      growth_cohort_styr = growth_cohort_styr,
      years_eval = if(use_cohort_growth) seq_len(growth_cohort_styr) else NULL
    )))

    # get sizeage transition based on fleet timing
    SizeAgeTrans_fish = tmp_growth$SizeAgeTrans_fish
    SizeAgeTrans_srv = tmp_growth$SizeAgeTrans_srv
    SizeAgeTrans_spawn = tmp_growth$SizeAgeTrans_spawn

    # get growth values
    mean_LAA_fish = tmp_growth$mean_LAA_fish
    sd_LAA_fish = tmp_growth$sd_LAA_fish
    mean_LAA_srv = tmp_growth$mean_LAA_srv
    sd_LAA_srv = tmp_growth$sd_LAA_srv
    mean_LAA_spawn = tmp_growth$mean_LAA_spawn
    sd_LAA_spawn = tmp_growth$sd_LAA_spawn

    # get growth parameters
    Linf = tmp_growth$Linf
    L_beg = tmp_growth$L_beg
    growth_pars_y = tmp_growth$growth_pars_y

    # get waa from fleet timing
    if(derive_waa == 1) {
      WAA = tmp_growth$WAA
      WAA_fish = tmp_growth$WAA_fish
      WAA_srv = tmp_growth$WAA_srv
    }

    # bin midpoints the weight-length relationship is read at
    growth_len_mid_vals = growth_len_mid(growth_len_lower)

  } # growth module

  ## Selectivity -------------------------------------------------------------
  if(fish_selex_type == 0) fish_selex_bins = ages # if age-based selectivity
  if(fish_selex_type == 1) fish_selex_bins = lens # if length-based selectivity
  if(ret_selex_type == 0) ret_selex_bins = ages # if age-based selectivity
  if(ret_selex_type == 1) ret_selex_bins = lens # if length-based selectivity
  if(srv_selex_type == 0) srv_selex_bins = ages # if age-based selectivity
  if(srv_selex_type == 1) srv_selex_bins = lens # if length-based selectivity

  # Total fishery selectivity
  tmp_fish_sel = Get_Selex_Array(
    selex_type = fish_selex_type,
    bins = fish_selex_bins,
    sel_blocks = fish_sel_blocks,
    sel_model = fish_sel_model,
    fixed_sel_pars = fish_fixed_sel_pars,
    cont_tv_sel = cont_tv_fish_sel,
    ln_seldevs = ln_fishsel_devs,
    use_fixed_sel = use_fixed_fish_sel,
    bin_devs = ln_fishsel_bin_devs,
    bin_dev_bins = fish_sel_bin_dev_bins,
    sel_norm_bins = fish_sel_norm_bins,
    sex_par_offset = fishsel_sex_par_offset,
    sex_scale_offset = fishsel_sex_scale_offset,
    sex_apical_offset = fishsel_sex_apical_offset,
    sex_scale = ln_fishsel_sex_scale,
    nselbins = fish_sel_bicubic_nselbins,
    dbnrml_raw = fish_dbnrml_raw,
    dbnrml_startbin = fish_dbnrml_startbin,
    sel_input = fish_sel_input,
    bicubic_Wbin = fish_sel_bicubic_Wbin,
    bicubic_Wyr = fish_sel_bicubic_Wyr,
    bicubic_binnodes = fish_sel_bicubic_binnodes,
    bicubic_yrnodes = fish_sel_bicubic_yrnodes,
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_proj_yrs_devs = n_proj_yrs_devs,
    n_seas = n_seas,
    n_ages = n_ages,
    n_lens = n_lens,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets
  )
  fish_sel = tmp_fish_sel$sel
  fish_sel_l = tmp_fish_sel$sel_l

  # Fishery retention selectivity
  tmp_ret_sel = Get_Selex_Array(
    selex_type = ret_selex_type,
    bins = ret_selex_bins,
    sel_blocks = ret_sel_blocks,
    sel_model = ret_sel_model,
    fixed_sel_pars = ret_fixed_sel_pars,
    cont_tv_sel = cont_tv_ret_sel,
    ln_seldevs = ln_retsel_devs,
    use_fixed_sel = use_fixed_ret_sel,
    bin_devs = ln_retsel_bin_devs,
    bin_dev_bins = ret_sel_bin_dev_bins,
    sel_norm_bins = ret_sel_norm_bins,
    sex_par_offset = retsel_sex_par_offset,
    sex_scale_offset = retsel_sex_scale_offset,
    sex_apical_offset = retsel_sex_apical_offset,
    sex_scale = ln_retsel_sex_scale,
    nselbins = ret_sel_bicubic_nselbins,
    dbnrml_raw = ret_dbnrml_raw,
    dbnrml_startbin = ret_dbnrml_startbin,
    sel_input = ret_sel_input,
    bicubic_Wbin = ret_sel_bicubic_Wbin,
    bicubic_Wyr = ret_sel_bicubic_Wyr,
    bicubic_binnodes = ret_sel_bicubic_binnodes,
    bicubic_yrnodes = ret_sel_bicubic_yrnodes,
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_proj_yrs_devs = n_proj_yrs_devs,
    n_seas = n_seas,
    n_ages = n_ages,
    n_lens = n_lens,
    n_sexes = n_sexes,
    n_fleets = n_fish_fleets
  )
  ret_sel = tmp_ret_sel$sel
  ret_sel_l = tmp_ret_sel$sel_l

  # Survey selectivity
  tmp_srv_sel = Get_Selex_Array(
    selex_type = srv_selex_type,
    bins = srv_selex_bins,
    sel_blocks = srv_sel_blocks,
    sel_model = srv_sel_model,
    fixed_sel_pars = srv_fixed_sel_pars,
    cont_tv_sel = cont_tv_srv_sel,
    ln_seldevs = ln_srvsel_devs,
    use_fixed_sel = use_fixed_srv_sel,
    bin_devs = ln_srvsel_bin_devs,
    bin_dev_bins = srv_sel_bin_dev_bins,
    sel_norm_bins = srv_sel_norm_bins,
    sex_par_offset = srvsel_sex_par_offset,
    sex_scale_offset = srvsel_sex_scale_offset,
    sex_apical_offset = srvsel_sex_apical_offset,
    sex_scale = ln_srvsel_sex_scale,
    nselbins = srv_sel_bicubic_nselbins,
    dbnrml_raw = srv_dbnrml_raw,
    dbnrml_startbin = srv_dbnrml_startbin,
    sel_input = srv_sel_input,
    bicubic_Wbin = srv_sel_bicubic_Wbin,
    bicubic_Wyr = srv_sel_bicubic_Wyr,
    bicubic_binnodes = srv_sel_bicubic_binnodes,
    bicubic_yrnodes = srv_sel_bicubic_yrnodes,
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_proj_yrs_devs = n_proj_yrs_devs,
    n_seas = n_seas,
    n_ages = n_ages,
    n_lens = n_lens,
    n_sexes = n_sexes,
    n_fleets = n_srv_fleets
  )
  srv_sel = tmp_srv_sel$sel
  srv_sel_l = tmp_srv_sel$sel_l

  ## Recruitment Transformations and Bias Ramp (Methot and Taylor) -----------
  ### Recruitment Parameters -------------------------------------------------
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

  # reshape virgin recruitment into a population by time block matrix
  if(is.null(dim(ln_global_R0))) ln_R0_mat = matrix(ln_global_R0, length(ln_global_R0), 1) # one block, single column
  else ln_R0_mat = ln_global_R0 # already population by block

  R0 = exp(ln_R0_mat[, R0_ref_block]) # reference block R0, read by the initial age structure, apportionment, priors and the stock-recruit scale

  # R0 in each year, taken from that year's block (only annual recruitment is blocked)
  R0_yr = matrix(0, n_pop, n_yrs) # container
  for(p in 1:n_pop) {
    for(y in 1:n_yrs) {
      R0_blk_idx = 1 # single block when no blocking is set up
      if(!is.null(R0_blocks)) R0_blk_idx = R0_blocks[1, y, p] # otherwise this year's block
      R0_yr[p,y] = exp(ln_R0_mat[p, R0_blk_idx]) # virgin recruitment in this year
    } # end y loop
  } # end p loop

  for(p in 1:n_pop) R0_r[p,] = R0[p] * rec_region_prop[p,] # apportion R0 across regions

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

  ### Bias ramp --------------------------------------------------------------
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

  # The initial age deviations are techincally based on years before the first model year; option to construct bias ramp on initial deviations as well, based on
  # the recruitment bias ramp values
  init_bias_ramp = rep(1, n_ages - 1)
  if(do_rec_bias_ramp == 1) {
    d_init = 1 - (1:(n_ages - 1)) # deviation index of the year each initial age was born
    init_bias_ramp = rep(0, n_ages - 1)
    r1 = which(d_init >= bias_year[1] & d_init < bias_year[2])
    r2 = which(d_init >= bias_year[2] & d_init < bias_year[3])
    r3 = which(d_init >= bias_year[3] & d_init < bias_year[4])
    if(length(r1) > 0) init_bias_ramp[r1] = (d_init[r1] - bias_year[1]) / (bias_year[2] - bias_year[1])
    if(length(r2) > 0) init_bias_ramp[r2] = 1
    if(length(r3) > 0) init_bias_ramp[r3] = 1 - ((d_init[r3] - bias_year[3]) / (bias_year[4] - bias_year[3]))
    init_bias_ramp = init_bias_ramp * max_bias_ramp_fct
  }

  ## Mortality at Age --------------------------------------------------------
  missing_catch = is.na(ObsCatch) # TRUE = aggregate catch is missing, not a true recorded zero

  ### Mortality Settings -----------------------------------------------------
  # settings compute_mortality_year reads, the same for every year it is called for
  mortality_args = list(
    growth_model = growth_model,
    derive_waa = derive_waa,
    fish_selex_type = fish_selex_type,
    ret_selex_type = ret_selex_type,
    srv_selex_type = srv_selex_type,
    fish_waa_selected = fish_waa_selected,
    srv_waa_selected = srv_waa_selected,
    fish_sel_l = fish_sel_l,
    ret_sel_l = ret_sel_l,
    srv_sel_l = srv_sel_l,
    wt_len_pars = wt_len_pars,
    growth_len_mid_vals = growth_len_mid_vals,
    UseCatch = UseCatch,
    UseCatch_pop = UseCatch_pop,
    missing_catch = missing_catch,
    UseCatchAA = UseCatchAA,
    UseCatchAA_pop = UseCatchAA_pop,
    use_catch_aa = use_catch_aa,
    catch_seas_agg = Catch_seas_Type,
    catch_pop_seas_agg = Catch_pop_seas_Type,
    catch_aa_seas_agg = CatchAA_seas_Type,
    catch_aa_pop_seas_agg = CatchAA_pop_seas_Type,
    ln_F_mean = ln_F_mean,
    ln_F_devs = ln_F_devs,
    logit_dmr_mean = logit_dmr_mean,
    logit_dmr_devs = logit_dmr_devs,
    SizeAgeTrans = SizeAgeTrans,
    natmort = natmort,
    seasdur = seasdur,
    n_pop = n_pop,
    n_regions = n_regions,
    n_seas = n_seas,
    n_ages = n_ages,
    n_sexes = n_sexes,
    n_fish_fleets = n_fish_fleets
  )

  ### Mortality Containers ---------------------------------------------------
  # arrays compute_mortality_year fills a year at a time, and that
  # update_cohort_growth_and_mortality keeps filling inside the population loop
  mortality_state = list(
    Fmort = Fmort,
    dmr = dmr,
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    ret_FAA = ret_FAA,
    disc_FAA = disc_FAA,
    tot_FAA = tot_FAA,
    ZAA = ZAA,
    WAA = WAA,
    SizeAgeTrans_fish = SizeAgeTrans_fish,
    SizeAgeTrans_srv = SizeAgeTrans_srv
  )

  # a growth model derives its own quantities year by year, so they travel in the same state
  if(growth_model != 0) {
    mortality_state = c(mortality_state, list(
      tmp_growth = tmp_growth,
      SizeAgeTrans_spawn = SizeAgeTrans_spawn,
      mean_LAA_fish = mean_LAA_fish,
      sd_LAA_fish = sd_LAA_fish,
      mean_LAA_srv = mean_LAA_srv,
      sd_LAA_srv = sd_LAA_srv,
      mean_LAA_spawn = mean_LAA_spawn,
      sd_LAA_spawn = sd_LAA_spawn,
      Linf = Linf,
      L_beg = L_beg,
      growth_pars_y = growth_pars_y
    ))
    if(derive_waa == 1) mortality_state = c(mortality_state, list(WAA_fish = WAA_fish, WAA_srv = WAA_srv))
  } # growth module

  ### Mortality by Year ------------------------------------------------------
  for(y in precomputed_mortality_yrs)
    mortality_state = do.call(compute_mortality_year, c(mortality_args, list(y = y, state = mortality_state))) # get mortality values

  # get mortality values
  Fmort = mortality_state$Fmort
  dmr = mortality_state$dmr
  fish_sel = mortality_state$fish_sel
  ret_sel = mortality_state$ret_sel
  ret_FAA = mortality_state$ret_FAA
  disc_FAA = mortality_state$disc_FAA
  tot_FAA = mortality_state$tot_FAA
  ZAA = mortality_state$ZAA
  WAA = mortality_state$WAA

  # get sizeage transition based on fleet timing
  SizeAgeTrans_fish = mortality_state$SizeAgeTrans_fish
  SizeAgeTrans_srv = mortality_state$SizeAgeTrans_srv

  if(growth_model != 0) {

    tmp_growth = mortality_state$tmp_growth
    SizeAgeTrans_spawn = mortality_state$SizeAgeTrans_spawn

    # get growth values
    mean_LAA_fish = mortality_state$mean_LAA_fish
    sd_LAA_fish = mortality_state$sd_LAA_fish
    mean_LAA_srv = mortality_state$mean_LAA_srv
    sd_LAA_srv = mortality_state$sd_LAA_srv
    mean_LAA_spawn = mortality_state$mean_LAA_spawn
    sd_LAA_spawn = mortality_state$sd_LAA_spawn

    # get growth parameters
    Linf = mortality_state$Linf
    L_beg = mortality_state$L_beg
    growth_pars_y = mortality_state$growth_pars_y

    # get waa from fleet timing
    if(derive_waa == 1) {
      WAA_fish = mortality_state$WAA_fish
      WAA_srv = mortality_state$WAA_srv
    }

  } # growth module


  # Population Dynamics -----------------------------------------------------
  ## Initial Age Structure ---------------------------------------------------
  # Get initial F for age structure
  catch_flag_base = array(UseCatch[,1,,], dim = c(n_regions, n_seas, n_fish_fleets))
  catch_flag_pop = apply(UseCatch_pop[,,1,,,drop = FALSE], c(2,4,5), max)
  catch_flag = pmax(catch_flag_base, catch_flag_pop)

  # a fleet reporting once a year is still fish in every season of that year, so the flag allows fishing instead of where the obs is
  for(f in seq_len(n_fish_fleets)) {
    if(Catch_seas_Type[f] == 1 || Catch_pop_seas_Type[f] == 1) {
      for(r in 1:n_regions) if(any(catch_flag[r,,f] > 0)) catch_flag[r,,f] = 1
    }
  } # end f loop

  # Set up how initial F is determined; == 0 use proportion of ln_F_mean, otherwise, use absolute value
  init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets))
  for(r in 1:n_regions) for(seas in 1:n_seas) for(f in 1:n_fish_fleets) {
    init_F[r,seas,f] = if(init_F_form == 0) RTMB::plogis(init_F_par[r,seas,f]) * exp(ln_F_mean[r,seas,f]) else exp(init_F_par[r,seas,f])
  }
  init_F = init_F * catch_flag

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
    natmort = array(natmort[,,1,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)), # natural mortality in first year
    init_F = init_F, # initial F applied
    fish_sel = array(fish_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # total fishery selectivity in first year
    R0_r = if(use_rinit == 0) R0_r else rinit_r, # regional mean or virgin recruitment
    sexratio = array(sexratio[,,1,], dim = c(n_pop, n_regions, n_sexes)), # sex ratio in first year
    Movement = array(Movement[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
    do_recruits_move = do_recruits_move, # whether recruits move
    ln_InitDevs = ln_InitDevs, # initial deviations
    ret_sel = array(ret_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # retained fishery selectivity in first year
    dmr = array(dmr[,1,,], dim = c(n_regions, n_seas, n_fish_fleets)),
    Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # rates in first year
    move_timing = move_timing,
    expm_nsub = move_expm_nsub
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
    natmort = array(natmort[,,1,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes)), # natural mortality in first year
    init_F = array(0, dim = c(n_regions, n_seas, n_fish_fleets)), # initial F applied (0 for unfished)
    fish_sel = array(fish_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # total fishery selectivity in first year
    R0_r = if(use_rinit == 0) R0_r else rinit_r, # regional mean or virgin recruitment
    sexratio = array(sexratio[,,1,], dim = c(n_pop, n_regions, n_sexes)), # sex ratio in first year
    Movement = array(Movement[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # movement in first year
    do_recruits_move = do_recruits_move, # whether recruits move
    ln_InitDevs = ln_InitDevs, # initial deviations
    ret_sel = array(ret_sel[,,1,,,,], dim = c(n_pop, n_regions, n_seas, n_ages, n_sexes, n_fish_fleets)), # retained fishery selectivity in first year
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)), # unfished
    Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # rates in first year
    move_timing = move_timing,
    expm_nsub = move_expm_nsub
  )

  # Input into model arrays (first year and season) - and add lognormal mean adjustment
  NAA[,,1,1,,] = Init_Fished_NAA
  NAA0[,,1,1,,] = Init_Unfished_NAA

  ## Cohort Growth and Mortality ---------------------------------------------
  # build out cohort influenced mortality, since mortality dynamics can change if large cohorts enter
  update_cohort_growth_and_mortality = NULL
  if(use_cohort_growth) {

    update_cohort_growth_and_mortality = function(y, NAA_y, state) {
      if(y >= growth_cohort_styr) {
        state$tmp_growth = do.call(Get_Growth_Year, c(growth_args, list(growth = state$tmp_growth, y = y, NAA_y = NAA_y)))
        state = growth_take_year(state, state$tmp_growth, y, derive_waa) # this year's growth quantity
        # only WAA_fish and WAA_srv need this year's growth. ZAA is recomputed as a side effect and only changes under length-based selectivity
        state = do.call(compute_mortality_year, c(mortality_args, list(y = y, state = state))) # get new selex after cohort growth, then compute mortality
      }

      list(state = state,
           ZAA_y = state$ZAA[,,y,,,, drop = FALSE],
           WAA_y = state$WAA[,,y,,,, drop = FALSE],
           MatAA_y = MatAA[,,y,,,, drop = FALSE])
    }
  }

  ## Population Projection ---------------------------------------------------

  # scale the stock-recruit curve is built at (based on R0, a separate R0 just for the SR curve, or rinit)
  if(sr_R0_spec == 0) sr_R0 = R0 # shared, virgin recruitment from the reference block
  if(sr_R0_spec == 1) sr_R0 = exp(ln_sr_R0) # est, a scale estimated for the curve alone
  if(sr_R0_spec == 2) sr_R0 = rinit # rinit, the initial equilibrium recruitment

  tmp_pop_proj = get_population_projection(

    # model dimensions
    n_pop = n_pop, # populations
    n_regions = n_regions, # regions
    n_seas = n_seas, # seasons
    n_ages = n_ages, # ages
    n_sexes = n_sexes, # sexes
    n_yrs = n_yrs, # years
    n_fish_fleets = n_fish_fleets, # fishery fleets
    n_est_rec_devs = n_est_rec_devs, # number of estimated recruitment deviations

    # recruitment
    rec_lag = rec_lag, # years between spawning and recruitment
    rec_model = rec_model, # mean recruitment, Beverton-Holt or Ricker
    rec_dd = rec_dd, # whether the curve is local or global density dependent
    R0 = R0, # reference block virgin recruitment
    R0_yr = R0_yr, # virgin recruitment in each year
    rec_region_prop = rec_region_prop, # recruitment proportions by region
    rec_seas_prop = rec_seas_prop, # recruitment proportions by season
    h_trans = h_trans, # steepness, bounded 0.2 to 1
    ln_RecDevs = ln_RecDevs, # recruitment deviations, log scale
    sexratio = sexratio, # recruitment sex ratio
    natal_region = natal_region, # natal region of each population
    SR_ref_yr = SR_ref_yr, # year the stock-recruit curve takes its biology from
    sr_penalty = sr_penalty, # stock-recruit penalty under mean recruitment
    sr_R0 = sr_R0, # scale that penalty's curve is built at

    # timing
    t_spawn = t_spawn, # fraction of the season elapsed at spawning
    spawn_seas = spawn_seas, # season spawning happens in
    seasdur = seasdur, # duration of each season, fraction of a year

    # biologicals
    WAA = WAA, # weight at age
    MatAA = MatAA, # maturity at age
    natmort = natmort, # natural mortality, rate per year in each season

    # movement
    Movement = Movement, # movement probabilities
    Mrate = Mrate, # movement rates, continuous time formulation
    move_timing = move_timing, # where movement sits relative to mortality
    expm_nsub = move_expm_nsub, # substeps in the matrix exponential
    stray_rate = stray_rate, # straying away from the natal region
    sgl_seas_spawning_movement = sgl_seas_spawning_movement, # movement back to spawn in a single season model
    do_recruits_move = do_recruits_move, # whether recruits move

    # mortality
    init_F = init_F, # fishing mortality generating the initial age structure
    fish_sel = fish_sel, # total fishery selectivity
    ret_sel = ret_sel, # retained fishery selectivity
    dmr = dmr, # discard mortality rate
    ZAA = ZAA, # total mortality at age

    # containers the projection fills
    NAA = NAA, # numbers at age
    NAA0 = NAA0, # unfished numbers at age
    NAA_bef = NAA_bef, # numbers at age before movement
    NAA_aft = NAA_aft, # numbers at age after movement
    Rec = Rec, # recruitment
    SSB = SSB, # spawning stock biomass
    Total_Biom = Total_Biom, # total biomass
    Dynamic_SSB0 = Dynamic_SSB0, # unfished spawning biomass with the estimated deviations
    eff_SSB = eff_SSB, # spawning biomass the curve reads under natal homing

    # cohort growth, updated year by year inside the projection
    growth_mortality_year_fn = update_cohort_growth_and_mortality, # builds one year of growth and mortality
    growth_mortality_state = mortality_state, # the arrays that function reads and writes

    # state-space numbers at age
    n_est_naa_re = n_est_naa_re, # number of estimated state cells
    ln_NAA = ln_NAA, # the state itself, log scale
    naa_re_ages = naa_re_ages, # ages the state runs over
    naa_re_yrs = naa_re_yrs, # years the state runs over
    naa_re_seas = naa_re_seas # seasons the state runs over
  )

  ### Cohort Growth and Mortality Outputs ------------------------------------
  # under cohort growth the projection builds growth and mortality one year at a time, so read them back out
  if(use_cohort_growth) {

    mortality_state = tmp_pop_proj$growth_mortality_state

    # get mortality values
    Fmort = mortality_state$Fmort
    dmr = mortality_state$dmr
    fish_sel = mortality_state$fish_sel
    ret_sel = mortality_state$ret_sel
    ret_FAA = mortality_state$ret_FAA
    disc_FAA = mortality_state$disc_FAA
    tot_FAA = mortality_state$tot_FAA
    ZAA = mortality_state$ZAA

    if(growth_model != 0) {

      tmp_growth = mortality_state$tmp_growth

      # get sizeage transition based on fleet timing
      SizeAgeTrans_fish = mortality_state$SizeAgeTrans_fish
      SizeAgeTrans_srv = mortality_state$SizeAgeTrans_srv
      SizeAgeTrans_spawn = mortality_state$SizeAgeTrans_spawn

      # get growth values
      mean_LAA_fish = mortality_state$mean_LAA_fish
      sd_LAA_fish = mortality_state$sd_LAA_fish
      mean_LAA_srv = mortality_state$mean_LAA_srv
      sd_LAA_srv = mortality_state$sd_LAA_srv
      mean_LAA_spawn = mortality_state$mean_LAA_spawn
      sd_LAA_spawn = mortality_state$sd_LAA_spawn

      # get growth parameters
      Linf = mortality_state$Linf
      L_beg = mortality_state$L_beg
      growth_pars_y = mortality_state$growth_pars_y

      # get waa from fleet timing
      if(derive_waa == 1) {
        WAA = mortality_state$WAA
        WAA_fish = mortality_state$WAA_fish
        WAA_srv = mortality_state$WAA_srv
      }

    } # growth module

  } # cohort growth

  # get numbers at age
  NAA = tmp_pop_proj$NAA
  NAA0 = tmp_pop_proj$NAA0 # unfished
  NAA_bef = tmp_pop_proj$NAA_bef # before movement
  NAA_aft = tmp_pop_proj$NAA_aft # after movement
  NAA_int = tmp_pop_proj$NAA_int # season-integrated abundance (move_timing == 2 only)
  NAA_pred = tmp_pop_proj$NAA_pred # deterministic prediction behind the state (NAA_re != none only)
  NAA_scalar = tmp_pop_proj$NAA_scalar # factor the state applied, one wherever it did not

  # get derived quantities
  Rec = tmp_pop_proj$Rec
  SSB = tmp_pop_proj$SSB
  Total_Biom = tmp_pop_proj$Total_Biom
  Dynamic_SSB0 = tmp_pop_proj$Dynamic_SSB0 # unfished spawning biomass with the estimated deviations
  eff_SSB = tmp_pop_proj$eff_SSB # effective spawning biomass the curve reads under natal homing
  Aggregated_SSB = tmp_pop_proj$Aggregated_SSB # summed across regions
  Dynamic_Aggregated_SSB0 = tmp_pop_proj$Dynamic_Aggregated_SSB0 # summed across regions
  SR_pred = tmp_pop_proj$SR_pred # stock-recruit prediction, only used under mean recruitment

  ## State-Space Numbers at Age Penalty --------------------------------------
  # penalized here rather than in the dynamics, so the deterministic prediction and the state are both available
  NAA_state_nLL = 0
  if(n_est_naa_re > 0) {

    sigmaNAA = array(exp(ln_sigmaNAA)[as.vector(naa_sigma_blocks)], dim = dim(naa_sigma_blocks)) # deviation sd, read from each cell's block

    NAA_state_nLL = Get_NAA_state_penalty(
      ln_NAA = ln_NAA, # the state itself, log scale
      NAA_pred = NAA_pred, # deterministic mortality and ageing prediction
      sigmaNAA = sigmaNAA, # deviation sd
      naa_re_ages = naa_re_ages, # ages the state runs over
      naa_re_yrs = naa_re_yrs, # years the state runs over
      naa_re_seas = naa_re_seas, # seasons the state runs over
      NAA_re = NAA_re, # process error form over ages and years
      NAA_pe_pars = NAA_pe_pars, # parameters of that form
      NAA_re_region = NAA_re_region, # correlation across regions
      NAA_region_corr_pars = NAA_region_corr_pars,
      NAA_re_pop = NAA_re_pop, # correlation across populations
      NAA_pop_corr_pars = NAA_pop_corr_pars,
      NAA_re_sex = NAA_re_sex, # correlation across sexes
      NAA_sex_corr_pars = NAA_sex_corr_pars,
      NAA_re_season = NAA_re_season, # correlation across seasons
      NAA_season_corr_pars = NAA_season_corr_pars,
      naa_re_where = naa_re_where # population by region cells the state runs over
    )

  } # end if estimating a state on numbers at age

  # Observation Models ------------------------------------------------------
  ## Fishery Observation Model -----------------------------------------------
  tmp_fish_obs = get_fishery_observation_model(
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fish_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    fish_q_blocks = fish_q_blocks,
    ln_fish_q = ln_fish_q,
    fish_q = fish_q,
    ret_FAA = ret_FAA,
    disc_FAA = disc_FAA,
    ZAA = ZAA,
    NAA = NAA,
    CAA = CAA,
    DAA = DAA,
    CAL = CAL,
    DAL = DAL,
    PredCatch = PredCatch,
    PredDiscard = PredDiscard,
    PredFishIdx = PredFishIdx,
    fit_lengths = fit_lengths,
    SizeAgeTrans = SizeAgeTrans,
    SizeAgeTrans_fish = SizeAgeTrans_fish,
    fish_len_comp_sel = fish_len_comp_sel,
    fish_selex_type = fish_selex_type,
    ret_selex_type = ret_selex_type,
    fish_sel_l = fish_sel_l,
    ret_sel_l = ret_sel_l,
    Fmort = Fmort,
    catch_units = catch_units,
    discard_units = discard_units,
    WAA_fish = WAA_fish,
    dmr = dmr,
    fish_idx_type = fish_idx_type,
    fish_sel = fish_sel,
    ret_sel = ret_sel,
    Mrate = Mrate,
    move_timing = move_timing,
    expm_nsub = move_expm_nsub,
    seasdur = seasdur,
    NAA_int = if(move_timing == 2) NAA_int else NULL, # reuse the dynamics' exponentials
    t_fish = t_fish,
    fish_idx_ages = fish_idx_ages,
    fish_q_type = fish_q_type,
    do_fish_q_cov = do_fish_q_cov,
    fish_q_cov = fish_q_cov,
    fish_q_coeff = fish_q_coeff,
    ObsFishIdx = ObsFishIdx,
    UseFishIdx = UseFishIdx,
    do_caal = do_caal,
    Fish_caal = Fish_caal,
    Fish_caal_discard = Fish_caal_discard
  )

  fish_q = tmp_fish_obs$fish_q; CAA = tmp_fish_obs$CAA; DAA = tmp_fish_obs$DAA
  CAL = tmp_fish_obs$CAL; DAL = tmp_fish_obs$DAL
  Fish_caal = tmp_fish_obs$Fish_caal; Fish_caal_discard = tmp_fish_obs$Fish_caal_discard
  PredCatch = tmp_fish_obs$PredCatch; PredDiscard = tmp_fish_obs$PredDiscard; PredFishIdx = tmp_fish_obs$PredFishIdx

  ## Survey Observation Model ------------------------------------------------

  # Get recruitment index - computed as an anomaly w/ a bias adjustment as an addition, rather than sutraction
  RecDev_anom = array(0, dim = dim(ln_RecDevs))
  if(any(srv_idx_type == 2)) {
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {
        sigma_idx = ifelse(n_pop == 1 && rec_dd == 0, r, natal_region[p])
        for(d in 1:n_est_rec_devs) {
          sigmaR_d = if(d < sigmaR_switch) exp(ln_sigmaR[1,p,sigma_idx]) else exp(ln_sigmaR[2,p,sigma_idx])
          RecDev_anom[p,r,d] = ln_RecDevs[p,r,d] + 0.5 * sigmaR_d^2 * bias_ramp[d]
        } # end d loop
      } # end r loop
    } # end p loop
  } # end if any fleet observes the recruitment deviations

  tmp_srv_obs = get_survey_observation_model(
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_srv_fleets = n_srv_fleets,
    n_sexes = n_sexes,
    srv_q_blocks = srv_q_blocks,
    ln_srv_q = ln_srv_q,
    srv_q = srv_q,
    do_srv_q_cov = do_srv_q_cov,
    srv_q_cov = srv_q_cov,
    srv_q_coeff = srv_q_coeff,
    srv_selex_type = srv_selex_type,
    srv_sel = srv_sel,
    srv_sel_l = srv_sel_l,
    SizeAgeTrans = SizeAgeTrans,
    SizeAgeTrans_srv = SizeAgeTrans_srv,
    srv_len_comp_sel = srv_len_comp_sel,
    NAA = NAA,
    ZAA = ZAA,
    t_srv = t_srv,
    SrvIAA = SrvIAA,
    fit_lengths = fit_lengths,
    SrvIAL = SrvIAL,
    srv_idx_type = srv_idx_type,
    WAA_srv = WAA_srv,
    PredSrvIdx = PredSrvIdx,
    Mrate = Mrate,
    move_timing = move_timing,
    expm_nsub = move_expm_nsub,
    seasdur = seasdur,
    srv_idx_ages = srv_idx_ages,
    srv_q_type = srv_q_type,
    ObsSrvIdx = ObsSrvIdx,
    UseSrvIdx = UseSrvIdx,
    RecDev_anom = RecDev_anom,
    do_caal = do_caal,
    Srv_caal = Srv_caal
  )

  srv_q = tmp_srv_obs$srv_q; srv_sel = tmp_srv_obs$srv_sel
  SrvIAA = tmp_srv_obs$SrvIAA; SrvIAL = tmp_srv_obs$SrvIAL; PredSrvIdx = tmp_srv_obs$PredSrvIdx
  Srv_caal = tmp_srv_obs$Srv_caal


  ## Conventional Tagging Observation Model ----------------------------------
  if(any(use_conv_fish_tagging == 1)) {

    tmp_tag_obs = get_tagging_observation_model(
      n_fish_fleets = n_fish_fleets,
      n_regions = n_regions,
      n_conv_tag_cohorts = n_conv_tag_cohorts,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_pop = n_pop,
      n_ages = n_ages,
      n_sexes = n_sexes,
      conv_tag_fish_reporting_blocks = conv_tag_fish_reporting_blocks,
      conv_tag_fish_reporting_pars = conv_tag_fish_reporting_pars,
      conv_tag_fish_reporting = conv_tag_fish_reporting,
      conv_tag_release_indicator = conv_tag_release_indicator,
      conv_tag_max_liberty = conv_tag_max_liberty,
      use_conv_fish_tagging = use_conv_fish_tagging,
      Fmort = Fmort,
      fish_sel = fish_sel,
      ret_sel = ret_sel,
      dmr = dmr,
      natmort = natmort,
      seasdur = seasdur,
      ln_conv_tag_shed = ln_conv_tag_shed,
      conv_tag_t_tagging = conv_tag_t_tagging,
      conv_tagged_fish = conv_tagged_fish,
      conv_fish_tag_attr = conv_fish_tag_attr,
      conv_tag_release_platform = conv_tag_release_platform,
      srv_sel = srv_sel,
      NAA_bef = NAA_bef,
      ln_init_conv_tag_mort = ln_init_conv_tag_mort,
      do_recruits_move = do_recruits_move,
      Movement = Movement,
      conv_tag_fish_avail = conv_tag_fish_avail,
      pred_conv_tag_fish_recap = pred_conv_tag_fish_recap,
      Mrate = Mrate,
      move_timing = move_timing,
      expm_nsub = move_expm_nsub,
      NAA_scalar = if(n_est_naa_re > 0) NAA_scalar else NULL
    )

    conv_tag_fish_reporting = tmp_tag_obs$conv_tag_fish_reporting
    conv_tag_fish_avail = tmp_tag_obs$conv_tag_fish_avail
    pred_conv_tag_fish_recap = tmp_tag_obs$pred_conv_tag_fish_recap

  } # end if for using tagging data


  # Likelihood Equations ----------------------------------------------------

  ## Shared Observation Error and Arrays -------------------------------------
  ### Index observation error ------------------------------------------------
  # standard deviations for every index, fishery and survey
  FishIdx_SD = build_idx_sd(ObsFishIdx_SE, ln_sigmaFishIdx, sigmaFishIdx_form)
  FishIdx_pop_SD = build_idx_sd(ObsFishIdx_pop_SE, ln_sigmaFishIdx_pop, sigmaFishIdx_pop_form)
  SrvIdx_SD = build_idx_sd(ObsSrvIdx_SE, ln_sigmaSrvIdx, sigmaSrvIdx_form)
  SrvIdx_pop_SD = build_idx_sd(ObsSrvIdx_pop_SE, ln_sigmaSrvIdx_pop, sigmaSrvIdx_pop_form)

  ### At-Age Observation Arrays ----------------------------------------------
  # model predictions the age-disaggregated data sources are summed from, shared by the fishery and survey blocks below
  at_age_arrays = list(CAA = CAA, # catch at age
                       DAA = DAA, # discards at age
                       SrvIAA = SrvIAA, # survey available numbers at age
                       WAA_fish = WAA_fish, # weight at age at the fishery
                       dmr = dmr, # discard mortality rate
                       catch_units = catch_units, # numbers or biomass
                       discard_units = discard_units) # numbers or biomass

  ## Fishery Likelihoods -----------------------------------------------------
  ### Retained Fishery Catches (Regional) ------------------------------------
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
                                               log(get_seas_pred(PredCatch, r, y, seas, f, Catch_seas_Type[f])),
                                               exp(ln_sigmaC[r,y,seas,f]), TRUE)
    }
  }

  ### Retained Fishery Catch at Age (Regional) -------------------------------
  # put observations on the scale this fleet's likelihood uses, then register them
  # so OSA residuals read the same vector the objective does
  ObsCatchAA = prep_at_age_obs(ObsCatchAA, UseCatchAA, CatchAA_LikeType)
  if(length(ObsCatchAA)) ObsCatchAA = RTMB::OBS(ObsCatchAA) # setup OSA

  caa = get_at_age_source_nLL(
    obs_t = ObsCatchAA, # observed catch at age
    use = UseCatchAA, # cells that are fit
    ln_sigma = ln_sigmaCAA, # observation error, log scale
    source = "catch", # prediction is summed from CAA
    pop = FALSE, # regional data source
    arrays = at_age_arrays, # model predictions
    obs_se = ObsCatchAA_SE, # reported standard errors
    sd_form = CatchAA_sigma_form, # how those standard errors enter
    like_type = CatchAA_LikeType, # lognormal or normal
    const = 0, # constant added inside the log
    corr_type = AgeObsCorr_catch, # correlation form across ages
    trans_rho = trans_rho_catch, # correlation across ages, unconstrained
    trans_rho_year = trans_rho_catch_year, # correlation across years, unconstrained
    us_pars = trans_rho_catch_us, # unstructured correlation parameters
    aa_type = CatchAA_Type, # dims the observations are split by
    seas_agg = CatchAA_seas_Type # whether this source is fit as a season total
  )

  CatchAA_nLL = caa$nLL # nLL
  PredCatchAA = caa$pred # predicted catch at age

  ### Retained Fishery Catches (Population-Specific) -------------------------
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
                                                     log(get_seas_pred_pop(PredCatch, p, r, y, seas, f, Catch_pop_seas_Type[f])),
                                                     exp(ln_sigmaC_pop[p,r,y,seas,f]), TRUE)
    }
  }


  ### Retained Fishery Catch at Age (Population-Specific) --------------------
  ObsCatchAA_pop = prep_at_age_obs(ObsCatchAA_pop, UseCatchAA_pop, CatchAA_pop_LikeType)
  if(length(ObsCatchAA_pop)) ObsCatchAA_pop = RTMB::OBS(ObsCatchAA_pop) # setup OSA

  caa_pop = get_at_age_source_nLL(
    obs_t = ObsCatchAA_pop, # observed catch at age
    use = UseCatchAA_pop, # cells that are fit
    ln_sigma = ln_sigmaCAA_pop, # observation error, log scale
    source = "catch", # prediction is summed from CAA
    pop = TRUE, # population-specific data source, never summed across populations
    arrays = at_age_arrays, # model predictions
    obs_se = ObsCatchAA_pop_SE, # reported standard errors
    sd_form = CatchAA_pop_sigma_form, # how those standard errors enter
    like_type = CatchAA_pop_LikeType, # lognormal or normal
    const = 0, # constant added inside the log
    corr_type = AgeObsCorr_catch_pop, # correlation form across ages
    trans_rho = trans_rho_catch_pop, # correlation across ages, unconstrained
    trans_rho_year = trans_rho_catch_pop_year, # correlation across years, unconstrained
    us_pars = trans_rho_catch_pop_us, # unstructured correlation parameters
    aa_type = CatchAA_pop_Type, # dims the observations are split by
    seas_agg = CatchAA_pop_seas_Type # whether this source is fit as a season total
  )

  CatchAA_pop_nLL = caa_pop$nLL # nLL
  PredCatchAA_pop = caa_pop$pred # predicted catch at age


  ### Discarded Fishery Discards (Regional) ----------------------------------
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
                                                 log(get_seas_pred(PredDiscard, r, y, seas, f, Discard_seas_Type[f])),
                                                 exp(ln_sigmaD[r,y,seas,f]), TRUE)
    }
  }


  ### Discarded Fishery Discard at Age (Regional) ----------------------------
  ObsDiscardAA = prep_at_age_obs(ObsDiscardAA, UseDiscardAA, DiscardAA_LikeType)
  if(length(ObsDiscardAA)) ObsDiscardAA = RTMB::OBS(ObsDiscardAA) # setup OSA

  daa = get_at_age_source_nLL(
    obs_t = ObsDiscardAA, # observed discards at age
    use = UseDiscardAA, # cells that are fit
    ln_sigma = ln_sigmaDAA, # observation error, log scale
    source = "discard", # prediction is summed from DAA
    pop = FALSE, # regional data source
    arrays = at_age_arrays, # model predictions
    obs_se = ObsDiscardAA_SE, # reported standard errors
    sd_form = DiscardAA_sigma_form, # how those standard errors enter
    like_type = DiscardAA_LikeType, # lognormal or normal
    const = 0, # constant added inside the log
    corr_type = AgeObsCorr_discard, # correlation form across ages
    trans_rho = trans_rho_discard, # correlation across ages, unconstrained
    trans_rho_year = trans_rho_discard_year, # correlation across years, unconstrained
    us_pars = trans_rho_discard_us, # unstructured correlation parameters
    aa_type = DiscardAA_Type, # dims the observations are split by
    seas_agg = DiscardAA_seas_Type # whether this source is fit as a season total
  )

  DiscardAA_nLL = daa$nLL # nLL
  PredDiscardAA = daa$pred # predicted discards at age


  ### Discarded Fishery Discards (Population-Specific) -----------------------
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
                                                       log(get_seas_pred_pop(PredDiscard, p, r, y, seas, f, Discard_pop_seas_Type[f])),
                                                       exp(ln_sigmaD_pop[p,r,y,seas,f]), TRUE)
    }
  }


  ### Discarded Fishery Discard at Age (Population-Specific) -----------------
  ObsDiscardAA_pop = prep_at_age_obs(ObsDiscardAA_pop, UseDiscardAA_pop, DiscardAA_pop_LikeType)
  if(length(ObsDiscardAA_pop)) ObsDiscardAA_pop = RTMB::OBS(ObsDiscardAA_pop) # setup OSA

  daa_pop = get_at_age_source_nLL(
    obs_t = ObsDiscardAA_pop, # observed discards at age
    use = UseDiscardAA_pop, # cells that are fit
    ln_sigma = ln_sigmaDAA_pop, # observation error, log scale
    source = "discard", # prediction is summed from DAA
    pop = TRUE, # population-specific data source, never summed across populations
    arrays = at_age_arrays, # model predictions
    obs_se = ObsDiscardAA_pop_SE, # reported standard errors
    sd_form = DiscardAA_pop_sigma_form, # how those standard errors enter
    like_type = DiscardAA_pop_LikeType, # lognormal or normal
    const = 0, # constant added inside the log
    corr_type = AgeObsCorr_discard_pop, # correlation form across ages
    trans_rho = trans_rho_discard_pop, # correlation across ages, unconstrained
    trans_rho_year = trans_rho_discard_pop_year, # correlation across years, unconstrained
    us_pars = trans_rho_discard_pop_us, # unstructured correlation parameters
    aa_type = DiscardAA_pop_Type, # dims the observations are split by
    seas_agg = DiscardAA_pop_seas_Type # whether this source is fit as a season total
  )

  DiscardAA_pop_nLL = daa_pop$nLL # nLL
  PredDiscardAA_pop = daa_pop$pred # predicted discards at age


  ### Retained Fishery Indices (Regional) ------------------------------------
  # reads ObsFishIdx as an array. the lognormal block below rebinds that name to a flat
  # log-scale vector, since OSA files residuals under the name they were registered on
  FishIdx_nLL = get_index_regional_nLL(
    nLL_arr = FishIdx_nLL, # container
    Use = UseFishIdx, # cells that are fit
    Obs = ObsFishIdx, # observed index
    Pred = PredFishIdx, # predicted index
    SD = FishIdx_SD, # index standard deviation
    LikeType = FishIdx_LikeType, # likelihood for each fleet
    Cov = FishIdx_Cov, # covariance, read only by a multivariate normal
    seas_Type = FishIdx_seas_Type, # whether each fleet is fit as a season total
    const = addtofishidx, # constant added inside the log
    n_fleets = n_fish_fleets
  )

  # only lognormal fleets register OSA observations, so drop the fleets already done above
  UseFishIdx_ln = UseFishIdx
  for(f in 1:n_fish_fleets) if(FishIdx_LikeType[f] != 0) UseFishIdx_ln[,,,f] = 0

  if(any(UseFishIdx_ln == 1)) { # setup OSA

    valid_idx_ir = which(UseFishIdx_ln == 1) # position of each fitted cell
    ObsFishIdx_map = arrayInd(valid_idx_ir, dim(UseFishIdx_ln)) # region, year, season and fleet of each
    ObsFishIdx = log(ObsFishIdx[valid_idx_ir] + addtofishidx) # observed index, log scale
    ObsFishIdx = RTMB::OBS(ObsFishIdx) # register for OSA residuals

    FishIdx_nLL = eval_index_osa_nLL(
      nLL_arr = FishIdx_nLL, # container
      obs_vec = ObsFishIdx, # registered observations, log scale
      obs_map = ObsFishIdx_map, # index of each observation
      Pred = PredFishIdx, # predicted index
      SD = FishIdx_SD, # index standard deviation
      seas_Type = FishIdx_seas_Type, # whether each fleet is fit as a season total
      const = addtofishidx, # constant added inside the log
      pop = FALSE
    )
  }

  ### Retained Fishery Indices (Population-Specific) -------------------------
  # reads ObsFishIdx_pop as an array, which the lognormal block below rebinds to a flat vector.
  # a multivariate normal covariance describes the regional series, so mvn fleets go there too
  FishIdx_pop_nLL = get_index_pop_nLL(
    nLL_arr = FishIdx_pop_nLL, # container
    Use = UseFishIdx_pop, # cells that are fit
    Obs = ObsFishIdx_pop, # observed index
    Pred = PredFishIdx, # predicted index
    SD = FishIdx_pop_SD, # index standard deviation
    LikeType = FishIdx_LikeType, # likelihood for each fleet
    seas_Type = FishIdx_pop_seas_Type, # whether each fleet is fit as a season total
    n_fleets = n_fish_fleets
  )

  # only lognormal fleets register OSA observations, so drop the fleets already done above
  UseFishIdx_pop_ln = UseFishIdx_pop
  for(f in 1:n_fish_fleets) if(FishIdx_LikeType[f] == 1) UseFishIdx_pop_ln[,,,,f] = 0

  if(any(UseFishIdx_pop_ln == 1)) { # setup OSA

    valid_idx_ip = which(UseFishIdx_pop_ln == 1) # position of each fitted cell
    ObsFishIdx_pop_map = arrayInd(valid_idx_ip, dim(UseFishIdx_pop_ln)) # population, region, year, season and fleet of each
    ObsFishIdx_pop = log(ObsFishIdx_pop[valid_idx_ip] + addtofishidx) # observed index, log scale
    ObsFishIdx_pop = RTMB::OBS(ObsFishIdx_pop) # register for OSA residuals

    FishIdx_pop_nLL = eval_index_osa_nLL(
      nLL_arr = FishIdx_pop_nLL, # container
      obs_vec = ObsFishIdx_pop, # registered observations, log scale
      obs_map = ObsFishIdx_pop_map, # index of each observation
      Pred = PredFishIdx, # predicted index
      SD = FishIdx_pop_SD, # index standard deviation
      seas_Type = FishIdx_pop_seas_Type, # whether each fleet is fit as a season total
      const = addtofishidx, # constant added inside the log
      pop = TRUE
    )
  }

  ### Retained Fishery Compositions (Region-Specific) ------------------------
  # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
  # looks the residuals up under that name, so each data source is registered here by name
  ObsFishAgeComps_osa = pack_comp_source_osa(
    ObsArr = ObsFishAgeComps, # observed age compositions
    ISSArr = ISS_FishAgeComps, # input sample size
    WtArr = Wt_FishAgeComps, # multinomial weight
    UseArr = UseFishAgeComps, # cells that are fit
    TypeMat = FishAgeComps_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_LikeType, # likelihood for each fleet
    comp_bins_spec = FishAgeComps_bins, # bins this data source is fit over
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    n_pop = n_pop,
    pop = FALSE, # regional data source
    addtocomp = addtocomp,
    do_internal_comp_osa = do_internal_comp_osa
  )

  ObsFishAgeComps_osa_discrete = ObsFishAgeComps_osa$discrete # fleets on a discrete likelihood
  ObsFishAgeComps_osa_continuous = ObsFishAgeComps_osa$continuous # fleets on a continuous likelihood
  if(!is.null(ObsFishAgeComps_osa_discrete)) ObsFishAgeComps_osa_discrete = RTMB::OBS(ObsFishAgeComps_osa_discrete)
  if(!is.null(ObsFishAgeComps_osa_continuous)) ObsFishAgeComps_osa_continuous = RTMB::OBS(ObsFishAgeComps_osa_continuous)

  FishAgeComps_nLL = get_comp_source_nLL(
    nLL_arr = FishAgeComps_nLL, # container
    ObsArr = ObsFishAgeComps, # observed age compositions
    ExpArr = CAA, # predicted catch at age
    UseArr = UseFishAgeComps, # cells that are fit
    ISSArr = ISS_FishAgeComps, # input sample size
    WtArr = Wt_FishAgeComps, # multinomial weight
    TypeMat = FishAgeComps_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_LikeType, # likelihood for each fleet
    lnThetaArr = ln_FishAge_theta, # dispersion, split compositions
    lnThetaAggVec = ln_FishAge_theta_agg, # dispersion, aggregated compositions
    LNcorrArr = FishAge_corr_pars, # correlation, split compositions
    LNcorrAggVec = FishAge_corr_pars_agg, # correlation, aggregated compositions
    age_or_len = 0, # age compositions
    n_model_bins = n_ages, # model bins
    comp_bins_spec = FishAgeComps_bins, # bins this data source is fit over
    AgeingErrorArr = AgeingError_fish, # ageing error
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    pop = FALSE, # regional data source
    addtocomp = addtocomp,
    comp_const_obs = comp_const_obs,
    do_internal_comp_osa = do_internal_comp_osa,
    tracked_discrete = ObsFishAgeComps_osa_discrete, # registered discrete observations
    tracked_continuous = ObsFishAgeComps_osa_continuous, # registered continuous observations
    seas_agg = FishAgeComps_seas_Type # whether this source is fit as a season total
  )

  if(fit_lengths == 1) {

    # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
    # looks the residuals up under that name, so each data source is registered here by name
    ObsFishLenComps_osa = pack_comp_source_osa(
      ObsArr = ObsFishLenComps, # observed length compositions
      ISSArr = ISS_FishLenComps, # input sample size
      WtArr = Wt_FishLenComps, # multinomial weight
      UseArr = UseFishLenComps, # cells that are fit
      TypeMat = FishLenComps_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_LikeType, # likelihood for each fleet
      comp_bins_spec = FishLenComps_bins, # bins this data source is fit over
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      n_pop = n_pop,
      pop = FALSE, # regional data source
      addtocomp = addtocomp,
      do_internal_comp_osa = do_internal_comp_osa
    )

    ObsFishLenComps_osa_discrete = ObsFishLenComps_osa$discrete # fleets on a discrete likelihood
    ObsFishLenComps_osa_continuous = ObsFishLenComps_osa$continuous # fleets on a continuous likelihood
    if(!is.null(ObsFishLenComps_osa_discrete)) ObsFishLenComps_osa_discrete = RTMB::OBS(ObsFishLenComps_osa_discrete)
    if(!is.null(ObsFishLenComps_osa_continuous)) ObsFishLenComps_osa_continuous = RTMB::OBS(ObsFishLenComps_osa_continuous)

    FishLenComps_nLL = get_comp_source_nLL(
      nLL_arr = FishLenComps_nLL, # container
      ObsArr = ObsFishLenComps, # observed length compositions
      ExpArr = CAL, # predicted catch at length
      UseArr = UseFishLenComps, # cells that are fit
      ISSArr = ISS_FishLenComps, # input sample size
      WtArr = Wt_FishLenComps, # multinomial weight
      TypeMat = FishLenComps_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_LikeType, # likelihood for each fleet
      lnThetaArr = ln_FishLen_theta, # dispersion, split compositions
      lnThetaAggVec = ln_FishLen_theta_agg, # dispersion, aggregated compositions
      LNcorrArr = FishLen_corr_pars, # correlation, split compositions
      LNcorrAggVec = FishLen_corr_pars_agg, # correlation, aggregated compositions
      age_or_len = 1, # length compositions
      n_model_bins = n_lens, # model bins
      comp_bins_spec = FishLenComps_bins, # bins this data source is fit over
      LenBinMap_fn = LenBinMap_fn, # model bin to observed bin map
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      pop = FALSE, # regional data source
      addtocomp = addtocomp,
      comp_const_obs = comp_const_obs,
      do_internal_comp_osa = do_internal_comp_osa,
      tracked_discrete = ObsFishLenComps_osa_discrete, # registered discrete observations
      tracked_continuous = ObsFishLenComps_osa_continuous, # registered continuous observations
      seas_agg = FishLenComps_seas_Type # whether this source is fit as a season total
    )

  } # end if fitting lengths

  ### Retained Fishery Compositions (Population-Specific) --------------------
  # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
  # looks the residuals up under that name, so each data source is registered here by name
  ObsFishAgeComps_pop_osa = pack_comp_source_osa(
    ObsArr = ObsFishAgeComps_pop, # observed age compositions
    ISSArr = ISS_FishAgeComps_pop, # input sample size
    WtArr = Wt_FishAgeComps_pop, # multinomial weight
    UseArr = UseFishAgeComps_pop, # cells that are fit
    TypeMat = FishAgeComps_pop_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_pop_LikeType, # likelihood for each fleet
    comp_bins_spec = FishAgeComps_pop_bins, # bins this data source is fit over
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    n_pop = n_pop,
    pop = TRUE, # population-specific data source
    addtocomp = addtocomp,
    do_internal_comp_osa = do_internal_comp_osa
  )

  ObsFishAgeComps_pop_osa_discrete = ObsFishAgeComps_pop_osa$discrete # fleets on a discrete likelihood
  ObsFishAgeComps_pop_osa_continuous = ObsFishAgeComps_pop_osa$continuous # fleets on a continuous likelihood
  if(!is.null(ObsFishAgeComps_pop_osa_discrete)) ObsFishAgeComps_pop_osa_discrete = RTMB::OBS(ObsFishAgeComps_pop_osa_discrete)
  if(!is.null(ObsFishAgeComps_pop_osa_continuous)) ObsFishAgeComps_pop_osa_continuous = RTMB::OBS(ObsFishAgeComps_pop_osa_continuous)

  FishAgeComps_pop_nLL = get_comp_source_nLL(
    nLL_arr = FishAgeComps_pop_nLL, # container
    ObsArr = ObsFishAgeComps_pop, # observed age compositions
    ExpArr = CAA, # predicted catch at age
    UseArr = UseFishAgeComps_pop, # cells that are fit
    ISSArr = ISS_FishAgeComps_pop, # input sample size
    WtArr = Wt_FishAgeComps_pop, # multinomial weight
    TypeMat = FishAgeComps_pop_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_pop_LikeType, # likelihood for each fleet
    lnThetaArr = ln_FishAge_pop_theta, # dispersion, split compositions
    lnThetaAggVec = ln_FishAge_pop_theta_agg, # dispersion, aggregated compositions
    LNcorrArr = FishAge_pop_corr_pars, # correlation, split compositions
    LNcorrAggVec = FishAge_pop_corr_pars_agg, # correlation, aggregated compositions
    age_or_len = 0, # age compositions
    n_model_bins = n_ages, # model bins
    comp_bins_spec = FishAgeComps_pop_bins, # bins this data source is fit over
    AgeingErrorArr = AgeingError_fish, # ageing error
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    pop = TRUE, # population-specific data source
    addtocomp = addtocomp,
    comp_const_obs = comp_const_obs,
    do_internal_comp_osa = do_internal_comp_osa,
    tracked_discrete = ObsFishAgeComps_pop_osa_discrete, # registered discrete observations
    tracked_continuous = ObsFishAgeComps_pop_osa_continuous, # registered continuous observations
    seas_agg = FishAgeComps_pop_seas_Type # whether this source is fit as a season total
  )

  if(fit_lengths == 1) {

    # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
    # looks the residuals up under that name, so each data source is registered here by name
    ObsFishLenComps_pop_osa = pack_comp_source_osa(
      ObsArr = ObsFishLenComps_pop, # observed length compositions
      ISSArr = ISS_FishLenComps_pop, # input sample size
      WtArr = Wt_FishLenComps_pop, # multinomial weight
      UseArr = UseFishLenComps_pop, # cells that are fit
      TypeMat = FishLenComps_pop_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_pop_LikeType, # likelihood for each fleet
      comp_bins_spec = FishLenComps_pop_bins, # bins this data source is fit over
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      n_pop = n_pop,
      pop = TRUE, # population-specific data source
      addtocomp = addtocomp,
      do_internal_comp_osa = do_internal_comp_osa
    )

    ObsFishLenComps_pop_osa_discrete = ObsFishLenComps_pop_osa$discrete # fleets on a discrete likelihood
    ObsFishLenComps_pop_osa_continuous = ObsFishLenComps_pop_osa$continuous # fleets on a continuous likelihood
    if(!is.null(ObsFishLenComps_pop_osa_discrete)) ObsFishLenComps_pop_osa_discrete = RTMB::OBS(ObsFishLenComps_pop_osa_discrete)
    if(!is.null(ObsFishLenComps_pop_osa_continuous)) ObsFishLenComps_pop_osa_continuous = RTMB::OBS(ObsFishLenComps_pop_osa_continuous)

    FishLenComps_pop_nLL = get_comp_source_nLL(
      nLL_arr = FishLenComps_pop_nLL, # container
      ObsArr = ObsFishLenComps_pop, # observed length compositions
      ExpArr = CAL, # predicted catch at length
      UseArr = UseFishLenComps_pop, # cells that are fit
      ISSArr = ISS_FishLenComps_pop, # input sample size
      WtArr = Wt_FishLenComps_pop, # multinomial weight
      TypeMat = FishLenComps_pop_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_pop_LikeType, # likelihood for each fleet
      lnThetaArr = ln_FishLen_pop_theta, # dispersion, split compositions
      lnThetaAggVec = ln_FishLen_pop_theta_agg, # dispersion, aggregated compositions
      LNcorrArr = FishLen_pop_corr_pars, # correlation, split compositions
      LNcorrAggVec = FishLen_pop_corr_pars_agg, # correlation, aggregated compositions
      age_or_len = 1, # length compositions
      n_model_bins = n_lens, # model bins
      comp_bins_spec = FishLenComps_pop_bins, # bins this data source is fit over
      LenBinMap_fn = LenBinMap_fn, # model bin to observed bin map
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      pop = TRUE, # population-specific data source
      addtocomp = addtocomp,
      comp_const_obs = comp_const_obs,
      do_internal_comp_osa = do_internal_comp_osa,
      tracked_discrete = ObsFishLenComps_pop_osa_discrete, # registered discrete observations
      tracked_continuous = ObsFishLenComps_pop_osa_continuous, # registered continuous observations
      seas_agg = FishLenComps_pop_seas_Type # whether this source is fit as a season total
    )

  } # end if fitting lengths

  ### Discarded Fishery Compositions (Region-Specific) -----------------------
  # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
  # looks the residuals up under that name, so each data source is registered here by name
  ObsFishAgeComps_discard_osa = pack_comp_source_osa(
    ObsArr = ObsFishAgeComps_discard, # observed age compositions
    ISSArr = ISS_FishAgeComps_discard, # input sample size
    WtArr = Wt_FishAgeComps_discard, # multinomial weight
    UseArr = UseFishAgeComps_discard, # cells that are fit
    TypeMat = FishAgeComps_discard_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_discard_LikeType, # likelihood for each fleet
    comp_bins_spec = FishAgeComps_discard_bins, # bins this data source is fit over
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    n_pop = n_pop,
    pop = FALSE, # regional data source
    addtocomp = addtocomp,
    do_internal_comp_osa = do_internal_comp_osa
  )

  ObsFishAgeComps_discard_osa_discrete = ObsFishAgeComps_discard_osa$discrete # fleets on a discrete likelihood
  ObsFishAgeComps_discard_osa_continuous = ObsFishAgeComps_discard_osa$continuous # fleets on a continuous likelihood
  if(!is.null(ObsFishAgeComps_discard_osa_discrete)) ObsFishAgeComps_discard_osa_discrete = RTMB::OBS(ObsFishAgeComps_discard_osa_discrete)
  if(!is.null(ObsFishAgeComps_discard_osa_continuous)) ObsFishAgeComps_discard_osa_continuous = RTMB::OBS(ObsFishAgeComps_discard_osa_continuous)

  FishAgeComps_discard_nLL = get_comp_source_nLL(
    nLL_arr = FishAgeComps_discard_nLL, # container
    ObsArr = ObsFishAgeComps_discard, # observed age compositions
    ExpArr = DAA, # predicted discards at age
    UseArr = UseFishAgeComps_discard, # cells that are fit
    ISSArr = ISS_FishAgeComps_discard, # input sample size
    WtArr = Wt_FishAgeComps_discard, # multinomial weight
    TypeMat = FishAgeComps_discard_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_discard_LikeType, # likelihood for each fleet
    lnThetaArr = ln_FishAge_discard_theta, # dispersion, split compositions
    lnThetaAggVec = ln_FishAge_discard_theta_agg, # dispersion, aggregated compositions
    LNcorrArr = FishAge_discard_corr_pars, # correlation, split compositions
    LNcorrAggVec = FishAge_discard_corr_pars_agg, # correlation, aggregated compositions
    age_or_len = 0, # age compositions
    n_model_bins = n_ages, # model bins
    comp_bins_spec = FishAgeComps_discard_bins, # bins this data source is fit over
    AgeingErrorArr = AgeingError_fish, # ageing error
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    pop = FALSE, # regional data source
    addtocomp = addtocomp,
    comp_const_obs = comp_const_obs,
    do_internal_comp_osa = do_internal_comp_osa,
    tracked_discrete = ObsFishAgeComps_discard_osa_discrete, # registered discrete observations
    tracked_continuous = ObsFishAgeComps_discard_osa_continuous, # registered continuous observations
    seas_agg = FishAgeComps_discard_seas_Type # whether this source is fit as a season total
  )

  if(fit_lengths == 1) {

    # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
    # looks the residuals up under that name, so each data source is registered here by name
    ObsFishLenComps_discard_osa = pack_comp_source_osa(
      ObsArr = ObsFishLenComps_discard, # observed length compositions
      ISSArr = ISS_FishLenComps_discard, # input sample size
      WtArr = Wt_FishLenComps_discard, # multinomial weight
      UseArr = UseFishLenComps_discard, # cells that are fit
      TypeMat = FishLenComps_discard_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_discard_LikeType, # likelihood for each fleet
      comp_bins_spec = FishLenComps_discard_bins, # bins this data source is fit over
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      n_pop = n_pop,
      pop = FALSE, # regional data source
      addtocomp = addtocomp,
      do_internal_comp_osa = do_internal_comp_osa
    )

    ObsFishLenComps_discard_osa_discrete = ObsFishLenComps_discard_osa$discrete # fleets on a discrete likelihood
    ObsFishLenComps_discard_osa_continuous = ObsFishLenComps_discard_osa$continuous # fleets on a continuous likelihood
    if(!is.null(ObsFishLenComps_discard_osa_discrete)) ObsFishLenComps_discard_osa_discrete = RTMB::OBS(ObsFishLenComps_discard_osa_discrete)
    if(!is.null(ObsFishLenComps_discard_osa_continuous)) ObsFishLenComps_discard_osa_continuous = RTMB::OBS(ObsFishLenComps_discard_osa_continuous)

    FishLenComps_discard_nLL = get_comp_source_nLL(
      nLL_arr = FishLenComps_discard_nLL, # container
      ObsArr = ObsFishLenComps_discard, # observed length compositions
      ExpArr = DAL, # predicted discards at length
      UseArr = UseFishLenComps_discard, # cells that are fit
      ISSArr = ISS_FishLenComps_discard, # input sample size
      WtArr = Wt_FishLenComps_discard, # multinomial weight
      TypeMat = FishLenComps_discard_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_discard_LikeType, # likelihood for each fleet
      lnThetaArr = ln_FishLen_discard_theta, # dispersion, split compositions
      lnThetaAggVec = ln_FishLen_discard_theta_agg, # dispersion, aggregated compositions
      LNcorrArr = FishLen_discard_corr_pars, # correlation, split compositions
      LNcorrAggVec = FishLen_discard_corr_pars_agg, # correlation, aggregated compositions
      age_or_len = 1, # length compositions
      n_model_bins = n_lens, # model bins
      comp_bins_spec = FishLenComps_discard_bins, # bins this data source is fit over
      LenBinMap_fn = LenBinMap_fn, # model bin to observed bin map
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      pop = FALSE, # regional data source
      addtocomp = addtocomp,
      comp_const_obs = comp_const_obs,
      do_internal_comp_osa = do_internal_comp_osa,
      tracked_discrete = ObsFishLenComps_discard_osa_discrete, # registered discrete observations
      tracked_continuous = ObsFishLenComps_discard_osa_continuous, # registered continuous observations
      seas_agg = FishLenComps_discard_seas_Type # whether this source is fit as a season total
    )

  } # end if fitting lengths

  ### Discarded Fishery Compositions (Population-Specific) -------------------
  # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
  # looks the residuals up under that name, so each data source is registered here by name
  ObsFishAgeComps_discard_pop_osa = pack_comp_source_osa(
    ObsArr = ObsFishAgeComps_discard_pop, # observed age compositions
    ISSArr = ISS_FishAgeComps_discard_pop, # input sample size
    WtArr = Wt_FishAgeComps_discard_pop, # multinomial weight
    UseArr = UseFishAgeComps_discard_pop, # cells that are fit
    TypeMat = FishAgeComps_discard_pop_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_discard_pop_LikeType, # likelihood for each fleet
    comp_bins_spec = FishAgeComps_discard_pop_bins, # bins this data source is fit over
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    n_pop = n_pop,
    pop = TRUE, # population-specific data source
    addtocomp = addtocomp,
    do_internal_comp_osa = do_internal_comp_osa
  )

  ObsFishAgeComps_discard_pop_osa_discrete = ObsFishAgeComps_discard_pop_osa$discrete # fleets on a discrete likelihood
  ObsFishAgeComps_discard_pop_osa_continuous = ObsFishAgeComps_discard_pop_osa$continuous # fleets on a continuous likelihood
  if(!is.null(ObsFishAgeComps_discard_pop_osa_discrete)) ObsFishAgeComps_discard_pop_osa_discrete = RTMB::OBS(ObsFishAgeComps_discard_pop_osa_discrete)
  if(!is.null(ObsFishAgeComps_discard_pop_osa_continuous)) ObsFishAgeComps_discard_pop_osa_continuous = RTMB::OBS(ObsFishAgeComps_discard_pop_osa_continuous)

  FishAgeComps_discard_pop_nLL = get_comp_source_nLL(
    nLL_arr = FishAgeComps_discard_pop_nLL, # container
    ObsArr = ObsFishAgeComps_discard_pop, # observed age compositions
    ExpArr = DAA, # predicted discards at age
    UseArr = UseFishAgeComps_discard_pop, # cells that are fit
    ISSArr = ISS_FishAgeComps_discard_pop, # input sample size
    WtArr = Wt_FishAgeComps_discard_pop, # multinomial weight
    TypeMat = FishAgeComps_discard_pop_Type, # how the compositions are split
    LikeTypeVec = FishAgeComps_discard_pop_LikeType, # likelihood for each fleet
    lnThetaArr = ln_FishAge_discard_pop_theta, # dispersion, split compositions
    lnThetaAggVec = ln_FishAge_discard_pop_theta_agg, # dispersion, aggregated compositions
    LNcorrArr = FishAge_discard_pop_corr_pars, # correlation, split compositions
    LNcorrAggVec = FishAge_discard_pop_corr_pars_agg, # correlation, aggregated compositions
    age_or_len = 0, # age compositions
    n_model_bins = n_ages, # model bins
    comp_bins_spec = FishAgeComps_discard_pop_bins, # bins this data source is fit over
    AgeingErrorArr = AgeingError_fish, # ageing error
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_fish_fleets,
    n_sexes = n_sexes,
    pop = TRUE, # population-specific data source
    addtocomp = addtocomp,
    comp_const_obs = comp_const_obs,
    do_internal_comp_osa = do_internal_comp_osa,
    tracked_discrete = ObsFishAgeComps_discard_pop_osa_discrete, # registered discrete observations
    tracked_continuous = ObsFishAgeComps_discard_pop_osa_continuous, # registered continuous observations
    seas_agg = FishAgeComps_discard_pop_seas_Type # whether this source is fit as a season total
  )

  if(fit_lengths == 1) {

    # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
    # looks the residuals up under that name, so each data source is registered here by name
    ObsFishLenComps_discard_pop_osa = pack_comp_source_osa(
      ObsArr = ObsFishLenComps_discard_pop, # observed length compositions
      ISSArr = ISS_FishLenComps_discard_pop, # input sample size
      WtArr = Wt_FishLenComps_discard_pop, # multinomial weight
      UseArr = UseFishLenComps_discard_pop, # cells that are fit
      TypeMat = FishLenComps_discard_pop_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_discard_pop_LikeType, # likelihood for each fleet
      comp_bins_spec = FishLenComps_discard_pop_bins, # bins this data source is fit over
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      n_pop = n_pop,
      pop = TRUE, # population-specific data source
      addtocomp = addtocomp,
      do_internal_comp_osa = do_internal_comp_osa
    )

    ObsFishLenComps_discard_pop_osa_discrete = ObsFishLenComps_discard_pop_osa$discrete # fleets on a discrete likelihood
    ObsFishLenComps_discard_pop_osa_continuous = ObsFishLenComps_discard_pop_osa$continuous # fleets on a continuous likelihood
    if(!is.null(ObsFishLenComps_discard_pop_osa_discrete)) ObsFishLenComps_discard_pop_osa_discrete = RTMB::OBS(ObsFishLenComps_discard_pop_osa_discrete)
    if(!is.null(ObsFishLenComps_discard_pop_osa_continuous)) ObsFishLenComps_discard_pop_osa_continuous = RTMB::OBS(ObsFishLenComps_discard_pop_osa_continuous)

    FishLenComps_discard_pop_nLL = get_comp_source_nLL(
      nLL_arr = FishLenComps_discard_pop_nLL, # container
      ObsArr = ObsFishLenComps_discard_pop, # observed length compositions
      ExpArr = DAL, # predicted discards at length
      UseArr = UseFishLenComps_discard_pop, # cells that are fit
      ISSArr = ISS_FishLenComps_discard_pop, # input sample size
      WtArr = Wt_FishLenComps_discard_pop, # multinomial weight
      TypeMat = FishLenComps_discard_pop_Type, # how the compositions are split
      LikeTypeVec = FishLenComps_discard_pop_LikeType, # likelihood for each fleet
      lnThetaArr = ln_FishLen_discard_pop_theta, # dispersion, split compositions
      lnThetaAggVec = ln_FishLen_discard_pop_theta_agg, # dispersion, aggregated compositions
      LNcorrArr = FishLen_discard_pop_corr_pars, # correlation, split compositions
      LNcorrAggVec = FishLen_discard_pop_corr_pars_agg, # correlation, aggregated compositions
      age_or_len = 1, # length compositions
      n_model_bins = n_lens, # model bins
      comp_bins_spec = FishLenComps_discard_pop_bins, # bins this data source is fit over
      LenBinMap_fn = LenBinMap_fn, # model bin to observed bin map
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_fish_fleets,
      n_sexes = n_sexes,
      pop = TRUE, # population-specific data source
      addtocomp = addtocomp,
      comp_const_obs = comp_const_obs,
      do_internal_comp_osa = do_internal_comp_osa,
      tracked_discrete = ObsFishLenComps_discard_pop_osa_discrete, # registered discrete observations
      tracked_continuous = ObsFishLenComps_discard_pop_osa_continuous, # registered continuous observations
      seas_agg = FishLenComps_discard_pop_seas_Type # whether this source is fit as a season total
    )

  } # end if fitting lengths
  ## Survey Likelihoods ------------------------------------------------------
  ### Survey Indices (Regional) ----------------------------------------------
  # reads ObsSrvIdx as an array. the lognormal block below rebinds that name to a flat
  # log-scale vector, since OSA files residuals under the name they were registered on
  SrvIdx_nLL = get_index_regional_nLL(
    nLL_arr = SrvIdx_nLL, # container
    Use = UseSrvIdx, # cells that are fit
    Obs = ObsSrvIdx, # observed index
    Pred = PredSrvIdx, # predicted index
    SD = SrvIdx_SD, # index standard deviation
    LikeType = SrvIdx_LikeType, # likelihood for each fleet
    Cov = SrvIdx_Cov, # covariance, read only by a multivariate normal
    seas_Type = SrvIdx_seas_Type, # whether each fleet is fit as a season total
    const = addtosrvidx, # constant added inside the log
    n_fleets = n_srv_fleets
  )

  # only lognormal fleets register OSA observations, so drop the fleets already done above
  UseSrvIdx_ln = UseSrvIdx
  for(f in 1:n_srv_fleets) if(SrvIdx_LikeType[f] != 0) UseSrvIdx_ln[,,,f] = 0

  if(any(UseSrvIdx_ln == 1)) { # setup OSA

    valid_idx_sr = which(UseSrvIdx_ln == 1) # position of each fitted cell
    ObsSrvIdx_map = arrayInd(valid_idx_sr, dim(UseSrvIdx_ln)) # region, year, season and fleet of each
    ObsSrvIdx = log(ObsSrvIdx[valid_idx_sr] + addtosrvidx) # observed index, log scale
    ObsSrvIdx = RTMB::OBS(ObsSrvIdx) # register for OSA residuals

    SrvIdx_nLL = eval_index_osa_nLL(
      nLL_arr = SrvIdx_nLL, # container
      obs_vec = ObsSrvIdx, # registered observations, log scale
      obs_map = ObsSrvIdx_map, # index of each observation
      Pred = PredSrvIdx, # predicted index
      SD = SrvIdx_SD, # index standard deviation
      seas_Type = SrvIdx_seas_Type, # whether each fleet is fit as a season total
      const = addtosrvidx, # constant added inside the log
      pop = FALSE
    )
  }

  ### Survey Index at Age (Regional) -----------------------------------------
  ObsSrvIdxAA = prep_at_age_obs(ObsSrvIdxAA, UseSrvIdxAA, SrvIdxAA_LikeType, addtosrvidx)
  if(length(ObsSrvIdxAA)) ObsSrvIdxAA = RTMB::OBS(ObsSrvIdxAA) # setup OSA

  siaa = get_at_age_source_nLL(
    obs_t = ObsSrvIdxAA, # observed survey index at age
    use = UseSrvIdxAA, # cells that are fit
    ln_sigma = ln_sigmaSrvIdxAA, # observation error, log scale
    source = "srv_index", # prediction is summed from SrvIAA
    pop = FALSE, # regional data source
    arrays = at_age_arrays, # model predictions
    obs_se = ObsSrvIdxAA_SE, # reported standard errors
    sd_form = SrvIdxAA_sigma_form, # how those standard errors enter
    like_type = SrvIdxAA_LikeType, # lognormal or normal
    const = addtosrvidx, # constant added inside the log
    corr_type = AgeObsCorr_srv_idx, # correlation form across ages
    trans_rho = trans_rho_srv_idx, # correlation across ages, unconstrained
    trans_rho_year = trans_rho_srv_idx_year, # correlation across years, unconstrained
    us_pars = trans_rho_srv_idx_us, # unstructured correlation parameters
    aa_type = SrvIdxAA_Type, # dims the observations are split by
    seas_agg = SrvIdxAA_seas_Type # whether this source is fit as a season total
  )

  SrvIdxAA_nLL = siaa$nLL # nLL
  PredSrvIdxAA = siaa$pred # predicted survey index at age


  ### Survey Indices (Population-Specific) -----------------------------------
  # reads ObsSrvIdx_pop as an array, which the lognormal block below rebinds to a flat vector.
  # a multivariate normal covariance describes the regional series, so mvn fleets go there too
  SrvIdx_pop_nLL = get_index_pop_nLL(
    nLL_arr = SrvIdx_pop_nLL, # container
    Use = UseSrvIdx_pop, # cells that are fit
    Obs = ObsSrvIdx_pop, # observed index
    Pred = PredSrvIdx, # predicted index
    SD = SrvIdx_pop_SD, # index standard deviation
    LikeType = SrvIdx_LikeType, # likelihood for each fleet
    seas_Type = SrvIdx_pop_seas_Type, # whether each fleet is fit as a season total
    n_fleets = n_srv_fleets
  )

  # only lognormal fleets register OSA observations, so drop the fleets already done above
  UseSrvIdx_pop_ln = UseSrvIdx_pop
  for(f in 1:n_srv_fleets) if(SrvIdx_LikeType[f] == 1) UseSrvIdx_pop_ln[,,,,f] = 0

  if(any(UseSrvIdx_pop_ln == 1)) { # setup OSA

    valid_idx_sp = which(UseSrvIdx_pop_ln == 1) # position of each fitted cell
    ObsSrvIdx_pop_map = arrayInd(valid_idx_sp, dim(UseSrvIdx_pop_ln)) # population, region, year, season and fleet of each
    ObsSrvIdx_pop = log(ObsSrvIdx_pop[valid_idx_sp] + addtosrvidx) # observed index, log scale
    ObsSrvIdx_pop = RTMB::OBS(ObsSrvIdx_pop) # register for OSA residuals

    SrvIdx_pop_nLL = eval_index_osa_nLL(
      nLL_arr = SrvIdx_pop_nLL, # container
      obs_vec = ObsSrvIdx_pop, # registered observations, log scale
      obs_map = ObsSrvIdx_pop_map, # index of each observation
      Pred = PredSrvIdx, # predicted index
      SD = SrvIdx_pop_SD, # index standard deviation
      seas_Type = SrvIdx_pop_seas_Type, # whether each fleet is fit as a season total
      const = addtosrvidx, # constant added inside the log
      pop = TRUE
    )
  }

  ### Survey Index at Age (Population-Specific) ------------------------------
  ObsSrvIdxAA_pop = prep_at_age_obs(ObsSrvIdxAA_pop, UseSrvIdxAA_pop, SrvIdxAA_pop_LikeType, addtosrvidx)
  if(length(ObsSrvIdxAA_pop)) ObsSrvIdxAA_pop = RTMB::OBS(ObsSrvIdxAA_pop) # setup OSA

  siaa_pop = get_at_age_source_nLL(
    obs_t = ObsSrvIdxAA_pop, # observed survey index at age
    use = UseSrvIdxAA_pop, # cells that are fit
    ln_sigma = ln_sigmaSrvIdxAA_pop, # observation error, log scale
    source = "srv_index", # prediction is summed from SrvIAA
    pop = TRUE, # population-specific data source, never summed across populations
    arrays = at_age_arrays, # model predictions
    obs_se = ObsSrvIdxAA_pop_SE, # reported standard errors
    sd_form = SrvIdxAA_pop_sigma_form, # how those standard errors enter
    like_type = SrvIdxAA_pop_LikeType, # lognormal or normal
    const = addtosrvidx, # constant added inside the log
    corr_type = AgeObsCorr_srv_idx_pop, # correlation form across ages
    trans_rho = trans_rho_srv_idx_pop, # correlation across ages, unconstrained
    trans_rho_year = trans_rho_srv_idx_pop_year, # correlation across years, unconstrained
    us_pars = trans_rho_srv_idx_pop_us, # unstructured correlation parameters
    aa_type = SrvIdxAA_pop_Type, # dims the observations are split by
    seas_agg = SrvIdxAA_pop_seas_Type # whether this source is fit as a season total
  )

  SrvIdxAA_pop_nLL = siaa_pop$nLL # nLL
  PredSrvIdxAA_pop = siaa_pop$pred # predicted survey index at age


  ### Survey Compositions (Region-Specific) ----------------------------------
  # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
  # looks the residuals up under that name, so each data source is registered here by name
  ObsSrvAgeComps_osa = pack_comp_source_osa(
    ObsArr = ObsSrvAgeComps, # observed age compositions
    ISSArr = ISS_SrvAgeComps, # input sample size
    WtArr = Wt_SrvAgeComps, # multinomial weight
    UseArr = UseSrvAgeComps, # cells that are fit
    TypeMat = SrvAgeComps_Type, # how the compositions are split
    LikeTypeVec = SrvAgeComps_LikeType, # likelihood for each fleet
    comp_bins_spec = SrvAgeComps_bins, # bins this data source is fit over
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_srv_fleets,
    n_sexes = n_sexes,
    n_pop = n_pop,
    pop = FALSE, # regional data source
    addtocomp = addtocomp,
    do_internal_comp_osa = do_internal_comp_osa
  )

  ObsSrvAgeComps_osa_discrete = ObsSrvAgeComps_osa$discrete # fleets on a discrete likelihood
  ObsSrvAgeComps_osa_continuous = ObsSrvAgeComps_osa$continuous # fleets on a continuous likelihood
  if(!is.null(ObsSrvAgeComps_osa_discrete)) ObsSrvAgeComps_osa_discrete = RTMB::OBS(ObsSrvAgeComps_osa_discrete)
  if(!is.null(ObsSrvAgeComps_osa_continuous)) ObsSrvAgeComps_osa_continuous = RTMB::OBS(ObsSrvAgeComps_osa_continuous)

  SrvAgeComps_nLL = get_comp_source_nLL(
    nLL_arr = SrvAgeComps_nLL, # container
    ObsArr = ObsSrvAgeComps, # observed age compositions
    ExpArr = SrvIAA, # predicted survey numbers at age
    UseArr = UseSrvAgeComps, # cells that are fit
    ISSArr = ISS_SrvAgeComps, # input sample size
    WtArr = Wt_SrvAgeComps, # multinomial weight
    TypeMat = SrvAgeComps_Type, # how the compositions are split
    LikeTypeVec = SrvAgeComps_LikeType, # likelihood for each fleet
    lnThetaArr = ln_SrvAge_theta, # dispersion, split compositions
    lnThetaAggVec = ln_SrvAge_theta_agg, # dispersion, aggregated compositions
    LNcorrArr = SrvAge_corr_pars, # correlation, split compositions
    LNcorrAggVec = SrvAge_corr_pars_agg, # correlation, aggregated compositions
    age_or_len = 0, # age compositions
    n_model_bins = n_ages, # model bins
    comp_bins_spec = SrvAgeComps_bins, # bins this data source is fit over
      AgeingErrorArr = AgeingError_srv, # ageing error
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_srv_fleets,
    n_sexes = n_sexes,
    pop = FALSE, # regional data source
    addtocomp = addtocomp,
    comp_const_obs = comp_const_obs,
    do_internal_comp_osa = do_internal_comp_osa,
    tracked_discrete = ObsSrvAgeComps_osa_discrete, # registered discrete observations
    tracked_continuous = ObsSrvAgeComps_osa_continuous, # registered continuous observations
    seas_agg = SrvAgeComps_seas_Type # whether this source is fit as a season total
  )

  if(fit_lengths == 1) {

    # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
    # looks the residuals up under that name, so each data source is registered here by name
    ObsSrvLenComps_osa = pack_comp_source_osa(
      ObsArr = ObsSrvLenComps, # observed length compositions
      ISSArr = ISS_SrvLenComps, # input sample size
      WtArr = Wt_SrvLenComps, # multinomial weight
      UseArr = UseSrvLenComps, # cells that are fit
      TypeMat = SrvLenComps_Type, # how the compositions are split
      LikeTypeVec = SrvLenComps_LikeType, # likelihood for each fleet
      comp_bins_spec = SrvLenComps_bins, # bins this data source is fit over
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_srv_fleets,
      n_sexes = n_sexes,
      n_pop = n_pop,
      pop = FALSE, # regional data source
      addtocomp = addtocomp,
      do_internal_comp_osa = do_internal_comp_osa
    )

    ObsSrvLenComps_osa_discrete = ObsSrvLenComps_osa$discrete # fleets on a discrete likelihood
    ObsSrvLenComps_osa_continuous = ObsSrvLenComps_osa$continuous # fleets on a continuous likelihood
    if(!is.null(ObsSrvLenComps_osa_discrete)) ObsSrvLenComps_osa_discrete = RTMB::OBS(ObsSrvLenComps_osa_discrete)
    if(!is.null(ObsSrvLenComps_osa_continuous)) ObsSrvLenComps_osa_continuous = RTMB::OBS(ObsSrvLenComps_osa_continuous)

    SrvLenComps_nLL = get_comp_source_nLL(
      nLL_arr = SrvLenComps_nLL, # container
      ObsArr = ObsSrvLenComps, # observed length compositions
      ExpArr = SrvIAL, # predicted survey numbers at length
      UseArr = UseSrvLenComps, # cells that are fit
      ISSArr = ISS_SrvLenComps, # input sample size
      WtArr = Wt_SrvLenComps, # multinomial weight
      TypeMat = SrvLenComps_Type, # how the compositions are split
      LikeTypeVec = SrvLenComps_LikeType, # likelihood for each fleet
      lnThetaArr = ln_SrvLen_theta, # dispersion, split compositions
      lnThetaAggVec = ln_SrvLen_theta_agg, # dispersion, aggregated compositions
      LNcorrArr = SrvLen_corr_pars, # correlation, split compositions
      LNcorrAggVec = SrvLen_corr_pars_agg, # correlation, aggregated compositions
      age_or_len = 1, # length compositions
      n_model_bins = n_lens, # model bins
      comp_bins_spec = SrvLenComps_bins, # bins this data source is fit over
      LenBinMap_fn = LenBinMap_fn, # model bin to observed bin map
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_srv_fleets,
      n_sexes = n_sexes,
      pop = FALSE, # regional data source
      addtocomp = addtocomp,
      comp_const_obs = comp_const_obs,
      do_internal_comp_osa = do_internal_comp_osa,
      tracked_discrete = ObsSrvLenComps_osa_discrete, # registered discrete observations
      tracked_continuous = ObsSrvLenComps_osa_continuous, # registered continuous observations
      seas_agg = SrvLenComps_seas_Type # whether this source is fit as a season total
    )

  } # end if fitting lengths

  ### Survey Compositions (Population-Specific) ------------------------------
  # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
  # looks the residuals up under that name, so each data source is registered here by name
  ObsSrvAgeComps_pop_osa = pack_comp_source_osa(
    ObsArr = ObsSrvAgeComps_pop, # observed age compositions
    ISSArr = ISS_SrvAgeComps_pop, # input sample size
    WtArr = Wt_SrvAgeComps_pop, # multinomial weight
    UseArr = UseSrvAgeComps_pop, # cells that are fit
    TypeMat = SrvAgeComps_pop_Type, # how the compositions are split
    LikeTypeVec = SrvAgeComps_pop_LikeType, # likelihood for each fleet
    comp_bins_spec = SrvAgeComps_pop_bins, # bins this data source is fit over
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_srv_fleets,
    n_sexes = n_sexes,
    n_pop = n_pop,
    pop = TRUE, # population-specific data source
    addtocomp = addtocomp,
    do_internal_comp_osa = do_internal_comp_osa
  )

  ObsSrvAgeComps_pop_osa_discrete = ObsSrvAgeComps_pop_osa$discrete # fleets on a discrete likelihood
  ObsSrvAgeComps_pop_osa_continuous = ObsSrvAgeComps_pop_osa$continuous # fleets on a continuous likelihood
  if(!is.null(ObsSrvAgeComps_pop_osa_discrete)) ObsSrvAgeComps_pop_osa_discrete = RTMB::OBS(ObsSrvAgeComps_pop_osa_discrete)
  if(!is.null(ObsSrvAgeComps_pop_osa_continuous)) ObsSrvAgeComps_pop_osa_continuous = RTMB::OBS(ObsSrvAgeComps_pop_osa_continuous)

  SrvAgeComps_pop_nLL = get_comp_source_nLL(
    nLL_arr = SrvAgeComps_pop_nLL, # container
    ObsArr = ObsSrvAgeComps_pop, # observed age compositions
    ExpArr = SrvIAA, # predicted survey numbers at age
    UseArr = UseSrvAgeComps_pop, # cells that are fit
    ISSArr = ISS_SrvAgeComps_pop, # input sample size
    WtArr = Wt_SrvAgeComps_pop, # multinomial weight
    TypeMat = SrvAgeComps_pop_Type, # how the compositions are split
    LikeTypeVec = SrvAgeComps_pop_LikeType, # likelihood for each fleet
    lnThetaArr = ln_SrvAge_pop_theta, # dispersion, split compositions
    lnThetaAggVec = ln_SrvAge_pop_theta_agg, # dispersion, aggregated compositions
    LNcorrArr = SrvAge_pop_corr_pars, # correlation, split compositions
    LNcorrAggVec = SrvAge_pop_corr_pars_agg, # correlation, aggregated compositions
    age_or_len = 0, # age compositions
    n_model_bins = n_ages, # model bins
    comp_bins_spec = SrvAgeComps_pop_bins, # bins this data source is fit over
      AgeingErrorArr = AgeingError_srv, # ageing error
    n_pop = n_pop,
    n_regions = n_regions,
    n_yrs = n_yrs,
    n_seas = n_seas,
    n_fleets = n_srv_fleets,
    n_sexes = n_sexes,
    pop = TRUE, # population-specific data source
    addtocomp = addtocomp,
    comp_const_obs = comp_const_obs,
    do_internal_comp_osa = do_internal_comp_osa,
    tracked_discrete = ObsSrvAgeComps_pop_osa_discrete, # registered discrete observations
    tracked_continuous = ObsSrvAgeComps_pop_osa_continuous, # registered continuous observations
    seas_agg = SrvAgeComps_pop_seas_Type # whether this source is fit as a season total
  )

  if(fit_lengths == 1) {

    # RTMB::OBS takes an observation's name from the variable it is called on, and get_osa
    # looks the residuals up under that name, so each data source is registered here by name
    ObsSrvLenComps_pop_osa = pack_comp_source_osa(
      ObsArr = ObsSrvLenComps_pop, # observed length compositions
      ISSArr = ISS_SrvLenComps_pop, # input sample size
      WtArr = Wt_SrvLenComps_pop, # multinomial weight
      UseArr = UseSrvLenComps_pop, # cells that are fit
      TypeMat = SrvLenComps_pop_Type, # how the compositions are split
      LikeTypeVec = SrvLenComps_pop_LikeType, # likelihood for each fleet
      comp_bins_spec = SrvLenComps_pop_bins, # bins this data source is fit over
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_srv_fleets,
      n_sexes = n_sexes,
      n_pop = n_pop,
      pop = TRUE, # population-specific data source
      addtocomp = addtocomp,
      do_internal_comp_osa = do_internal_comp_osa
    )

    ObsSrvLenComps_pop_osa_discrete = ObsSrvLenComps_pop_osa$discrete # fleets on a discrete likelihood
    ObsSrvLenComps_pop_osa_continuous = ObsSrvLenComps_pop_osa$continuous # fleets on a continuous likelihood
    if(!is.null(ObsSrvLenComps_pop_osa_discrete)) ObsSrvLenComps_pop_osa_discrete = RTMB::OBS(ObsSrvLenComps_pop_osa_discrete)
    if(!is.null(ObsSrvLenComps_pop_osa_continuous)) ObsSrvLenComps_pop_osa_continuous = RTMB::OBS(ObsSrvLenComps_pop_osa_continuous)

    SrvLenComps_pop_nLL = get_comp_source_nLL(
      nLL_arr = SrvLenComps_pop_nLL, # container
      ObsArr = ObsSrvLenComps_pop, # observed length compositions
      ExpArr = SrvIAL, # predicted survey numbers at length
      UseArr = UseSrvLenComps_pop, # cells that are fit
      ISSArr = ISS_SrvLenComps_pop, # input sample size
      WtArr = Wt_SrvLenComps_pop, # multinomial weight
      TypeMat = SrvLenComps_pop_Type, # how the compositions are split
      LikeTypeVec = SrvLenComps_pop_LikeType, # likelihood for each fleet
      lnThetaArr = ln_SrvLen_pop_theta, # dispersion, split compositions
      lnThetaAggVec = ln_SrvLen_pop_theta_agg, # dispersion, aggregated compositions
      LNcorrArr = SrvLen_pop_corr_pars, # correlation, split compositions
      LNcorrAggVec = SrvLen_pop_corr_pars_agg, # correlation, aggregated compositions
      age_or_len = 1, # length compositions
      n_model_bins = n_lens, # model bins
      comp_bins_spec = SrvLenComps_pop_bins, # bins this data source is fit over
      LenBinMap_fn = LenBinMap_fn, # model bin to observed bin map
      n_pop = n_pop,
      n_regions = n_regions,
      n_yrs = n_yrs,
      n_seas = n_seas,
      n_fleets = n_srv_fleets,
      n_sexes = n_sexes,
      pop = TRUE, # population-specific data source
      addtocomp = addtocomp,
      comp_const_obs = comp_const_obs,
      do_internal_comp_osa = do_internal_comp_osa,
      tracked_discrete = ObsSrvLenComps_pop_osa_discrete, # registered discrete observations
      tracked_continuous = ObsSrvLenComps_pop_osa_continuous, # registered continuous observations
      seas_agg = SrvLenComps_pop_seas_Type # whether this source is fit as a season total
    )

  } # end if fitting lengths
  ## Conditional Age-at-Length Likelihoods -----------------------------------
  # population summing helpers for the conditional age-at-length data sources live in
  # model_lik_caal.R, alongside the likelihood that uses them
  caal_exp = function(arr, y, seas, f) caal_sum_pop(arr, y, seas, f, n_pop, n_regions, n_lens, n_ages, n_sexes)
  caal_exp_len = function(arr, y, seas, l, f) caal_sum_pop_len(arr, y, seas, l, f, n_pop, n_regions, n_ages, n_sexes)


  ### Fishery Conditional Age-at-Length --------------------------------------
  if(do_fish_caal) {
    if(do_internal_comp_osa == FALSE) {

      for(y in 1:n_yrs) {
        for(f in 1:n_fish_fleets) {
          for(seas in 1:n_seas) {
            if(sum(UseFish_caal[,y,seas,,f]) >= 1) {
              Fish_caal_nLL[,y,seas,,,f] = Get_CAAL_Likelihoods(
                comp_const_obs = comp_const_obs,
                Exp = caal_exp(Fish_caal, y, seas, f),
                Obs = ObsFish_caal[,y,seas,,,,f],
                ISS = ISS_Fish_caal[,y,seas,,,f],
                Wt_Mltnml = Wt_Fish_caal[,y,seas,,,f],
                ln_theta = ln_Fish_caal_theta[,,f],
                ln_theta_agg = ln_Fish_caal_theta_agg[f],
                Comp_Type = Fish_caal_Type[y,f],
                Likelihood_Type = Fish_caal_LikeType[f],
                n_regions = n_regions,
                n_lens = n_lens,
                n_model_bins = n_ages,
                n_obs_bins = dim(ObsFish_caal)[5],
                n_sexes = n_sexes,
                AgeingError = AgeingError_fish[y,,,f],
                use = UseFish_caal[,y,seas,,f],
                addtocomp = addtocomp,
                comp_bins = fleet_bins_or_null(Fish_caal_bins, f)
              )
            } # if we have fishery caal
          } # end seas loop
        } # end f loop
      } # end y loop

    } else {

      # Using OSA compositions
      ObsFish_caal_osa = pack_caal_osa(
        ObsArr = ObsFish_caal,
        ISSArr = ISS_Fish_caal,
        WtArr = Wt_Fish_caal,
        UseArr = UseFish_caal,
        TypeMat = Fish_caal_Type,
        LikeTypeVec = Fish_caal_LikeType,
        n_yrs = n_yrs,
        n_seas = n_seas,
        n_lens = n_lens,
        n_fleets = n_fish_fleets,
        n_sexes = n_sexes,
        addtocomp = addtocomp,
        BinsArr = bins_or_null(Fish_caal_bins)
      )

      if(!is.null(ObsFish_caal_osa)) {
        ObsFish_caal_osa = RTMB::OBS(ObsFish_caal_osa)
        Fish_caal_nLL = eval_caal_osa(
          nLL_arr = Fish_caal_nLL,
          tracked = ObsFish_caal_osa,
          ExpArrFn = function(y, seas, l, f) caal_exp_len(Fish_caal, y, seas, l, f),
          UseArr = UseFish_caal,
          TypeMat = Fish_caal_Type,
          LikeTypeVec = Fish_caal_LikeType,
          ISSArr = ISS_Fish_caal,
          lnThetaArr = ln_Fish_caal_theta,
          lnThetaAggVec = ln_Fish_caal_theta_agg,
          n_regions = n_regions,
          n_yrs = n_yrs,
          n_seas = n_seas,
          n_lens = n_lens,
          n_fleets = n_fish_fleets,
          n_sexes = n_sexes,
          n_model_bins = n_ages,
          n_obs_bins = dim(ObsFish_caal)[5],
          AgeingErrorFn = function(y, f) AgeingError_fish[y,,,f],
          addtocomp = addtocomp,
          BinsArr = bins_or_null(Fish_caal_bins)
        )
      }

    } # end osa switch
  } # end fishery caal

  ### Survey Conditional Age-at-Length ---------------------------------------
  if(do_srv_caal) {
    if(do_internal_comp_osa == FALSE) {

      for(y in 1:n_yrs) {
        for(sf in 1:n_srv_fleets) {
          for(seas in 1:n_seas) {
            if(sum(UseSrv_caal[,y,seas,,sf]) >= 1) {
              Srv_caal_nLL[,y,seas,,,sf] = Get_CAAL_Likelihoods(
                comp_const_obs = comp_const_obs,
                Exp = caal_exp(Srv_caal, y, seas, sf),
                Obs = ObsSrv_caal[,y,seas,,,,sf],
                ISS = ISS_Srv_caal[,y,seas,,,sf],
                Wt_Mltnml = Wt_Srv_caal[,y,seas,,,sf],
                ln_theta = ln_Srv_caal_theta[,,sf],
                ln_theta_agg = ln_Srv_caal_theta_agg[sf],
                Comp_Type = Srv_caal_Type[y,sf],
                Likelihood_Type = Srv_caal_LikeType[sf],
                n_regions = n_regions,
                n_lens = n_lens,
                n_model_bins = n_ages,
                n_obs_bins = dim(ObsSrv_caal)[5],
                n_sexes = n_sexes,
                AgeingError = AgeingError_srv[y,,,sf],
                use = UseSrv_caal[,y,seas,,sf],
                addtocomp = addtocomp,
                comp_bins = fleet_bins_or_null(Srv_caal_bins, sf)
              )
            } # if we have survey caal
          } # end seas loop
        } # end sf loop
      } # end y loop

    } else {

      # Using OSA compositions
      ObsSrv_caal_osa = pack_caal_osa(
        ObsArr = ObsSrv_caal,
        ISSArr = ISS_Srv_caal,
        WtArr = Wt_Srv_caal,
        UseArr = UseSrv_caal,
        TypeMat = Srv_caal_Type,
        LikeTypeVec = Srv_caal_LikeType,
        n_yrs = n_yrs,
        n_seas = n_seas,
        n_lens = n_lens,
        n_fleets = n_srv_fleets,
        n_sexes = n_sexes,
        addtocomp = addtocomp,
        BinsArr = bins_or_null(Srv_caal_bins)
      )

      if(!is.null(ObsSrv_caal_osa)) {
        ObsSrv_caal_osa = RTMB::OBS(ObsSrv_caal_osa)
        Srv_caal_nLL = eval_caal_osa(
          nLL_arr = Srv_caal_nLL,
          tracked = ObsSrv_caal_osa,
          ExpArrFn = function(y, seas, l, sf) caal_exp_len(Srv_caal, y, seas, l, sf),
          UseArr = UseSrv_caal,
          TypeMat = Srv_caal_Type,
          LikeTypeVec = Srv_caal_LikeType,
          ISSArr = ISS_Srv_caal,
          lnThetaArr = ln_Srv_caal_theta,
          lnThetaAggVec = ln_Srv_caal_theta_agg,
          n_regions = n_regions,
          n_yrs = n_yrs,
          n_seas = n_seas,
          n_lens = n_lens,
          n_fleets = n_srv_fleets,
          n_sexes = n_sexes,
          n_model_bins = n_ages,
          n_obs_bins = dim(ObsSrv_caal)[5],
          AgeingErrorFn = function(y, f) AgeingError_srv[y,,,f],
          addtocomp = addtocomp,
          BinsArr = bins_or_null(Srv_caal_bins)
        )
      }

    } # end osa switch
  } # end survey caal

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

  # Priors and Penalties ----------------------------------------------------
  ## Mortality Deviations ----------------------------------------------------
  ### Fishing Mortality (Penalty) --------------------------------------------
  if(Use_F_pen == 1) {
    Fmort_nLL = Get_Fdev_PE_loglik(
      PE_model = Fdev_model,
      ln_sigmaF = ln_sigmaF,
      Fdev_rho = Fdev_rho,
      ln_F_devs = ln_F_devs,
      map_ln_F_devs = map_ln_F_devs,
      Fdev_pen_center = Fdev_pen_center
    )
  } #  if using fishing mortality penalty

  ### Discard Mortality Rate (Penalty) ---------------------------------------
  if(Use_dmr_pen == 1) {
    dmr_nLL = get_dmr_penalty(
      logit_dmr_devs = logit_dmr_devs,
      ln_sigma_dmr = ln_sigma_dmr,
      map_logit_dmr_devs = map_logit_dmr_devs,
      n_fish_fleets = n_fish_fleets,
      n_yrs = n_yrs,
      n_regions = n_regions,
      n_seas = n_seas
    )
  } #  if using discard mortality rate penalty


  ## Selectivity -------------------------------------------------------------
  ### Selectivity (Penalty) --------------------------------------------------
  for(r in 1:n_regions) {

    for(f in 1:n_fish_fleets) {

      # Total Fishery Selectivity Deviations
      if(cont_tv_fish_sel[r,f] > 0 && fishsel_pe_wt[f] != 0) {

        sel_nLL = sel_nLL + - fishsel_pe_wt[f] * Get_PE_loglik(PE_model = cont_tv_fish_sel[r,f], # process error model
                                                PE_pars = fishsel_pe_pars[r,,,f, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_fishsel_devs[r,,,,f, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_fishsel_devs[r,,,,f, drop = FALSE],
                                                map_sel_devs_full = map_ln_fishsel_devs[,,,,f, drop = FALSE],
                                                min_sel_devs_shared_bins = fishsel_devs_min_shared_bins,
                                                rw_init_sigma = fishsel_rw_init_sigma[f]

        )
      } # end if

      # Retained Fishery Selectivity Deviations
      if(cont_tv_ret_sel[r,f] > 0 && retsel_pe_wt[f] != 0) {

        sel_nLL = sel_nLL + - retsel_pe_wt[f] * Get_PE_loglik(PE_model = cont_tv_ret_sel[r,f], # process error model
                                                PE_pars = retsel_pe_pars[r,,,f, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_retsel_devs[r,,,,f, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_retsel_devs[r,,,,f, drop = FALSE],
                                                map_sel_devs_full = map_ln_retsel_devs[,,,,f, drop = FALSE],
                                                min_sel_devs_shared_bins = retsel_devs_min_shared_bins,
                                                rw_init_sigma = retsel_rw_init_sigma[f]

        )
      } # end if
    } # end f loop

    # Survey Selectivity Deviations
    for(sf in 1:n_srv_fleets) {

      if(cont_tv_srv_sel[r,sf] > 0 && srvsel_pe_wt[sf] != 0) {

        sel_nLL = sel_nLL + - srvsel_pe_wt[sf] * Get_PE_loglik(PE_model = cont_tv_srv_sel[r,sf], # process error model
                                                PE_pars = srvsel_pe_pars[r,,,sf, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_srvsel_devs[r,,,,sf, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_srvsel_devs[r,,,,sf, drop = FALSE],
                                                map_sel_devs_full = map_ln_srvsel_devs[,,,,sf, drop = FALSE],
                                                min_sel_devs_shared_bins = srvsel_devs_min_shared_bins,
                                                rw_init_sigma = srvsel_rw_init_sigma[sf]

        )
      } # end if

    } # end sf loop
  } # end r loop


  ### Bin-Override Selectivity Deviations (Process Error) --------------------

  # Same process as above but overide bins and allow PE there
  for(r in 1:n_regions) {
    for(f in 1:n_fish_fleets) {
      if(cont_tv_fishsel_bin_devs[f] > 0) {
        sel_nLL = sel_nLL + - Get_PE_loglik(PE_model = cont_tv_fishsel_bin_devs[f],
                                            PE_pars = fishsel_bin_devs_pe_pars[r,,,f, drop = FALSE],
                                            ln_devs = ln_fishsel_bin_devs[r,,,,f, drop = FALSE],
                                            map_sel_devs = map_ln_fishsel_bin_devs[r,,,,f, drop = FALSE],
                                            map_sel_devs_full = map_ln_fishsel_bin_devs[,,,,f, drop = FALSE],
                                            min_sel_devs_shared_bins = 1:dim(ln_fishsel_bin_devs)[3],
                                            rw_init_sigma = fishsel_bin_devs_rw_init_sigma[f])
      } # end if
      if(cont_tv_retsel_bin_devs[f] > 0) {
        sel_nLL = sel_nLL + - Get_PE_loglik(PE_model = cont_tv_retsel_bin_devs[f],
                                             PE_pars = retsel_bin_devs_pe_pars[r,,,f, drop = FALSE],
                                             ln_devs = ln_retsel_bin_devs[r,,,,f, drop = FALSE],
                                             map_sel_devs = map_ln_retsel_bin_devs[r,,,,f, drop = FALSE],
                                            map_sel_devs_full = map_ln_retsel_bin_devs[,,,,f, drop = FALSE],
                                             min_sel_devs_shared_bins = 1:dim(ln_retsel_bin_devs)[3],
                                             rw_init_sigma = retsel_bin_devs_rw_init_sigma[f])
      } # end if
    } # end f loop
    for(sf in 1:n_srv_fleets) {
      if(cont_tv_srvsel_bin_devs[sf] > 0) {
        sel_nLL = sel_nLL + - Get_PE_loglik(PE_model = cont_tv_srvsel_bin_devs[sf],
                                             PE_pars = srvsel_bin_devs_pe_pars[r,,,sf, drop = FALSE],
                                             ln_devs = ln_srvsel_bin_devs[r,,,,sf, drop = FALSE],
                                             map_sel_devs = map_ln_srvsel_bin_devs[r,,,,sf, drop = FALSE],
                                            map_sel_devs_full = map_ln_srvsel_bin_devs[,,,,sf, drop = FALSE],
                                             min_sel_devs_shared_bins = 1:dim(ln_srvsel_bin_devs)[3],
                                             rw_init_sigma = srvsel_bin_devs_rw_init_sigma[sf])
      } # end if
    } # end sf loop
  } # end r loop

  ### Selectivity Smoothness (Penalty) ---------------------------------------
  smooth_pen_terms = c("smooth_bin_curve", "smooth_bin_diff", "smooth_yr_diff", "smooth_yr_curve", "smooth_dome", "smooth_mean_center")
  for(r in 1:n_regions) {

    for(f in 1:n_fish_fleets) {

      # If Bicubic spline
      bicubic_yrs = which(fish_sel_model[r,,f] == 8)
      has_bicubic = length(bicubic_yrs) > 0
      has_nonzero_pen = any(unlist(lapply(smooth_pen_terms, function(data_name) safe_extract(fish_sel_pen_wts[[f]], data_name))) != 0)

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
                                                         wt_bin_curve = safe_extract(fish_sel_pen_wts[[f]], "smooth_bin_curve"),
                                                         wt_bin_diff = safe_extract(fish_sel_pen_wts[[f]], "smooth_bin_diff"),
                                                         wt_yr_diff = safe_extract(fish_sel_pen_wts[[f]], "smooth_yr_diff"),
                                                         wt_yr_curve = safe_extract(fish_sel_pen_wts[[f]], "smooth_yr_curve"),
                                                         wt_dome = safe_extract(fish_sel_pen_wts[[f]], "smooth_dome"),
                                                         wt_mean_center = safe_extract(fish_sel_pen_wts[[f]], "smooth_mean_center"),
                                                         normalize = fish_sel_pen_wts[[f]]$normalize,
                                                         bin_range = fish_sel_pen_wts[[f]]$bin_range,
                                                         yr_diff_ref = fish_sel_pen_wts[[f]]$yr_diff_ref)
      } # end if
    } # end f loop

    for(f in 1:n_fish_fleets) {

      # If Bicubic spline (retention)
      bicubic_yrs = which(ret_sel_model[r,,f] == 8)
      has_bicubic = length(bicubic_yrs) > 0
      has_nonzero_pen = any(unlist(lapply(smooth_pen_terms, function(data_name) safe_extract(ret_sel_pen_wts[[f]], data_name))) != 0)

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
                                                         wt_bin_curve = safe_extract(ret_sel_pen_wts[[f]], "smooth_bin_curve"),
                                                         wt_bin_diff = safe_extract(ret_sel_pen_wts[[f]], "smooth_bin_diff"),
                                                         wt_yr_diff = safe_extract(ret_sel_pen_wts[[f]], "smooth_yr_diff"),
                                                         wt_yr_curve = safe_extract(ret_sel_pen_wts[[f]], "smooth_yr_curve"),
                                                         wt_dome = safe_extract(ret_sel_pen_wts[[f]], "smooth_dome"),
                                                         wt_mean_center = safe_extract(ret_sel_pen_wts[[f]], "smooth_mean_center"),
                                                         normalize = ret_sel_pen_wts[[f]]$normalize,
                                                         bin_range = ret_sel_pen_wts[[f]]$bin_range,
                                                         yr_diff_ref = ret_sel_pen_wts[[f]]$yr_diff_ref)
      } # end if
    } # end f loop

    for(sf in 1:n_srv_fleets) {

      # if bicubic
      bicubic_yrs = which(srv_sel_model[r,,sf] == 8)
      has_bicubic = length(bicubic_yrs) > 0
      has_nonzero_pen = any(unlist(lapply(smooth_pen_terms, function(data_name) safe_extract(srv_sel_pen_wts[[sf]], data_name))) != 0)
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
                                                         wt_bin_curve = safe_extract(srv_sel_pen_wts[[sf]], "smooth_bin_curve"),
                                                         wt_bin_diff = safe_extract(srv_sel_pen_wts[[sf]], "smooth_bin_diff"),
                                                         wt_yr_diff = safe_extract(srv_sel_pen_wts[[sf]], "smooth_yr_diff"),
                                                         wt_yr_curve = safe_extract(srv_sel_pen_wts[[sf]], "smooth_yr_curve"),
                                                         wt_dome = safe_extract(srv_sel_pen_wts[[sf]], "smooth_dome"),
                                                         wt_mean_center = safe_extract(srv_sel_pen_wts[[sf]], "smooth_mean_center"),
                                                         normalize = srv_sel_pen_wts[[sf]]$normalize,
                                                         bin_range = srv_sel_pen_wts[[sf]]$bin_range,
                                                         yr_diff_ref = srv_sel_pen_wts[[sf]]$yr_diff_ref)
      } # end if
    } # end sf loop

  } # end r loop


  ### Selectivity (Prior) ----------------------------------------------------
  # Total Fishery selectivity
  if(Use_fish_selex_prior == 1) sel_nLL = sel_nLL + get_selex_prior(fish_selex_prior, fish_fixed_sel_pars, fish_sel, fish_sel_l, fish_selex_type, fish_sel_blocks)

  # Retained Fishery selectivity
  if(Use_ret_selex_prior == 1) sel_nLL = sel_nLL + get_selex_prior(ret_selex_prior, ret_fixed_sel_pars, ret_sel, ret_sel_l, ret_selex_type, ret_sel_blocks)

  # Survey selectivity
  if(Use_srv_selex_prior == 1) sel_nLL = sel_nLL + get_selex_prior(srv_selex_prior, srv_fixed_sel_pars, srv_sel, srv_sel_l, srv_selex_type, srv_sel_blocks)


  ### Selectivity Parameter Centering (Penalty) ------------------------------
  if(Use_fish_selex_penalty == 1) sel_nLL = sel_nLL + get_selex_fixed_penalty(fish_selex_penalty, fish_fixed_sel_pars)
  if(Use_ret_selex_penalty == 1) sel_nLL = sel_nLL + get_selex_fixed_penalty(ret_selex_penalty, ret_fixed_sel_pars)
  if(Use_srv_selex_penalty == 1) sel_nLL = sel_nLL + get_selex_fixed_penalty(srv_selex_penalty, srv_fixed_sel_pars)

  ## Recruitment and Growth --------------------------------------------------
  ### Recruitment (Penalty) --------------------------------------------------
  tmp_rec_pen = get_recruitment_penalty(
    n_pop = n_pop,
    n_regions = n_regions,
    n_ages = n_ages,
    n_est_rec_devs = n_est_rec_devs,
    rec_region_prop_spec = rec_region_prop_spec,
    rec_region_prop = rec_region_prop,
    equil_init_age_strc = equil_init_age_strc,
    ln_InitDevs = ln_InitDevs,
    init_age_devs_shared = init_age_devs_shared,
    ln_sigmaR = ln_sigmaR,
    bias_ramp = bias_ramp,
    sigmaR_switch = sigmaR_switch,
    ln_RecDevs = ln_RecDevs,
    sigmaR2_early = sigmaR2_early,
    sigmaR2_late = sigmaR2_late,
    do_rec_bias_ramp = do_rec_bias_ramp,
    map_ln_RecDevs = map_ln_RecDevs,
    RecDevs_model = RecDevs_model,
    RecDevs_rho = RecDevs_rho,
    RecDevs_rw_init_sigma = RecDevs_rw_init_sigma,
    RecDevs_pen_center = RecDevs_pen_center,
    InitDevs_pen_center = InitDevs_pen_center,
    init_devs_pen_use = init_devs_pen_use,
    init_bias_ramp = init_bias_ramp,
    map_ln_InitDevs = map_ln_InitDevs,
    Use_init_sex_pen = Use_init_sex_pen,
    ln_sigma_init_sex = ln_sigma_init_sex
  )
  Init_Rec_nLL = tmp_rec_pen$Init_Rec_nLL
  Init_Sex_nLL = tmp_rec_pen$Init_Sex_nLL
  Rec_nLL = tmp_rec_pen$Rec_nLL

  ### Recruitment Level (Penalty) --------------------------------------------
  # A second, separate statement about recruitment: the stock-recruit penalty
  # constrains the residuals about the curve, this one constrains the series.
  Rec_level_nLL = array(0, dim = dim(Rec))
  if(Use_rec_level_pen == 1) {
    Rec_level_nLL = get_rec_level_penalty(Rec, exp(ln_sigma_rec_level), rec_level_pen_center, if(all(rec_level_pen_yrs == 1)) NULL else which(rec_level_pen_yrs == 1))
  }

  ### Stock-Recruit Residual (Penalty) ---------------------------------------
  # only reachable under mean recruitment; the curve penalizes the recruitment
  # series without generating it, so the deviations stay free
  SR_pen_nLL = array(0, dim = dim(Rec))
  if(sr_penalty > 0) {
    SR_pen_nLL = get_sr_penalty(Rec, SR_pred, exp(ln_sigma_sr_pen),
                                if(all(sr_pen_yrs == 1)) NULL else which(sr_pen_yrs == 1))
  }

  ### Growth (Process Error) -------------------------------------------------
  # Growth deviations either on parameters or semi-parametric
  if(growth_model != 0) {
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {

        # growth deviations on parameters
        if(!is.null(growth_tv_model)) for(k in which(growth_tv_model > 0)) {
          growth_tv_nLL = growth_tv_nLL - Get_PE_loglik(
            PE_model = growth_tv_model[k],
            PE_pars = array(growth_pe_pars[p,r,k,,1], dim = c(1, 1, n_sexes, 1)),
            ln_devs = array(ln_growth_devs[p,r,,k,], dim = c(1, n_yrs, 1, n_sexes, 1)),
            map_sel_devs = array(map_ln_growth_devs[p,r,,k,], dim = c(1, n_yrs, 1, n_sexes)),
            map_sel_devs_full = array(map_ln_growth_devs[,,,k,], dim = c(n_pop * n_regions, n_yrs, 1, n_sexes)),
            min_sel_devs_shared_bins = 1,
            rw_init_sigma = growth_rw_init_sigma
          )
        } # end k loop

        # growth deviations on parametric form
        if(growth_semipar > 0) {
          growth_semipar_nLL = growth_semipar_nLL - Get_PE_loglik(
            PE_model = growth_semipar,
            PE_pars = array(growth_pe_pars[p,r,,,2], dim = c(1, dim(growth_pe_pars)[3], n_sexes, 1)),
            ln_devs = array(ln_growth_semipar_devs[p,r,,,], dim = c(1, n_yrs, n_ages, n_sexes, 1)),
            map_sel_devs = array(map_ln_growth_semipar_devs[p,r,,,], dim = c(1, n_yrs, n_ages, n_sexes)),
            map_sel_devs_full = array(map_ln_growth_semipar_devs, dim = c(n_pop * n_regions, n_yrs, n_ages, n_sexes)),
            min_sel_devs_shared_bins = if(is.null(growth_semipar_bins)) 1:n_ages else growth_semipar_bins,
            rw_init_sigma = growth_rw_init_sigma
          )
        }

      } # end r loop
    } # end p loop
  } # end growth process error

  ### Initial Recruitment Offset (Penalty) -----------------------------------
  # The initial equilibrium recruitment against the recruitment level
  if(Use_rinit_pen == 1 && use_rinit == 1) {
    for(p in 1:n_pop) rinit_nLL = rinit_nLL - RTMB::dnorm(ln_rinit[p] - ln_R0_mat[p, R0_ref_block], 0, rinit_pen_sd, TRUE)
  }

  ## Catchability and Biology Priors -----------------------------------------
  ### Fishery Catchability (Prior) -------------------------------------------
  if(Use_fish_q_prior == 1) fish_q_nLL = fish_q_nLL + get_q_prior(fish_q_prior, ln_fish_q)

  ### Survey Catchability (Prior) --------------------------------------------
  if(Use_srv_q_prior == 1) srv_q_nLL = srv_q_nLL + get_q_prior(srv_q_prior, ln_srv_q)

  ### Natural Mortality (Prior) ----------------------------------------------
  if(Use_M_prior == 1) M_nLL = M_nLL + get_natmort_prior(M_prior, ln_M, M_blocks)

  ### Steepness (Prior) ------------------------------------------------------
  if(Use_h_prior == 1) h_nLL = h_nLL + get_steepness_prior(h_prior, h_trans)

  ## Movement ----------------------------------------------------------------
  ### Movement Rates (Penalty) -----------------------------------------------
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

  ### Movement Rates (Prior) -------------------------------------------------
  # For CTMC movement the prior is placed on the annual fractions exp(Q), so that alpha
  # means the same thing regardless of how the year is divided into seasons.
  if(Use_Movement_Prior == 1) Movement_nLL = Movement_nLL + get_movement_dirichlet_prior(Movement_prior, Movement, Mrate)

  ## Recruitment Apportionment and Tagging -----------------------------------
  ### Recruitment R0 and Proportions (Prior) ---------------------------------
  # Regional/seasonal apportionment and stray rate priors all feed rec_prop_nLL
  rec_prop_nLL = get_recruitment_proportion_priors(
    use_rec_region_prop_prior = use_rec_region_prop_prior,
    rec_region_prop_prior = rec_region_prop_prior,
    rec_region_prop = rec_region_prop,
    use_rec_seas_prop_prior = use_rec_seas_prop_prior,
    use_fixed_rec_seas_prop = use_fixed_rec_seas_prop,
    rec_seas_prop_prior = rec_seas_prop_prior,
    rec_seas_prop = rec_seas_prop,
    rec_lag = rec_lag,
    spawn_seas = spawn_seas,
    n_seas = n_seas,
    use_stray_rate_prior = use_stray_rate_prior,
    stray_rate_prior = stray_rate_prior,
    stray_rate_pars = stray_rate_pars
  )

  if(use_r0_prior == 1) R0_nLL = get_r0_prior(r0_prior, ln_R0_mat[, R0_ref_block]) # recruitment R0 (global scalar prior, unweighted by Wt_Rec)

  ### Tag Reporting Rate (Prior) ---------------------------------------------
  if(use_conv_tag_fishrep_prior == 1) TagRep_nLL = TagRep_nLL + get_tagrep_prior(conv_tag_fishrep_prior, conv_tag_fish_reporting_pars)

  # Joint Negative Log Likelihood -------------------------------------------
  jnLL = sum(Wt_Catch * Catch_nLL) +             # Aggregated catch likelihoods
    sum(Wt_Catch * apply(CatchAA_nLL, c(1,2,3,6), sum)) + # Catch at age likelihoods
    sum(Wt_Catch_pop * apply(CatchAA_pop_nLL, c(1,2,3,4,7), sum)) + # Population-specific catch at age
    sum(Wt_Discard * apply(DiscardAA_nLL, c(1,2,3,6), sum)) +   # Discard at age likelihoods
    sum(Wt_Discard_pop * apply(DiscardAA_pop_nLL, c(1,2,3,4,7), sum)) + # Population-specific discard at age
    sum(Wt_SrvIdx * apply(SrvIdxAA_nLL, c(1,2,3,6), sum)) +     # Survey index at age likelihoods
    sum(Wt_SrvIdx_pop * apply(SrvIdxAA_pop_nLL, c(1,2,3,4,7), sum)) + # Population-specific survey index at age
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
    sum(Fish_caal_nLL) +                      # Fishery conditional age-at-length likelihoods
    sum(Srv_caal_nLL) +                       # Survey conditional age-at-length likelihoods
    (Wt_Tagging * sum(conv_fish_tag_nLL)) +   # Tagging likelihood
    (Wt_F * sum(Fmort_nLL)) +                 # Fishing mortality penalty
    (Wt_D * sum(dmr_nLL)) +                   # Discard mortality rate penalty
    sum(Wt_Rec * Rec_nLL) +                   # Recruitment penalty
    sum(Wt_Init_Rec * Init_Rec_nLL) +         # Initial age penalty
    sum(Init_Sex_nLL) +                       # Initial age deviations, tie between sexes
    sum(Rec_level_nLL) +                      # Recruitment level penalty
    sum(SR_pen_nLL) +                         # Stock-recruit residual penalty (mean recruitment only)
    NAA_state_nLL +                            # State-space numbers at age penalty
    sel_nLL +                                  # Selectivity penalty
    M_nLL +                                    # Natural mortality prior
    R0_nLL +                                   # Global R0 prior
    h_nLL +                                    # Steepness prior
    Movement_nLL +                             # Movement prior
    TagRep_nLL +                               # Tag reporting rate prior
    fish_q_nLL +                               # Fishery q prior
    srv_q_nLL +                                # Survey q prior
    rec_prop_nLL +                             # Recruitment proportion prior
    growth_tv_nLL +                            # Time-varying growth process error
    growth_semipar_nLL +                       # Semi-parametric growth process error
    rinit_nLL                                  # Initial recruitment offset penalty

  # Report Section ----------------------------------------------------------
  ## Biological Processes ----------------------------------------------------
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
  RTMB::REPORT(NAA_int)
  RTMB::REPORT(NAA_pred)
  RTMB::REPORT(NAA_scalar)
  RTMB::REPORT(NAA_state_nLL)
  RTMB::REPORT(ZAA)
  RTMB::REPORT(natmort)
  RTMB::REPORT(bias_ramp)
  RTMB::REPORT(init_bias_ramp)
  RTMB::REPORT(Movement)
  RTMB::REPORT(sgl_seas_spawning_movement)
  RTMB::REPORT(Mrate)

  ## Fishery Processes -------------------------------------------------------
  RTMB::REPORT(init_F)
  RTMB::REPORT(ln_sigmaC)
  RTMB::REPORT(ln_sigmaCAA)
  RTMB::REPORT(trans_rho_catch)
  RTMB::REPORT(trans_rho_discard)
  RTMB::REPORT(trans_rho_srv_idx)
  RTMB::REPORT(CatchAA_nLL)
  RTMB::REPORT(PredCatchAA)
  RTMB::REPORT(PredDiscardAA)
  RTMB::REPORT(PredSrvIdxAA)
  RTMB::REPORT(PredCatchAA_pop)
  RTMB::REPORT(PredDiscardAA_pop)
  RTMB::REPORT(PredSrvIdxAA_pop)
  RTMB::REPORT(DiscardAA_nLL)
  RTMB::REPORT(DiscardAA_pop_nLL)
  RTMB::REPORT(CatchAA_pop_nLL)
  RTMB::REPORT(SrvIdxAA_nLL)
  RTMB::REPORT(SrvIdxAA_pop_nLL)
  RTMB::REPORT(ln_sigmaSrvIdxAA)
  RTMB::REPORT(ln_sigmaFishIdx)
  RTMB::REPORT(ln_sigmaSrvIdx)
  RTMB::REPORT(ln_sigmaFishIdx_pop)
  RTMB::REPORT(ln_sigmaSrvIdx_pop)
  RTMB::REPORT(FishIdx_SD)
  RTMB::REPORT(SrvIdx_SD)
  RTMB::REPORT(FishIdx_pop_SD)
  RTMB::REPORT(SrvIdx_pop_SD)
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
  if(do_caal == 1) {
    RTMB::REPORT(Fish_caal)
    RTMB::REPORT(Fish_caal_discard)
  }
  RTMB::REPORT(PredCatch)
  RTMB::REPORT(PredDiscard)
  RTMB::REPORT(PredFishIdx)
  RTMB::REPORT(fish_sel)
  RTMB::REPORT(ret_sel)
  RTMB::REPORT(fish_q)

  ## Survey Processes --------------------------------------------------------
  RTMB::REPORT(PredSrvIdx)
  RTMB::REPORT(srv_sel)
  RTMB::REPORT(srv_q)
  RTMB::REPORT(SrvIAA)
  RTMB::REPORT(SrvIAL)
  if(do_caal == 1) RTMB::REPORT(Srv_caal)
  RTMB::REPORT(RecDev_anom)


  # per-fleet keys, growth-derived or fixed data, reported whenever either source supplied one,
  # so a fixed SizeAgeTrans_fish shows up the way a supplied WAA_fish already does
  if(!is.null(SizeAgeTrans_fish)) RTMB::REPORT(SizeAgeTrans_fish)
  if(!is.null(SizeAgeTrans_srv)) RTMB::REPORT(SizeAgeTrans_srv)

  # Growth module
  if(growth_model != 0) {
    RTMB::REPORT(mean_LAA_fish)
    RTMB::REPORT(sd_LAA_fish)
    RTMB::REPORT(mean_LAA_srv)
    RTMB::REPORT(sd_LAA_srv)
    RTMB::REPORT(mean_LAA_spawn)
    RTMB::REPORT(sd_LAA_spawn)
    RTMB::REPORT(Linf)
    RTMB::REPORT(SizeAgeTrans_spawn)
    RTMB::REPORT(L_beg)
    RTMB::REPORT(growth_pars_y)
    if(derive_waa == 1) {
      RTMB::REPORT(WAA)
      RTMB::REPORT(WAA_fish)
      RTMB::REPORT(WAA_srv)
    }
  }

  # Report length-based selectivity
  if(fish_selex_type == 1) RTMB::REPORT(fish_sel_l)
  if(ret_selex_type == 1) RTMB::REPORT(ret_sel_l)
  if(srv_selex_type == 1) RTMB::REPORT(srv_sel_l)

  ## Tagging -----------------------------------------------------------------
  if(any(use_conv_fish_tagging == 1)) {
    RTMB::REPORT(pred_conv_tag_fish_recap)
    RTMB::REPORT(conv_tag_fish_avail)
    RTMB::REPORT(conv_tag_fish_reporting)
  }

  ## Parameter Deviations ----------------------------------------------------
  RTMB::REPORT(ln_RecDevs)
  RTMB::REPORT(move_devs)
  RTMB::REPORT(ln_fishsel_devs)
  RTMB::REPORT(ln_srvsel_devs)

  ## Aggregated Likelihoods --------------------------------------------------
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

  ## Population-specific Likelihoods -----------------------------------------
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
  if(do_fish_caal) RTMB::REPORT(Fish_caal_nLL)
  if(do_srv_caal) RTMB::REPORT(Srv_caal_nLL)

  ## Penalties and Priors ----------------------------------------------------
  RTMB::REPORT(M_nLL)
  RTMB::REPORT(Fmort_nLL)
  RTMB::REPORT(dmr_nLL)
  RTMB::REPORT(Rec_nLL)
  RTMB::REPORT(Init_Rec_nLL)
  RTMB::REPORT(Init_Sex_nLL)
  RTMB::REPORT(Rec_level_nLL)
  RTMB::REPORT(SR_pen_nLL)
  RTMB::REPORT(SR_pred)
  RTMB::REPORT(sr_R0)
  RTMB::REPORT(conv_fish_tag_nLL)
  RTMB::REPORT(h_nLL)
  RTMB::REPORT(R0_nLL)
  RTMB::REPORT(fish_q_nLL)
  RTMB::REPORT(sel_nLL)
  RTMB::REPORT(R0_yr)
  RTMB::REPORT(growth_tv_nLL)
  RTMB::REPORT(growth_semipar_nLL)
  RTMB::REPORT(rinit_nLL)
  RTMB::REPORT(srv_q_nLL)
  RTMB::REPORT(Movement_nLL)
  RTMB::REPORT(TagRep_nLL)
  RTMB::REPORT(rec_prop_nLL)
  RTMB::REPORT(jnLL)

  ## Derived Quantities ------------------------------------------------------
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
