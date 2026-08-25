# Stage 2 of 3: objective function
#
# The one RTMB objective function, evaluated by fit_model. There is no per
# configuration branching at the file level: SPoRC_rtmb always runs end to end
# and is regulated by the switches set during setup (rec_model, move_type, the
# likelihood and OSA flags). It works through parameter transforms, mortality,
# recruitment, initial age structure, population projection, the observation
# models, the likelihood components, and finally priors and penalties.
#
# Almost none of the array arithmetic happens inline. Each section delegates to a
# module in one of the model_*.R files, so this file reads as the order of
# operations rather than as the arithmetic itself.

# Development history of the objective function, kept as provenance for the
# assessment. Current behaviour is documented at each section below, not here.

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

#' Fill in defaults for input lists built by older versions of SPoRC
#'
#' Assigns any data or parameter objects that a current \code{SPoRC_rtmb} call
#' expects but that older input lists predate, so previously built objects keep
#' evaluating unchanged. Values are written into \code{env} only when absent and
#' are never overwritten.
#'
#' @param env Environment holding the unpacked data and parameters, i.e. the
#'   \code{SPoRC_rtmb} frame after \code{RTMB::getAll}. Defaults to the caller.
#'
#' @return \code{NULL}, invisibly. Called for its side effects on \code{env}.
#'
#' @keywords internal
maintain_backwards_compatibility <- function(env = parent.frame()) {

  has <- function(x) exists(x, envir = env, inherits = FALSE)
  set <- function(x, value) assign(x, value, envir = env)

  # Movement timing options.
  if(!has("move_timing")) set("move_timing", 0)
  if(!has("comp_const_obs")) set("comp_const_obs", 1)

  # CAAL
  if(!has("do_caal")) set("do_caal", 0)

  # Growth stuff
  if(!has("growth_model")) set("growth_model", 0)
  if(!has("derive_waa")) set("derive_waa", 0)
  if(!has("wt_len_pars")) set("wt_len_pars", NULL)
  if(!has("growth_len_mid_vals")) set("growth_len_mid_vals", NULL)
  if(!has("growth_tv_model")) set("growth_tv_model", NULL)
  if(!has("growth_tv_link")) set("growth_tv_link", 0)
  if(!has("growth_par_bounds")) set("growth_par_bounds", NULL)
  if(!has("growth_tv_type")) set("growth_tv_type", 0)
  if(!has("growth_cohort_styr")) set("growth_cohort_styr", 1)
  if(!has("growth_rw_init_sigma")) set("growth_rw_init_sigma", 5)
  if(!has("growth_semipar")) set("growth_semipar", 0)
  if(!has("growth_L2_asymptote")) set("growth_L2_asymptote", 0)
  if(!has("ln_growth_semipar_devs")) set("ln_growth_semipar_devs", NULL)
  if(!has("growth_pe_pars")) set("growth_pe_pars", NULL)
  if(!has("map_ln_growth_semipar_devs")) set("map_ln_growth_semipar_devs", NULL)
  if(!has("growth_semipar_bins")) set("growth_semipar_bins", NULL)
  if(!has("ln_growth_devs")) set("ln_growth_devs", NULL)
  if(!has("map_ln_growth_devs")) set("map_ln_growth_devs", NULL)
  if(!has("fish_len_comp_sel")) set("fish_len_comp_sel", rep(0, get("n_fish_fleets", envir = env)))
  if(!has("srv_len_comp_sel")) set("srv_len_comp_sel", rep(0, get("n_srv_fleets", envir = env)))
  if(!has("fish_waa_selected")) set("fish_waa_selected", rep(0, get("n_fish_fleets", envir = env)))
  if(!has("srv_waa_selected")) set("srv_waa_selected", rep(0, get("n_srv_fleets", envir = env)))
  if(!has("SizeAgeTrans_fish")) set("SizeAgeTrans_fish", NULL)
  if(!has("SizeAgeTrans_srv")) set("SizeAgeTrans_srv", NULL)
  if(!has("LenBinMap")) set("LenBinMap", NULL)
  if(!has("Use_rinit_pen")) set("Use_rinit_pen", 0)
  if(!has("rinit_pen_sd")) set("rinit_pen_sd", 1)
  if(!has("UseFish_caal")) set("UseFish_caal", NULL)
  if(!has("UseSrv_caal")) set("UseSrv_caal", NULL)

  # The conditional age-at-length weights live in Setup_Mod_Weighting; a list
  # without them weights every row at one.
  caal_wt_dim <- function(n_fleets) c(get("n_regions", envir = env), length(get("years", envir = env)), get("n_seas", envir = env),
                                      length(get("lens", envir = env)), get("n_sexes", envir = env), n_fleets)
  if(!has("Wt_Fish_caal")) set("Wt_Fish_caal", array(1, dim = caal_wt_dim(get("n_fish_fleets", envir = env))))
  if(!has("Wt_Srv_caal")) set("Wt_Srv_caal", array(1, dim = caal_wt_dim(get("n_srv_fleets", envir = env))))

  # Reference year for the biological inputs to unfished spawning biomass per
  # recruit. The first year is what the model always used, so that stays the default.
  if(!has("SR_ref_yr")) set("SR_ref_yr", 1)
  if(!has("ctmc_scale_by_seasdur")) set("ctmc_scale_by_seasdur", 0)

  # Initialization of F. Older input lists carried init_F_prop, a fixed proportion of
  # the mean F held as data; it is now the init_F_par parameter plus init_F_form,
  # so a stored proportion maps onto the logit scale of the "prop" form.
  if(!has("init_F_form")) set("init_F_form", 0)
  if(!has("init_F_par")) {
    init_F_dim <- c(get("n_regions", envir = env), get("n_seas", envir = env), get("n_fish_fleets", envir = env))
    init_F_prop <- if(has("init_F_prop")) get("init_F_prop", envir = env) else array(0, dim = init_F_dim)
    set("init_F_par", array(stats::qlogis(pmin(pmax(init_F_prop, 1e-10), 1 - 1e-10)), dim = init_F_dim))
  }

  # Deviation maps mirrored into the data lists. Without the initial age map the
  # penalty falls on every cell, which is what lists built before it carried did.
  if(!has("map_ln_InitDevs")) set("map_ln_InitDevs", NULL)
  if(!has("map_ln_F_devs") || !has("map_logit_dmr_devs")) {
    UseCatch <- get("UseCatch", envir = env)
    has_catch <- UseCatch == 1 |
      apply(get("UseCatch_pop", envir = env) == 1, c(2,3,4,5), any) |
      is.na(get("ObsCatch", envir = env))
    legacy_map <- array(NA_real_, dim = dim(UseCatch))
    legacy_map[has_catch] <- seq_len(sum(has_catch))
    if(!has("map_ln_F_devs")) set("map_ln_F_devs", legacy_map)
    if(!has("map_logit_dmr_devs")) set("map_logit_dmr_devs", legacy_map)
  }

  # Data lists built before the recruitment map mirror existed penalize every
  # deviation, so default to that rather than silently changing their objective
  if(!has("map_ln_RecDevs")) {
    set("map_ln_RecDevs", array(1, dim = dim(get("ln_RecDevs", envir = env))))
  }

  # Fishery index timing. Older input lists have no t_fish, and the model used to
  # form the index from start-of-season numbers, so default to zero rather than
  # shifting their predicted index.
  if(!has("t_fish")) {
    set("t_fish", array(0, dim = c(get("n_regions", envir = env),
                                   get("n_seas", envir = env),
                                   get("n_fish_fleets", envir = env))))
  }

  # Index age selection, catchability solving, and index error structure. Older
  # input lists sum the index over all ages, estimate every catchability, and
  # use a lognormal index likelihood, so default to that.
  n_ages_bc <- length(get("ages", envir = env))
  n_srv_bc <- get("n_srv_fleets", envir = env)
  n_fish_bc <- get("n_fish_fleets", envir = env)
  if(!has("srv_idx_ages")) set("srv_idx_ages", array(1, dim = c(n_ages_bc, n_srv_bc)))
  if(!has("fish_idx_ages")) set("fish_idx_ages", array(1, dim = c(n_ages_bc, n_fish_bc)))
  if(!has("srv_q_type")) set("srv_q_type", rep(0, n_srv_bc))

  # Fishery catchability gained the same analytic and covariate machinery the
  # survey already had. Absent on older input lists, where every fishery q was
  # estimated with no covariates, which is what these defaults reproduce.
  n_fish_bc <- get("n_fish_fleets", envir = env)
  if(!has("fish_q_type")) set("fish_q_type", rep(0, n_fish_bc))
  if(!has("do_fish_q_cov")) set("do_fish_q_cov", 0)
  if(!has("fish_q_cov")) set("fish_q_cov", array(0, dim = c(get("n_regions", envir = env), length(get("years", envir = env)), n_fish_bc, 1)))
  if(!has("fish_q_coeff")) set("fish_q_coeff", array(0, dim = c(get("n_regions", envir = env), n_fish_bc, 1)))
  if(!has("SrvIdx_LikeType")) set("SrvIdx_LikeType", rep(0, n_srv_bc))
  if(!has("FishIdx_LikeType")) set("FishIdx_LikeType", rep(0, n_fish_bc))
  if(!has("SrvIdx_Cov")) set("SrvIdx_Cov", vector("list", n_srv_bc))
  if(!has("FishIdx_Cov")) set("FishIdx_Cov", vector("list", n_fish_bc))

  # Deviation penalties centred on a fixed prior mean unless asked otherwise.
  if(!has("Fdev_pen_center")) set("Fdev_pen_center", 0)
  if(!has("RecDevs_pen_center")) set("RecDevs_pen_center", 0)
  if(!has("InitDevs_pen_center")) set("InitDevs_pen_center", 0)

  # Only read by the initial-age penalty's shared-subset case
  # (equil_init_age_strc == 3), which setup now stores in data; a list without
  # it cannot be using that case, so NULL is never evaluated.
  if(!has("init_age_devs_shared")) set("init_age_devs_shared", NULL)

  # The initial age penalty used to share Wt_Rec, which only worked because both
  # were scalars applied outside the sum.
  if(!has("Wt_Init_Rec")) set("Wt_Init_Rec", get("Wt_Rec", envir = env))
  # An array weight from before the sex dimension is one weight per age; it repeats across sexes so it conforms with the sex-dimensioned penalty array
  wt_init_bc <- get("Wt_Init_Rec", envir = env)
  if(length(wt_init_bc) > 1 && length(dim(wt_init_bc)) == 3) set("Wt_Init_Rec", array(rep(wt_init_bc, get("n_sexes", envir = env)), dim = c(dim(wt_init_bc), get("n_sexes", envir = env))))

  # The recruitment level penalty is off unless asked for.
  if(!has("Use_rec_level_pen")) set("Use_rec_level_pen", 0)
  if(!has("ln_sigma_rec_level")) set("ln_sigma_rec_level", 0)

  # The between-sex likelihood on initial age deviations is off unless specified
  if(!has("Use_init_sex_pen")) set("Use_init_sex_pen", 0)
  if(!has("ln_sigma_init_sex")) set("ln_sigma_init_sex", 0)
  if(!has("rec_level_pen_center")) set("rec_level_pen_center", 1)
  if(!has("rec_level_pen_yrs")) set("rec_level_pen_yrs", rep(1, length(get("years", envir = env))))

  # The stock-recruit penalty under mean recruitment is off unless asked for.
  if(!has("sr_penalty")) set("sr_penalty", 0)
  if(!has("sr_R0_spec")) set("sr_R0_spec", 0)
  if(!has("ln_sigma_sr_pen")) set("ln_sigma_sr_pen", 0)
  if(!has("sr_pen_yrs")) set("sr_pen_yrs", rep(1, length(get("years", envir = env))))

  # Selectivity process error weights. Older input lists always applied it.
  if(!has("fishsel_pe_wt")) set("fishsel_pe_wt", rep(1, n_fish_bc))
  if(!has("retsel_pe_wt")) set("retsel_pe_wt", rep(1, n_fish_bc))
  if(!has("srvsel_pe_wt")) set("srvsel_pe_wt", rep(1, n_srv_bc))
  if(!has("fishsel_rw_init_sigma")) set("fishsel_rw_init_sigma", rep(5, n_fish_bc))
  if(!has("retsel_rw_init_sigma")) set("retsel_rw_init_sigma", rep(5, n_fish_bc))
  if(!has("srvsel_rw_init_sigma")) set("srvsel_rw_init_sigma", rep(5, n_srv_bc))
  if(!has("fishsel_bin_devs_rw_init_sigma")) set("fishsel_bin_devs_rw_init_sigma", rep(5, n_fish_bc))
  if(!has("retsel_bin_devs_rw_init_sigma")) set("retsel_bin_devs_rw_init_sigma", rep(5, n_fish_bc))
  if(!has("srvsel_bin_devs_rw_init_sigma")) set("srvsel_bin_devs_rw_init_sigma", rep(5, n_srv_bc))

  # Composition bin ranges. Older input lists fit every bin.
  if(!has("FishAgeComps_bins")) set("FishAgeComps_bins", array(1, dim = c(n_ages_bc, n_fish_bc)))
  if(!has("SrvAgeComps_bins")) set("SrvAgeComps_bins", array(1, dim = c(n_ages_bc, n_srv_bc)))

  # Bin-override selectivity deviations. Older input lists have none
  for(pre in c("fish", "ret", "srv")) {
    n_fl_bc <- if(pre == "srv") n_srv_bc else n_fish_bc
    if(!has(paste0(pre, "_sel_bin_dev_bins"))) set(paste0(pre, "_sel_bin_dev_bins"), array(0, dim = c(n_ages_bc, n_fl_bc)))
    if(!has(paste0(pre, "_sel_norm_bins"))) set(paste0(pre, "_sel_norm_bins"), array(1, dim = c(n_ages_bc, n_fl_bc)))
    if(!has(paste0("cont_tv_", pre, "sel_bin_devs"))) set(paste0("cont_tv_", pre, "sel_bin_devs"), rep(0, n_fl_bc))
    if(!has(paste0("ln_", pre, "sel_bin_devs"))) set(paste0("ln_", pre, "sel_bin_devs"), array(0, dim = c(get("n_regions", envir = env), length(get("years", envir = env)) + get("n_proj_yrs_devs", envir = env), n_ages_bc, get("n_sexes", envir = env), n_fl_bc)))
    if(!has(paste0(pre, "sel_bin_devs_pe_pars"))) set(paste0(pre, "sel_bin_devs_pe_pars"), array(0, dim = c(get("n_regions", envir = env), n_ages_bc, get("n_sexes", envir = env), n_fl_bc)))
    if(!has(paste0("map_ln_", pre, "sel_bin_devs"))) set(paste0("map_ln_", pre, "sel_bin_devs"), array(NA_real_, dim = dim(get(paste0("ln_", pre, "sel_bin_devs"), envir = env))))
  } # end pre loop

  # Selectivity parameter centering penalties. The flag guards every reference to
  # the table, so older input lists need only the flag.
  if(!has("Use_fish_selex_penalty")) set("Use_fish_selex_penalty", 0)
  if(!has("Use_ret_selex_penalty")) set("Use_ret_selex_penalty", 0)
  if(!has("Use_srv_selex_penalty")) set("Use_srv_selex_penalty", 0)

  # Bicubic bookkeeping arrays, which the parametric plateau now also reads.
  # Input lists predating bicubic selectivity have none of them; all-zero means
  # no plateau and no bicubic block anywhere, which is what those lists were.
  for(pre_arr in c("fish", "ret")) for(suf in c("binnodes", "yrnodes", "selstyr", "nselbins")) {
    nm_arr <- paste0(pre_arr, "_sel_bicubic_", suf)
    if(!has(nm_arr)) set(nm_arr, array(0, dim = c(get("n_regions", envir = env), length(get("years", envir = env)), n_fish_bc)))
  }
  for(suf in c("binnodes", "yrnodes", "selstyr", "nselbins")) {
    nm_arr <- paste0("srv_sel_bicubic_", suf)
    if(!has(nm_arr)) set(nm_arr, array(0, dim = c(get("n_regions", envir = env), length(get("years", envir = env)), n_srv_bc)))
  }

  # Sex offsets on selectivity. Older input lists have neither the flags nor the
  # scale parameters, and both defaults reproduce sex-independent selectivity.
  # Double normal limbs are anchored at the end bins unless a model asked for the
  # raw form, which no older input list did.
  if(!has("fish_dbnrml_raw")) set("fish_dbnrml_raw", array(0, dim = c(n_fish_bc, 2)))
  if(!has("ret_dbnrml_raw")) set("ret_dbnrml_raw", array(0, dim = c(n_fish_bc, 2)))
  if(!has("srv_dbnrml_raw")) set("srv_dbnrml_raw", array(0, dim = c(n_srv_bc, 2)))
  if(!has("fish_dbnrml_startbin")) set("fish_dbnrml_startbin", rep(1, n_fish_bc))
  if(!has("ret_dbnrml_startbin")) set("ret_dbnrml_startbin", rep(1, n_fish_bc))
  if(!has("srv_dbnrml_startbin")) set("srv_dbnrml_startbin", rep(1, n_srv_bc))
  if(!has("fishsel_sex_par_offset")) set("fishsel_sex_par_offset", rep(0, n_fish_bc))
  if(!has("srvsel_sex_par_offset")) set("srvsel_sex_par_offset", rep(0, n_srv_bc))
  if(!has("fishsel_sex_scale_offset")) set("fishsel_sex_scale_offset", rep(0, n_fish_bc))
  if(!has("srvsel_sex_scale_offset")) set("srvsel_sex_scale_offset", rep(0, n_srv_bc))
  if(!has("fishsel_sex_apical_offset")) set("fishsel_sex_apical_offset", rep(0, n_fish_bc))
  if(!has("srvsel_sex_apical_offset")) set("srvsel_sex_apical_offset", rep(0, n_srv_bc))
  if(!has("ln_fishsel_sex_scale")) set("ln_fishsel_sex_scale", array(0, dim = c(get("n_regions", envir = env), dim(get("fish_fixed_sel_pars", envir = env))[3], get("n_sexes", envir = env), n_fish_bc)))
  if(!has("ln_srvsel_sex_scale")) set("ln_srvsel_sex_scale", array(0, dim = c(get("n_regions", envir = env), dim(get("srv_fixed_sel_pars", envir = env))[3], get("n_sexes", envir = env), n_srv_bc)))
  if(!has("retsel_sex_par_offset")) set("retsel_sex_par_offset", rep(0, n_fish_bc))
  if(!has("retsel_sex_scale_offset")) set("retsel_sex_scale_offset", rep(0, n_fish_bc))
  if(!has("retsel_sex_apical_offset")) set("retsel_sex_apical_offset", rep(0, n_fish_bc))
  if(!has("ln_retsel_sex_scale")) set("ln_retsel_sex_scale", array(0, dim = c(get("n_regions", envir = env), dim(get("ret_fixed_sel_pars", envir = env))[3], get("n_sexes", envir = env), n_fish_bc)))

  # Initial age deviations gained a sex dimension. An older 3-D array is one
  # shared curve, so it broadcasts across sexes and only the first sex's copy is
  # penalized, which is exactly what the model computed before.
  if(length(dim(get("ln_InitDevs", envir = env))) == 3) {
    init3 <- get("ln_InitDevs", envir = env)
    set("ln_InitDevs", array(rep(init3, get("n_sexes", envir = env)), dim = c(dim(init3), get("n_sexes", envir = env))))
  }
  if(!has("init_devs_pen_use")) {
    pen_use <- array(0, dim = dim(get("ln_InitDevs", envir = env)))
    pen_use[,,,1] <- 1
    set("init_devs_pen_use", pen_use)
  }

  # Selectivity penalty weights are now one specification per fleet. Older input
  # lists hold a single named vector shared by every fleet, so replicate it.
  for(nm in c("fish_sel_pen_wts", "ret_sel_pen_wts", "srv_sel_pen_wts")) {
    if(!has(nm)) next
    spec <- get(nm, envir = env)
    if(is.null(names(spec))) {
      # Already per fleet, but predates the per-term normalize switch
      spec <- lapply(spec, function(s) { if(is.null(s$normalize)) s$normalize <- TRUE; s })
      set(nm, spec)
      next
    }
    n_fleets_bc <- if(nm == "srv_sel_pen_wts") n_srv_bc else n_fish_bc
    spec <- as.list(spec)
    spec$normalize <- TRUE
    set(nm, rep(list(spec), n_fleets_bc))
  } # end nm loop

  invisible(NULL)
}


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
  n_ages = length(ages) # number of ages
  n_yrs = length(years) # number of years
  n_lens = length(lens) # number of lengths
  # length compositions recorded on coarser bins than the model carries are
  # mapped through LenBinMap inside the likelihood, the same way an ageing error
  # matrix maps model ages onto observed bins; NA leaves the bins as they are
  LenBinMap_lik = if(is.null(LenBinMap)) NA else LenBinMap
  LenBinMap_fn = if(is.null(LenBinMap)) NULL else function(y) LenBinMap

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
  Fish_caal = if(do_caal == 1) array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_fish_fleets)) else NULL # Retained catch at length and age
  Fish_caal_discard = if(do_caal == 1) array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_fish_fleets)) else NULL # Discarded catch at length and age
  PredCatch = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted retained catch (can be abundance or biomass)
  PredDiscard = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted discarded catch (can be abundance, biomass, or abdunance or biomass fraction of retained catch)
  PredFishIdx = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_fish_fleets)) # Predicted fishery index
  fish_q = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery catchability

  # Survey Processes
  SrvIAA = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes, n_srv_fleets)) # Survey index at age
  SrvIAL = array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets)) # Survey index at length
  Srv_caal = if(do_caal == 1) array(data = 0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_lens, n_ages, n_sexes, n_srv_fleets)) else NULL # Survey index at length and age
  PredSrvIdx = array(0, dim = c(n_pop, n_regions, n_yrs, n_seas, n_srv_fleets)) # Predicted survey index
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
  Fish_caal_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_fish_fleets)) # Fishery Conditional Age-at-Length Likelihoods
  Srv_caal_nLL = array(0, dim = c(n_regions, n_yrs, n_seas, n_lens, n_sexes, n_srv_fleets)) # Survey Conditional Age-at-Length Likelihoods
  # Conditional age-at-length is on only where use flags say so
  do_fish_caal = !is.null(UseFish_caal) && any(UseFish_caal == 1)
  do_srv_caal = !is.null(UseSrv_caal) && any(UseSrv_caal == 1)
  if((do_fish_caal || do_srv_caal) && do_caal != 1)
    stop("Conditional age-at-length data are being fit, but do_caal is 0, so the joint arrays at length and age were never built. Set do_caal = 1 in Setup_Mod_Biologicals.")

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
    ctmc_diffusion_bounds = ctmc_diffusion_bounds,
    seasdur = seasdur,
    ctmc_scale_by_seasdur = ctmc_scale_by_seasdur
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

  ## Growth -------------------------------------------------------------------

  use_cohort_growth = growth_model != 0 && growth_tv_type == 1 # whether using growth cohort mode or not
  precomputed_mortality_yrs = if(use_cohort_growth) seq_len(max(1, growth_cohort_styr - 1)) else 1:n_yrs

  if(growth_model != 0) {

    # grab growth arguments into environment
    growth_args = growth_args_from_model()

    # get growth function
    tmp_growth = do.call(Get_Growth, c(growth_args, list(
      n_yrs = n_yrs, n_fish_fleets = n_fish_fleets, n_srv_fleets = n_srv_fleets,
      growth_tv_type = growth_tv_type, growth_cohort_styr = growth_cohort_styr,
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

  ## Selectivity ------------------------------------------------------------
  if(fish_selex_type == 0) fish_selex_bins = ages # if age-based selectivity
  if(fish_selex_type == 1) fish_selex_bins = lens # if length-based selectivity
  if(ret_selex_type == 0) ret_selex_bins = ages # if age-based selectivity
  if(ret_selex_type == 1) ret_selex_bins = lens # if length-based selectivity
  if(srv_selex_type == 0) srv_selex_bins = ages # if age-based selectivity
  if(srv_selex_type == 1) srv_selex_bins = lens # if length-based selectivity

  # Total fishery selectivity
  tmp_fish_sel = Get_Selex_Array(selex_type = fish_selex_type, bins = fish_selex_bins,
                                 sel_blocks = fish_sel_blocks, sel_model = fish_sel_model,
                                 fixed_sel_pars = fish_fixed_sel_pars, cont_tv_sel = cont_tv_fish_sel,
                                 ln_seldevs = ln_fishsel_devs, use_fixed_sel = use_fixed_fish_sel,
                                 bin_devs = ln_fishsel_bin_devs, bin_dev_bins = fish_sel_bin_dev_bins,
                                 sel_norm_bins = fish_sel_norm_bins,
                                 sex_par_offset = fishsel_sex_par_offset, sex_scale_offset = fishsel_sex_scale_offset, sex_apical_offset = fishsel_sex_apical_offset,
                                 sex_scale = ln_fishsel_sex_scale,
                                 nselbins = fish_sel_bicubic_nselbins,
                                 dbnrml_raw = fish_dbnrml_raw,
                                 dbnrml_startbin = fish_dbnrml_startbin,
                                 sel_input = fish_sel_input,
                                 bicubic_Wbin = fish_sel_bicubic_Wbin, bicubic_Wyr = fish_sel_bicubic_Wyr,
                                 bicubic_binnodes = fish_sel_bicubic_binnodes, bicubic_yrnodes = fish_sel_bicubic_yrnodes,
                                 n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_proj_yrs_devs = n_proj_yrs_devs,
                                 n_seas = n_seas, n_ages = n_ages, n_lens = n_lens, n_sexes = n_sexes,
                                 n_fleets = n_fish_fleets)
  fish_sel = tmp_fish_sel$sel
  fish_sel_l = tmp_fish_sel$sel_l

  # Fishery retention selectivity
  tmp_ret_sel = Get_Selex_Array(selex_type = ret_selex_type, bins = ret_selex_bins,
                                sel_blocks = ret_sel_blocks, sel_model = ret_sel_model,
                                fixed_sel_pars = ret_fixed_sel_pars, cont_tv_sel = cont_tv_ret_sel,
                                ln_seldevs = ln_retsel_devs, use_fixed_sel = use_fixed_ret_sel,
                                bin_devs = ln_retsel_bin_devs, bin_dev_bins = ret_sel_bin_dev_bins,
                                 sel_norm_bins = ret_sel_norm_bins,
                                sex_par_offset = retsel_sex_par_offset, sex_scale_offset = retsel_sex_scale_offset, sex_apical_offset = retsel_sex_apical_offset,
                                sex_scale = ln_retsel_sex_scale,
                                nselbins = ret_sel_bicubic_nselbins,
                                dbnrml_raw = ret_dbnrml_raw,
                                dbnrml_startbin = ret_dbnrml_startbin,
                                sel_input = ret_sel_input,
                                bicubic_Wbin = ret_sel_bicubic_Wbin, bicubic_Wyr = ret_sel_bicubic_Wyr,
                                bicubic_binnodes = ret_sel_bicubic_binnodes, bicubic_yrnodes = ret_sel_bicubic_yrnodes,
                                n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_proj_yrs_devs = n_proj_yrs_devs,
                                n_seas = n_seas, n_ages = n_ages, n_lens = n_lens, n_sexes = n_sexes,
                                n_fleets = n_fish_fleets)
  ret_sel = tmp_ret_sel$sel
  ret_sel_l = tmp_ret_sel$sel_l

  # Survey selectivity
  tmp_srv_sel = Get_Selex_Array(selex_type = srv_selex_type, bins = srv_selex_bins,
                                sel_blocks = srv_sel_blocks, sel_model = srv_sel_model,
                                fixed_sel_pars = srv_fixed_sel_pars, cont_tv_sel = cont_tv_srv_sel,
                                ln_seldevs = ln_srvsel_devs, use_fixed_sel = use_fixed_srv_sel,
                                bin_devs = ln_srvsel_bin_devs, bin_dev_bins = srv_sel_bin_dev_bins,
                                 sel_norm_bins = srv_sel_norm_bins,
                                sex_par_offset = srvsel_sex_par_offset, sex_scale_offset = srvsel_sex_scale_offset, sex_apical_offset = srvsel_sex_apical_offset,
                                sex_scale = ln_srvsel_sex_scale,
                                nselbins = srv_sel_bicubic_nselbins,
                                dbnrml_raw = srv_dbnrml_raw,
                                dbnrml_startbin = srv_dbnrml_startbin,
                                sel_input = srv_sel_input,
                                bicubic_Wbin = srv_sel_bicubic_Wbin, bicubic_Wyr = srv_sel_bicubic_Wyr,
                                bicubic_binnodes = srv_sel_bicubic_binnodes, bicubic_yrnodes = srv_sel_bicubic_yrnodes,
                                n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_proj_yrs_devs = n_proj_yrs_devs,
                                n_seas = n_seas, n_ages = n_ages, n_lens = n_lens, n_sexes = n_sexes,
                                n_fleets = n_srv_fleets)
  srv_sel = tmp_srv_sel$sel
  srv_sel_l = tmp_srv_sel$sel_l

  ## Mortality ---------------------------------------------------------------
  missing_catch = is.na(ObsCatch) # TRUE = aggregate catch is missing, not a true recorded zero

  # get mortality arguments for the model
  mortality_args = mortality_args_from_model()

  # Values needed for mortality function, and later update_cohort_growth_and_mortality
  mortality_state_fields = c("Fmort", "dmr", "fish_sel", "ret_sel", "ret_FAA", "disc_FAA", "tot_FAA", "ZAA", "WAA") # mortality states
  if(growth_model == 0) mortality_state_fields = c(mortality_state_fields, "SizeAgeTrans_fish", "SizeAgeTrans_srv") # add in sizeage if not using a growth model
  if(growth_model != 0) { # add in additional states if using growth model
    mortality_state_fields = c(mortality_state_fields, "tmp_growth", "SizeAgeTrans_fish", "SizeAgeTrans_srv", "SizeAgeTrans_spawn",
                          "mean_LAA_fish", "sd_LAA_fish", "mean_LAA_srv", "sd_LAA_srv", "mean_LAA_spawn", "sd_LAA_spawn",
                          "Linf", "L_beg", "growth_pars_y")
    if(derive_waa == 1) mortality_state_fields = c(mortality_state_fields, "WAA_fish", "WAA_srv")
  } # growth module
  mortality_state = mget(mortality_state_fields, envir = environment()) # output values to be easily referenced

  for(y in precomputed_mortality_yrs) mortality_state = do.call(compute_mortality_year, c(mortality_args, list(y = y, st = mortality_state))) # get mortality values
  list2env(mortality_state, envir = environment()) # update mortality values in environment

  ## Growth x Mortality (cohort growth) --------------------------------------
  # build out cohort influenced mortality, since mortality dynamics can change if large cohorts enter
  update_cohort_growth_and_mortality = NULL
  if(use_cohort_growth) {

    update_cohort_growth_and_mortality = function(y, NAA_y, st) {
      if(y >= growth_cohort_styr) {
        st$tmp_growth = do.call(Get_Growth_Year, c(growth_args, list(growth = st$tmp_growth, y = y, NAA_y = NAA_y)))
        st = growth_take_year(st, st$tmp_growth, y, derive_waa) # this year's growth quantity
        # Only WAA_fish/WAA_srv genuinely need this year's growth; ZAA is
        # recomputed as a side effect and only actually changes if selectivity
        # is length-based (age-based selectivity reproduces the same ZAA).
        st = do.call(compute_mortality_year, c(mortality_args, list(y = y, st = st))) # get new selex after cohort growth, then compute mortality
      }

      list(state = st,
           ZAA_y = st$ZAA[,,y,,,, drop = FALSE],
           WAA_y = st$WAA[,,y,,,, drop = FALSE],
           MatAA_y = MatAA[,,y,,,, drop = FALSE])
    }
  }

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

  ## Initial Age Structure ---------------------------------------------------
  # Get initial F for age structure
  catch_flag_base = array(UseCatch[,1,,], dim = c(n_regions, n_seas, n_fish_fleets))
  catch_flag_pop = apply(UseCatch_pop[,,1,,,drop = FALSE], c(2,4,5), max)
  catch_flag = pmax(catch_flag_base, catch_flag_pop)

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
    natmort = array(natmort[,,1,,], dim = c(n_pop, n_regions, n_ages, n_sexes)), # natural mortality in first year
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
    move_timing = move_timing
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
    dmr = array(0, dim = c(n_regions, n_seas, n_fish_fleets)), # unfished
    Mrate = if(is.null(Mrate)) NULL else array(Mrate[,,,1,,,], dim = c(n_pop, n_regions, n_regions, n_seas, n_ages, n_sexes)), # rates in first year
    move_timing = move_timing
  )

  # Input into model arrays (first year and season) - and add lognormal mean adjustment
  NAA[,,1,1,,] = Init_Fished_NAA
  NAA0[,,1,1,,] = Init_Unfished_NAA

  ## Population Projection ---------------------------------------------------

  sr_R0 = if(sr_R0_spec == 1) exp(ln_sr_R0) else if(sr_R0_spec == 2) rinit else R0

  tmp_pop_proj = get_population_projection(
    n_pop = n_pop, n_regions = n_regions, n_seas = n_seas, n_ages = n_ages, n_sexes = n_sexes,
    n_yrs = n_yrs, n_fish_fleets = n_fish_fleets, n_est_rec_devs = n_est_rec_devs,
    rec_lag = rec_lag, rec_model = rec_model, rec_dd = rec_dd, R0 = R0,
    rec_region_prop = rec_region_prop, rec_seas_prop = rec_seas_prop, h_trans = h_trans,
    natal_region = natal_region, t_spawn = t_spawn, spawn_seas = spawn_seas, seasdur = seasdur,
    init_F = init_F, ln_RecDevs = ln_RecDevs,
    sexratio = sexratio, WAA = WAA, MatAA = MatAA, natmort = natmort, Movement = Movement,
    stray_rate = stray_rate, sgl_seas_spawning_movement = sgl_seas_spawning_movement,
    do_recruits_move = do_recruits_move, fish_sel = fish_sel, ret_sel = ret_sel, dmr = dmr, ZAA = ZAA,
    NAA = NAA, NAA0 = NAA0, NAA_bef = NAA_bef, NAA_aft = NAA_aft, Rec = Rec,
    SSB = SSB, Total_Biom = Total_Biom, Dynamic_SSB0 = Dynamic_SSB0, eff_SSB = eff_SSB,
    Mrate = Mrate, move_timing = move_timing, SR_ref_yr = SR_ref_yr,
    sr_penalty = sr_penalty,
    sr_R0 = sr_R0,
    growth_mortality_year_fn = update_cohort_growth_and_mortality,
    growth_mortality_state = mortality_state
  )

  # the years update_cohort_growth_and_mortality built come back through the state it carried
  if(use_cohort_growth) {
    mortality_state = tmp_pop_proj$growth_mortality_state
    Fmort = mortality_state$Fmort; dmr = mortality_state$dmr; fish_sel = mortality_state$fish_sel
    ret_sel = mortality_state$ret_sel; ret_FAA = mortality_state$ret_FAA; disc_FAA = mortality_state$disc_FAA
    tot_FAA = mortality_state$tot_FAA; ZAA = mortality_state$ZAA
    if(growth_model != 0) {
      tmp_growth = mortality_state$tmp_growth
      SizeAgeTrans_fish = mortality_state$SizeAgeTrans_fish; SizeAgeTrans_srv = mortality_state$SizeAgeTrans_srv
      SizeAgeTrans_spawn = mortality_state$SizeAgeTrans_spawn
      mean_LAA_fish = mortality_state$mean_LAA_fish; sd_LAA_fish = mortality_state$sd_LAA_fish
      mean_LAA_srv = mortality_state$mean_LAA_srv; sd_LAA_srv = mortality_state$sd_LAA_srv
      mean_LAA_spawn = mortality_state$mean_LAA_spawn; sd_LAA_spawn = mortality_state$sd_LAA_spawn
      Linf = mortality_state$Linf; L_beg = mortality_state$L_beg; growth_pars_y = mortality_state$growth_pars_y
      if(derive_waa == 1) { WAA = mortality_state$WAA; WAA_fish = mortality_state$WAA_fish; WAA_srv = mortality_state$WAA_srv }
    } # growth module
  } # cohort growth

  NAA = tmp_pop_proj$NAA; NAA0 = tmp_pop_proj$NAA0; NAA_bef = tmp_pop_proj$NAA_bef; NAA_aft = tmp_pop_proj$NAA_aft
  NAA_int = tmp_pop_proj$NAA_int # season-integrated abundance (move_timing == 2 only)
  Rec = tmp_pop_proj$Rec; SSB = tmp_pop_proj$SSB; Total_Biom = tmp_pop_proj$Total_Biom
  Dynamic_SSB0 = tmp_pop_proj$Dynamic_SSB0; eff_SSB = tmp_pop_proj$eff_SSB
  Aggregated_SSB = tmp_pop_proj$Aggregated_SSB; Dynamic_Aggregated_SSB0 = tmp_pop_proj$Dynamic_Aggregated_SSB0
  SR_pred = tmp_pop_proj$SR_pred

  ## Fishery Observation Model -----------------------------------------------
  tmp_fish_obs = get_fishery_observation_model(
    n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_fish_fleets = n_fish_fleets, n_sexes = n_sexes,
    fish_q_blocks = fish_q_blocks, ln_fish_q = ln_fish_q, fish_q = fish_q,
    ret_FAA = ret_FAA, disc_FAA = disc_FAA, ZAA = ZAA, NAA = NAA,
    CAA = CAA, DAA = DAA, CAL = CAL, DAL = DAL, PredCatch = PredCatch, PredDiscard = PredDiscard, PredFishIdx = PredFishIdx,
    fit_lengths = fit_lengths, SizeAgeTrans = SizeAgeTrans,
    SizeAgeTrans_fish = SizeAgeTrans_fish,
    fish_len_comp_sel = fish_len_comp_sel, fish_selex_type = fish_selex_type, ret_selex_type = ret_selex_type,
    fish_sel_l = fish_sel_l, ret_sel_l = ret_sel_l, Fmort = Fmort,
    catch_units = catch_units, discard_units = discard_units, WAA_fish = WAA_fish, dmr = dmr,
    fish_idx_type = fish_idx_type, fish_sel = fish_sel, ret_sel = ret_sel,
    Mrate = Mrate, move_timing = move_timing, seasdur = seasdur,
    NAA_int = if(move_timing == 2) NAA_int else NULL, # reuse the dynamics' exponentials
    t_fish = t_fish, fish_idx_ages = fish_idx_ages,
    fish_q_type = fish_q_type, do_fish_q_cov = do_fish_q_cov,
    fish_q_cov = fish_q_cov, fish_q_coeff = fish_q_coeff,
    ObsFishIdx = ObsFishIdx, UseFishIdx = UseFishIdx,
    do_caal = do_caal, Fish_caal = Fish_caal, Fish_caal_discard = Fish_caal_discard
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
    n_pop = n_pop, n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_srv_fleets = n_srv_fleets, n_sexes = n_sexes,
    srv_q_blocks = srv_q_blocks, ln_srv_q = ln_srv_q, srv_q = srv_q, do_srv_q_cov = do_srv_q_cov, srv_q_cov = srv_q_cov, srv_q_coeff = srv_q_coeff,
    srv_selex_type = srv_selex_type, srv_sel = srv_sel, srv_sel_l = srv_sel_l, SizeAgeTrans = SizeAgeTrans,
    SizeAgeTrans_srv = SizeAgeTrans_srv,
    srv_len_comp_sel = srv_len_comp_sel,
    NAA = NAA, ZAA = ZAA, t_srv = t_srv, SrvIAA = SrvIAA, fit_lengths = fit_lengths, SrvIAL = SrvIAL,
    srv_idx_type = srv_idx_type, WAA_srv = WAA_srv, PredSrvIdx = PredSrvIdx,
    Mrate = Mrate, move_timing = move_timing, seasdur = seasdur,
    srv_idx_ages = srv_idx_ages, srv_q_type = srv_q_type,
    ObsSrvIdx = ObsSrvIdx, UseSrvIdx = UseSrvIdx,
    RecDev_anom = RecDev_anom,
    do_caal = do_caal, Srv_caal = Srv_caal
  )

  srv_q = tmp_srv_obs$srv_q; srv_sel = tmp_srv_obs$srv_sel
  SrvIAA = tmp_srv_obs$SrvIAA; SrvIAL = tmp_srv_obs$SrvIAL; PredSrvIdx = tmp_srv_obs$PredSrvIdx
  Srv_caal = tmp_srv_obs$Srv_caal


  ## Conventional Tagging Observation Model -----------------------------------------------
  if(any(use_conv_fish_tagging == 1)) {

    tmp_tag_obs = get_tagging_observation_model(
      n_fish_fleets = n_fish_fleets, n_regions = n_regions, n_conv_tag_cohorts = n_conv_tag_cohorts,
      n_yrs = n_yrs, n_seas = n_seas, n_pop = n_pop, n_ages = n_ages, n_sexes = n_sexes,
      conv_tag_fish_reporting_blocks = conv_tag_fish_reporting_blocks,
      conv_tag_fish_reporting_pars = conv_tag_fish_reporting_pars, conv_tag_fish_reporting = conv_tag_fish_reporting,
      conv_tag_release_indicator = conv_tag_release_indicator, conv_tag_max_liberty = conv_tag_max_liberty,
      use_conv_fish_tagging = use_conv_fish_tagging,
      Fmort = Fmort, fish_sel = fish_sel, ret_sel = ret_sel, dmr = dmr, natmort = natmort, seasdur = seasdur,
      ln_conv_tag_shed = ln_conv_tag_shed, conv_tag_t_tagging = conv_tag_t_tagging, conv_tagged_fish = conv_tagged_fish,
      conv_fish_tag_attr = conv_fish_tag_attr, conv_tag_release_platform = conv_tag_release_platform,
      srv_sel = srv_sel, NAA_bef = NAA_bef, ln_init_conv_tag_mort = ln_init_conv_tag_mort,
      do_recruits_move = do_recruits_move, Movement = Movement,
      conv_tag_fish_avail = conv_tag_fish_avail, pred_conv_tag_fish_recap = pred_conv_tag_fish_recap,
      Mrate = Mrate, move_timing = move_timing
    )

    conv_tag_fish_reporting = tmp_tag_obs$conv_tag_fish_reporting
    conv_tag_fish_avail = tmp_tag_obs$conv_tag_fish_avail
    pred_conv_tag_fish_recap = tmp_tag_obs$pred_conv_tag_fish_recap

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
  # Non-lognormal fleets go first, while ObsFishIdx is still untransformed.
  for(f in 1:n_fish_fleets) {
    if(FishIdx_LikeType[f] == 0) next
    use_f = array(UseFishIdx[,,,f], dim = dim(UseFishIdx)[1:3])
    if(!any(use_f == 1)) next

    obs_pos_f = which(use_f == 1)
    obs_map_f = arrayInd(obs_pos_f, dim(use_f))
    obs_vec_f = se_vec_f = rep(0, length(obs_pos_f))
    pred_vec_f = rep(0, length(obs_pos_f))

    for(i in seq_along(obs_pos_f)) {
      r    = obs_map_f[i, 1]
      y    = obs_map_f[i, 2]
      seas = obs_map_f[i, 3]
      obs_vec_f[i] = ObsFishIdx[r,y,seas,f]
      se_vec_f[i] = ObsFishIdx_SE[r,y,seas,f]
      pred_vec_f[i] = sum(PredFishIdx[,r,y,seas,f])
    } # end i loop

    tmp_fidx_nLL = get_index_nLL(obs_vec_f, pred_vec_f, se_vec_f, FishIdx_LikeType[f],
                                 FishIdx_Cov[[f]], addtofishidx)

    # input into likelihoods
    for(i in seq_along(obs_pos_f)) FishIdx_nLL[obs_map_f[i,1], obs_map_f[i,2], obs_map_f[i,3], f] = tmp_fidx_nLL[i]

  } # end f loop

  UseFishIdx_ln = UseFishIdx
  for(f in 1:n_fish_fleets) if(FishIdx_LikeType[f] != 0) UseFishIdx_ln[,,,f] = 0

  if(any(UseFishIdx_ln == 1)) { # setup OSA

    valid_idx_ir = which(UseFishIdx_ln == 1)
    ObsFishIdx_map = arrayInd(valid_idx_ir, dim(UseFishIdx_ln))
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
  # Arithmetic-scale normal fleets go first, while ObsFishIdx_pop is still
  # untransformed, and register no OSA observations, mirroring the regional
  # block. A multivariate normal covariance describes the regional series and
  # says nothing about the population stream, so mvn fleets keep lognormal
  # error here.
  for(f in 1:n_fish_fleets) {
    if(FishIdx_LikeType[f] != 1) next
    use_fp = array(UseFishIdx_pop[,,,,f], dim = dim(UseFishIdx_pop)[1:4])
    if(!any(use_fp == 1)) next

    obs_pos_fp = which(use_fp == 1)
    obs_map_fp = arrayInd(obs_pos_fp, dim(use_fp))

    for(i in seq_along(obs_pos_fp)) {
      p    = obs_map_fp[i, 1]
      r    = obs_map_fp[i, 2]
      y    = obs_map_fp[i, 3]
      seas = obs_map_fp[i, 4]

      FishIdx_pop_nLL[p,r,y,seas,f] = -1 * RTMB::dnorm(ObsFishIdx_pop[p,r,y,seas,f],
                                                       PredFishIdx[p,r,y,seas,f],
                                                       ObsFishIdx_pop_SE[p,r,y,seas,f], TRUE)
    }
  } # end f loop

  UseFishIdx_pop_ln = UseFishIdx_pop
  for(f in 1:n_fish_fleets) if(FishIdx_LikeType[f] == 1) UseFishIdx_pop_ln[,,,,f] = 0

  if(any(UseFishIdx_pop_ln == 1)) { # setup OSA

    valid_idx_ip = which(UseFishIdx_pop_ln == 1)
    ObsFishIdx_pop_map = arrayInd(valid_idx_ip, dim(UseFishIdx_pop_ln))
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
              comp_const_obs = comp_const_obs,
              Exp = apply(CAA[,,y,seas,,,f, drop = FALSE], 2:7, sum),
              Obs = ObsFishAgeComps[,y,seas,,,f], ISS = ISS_FishAgeComps[,y,seas,,f], Wt_Mltnml = Wt_FishAgeComps[,y,seas,,f],
              Comp_Type = FishAgeComps_Type[y,f], Likelihood_Type = FishAgeComps_LikeType[f], ln_theta = ln_FishAge_theta[,,f],
              ln_theta_agg = ln_FishAge_theta_agg[f], LN_corr_pars = FishAge_corr_pars[,,f,], LN_corr_pars_agg = FishAge_corr_pars_agg[f],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError[y,,], use = UseFishAgeComps[,y,seas,f],
              n_model_bins = n_ages, n_obs_bins = dim(ObsFishAgeComps)[4], addtocomp = addtocomp,
              comp_bins = if(all(FishAgeComps_bins[,f] == 1)) NULL else which(FishAgeComps_bins[,f] == 1)
            )
          } # if we have fishery age comps

          # Fishery Length Compositions
          if(sum(UseFishLenComps[,y,seas,f]) >= 1 && fit_lengths == 1) {
            FishLenComps_nLL[,y,seas,,f] = Get_Comp_Likelihoods(
              comp_const_obs = comp_const_obs,
              Exp = apply(CAL[,,y,seas,,,f, drop = FALSE], 2:7, sum), Obs = ObsFishLenComps[,y,seas,,,f],
              ISS = ISS_FishLenComps[,y,seas,,f], Wt_Mltnml = Wt_FishLenComps[,y,seas,,f], Comp_Type = FishLenComps_Type[y,f],
              Likelihood_Type = FishLenComps_LikeType[f], ln_theta = ln_FishLen_theta[,,f], ln_theta_agg = ln_FishLen_theta_agg[f],
              LN_corr_pars = FishLen_corr_pars[,,f,], LN_corr_pars_agg = FishLen_corr_pars_agg[f],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = LenBinMap_lik,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
                comp_const_obs = comp_const_obs,
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
                comp_const_obs = comp_const_obs,
                Exp = CAL[p,,y,seas,,,f], Obs = ObsFishLenComps_pop[p,,y,seas,,,f], ISS = ISS_FishLenComps_pop[p,,y,seas,,f],
                Wt_Mltnml = Wt_FishLenComps_pop[p,,y,seas,,f], Comp_Type = FishLenComps_pop_Type[y,f],
                Likelihood_Type = FishLenComps_pop_LikeType[f], ln_theta = ln_FishLen_pop_theta[p,,,f],
                ln_theta_agg = ln_FishLen_pop_theta_agg[p,f], LN_corr_pars = FishLen_pop_corr_pars[p,,,f,],
                LN_corr_pars_agg = FishLen_pop_corr_pars_agg[p,f], n_regions = n_regions,
                n_sexes = n_sexes, age_or_len = 1, AgeingError = LenBinMap_lik, use = UseFishLenComps_pop[p,,y,seas,f],
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
              comp_const_obs = comp_const_obs,
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
              comp_const_obs = comp_const_obs,
              Exp = apply(DAL[,,y,seas,,,f, drop = FALSE], 2:7, sum), Obs = ObsFishLenComps_discard[,y,seas,,,f],
              ISS = ISS_FishLenComps_discard[,y,seas,,f], Wt_Mltnml = Wt_FishLenComps_discard[,y,seas,,f],
              Comp_Type = FishLenComps_discard_Type[y,f], Likelihood_Type = FishLenComps_discard_LikeType[f],
              ln_theta = ln_FishLen_discard_theta[,,f], ln_theta_agg = ln_FishLen_discard_theta_agg[f],
              LN_corr_pars = FishLen_discard_corr_pars[,,f,], LN_corr_pars_agg = FishLen_discard_corr_pars_agg[f],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = LenBinMap_lik,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
                comp_const_obs = comp_const_obs,
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
                comp_const_obs = comp_const_obs,
                Exp = DAL[p,,y,seas,,,f], Obs = ObsFishLenComps_discard_pop[p,,y,seas,,,f], ISS = ISS_FishLenComps_discard_pop[p,,y,seas,,f],
                Wt_Mltnml = Wt_FishLenComps_discard_pop[p,,y,seas,,f], Comp_Type = FishLenComps_discard_pop_Type[y,f],
                Likelihood_Type = FishLenComps_discard_pop_LikeType[f], ln_theta = ln_FishLen_discard_pop_theta[p,,,f],
                ln_theta_agg = ln_FishLen_discard_pop_theta_agg[p,f], LN_corr_pars = FishLen_discard_pop_corr_pars[p,,,f,],
                LN_corr_pars_agg = FishLen_discard_pop_corr_pars_agg[p,f], n_regions = n_regions, n_sexes = n_sexes,
                age_or_len = 1, AgeingError = LenBinMap_lik, use = UseFishLenComps_discard_pop[p,,y,seas,f], n_model_bins = n_lens,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsFishLenComps_discard_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }
  }

  ## Survey Likelihoods ------------------------------------------------------
  ### Survey Indices (Regional) ---------------------------------------------------------
  # Fleets on an arithmetic-scale normal or a multivariate normal are evaluated
  # first, while ObsSrvIdx still holds untransformed observations. They register
  # no OSA observations, so the lognormal block below keeps exactly the
  # observation vector it has always had.
  for(sf in 1:n_srv_fleets) {
    if(SrvIdx_LikeType[sf] == 0) next
    use_sf = array(UseSrvIdx[,,,sf], dim = dim(UseSrvIdx)[1:3])
    if(!any(use_sf == 1)) next

    obs_pos = which(use_sf == 1)
    obs_map = arrayInd(obs_pos, dim(use_sf))
    obs_vec = se_vec = rep(0, length(obs_pos))
    pred_vec = rep(0, length(obs_pos))

    for(i in seq_along(obs_pos)) {
      r    = obs_map[i, 1]
      y    = obs_map[i, 2]
      seas = obs_map[i, 3]
      obs_vec[i] = ObsSrvIdx[r,y,seas,sf]
      se_vec[i] = ObsSrvIdx_SE[r,y,seas,sf]
      pred_vec[i] = sum(PredSrvIdx[,r,y,seas,sf])
    } # end i loop

    # A multivariate normal cannot be split across years, so its total lands in
    # the first observation's cell and the rest stay zero.
    tmp_idx_nLL = get_index_nLL(obs_vec, pred_vec, se_vec, SrvIdx_LikeType[sf],
                                SrvIdx_Cov[[sf]], addtosrvidx)

    for(i in seq_along(obs_pos)) SrvIdx_nLL[obs_map[i,1], obs_map[i,2], obs_map[i,3], sf] = tmp_idx_nLL[i]

  } # end sf loop

  UseSrvIdx_ln = UseSrvIdx
  for(sf in 1:n_srv_fleets) if(SrvIdx_LikeType[sf] != 0) UseSrvIdx_ln[,,,sf] = 0

  if(any(UseSrvIdx_ln == 1)) { # setup OSA
    valid_idx_sr = which(UseSrvIdx_ln == 1)
    ObsSrvIdx_map = arrayInd(valid_idx_sr, dim(UseSrvIdx_ln))
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
  # Arithmetic-scale normal fleets go first, while ObsSrvIdx_pop is still
  # untransformed, and register no OSA observations, mirroring the regional
  # block. A multivariate normal covariance describes the regional series and
  # says nothing about the population stream, so mvn fleets keep lognormal
  # error here.
  for(sf in 1:n_srv_fleets) {
    if(SrvIdx_LikeType[sf] != 1) next
    use_sp = array(UseSrvIdx_pop[,,,,sf], dim = dim(UseSrvIdx_pop)[1:4])
    if(!any(use_sp == 1)) next

    obs_pos_sp = which(use_sp == 1)
    obs_map_sp = arrayInd(obs_pos_sp, dim(use_sp))

    for(i in seq_along(obs_pos_sp)) {
      p    = obs_map_sp[i, 1]
      r    = obs_map_sp[i, 2]
      y    = obs_map_sp[i, 3]
      seas = obs_map_sp[i, 4]

      SrvIdx_pop_nLL[p,r,y,seas,sf] = -1 * RTMB::dnorm(ObsSrvIdx_pop[p,r,y,seas,sf],
                                                       PredSrvIdx[p,r,y,seas,sf],
                                                       ObsSrvIdx_pop_SE[p,r,y,seas,sf], TRUE)
    }
  } # end sf loop

  UseSrvIdx_pop_ln = UseSrvIdx_pop
  for(sf in 1:n_srv_fleets) if(SrvIdx_LikeType[sf] == 1) UseSrvIdx_pop_ln[,,,,sf] = 0

  if(any(UseSrvIdx_pop_ln == 1)) { # setup OSA
    valid_idx_sp = which(UseSrvIdx_pop_ln == 1)
    ObsSrvIdx_pop_map = arrayInd(valid_idx_sp, dim(UseSrvIdx_pop_ln))
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
              comp_const_obs = comp_const_obs,
              Exp = apply(SrvIAA[,,y,seas,,,sf, drop = FALSE], 2:7, sum), Obs = ObsSrvAgeComps[,y,seas,,,sf],
              ISS = ISS_SrvAgeComps[,y,seas,,sf], Wt_Mltnml = Wt_SrvAgeComps[,y,seas,,sf],
              Comp_Type = SrvAgeComps_Type[y,sf], Likelihood_Type = SrvAgeComps_LikeType[sf],
              ln_theta = ln_SrvAge_theta[,,sf], ln_theta_agg = ln_SrvAge_theta_agg[sf],
              LN_corr_pars = SrvAge_corr_pars[,,sf,], LN_corr_pars_agg = SrvAge_corr_pars_agg[sf],
              n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError[y,,],
              use = UseSrvAgeComps[,y,seas,sf], n_model_bins = n_ages, n_obs_bins = dim(ObsSrvAgeComps)[4], addtocomp = addtocomp,
              comp_bins = if(all(SrvAgeComps_bins[,sf] == 1)) NULL else which(SrvAgeComps_bins[,sf] == 1)
            )
          }

          # Survey Length Compositions
          if(sum(UseSrvLenComps[,y,seas,sf]) >= 1 && fit_lengths == 1) {
            SrvLenComps_nLL[,y,seas,,sf] = Get_Comp_Likelihoods(
              comp_const_obs = comp_const_obs,
              Exp = apply(SrvIAL[,,y,seas,,,sf, drop = FALSE], 2:7, sum), Obs = ObsSrvLenComps[,y,seas,,,sf],
              ISS = ISS_SrvLenComps[,y,seas,,sf], Wt_Mltnml = Wt_SrvLenComps[,y,seas,,sf], Comp_Type = SrvLenComps_Type[y,sf],
              Likelihood_Type = SrvLenComps_LikeType[sf], ln_theta = ln_SrvLen_theta[,,sf], ln_theta_agg = ln_SrvLen_theta_agg[sf],
              LN_corr_pars = SrvLen_corr_pars[,,sf,], LN_corr_pars_agg = SrvLen_corr_pars_agg[sf], n_regions = n_regions,
              n_sexes = n_sexes, age_or_len = 1, AgeingError = LenBinMap_lik, use = UseSrvLenComps[,y,seas,sf],
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
                comp_const_obs = comp_const_obs,
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
                comp_const_obs = comp_const_obs,
                Exp = SrvIAL[p,,y,seas,,,sf], Obs = ObsSrvLenComps_pop[p,,y,seas,,,sf],
                ISS = ISS_SrvLenComps_pop[p,,y,seas,,sf], Wt_Mltnml = Wt_SrvLenComps_pop[p,,y,seas,,sf],
                Comp_Type = SrvLenComps_pop_Type[y,sf], Likelihood_Type = SrvLenComps_pop_LikeType[sf],
                ln_theta = ln_SrvLen_pop_theta[p,,,sf], ln_theta_agg = ln_SrvLen_pop_theta_agg[p,sf],
                LN_corr_pars = SrvLen_pop_corr_pars[p,,,sf,], LN_corr_pars_agg = SrvLen_pop_corr_pars_agg[p,sf],
                n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = LenBinMap_lik,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
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
          AgeingErrorFn = LenBinMap_fn, addtocomp = addtocomp,
          family = "continuous", zero_init = is.null(ObsSrvLenComps_pop_osa_discrete), pop = TRUE, n_pop = n_pop
        )
      }
    }
  }

  ## Conditional Age-at-Length Likelihoods ------------------------------------
  # helpers to sum across population arrays
  caal_exp = function(arr, y, seas, f) {
    e = arr[,,y,seas,,,,f, drop = FALSE]
    if(n_pop == 1) dim(e) = c(n_regions, n_lens, n_ages, n_sexes) else e = apply(e, c(2,5,6,7), sum)
    return(e)
  }
  caal_exp_len = function(arr, y, seas, l, f) {
    e = arr[,,y,seas,l,,,f, drop = FALSE]
    if(n_pop == 1) dim(e) = c(n_regions, n_ages, n_sexes) else e = apply(e, c(2,6,7), sum)
    return(e)
  }

  ### Fishery Conditional Age-at-Length ---------------------------------------
  if(do_fish_caal) {
    if(do_internal_comp_osa == FALSE) {

      for(y in 1:n_yrs) {
        for(f in 1:n_fish_fleets) {
          for(seas in 1:n_seas) {
            if(sum(UseFish_caal[,y,seas,,f]) >= 1) {
              Fish_caal_nLL[,y,seas,,,f] = Get_CAAL_Likelihoods(
                comp_const_obs = comp_const_obs,
                Exp = caal_exp(Fish_caal, y, seas, f), Obs = ObsFish_caal[,y,seas,,,,f],
                ISS = ISS_Fish_caal[,y,seas,,,f], Wt_Mltnml = Wt_Fish_caal[,y,seas,,,f],
                ln_theta = ln_Fish_caal_theta[,,f], ln_theta_agg = ln_Fish_caal_theta_agg[f],
                Comp_Type = Fish_caal_Type[y,f], Likelihood_Type = Fish_caal_LikeType[f],
                n_regions = n_regions, n_lens = n_lens, n_model_bins = n_ages,
                n_obs_bins = dim(ObsFish_caal)[5], n_sexes = n_sexes,
                AgeingError = AgeingError[y,,], use = UseFish_caal[,y,seas,,f],
                addtocomp = addtocomp
              )
            } # if we have fishery caal
          } # end seas loop
        } # end f loop
      } # end y loop

    } else {

      # Using OSA compositions
      ObsFish_caal_osa = pack_caal_osa(
        ObsArr = ObsFish_caal, ISSArr = ISS_Fish_caal, WtArr = Wt_Fish_caal,
        UseArr = UseFish_caal, TypeMat = Fish_caal_Type, LikeTypeVec = Fish_caal_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_lens = n_lens, n_fleets = n_fish_fleets,
        n_sexes = n_sexes, addtocomp = addtocomp
      )

      if(!is.null(ObsFish_caal_osa)) {
        ObsFish_caal_osa = RTMB::OBS(ObsFish_caal_osa)
        Fish_caal_nLL = eval_caal_osa(
          nLL_arr = Fish_caal_nLL, tracked = ObsFish_caal_osa,
          ExpArrFn = function(y, seas, l, f) caal_exp_len(Fish_caal, y, seas, l, f),
          UseArr = UseFish_caal, TypeMat = Fish_caal_Type, LikeTypeVec = Fish_caal_LikeType,
          ISSArr = ISS_Fish_caal, lnThetaArr = ln_Fish_caal_theta, lnThetaAggVec = ln_Fish_caal_theta_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_lens = n_lens,
          n_fleets = n_fish_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsFish_caal)[5],
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp
        )
      }

    } # end osa switch
  } # end fishery caal

  ### Survey Conditional Age-at-Length ----------------------------------------
  if(do_srv_caal) {
    if(do_internal_comp_osa == FALSE) {

      for(y in 1:n_yrs) {
        for(sf in 1:n_srv_fleets) {
          for(seas in 1:n_seas) {
            if(sum(UseSrv_caal[,y,seas,,sf]) >= 1) {
              Srv_caal_nLL[,y,seas,,,sf] = Get_CAAL_Likelihoods(
                comp_const_obs = comp_const_obs,
                Exp = caal_exp(Srv_caal, y, seas, sf), Obs = ObsSrv_caal[,y,seas,,,,sf],
                ISS = ISS_Srv_caal[,y,seas,,,sf], Wt_Mltnml = Wt_Srv_caal[,y,seas,,,sf],
                ln_theta = ln_Srv_caal_theta[,,sf], ln_theta_agg = ln_Srv_caal_theta_agg[sf],
                Comp_Type = Srv_caal_Type[y,sf], Likelihood_Type = Srv_caal_LikeType[sf],
                n_regions = n_regions, n_lens = n_lens, n_model_bins = n_ages,
                n_obs_bins = dim(ObsSrv_caal)[5], n_sexes = n_sexes,
                AgeingError = AgeingError[y,,], use = UseSrv_caal[,y,seas,,sf],
                addtocomp = addtocomp
              )
            } # if we have survey caal
          } # end seas loop
        } # end sf loop
      } # end y loop

    } else {

      # Using OSA compositions
      ObsSrv_caal_osa = pack_caal_osa(
        ObsArr = ObsSrv_caal, ISSArr = ISS_Srv_caal, WtArr = Wt_Srv_caal,
        UseArr = UseSrv_caal, TypeMat = Srv_caal_Type, LikeTypeVec = Srv_caal_LikeType,
        n_yrs = n_yrs, n_seas = n_seas, n_lens = n_lens, n_fleets = n_srv_fleets,
        n_sexes = n_sexes, addtocomp = addtocomp
      )

      if(!is.null(ObsSrv_caal_osa)) {
        ObsSrv_caal_osa = RTMB::OBS(ObsSrv_caal_osa)
        Srv_caal_nLL = eval_caal_osa(
          nLL_arr = Srv_caal_nLL, tracked = ObsSrv_caal_osa,
          ExpArrFn = function(y, seas, l, sf) caal_exp_len(Srv_caal, y, seas, l, sf),
          UseArr = UseSrv_caal, TypeMat = Srv_caal_Type, LikeTypeVec = Srv_caal_LikeType,
          ISSArr = ISS_Srv_caal, lnThetaArr = ln_Srv_caal_theta, lnThetaAggVec = ln_Srv_caal_theta_agg,
          n_regions = n_regions, n_yrs = n_yrs, n_seas = n_seas, n_lens = n_lens,
          n_fleets = n_srv_fleets, n_sexes = n_sexes,
          n_model_bins = n_ages, n_obs_bins = dim(ObsSrv_caal)[5],
          AgeingErrorFn = function(y) AgeingError[y,,], addtocomp = addtocomp
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

  ## Priors and Penalties ----------------------------------------------------
  ### Fishing Mortality (Penalty) ---------------------------------------------
  if(Use_F_pen == 1) {
    Fmort_nLL = Get_Fdev_PE_loglik(PE_model = Fdev_model, ln_sigmaF = ln_sigmaF, Fdev_rho = Fdev_rho,
                                   ln_F_devs = ln_F_devs, map_ln_F_devs = map_ln_F_devs, Fdev_pen_center = Fdev_pen_center)
  } #  if using fishing mortality penalty

  ### Discard Mortality Rate (Penalty) ---------------------------------------------
  if(Use_dmr_pen == 1) {
    dmr_nLL = get_dmr_penalty(logit_dmr_devs = logit_dmr_devs, ln_sigma_dmr = ln_sigma_dmr,
                              map_logit_dmr_devs = map_logit_dmr_devs,
                              n_fish_fleets = n_fish_fleets, n_yrs = n_yrs, n_regions = n_regions, n_seas = n_seas)
  } #  if using discard mortality rate penalty


  ### Selectivity (Penalty) ---------------------------------------------------
  for(r in 1:n_regions) {

    for(f in 1:n_fish_fleets) {

      # Total Fishery Selectivity Deviations
      if(cont_tv_fish_sel[r,f] > 0 && fishsel_pe_wt[f] != 0) {

        sel_nLL = sel_nLL + - fishsel_pe_wt[f] * Get_PE_loglik(PE_model = cont_tv_fish_sel[r,f], # process error model
                                                PE_pars = fishsel_pe_pars[r,,,f, drop = FALSE], # process error parameters for a given fleet (correlaiton and sigmas)
                                                ln_devs = ln_fishsel_devs[r,,,,f, drop = FALSE], # extract out process error deviations for a given fleet
                                                map_sel_devs = map_ln_fishsel_devs[r,,,,f, drop = FALSE],
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
                                                min_sel_devs_shared_bins = srvsel_devs_min_shared_bins,
                                                rw_init_sigma = srvsel_rw_init_sigma[sf]

        )
      } # end if

    } # end sf loop
  } # end r loop


  ### Bin-Override Selectivity Deviations (Process Error) ---------------------
  # The bin dimension of these arrays is n_ages for age-based selectivity and
  # n_lens for length-based, so the reference bins are taken from the array
  # itself rather than assumed to be ages. Inert while setup_sel_bin_devs
  # restricts these to none/iid/rw, since min_sel_devs_shared_bins is only read
  # by the 3D separable process errors, but wrong the moment that is relaxed.
  # The override deviations carry their own process error, so a bin that is free
  # of the functional form is still smoothed in time when asked to be.
  for(r in 1:n_regions) {
    for(f in 1:n_fish_fleets) {
      if(cont_tv_fishsel_bin_devs[f] > 0) {
        sel_nLL = sel_nLL + - Get_PE_loglik(PE_model = cont_tv_fishsel_bin_devs[f],
                                                PE_pars = fishsel_bin_devs_pe_pars[r,,,f, drop = FALSE],
                                                ln_devs = ln_fishsel_bin_devs[r,,,,f, drop = FALSE],
                                                map_sel_devs = map_ln_fishsel_bin_devs[r,,,,f, drop = FALSE],
                                                min_sel_devs_shared_bins = 1:dim(ln_fishsel_bin_devs)[3],
                                                rw_init_sigma = fishsel_bin_devs_rw_init_sigma[f])
      } # end if
      if(cont_tv_retsel_bin_devs[f] > 0) {
        sel_nLL = sel_nLL + - Get_PE_loglik(PE_model = cont_tv_retsel_bin_devs[f],
                                                PE_pars = retsel_bin_devs_pe_pars[r,,,f, drop = FALSE],
                                                ln_devs = ln_retsel_bin_devs[r,,,,f, drop = FALSE],
                                                map_sel_devs = map_ln_retsel_bin_devs[r,,,,f, drop = FALSE],
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
                                                min_sel_devs_shared_bins = 1:dim(ln_srvsel_bin_devs)[3],
                                                rw_init_sigma = srvsel_bin_devs_rw_init_sigma[sf])
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
      has_nonzero_pen = any(unlist(lapply(smooth_pen_terms, function(nm) safe_extract(fish_sel_pen_wts[[f]], nm))) != 0)

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
      has_nonzero_pen = any(unlist(lapply(smooth_pen_terms, function(nm) safe_extract(ret_sel_pen_wts[[f]], nm))) != 0)

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
      has_nonzero_pen = any(unlist(lapply(smooth_pen_terms, function(nm) safe_extract(srv_sel_pen_wts[[sf]], nm))) != 0)
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


  ### Selectivity (Prior) -----------------------------------------------------
  # Each row of a prior table constrains either one fixed selectivity parameter
  # (type "par", the default) or the realized selectivity value at one bin
  # (type "value"), per get_selex_prior.
  # Total Fishery selectivity
  if(Use_fish_selex_prior == 1) sel_nLL = sel_nLL + get_selex_prior(fish_selex_prior, fish_fixed_sel_pars, fish_sel, fish_sel_l, fish_selex_type, fish_sel_blocks)

  # Retained Fishery selectivity
  if(Use_ret_selex_prior == 1) sel_nLL = sel_nLL + get_selex_prior(ret_selex_prior, ret_fixed_sel_pars, ret_sel, ret_sel_l, ret_selex_type, ret_sel_blocks)

  # Survey selectivity
  if(Use_srv_selex_prior == 1) sel_nLL = sel_nLL + get_selex_prior(srv_selex_prior, srv_fixed_sel_pars, srv_sel, srv_sel_l, srv_selex_type, srv_sel_blocks)

  ### Selectivity Parameter Centering (Penalty) -------------------------------
  if(Use_fish_selex_penalty == 1) sel_nLL = sel_nLL + get_selex_fixed_penalty(fish_selex_penalty, fish_fixed_sel_pars)
  if(Use_ret_selex_penalty == 1) sel_nLL = sel_nLL + get_selex_fixed_penalty(ret_selex_penalty, ret_fixed_sel_pars)
  if(Use_srv_selex_penalty == 1) sel_nLL = sel_nLL + get_selex_fixed_penalty(srv_selex_penalty, srv_fixed_sel_pars)

  ### Recruitment (Penalty) ----------------------------------------------------
  tmp_rec_pen = get_recruitment_penalty(n_pop = n_pop, n_regions = n_regions, n_ages = n_ages,
                                        n_est_rec_devs = n_est_rec_devs, rec_dd = rec_dd, natal_region = natal_region,
                                        rec_region_prop_spec = rec_region_prop_spec, rec_region_prop = rec_region_prop,
                                        equil_init_age_strc = equil_init_age_strc, ln_InitDevs = ln_InitDevs,
                                        init_age_devs_shared = init_age_devs_shared, ln_sigmaR = ln_sigmaR,
                                        bias_ramp = bias_ramp, sigmaR_switch = sigmaR_switch, ln_RecDevs = ln_RecDevs,
                                        sigmaR2_early = sigmaR2_early, sigmaR2_late = sigmaR2_late,
                                        do_rec_bias_ramp = do_rec_bias_ramp,
                                        map_ln_RecDevs = map_ln_RecDevs,
                                        RecDevs_pen_center = RecDevs_pen_center,
                                        InitDevs_pen_center = InitDevs_pen_center,
                                        init_devs_pen_use = init_devs_pen_use,
                                        init_bias_ramp = init_bias_ramp,
                                        map_ln_InitDevs = map_ln_InitDevs,
                                        Use_init_sex_pen = Use_init_sex_pen, ln_sigma_init_sex = ln_sigma_init_sex)
  Init_Rec_nLL = tmp_rec_pen$Init_Rec_nLL
  Init_Sex_nLL = tmp_rec_pen$Init_Sex_nLL
  Rec_nLL = tmp_rec_pen$Rec_nLL

  ### Recruitment Level (Penalty) ---------------------------------------------
  # A second, separate statement about recruitment: the stock-recruit penalty
  # constrains the residuals about the curve, this one constrains the series.
  Rec_level_nLL = array(0, dim = dim(Rec))
  if(Use_rec_level_pen == 1) {
    Rec_level_nLL = get_rec_level_penalty(Rec, exp(ln_sigma_rec_level), rec_level_pen_center,
                                          if(all(rec_level_pen_yrs == 1)) NULL else which(rec_level_pen_yrs == 1))
  }

  ### Stock-Recruit Residual (Penalty) ----------------------------------------
  # Only reachable under mean recruitment. The curve scores the recruitment
  # series without generating it, so the deviations stay free and a weakly
  # determined relationship informs the series rather than dictating it.
  SR_pen_nLL = array(0, dim = dim(Rec))
  if(sr_penalty > 0) {
    SR_pen_nLL = get_sr_penalty(Rec, SR_pred, exp(ln_sigma_sr_pen),
                                if(all(sr_pen_yrs == 1)) NULL else which(sr_pen_yrs == 1))
  }

  ### Growth (Process Error) --------------------------------------------------
  # Growth carries two deviation surfaces, both scored through the same process
  # error machinery the selectivity deviations use. A time-varying growth
  # PARAMETER has one series per parameter, a surface one column wide, iid or a
  # random walk. Semi-parametric growth has one surface over years and ages,
  # which may also be a 3D GMRF or a separable 2D AR(1).
  if(growth_model != 0) {
    for(p in 1:n_pop) {
      for(r in 1:n_regions) {

        if(!is.null(growth_tv_model)) for(k in which(growth_tv_model > 0)) {
          growth_tv_nLL = growth_tv_nLL - Get_PE_loglik(
            PE_model = growth_tv_model[k],
            PE_pars = array(growth_pe_pars[p,r,k,,1], dim = c(1, 1, n_sexes, 1)),
            ln_devs = array(ln_growth_devs[p,r,,k,], dim = c(1, n_yrs, 1, n_sexes, 1)),
            map_sel_devs = array(map_ln_growth_devs[p,r,,k,], dim = c(1, n_yrs, 1, n_sexes)),
            min_sel_devs_shared_bins = 1,
            rw_init_sigma = growth_rw_init_sigma
          )
        } # end k loop

        if(growth_semipar > 0) {
          growth_semipar_nLL = growth_semipar_nLL - Get_PE_loglik(
            PE_model = growth_semipar,
            PE_pars = array(growth_pe_pars[p,r,,,2], dim = c(1, dim(growth_pe_pars)[3], n_sexes, 1)),
            ln_devs = array(ln_growth_semipar_devs[p,r,,,], dim = c(1, n_yrs, n_ages, n_sexes, 1)),
            map_sel_devs = array(map_ln_growth_semipar_devs[p,r,,,], dim = c(1, n_yrs, n_ages, n_sexes)),
            min_sel_devs_shared_bins = if(is.null(growth_semipar_bins)) 1:n_ages else growth_semipar_bins,
            rw_init_sigma = growth_rw_init_sigma
          )
        }

      } # end r loop
    } # end p loop
  } # end growth process error

  ### Initial Recruitment Offset (Penalty) ---------------------------------------
  # The initial equilibrium recruitment against the recruitment level, as a
  # normal on the log ratio; an equilibrium recruitment averages over several
  # years of deviations, so its sd is typically sigmaR shrunk by that span
  if(Use_rinit_pen == 1 && use_rinit == 1) {
    for(p in 1:n_pop) rinit_nLL = rinit_nLL - RTMB::dnorm(ln_rinit[p] - ln_global_R0[p], 0, rinit_pen_sd, TRUE)
  }

  ### Fishery Catchability (Prior) -----------------------------------------------
  if(Use_fish_q_prior == 1) fish_q_nLL = fish_q_nLL + get_q_prior(fish_q_prior, ln_fish_q)

  ### Survey Catchability (Prior) -----------------------------------------------
  if(Use_srv_q_prior == 1) srv_q_nLL = srv_q_nLL + get_q_prior(srv_q_prior, ln_srv_q)

  ### Natural Mortality (Prior) -----------------------------------------------
  if(Use_M_prior == 1) M_nLL = M_nLL + get_natmort_prior(M_prior, ln_M, M_blocks)

  ### Steepness (Prior) -----------------------------------------------
  if(Use_h_prior == 1) h_nLL = h_nLL + get_steepness_prior(h_prior, h_trans)


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
  # For CTMC movement the prior is placed on the annual fractions exp(Q), so that alpha
  # means the same thing regardless of how the year is divided into seasons.
  if(Use_Movement_Prior == 1) Movement_nLL = Movement_nLL + get_movement_dirichlet_prior(Movement_prior, Movement, Mrate)

  ### Recruitment R0 and Proportions (Prior) -----------------------------------------
  # Regional/seasonal apportionment and stray rate priors all feed rec_prop_nLL
  rec_prop_nLL = get_recruitment_proportion_priors(
    use_rec_region_prop_prior = use_rec_region_prop_prior, rec_region_prop_prior = rec_region_prop_prior, rec_region_prop = rec_region_prop,
    use_rec_seas_prop_prior = use_rec_seas_prop_prior, use_fixed_rec_seas_prop = use_fixed_rec_seas_prop, rec_seas_prop_prior = rec_seas_prop_prior,
    rec_seas_prop = rec_seas_prop, rec_lag = rec_lag, spawn_seas = spawn_seas, n_seas = n_seas,
    use_stray_rate_prior = use_stray_rate_prior, stray_rate_prior = stray_rate_prior, stray_rate_pars = stray_rate_pars
  )

  if(use_r0_prior == 1) R0_nLL = get_r0_prior(r0_prior, ln_global_R0) # recruitment R0 (global scalar prior, unweighted by Wt_Rec)

  ### Tag Reporting Rate (Prior) --------------------------------------------
  if(use_conv_tag_fishrep_prior == 1) TagRep_nLL = TagRep_nLL + get_tagrep_prior(conv_tag_fishrep_prior, conv_tag_fish_reporting_pars)

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
  RTMB::REPORT(NAA_int)
  RTMB::REPORT(ZAA)
  RTMB::REPORT(natmort)
  RTMB::REPORT(bias_ramp)
  RTMB::REPORT(init_bias_ramp)
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

  # Survey Processes
  RTMB::REPORT(PredSrvIdx)
  RTMB::REPORT(srv_sel)
  RTMB::REPORT(srv_q)
  RTMB::REPORT(SrvIAA)
  RTMB::REPORT(SrvIAL)
  if(do_caal == 1) RTMB::REPORT(Srv_caal)
  RTMB::REPORT(RecDev_anom)


  # Per-fleet keys: growth-derived or fixed data, reported whenever either
  # source supplied one, so a fixed SizeAgeTrans_fish/SizeAgeTrans_srv shows up
  # in the report the same way a supplied WAA_fish/WAA_srv already does
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
  if(do_fish_caal) RTMB::REPORT(Fish_caal_nLL)
  if(do_srv_caal) RTMB::REPORT(Srv_caal_nLL)

  # Penalties and priors
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
  RTMB::REPORT(growth_tv_nLL)
  RTMB::REPORT(growth_semipar_nLL)
  RTMB::REPORT(rinit_nLL)
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
