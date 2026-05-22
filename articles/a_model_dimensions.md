# Description of Model and Data Dimensions

The following tables describes all elements contained within
`input_list$data`, which is generated using the `SPoRC::Setup_x`
functions. Note that they do not detail how the `SPoRC::Setup_x`
functions should be utilized. Rather, they detail the components of
`input_list$data`. Thus, for further details on how arguments for
`SPoRC::Setup_x` functions should be defined, users should refer to the
function documentation. Furthermore, note that when `n_sexes > 1`, the
first dimension will always be females and the second dimension will
always be males.

## Data Inputs for Defining Model Dimensions

| Name | Description |
|----|----|
| `years` | Vector specifying number of years to model |
| `ages` | Vector specifying number of ages to model |
| `lens` | Vector specifying number of lengths to model |
| `n_regions` | Value specifying number of regions to model |
| `n_sexes` | Value specifying number of sexes to model |
| `n_fish_fleets` | Value specifying number of fishery fleets to model |
| `n_srv_fleets` | Value specifying number of survey fleets to model |
| `n_proj_yrs_devs` | Number of years to project deviations / random effect parameters forward |

## Data Inputs for Defining Recruitment Processes

| Name | Description |
|----|----|
| `rec_lag` | Value specifying the delay between spawning and when recruits enter the population. For example, if recruits enter the population as age 2, rec_lag would be specified as 2, such that the spawning biomass from year – 2 produces these recruits |
| `Use_h_prior` | Value specifying whether steepness priors are used. 0: Don’t use priors, 1: Use Priors. Steepness priors are bounded between 0.2 and 1 with a scaled beta penalty. |
| `h_prior` | Data frame specifying prior distributions for steepness parameters. Must include columns: `region` (region index), `mu` (mean of the prior in normal space), and `sd` (standard deviation of the prior in normal space). For each row, a beta distribution is scaled to the interval \[0.2, 1\], and the corresponding `h_trans` value is transformed and penalized using the log-density from the scaled beta distribution. |
| `do_rec_bias_ramp` | Value specifying whether or not the Methot and Taylor recruitment bias ramp is conducted. 0: Don’t do bias ramp, 1: Do bias ramp |
| `max_bias_ramp_fct` | Value specifying the maximum bias correction factor applied to the recruitment bias ramp. |
| `bias_year` | Vector of years, specifying when to change the bias ramp. Can be specified with an NA if do_rec_bias_ramp = 0 |
| `sigmaR_switch` | Value specifying when to transition between using an early period sigma and late period sigma R for penalizing initial age deviations and recruitment deviations. Specify as 0 if sigmaR early is equal to sigmaR late. |
| `init_age_strc` | Value specifying how the population age structure should be initialized. 0: Initialize via iteration to equilibrium, 1: Initialize with geometric series solution. |
| `equil_init_age_strc` | Value specifying how initial age deviations arise. 0 == equilibrium, 1 == stochastic except for plus group, 2 == stochastic for all ages (including the plus group) |
| `init_F_prop` | Value specifying the proportion of fishing mortality to apply to the initial age deviations relative to the mean fishing mortality parameter for fishery fleet 1 (the dominant fleet) |
| `rec_model` | Value specifying the recruitment model. 0 == Mean Recruitment, 1 == Beverton-Holt with steepness parameterization |
| `rec_dd` | Value specifying the recruitment density dependence (only used when there is a stock recruitment relationship) 0 == local density, 1 == global density, 999 == no density dependent stock recruitment form is used |
| `t_spawn` | Fraction of year in which spawning occurs |
| `Use_Rec_prop_Prior` | Value specifying whether recruitment proportion priors are used. 0: Don’t use priors, 1: Use Priors. Recruitment proportion priors are applied using a Dirichlet distribution. |
| `Rec_prop_prior` | Prior values defined for recruitment proportions if used (i.e., the Dirichlet concentration parameters) dimensioned by `n_regions` |
| `sexratio_blocks` | Matrix specifying sexratio blocks (`n_regions x n_years`) |

## Data Inputs for Defining Biological Processes

| Name | Description |
|----|----|
| `WAA` | Weight-at-age values dimensioned by n_regions, n_years, n_ages, n_sexes |
| `WAA_fish` | Weight-at-age values for the fishery dimensioned by n_regions, n_years, n_ages, n_sexes, n_fish_fleets |
| `WAA_srv` | Weight-at-age values for the survey dimensioned by n_regions, n_years, n_ages, n_sexes, n_srv_fleets |
| `MatAA` | Maturity-at-age values dimensioned by n_regions, n_years, n_ages, n_sexes |
| `AgeingError` | Ageing error dimensioned by number of years, number of modelled ages, number of observed composition ages, where rows across modelled ages sum to 1 |
| `fit_lengths` | Value describing whether or not to fit length composition data. 0: Don’t fit lengths, 1: Fit lengths |
| `SizeAgeTrans` | Size-age transition matrix dimensioned by n_regions, n_years, n_lens, n_ages, n_sexes. Can be specified as NA if length compositions are not fit |
| `addtocomp` | Constant to add to all composition data. |
| `addtofishidx` | Constant to add to all fishery index data. |
| `addtosrvidx` | Constant to add to all survey index data. |
| `addtotag` | Constant to add to all tagging data. |
| `Use_M_prior` | Value specifying whether natural mortality priors are used. 0: Don’t use prior, 1: Use prior |
| `M_prior` | A dataframe specifying how natural mortality priors should be applied. |
| `Selex_Type` | Value specifying whether age or length-based selectivity is used. 0: Age-based, 1: Length-based |
| `use_fixed_natmort` | Value specifying whether a fixed mortality array is used. 0: Don’t use and estimate M, 1: Use |
| `Fixed_natmort` | Natural mortality array dimensioned by n_regions, n_years, n_ages, and n_sexes. |
| `M_blocks` | Array that specifies the dimensions in which natural mortality is blocked (n_regions, n_years, n_ages, and n_sexes) |

## Data Inputs for Defining Movement Processes

| Name | Description |
|----|----|
| `do_recruits_move` | Value specifying whether recruits are allowed to move. 0: Recruits don’t move, 1: Recruits move |
| `use_fixed_movement` | Value specifying whether or not to use a fixed movement matrix. 0: Don’t use fixed movement, 1: Use fixed movement |
| `Fixed_Movement` | Fixed movement matrix dimensioned by n_regions, n_regions, n_years, n_ages, n_sexes |
| `Use_Movement_Prior` | Value specifying whether or not to use movement priors. 0: Don’t use movement prior, 1: Use movement prior |
| `Movement_prior` | Data frame with columns representing origin region (`region_from`), representative age (`age`), sex (`sex`), and a list column (`alpha`) containing Dirichlet prior parameters. Each element of `alpha` is a numeric vector of length *n_regions*, specifying the prior concentration (`alpha` values) for movement among destination regions. |
| `cont_vary_movement` | Integer indicating whether movement is continuously varying across regions, years, and ages (0 = none, 1 = iid deviations) |
| `map_move_devs` | Array dimensioned by n_regions, n_regions - 1, n_years, n_ages, and n_sexes specifying which movement parameters are shared and mapped off |
| `move_type` | Integer indicating movement model type: 0 = unstructured Markov (multinomial logit), 1 = Continuous Time Markov Chain (CTMC) |
| `ctmc_move_dat` | Data frame with CTMC covariates (`regions`, `years`, `ages`, `sexes`, plus formula variables) used to build design matrices for diffusion and preference when `move_type == 1`. Can include projection years with projected covariate values. |
| `diffusion_formula` | R formula specifying diffusion covariates for CTMC movement (when `move_type == 1`) |
| `preference_formula` | R formula specifying preference covariates for CTMC movement (when `move_type == 1`) |
| `adjacency_mat` | Square adjacency matrix (`n_regions` × `n_regions`) defining allowed transitions between regions for CTMC movement |
| `adjacency_collapsed` | Collasped adjacency matrix (`n_regions` × `n_regions - 1`) defining allowed transitions between regions for CTMC movement (but excluding diagonals) |
| `area_r` | Numeric vector of region areas (length `n_regions`) used to scale diffusion rates in CTMC movement |
| `ctmc_diffusion_bounds` | Integer flag (0 or 1) indicating whether diffusion bounds are applied to ensure the generator matrix is Metzler (1 = bounds enforced) |

## Data Inputs for Defining Tagging Processes

| Name | Description |
|----|----|
| `UseTagging` | Value indicating whether to use tagging data. 0: Don’t use tagging data, 1: Use tagging data |
| `tag_release_indicator` | Matrix dimensioned by n_tag_cohorts, 2. The 2 columns in this matrix represent the tag region (column 1) for a given cohort (rows of the matrix) and the tag year of the cohort (column 2). |
| `n_tag_cohorts` | Value specifying the number of tag cohorts released |
| `max_tag_liberty` | Value specifying the maximum tag liberty to track cohorts |
| `Tagged_Fish` | Array dimensioned by n_tag_cohorts, n_ages, n_sexes for the number of tagged fish from a given cohort |
| `Obs_Tag_Recap` | Array dimensioned by max_tag_liberty, n_tag_cohorts, n_regions, n_ages, n_sexes of individuals that are recaptured. If no age or sex information is available to distinguish recaptured individuals, users can initialize the array with 0s and the recaptured individuals should be input into the first dimension of n_ages and the first dimension of n_sexes. This is then coupled with move_age_tag_pool and move_sex_tag_pool (see last two rows in this table for an example) to ensure that the predicted recaptures are summed across ages and sexes when fitting to the observed recaptures for a given likelihood. |
| `Tag_LikeType` | Value specifying the tag likelihood to use. 0: Poisson, 1: Negative Binomial, 2: Multinomial release conditioned, 3: Multinomial recaptured conditioned |
| `mixing_period` | Value specifying the mixing period for tag cohorts. This still requires the full history of the tag cohort to be input, but it ignores the mixing period in the likelihood calculations. |
| `t_tagging` | Fraction of year in which tagging occurs. Specify at 1 if tagging happens at the start of the year such that mortality is not discounted |
| `Use_TagRep_Prior` | Value specifying whether or not a tag reporting rate prior should be used. 0: Don’t use prior, 1: Use prior |
| `TagRep_Prior` | Data frame containing prior specifications for tag reporting parameters. Must include columns: `region` (region index), `block` (time block index), `mu` (numeric mean for tag reporting prior in normal space; use `NA` if symmetric beta is used), `sd` (numeric standard deviation for the prior in normal space), and `type` (0 = symmetric beta, 1 = regular beta). Each row specifies a beta prior for one tag reporting parameter. Only parameters with rows in this data frame will have priors applied. |
| `tag_selex` | Value specifying how tag selectivity is specified. 0 == Uniform with F resulting from Fleet 1, 1 == Sex Aggregated Selectivity with F resulting from Fleet 1, 2 == Sex-Specific Selectivity with F resulting from Fleet 1, 3 == Uniform with F resulting from a weighted sum of all fleets, 4 == Sex Aggregated Selectivity with F resulting from a weighted sum of all fleets, 5 == Sex-Specific Selectivity with F resulting from a weighted sum of all fleets |
| `tag_natmort` | Value specifying how tag natural mortality is specified. 0 == Natural mortality aggregated across ages and sexes, 1 == Natural mortality is age-specific but sex-aggregated, 2 == Natural mortality is age aggregated but sex-specific, 3 == Natural mortality is age and sex-specific |
| `move_age_tag_pool` | Data list object detailing how tag recaptures along the age axis should be pooled. If we have age information for tagging data, but are fitting two age blocks (10 total ages) and want to pool these two age blocks together when fitting to reduce computational cost, this can be specified as: `move_age_tag_pool <- list(c(1:5), c(6:10))` where the first element in the list specifies the ages to sum or pool across when, while the second element in the list specifies the age groups to pool across when fitting. Conversely, if we don’t have an age information, this list would be specified as: `move_age_tag_pool <- list(c(1:30))` where we are summing across all ages for the observed and predicted tag recaptures when fitting to these data |
| `move_sex_tag_pool` | Data list object detailing how tag recaptures along the sex axis should be pooled. If we have sex information for tagging data and are estimating sex-specific movement (2 sexes in this case), this would be specified as: `move_sex_tag_pool <- list(1, 2)` where the first element in the list specifies the first sex and the second element in the list specifies the second sex (i.e., no pooling in this case). Conversely, if we don’t have sex information, this list would be specified as: `move_sex_tag_pool <- list(c(1:2))` where we are summing across all sexes for the observed and predicted tag recaptures when fitting to these data |
| `Tag_Reporting_blocks` | Array dimensioned by n_regions and n_years indicating the regions and years for which a block is specified |

## Data Inputs for Defining Catch and Fishing Mortality Processes

| Name | Description |
|----|----|
| `ObsCatch` | Observed catch data dimensioned by n_regions, n_years, n_fish_fleets |
| `Catch_Type` | Matrix dimensioned by n_years, n_fish_fleets describing when and for what fishery fleet catch is aggregated across regions (0) or spatially explicit (1). For example, in a case where we have 10 years and 1 fishery fleet, where fishery fleet 1 is spatially aggregated in the first 5 years and fishery fleet 2 is spatially explicit, this would be structured as: `Catch_Type <- array(c(rep(0, 5), rep(1, 5)), dim = c(n_years, n_fish_fleets))` |
| `catch_units` | Array dimensioned by n_regions, n_fish_fleets describing catch units. 0 == Abundance, 1 == Biomass (default) |
| `UseCatch` | Array dimensioned by n_regions, n_years, n_fish_fleets describing whether to fit to catch data in a given year for a given fleet. 0: Don’t fit catch, 1: Fit catch |
| `est_all_regional_F` | Value specifying whether regional fishing mortality deviations are all estimated or only a subset are. 0: A subset of fishing mortality deviations are spatially aggregated, 1: All fishing mortality deviations are spatially explicit, irrespective of whether catch is aggregated across regions in certain periods |
| `Use_F_pen` | Value specifying whether to use a fishing mortality penalty to regularize fishing mortality deviations. 0: Don’t use regularity penalty, 1: Use regularity penalty |

## Data Inputs for Defining Fishery Indices and Compositions

| Name | Description |
|----|----|
| `ObsFishIdx` | Fishery index dimensioned by n_regions, n_years, n_fish_fleets |
| `ObsFishIdx_SE` | Fishery index standard errors dimensioned by n_regions, n_years, n_fish_fleets |
| `UseFishIdx` | Array dimensioned by n_regions, n_years, n_fish_fleets describing whether to fit to fishery index in a given year for a given fleet. 0: Don’t fit fishery index, 1: Fit fishery index |
| `fish_idx_type` | Matrix dimensioned by n_regions, n_fish_fleets specifying the index type for a given region and fishery fleet. 0: Abundance index, 1: Biomass index (uses `WAA` for calculations), 999: None Available |
| `ObsFishAgeComps` | Observed fishery age compositions dimensioned by n_regions, n_years, n_ages, n_sexes, n_fish_fleets. This can be either input as proportions or as numbers, since these are normalized again within the model |
| `UseFishAgeComps` | Array dimensioned by n_regions, n_years, n_fish_fleets specifying whether or not to fit to fishery age compositions. 0: Don’t fit fishery ages, 1: Fit fishery ages |
| `ISS_FishAgeComps` | Array dimensioned by n_regions, n_years, n_fish_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood. This is used in conjunction with `Wt_FishAgeComps` or `ln_FishAge_theta` if a Dirichlet-multinomial is used |
| `Wt_FishAgeComps` | Array dimensioned by n_regions, n_sexes, n_fish_fleets specifying a multinomial weight to apply to fishery age compositions, which should ideally be derived using Francis re-weighting methods |
| `ObsFishLenComps` | Observed fishery length compositions dimensioned by n_regions, n_years, n_lens, n_sexes, n_fish_fleets. This can be either input as proportions or as numbers, since these are normalized again within the model |
| `UseFishLenComps` | Array dimensioned by n_regions, n_years, n_fish_fleets specifying whether or not to fit to fishery length compositions. 0: Don’t fit fishery lengths, 1: Fit fishery lengths |
| `ISS_FishLenComps` | Array dimensioned by n_regions, n_years, n_fish_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood. This is used in conjunction with `Wt_FishLenComps` or `ln_FishLen_theta` if a Dirichlet-multinomial is used |
| `FishAgeComps_LikeType` | Vector dimensioned by n_fish_fleets specifying the likelihood to use for fitting fishery age compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |
| `FishLenComps_LikeType` | Vector dimensioned by n_fish_fleets specifying the likelihood to use for fitting fishery length compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |
| `FishAgeComps_Type` | Matrix dimensioned by n_years, n_fish_fleets specifying how age composition data should be structured and fit to. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |
| `FishLenComps_Type` | Matrix dimensioned by n_years, n_fish_fleets specifying how length composition data should be structured and fit to. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |

## Data Inputs for Defining Survey Indices and Compositions

| Name | Description |
|----|----|
| `ObsSrvIdx` | Survey index dimensioned by n_regions, n_years, n_srv_fleets |
| `ObsSrvIdx_SE` | Survey index standard errors dimensioned by n_regions, n_years, n_srv_fleets |
| `UseSrvIdx` | Array dimensioned by n_regions, n_years, n_srv_fleets describing whether to fit to survey index in a given year for a given fleet. 0: Don’t fit survey index, 1: Fit survey index |
| `srv_idx_type` | Matrix dimensioned by n_regions, n_srv_fleets specifying the index type for a given region and survey fleet. 0: Abundance index, 1: Biomass index (uses `WAA` for calculations), 999: None Available. |
| `ObsSrvAgeComps` | Observed survey age compositions dimensioned by n_regions, n_years, n_ages, n_sexes, n_srv_fleets. This can be either input as proportions or as numbers, since these are normalized again within the model |
| `UseSrvAgeComps` | Array dimensioned by n_regions, n_years, n_srv_fleets specifying whether or not to fit to survey age compositions. 0: Don’t fit survey ages, 1: Fit survey ages |
| `ISS_SrvAgeComps` | Array dimensioned by n_regions, n_years, n_srv_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood. This is used in conjunction with `Wt_SrvAgeComps` or `ln_SrvAge_theta` if a Dirichlet-multinomial is used |
| `Wt_SrvAgeComps` | Array dimensioned by n_regions, n_sexes, n_srv_fleets specifying a multinomial weight to apply to survey age compositions, which should ideally be derived using Francis re-weighting methods |
| `ObsSrvLenComps` | Observed survey length compositions dimensioned by n_regions, n_years, n_lens, n_sexes, n_srv_fleets. This can be either input as proportions or as numbers, since these are normalized again within the model |
| `UseSrvLenComps` | Array dimensioned by n_regions, n_years, n_srv_fleets specifying whether or not to fit to survey length compositions. 0: Don’t fit survey lengths, 1: Fit survey lengths |
| `ISS_SrvLenComps` | Array dimensioned by n_regions, n_years, n_srv_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood. This is used in conjunction with `Wt_SrvLenComps` or `ln_SrvLen_theta` if a Dirichlet-multinomial is used |
| `SrvAgeComps_LikeType` | Vector dimensioned by n_srv_fleets specifying the likelihood to use for fitting survey age compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |
| `SrvLenComps_LikeType` | Vector dimensioned by n_srv_fleets specifying the likelihood to use for fitting survey length compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |
| `SrvAgeComps_Type` | Matrix dimensioned by n_years, n_srv_fleets specifying how age composition data should be structured and fit to. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |
| `SrvLenComps_Type` | Matrix dimensioned by n_years, n_srv_fleets specifying how length composition data should be structured and fit to. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. More options and details can be found in `Get_Comp_Likelihoods.R` |
| `t_srv` | Matrix dimensioned by n_years, n_srv_fleets specifying survey timing as a proportion of the year. |

## Data Inputs for Defining Fishery Selectivity and Catchability

| Name | Description |
|----|----|
| `cont_tv_fish_sel` | Matrix dimensioned by n_regions, n_fish_fleets specifying whether or not to do continuous time-varying selectivity, and the type of continuous time-varying selectivity to do. 0: Don’t do continuous time-varying selectivity, 1: iid continuous time-varying selectivity, 2: random walk time-varying selectivity, 3: 3d gmrf continuous time-varying selectivity with marginal variance (semi-parametric), 4: 3d gmrf continuous time-varying selectivity with conditional variance (semi-parametric). More options and details on this can be found in `Get_PE_loglik.R`, `Get_3d_precision.R,` and `Get_Selex.R` |
| `fish_sel_blocks` | Array dimensioned by n_regions, n_years, n_fish_fleets specifying when a new selectivity might be specified. Unique numbers denote when a given fishery selectivity block should be specified. For example, if we have 1 region, 1 fishery fleet and 10 years, and we have two evenly spaced time blocks, this would be specified as: `fish_sel_blocks <- array(c(rep(0,5),rep(1,5)), dim = c(n_regions, n_years, n_fish_fleets))` where each unique number represents the selectivity block pattern to apply in that given period. |
| `fish_sel_model` | Array dimensioned by n_regions, n_years, n_fish_fleets specifying the selectivity pattern for a given region, year, and fishery fleet. Several selectivity forms are currently available. 0: Logistic with a50 and slope, 1: Gamma dome-shaped, 2: Power function, 3: Logistic with a50 and a95, 4: Double Normal with 6 parameters. More options and details can be found in `Get_Selex.R` and the model equations vignette. |
| `fish_q_blocks` | Array dimensioned by n_regions, n_years, n_fish_fleets specifying when a catchability block might be specified. Unique numbers denote when a new catchability block should be specified. For example, if we have 1 region, 1 fishery fleet and 10 years, and we have two evenly spaced catchability blocks, this would be specified as: `fish_q_blocks <- array(c(rep(0,5),rep(1,5)), dim = c(n_regions, n_years, n_fish_fleets))` where each unique number represents the unique catchability estimate to apply in that given period. |
| `Use_fish_q_prior` | Fishery catchability prior indicator 0 == don’t use, 1 == use |
| `fish_q_prior` | Data frame containing prior specifications for fishery catchability parameters. Must include columns: `region` (region index), `fleet` (fleet index), `block` (time block index), `mu` (prior mean on natural scale), and `sd` (prior standard deviation on log scale). Each row specifies a log-normal prior N(log(mu), sd) for one catchability parameter. |
| `map_ln_fishsel_devs` | Array of values dimensioned by n_regions, n_years, n_ages, n_sexes, and n_fish_fleets indicating which values are mapped off |
| `Use_fish_selex_prior` | Integer (0 or 1). Flag to enable/disable fishery selectivity priors. When set to 1, applies log-normal priors to fishery selectivity parameters as specified in `fish_selex_prior`. When set to 0, no priors are applied. |
| `fish_selex_prior` | Data frame containing prior specifications for fishery selectivity parameters. Must include columns: `region` (region index), `fleet` (fleet index), `block` (time block index), `sex` (sex index), `par` (parameter index), `mu` (prior mean on natural scale), and `sd` (prior standard deviation on log scale). Each row specifies a log-normal prior N(log(mu), sd) for one selectivity parameter. |
| `cont_tv_fish_sel_penalty` | Boolean on whether fishery selectivity penalties are applied to continous time-varying selectivity |

## Data Inputs for Defining Survey Selectivity and Catchability

| Name | Description |
|----|----|
| `cont_tv_srv_sel` | Matrix dimensioned by n_regions, n_srv_fleets specifying whether or not to do continuous time-varying selectivity, and the type of continuous time-varying selectivity to do. 0: Don’t do continuous time-varying selectivity, 1: iid continuous time-varying selectivity, 2: random walk time-varying selectivity, 3: 3d gmrf continuous time-varying selectivity with marginal variance (semi-parametric), 4: 3d gmrf continuous time-varying selectivity with conditional variance (semi-parametric). More options and details on this can be found in `Get_PE_loglik.R`, `Get_3d_precision.R,` and `Get_Selex.R` |
| `srv_sel_blocks` | Array dimensioned by n_regions, n_years, n_srv_fleets specifying when a new selectivity might be specified. Unique numbers denote when a given survey selectivity block should be specified. For example, if we have 1 region, 1 survey fleet and 10 years, and we have two evenly spaced time blocks, this would be specified as: `srv_sel_blocks <- array(c(rep(0,5),rep(1,5)), dim = c(n_regions, n_years, n_srv_fleets))` where each unique number represents the selectivity block pattern to apply in that given period. |
| `srv_sel_model` | Array dimensioned by n_regions, n_years, n_srv_fleets specifying the selectivity pattern for a given region, year, and survey fleet. Several selectivity forms are currently available. 0: Logistic with a50 and slope, 1: Gamma dome-shaped, 2: Power function, 3: Logistic with a50 and a95, 4: Double Normal with 6 parameters. More options and details can be found in `Get_Selex.R` and the model equations vignette. |
| `srv_q_blocks` | Array dimensioned by n_regions, n_years, n_srv_fleets specifying when a catchability block might be specified. Unique numbers denote when a new catchability block should be specified. For example, if we have 1 region, 1 survey fleet and 10 years, and we have two evenly spaced catchability blocks, this would be specified as: `srv_q_blocks <- array(c(rep(0,5),rep(1,5)), dim = c(n_regions, n_years, n_srv_fleets))` where each unique number represents the unique catchability estimate to apply in that given period. |
| `Use_srv_q_prior` | Survey catchability prior indicator 0 == don’t use, 1 == use |
| `srv_q_prior` | Data frame containing prior specifications for survey catchability parameters. Must include columns: `region` (region index), `fleet` (fleet index), `block` (time block index), `mu` (prior mean on natural scale), and `sd` (prior standard deviation on log scale). Each row specifies a log-normal prior N(log(mu), sd) for one catchability parameter. |
| `map_ln_srvsel_devs` | Array of values dimensioned by n_regions, n_years, n_ages, n_sexes, and n_srv_fleets indicating which values are mapped off |
| `do_srv_q_cov` | Numeric value indicating whether a survey catchability covariate is used in the model. 0 == Not used, 1 == Used |
| `srv_q_cov` | Array dimensioned by `n_regions`, `n_years`, `n_srv_fleets`, `n_covariates` representing the covariates used to compute survey catchability. |
| `Use_srv_selex_prior` | Integer (0 or 1). Flag to enable/disable survey selectivity priors. When set to 1, applies log-normal priors to survey selectivity parameters as specified in `srv_selex_prior`. When set to 0, no priors are applied. |
| `srv_selex_prior` | Data frame containing prior specifications for survey selectivity parameters. Must include columns: `region` (region index), `fleet` (fleet index), `block` (time block index), `sex` (sex index), `par` (parameter index), `mu` (prior mean on natural scale), and `sd` (prior standard deviation on log scale). Each row specifies a log-normal prior N(log(mu), sd) for one selectivity parameter. |
| `cont_tv_srv_sel_penalty` | Boolean on whether survey selectivity penalties are applied to continous time-varying selectivity |

## Data Inputs for Defining Model Weighting

| Name | Description |
|----|----|
| `Wt_Catch` | Weight applied to fishery catch likelihoods. Either a numeric scalar or an array of scalars dimensioned by `n_regions`, `n_years`, `n_fish_fleets`. |
| `Wt_FishIdx` | Weight applied to fishery index likelihoods. Either a numeric scalar or an array of scalars dimensioned by `n_regions`, `n_years`, `n_fish_fleets`. |
| `Wt_SrvIdx` | Weight applied to survey index likelihoods. Either a numeric scalar or an array of scalars dimensioned by `n_regions`, `n_years`, `n_srv_fleets`. |
| `Wt_Rec` | Weight applied to initial age deviations and recruitment deviations. |
| `Wt_F` | Weight applied to fishing mortality deviations. |
| `Wt_Tagging` | Weight applied to tagging data. |
| `Wt_FishAgeComps` | Array dimensioned by n_regions, n_years, n_sexes, n_fish_fleets specifying a multinomial weight to apply to fishery age compositions, which should ideally be derived using Francis re-weighting methods |
| `Wt_FishLenComps` | Array dimensioned by n_regions, n_years, n_sexes, n_fish_fleets specifying a multinomial weight to apply to fishery length compositions, which should ideally be derived using Francis re-weighting methods |
| `Wt_SrvAgeComps` | Array dimensioned by n_regions, n_years, n_sexes, n_srv_fleets specifying a multinomial weight to apply to survey age compositions, which should ideally be derived using Francis re-weighting methods |
| `Wt_SrvLenComps` | Array dimensioned by n_regions, n_years, n_sexes, n_srv_fleets specifying a multinomial weight to apply to survey length compositions, which should ideally be derived using Francis re-weighting methods |
