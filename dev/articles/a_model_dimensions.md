# Description of Model and Data Dimensions

The following tables describe all elements contained within
`input_list$data`, which is generated using the `SPoRC::Setup_x`
functions. Note that they do not detail how the `SPoRC::Setup_x`
functions should be utilized. Rather, they detail the components of
`input_list$data`. Thus, for further details on how arguments for
`SPoRC::Setup_x` functions should be defined, users should refer to the
function documentation. Furthermore, note that when `n_sexes > 1`, the
first dimension will always be females and the second dimension will
always be males. Similarly, when `n_pop > 1`, populations are indexed in
the order they are defined.

## Data Inputs for Defining Model Dimensions

| Name | Description |
|----|----|
| years | Vector specifying number of years to model |
| ages | Vector specifying number of ages to model |
| lens | Vector specifying number of lengths to model |
| n_pop | Value specifying number of populations to model |
| n_regions | Value specifying number of regions to model |
| n_sexes | Value specifying number of sexes to model |
| n_fish_fleets | Value specifying number of fishery fleets to model |
| n_srv_fleets | Value specifying number of survey fleets to model |
| n_seas | Value specifying number of seasons to model within a year |
| seasdur | Numeric vector of length n_seas specifying the duration of each season as a fraction of the year (must sum to 1) |
| spawn_seas | Integer specifying which season spawning occurs in |
| natal_region | Integer vector of length n_pop specifying the natal region for each population. Used when n_pop \> 1 to define population-specific recruitment and density dependence |
| n_proj_yrs_devs | Number of years to project deviations / random effect parameters forward |
| do_internal_comp_osa | Boolean. Whether or not internal composition OSAs are planned to be used |
| do_internal_conv_tag_osa | Boolean. Whether or not internal tagging OSAs are planned to be used |

## Data Inputs for Defining Recruitment Processes

| Name | Description |
|----|----|
| rec_lag | Value specifying the delay between spawning and when recruits enter the population. For example, if recruits enter the population as age 2, rec_lag would be specified as 2, such that the spawning biomass from year, 2 produces these recruits. A special case, rec_lag = 0 (age-0 recruitment), uses the *same* year’s own spawning biomass instead of a prior year’s; because that year’s SSB isn’t known until the spawning season is reached, recruits may only enter in the spawning season itself or later that same year |
| Use_h_prior | Value specifying whether steepness priors are used. 0: Don’t use priors, 1: Use Priors. Steepness priors are bounded between 0.2 and 1 with a scaled beta penalty |
| h_prior | Data frame specifying prior distributions for steepness parameters. Must include columns: pop (population index), region (region index), mu (mean of the prior in normal space), and sd (standard deviation of the prior in normal space). Optional lb and ub columns set the support of the scaled beta (default 0.2 and 1). For each row, a beta distribution is scaled to that interval, and the corresponding h_trans value is transformed and penalized using the log-density from the scaled beta distribution |
| do_rec_bias_ramp | Value specifying whether or not the Methot and Taylor recruitment bias ramp is conducted. 0: Don’t do bias ramp, 1: Do bias ramp |
| max_bias_ramp_fct | Value specifying the maximum bias correction factor applied to the recruitment bias ramp |
| bias_year | Vector of years specifying when to change the bias ramp. Can be specified with an NA if do_rec_bias_ramp = 0 |
| sigmaR_switch | Value specifying when to transition between using an early period sigma and late period sigma R for penalizing initial age deviations and recruitment deviations. Specify as 0 if sigmaR early is equal to sigmaR late |
| init_age_strc | Value specifying how the population age structure should be initialized. 0: Iterate to equilibrium, 1: Scalar geometric series (no movement), 2: Matrix geometric series with movement (default), 3: Matrix approach with a scalar plus group, 4: Free (no equilibrium; ln_InitDevs are the initial log numbers at age for ages 2 and older, apportioned by sex ratio) |
| equil_init_age_strc | Value specifying how initial age deviations arise. 0 (`"equil"`) == deterministic equilibrium; 1 (`"stoch_no_plus"`) == stochastic for all ages except the plus group; 2 (`"stoch_all"`) == stochastic for all ages including the plus group; 3 (`"stoch_shared_ages"`) == stochastic with user-defined age sharing via `init_age_devs_shared` |
| init_F_par | Array specifying the initialization fishing mortality for regions, seasons, and fishery fleets. Interpreted as a proportion of the mean fishing mortality (inverse-logit scale) or as an absolute rate (log scale), per `init_F_form`; `init_F_spec` sets whether it is estimated. Supersedes the legacy `init_F_prop`, which is still accepted |
| rec_model | Value specifying the recruitment model. 0 == Mean Recruitment, 1 == Beverton-Holt with steepness parameterization, 2 == Ricker in depletion form with steepness mapped through the Beverton-Holt compensation ratio (see [`vignette("c_model_equations")`](https://chengmatt.github.io/SPoRC/dev/articles/c_model_equations.md)) |
| SR_ref_yr | Year index (not a calendar year) supplying the biological inputs (WAA, maturity, natural mortality, movement) to unfished spawning biomass per recruit, and therefore to S0 and the scale of the stock-recruit curve. Default 1 (first model year); set to the terminal year index to condition the curve on terminal biologicals. Ignored under mean recruitment |
| RecDevs_pen_center | Value specifying where the recruitment deviation penalty is centered. 0 (`"fixed"`): the asserted prior mean (zero or the bias-corrected offset), constraining both level and spread; 1 (`"own_mean"`): the mean of the estimated deviations, penalizing only their spread. Cannot be combined with do_rec_bias_ramp = 1 |
| InitDevs_pen_center | Value specifying where the initial age deviation penalty is centered, with the same 0/1 coding as RecDevs_pen_center. Under `"own_mean"` with sex-specific deviations the mean pools every penalized cell across sexes |
| init_devs_pen_use | Array dimensioned by n_pop, n_regions, n_ages - 1, n_sexes of 0/1 flags naming which cells of `ln_InitDevs` the initial age penalty scores, so a curve shared across sexes (`InitDevs_sex_spec = "est_shared_s"`) is penalized once and sex-specific curves (`"est_all"`) each are |
| Use_init_sex_pen | Value specifying whether each later sex’s initial age deviations are tied to the first sex’s by a Gaussian on their difference. 0: Don’t use (default), 1: Use; requires sex-specific deviations |
| ln_sigma_init_sex | Log-scale standard deviation of that tie. A sum of squares with weight w corresponds to sigma = 1/sqrt(2w) |
| Use_rec_level_pen | Value specifying whether a penalty is applied to the log recruitment series itself, separately from the deviation penalty. 0: Don’t use (default), 1: Use |
| ln_sigma_rec_level | Log-scale standard deviation of the recruitment level penalty. A sum of squares with weight w corresponds to sigma = 1/sqrt(2w) |
| rec_level_pen_center | Value specifying where the recruitment level penalty is centered. 0 (`"fixed"`): zero; 1 (`"own_mean"`, default): the mean of the log recruitment series, so only its variability is penalized |
| rec_level_pen_yrs | Vector of length n_years with 1 marking the years the recruitment level penalty applies over and 0 elsewhere. All ones by default |
| rec_dd | Value specifying the recruitment density dependence (only used when there is a stock-recruitment relationship). 0 == local density dependence (region-specific SSB drives regional recruitment), 1 == global density dependence (summed SSB across regions drives recruitment), 999 == no density-dependent stock-recruitment form is used |
| rec_region_prop_spec | Integer specifying how recruitment regional apportionment is handled when n_pop \> 1. 0 == recruitment dispersal (recruits can be distributed across regions via estimated proportions), 1 == strict natal homing (recruits are assigned entirely to each population’s natal region) |
| t_spawn | Fraction of year in which spawning occurs |
| use_fixed_stray_rate | Integer specifying whether stray rates are supplied as a fixed external array (1) or estimated as model parameters via stray_rate_pars (0). Default 1. Only relevant when n_pop \> 1 |
| fixed_stray_rate | Array dimensioned by n_pop, n_years specifying fixed stray rate values used when use_fixed_stray_rate = 1. Values should be in \[0, 1\]. Ignored when use_fixed_stray_rate = 0 |
| stray_rate_blocks | Array dimensioned by n_pop, n_years specifying the time block index for stray rate parameters for each population and year. Unique integer values denote distinct stray rate parameter blocks |
| use_stray_rate_prior | Integer specifying whether a Beta prior is placed on estimated stray rate parameters. 0: Don’t use prior, 1: Use prior. Only relevant when use_fixed_stray_rate = 0 and n_pop \> 1. Prior contributions accumulate into rec_prop_nLL |
| stray_rate_prior | Data frame specifying Beta prior parameters for stray rates. Must include columns: pop (population index), block (time block index matching stray_rate_blocks), mu (prior mean in (0,1)), and sd (prior standard deviation). One row per population × block combination when stray_rate_spec = “est_all”; one row per block only when stray_rate_spec = “est_shared_p” |
| use_rec_region_prop_prior | Value specifying whether a Dirichlet prior is placed on recruitment regional apportionment proportions. 0: Don’t use prior, 1: Use prior |
| rec_region_prop_prior | Data frame specifying Dirichlet prior parameters for recruitment regional apportionment. Must include columns: pop (population index) and alpha (list column containing Dirichlet concentration parameter vectors of length n_regions) |
| use_fixed_rec_seas_prop | Integer specifying whether seasonal recruitment apportionment proportions are fixed (1) or estimated (0) |
| fixed_rec_seas_prop | Array dimensioned by n_pop, n_seas specifying fixed seasonal recruitment proportions when use_fixed_rec_seas_prop = 1. Ignored otherwise |
| use_rec_seas_prop_prior | Value specifying whether a Dirichlet prior is placed on estimated seasonal recruitment apportionment proportions. 0: Don’t use prior, 1: Use prior. Only relevant when use_fixed_rec_seas_prop = 0 |
| rec_seas_prop_prior | Data frame specifying Dirichlet prior parameters for seasonal recruitment apportionment. Must include columns: pop (population index) and alpha (list column containing Dirichlet concentration parameter vectors of length n_seas) |
| sexratio_blocks | Array specifying sex-ratio blocks dimensioned by n_pop, n_regions, n_years |
| use_rinit | Integer specifying whether a separate initial recruitment scalar is used to initialize the population independently of ln_global_R0. 0: Population initialized using `ln_global_R0` (default), 1: Population initialized using `ln_rinit`, with `ln_global_R0` governing only the recruitment relationship |
| init_age_devs_shared | Integer vector of length `n_ages - 1` specifying explicit age-sharing for `ln_InitDevs`. Positions with the same value share a single estimated parameter (e.g. `c(1:42, rep(42, 9))` shares the last 9 ages with age 42, giving 42 free parameters). Required when `equil_init_age_strc = 3`; `NULL` (default) uses standard behavior. |
| use_r0_prior | Integer specifying whether a lognormal prior is placed on R0. 0: Don’t use prior (default), 1: Use prior |
| r0_prior | Data frame specifying lognormal prior parameters for R0. Must include columns: pop (population index), mu (prior mean on the natural scale), and sd (prior standard deviation on the log scale) |
| map_ln_RecDevs | Array dimensioned by n_pop, n_regions, n_years mirroring the ln_RecDevs factor map: an estimation index where a deviation is estimated, NA where it is fixed (mapped off). Built by Setup_Mod_Rec and refreshed by fit_model from the map handed to RTMB::MakeADFun. The recruitment deviation penalty is evaluated only where this is not NA |

## Data Inputs for Defining Biological Processes

| Name | Description |
|----|----|
| WAA | Weight-at-age values dimensioned by n_pop, n_regions, n_years, n_seas, n_ages, n_sexes |
| WAA_fish | Fishery weight-at-age values dimensioned by n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets |
| WAA_srv | Survey weight-at-age values dimensioned by n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets |
| MatAA | Maturity-at-age values dimensioned by n_pop, n_regions, n_years, n_seas, n_ages, n_sexes |
| AgeingError | Ageing error matrix dimensioned by n_years, n_modelled_ages, n_observed_ages, where rows across modeled ages sum to 1 |
| fit_lengths | Value describing whether or not to fit length composition data. 0: Don’t fit lengths, 1: Fit lengths |
| SizeAgeTrans | Size-age transition matrix dimensioned by n_pop, n_regions, n_years, n_seas, n_lens, n_ages, n_sexes. Can be specified as NA if length compositions are not fit |
| LenBinMap | Matrix dimensioned by n_lens, n_observed_lens mapping the model’s length bins onto the bins the compositions are recorded on, each row summing to 1. NULL when the two coincide; when supplied, the observed length composition arrays are dimensioned by n_observed_lens rather than n_lens, and expected compositions are mapped through it inside the likelihood the way AgeingError maps ages |
| ln_growth_pars | Parameter array dimensioned by n_pop, n_regions, n_sexes, n_growth_pars of log growth parameters in the order L1, L2, K, CV1, CV2, and rho under the Richards form |
| ln_growth_devs | Parameter array dimensioned by n_pop, n_regions, n_years, n_growth_pars, n_sexes of deviations on the growth parameters, zero and mapped off for parameters that do not vary and for years outside growth_tv_years |
| map_ln_growth_devs | Array dimensioned like ln_growth_devs indicating which growth parameter deviations are fixed (mapped off), read by the process error penalty |
| ln_growth_semipar_devs | Parameter array dimensioned by n_pop, n_regions, n_years, n_ages, n_sexes of semi-parametric deviations on mean length at age, multiplying the parametric curve |
| growth_pe_pars | Parameter array dimensioned by n_pop, n_regions, max(4, n_ages, n_growth_pars), n_sexes, 2 of process error hyperparameters for both growth deviation streams. The first stream holds one log process error standard deviation per growth parameter for the time-varying deviations. The second holds the semi-parametric surface’s correlations by age, year and cohort in slots one to three and the log scale in slot four for the correlated forms; the iid and random walk forms instead read one log sigma per age. Slots a form does not read are mapped off |
| map_ln_growth_semipar_devs | Array dimensioned like ln_growth_semipar_devs indicating which semi-parametric growth deviations are fixed (mapped off) |
| growth_semipar_bins | Integer vector of the age indices the semi-parametric surface is evaluated over, used to subset the age dimension when evaluating GMRF or 2D AR(1) likelihoods. Defaults to every age |
| addtocomp | Constant to add to all composition data |
| comp_const_obs | Value specifying whether addtocomp is also added to the observed proportions used as multinomial weights. 1 (default): added to both observed and expected, the long-standing behavior; 0: added only inside the logarithms, a convention several existing assessments use |
| addtofishidx | Constant to add to all fishery index data |
| addtosrvidx | Constant to add to all survey index data |
| addtotag | Constant to add to all tagging data |
| Use_M_prior | Value specifying whether natural mortality priors are used. 0: Don’t use prior, 1: Use prior |
| M_prior | A data frame specifying how natural mortality priors should be applied. Must include columns: popblk (population block index), regionblk (region block index), yearblk (year block index), ageblk (age block index), sexblk (sex block index), mu (prior mean on natural scale), and sd (prior standard deviation on log scale) |
| Fixed_natmort | Natural mortality array dimensioned by n_pop, n_regions, n_years, n_ages, n_sexes |
| M_blocks | Array that specifies the natural mortality blocking structure, dimensioned by n_pop, n_regions, n_years, n_ages, n_sexes. Unique integer values denote distinct natural mortality parameter blocks |

## Data Inputs for Defining Spawning and Multi-Population Processes

| Name | Description |
|----|----|
| sgl_seas_spawning_movement | Movement matrix applied during the spawning event when n_seas == 1 and n_pop \> 1, representing fish returning to natal areas to spawn. Dimensioned by n_pop, n_regions, n_regions, n_years, n_ages, n_sexes. Ignored when n_seas \> 1 or n_pop == 1 |

## Data Inputs for Defining Movement Processes

| Name | Description |
|----|----|
| do_recruits_move | Value specifying whether recruits are allowed to move. 0: Recruits don’t move, 1: Recruits move |
| use_fixed_movement | Value specifying whether or not to use a fixed movement matrix. 0: Don’t use fixed movement, 1: Use fixed movement |
| Fixed_Movement | Fixed movement matrix dimensioned by n_pop, n_regions, n_regions, n_years, n_seas, n_ages, n_sexes |
| Use_Movement_Prior | Value specifying whether or not to use movement priors. 0: Don’t use movement prior, 1: Use movement prior |
| Movement_prior | Data frame with columns representing population (pop), origin region (region_from), year (year), season (seas), representative age (age), sex (sex), and a list column (alpha) containing Dirichlet prior parameters. Each element of alpha is a numeric vector of length n_regions, specifying the prior concentration values for movement among destination regions |
| cont_vary_movement | Integer indicating whether movement is continuously varying across regions, years, and ages. 0 = none, 1 = iid deviations |
| map_move_devs | Array dimensioned by n_regions, n_regions - 1, n_years, n_ages, and n_sexes specifying which movement parameters are shared and mapped off |
| move_type | Integer indicating movement model type. 0 = unstructured Markov (multinomial logit), 1 = Continuous Time Markov Chain (CTMC) |
| ctmc_move_dat | Data frame with CTMC covariates (regions, years, ages, sexes, plus formula variables) used to build design matrices for diffusion and preference when move_type == 1. Can include projection years with projected covariate values |
| diffusion_formula | R formula specifying diffusion covariates for CTMC movement (when move_type == 1) |
| preference_formula | R formula specifying preference covariates for CTMC movement (when move_type == 1) |
| adjacency_mat | Square adjacency matrix (n_regions × n_regions) defining allowed transitions between regions for CTMC movement |
| adjacency_collapsed | Collapsed adjacency matrix (n_regions × n_regions - 1) defining allowed transitions between regions for CTMC movement, excluding diagonals |
| area_r | Numeric vector of region areas (length n_regions) used to scale diffusion rates in CTMC movement |
| ctmc_diffusion_bounds | How the generator is kept Metzler: `"none"` (0), `"softplus"` (1), or `"upwind"` (2), the discontinuous Galerkin flux, which carries diffusion whole and reads no eps |
| ctmc_diffusion_eps | Positive numeric softplus width used by `ctmc_diffusion_bounds = "softplus"` (default 0.1); the flow on an edge where taxis exactly cancels diffusion is eps \* log(2), so it is a floor on exchange and not only a smoothing width |

## Data Inputs for Defining Tagging Processes

| Name | Description |
|----|----|
| use_conv_fish_tagging | Integer vector of length n_fish_fleets indicating whether conventional fishery tagging data are used for each fleet. 0: Don’t use tagging data, 1: Use tagging data |
| conv_tag_release_indicator | Matrix dimensioned by n_conv_tag_cohorts, 3. Columns represent the tag release region (column 1), tag release year (column 2), and tag release season (column 3) for each cohort |
| n_conv_tag_cohorts | Value specifying the number of conventional tag cohorts released |
| conv_tag_max_liberty | Value specifying the maximum years at liberty to track tag cohorts |
| conv_tagged_fish | Array specifying the number of tagged fish for a given cohort, dimensioned by n_conv_tag_cohorts, n_pop, n_ages, n_sexes |
| obs_conv_tag_fish_recap | Array of observed conventional tag recaptures dimensioned by conv_tag_max_liberty, n_seas, n_conv_tag_cohorts, n_pop, n_regions, n_ages, n_sexes, n_fish_fleets. If no age or sex information is available, the array can be initialized with 0s and recaptured individuals input into the first age and first sex dimension, coupled with conv_tag_age_pool and conv_tag_sex_pool to sum across those dimensions when computing the likelihood |
| conv_fish_tag_like | Integer specifying the tag likelihood to use. 0: Poisson, 1: Negative Binomial, 2: Multinomial release-conditioned, 3: Multinomial recapture-conditioned, 4: Dirichlet-Multinomial release-conditioned, 5: Dirichlet-Multinomial recapture-conditioned |
| conv_tag_mixing_period | Value specifying the mixing period for tag cohorts in seasonal (or annual if model is annual) time steps. Tag observations within the mixing period are excluded from the likelihood |
| conv_tag_t_tagging | Fraction of the season in which tagging occurs. Specify as 1 if tagging happens at the start of the season such that mortality is not discounted |
| use_conv_tag_fishrep_prior | Value specifying whether or not a tag reporting rate prior should be used. 0: Don’t use prior, 1: Use prior |
| conv_tag_fishrep_prior | Data frame containing prior specifications for conventional fishery tag reporting parameters. Must include columns: region (region index), block (time block index), fleet (fleet index), mu (numeric mean for reporting rate prior in normal space; use NA if symmetric beta is used), sd (prior standard deviation in normal space), and type (0 = symmetric beta, 1 = regular beta) |
| conv_tag_pop_pool | List specifying how tag recaptures along the population axis should be pooled when computing the likelihood. Each element of the list contains the population indices to sum across for a given pooled group |
| conv_tag_age_pool | List specifying how tag recaptures along the age axis should be pooled when computing the likelihood. Each element of the list contains the age indices to sum across for a given pooled group. For example, if no age information is available: list(c(1:n_ages)) sums all ages together |
| conv_tag_sex_pool | List specifying how tag recaptures along the sex axis should be pooled when computing the likelihood. Each element of the list contains the sex indices to sum across for a given pooled group. For example, if no sex information is available: list(c(1:n_sexes)) sums both sexes together |
| conv_tag_fish_reporting_blocks | Array dimensioned by n_regions, n_years, n_fish_fleets indicating the time block index for tag reporting rate parameters in each region and fleet |
| conv_fish_tag_attr | Data object specifying how tagged fish attributes (e.g., age, sex, population) are apportioned at release when full individual-level data are unavailable. Passed to the internal release_conv_tag_attr function |
| conv_tag_release_platform | Matrix dimensioned by n_conv_tag_cohorts specifying the release platform (e.g., population, fishery or survey fleet) used to attribute selectivity and effort to tagged fish at the time of release |

## Data Inputs for Defining Catch and Fishing Mortality Processes

| Name | Description |
|----|----|
| ObsCatch | Observed catch data dimensioned by n_regions, n_years, n_seas, n_fish_fleets. Used in region-aggregated catch likelihoods |
| ObsCatch_pop | Population-specific observed catch data dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets. Used in population-specific catch likelihoods |
| catch_units | Array dimensioned by n_fish_fleets describing catch units. 0 == Abundance, 1 == Biomass (default) |
| UseCatch | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets describing whether to fit to catch data in a given year, season, and fleet. 0: Don’t fit catch, 1: Fit catch |
| UseCatch_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets describing whether to fit to population-specific catch data. 0: Don’t fit catch, 1: Fit catch |
| Use_F_pen | Value specifying whether to use a fishing mortality penalty to regularize fishing mortality deviations. 0: Don’t use regularity penalty, 1: Use regularity penalty. Which cells are penalized is read off map_ln_F_devs |
| Fdev_pen_center | Value specifying where the iid fishing mortality deviation penalty is centered. 0 (`"fixed"`): zero, constraining both level and spread; 1 (`"own_mean"`): the mean of the estimated deviations, penalizing only their spread |
| ObsDiscard | Observed discard data dimensioned by n_regions, n_years, n_seas, n_fish_fleets. Used in region-aggregated discard likelihoods |
| discard_units | Array dimensioned by n_fish_fleets describing discard units. 0 == Abundance, 1 == Biomass |
| UseDiscard | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets describing whether to fit to discard data in a given year, season, and fleet. 0: Don’t fit discard, 1: Fit discard |
| UseDiscard_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets describing whether to fit to population-specific discard data. 0: Don’t fit discard, 1: Fit discard |
| Use_dmr_pen | Value specifying whether to use a discard mortality rate penalty. 0: Don’t use penalty, 1: Use penalty. Which cells are penalized is read off map_logit_dmr_devs |
| map_ln_F_devs | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets mirroring the ln_F_devs factor map: an estimation index where a deviation is estimated, NA where it is fixed (mapped off). Built by Setup_Mod_Catch_and_F and refreshed by fit_model from the map handed to RTMB::MakeADFun. The fishing mortality deviation penalty is evaluated only where this is not NA |
| map_logit_dmr_devs | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets mirroring the logit_dmr_devs factor map, in the same form as map_ln_F_devs. The discard mortality rate deviation penalty is evaluated only where this is not NA |

## Data Inputs for Defining Fishery Indices and Compositions

Region-aggregated observations (summed across populations) and
population-specific observations are handled separately.
Region-aggregated inputs use dimensions beginning with n_regions;
population-specific inputs use dimensions beginning with n_pop,
n_regions.

| Name | Description |
|----|----|
| ObsFishIdx | Fishery index dimensioned by n_regions, n_years, n_seas, n_fish_fleets |
| ObsFishIdx_SE | Fishery index standard errors dimensioned by n_regions, n_years, n_seas, n_fish_fleets |
| UseFishIdx | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets describing whether to fit to the fishery index in a given year, season, and fleet. 0: Don’t fit fishery index, 1: Fit fishery index |
| fish_idx_type | Matrix dimensioned by n_fish_fleets specifying the index type for a given fishery fleet. 0: Abundance index, 1: Biomass index (uses WAA_fish for calculations), 999: None Available |
| t_fish | Array dimensioned by n_regions, n_seas, n_fish_fleets specifying fishery index timing as a proportion of the season. Numbers at age are decayed by exp(-t_fish \* ZAA) before the index is formed, the same convention t_srv uses for surveys. Defaults to 0 (start of season); use 0.5 for a mid-season index. Set via t_fish in Setup_Mod_FishIdx_and_Comps, and via the matching argument to Setup_Sim_Fishing on the simulation side |
| fish_idx_ages | Array dimensioned by n_ages, n_fish_fleets of 0/1 weights selecting which ages contribute to each fleet’s index total. All ones by default. Applies to the index sum only; selectivity, catch, and compositions are unaffected |
| FishIdx_LikeType | Vector dimensioned by n_fish_fleets specifying the index error structure. 0: Lognormal (default; standard errors on the log scale), 1: Normal on the arithmetic scale, 2: Multivariate normal on the arithmetic scale with a fixed covariance from FishIdx_Cov. One-step-ahead residuals are only available for lognormal fleets |
| FishIdx_Cov | List with one element per fishery fleet holding the fixed covariance matrix for fleets using the multivariate normal, NULL otherwise. Each matrix is square with one row per fitted observation, ordered as the observations appear when scanning that fleet’s UseFishIdx slice in array order. Validated at setup for symmetry and positive definiteness |
| ObsCatchAA, ObsDiscardAA | Age-disaggregated retained catch and discards, each dimensioned by n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets. A fleet fits these or the aggregated form with compositions, never both. The sex margin is required whatever the fleet reports: a stream summed over sexes carries its observation in sex slot one, and an array a dimension short is refused |
| UseCatchAA, UseDiscardAA | Arrays shaped like their observations, 1 where an observation is fit. A cell not fit for catch at age is also not fished, the way UseCatch governs closures for the aggregated stream |
| ObsCatchAA_SE, ObsDiscardAA_SE | Reported standard errors shaped like their observations, read only when the stream’s sigma_form asks for them |
| CatchAA_Type, DiscardAA_Type | Matrices dimensioned by n_years, n_fish_fleets naming the margins each fleet reports separately, which may change between years. 0: summed over regions and sexes, 1: split by region and summed over sexes (default), 2: summed over regions and split by sex, 3: split by both |
| CatchAA_LikeType, DiscardAA_LikeType | Vectors dimensioned by n_fish_fleets. 0: Lognormal (default), 1: Normal on the arithmetic scale |
| CatchAA_sigma_form, DiscardAA_sigma_form | Vectors dimensioned by n_fish_fleets naming where the observation error comes from. 0: the estimated parameter alone (default), 1: the reported standard errors alone, 2: additive, 3: in quadrature |
| AgeObsCorr_catch, AgeObsCorr_discard | Vectors dimensioned by n_fish_fleets giving the correlation across ages within a cell. 0: iid, 1: AR(1) in age distance, 2: unstructured, 3: separable AR(1) over ages and years. Each has a \_pop counterpart, since the population-specific streams carry their own |
| ObsCatchAA_pop, UseCatchAA_pop, and the rest | Population-specific counterparts of every row above, with a leading n_pop dimension on the arrays |
| FishAgeComps_bins | Array dimensioned by n_ages, n_fish_fleets of 0/1 weights selecting which observed bins (after ageing error) each fleet’s age compositions are fit over. All ones by default. Observed and expected compositions are subset and renormalized within the selected bins; bins outside are left out of the likelihood |
| ObsFishIdx_pop | Population-specific fishery index dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets |
| ObsFishIdx_pop_SE | Population-specific fishery index standard errors dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets |
| UseFishIdx_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets describing whether to fit to population-specific fishery indices. 0: Don’t fit, 1: Fit |
| ObsFishAgeComps | Observed fishery age compositions dimensioned by n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets. Can be input as proportions or numbers, as these are normalized within the model |
| UseFishAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets specifying whether or not to fit to fishery age compositions. 0: Don’t fit, 1: Fit |
| ISS_FishAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood |
| Wt_FishAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight to apply to fishery age compositions, ideally derived using Francis re-weighting |
| ObsFishAgeComps_pop | Population-specific observed fishery age compositions dimensioned by n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets. Can be input as proportions or numbers |
| UseFishAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets specifying whether to fit to population-specific fishery age compositions. 0: Don’t fit, 1: Fit |
| ISS_FishAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying the input sample size for population-specific fishery age composition likelihoods |
| Wt_FishAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight to apply to population-specific fishery age compositions |
| ObsFishLenComps | Observed fishery length compositions dimensioned by n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets. Can be input as proportions or numbers, as these are normalized within the model |
| UseFishLenComps | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets specifying whether or not to fit to fishery length compositions. 0: Don’t fit, 1: Fit |
| ISS_FishLenComps | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood |
| Wt_FishLenComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight to apply to fishery length compositions, ideally derived using Francis re-weighting |
| ObsFishLenComps_pop | Population-specific observed fishery length compositions dimensioned by n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets. Can be input as proportions or numbers |
| UseFishLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets specifying whether to fit to population-specific fishery length compositions. 0: Don’t fit, 1: Fit |
| ObsFish_caal | Observed fishery conditional age-at-length dimensioned by n_regions, n_years, n_seas, n_lens, n_ages, n_sexes, n_fish_fleets: the ages of the fish aged from each length bin. Can be input as proportions or numbers |
| UseFish_caal | Array dimensioned by n_regions, n_years, n_seas, n_lens, n_fish_fleets specifying whether to fit each length bin’s age composition. 0: Don’t fit, 1: Fit |
| ISS_Fish_caal | Array dimensioned by n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets specifying the input sample size of each length bin, the number of fish aged from it |
| Wt_Fish_caal | Array dimensioned by n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets specifying a multinomial weight to apply to fishery conditional age-at-length, set in [`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md) |
| ObsSrv_caal, UseSrv_caal, ISS_Srv_caal, Wt_Srv_caal | The survey counterparts, with n_srv_fleets as the last dimension |
| ISS_FishLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets specifying the input sample size for population-specific fishery length composition likelihoods |
| Wt_FishLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight to apply to population-specific fishery length compositions |
| FishAgeComps_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for region-aggregated fishery age compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| FishAgeComps_pop_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for population-specific fishery age compositions. Same options as FishAgeComps_LikeType |
| FishLenComps_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for region-aggregated fishery length compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| FishLenComps_pop_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for population-specific fishery length compositions. Same options as FishLenComps_LikeType |
| FishAgeComps_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how region-aggregated age composition data should be structured. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| FishAgeComps_pop_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how population-specific age composition data should be structured. Same options as FishAgeComps_Type |
| FishLenComps_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how region-aggregated length composition data should be structured. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| FishLenComps_pop_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how population-specific length composition data should be structured. Same options as FishLenComps_Type |
| ObsFishAgeComps_discard | Observed discard age compositions dimensioned by n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets. Can be input as proportions or numbers, as these are normalized within the model |
| UseFishAgeComps_discard | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets specifying whether or not to fit to discard age compositions. 0: Don’t fit, 1: Fit |
| ISS_FishAgeComps_discard | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial discard age composition likelihood |
| ObsFishLenComps_discard | Observed discard length compositions dimensioned by n_regions, n_years, n_seas, n_lens, n_sexes, n_fish_fleets. Can be input as proportions or numbers, as these are normalized within the model |
| UseFishLenComps_discard | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets specifying whether or not to fit to discard length compositions. 0: Don’t fit, 1: Fit |
| ISS_FishLenComps_discard | Array dimensioned by n_regions, n_years, n_seas, n_fish_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial discard length composition likelihood |
| ObsFishAgeComps_discard_pop | Population-specific observed discard age compositions dimensioned by n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets. Can be input as proportions or numbers |
| UseFishAgeComps_discard_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets specifying whether to fit to population-specific discard age compositions. 0: Don’t fit, 1: Fit |
| ISS_FishAgeComps_discard_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying the input sample size for population-specific discard age composition likelihoods |
| UseFishLenComps_discard_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets specifying whether to fit to population-specific discard length compositions. 0: Don’t fit, 1: Fit |
| ISS_FishLenComps_discard_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets specifying the input sample size for population-specific discard length composition likelihoods |
| FishAgeComps_discard_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for region-aggregated discard age compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available |
| FishLenComps_discard_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for region-aggregated discard length compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available |
| FishAgeComps_discard_pop_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for population-specific discard age compositions. Same options as FishAgeComps_discard_LikeType |
| FishLenComps_discard_pop_LikeType | Vector dimensioned by n_fish_fleets specifying the likelihood for population-specific discard length compositions. Same options as FishLenComps_discard_LikeType |
| FishAgeComps_discard_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how region-aggregated discard age composition data should be structured. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available |
| FishLenComps_discard_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how region-aggregated discard length composition data should be structured. Same options as FishAgeComps_discard_Type |
| FishAgeComps_discard_pop_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how population-specific discard age composition data should be structured. Same options as FishAgeComps_discard_Type |
| FishLenComps_discard_pop_Type | Matrix dimensioned by n_years, n_fish_fleets specifying how population-specific discard length composition data should be structured. Same options as FishAgeComps_discard_Type |

## Data Inputs for Defining Survey Indices and Compositions

Region-aggregated observations (summed across populations) and
population-specific observations are handled separately.
Region-aggregated inputs use dimensions beginning with n_regions;
population-specific inputs use dimensions beginning with n_pop,
n_regions.

| Name | Description |
|----|----|
| ObsSrvIdx | Survey index dimensioned by n_regions, n_years, n_seas, n_srv_fleets |
| ObsSrvIdx_SE | Survey index standard errors dimensioned by n_regions, n_years, n_seas, n_srv_fleets |
| UseSrvIdx | Array dimensioned by n_regions, n_years, n_seas, n_srv_fleets describing whether to fit to the survey index in a given year, season, and fleet. 0: Don’t fit, 1: Fit |
| srv_idx_type | Matrix dimensioned by n_srv_fleets specifying the index type for a given survey fleet. 0: Abundance index, 1: Biomass index (uses WAA_srv for calculations), 999: None Available |
| srv_idx_ages | Array dimensioned by n_ages, n_srv_fleets of 0/1 weights selecting which ages contribute to each fleet’s index total. All ones by default. Restricting a fleet to a single age turns it into an index of that age alone; the fleet’s compositions keep using the full age range |
| SrvIdx_LikeType | Vector dimensioned by n_srv_fleets specifying the index error structure. 0: Lognormal (default), 1: Normal on the arithmetic scale, 2: Multivariate normal with a fixed covariance from SrvIdx_Cov. One-step-ahead residuals are only available for lognormal fleets |
| SrvIdx_Cov | List with one element per survey fleet holding the fixed covariance matrix for fleets using the multivariate normal, NULL otherwise. Same conventions as FishIdx_Cov |
| ObsSrvIdxAA | Age-disaggregated survey index dimensioned by n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets. A fleet fits this or the aggregated index, never both. Same shape rules as ObsCatchAA |
| UseSrvIdxAA | Array shaped like ObsSrvIdxAA, 1 where an observation is fit |
| ObsSrvIdxAA_SE | Reported standard errors shaped like ObsSrvIdxAA, read only when SrvIdxAA_sigma_form asks for them |
| SrvIdxAA_Type | Matrix dimensioned by n_years, n_srv_fleets, following the same codes as its fishery counterpart |
| SrvIdxAA_LikeType, SrvIdxAA_sigma_form | Vectors dimensioned by n_srv_fleets, following the same codes as their fishery counterparts |
| AgeObsCorr_srv_idx | Vector dimensioned by n_srv_fleets giving the correlation across ages, with the same codes as AgeObsCorr_catch and a \_pop counterpart |
| ObsSrvIdxAA_pop, UseSrvIdxAA_pop, ObsSrvIdxAA_pop_SE | Population-specific counterparts, with a leading n_pop dimension |
| SrvAgeComps_bins | Array dimensioned by n_ages, n_srv_fleets of 0/1 weights selecting which observed bins each fleet’s age compositions are fit over. Same conventions as FishAgeComps_bins |
| ObsSrvIdx_pop | Population-specific survey index dimensioned by n_pop, n_regions, n_years, n_seas, n_srv_fleets |
| ObsSrvIdx_pop_SE | Population-specific survey index standard errors dimensioned by n_pop, n_regions, n_years, n_seas, n_srv_fleets |
| UseSrvIdx_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_srv_fleets describing whether to fit to population-specific survey indices. 0: Don’t fit, 1: Fit |
| ObsSrvAgeComps | Observed survey age compositions dimensioned by n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets. Can be input as proportions or numbers, as these are normalized within the model |
| UseSrvAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_srv_fleets specifying whether or not to fit to survey age compositions. 0: Don’t fit, 1: Fit |
| ISS_SrvAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood |
| Wt_SrvAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight to apply to survey age compositions, ideally derived using Francis re-weighting |
| ObsSrvAgeComps_pop | Population-specific observed survey age compositions dimensioned by n_pop, n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets. Can be input as proportions or numbers |
| UseSrvAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_srv_fleets specifying whether to fit to population-specific survey age compositions. 0: Don’t fit, 1: Fit |
| ISS_SrvAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying the input sample size for population-specific survey age composition likelihoods |
| Wt_SrvAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight to apply to population-specific survey age compositions |
| ObsSrvLenComps | Observed survey length compositions dimensioned by n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets. Can be input as proportions or numbers, as these are normalized within the model |
| UseSrvLenComps | Array dimensioned by n_regions, n_years, n_seas, n_srv_fleets specifying whether or not to fit to survey length compositions. 0: Don’t fit, 1: Fit |
| ISS_SrvLenComps | Array dimensioned by n_regions, n_years, n_seas, n_srv_fleets specifying the input sample size for a multinomial or Dirichlet-multinomial likelihood |
| Wt_SrvLenComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight to apply to survey length compositions, ideally derived using Francis re-weighting |
| ObsSrvLenComps_pop | Population-specific observed survey length compositions dimensioned by n_pop, n_regions, n_years, n_seas, n_lens, n_sexes, n_srv_fleets. Can be input as proportions or numbers |
| UseSrvLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_srv_fleets specifying whether to fit to population-specific survey length compositions. 0: Don’t fit, 1: Fit |
| ISS_SrvLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_srv_fleets specifying the input sample size for population-specific survey length composition likelihoods |
| Wt_SrvLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight to apply to population-specific survey length compositions |
| SrvAgeComps_LikeType | Vector dimensioned by n_srv_fleets specifying the likelihood for region-aggregated survey age compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| SrvAgeComps_pop_LikeType | Vector dimensioned by n_srv_fleets specifying the likelihood for population-specific survey age compositions. Same options as SrvAgeComps_LikeType |
| SrvLenComps_LikeType | Vector dimensioned by n_srv_fleets specifying the likelihood for region-aggregated survey length compositions. 0: Multinomial, 1: Dirichlet-multinomial, 2: Logistic-normal with iid covariance, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| SrvLenComps_pop_LikeType | Vector dimensioned by n_srv_fleets specifying the likelihood for population-specific survey length compositions. Same options as SrvLenComps_LikeType |
| SrvAgeComps_Type | Matrix dimensioned by n_years, n_srv_fleets specifying how region-aggregated age composition data should be structured. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| SrvAgeComps_pop_Type | Matrix dimensioned by n_years, n_srv_fleets specifying how population-specific survey age composition data should be structured. Same options as SrvAgeComps_Type |
| SrvLenComps_Type | Matrix dimensioned by n_years, n_srv_fleets specifying how region-aggregated length composition data should be structured. 0: Aggregated across sexes and regions, 1: Split by sexes and regions, 2: Joint by sex but split by region, 999: None Available. Further options in Get_Comp_Likelihoods.R |
| SrvLenComps_pop_Type | Matrix dimensioned by n_years, n_srv_fleets specifying how population-specific survey length composition data should be structured. Same options as SrvLenComps_Type |
| t_srv | Array dimensioned by n_regions, n_seas, n_srv_fleets specifying survey timing as a proportion of the season |
| do_srv_q_cov | Numeric value indicating whether a survey catchability covariate is used. 0 == Not used, 1 == Used |
| srv_q_cov | Array dimensioned by n_regions, n_years, n_srv_fleets, n_covariates representing the covariates used to compute survey catchability |

## Data Inputs for Defining Fishery Selectivity and Catchability

| Name | Description |
|----|----|
| cont_tv_fish_sel | Matrix dimensioned by n_regions, n_fish_fleets specifying whether and how continuous time-varying selectivity is applied. 0: None, 1: iid deviations, 2: random walk, 3: 3D GMRF with marginal variance (semi-parametric), 4: 3D GMRF with conditional variance (semi-parametric). Further details in model_priors_penalties.R, model_precision.R, and model_selectivity.R |
| fish_sel_blocks | Array dimensioned by n_regions, n_years, n_fish_fleets specifying selectivity time blocks. Unique integers denote distinct selectivity parameter blocks |
| fish_sel_model | Array dimensioned by n_regions, n_years, n_fish_fleets specifying the selectivity functional form. 0: Logistic (a50 and slope), 1: Gamma dome-shaped, 2: Power function, 3: Logistic (a50 and a95), 4: Double Normal (6 parameters), 5: Non-parametric (logit scale), 6/7: Logistic with asymptote, 8: Bicubic spline, 9: Non-parametric on the log scale, standardized within each year. Further details in model_selectivity.R |
| fish_sel_bicubic_binnodes, fish_sel_bicubic_yrnodes | Arrays dimensioned by n_regions, n_years, n_fish_fleets giving the number of bin/year spline nodes where fish_sel_model == 8; 0 elsewhere |
| fish_sel_bicubic_selstyr, fish_sel_bicubic_nselbins | Arrays dimensioned by n_regions, n_years, n_fish_fleets giving the optional fitted-region restrictions (start year, number of bins); 0 means unrestricted. The start year applies to bicubic blocks only; the number of bins (the `_NSelBins_<n>` suffix) applies to every form, holding bins beyond n at bin n’s value |
| fish_sel_bicubic_Wbin, fish_sel_bicubic_Wyr | Precomputed natural-cubic-spline interpolation weight matrices mapping bin-node/year-node values onto every bin/year, for bicubic blocks |
| fish_selex_type | Integer specifying whether fishery selectivity is age-based (0) or length-based (1) |
| use_fixed_fish_sel | Integer specifying whether fishery selectivity is fixed externally (1) or estimated (0) |
| fish_q_blocks | Array dimensioned by n_regions, n_years, n_fish_fleets specifying catchability time blocks. Unique integers denote distinct catchability parameter blocks |
| fish_q_type | Vector dimensioned by n_fish_fleets specifying how fishery catchability is obtained. Same codes and behavior as srv_q_type. 0: Estimated as exp(ln_fish_q) (default), 1: Solved analytically as the ratio of mean observed to mean predicted index, 2: Solved analytically on the log scale as exp(mean(log(obs) - log(pred))). Analytic fleets have their ln_fish_q fixed automatically, ignore block structure, and cannot carry catchability covariates or priors |
| Use_fish_q_prior | Fishery catchability prior indicator. 0 == don’t use, 1 == use |
| fish_q_prior | Data frame containing prior specifications for fishery catchability parameters. Must include columns: region, fleet, block, mu (prior mean on natural scale), and sd (prior standard deviation on log scale). Each row specifies a log-normal prior for one catchability parameter |
| map_ln_fishsel_devs | Array dimensioned by n_regions, n_years, n_ages, n_sexes, n_fish_fleets indicating which continuous time-varying selectivity deviation values are fixed (mapped off) |
| Use_fish_selex_prior | Integer (0 or 1). Flag to enable log-normal priors on fishery selectivity parameters as specified in fish_selex_prior |
| fish_selex_prior | Data frame containing prior specifications for fishery selectivity parameters. Must include columns: region, fleet, block, sex, par (parameter index), mu (prior mean on natural scale), and sd (prior standard deviation on log scale). Each row specifies a log-normal prior for one selectivity parameter |
| fishsel_devs_min_shared_bins | Integer vector specifying the reference (minimum) bin index within each shared deviation group, used to subset the bin dimension when evaluating GMRF or 2D AR(1) likelihoods (PE models 3-5). Defaults to 1:n_ages when no bin sharing is specified |
| fish_sel_pen_wts | List with one element per fishery fleet, each a named list of smoothness penalty weights (smooth_bin_curve, smooth_bin_diff, smooth_yr_diff, smooth_yr_curve, smooth_dome, smooth_mean_center; unset names default to 0) evaluated on the realized fishery selectivity surface, plus optional bin_range (bins the penalties act over), normalize (whether weights are divided by the number of penalized bins/years), and yr_diff_ref (reference log-selectivity anchoring the first penalized year of smooth_yr_diff). Each weight may be a scalar or a per-year vector. Built by Setup_Mod_Weighting from a single shared specification or an unnamed per-fleet list |
| fishsel_pe_wt | Vector dimensioned by n_fish_fleets multiplying each fleet’s selectivity process error likelihood. Default 1; 0 removes the distributional penalty while the deviations remain estimated |
| fishsel_rw_init_sigma | Vector dimensioned by n_fish_fleets giving the standard deviation on the first year of a random walk deviation series. Default 5 (first year effectively free); NA starts the walk at zero under the walk’s own estimated sigma |
| fish_sel_bin_dev_bins | Array dimensioned by n_ages (or n_lens), n_fish_fleets of 0/1 flags marking the bins whose selectivity is overridden by exp(ln_fishsel_bin_devs) rather than taken from the functional form. All zeros by default |
| cont_tv_fishsel_bin_devs | Vector dimensioned by n_fish_fleets specifying process error on the bin-override deviations. 0: None, 1: iid, 2: random walk |
| fishsel_sex_par_offset | Vector dimensioned by n_fish_fleets of 0/1 flags. 1 reads the stored fixed-effect slots of every sex beyond the first as additive offsets on the first sex’s transformed parameters (`fish_sel_sex_offset = "par"` or `"par_scale"`) |
| fishsel_sex_scale_offset | Vector dimensioned by n_fish_fleets of 0/1 flags. 1 multiplies the realized curve of every sex beyond the first by `exp(ln_fishsel_sex_scale)` (`fish_sel_sex_offset = "scale"` or `"par_scale"`) |
| fishsel_bin_devs_rw_init_sigma | Vector dimensioned by n_fish_fleets giving the first-year standard deviation for random walk bin-override deviations, with the same conventions as fishsel_rw_init_sigma |
| map_ln_fishsel_bin_devs | Array mirroring the ln_fishsel_bin_devs factor map: an estimation index where a deviation is estimated, NA where fixed. Only overridden bins carry estimated deviations |
| Use_fish_selex_penalty | Integer (0 or 1). Flag to enable the centering penalty on sets of fishery selectivity fixed-effect parameters specified in fish_selex_penalty |
| fish_selex_penalty | Data frame with columns region, fleet, block, sex, par (a single index or a list column of integer vectors naming a set), and wt. Each row penalizes wt \* (log(mean(exp(pars))))^2, pushing the set’s average selectivity toward one. Intended for log-scale parameter sets such as the log-scale non-parametric form |

## Data Inputs for Defining Retention Selectivity

| Name | Description |
|----|----|
| cont_tv_ret_sel | Matrix dimensioned by n_regions, n_fish_fleets specifying whether and how continuous time-varying retention selectivity is applied. Same options as cont_tv_fish_sel |
| ret_sel_blocks | Array dimensioned by n_regions, n_years, n_fish_fleets specifying retention selectivity time blocks. Unique integers denote distinct retention selectivity parameter blocks |
| ret_sel_model | Array dimensioned by n_regions, n_years, n_fish_fleets specifying the retention selectivity functional form. Same options as fish_sel_model, including bicubic (8) |
| ret_sel_bicubic_binnodes, ret_sel_bicubic_yrnodes, ret_sel_bicubic_selstyr, ret_sel_bicubic_nselbins, ret_sel_bicubic_Wbin, ret_sel_bicubic_Wyr | Same meaning as their fish_sel_bicubic\_\* counterparts, for retention |
| ret_selex_type | Integer specifying whether retention selectivity is age-based (0) or length-based (1) |
| use_fixed_ret_sel | Integer specifying whether retention selectivity is fixed externally (1) or estimated (0) |
| ret_sel_input | Array specifying fixed retention selectivity values when use_fixed_ret_sel = 1 |
| Use_ret_selex_prior | Integer (0 or 1). Flag to enable log-normal priors on retention selectivity parameters |
| retsel_devs_min_shared_bins | Integer vector specifying the reference (minimum) bin index within each shared deviation group for retention selectivity, used when evaluating GMRF or 2D AR(1) likelihoods. Defaults to 1:n_ages when no bin sharing is specified |
| map_ln_retsel_devs | Array indicating which continuous time-varying retention selectivity deviation values are fixed (mapped off) |
| ret_sel_pen_wts | Same format as fish_sel_pen_wts, evaluated on the realized retention selectivity surface |
| retsel_pe_wt, retsel_rw_init_sigma | Same meaning as their fishsel counterparts, for retention selectivity |
| ret_sel_bin_dev_bins, cont_tv_retsel_bin_devs, retsel_bin_devs_rw_init_sigma, map_ln_retsel_bin_devs | Bin-override deviation controls with the same meaning as their fishsel counterparts, for retention selectivity |
| retsel_sex_par_offset, retsel_sex_scale_offset | Sex offset flags with the same meaning as their fishsel counterparts, for retention selectivity (`ret_sel_sex_offset`) |
| Use_ret_selex_penalty, ret_selex_penalty | Centering penalty flag and specification table with the same meaning as their fish counterparts, for retention selectivity |

## Data Inputs for Defining Survey Selectivity and Catchability

| Name | Description |
|----|----|
| cont_tv_srv_sel | Matrix dimensioned by n_regions, n_srv_fleets specifying whether and how continuous time-varying selectivity is applied. 0: None, 1: iid deviations, 2: random walk, 3: 3D GMRF with marginal variance (semi-parametric), 4: 3D GMRF with conditional variance (semi-parametric). Further details in model_priors_penalties.R, model_precision.R, and model_selectivity.R |
| srv_sel_blocks | Array dimensioned by n_regions, n_years, n_srv_fleets specifying selectivity time blocks. Unique integers denote distinct selectivity parameter blocks |
| srv_sel_model | Array dimensioned by n_regions, n_years, n_srv_fleets specifying the selectivity functional form. Same options as fish_sel_model, including bicubic (8). Further details in model_selectivity.R |
| srv_sel_bicubic_binnodes, srv_sel_bicubic_yrnodes, srv_sel_bicubic_selstyr, srv_sel_bicubic_nselbins, srv_sel_bicubic_Wbin, srv_sel_bicubic_Wyr | Same meaning as their fish_sel_bicubic\_\* counterparts, for survey selectivity |
| srv_selex_type | Integer specifying whether survey selectivity is age-based (0) or length-based (1) |
| use_fixed_srv_sel | Integer specifying whether survey selectivity is fixed externally (1) or estimated (0) |
| srv_q_blocks | Array dimensioned by n_regions, n_years, n_srv_fleets specifying catchability time blocks. Unique integers denote distinct catchability parameter blocks |
| Use_srv_q_prior | Survey catchability prior indicator. 0 == don’t use, 1 == use |
| srv_q_prior | Data frame containing prior specifications for survey catchability parameters. Must include columns: region, fleet, block, mu (prior mean on natural scale), and sd (prior standard deviation on log scale). Each row specifies a log-normal prior for one catchability parameter |
| map_ln_srvsel_devs | Array dimensioned by n_regions, n_years, n_ages, n_sexes, n_srv_fleets indicating which continuous time-varying selectivity deviation values are fixed (mapped off) |
| Use_srv_selex_prior | Integer (0 or 1). Flag to enable log-normal priors on survey selectivity parameters as specified in srv_selex_prior |
| srv_selex_prior | Data frame containing prior specifications for survey selectivity parameters. Must include columns: region, fleet, block, sex, par (parameter index), mu (prior mean on natural scale), and sd (prior standard deviation on log scale). Each row specifies a log-normal prior for one selectivity parameter |
| srvsel_devs_min_shared_bins | Integer vector specifying the reference (minimum) bin index within each shared deviation group, used to subset the bin dimension when evaluating GMRF or 2D AR(1) likelihoods (PE models 3-5). Defaults to 1:n_ages when no bin sharing is specified |
| srv_sel_pen_wts | Same format as fish_sel_pen_wts, evaluated on the realized survey selectivity surface |
| srv_q_type | Vector dimensioned by n_srv_fleets specifying how survey catchability is obtained. 0: Estimated as exp(ln_srv_q) (default), 1: Solved analytically as the ratio of mean observed to mean predicted index, 2: Solved analytically on the log scale as exp(mean(log(obs) - log(pred))). Analytic fleets have their ln_srv_q fixed automatically, ignore block structure, and cannot carry catchability covariates or priors |
| srvsel_pe_wt, srvsel_rw_init_sigma | Same meaning as their fishsel counterparts, for survey selectivity |
| srv_sel_bin_dev_bins, cont_tv_srvsel_bin_devs, srvsel_bin_devs_rw_init_sigma, map_ln_srvsel_bin_devs | Bin-override deviation controls with the same meaning as their fishsel counterparts, for survey selectivity |
| srvsel_sex_par_offset, srvsel_sex_scale_offset | Sex offset flags with the same meaning as their fishsel counterparts, for survey selectivity (`srv_sel_sex_offset`) |
| Use_srv_selex_penalty, srv_selex_penalty | Centering penalty flag and specification table with the same meaning as their fish counterparts, for survey selectivity |

## Data Inputs for Defining Model Weighting

| Name | Description |
|----|----|
| Wt_Catch | Weight applied to region-aggregated fishery catch likelihoods. Either a numeric scalar or an array dimensioned by n_regions, n_years, n_seas, n_fish_fleets |
| Wt_Catch_pop | Weight applied to population-specific fishery catch likelihoods. Either a numeric scalar or an array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets |
| Wt_FishIdx | Weight applied to region-aggregated fishery index likelihoods. Either a numeric scalar or an array dimensioned by n_regions, n_years, n_seas, n_fish_fleets |
| Wt_FishIdx_pop | Weight applied to population-specific fishery index likelihoods. Either a numeric scalar or an array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets |
| Wt_SrvIdx | Weight applied to region-aggregated survey index likelihoods. Either a numeric scalar or an array dimensioned by n_regions, n_years, n_seas, n_srv_fleets |
| Wt_SrvIdx_pop | Weight applied to population-specific survey index likelihoods. Either a numeric scalar or an array dimensioned by n_pop, n_regions, n_years, n_seas, n_srv_fleets |
| Wt_Rec | Weight applied to the recruitment deviation penalty. Either a scalar or an array dimensioned by n_pop, n_regions, and the third dimension of ln_RecDevs (which dont_est_recdev_last and n_proj_yrs_devs both change) for per-deviation weighting. A weight of zero excludes a deviation from the penalty while it remains estimated |
| Wt_Init_Rec | Weight applied to the initial age deviation penalty. Either a scalar or an array dimensioned by n_pop, n_regions, n_ages - 1. Defaults to a scalar Wt_Rec; must be supplied explicitly when Wt_Rec is an array, since the two penalties are dimensioned differently |
| Wt_F | Weight applied to fishing mortality deviations |
| Wt_D | Weight applied to discard mortality rate deviations |
| Wt_Tagging | Weight applied to tagging data likelihoods |
| Wt_Discard | Weight applied to region-aggregated discard likelihoods. Either a numeric scalar or an array dimensioned by n_regions, n_years, n_seas, n_fish_fleets |
| Wt_Discard_pop | Weight applied to population-specific discard likelihoods. Either a numeric scalar or an array dimensioned by n_pop, n_regions, n_years, n_seas, n_fish_fleets |
| Wt_FishAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to region-aggregated fishery age compositions, ideally derived using Francis re-weighting |
| Wt_FishAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to population-specific fishery age compositions |
| Wt_FishLenComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to region-aggregated fishery length compositions, ideally derived using Francis re-weighting |
| Wt_FishLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to population-specific fishery length compositions |
| Wt_SrvAgeComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight applied to region-aggregated survey age compositions, ideally derived using Francis re-weighting |
| Wt_SrvAgeComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight applied to population-specific survey age compositions |
| Wt_SrvLenComps | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight applied to region-aggregated survey length compositions, ideally derived using Francis re-weighting |
| Wt_SrvLenComps_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_srv_fleets specifying a multinomial weight applied to population-specific survey length compositions |
| Wt_FishAgeComps_discard | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to region-aggregated discard age compositions |
| Wt_FishAgeComps_discard_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to population-specific discard age compositions |
| Wt_FishLenComps_discard | Array dimensioned by n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to region-aggregated discard length compositions |
| Wt_FishLenComps_discard_pop | Array dimensioned by n_pop, n_regions, n_years, n_seas, n_sexes, n_fish_fleets specifying a multinomial weight applied to population-specific discard length compositions |
