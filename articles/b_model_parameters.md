# Description of Model Parameters

The following table describes all elements contained within
`input_list$pars`, which are generated using the `SPoRC::Setup_x`
functions. Note that when `n_sexes > 1`, the first dimension will always
be females and the second dimension will always be males.

## Fishery Removals

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_sigmaC` | `n_regions, n_years, n_fish_fleets` | Standard deviation for catch likelihood | 0.01 |
| `ln_sigmaC_agg` | `n_years, n_fish_fleets` | Standard deviation for region-aggregated catch likelihood | 0.01 |
| `ln_sigmaF` | `n_regions, n_fish_fleets` | Standard deviation for fishing mortality process error | 0 |
| `ln_sigmaF_agg` | `n_fish_fleets` | Standard deviation for region-aggregated fishing mortality process error | 0 |
| `ln_F_mean` | `n_regions, n_fish_fleets` | Log-scale mean fishing mortality rate | 0.1 |
| `ln_F_devs` | `n_regions, n_years, n_fish_fleets` | Log-scale annual fishing mortality deviations (can be specified as random effects) | 0 |
| `ln_F_mean_AggCatch` | `n_fish_fleets` | Log-scale mean fishing mortality for region-aggregated catch | 0.1 |
| `ln_F_devs_AggCatch` | `n_Catch_Type_years, n_fish_fleets` | Log-scale fishing mortality deviations for region-aggregated catch (can be specified as random effects) | 0 |

## Fishery Selectivity

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_fish_fixed_sel_pars` | `n_regions, max_fish_pars, max_fishsel_blks, n_sexes, n_fish_fleets` | Log-scale selectivity curve parameters | 0 |
| `fishsel_pe_pars` | `n_regions, 4, n_sexes, n_fish_fleets` | Process error parameters for time-varying selectivity | 0 |
| `ln_fishsel_devs` | `n_regions, n_years + n_proj_yrs_devs, n_bins, n_sexes, n_fish_fleets` | Log-scale selectivity deviations w/ projection years (n_bins = n_ages or n_lens depending on selectivity type) (can be specified as random effects) | 0 |

## Fishery Catchability

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_fish_q` | `n_regions, max_fishq_blks, n_fish_fleets` | Log-scale fishery catchability coefficient | 0 |

## Fishery Composition Likelihoods

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_FishAge_theta` | `n_regions, n_fish_fleets, n_sexes` | Log-scale overdispersion parameter for fishery age compositions | 0 |
| `ln_FishAge_theta_agg` | `n_fish_fleets` | Log-scale overdispersion for region-aggregated fishery age compositions | 0 |
| `FishAge_corr_pars` | `n_regions, n_sexes, n_fish_fleets, 2` | Correlation parameters for logistic-normal fishery age composition likelihood | 0.01 |
| `FishAge_corr_pars_agg` | `n_fish_fleets` | Correlation parameters for region-aggregated logistic-normal fishery age compositions | 0.01 |
| `ln_FishLen_theta` | `n_regions, n_fish_fleets, n_sexes` | Log-scale overdispersion parameter for fishery length compositions | 0 |
| `ln_FishLen_theta_agg` | `n_fish_fleets` | Log-scale overdispersion for region-aggregated fishery length compositions | 0 |
| `FishLen_corr_pars` | `n_regions, n_sexes, n_fish_fleets, 2` | Correlation parameters for logistic-normal fishery length composition likelihood. The first parameter represents bin correlations, while the second parameter represents sex correlations (if specified) | 0.01 |
| `FishLen_corr_pars_agg` | `n_fish_fleets` | Correlation parameters for region- and sex-aggregated logistic-normal fishery length compositions. The correlation represents bin correlations. | 0.01 |

## Survey Selectivity

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_srv_fixed_sel_pars` | `n_regions, max_srv_pars, max_srv_blks, n_sexes, n_srv_fleets` | Log-scale selectivity curve parameters | 0 |
| `srvsel_pe_pars` | `n_regions, 4, n_sexes, n_srv_fleets` | Process error parameters for time-varying selectivity | 0 |
| `ln_srvsel_devs` | `n_regions, n_years + n_proj_yrs_devs, n_bins, n_sexes, n_srv_fleets` | Log-scale selectivity deviations w/ projection years (n_bins = n_ages or n_lens depending on selectivity type) (can be specified as random effects) | 0 |

## Survey Catchability

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_srv_q` | `n_regions, max_srvq_blks, n_srv_fleets` | Log-scale survey catchability coefficient | 0 |
| `srv_q_coeff` | `n_regions, n_srv_fleets, n_covariates` | Regression coefficients for covariate-based survey catchability | 0 |

## Survey Composition Likelihoods

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_SrvAge_theta` | `n_regions, n_srv_fleets, n_sexes` | Log-scale overdispersion parameter for survey age compositions | 0 |
| `ln_SrvAge_theta_agg` | `n_srv_fleets` | Log-scale overdispersion for region-aggregated survey age compositions | 0 |
| `SrvAge_corr_pars` | `n_regions, n_sexes, n_srv_fleets, 2` | Correlation parameters for logistic-normal survey age composition likelihood | 0.01 |
| `SrvAge_corr_pars_agg` | `n_srv_fleets` | Correlation parameters for region-aggregated logistic-normal survey age compositions | 0.01 |
| `ln_SrvLen_theta` | `n_regions, n_srv_fleets, n_sexes` | Log-scale overdispersion parameter for survey length compositions | 0 |
| `ln_SrvLen_theta_agg` | `n_srv_fleets` | Log-scale overdispersion for region-aggregated survey length compositions | 0 |
| `SrvLen_corr_pars` | `n_regions, n_sexes, n_srv_fleets, 2` | Correlation parameters for logistic-normal survey length composition likelihood. The first parameter represents bin correlations, while the second parameter represents sex correlations (if specified) | 0.01 |
| `SrvLen_corr_pars_agg` | `n_srv_fleets` | Correlation parameters for region- and sex-aggregated logistic-normal survey length compositions. The correlation represents bin correlations. | 0.01 |

## Natural Mortality

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_M` | `n_regions_block, n_years_block, n_ages_block, n_sexes_block` | Log-scale natural mortality rate (dimensioned by user-defined blocks) | 0.5 |

## Recruitment

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_global_R0` | `1` | Log-scale unfished/mean recruitment | 15 |
| `Rec_prop` | `n_regions - 1` | Regional recruitment proportions (simplex parameterization) | 0 |
| `steepness_h` | `n_regions` | Steepness of stock-recruitment relationship (bounded logit scale) | 0 |
| `ln_sigmaR` | `2` | Log-scale standard deviations for initial age deviations and recruitment deviations | 0 |
| `ln_InitDevs` | `n_regions, n_ages - 2` | Log-scale initial age structure deviations (can be specified as random effects) | 0 |
| `ln_RecDevs` | `n_regions, n_recdev_years + n_proj_yrs_devs` | Log-scale annual recruitment deviations w/ projection years (can be specified as random effects) | 0 |
| `sexratio_pars` | `n_regions, n_blocks` | Logit-scale sex ratio parameters (defaults to 50:50 when n_sexes = 1) | 0 |

## Movement

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `move_pars` | `n_regions, n_regions - 1, n_years, n_ages, n_sexes` | Logit-scale movement probabilities (unstructured Markov; move_type = 0) | 0 |
| `log_move_diffusion_pars` | (Vector; Depends on formula) | Log-scale diffusion parameters for CTMC movement (move_type = 1) | log(0.1) |
| `move_preference_pars` | (Vector; Depends on formula) | Preference parameters for CTMC movement (move_type = 1) | 0 |
| `move_devs` | `n_regions, n_regions - 1, n_years + n_proj_yrs_devs, n_ages, n_sexes` | Log Movement deviations w/ projection years (can be specified as random effects) | 0 |
| `move_pe_pars` | `n_regions, max(4, n_ages), n_sexes` | Process error parameters for continuously varying movement | 0 |

## Tagging

| Parameter | Dimension | Description | Default |
|----|----|----|----|
| `ln_Init_Tag_Mort` | `1` | Log-scale initial tagging mortality rate | 1e-50 |
| `ln_Tag_Shed` | `1` | Log-scale chronic tag shedding rate | 1e-50 |
| `Tag_Reporting_Pars` | `n_regions, n_blocks` | Logit-scale tag reporting rates | 0 |
| `ln_tag_theta` | `1` | Log-scale overdispersion parameter for tag recapture likelihood | 0 |
