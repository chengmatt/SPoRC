# Description of Model Report

The following tables describes all elements contained within `obj$rep`,
which is generated after a model is run. Note that when `n_sexes > 1`,
the first dimension will always be females and the second dimension will
always be males.

## Biological Processes

| Name | Description | Dimensions |
|----|----|----|
| `R0` | Global mean or virgin recruitment parameter | scalar |
| `Rec_trans_prop` | Recruitment apportionment parameters | `n_regions` |
| `sexratio` | Recruitment sex ratio | `n_regions × n_years × n_sexes` |
| `h_trans` | Steepness parameter by region; ignored if mean recruitment specified | `n_regions` |
| `NAA` | Numbers-at-age | `n_regions × n_years × n_ages × n_sexes` |
| `NAA0` | Unfished Numbers-at-age | `n_regions × n_years × n_ages × n_sexes` |
| `ZAA` | Total instantaneous mortality | `n_regions × n_years × n_ages × n_sexes` |
| `natmort` | Instantaneous natural mortality | `n_regions × n_years × n_ages × n_sexes` |
| `bias_ramp` | Vector of bias ramp values | `n_years` |
| `Movement` | Movement array | `n_regions × n_regions × (n_years + n_proj_yrs_devs) × n_ages × n_sexes` |
| `Mrate` | Instantaneous movement rate matrix (CTMC movement only; NULL otherwise) | `n_regions × n_regions × (n_years + n_proj_yrs_devs) × n_ages × n_sexes` |

## Fishery Processes

| Name | Description | Dimensions |
|----|----|----|
| `init_F` | Initial fishing mortality applied to equilibrium age structure | `n_regions × n_ages × n_sexes` |
| `ln_sigmaC` | Log standard deviation of regional catches | `n_regions × n_years × n_fish_fleets` |
| `ln_sigmaC_agg` | Log standard deviation for aggregated catch | `n_years × n_fish_fleets` |
| `Fmort` | Fishing mortality rate | `n_regions × n_years × n_fish_fleets` |
| `FAA` | Fishing mortality at age | `n_regions × n_years × n_ages × n_sexes × n_fish_fleets` |
| `CAA` | Catch-at-age | `n_regions × n_years × n_ages × n_sexes × n_fish_fleets` |
| `CAL` | Catch-at-length | `n_regions × n_years × n_lens × n_fish_fleets` |
| `PredCatch` | Predicted catch | `n_regions × n_years × n_fish_fleets` |
| `PredFishIdx` | Predicted fishery index | `n_regions × n_years × n_fish_fleets` |
| `fish_sel` | Fishery age-selectivity | `n_regions × (n_years + n_proj_yrs_devs) × n_ages × n_sexes × n_fish_fleets` |
| `fish_sel_l` | Fishery length-selectivity | `n_regions × (n_years + n_proj_yrs_devs) × n_lens × n_sexes × n_fish_fleets` |
| `fish_q` | Fishery catchability | `n_regions × n_years × n_fish_fleets` |

## Survey Processes

| Name | Description | Dimensions |
|----|----|----|
| `PredSrvIdx` | Predicted survey index | `n_regions × n_years × n_srv_fleets` |
| `srv_sel` | Survey age-selectivity | `n_regions × (n_years + n_proj_yrs_devs) × n_ages × n_srv_fleets × n_sexes` |
| `srv_sel_l` | Survey length-selectivity | `n_regions × (n_years + n_proj_yrs_devs) × n_lens × n_srv_fleets × n_sexes` |
| `srv_q` | Survey catchability | `n_regions × n_years × n_srv_fleets` |
| `SrvIAA` | Survey index-at-age | `n_regions × n_years × n_ages × n_sexes × n_srv_fleets` |
| `SrvIAL` | Survey index-at-length | `n_regions × n_years × n_lens × n_sexes × n_srv_fleets` |

## Tagging Processes

| Name | Description | Dimensions |
|----|----|----|
| `Pred_Tag_Recap` | Predicted tag recaptures | `n_tag_liberty × n_years × n_regions × n_ages × n_sexes` |
| `Tags_Avail` | Available tags | `(n_tag_liberty + 1) × n_years × n_regions × n_ages × n_sexes` |
| `Tag_Reporting` | Tag reporting rate | `n_regions × n_years` |

## Likelihoods

| Name | Description | Dimensions |
|----|----|----|
| `Catch_nLL` | Negative log-likelihood for catch | `n_regions × n_years × n_fish_fleets` |
| `FishIdx_nLL` | Negative log-likelihood for fishery index | `n_regions × n_years × n_fish_fleets` |
| `SrvIdx_nLL` | Negative log-likelihood for survey index | `n_regions × n_years × n_srv_fleets` |
| `FishAgeComps_nLL` | Negative log-likelihood for fish age composition | `n_regions × n_years × n_sexes × n_fish_fleets` |
| `SrvAgeComps_nLL` | Negative log-likelihood for survey age composition | `n_regions × n_years × n_sexes × n_srv_fleets` |
| `FishLenComps_nLL` | Negative log-likelihood for fish length composition | `n_regions × n_years × n_sexes × n_fish_fleets` |
| `SrvLenComps_nLL` | Negative log-likelihood for survey length composition | `n_regions × n_years × n_sexes × n_srv_fleets` |
| `M_nLL` | Negative log-likelihood for natural mortality | scalar |
| `Fmort_nLL` | Negative log-likelihood for fishing mortality | `n_regions × n_years × n_fish_fleets` |
| `Rec_nLL` | Negative log-likelihood for recruitment | `n_regions × n_years` |
| `Init_Rec_nLL` | Negative log-likelihood for initial recruitment | `n_regions × (n_ages - 2)` |
| `Tag_nLL` | Negative log-likelihood for tag recapture | `n_tag_liberty × n_years × n_regions × n_ages × n_sexes` |
| `h_nLL` | Negative log-likelihood for steepness | scalar |
| `fish_q_nLL` | Negative log-likelihood for fishery catchability | scalar |
| `sel_nLL` | Negative log-likelihood for selectivity | scalar |
| `srv_q_nLL` | Negative log-likelihood for survey catchability | scalar |
| `Movement_nLL` | Negative log-likelihood for movement | scalar |
| `TagRep_nLL` | Negative log-likelihood for tag reporting | scalar |
| `Rec_prop_nLL` | Negative log-likelihood for recruitment proportions | scalar |
| `jnLL` | Joint negative log-likelihood | scalar |

## Parameter Deviations

| Name | Description | Dimensions |
|----|----|----|
| `ln_RecDevs` | Log recruitment deviations by region and year (including projected years). | `n_regions × (n_years + n_proj_yrs_devs)` |
| `move_devs` | Movement deviations among regions, by year, age, and sex. | `n_regions × (n_regions - 1) × (n_years + n_proj_yrs_devs) × n_ages × n_sexes` |
| `ln_fishsel_devs` | Log-scale deviations in fishery selectivity-at-length or age, by fleet. | `n_regions × (n_years + n_proj_yrs_devs) × bins × n_sexes × n_fish_fleets` |
| `ln_srvsel_devs` | Log-scale deviations in survey selectivity-at-length or age, by fleet. | `n_regions × (n_years + n_proj_yrs_devs) × bins × n_sexes × n_srv_fleets` |

## Derived Quantities (also reported in ADREPORT)

| Name | Description | Dimensions |
|----|----|----|
| `Total_Biom` | Total biomass | `n_regions × n_years` |
| `SSB` | Spawning stock biomass | `n_regions × n_years` |
| `Aggregated_SSB` | Aggregated spawning stock biomass | `n_years` |
| `Dynamic_SSB0` | Unfished spawning stock biomass | `n_regions × n_years` |
| `Dynamic_Aggregated_SSB0` | Unfished aggregated spawning stock biomass | `n_years` |
| `Rec` | Recruitment | `n_regions × n_years` |
