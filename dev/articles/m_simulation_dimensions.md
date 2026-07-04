# Description of Simulation Dimensions

The following tables describe key elements generated from simulations
conducted in `SPoRC`, which can be input into an estimation framework
and utilized for comparison.

## Simulation Settings

| Element | Dimension | Description |
|----|----|----|
| `n_sims` | `1` | Number of simulations |
| `n_pop` | `1` | Number of populations |
| `n_regions` | `1` | Number of regions |
| `n_years` | `1` | Number of years (alias for `n_yrs`) |
| `n_yrs` | `1` | Number of years |
| `n_ages` | `1` | Number of ages |
| `n_lens` | `1` | Number of length bins |
| `n_sexes` | `1` | Number of sexes |
| `n_seas` | `1` | Number of seasons |
| `n_fish_fleets` | `1` | Number of fishery fleets |
| `n_srv_fleets` | `1` | Number of survey fleets |
| `seasdur` | `n_seas` | Duration of each season as a fraction of the year |
| `spawn_seas` | `1` | Season in which spawning occurs |
| `natal_region` | `n_pop` | Natal region index for each population |

## Population, Recruitment, and Biological Processes

| Element | Dimension | Description |
|----|----|----|
| `sexratio` | `n_pop × n_regions × n_years × n_sexes × n_sims` | Sex ratio (female proportion) by population, region, year, and sex |
| `R0` | `n_pop × n_regions × n_years × n_sims` | Regional unfished or mean recruitment by population |
| `Rec` | `n_pop × n_regions × n_years × n_sims` | Realized recruitment by population and region |
| `rec_seas_prop` | `n_pop × n_seas × n_sims` | Seasonal recruitment apportionment proportions by population |
| `stray_rate` | `n_pop × n_years × n_sims` | Proportion of each population’s spawning biomass contributing to non-natal regions |
| `NAA` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Numbers-at-age (fished) |
| `NAA0` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Numbers-at-age (unfished) |
| `NAA_bef` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Numbers-at-age before movement is applied |
| `NAA_aft` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Numbers-at-age after movement is applied |
| `Dynamic_SSB0` | `n_pop × n_regions × n_years × n_sims` | Dynamic unfished spawning stock biomass by population and region |
| `SSB` | `n_pop × n_regions × n_years × n_sims` | Spawning stock biomass by population and region |
| `eff_SSB` | `n_pop × n_years × n_sims` | Effective SSB at each population’s natal region, including stray contributions |
| `Total_Biom` | `n_pop × n_regions × n_years × n_sims` | Total biomass by population and region |
| `ln_RecDevs` | `n_pop × n_regions × n_years × n_sims` | Log-scale recruitment deviations |
| `ln_InitDevs` | `n_pop × n_regions × (n_ages - 1) × n_sims` | Log-scale initial age deviations |
| `natmort` | `n_pop × n_regions × n_years × n_ages × n_sexes × n_sims` | Natural mortality-at-age |
| `ZAA` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Total instantaneous mortality (natural + retained F + dead discard F) |
| `WAA` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Mean weight-at-age |
| `WAA_fish` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims` | Fishery weight-at-age |
| `WAA_srv` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims` | Survey weight-at-age |
| `MatAA` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Proportion mature-at-age |
| `SizeAgeTrans` | `n_pop × n_regions × n_years × n_seas × n_lens × n_ages × n_sexes × n_sims` | Size–age transition array |
| `AgeingError` | `n_years × n_model_ages × n_obs_ages × n_sims` | Ageing error transition matrix |
| `Movement` | `n_pop × n_regions × n_regions × n_years × n_seas × n_ages × n_sexes × n_sims` | Movement (transition) matrix |
| `sgl_seas_spawning_movement` | `n_pop × n_regions × n_regions × n_years × n_ages × n_sexes × n_sims` | Spawning movement applied when `n_seas == 1` and `n_pop > 1` |

## Fishery Processes

### Fishing Mortality and Selectivity

| Element | Dimension | Description |
|----|----|----|
| `init_F` | `n_seas` | Initial fishing mortality rate |
| `Fmort` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Annual fishing mortality rate by region, season, and fleet |
| `fish_sel` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims` | Total fishery selectivity-at-age |
| `ret_sel` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims` | Retention selectivity-at-age (fraction of selected fish that are retained) |
| `dmr` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Discard mortality rate (fraction of discarded fish that die); used to compute dead discard fishing mortality |
| `fish_q` | `n_regions × n_years × n_fish_fleets × n_sims` | Fishery catchability |

### Retained Catch

| Element | Dimension | Description |
|----|----|----|
| `ln_sigmaC` | `n_regions × n_years × n_seas × n_fish_fleets` | Log standard deviation for regional retained catch likelihood |
| `ln_sigmaC_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets` | Log standard deviation for population-specific retained catch likelihood |
| `CAA` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims` | Retained catch-at-age (Baranov equation, retained component only) |
| `CAL` | `n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims` | Retained catch-at-length |
| `TrueCatch` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | True regional retained catch (abundance or biomass) |
| `ObsCatch` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Observed regional retained catch (with lognormal error) |
| `TrueCatch_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets × n_sims` | True population-specific retained catch |
| `ObsCatch_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Observed population-specific retained catch (with lognormal error) |

### Discards

| Element | Dimension | Description |
|----|----|----|
| `ln_sigmaD` | `n_regions × n_years × n_seas × n_fish_fleets` | Log standard deviation for regional discard likelihood |
| `ln_sigmaD_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets` | Log standard deviation for population-specific discard likelihood |
| `DAA` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets × n_sims` | Dead discard catch-at-age (Baranov equation, dead discard component: `F × fish_sel × (1 - ret_sel) × dmr / Z`) |
| `DAL` | `n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims` | Dead discard catch-at-length |
| `TrueDiscard` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | True regional discard index (units determined by `discard_units`: abundance, biomass, or fraction) |
| `ObsDiscard` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Observed regional discard index (with lognormal error) |
| `TrueDiscard_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets × n_sims` | True population-specific discard index |
| `ObsDiscard_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Observed population-specific discard index (with lognormal error) |

### Fishery Indices

| Element | Dimension | Description |
|----|----|----|
| `TrueFishIdx` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | True regional fishery index (abundance or biomass) |
| `ObsFishIdx` | `n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Observed regional fishery index (with lognormal error) |
| `ObsFishIdx_SE` | `n_regions × n_years × n_seas × n_fish_fleets` | Standard error of observed regional fishery index |
| `TrueFishIdx_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets × n_sims` | True population-specific fishery index |
| `ObsFishIdx_pop` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets × n_sims` | Observed population-specific fishery index |
| `ObsFishIdx_pop_SE` | `n_pop × n_regions × n_years × n_seas × n_fish_fleets` | Standard error of population-specific fishery index |

### Retained Fishery Compositions

| Element | Dimension | Description |
|----|----|----|
| `ObsFishAgeComps` | `n_regions × n_years × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims` | Observed regional retained fishery age compositions |
| `ObsFishLenComps` | `n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims` | Observed regional retained fishery length compositions |
| `ISS_FishAgeComps` | `n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for regional retained fishery age compositions |
| `ISS_FishLenComps` | `n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for regional retained fishery length compositions |
| `ObsFishAgeComps_pop` | `n_pop × n_regions × n_years × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims` | Observed population-specific retained fishery age compositions |
| `ObsFishLenComps_pop` | `n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims` | Observed population-specific retained fishery length compositions |
| `ISS_FishAgeComps_pop` | `n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for population-specific retained fishery age comps |
| `ISS_FishLenComps_pop` | `n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for population-specific retained fishery length comps |

### Discard Fishery Compositions

| Element | Dimension | Description |
|----|----|----|
| `ObsFishAgeComps_discard` | `n_regions × n_years × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims` | Observed regional discard fishery age compositions |
| `ObsFishLenComps_discard` | `n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims` | Observed regional discard fishery length compositions |
| `ISS_FishAgeComps_discard` | `n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for regional discard fishery age compositions |
| `ISS_FishLenComps_discard` | `n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for regional discard fishery length compositions |
| `ObsFishAgeComps_discard_pop` | `n_pop × n_regions × n_years × n_seas × n_obs_ages × n_sexes × n_fish_fleets × n_sims` | Observed population-specific discard fishery age compositions |
| `ObsFishLenComps_discard_pop` | `n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets × n_sims` | Observed population-specific discard fishery length compositions |
| `ISS_FishAgeComps_discard_pop` | `n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for population-specific discard fishery age comps |
| `ISS_FishLenComps_discard_pop` | `n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets × n_sims` | Input sample size for population-specific discard fishery length comps |

## Survey Processes

### Survey Selectivity and Indices

| Element | Dimension | Description |
|----|----|----|
| `srv_sel` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims` | Survey selectivity-at-age |
| `srv_q` | `n_regions × n_years × n_srv_fleets × n_sims` | Survey catchability |
| `TrueSrvIdx` | `n_regions × n_years × n_seas × n_srv_fleets × n_sims` | True regional survey index (abundance or biomass) |
| `ObsSrvIdx` | `n_regions × n_years × n_seas × n_srv_fleets × n_sims` | Observed regional survey index (with lognormal error) |
| `ObsSrvIdx_SE` | `n_regions × n_years × n_seas × n_srv_fleets` | Standard error of observed regional survey index |
| `TrueSrvIdx_pop` | `n_pop × n_regions × n_years × n_seas × n_srv_fleets × n_sims` | True population-specific survey index |
| `ObsSrvIdx_pop` | `n_pop × n_regions × n_years × n_seas × n_srv_fleets × n_sims` | Observed population-specific survey index (with lognormal error) |
| `ObsSrvIdx_pop_SE` | `n_pop × n_regions × n_years × n_seas × n_srv_fleets` | Standard error of population-specific survey index |

### Survey Compositions

| Element | Dimension | Description |
|----|----|----|
| `SrvIAA` | `n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets × n_sims` | Survey index-at-age |
| `SrvIAL` | `n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims` | Survey index-at-length |
| `ObsSrvAgeComps` | `n_regions × n_years × n_seas × n_obs_ages × n_sexes × n_srv_fleets × n_sims` | Observed regional survey age compositions |
| `ObsSrvLenComps` | `n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims` | Observed regional survey length compositions |
| `ISS_SrvAgeComps` | `n_regions × n_years × n_seas × n_sexes × n_srv_fleets × n_sims` | Input sample size for regional survey age compositions |
| `ISS_SrvLenComps` | `n_regions × n_years × n_seas × n_sexes × n_srv_fleets × n_sims` | Input sample size for regional survey length compositions |
| `ObsSrvAgeComps_pop` | `n_pop × n_regions × n_years × n_seas × n_obs_ages × n_sexes × n_srv_fleets × n_sims` | Observed population-specific survey age compositions |
| `ObsSrvLenComps_pop` | `n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets × n_sims` | Observed population-specific survey length compositions |
| `ISS_SrvAgeComps_pop` | `n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets × n_sims` | Input sample size for population-specific survey age compositions |
| `ISS_SrvLenComps_pop` | `n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets × n_sims` | Input sample size for population-specific survey length compositions |

## Tagging

| Element | Dimension | Description |
|----|----|----|
| `use_conv_fish_tagging` | `n_fish_fleets` | Indicator vector for whether conventional tagging is simulated for each fishery fleet |
| `conv_tag_release_indicator` | `n_conv_tag_cohorts × 3` (col 1 = release region, col 2 = release year, col 3 = release season) | Indicator for tag release events |
| `conv_tag_fish_reporting` | `n_regions × n_years × n_fish_fleets × n_sims` | Tag reporting rate by region, year, and fleet |
| `conv_tagged_fish` | `n_conv_tag_cohorts × n_pop × n_ages × n_sexes × n_sims` | Number of tagged fish released per cohort, population, age, and sex |
| `conv_tagged_fish_attr` | `n_conv_tag_cohorts × n_pop × n_ages × n_sexes × n_sims` | Tagged fish release array marginalised to the attended likelihood resolution (via `conv_fish_tag_attr`) |
| `ln_init_conv_tag_mort` | `1` | Log-scale initial tag-induced mortality |
| `ln_conv_tag_shed` | `1` | Log-scale chronic tag shedding rate |
| `conv_tag_fish_avail` | `(conv_tag_max_liberty + 1) × n_seas × n_conv_tag_cohorts × n_pop × n_regions × n_ages × n_sexes × n_sims` | Tagged fish available for recapture at each time step |
| `pred_conv_tag_fish_recap` | `conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop × n_regions × n_ages × n_sexes × n_fish_fleets × n_sims` | Predicted conventional tag recaptures (retained catch only) |
| `obs_conv_tag_fish_recap` | `conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop × n_regions × n_ages × n_sexes × n_fish_fleets × n_sims` | Observed (simulated) conventional tag recaptures |
