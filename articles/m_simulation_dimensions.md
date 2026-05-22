# Description of Simulation Dimensions

The following tables describes key elements generated from simulations
conducted in `SPoRC`, for which can be input into an estimation
framework and utilized for comparison.

## Simulation Settings

| Element         | Dimension | Description              |
|-----------------|-----------|--------------------------|
| `n_sims`        | `1`       | Number of simulations    |
| `n_regions`     | `1`       | Number of regions        |
| `n_years`       | `1`       | Number of years          |
| `n_ages`        | `1`       | Number of ages           |
| `n_lens`        | `1`       | Number of length bins    |
| `n_sexes`       | `1`       | Number of sexes          |
| `n_fish_fleets` | `1`       | Number of fishery fleets |
| `n_srv_fleets`  | `1`       | Number of survey fleets  |

## Population, Recruitment, and Biological Processes

| Element | Dimension | Description |
|----|----|----|
| `sexratio` | `n_regions, n_years, n_sexes, n_sims` | Sex ratio values |
| `R0` | `n_regions, n_years, n_sims` | Regional Unfished or mean recruitment |
| `Rec` | `n_regions, n_years, n_sims` | Realized recruitment |
| `NAA` | `n_regions, n_ages, n_years, n_sexes, n_sims` | Numbers-at-age |
| `NAA0` | `n_regions, n_ages, n_years, n_sexes, n_sims` | Dynamic unfished numbers-at-age |
| `Dynamic_SSB0` | `n_regions, n_years, n_sims` | Dynamic unfished spawning stock biomass |
| `SSB` | `n_regions, n_years, n_sims` | Spawning stock biomass |
| `Total_Biom` | `n_regions, n_years, n_sims` | Total biomass |
| `ln_RecDevs` | `n_regions, n_years, n_sims` | Log-scale recruitment deviations |
| `ln_InitDevs` | `n_regions, n_ages - 1, n_sims` | Log-scale initial age deviations |
| `natmort` | `n_regions, n_ages, n_sexes, n_sims` | Natural mortality-at-age |
| `ZAA` | `n_regions, n_years, n_ages, n_sexes, n_sims` | Total mortality |
| `WAA` | `n_regions, n_years, n_ages, n_sexes, n_sims` | Mean weight-at-age |
| `WAA_fish` | `n_regions, n_years, n_ages, n_sexes, n_fish_fleets, n_sims` | Weight-at-age for fishery |
| `WAA_srv` | `n_regions, n_years, n_ages, n_sexes, n_srv_fleets, n_sims` | Weight-at-age for survey |
| `MatAA` | `n_regions, n_years, n_ages, n_sexes, n_sims` | Proportion mature-at-age |
| `SizeAgeTrans` | `n_regions, n_years, n_lens, n_ages, n_sexes, n_sims` | Size–age transition array |
| `AgeingError` | `n_years, n_model_ages, n_obs_ages, n_sims` | Ageing error array |
| `Movement` | `n_regions, n_regions, n_years, n_ages, n_sexes, n_sims` | Movement array |

## Fishery Processes

| Element | Dimension | Description |
|----|----|----|
| `init_F` | `n_regions, n_fish_fleets` | Initial fishing mortality rate |
| `Fmort` | `n_regions, n_years, n_fish_fleets, n_sims` | Annual fishing mortality |
| `ln_sigmaC` | `n_regions, n_years, n_fish_fleets` | Standard deviation for catch likelihood |
| `fish_sel` | `n_regions, n_years, n_ages, n_sexes, n_fish_fleets, n_sims` | Fishery selectivity-at-age |
| `fish_q` | `n_regions, n_years, n_fish_fleets, n_sims` | Fishery catchability |
| `TrueCatch` | `n_regions, n_years, n_fish_fleets, n_sims` | True catch |
| `ObsCatch` | `n_regions, n_years, n_fish_fleets, n_sims` | Observed catch |
| `ObsFishIdx` | `n_regions, n_years, n_fish_fleets, n_sims` | Observed fishery index |
| `TrueFishIdx` | `n_regions, n_years, n_fish_fleets, n_sims` | True fishery index |
| `ObsFishIdx_SE` | `n_regions, n_years, n_fish_fleets, n_sims` | Standard error of observed index |
| `CAA` | `n_regions, n_years, n_ages, n_sexes, n_fish_fleets, n_sims` | Catch-at-age |
| `CAL` | `n_regions, n_years, n_lens, n_sexes, n_fish_fleets, n_sims` | Catch-at-length |
| `ObsFishAgeComps` | `n_regions, n_years, n_ages, n_sexes, n_fish_fleets, n_sims` | Observed fishery age comps |
| `ObsFishLenComps` | `n_regions, n_years, n_lens, n_sexes, n_fish_fleets, n_sims` | Observed fishery length comps |
| `ISS_FishAgeComps` | `n_regions, n_years, n_fish_fleets, n_sexes, n_sims` | Input sample size for fishery age comps |
| `ISS_FishLenComps` | `n_regions, n_years, n_fish_fleets, n_sexes, n_sims` | Input sample size for fishery length comps |

## Survey Processes

| Element | Dimension | Description |
|----|----|----|
| `srv_sel` | `n_regions, n_years, n_ages, n_sexes, n_srv_fleets, n_sims` | Survey selectivity-at-age |
| `srv_q` | `n_regions, n_years, n_srv_fleets, n_sims` | Survey catchability |
| `ObsSrvIdx` | `n_regions, n_years, n_srv_fleets, n_sims` | Observed survey index |
| `TrueSrvIdx` | `n_regions, n_years, n_srv_fleets, n_sims` | True survey index |
| `ObsSrvIdx_SE` | `n_regions, n_years, n_srv_fleets, n_sims` | Standard error of observed index |
| `Srv_IAA` | `n_regions, n_years, n_ages, n_sexes, n_srv_fleets, n_sims` | Survey index-at-age |
| `Srv_IAL` | `n_regions, n_years, n_lens, n_sexes, n_srv_fleets, n_sims` | Survey index-at-length |
| `ObsSrvAgeComps` | `n_regions, n_years, n_ages, n_sexes, n_srv_fleets, n_sims` | Observed survey age comps |
| `ObsSrvLenComps` | `n_regions, n_years, n_lens, n_sexes, n_srv_fleets, n_sims` | Observed survey length comps |
| `ISS_SrvAgeComps` | `n_regions, n_years, n_srv_fleets, n_sexes, n_sims` | Input sample size for survey age comps |
| `ISS_SrvLenComps` | `n_regions, n_years, n_srv_fleets, n_sexes, n_sims` | Input sample size for survey length comps |

## Tagging

| Element | Dimension | Description |
|----|----|----|
| `tag_release_indicator` | `n_tag_cohorts, 2 (Column 1 = tag region, Column 2 = tag year)` | Indicator for tag release events |
| `Tag_Reporting` | `n_regions, n_years, n_sims` | Tag reporting rate |
| `Tagged_Fish` | `n_regions, n_years, n_ages, n_sexes, n_sims` | Number of tagged fish |
| `ln_Init_Tag_Mort` | `1` | Initial tag-induced mortality |
| `ln_Tag_Shed` | `1` | Tag shedding rate |
| `Tag_Avail` | `n_regions, n_years, n_ages, n_sexes, n_sims` | Available tagged fish for recapture |
| `UseTagging` | `1` | Indicator if tagging is simulated |
| `Pred_Tag_Recap` | `n_regions, n_years, n_ages, n_sexes, n_sims` | True tag recaptures |
| `Obs_Tag_Recap` | `n_regions, n_years, n_ages, n_sexes, n_sims` | Observed tag recaptures |
