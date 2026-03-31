# Setup values and dimensions of fishing processes

Setup values and dimensions of fishing processes

## Usage

``` r
Setup_Sim_Fishing(
  sim_list,
  ln_sigmaC = array(log(0.02), dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_fish_fleets)),
  catch_units = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
  init_F_val = 0,
  Fmort_input = array(0.1, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_fish_fleets, sim_list$n_sims)),
  fish_sel_input,
  fish_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_fish_fleets, sim_list$n_sims)),
  ObsFishIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_fish_fleets)),
  fish_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_fish_fleets)),
  comp_fishage_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
  ln_FishAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_FishAge_theta_agg = rep(log(1), sim_list$n_fish_fleets),
  FishAge_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
  FishAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets, 2)),
  FishAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets)),
  comp_fishlen_like = rep(0, sim_list$n_fish_fleets),
  ISS_FishLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_sexes, sim_list$n_fish_fleets, sim_list$n_sims)),
  ln_FishLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets)),
  ln_FishLen_theta_agg = rep(log(1), sim_list$n_fish_fleets),
  FishLen_corr_pars_agg = rep(0.01, sim_list$n_fish_fleets),
  FishLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_fish_fleets, 2)),
  FishLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_fish_fleets))
)
```

## Arguments

- sim_list:

  Simulation list object from \`Setup_Sim_Dim()\`

- ln_sigmaC:

  Observation error for catch \[n_regions × n_yrs × n_fish_fleets\]
  (default: \`log(0.02)\`)

- catch_units:

  Units of catch - Array \[n_regions × n_fish_fleets\]

  - `0`: Abundance

  - `1`: Biomass (default)

- init_F_val:

  Initial fishing mortality value (default: \`0\`)

- Fmort_input:

  Fishing mortality input array \[n_regions × n_yrs × n_fish_fleets ×
  n_sims\] (default: \`0.1\`)

- fish_sel_input:

  Fishery selectivity array \[n_regions × n_yrs × n_ages × n_sexes ×
  n_fish_fleets × n_sims\] (no default, must be provided)

- fish_q_input:

  Fishery catchability array \[n_regions × n_yrs × n_fish_fleets ×
  n_sims\] (default: \`1\`)

- ObsFishIdx_SE:

  Observation error of fishery index \[n_regions × n_yrs ×
  n_fish_fleets\] (default: \`0.2\`)

- fish_idx_type:

  Array of index types \[n_regions x n_fish_fleets\] (default: all \`1\`
  = biomass index)

  - `0`: Abundance index

  - `1`: Biomass index

- comp_fishage_like:

  Vector \[n_fish_fleets\] specifying likelihood for simulating age
  comps (default: all \`0\` = multinomial)

  - `0`: Multinomial

  - `1`: Dirichlet-Multinomial

  - `2`: Logistic Normal iid

  - `3`: Logistic Normal 1dar1

  - `4`: Logistic Normal 2d correlation (constant by sex, 1dar1 by age)

- ISS_FishAgeComps:

  Input sample sizes \[n_regions × n_yrs × n_sexes × n_fish_fleets ×
  n_sims\] (default: \`100\`)

- ln_FishAge_theta:

  Overdispersion parameters \[n_regions × n_sexes × n_fish_fleets\]
  (default: \`log(1)\`)

- ln_FishAge_theta_agg:

  Overdispersion parameters for aggregated comps \[n_fish_fleets\]
  (default: \`log(1)\`)

- FishAge_corr_pars_agg:

  Correlation parameters (agg.) for options 3–4 \[n_fish_fleets\]
  (default: \`0.01\`)

- FishAge_corr_pars:

  Correlation parameters \[n_regions × n_sexes × n_fish_fleets x 2\]
  (default: \`0.01\`)

- FishAgeComps_Type:

  Array \[n_yrs × n_fish_fleets\] (default: \`2\` = joint by sex, split
  by region)

  - `0`: Aggregated

  - `1`: Split by sex and region

  - `2`: Joint by sex, split by region

  - `999`: Not simulated

- comp_fishlen_like:

  Vector \[n_fish_fleets\] specifying likelihood for simulating length
  comps (default: all \`0\` = multinomial)

  - `0`: Multinomial

  - `1`: Dirichlet-Multinomial

  - `2`: Logistic Normal iid

  - `3`: Logistic Normal 1dar1

  - `4`: Logistic Normal 2d correlation (constant by sex, 1dar1 by
    length)

- ISS_FishLenComps:

  Input sample sizes \[n_regions × n_yrs × n_sexes × n_fish_fleets ×
  n_sims\] (default: \`100\`)

- ln_FishLen_theta:

  Overdispersion parameters \[n_regions × n_sexes × n_fish_fleets x 2\]
  (default: \`log(1)\`)

- ln_FishLen_theta_agg:

  Overdispersion parameters for aggregated comps \[n_fish_fleets\]
  (default: \`log(1)\`)

- FishLen_corr_pars_agg:

  Correlation parameters (agg.) for options 3–4 \[n_fish_fleets\]
  (default: \`0.01\`)

- FishLen_corr_pars:

  Correlation parameters \[n_regions × n_sexes × n_fish_fleets\]
  (default: \`0.01\`)

- FishLenComps_Type:

  Array \[n_yrs × n_fish_fleets\] (default: \`2\` = joint by sex, split
  by region)

  - `0`: Aggregated

  - `1`: Split by sex and region

  - `2`: Joint by sex, split by region

  - `999`: Not simulated

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
