# Setup values for survey parameterization

Setup values for survey parameterization

## Usage

``` r
Setup_Sim_Survey(
  ObsSrvIdx_SE = array(0.2, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_srv_fleets)),
  sim_list,
  srv_sel_input,
  srv_q_input = array(1, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_srv_fleets, sim_list$n_sims)),
  t_srv = array(0, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)),
  srv_idx_type = array(1, dim = c(sim_list$n_regions, sim_list$n_srv_fleets)),
  comp_srvage_like = rep(0, sim_list$n_srv_fleets),
  ISS_SrvAgeComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
  ln_SrvAge_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets)),
  ln_SrvAge_theta_agg = rep(log(1), sim_list$n_srv_fleets),
  SrvAge_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
  SrvAge_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets, 2)),
  SrvAgeComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets)),
  comp_srvlen_like = rep(0, sim_list$n_srv_fleets),
  ISS_SrvLenComps = array(100, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_sexes, sim_list$n_srv_fleets, sim_list$n_sims)),
  ln_SrvLen_theta = array(log(1), dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets)),
  ln_SrvLen_theta_agg = rep(log(1), sim_list$n_srv_fleets),
  SrvLen_corr_pars_agg = rep(0.01, sim_list$n_srv_fleets),
  SrvLen_corr_pars = array(0.01, dim = c(sim_list$n_regions, sim_list$n_sexes,
    sim_list$n_srv_fleets, 2)),
  SrvLenComps_Type = array(2, dim = c(sim_list$n_yrs, sim_list$n_srv_fleets))
)
```

## Arguments

- ObsSrvIdx_SE:

  Survey index observation error \[n_regions × n_yrs × n_srv_fleets\]
  (default: \`0.2\`)

- sim_list:

  Simulation list object from \`Setup_Sim_Dim()\`

- srv_sel_input:

  Survey selectivity array \[n_regions × n_yrs × n_ages × n_sexes ×
  n_srv_fleets × n_sims\] (no default, must be provided)

- srv_q_input:

  Survey catchability array \[n_regions × n_yrs × n_srv_fleets ×
  n_sims\] (default: \`1\`)

- t_srv:

  Survey timing fraction \[n_regions × n_srv_fleets\] (default: \`0\`)

- srv_idx_type:

  Array of index types \[n_regions x n_srv_fleets\] (default: all \`1\`
  = biomass index)

  - `0`: Abundance index

  - `1`: Biomass index

- comp_srvage_like:

  Vector \[n_srv_fleets\] specifying likelihood for simulating age comps
  (default: all \`0\` = multinomial)

  - `0`: Multinomial

  - `1`: Dirichlet-Multinomial

  - `2`: Logistic Normal iid

  - `3`: Logistic Normal 1dar1

  - `4`: Logistic Normal 2d correlation (constant by sex, 1dar1 by age)

- ISS_SrvAgeComps:

  Input sample sizes \[n_regions × n_yrs × n_sexes × n_srv_fleets ×
  n_sims\] (default: \`100\`)

- ln_SrvAge_theta:

  Overdispersion parameters \[n_regions × n_sexes × n_srv_fleets\]
  (default: \`log(1)\`)

- ln_SrvAge_theta_agg:

  Overdispersion parameters for aggregated comps \[n_srv_fleets\]
  (default: \`log(1)\`)

- SrvAge_corr_pars_agg:

  Correlation parameters (agg.) for options 3–4 \[n_srv_fleets\]
  (default: \`0.01\`)

- SrvAge_corr_pars:

  Correlation parameters \[n_regions × n_sexes × n_srv_fleets x 2\]
  (default: \`0.01\`)

- SrvAgeComps_Type:

  Array \[n_yrs × n_srv_fleets\] (default: \`2\` = joint by sex, split
  by region)

  - `0`: Aggregated

  - `1`: Split by sex and region

  - `2`: Joint by sex, split by region

  - `999`: Not simulated

- comp_srvlen_like:

  Vector \[n_srv_fleets\] specifying likelihood for simulating length
  comps (default: all \`0\` = multinomial)

  - `0`: Multinomial

  - `1`: Dirichlet-Multinomial

  - `2`: Logistic Normal iid

  - `3`: Logistic Normal 1dar1

  - `4`: Logistic Normal 2d correlation (constant by sex, 1dar1 by
    length)

- ISS_SrvLenComps:

  Input sample sizes \[n_regions × n_yrs × n_sexes × n_srv_fleets ×
  n_sims\] (default: \`100\`)

- ln_SrvLen_theta:

  Overdispersion parameters \[n_regions × n_sexes × n_srv_fleets\]
  (default: \`log(1)\`)

- ln_SrvLen_theta_agg:

  Overdispersion parameters for aggregated comps \[n_srv_fleets\]
  (default: \`log(1)\`)

- SrvLen_corr_pars_agg:

  Correlation parameters (agg.) for options 3–4 \[n_srv_fleets\]
  (default: \`0.01\`)

- SrvLen_corr_pars:

  Correlation parameters \[n_regions × n_sexes × n_srv_fleets x 2\]
  (default: \`0.01\`)

- SrvLenComps_Type:

  Array \[n_yrs × n_srv_fleets\] (default: \`2\` = joint by sex, split
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
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
