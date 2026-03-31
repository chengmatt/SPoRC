# Set up simulated tagging dynamics

Set up simulated tagging dynamics

## Usage

``` r
Setup_Sim_Tagging(
  n_tags = NULL,
  n_tags_rel_input = NULL,
  UseTagging = 0,
  max_liberty = sim_list$n_ages/2,
  tag_release_indicator = expand.grid(regions = 1:sim_list$n_regions, tag_years =
    1:sim_list$n_yrs),
  t_tagging = 0,
  ln_Init_Tag_Mort = log(1e-05),
  ln_Tag_Shed = log(1e-05),
  Tag_Reporting_input = array(0.5, dim = c(sim_list$n_regions, sim_list$n_yrs,
    sim_list$n_sims)),
  tag_selex = 5,
  tag_natmort = 3,
  tag_like = 0,
  ln_tag_theta = log(1),
  sim_list
)
```

## Arguments

- n_tags:

  Number of tags to release in a given year (scalar, default = NULL)

- n_tags_rel_input:

  Number of tag releases by tag cohort length (default = NULL)

- UseTagging:

  Boolean to use tagging (default = 0):

  - 0: Do not simulate tagging

  - 1: Simulate tagging

- max_liberty:

  Maximum liberty to track cohorts (default = sim_list\$n_ages / 2)

- tag_release_indicator:

  Tag release indicator \[regions × tag_years\] (default = all
  combinations of regions × years via \`expand.grid\`)

- t_tagging:

  Time of tagging (e.g., start year == 0, mid year == 0.5; default = 0)

- ln_Init_Tag_Mort:

  Log initial tag-induced mortality (default = log(1e-5))

- ln_Tag_Shed:

  Log chronic tag shedding rate (default = log(1e-5))

- Tag_Reporting_input:

  Tag reporting input \[n_regions × n_yrs × n_sims\] (default = 0.5)

- tag_selex:

  Tag selectivity type (integer, default = 5):

  - 0: Uniform_DomFleet

  - 1: SexAgg_DomFleet

  - 2: SexSp_DomFleet

  - 3: Uniform_AllFleet

  - 4: SexAgg_AllFleet

  - 5: SexSp_AllFleet

- tag_natmort:

  Tag natural mortality type (integer, default = 3):

  - 0: AgeAgg_SexAgg

  - 1: AgeSp_SexAgg

  - 2: AgeAgg_SexSp

  - 3: AgeSp_SexSp

- tag_like:

  Tag likelihood type (integer, default = 0):

  - 0: Poisson

  - 1: NegBin

  - 2: Multinomial_Release

  - 3: Multinomial_Recapture

  - 4: Dirichlet-Multinomial_Release

  - 5: Dirichlet-Multinomial_Recapture

- ln_tag_theta:

  Scalar in log space describing tag likelihood overdispersion (default
  = log(1))

- sim_list:

  Simulation list (required)

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_data_to_SPoRC()`](https://chengmatt.github.io/SPoRC/reference/simulation_data_to_SPoRC.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
