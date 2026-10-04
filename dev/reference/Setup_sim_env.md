# Construct and populate a simulation execution environment

Creates a new R environment populated with all objects from `sim_list`
and binds the SPoRC simulation functions required by
[`run_annual_cycle`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md).
Keeping the simulation state in its own environment lets the annual
cycle helpers update shared arrays in place from inside
[`with()`](https://rdrr.io/r/base/with.html) blocks, without threading
them through every call.

## Usage

``` r
Setup_sim_env(sim_list)
```

## Arguments

- sim_list:

  Named list returned by
  [`Setup_Sim_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)
  (or the last upstream setup function called). All elements are copied
  into the new environment via `list2env`.

## Value

A new environment (parent = calling frame) containing every element of
`sim_list` as a named object, plus bound references to the following
SPoRC simulation functions: `generate_initial_age_structure`,
`generate_recruitment`, `apply_pop_dy`, `compute_biom_y_sim`,
`generate_fishery_catch_comp_idx`, `generate_survey_comp_idx`,
`release_conv_tags`, `generate_fishery_conv_tags_recap`,
`Get_Det_Recruitment`, `Get_Init_NAA`, `predict_sim_fish_iss_fmort`,
`rho_trans`, `simulate_comps`, `simulate_conv_tag_fish_recaptures`,
`draw_index_obs`, `resolve_idx_factor`.

## Details

Each helper first stores the environment inside itself as `sim_env`, the
same object under a second name, so its
[`with()`](https://rdrr.io/r/base/with.html) block writes to the
environment it was given whatever the script calls it; see
[`vignette("architecture")`](https://chengmatt.github.io/SPoRC/dev/articles/architecture.md)
for why.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_NAA_state()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_NAA_state.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim_env <- Setup_sim_env(sim_list)
} # }
```
