# Simulate a static (open-loop) spatial age- and sex-structured population

Runs a complete multi-replicate operating model simulation with no
feedback between the population and the harvest control rule (i.e.,
fishing mortality is fixed as supplied in `sim_list`). Calls
[`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
to create an isolated execution environment and then iterates
[`run_annual_cycle`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md)
over all years and simulation replicates. All simulation outputs are
collected from the environment and returned as a named list. Optionally
writes the output to an RDS file.

## Usage

``` r
Simulate_Pop_Static(sim_list, output_path = NULL)
```

## Arguments

- sim_list:

  Simulation list returned by the last upstream setup function
  (typically
  [`Setup_Sim_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md)
  or
  [`Setup_Sim_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md)).

- output_path:

  Character string. File path for saving the output list as an RDS file
  via `saveRDS`. If `NULL` (default), no file is written.

## Value

A named list containing all simulation outputs, including (among
others): `NAA`, `NAA0`, `SSB`, `Dynamic_SSB0`, `eff_SSB`, `Rec`,
`ln_RecDevs`, `ln_InitDevs`, `ZAA`, `TrueCatch`, `ObsCatch`,
`TrueCatch_pop`, `ObsCatch_pop`, `CAA`, `CAL`, `ObsFishAgeComps`,
`ObsFishAgeComps_pop`, `ObsFishLenComps`, `ObsFishLenComps_pop`,
`ObsFishIdx`, `TrueFishIdx`, `ObsFishIdx_pop`, `TrueFishIdx_pop`,
`SrvIAA`, `SrvIAL`, `ObsSrvAgeComps`, `ObsSrvAgeComps_pop`,
`ObsSrvLenComps`, `ObsSrvLenComps_pop`, `ObsSrvIdx`, `TrueSrvIdx`,
`ObsSrvIdx_pop`, `TrueSrvIdx_pop`, `conv_tagged_fish`,
`conv_tagged_fish_attr`, `conv_tag_fish_avail`,
`pred_conv_tag_fish_recap`, `obs_conv_tag_fish_recap`, and key dimension
scalars (`n_regions`, `n_pop`, `n_yrs`, `n_ages`, etc.). Note that
`n_years` and `n_yrs` are both present for backwards compatibility.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/dev/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
