# Extract simulation data into SPoRC format

This function subsets and reshapes biological, tagging, fishery, and
survey data from a simulation environment for use in SPoRC analyses.

## Usage

``` r
simulation_data_to_SPoRC(sim_env, y, sim)
```

## Arguments

- sim_env:

  A simulation environment / object (list or environment) containing
  arrays of biological quantities, tagging information, fishery data,
  and survey data.

- y:

  Integer. Number of years to retain (subset from \`1:y\`).

- sim:

  Integer. Simulation replicate index to extract.

## Value

A named list with the following elements:

- WAA:

  Weight-at-age array \[region × year × age × sex\].

- MatAA:

  Maturity-at-age array \[region × year × age × sex\].

- SizeAgeTrans:

  Size–age transition array \[region × year × length × age × sex\].

- AgeingError:

  Ageing error matrix \[year × age × error × sim\].

- tag_release_indicator:

  Tag release indicators (or \`NULL\` if tagging not used).

- Obs_Tag_Recap:

  Observed tag recapture array (or \`NULL\`).

- Tagged_Fish:

  Tagged fish counts (or \`NULL\`).

- n_tag_cohorts:

  Number of tag release cohorts (or \`NULL\`).

- ObsCatch:

  Observed fishery catch array \[region × year × fleet\].

- ln_sigmaC:

  Log Fishery Catch SD \[region × year × fleet\].

- UseCatch:

  Binary indicator array for catch data availability.

- ObsFishIdx:

  Observed fishery index array \[region × year × fleet\].

- ObsFishIdx_SE:

  Standard error for fishery index array.

- UseFishIdx:

  Binary indicator array for fishery indices.

- ObsFishAgeComps:

  Observed fishery age composition array.

- ObsFishLenComps:

  Observed fishery length composition array.

- ISS_FishAgeComps:

  Implied sample sizes for fishery age compositions.

- ISS_FishLenComps:

  Implied sample sizes for fishery length compositions.

- UseFishAgeComps:

  Binary indicator array for fishery age comps.

- UseFishLenComps:

  Binary indicator array for fishery length comps.

- ObsSrvIdx:

  Observed survey index array \[region × year × fleet\].

- ObsSrvIdx_SE:

  Standard error for survey index array.

- UseSrvIdx:

  Binary indicator array for survey indices.

- ObsSrvAgeComps:

  Observed survey age composition array.

- ObsSrvLenComps:

  Observed survey length composition array.

- ISS_SrvAgeComps:

  Implied sample sizes for survey age compositions.

- ISS_SrvLenComps:

  Implied sample sizes for survey length compositions.

- UseSrvAgeComps:

  Binary indicator array for survey age comps.

- UseSrvLenComps:

  Binary indicator array for survey length comps.

## See also

Other Simulation Setup:
[`Setup_Sim_Biologicals()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Biologicals.md),
[`Setup_Sim_Containers()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Containers.md),
[`Setup_Sim_Dim()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Dim.md),
[`Setup_Sim_Fishing()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Fishing.md),
[`Setup_Sim_Rec()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Rec.md),
[`Setup_Sim_Survey()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Survey.md),
[`Setup_Sim_Tagging()`](https://chengmatt.github.io/SPoRC/reference/Setup_Sim_Tagging.md),
[`Setup_sim_env()`](https://chengmatt.github.io/SPoRC/reference/Setup_sim_env.md),
[`Simulate_Pop_Static()`](https://chengmatt.github.io/SPoRC/reference/Simulate_Pop_Static.md),
[`run_annual_cycle()`](https://chengmatt.github.io/SPoRC/reference/run_annual_cycle.md),
[`simulation_self_test()`](https://chengmatt.github.io/SPoRC/reference/simulation_self_test.md)
