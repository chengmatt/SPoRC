# Extract simulation outputs into SPoRC estimation model format

Subsets and reshapes the biological, tagging, fishery and survey arrays
of a simulation environment to years `1:y` and replicate `sim`, giving a
list ready for
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md)
and
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md).
The `Use*` flags are derived from the extracted observations, 1 where a
value is present and positive.

## Usage

``` r
simulation_data_to_SPoRC(sim_env, y, sim)
```

## Arguments

- sim_env:

  Simulation environment or list, from
  [`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md)
  or
  [`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md),
  holding the operating model arrays.

- y:

  Integer. Last year to include; years `1:y` are kept.

- sim:

  Integer. Simulation replicate to extract.

## Value

Named list, every array with `y` in its year dim. The biological
elements are `WAA`
`[n_pop x n_regions x y x n_seas x n_ages x n_sexes]`, `WAA_fish` and
`WAA_srv` with a trailing fleet dim, `MatAA`, `SizeAgeTrans`,
`AgeingError` `[y x n_ages x n_obs_ages]`, whose columns are the
observed ages the model ages are read onto, and `AgeingError_fish` and
`AgeingError_srv`, `NULL` when the fleets share one matrix. `LenBinMap`
is the map onto the length bins the length compositions are recorded on,
`NULL` when those are the model's bins, and goes to
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md)
with the rest, as does `CAAL_LenBinMap`, the model bins each
age-at-length row covers. The tagging elements are
`use_conv_fish_tagging`, `conv_tag_release_indicator`,
`obs_conv_tag_fish_recap`, `conv_tagged_fish`, `conv_tagged_fish_attr`
and `n_tag_cohorts`.

Each fishery and survey data source contributes its observations, its
observation error or input sample sizes, and its use flag: `ObsCatch`
with `ln_sigmaC` and `UseCatch`, `ObsDiscard` with `ln_sigmaD` and
`UseDiscard`, `ObsFishIdx` and `ObsSrvIdx` with their `_SE` and use
arrays, and the four composition streams (fishery and survey age and
length) with their `ISS_*` and `Use*` arrays. Each has a `_pop`
counterpart, and the fishery compositions also have `_discard` and
`_discard_pop` ones. The length composition elements and their input
sample sizes are `NULL` when no size-age transition matrix is present.

The at-age data sources, `ObsCatchAA`, `ObsDiscardAA` and `ObsSrvIdxAA`
and their `_pop` counterparts, come with their `_SE` arrays and the
operating model's own use flags, since a normal likelihood can draw a
value at or below zero. Each is `NULL` when the operating model draws
none of it.

## Details

Population-specific arrays are extracted when `sim_env` holds them, with
their own flags derived the same way. The length composition outputs and
`SizeAgeTrans` are `NULL` when the environment holds no size-age
transition matrix, and the tagging outputs are `NULL` under
`use_conv_fish_tagging = 0`; otherwise only cohorts released in `1:y`
are kept.

## See also

[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Tagging`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Simulate_Pop_Static`](https://chengmatt.github.io/SPoRC/dev/reference/Simulate_Pop_Static.md),
[`Setup_sim_env`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_sim_env.md)
