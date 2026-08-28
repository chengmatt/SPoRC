# Initialize model dimension settings

Creates the foundational `input_list` object used throughout the
estimation model setup. All downstream configuration functions
([`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
etc.) expect an `input_list` produced by this function. Dimension
vectors and scalars are stored in `$data`, with empty `$par` and `$map`
sublists ready to be populated by subsequent setup calls. Optionally, a
`$config` sublist can be stored to retain model configuration options.

## Usage

``` r
Setup_Mod_Dim(
  years,
  ages,
  lens,
  n_pop = 1,
  natal_region = NULL,
  n_seas = 1,
  seasdur = if (n_seas == 1) 1 else rep(1/n_seas, n_seas),
  n_regions,
  n_sexes,
  n_fish_fleets,
  n_srv_fleets,
  n_proj_yrs_devs = 0,
  do_internal_comp_osa = FALSE,
  do_internal_conv_tag_osa = FALSE,
  verbose = FALSE,
  store_config = FALSE
)
```

## Arguments

- years:

  Numeric vector of calendar years included in the model (e.g.,
  `1990:2024`). The length of this vector determines `n_years`
  throughout the model.

- ages:

  Numeric vector of modeled age classes (e.g., `2:31` for a model
  spanning ages 2-31). The final element is treated as a plus-group.

- lens:

  Numeric vector of length bin midpoints. Set to `NULL` when length data
  are not modeled; a scalar placeholder of `1` is stored internally in
  that case.

- n_pop:

  Positive integer. Number of distinct populations (default `1`).
  Populations can have independent stock-recruit relationships and natal
  regions but share the spatial domain defined by `n_regions`.

- natal_region:

  Integer vector of length `n_pop` mapping each population to its natal
  region (1-indexed). Defaults are applied when `NULL`:

  - `n_regions == 1`: all populations assigned to region 1.

  - `n_pop == n_regions`: one-to-one mapping (`1:n_pop`).

  - `n_pop == 1`: assigned to region 1.

  Must be supplied explicitly when `n_pop > 1`, `n_pop != n_regions`,
  and `n_regions > 1`, as no sensible default exists when populations
  must share a natal region.

- n_seas:

  Positive integer. Number of seasons within each year (default `1`).

- seasdur:

  Numeric vector of length `n_seas` giving the duration of each season
  as a fraction of a year (values should sum to 1). Defaults to `1` for
  a single season or `rep(1 / n_seas, n_seas)` for equal-length seasons
  when `n_seas > 1`.

- n_regions:

  Positive integer. Number of spatial regions.

- n_sexes:

  Integer. Number of sexes; must be `1` (sex-aggregated) or `2`
  (sex-structured).

- n_fish_fleets:

  Positive integer. Number of fishery fleets.

- n_srv_fleets:

  Positive integer. Number of survey fleets.

- n_proj_yrs_devs:

  Non-negative integer. Number of projection years for which deviation
  parameters (`ln_RecDevs`, `move_devs`, `ln_fishsel_devs`,
  `ln_srvsel_devs`) are allocated. Set to `0` (default) when projections
  beyond the assessment period are not required.

- do_internal_comp_osa:

  Logical. If `TRUE`, allows OSA residuals for composition datasets.
  Default `FALSE`.

- do_internal_conv_tag_osa:

  Logical. If `TRUE`, allows OSA residuals for tagging datasets. Default
  `FALSE`.

- verbose:

  Logical. If `TRUE`, prints a summary of all dimension settings to the
  console via [`message()`](https://rdrr.io/r/base/message.html) after
  setup. Default `FALSE`.

- store_config:

  Logical. If `TRUE`, stores configuration of the model. Default is
  `FALSE`

## Value

A named list (`input_list`) with three sublists:

- `$data`:

  Dimension scalars and vectors stored for use by the TMB/RTMB objective
  function: `years`, `ages`, `lens`, `n_pop`, `n_regions`,
  `natal_region`, `n_seas`, `seasdur`, `n_sexes`, `n_fish_fleets`,
  `n_srv_fleets`, `n_proj_yrs_devs`.

- `$par`:

  Empty list to be populated with parameter starting values by
  downstream setup functions.

- `$map`:

  Empty list to be populated with TMB/RTMB factor maps by downstream
  setup functions.

The top-level `$verbose` flag is also set and respected by all
subsequent setup functions.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
