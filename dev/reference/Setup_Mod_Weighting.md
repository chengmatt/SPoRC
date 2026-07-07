# Set likelihood and penalty weights for the estimation model

Assigns lambda (\\\lambda\\) multipliers to each likelihood component
and penalty term in the TMB/RTMB objective function. Weights scale the
relative contribution of each data source during estimation and can be
used to down-weight noisy data, implement iterative reweighting schemes
(e.g., Francis), or disable a component entirely by setting its weight
to `1`. Must be called after all data setup functions.

## Usage

``` r
Setup_Mod_Weighting(
  input_list,
  Wt_Catch = 1,
  Wt_FishIdx = 1,
  Wt_SrvIdx = 1,
  Wt_Catch_pop = 1,
  Wt_FishIdx_pop = 1,
  Wt_SrvIdx_pop = 1,
  Wt_Rec = 1,
  Wt_F = 1,
  Wt_Tagging = 1,
  Wt_FishAgeComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_FishLenComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_SrvLenComps = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_FishAgeComps_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_SrvAgeComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_FishLenComps_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_SrvLenComps_pop = array(1, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_srv_fleets)),
  Wt_Discard = 1,
  Wt_Discard_pop = 1,
  Wt_D = 1,
  Wt_FishAgeComps_discard = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_FishLenComps_discard = array(1, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_sexes,
    input_list$data$n_fish_fleets)),
  Wt_FishAgeComps_discard_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  Wt_FishLenComps_discard_pop = array(1, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  fish_sel_pen_wts = NULL,
  ret_sel_pen_wts = NULL,
  srv_sel_pen_wts = NULL
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions.

- Wt_Catch:

  Weight applied to the catch likelihood. Either a scalar applied
  uniformly across all fleets, regions, years, and seasons, or a numeric
  array `[n_regions × n_years × n_seas × n_fish_fleets]` for fleet- or
  time-specific weighting. Default `1`.

- Wt_FishIdx:

  Weight applied to the fishery index likelihood. Accepts the same
  scalar or array format as `Wt_Catch`, dimensioned
  `[n_regions × n_years × n_seas × n_fish_fleets]`. Default `1`.

- Wt_SrvIdx:

  Weight applied to the survey index likelihood. Accepts the same scalar
  or array format, dimensioned
  `[n_regions × n_years × n_seas × n_srv_fleets]`. Default `1`.

- Wt_Catch_pop:

  Weight applied to the population-specific catch likelihood. Either a
  scalar applied uniformly or a numeric array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`. Default `1`.

- Wt_FishIdx_pop:

  Weight applied to the population-specific fishery index likelihood.
  Same scalar or array format as `Wt_Catch_pop`, dimensioned
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`. Default `1`.

- Wt_SrvIdx_pop:

  Weight applied to the population-specific survey index likelihood.
  Same scalar or array format as `Wt_Catch_pop`, dimensioned
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`. Default `1`.

- Wt_Rec:

  Scalar weight applied to the recruitment deviation penalty
  (`ln_RecDevs`). Default `1`.

- Wt_F:

  Scalar weight applied to the fishing mortality deviation penalty
  (`ln_F_devs`). Default `1`.

- Wt_Tagging:

  Scalar weight applied to the tag-recovery likelihood. Default `1`.

- Wt_FishAgeComps:

  Weight applied to the fishery age composition likelihood. Either a
  scalar or a numeric array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. Default:
  array of `1`s.

- Wt_SrvAgeComps:

  Weight applied to the survey age composition likelihood. Either a
  scalar or a numeric array
  `[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`. Default:
  array of `1`s.

- Wt_FishLenComps:

  Weight applied to the fishery length composition likelihood. Same
  format as `Wt_FishAgeComps`,
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. Default:
  array of `1`s.

- Wt_SrvLenComps:

  Weight applied to the survey length composition likelihood. Same
  format as `Wt_SrvAgeComps`,
  `[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`. Default:
  array of `1`s.

- Wt_FishAgeComps_pop:

  Weight applied to the population-specific fishery age composition
  likelihood. Either a scalar or a numeric array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  Default: array of `1`s.

- Wt_SrvAgeComps_pop:

  Weight applied to the population-specific survey age composition
  likelihood. Either a scalar or a numeric array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`.
  Default: array of `1`s.

- Wt_FishLenComps_pop:

  Weight applied to the population-specific fishery length composition
  likelihood. Same format as `Wt_FishAgeComps_pop`,
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  Default: array of `1`s.

- Wt_SrvLenComps_pop:

  Weight applied to the population-specific survey length composition
  likelihood. Same format as `Wt_SrvAgeComps_pop`,
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`.
  Default: array of `1`s.

- Wt_Discard:

  Weight applied to the aggregated discard amount or fraction
  likelihood. Either a scalar applied uniformly or a numeric array
  `[n_regions × n_years × n_seas × n_fish_fleets]`. Default `1`.

- Wt_Discard_pop:

  Weight applied to the population-specific discard likelihood. Either a
  scalar or a numeric array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`. Default `1`.

- Wt_D:

  Scalar weight applied to the discard mortality rate deviation penalty
  (`logit_dmr_devs`). Default `1`.

- Wt_FishAgeComps_discard:

  Weight applied to the discard fishery age composition likelihood.
  Either a scalar or a numeric array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. Default:
  array of `1`s.

- Wt_FishLenComps_discard:

  Weight applied to the discard fishery length composition likelihood.
  Same format as `Wt_FishAgeComps_discard`,
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. Default:
  array of `1`s.

- Wt_FishAgeComps_discard_pop:

  Weight applied to the population-specific discard fishery age
  composition likelihood. Either a scalar or a numeric array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  Default: array of `1`s.

- Wt_FishLenComps_discard_pop:

  Weight applied to the population-specific discard fishery length
  composition likelihood. Same format as `Wt_FishAgeComps_discard_pop`,
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  Default: array of `1`s.

- fish_sel_pen_wts:

  `NULL` (default), or a named numeric vector/list with independent
  weights for any subset of six selectivity smoothness penalty terms
  (see
  [`resolve_sel_pen_wts`](https://chengmatt.github.io/SPoRC/dev/reference/resolve_sel_pen_wts.md)
  and
  [`Get_Selex_Smoothness_Penalty`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Selex_Smoothness_Penalty.md)),
  evaluated directly on the fleet's realized selectivity-at-bin-at-year
  surface and so applicable to any selectivity functional form:

  `"smooth_bin_curve"`

  :   Second-difference (curvature) penalty across bins.

  `"smooth_bin_diff"`

  :   Unconditional first-difference penalty across bins.

  `"smooth_yr_diff"`

  :   First-difference penalty across years.

  `"smooth_yr_curve"`

  :   Second-difference (curvature) penalty across years.

  `"smooth_dome"`

  :   Dome-shape (non-monotonicity) penalty across bins.

  `"smooth_mean_center"`

  :   Per-year mean-centering regularization.

  Any name not supplied defaults to `0` (off). Must be called after
  `Setup_Mod_Fishsel_and_Q`.

- ret_sel_pen_wts:

  Same format as `fish_sel_pen_wts`, for the retained fishery
  selectivity penalty.

- srv_sel_pen_wts:

  Same format as `fish_sel_pen_wts`, for the survey selectivity penalty.
  Must be called after `Setup_Mod_Srvsel_and_Q`.

## Value

The input `input_list` with all weight values stored in `$data` under
their respective names (`Wt_Catch`, `Wt_FishIdx`, `Wt_SrvIdx`,
`Wt_Catch_pop`, `Wt_FishIdx_pop`, `Wt_SrvIdx_pop`, `Wt_FishAgeComps`,
`Wt_SrvAgeComps`, `Wt_FishLenComps`, `Wt_SrvLenComps`,
`Wt_FishAgeComps_pop`, `Wt_SrvAgeComps_pop`, `Wt_FishLenComps_pop`,
`Wt_SrvLenComps_pop`, `Wt_Rec`, `Wt_F`, `Wt_Tagging`, `Wt_Discard`,
`Wt_Discard_pop`, `Wt_D`, `Wt_FishAgeComps_discard`,
`Wt_FishLenComps_discard`, `Wt_FishAgeComps_discard_pop`,
`Wt_FishLenComps_discard_pop`).

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_FishIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md)
