# Set up discards, fishery index, age composition, and length composition inputs

Populates `input_list` with observed fishery indices, age compositions,
and length compositions (both pooled and population-specific) along with
their usage indicators, likelihood types, composition structure types,
input sample sizes, and overdispersion and correlation parameter
starting values and mappings. Must be called after
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

## Usage

``` r
Setup_Mod_FishIdx_and_Comps(
  input_list,
  ObsFishIdx,
  ObsFishIdx_SE,
  ObsFishIdx_pop = NULL,
  ObsFishIdx_pop_SE = NULL,
  UseFishIdx_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  fish_idx_type,
  t_fish = array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  UseFishIdx,
  ObsFishAgeComps,
  UseFishAgeComps,
  ISS_FishAgeComps = NULL,
  ObsFishLenComps,
  UseFishLenComps,
  ISS_FishLenComps = NULL,
  FishAgeComps_LikeType,
  FishLenComps_LikeType,
  FishAgeComps_Type,
  FishLenComps_Type,
  ObsFishAgeComps_pop = NULL,
  UseFishAgeComps_pop = array(0, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ISS_FishAgeComps_pop = NULL,
  ObsFishLenComps_pop = NULL,
  UseFishLenComps_pop = array(0, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ISS_FishLenComps_pop = NULL,
  FishAgeComps_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
  FishLenComps_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
  FishAgeComps_pop_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_fish_fleets, sep = ""),
  FishLenComps_pop_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_fish_fleets, sep = ""),
  fish_idx_ages = NULL,
  FishAgeComps_bins = NULL,
  FishIdx_LikeType = rep("lognormal", input_list$data$n_fish_fleets),
  FishIdx_Cov = NULL,
  ObsFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages),
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  UseFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ISS_FishAgeComps_discard = NULL,
  ObsFishLenComps_discard = array(0, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, length(input_list$data$lens),
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  UseFishLenComps_discard = array(0, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ISS_FishLenComps_discard = NULL,
  FishAgeComps_discard_LikeType = rep("none", input_list$data$n_fish_fleets),
  FishLenComps_discard_LikeType = rep("none", input_list$data$n_fish_fleets),
  FishAgeComps_discard_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_fish_fleets, sep = ""),
  FishLenComps_discard_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_fish_fleets, sep = ""),
  ObsFishAgeComps_discard_pop = NULL,
  UseFishAgeComps_discard_pop = array(0, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ISS_FishAgeComps_discard_pop = NULL,
  ObsFishLenComps_discard_pop = NULL,
  UseFishLenComps_discard_pop = array(0, dim = c(input_list$data$n_pop,
    input_list$data$n_regions, length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ISS_FishLenComps_discard_pop = NULL,
  FishAgeComps_discard_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
  FishLenComps_discard_pop_LikeType = rep("none", input_list$data$n_fish_fleets),
  FishAgeComps_discard_pop_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_fish_fleets, sep = ""),
  FishLenComps_discard_pop_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_fish_fleets, sep = ""),
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions.

- ObsFishIdx:

  Observed fishery CPUE or biomass index array
  `[n_regions × n_years × n_seas × n_fish_fleets]`.

- ObsFishIdx_SE:

  Standard errors of `ObsFishIdx` on the log scale, same dimensions as
  `ObsFishIdx`.

- ObsFishIdx_pop:

  Observed population-specific fishery index array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.

- ObsFishIdx_pop_SE:

  Lognormal standard errors for `ObsFishIdx_pop`, same dimensions
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.

- UseFishIdx_pop:

  Binary indicator array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`. `1` =
  include population-specific index in likelihood; `0` = exclude.
  Default: all zeros.

- fish_idx_type:

  Character vector of length `n_fish_fleets` specifying the index type
  for each fleet. `"biom"` = biomass; `"abd"` = abundance; `"none"` = no
  index for this fleet.

- t_fish:

  Array `[n_regions x n_seas x n_fish_fleets]` giving the fishery index
  timing: the fraction of the season elapsed when each index is
  observed. Numbers at age are decayed by `exp(-t_fish * ZAA)` before
  the index is formed, the same convention `t_srv` uses for surveys.
  Defaults to `0` (start of season), which is what the model did before
  this argument existed; set `0.5` for a mid-season index.

- UseFishIdx:

  Binary indicator array
  `[n_regions × n_years × n_seas × n_fish_fleets]`. `1` = include index
  in the likelihood; `0` = exclude.

- ObsFishAgeComps:

  Observed fishery age composition array
  `[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`.
  Values may be raw counts or proportions; if proportions, supply
  `ISS_FishAgeComps` explicitly.

- UseFishAgeComps:

  Binary indicator array
  `[n_regions × n_years × n_seas × n_fish_fleets]`. `1` = fit age
  compositions; `0` = exclude.

- ISS_FishAgeComps:

  Input sample size array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If `NULL`
  (default), computed automatically by summing `ObsFishAgeComps` within
  each year-fleet-season-region cell according to `FishAgeComps_Type`.

- ObsFishLenComps:

  Observed fishery length composition array
  `[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.
  Only required when `input_list$data$fit_lengths == 1`.

- UseFishLenComps:

  Binary indicator array
  `[n_regions × n_years × n_seas × n_fish_fleets]`. `1` = fit length
  compositions; `0` = exclude.

- ISS_FishLenComps:

  Input sample size array for length compositions
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If `NULL`
  (default), derived automatically from `ObsFishLenComps`.

- FishAgeComps_LikeType:

  Character vector of length `n_fish_fleets` specifying the likelihood
  for fishery age compositions. Options: `"Multinomial"`,
  `"Dirichlet-Multinomial"`, `"iid-Logistic-Normal"`,
  `"1d-Logistic-Normal"`, `"2d-Logistic-Normal"`, `"none"`.

- FishLenComps_LikeType:

  Same as `FishAgeComps_LikeType` but for length compositions.

- FishAgeComps_Type:

  Character vector defining the age composition structure (aggregation
  level) for each fleet and time period. Each element must follow the
  format `"<type>_Year_<start>-<end>_Fleet_<f>"` or
  `"<type>_Year_<start>-terminal_Fleet_<f>"`. Valid types:

  `"agg"`

  :   Aggregated across regions and sexes (incompatible with
      `"2d-Logistic-Normal"`).

  `"spltRspltS"`

  :   Split by region and sex.

  `"spltRjntS"`

  :   Split by region, summed jointly across sexes.

  `"none"`

  :   No composition data for this fleet and period.

  Example:
  `c("spltRjntS_Year_1-10_Fleet_1", "agg_Year_11-terminal_Fleet_1")`.

- FishLenComps_Type:

  Same format and options as `FishAgeComps_Type` but applied to length
  compositions.

- ObsFishAgeComps_pop:

  Observed population-specific fishery age composition array
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`.
  Required when any element of `UseFishAgeComps_pop` is `1`.

- UseFishAgeComps_pop:

  Binary indicator array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`. `1` = fit
  population-specific age compositions; `0` = exclude. Default: all
  zeros.

- ISS_FishAgeComps_pop:

  Input sample size array for population-specific age compositions
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If
  `NULL` (default), computed automatically by summing
  `ObsFishAgeComps_pop` within each population-year-fleet-season-region
  cell according to `FishAgeComps_pop_Type`.

- ObsFishLenComps_pop:

  Observed population-specific fishery length composition array
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.
  Required when `input_list$data$fit_lengths == 1` and any element of
  `UseFishLenComps_pop` is `1`.

- UseFishLenComps_pop:

  Binary indicator array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`. `1` = fit
  population-specific length compositions; `0` = exclude. Default: all
  zeros.

- ISS_FishLenComps_pop:

  Input sample size array for population-specific length compositions
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If
  `NULL` (default), derived automatically from `ObsFishLenComps_pop`.

- FishAgeComps_pop_LikeType:

  Character vector of length `n_fish_fleets` specifying the likelihood
  for population-specific fishery age compositions. Same options as
  `FishAgeComps_LikeType`. Default: `"none"` for all fleets.

- FishLenComps_pop_LikeType:

  Character vector of length `n_fish_fleets` specifying the likelihood
  for population-specific fishery length compositions. Same options as
  `FishLenComps_LikeType`. Default: `"none"` for all fleets.

- FishAgeComps_pop_Type:

  Character vector defining the composition structure for
  population-specific age compositions. Same format and options as
  `FishAgeComps_Type`. Default: `"none"` for all fleets across all
  years.

- FishLenComps_pop_Type:

  Character vector defining the composition structure for
  population-specific length compositions. Same format and options as
  `FishLenComps_Type`. Default: `"none"` for all fleets across all
  years.

- fish_idx_ages:

  Per-fleet selection of which ages contribute to the index total.
  Either a list with one element per fishery fleet, where each element
  is a vector of ages or `NULL` for all ages, or an array
  `[n_ages x n_fish_fleets]` of 0/1 weights. Default `NULL` uses every
  age for every fleet. The fleet's compositions are unaffected.

- FishAgeComps_bins:

  Which age bins each fishery fleet's age composition is fitted over.
  Supply a list with one element per fleet, each a vector of age indices
  or `NULL` for all ages, or an `[n_ages x n_fish_fleets]` array of 0/1
  weights. Both observed and expected compositions are restricted to the
  named bins and renormalized within them, so excluded bins are left out
  of the likelihood rather than being forced to be explained; this is
  how a fleet that only ages part of its age range is fitted. Indices
  refer to observed bins, that is after any ageing error has mapped
  model ages onto observed ones. Every fleet must retain at least one
  bin. Default `NULL`, which fits all ages for all fleets.

- FishIdx_LikeType:

  Character vector `[n_fish_fleets]` giving the error structure of each
  fishery index. Options are `"lognormal"` (default, the observation
  standard errors are on the log scale), `"normal"` (arithmetic scale),
  and `"mvn"` (multivariate normal on the arithmetic scale using a fixed
  covariance supplied through `FishIdx_Cov`). One-step-ahead residuals
  are available only for lognormal fleets. A fleet's population-specific
  index stream follows the same choice for `"lognormal"` and `"normal"`,
  but stays lognormal under `"mvn"`, whose covariance describes the
  regional series only.

- FishIdx_Cov:

  List with one element per fishery fleet holding the fixed covariance
  matrix for fleets using `"mvn"`, and `NULL` otherwise. Each matrix
  must be square with one row per observation the fleet fits, ordered as
  the observations appear when scanning that fleet's `UseFishIdx` slice
  in array order.

- ObsFishAgeComps_discard:

  Observed fishery age composition from discards
  `[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`.
  Structure must match `ObsFishAgeComps`.

- UseFishAgeComps_discard:

  Binary indicator array for discard age compositions
  `[n_regions × n_years × n_seas × n_fish_fleets]`. `1` = include
  discard age compositions in likelihood; `0` = exclude.

- ISS_FishAgeComps_discard:

  Input sample size array for discard age compositions
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If `NULL`,
  derived automatically from `ObsFishAgeComps_discard` using
  `FishAgeComps_discard_Type`.

- ObsFishLenComps_discard:

  Observed fishery length composition from discards
  `[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.
  Required if `input_list$data$fit_lengths == 1`.

- UseFishLenComps_discard:

  Binary indicator array for discard length compositions
  `[n_regions × n_years × n_seas × n_fish_fleets]`. `1` = include
  discard length compositions in likelihood; `0` = exclude.

- ISS_FishLenComps_discard:

  Input sample size array for discard length compositions
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If `NULL`,
  derived automatically from `ObsFishLenComps_discard`.

- FishAgeComps_discard_LikeType:

  Character vector of length `n_fish_fleets` specifying likelihood type
  for discard age compositions. Options:

  `"Multinomial"`

  :   Standard multinomial likelihood

  `"Dirichlet-Multinomial"`

  :   Overdispersed multinomial

  `"iid-Logistic-Normal"`

  :   Independent logistic-normal

  `"1d-Logistic-Normal"`

  :   1D correlated logistic-normal

  `"2d-Logistic-Normal"`

  :   2D correlated logistic-normal

  `"none"`

  :   No discard age composition likelihood

- FishLenComps_discard_LikeType:

  Same specification as `FishAgeComps_discard_LikeType`, but for discard
  length compositions.

- FishAgeComps_discard_Type:

  Character vector defining discard age composition structure by fleet
  and year block. Format: `"<type>_Year_<start>-<end>_Fleet_<f>"` or
  `"<type>_Year_<start>-terminal_Fleet_<f>"`. Valid types:

  `"agg"`

  :   Aggregated across regions and sexes

  `"spltRspltS"`

  :   Split by region and sex

  `"spltRjntS"`

  :   Split by region, joint across sexes

  `"none"`

  :   No discard age composition

- FishLenComps_discard_Type:

  Same format and options as `FishAgeComps_discard_Type`, applied to
  discard length compositions.

- ObsFishAgeComps_discard_pop:

  Observed population-specific discard age composition array
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`.

- UseFishAgeComps_discard_pop:

  Binary indicator array for population-specific discard age
  compositions `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.

- ISS_FishAgeComps_discard_pop:

  Input sample size array for population-specific discard age
  compositions
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If
  `NULL`, computed from `ObsFishAgeComps_discard_pop`.

- ObsFishLenComps_discard_pop:

  Observed population-specific discard length composition array
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.

- UseFishLenComps_discard_pop:

  Binary indicator array for population-specific discard length
  compositions `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.

- ISS_FishLenComps_discard_pop:

  Input sample size array for population-specific discard length
  compositions
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. If
  `NULL`, derived from `ObsFishLenComps_discard_pop`.

- FishAgeComps_discard_pop_LikeType:

  Character vector of length `n_fish_fleets` specifying likelihood type
  for population-specific discard age compositions. Same options as
  `FishAgeComps_discard_LikeType`.

- FishLenComps_discard_pop_LikeType:

  Same as above but for discard length compositions.

- FishAgeComps_discard_pop_Type:

  Character vector defining structure for population-specific discard
  age compositions. Same format as `FishAgeComps_discard_Type`.

- FishLenComps_discard_pop_Type:

  Character vector defining structure for population-specific discard
  length compositions. Same format as `FishLenComps_discard_Type`.

- ...:

  Optional starting value overrides for overdispersion and correlation
  parameters.

## Value

The input `input_list` with `$data`, `$par`, and `$map` updated with all
fishery index and composition fields, including pooled and
population-specific observed arrays, computed or supplied ISS arrays,
integer-coded likelihood and composition type matrices, overdispersion
parameters, and their factor maps.

## Details

When `ISS_FishAgeComps`, `ISS_FishLenComps`, `ISS_FishAgeComps_pop`, or
`ISS_FishLenComps_pop` are `NULL`, input sample sizes are derived
automatically by summing the observed composition arrays within each
year-fleet-season-region cell, consistent with the specified composition
type.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Catch_and_F()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
[`Setup_Mod_Discard_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Discard_Comps.md),
[`Setup_Mod_Fishsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Fishsel_and_Q.md),
[`Setup_Mod_Movement()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Movement.md),
[`Setup_Mod_Rec()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md),
[`Setup_Mod_Retsel()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Retsel.md),
[`Setup_Mod_SrvIdx_and_Comps()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md),
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
