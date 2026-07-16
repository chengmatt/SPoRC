# Set up observed survey indices and composition data

Ingests observed survey index, age composition, and length composition
data (both pooled and population-specific) into `input_list$data`,
initialises overdispersion and correlation starting values in
`input_list$par`, and constructs parameter maps via
[`do_SrvAge_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_theta_mapping.md),
[`do_SrvLen_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvLen_theta_mapping.md),
[`do_SrvAge_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvAge_corr_pars_mapping.md),
and
[`do_SrvLen_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_SrvLen_corr_pars_mapping.md).
When `ISS_SrvAgeComps`, `ISS_SrvLenComps`, `ISS_SrvAgeComps_pop`, or
`ISS_SrvLenComps_pop` is `NULL`, input sample sizes are derived
automatically by summing observed composition counts across the
appropriate dimensions each year. Must be called after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
and before model compilation.

## Usage

``` r
Setup_Mod_SrvIdx_and_Comps(
  input_list,
  ObsSrvIdx,
  ObsSrvIdx_SE,
  UseSrvIdx,
  ObsSrvIdx_pop = NULL,
  ObsSrvIdx_pop_SE = NULL,
  UseSrvIdx_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
  srv_idx_type,
  ObsSrvAgeComps,
  UseSrvAgeComps,
  ObsSrvLenComps,
  UseSrvLenComps,
  ISS_SrvAgeComps = NULL,
  ISS_SrvLenComps = NULL,
  SrvAgeComps_LikeType,
  SrvLenComps_LikeType,
  SrvAgeComps_Type,
  SrvLenComps_Type,
  ObsSrvAgeComps_pop = NULL,
  UseSrvAgeComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
  ISS_SrvAgeComps_pop = NULL,
  ObsSrvLenComps_pop = NULL,
  UseSrvLenComps_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, input_list$data$n_srv_fleets)),
  ISS_SrvLenComps_pop = NULL,
  SrvAgeComps_pop_LikeType = rep("none", input_list$data$n_srv_fleets),
  SrvLenComps_pop_LikeType = rep("none", input_list$data$n_srv_fleets),
  SrvAgeComps_pop_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_srv_fleets, sep = ""),
  SrvLenComps_pop_Type = paste("none_Year_1-terminal_Fleet_",
    1:input_list$data$n_srv_fleets, sep = ""),
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions.

- ObsSrvIdx:

  Observed survey index array
  `[n_regions × n_years × n_seas × n_srv_fleets]`.

- ObsSrvIdx_SE:

  Lognormal standard errors for `ObsSrvIdx`, same dimensions
  `[n_regions × n_years × n_seas × n_srv_fleets]`.

- UseSrvIdx:

  Binary indicator array
  `[n_regions × n_years × n_seas × n_srv_fleets]`. `1` = include in
  likelihood; `0` = exclude.

- ObsSrvIdx_pop:

  Observed population-specific survey index array
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`.

- ObsSrvIdx_pop_SE:

  Lognormal standard errors for `ObsSrvIdx_pop`, same dimensions
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`.

- UseSrvIdx_pop:

  Binary indicator array
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`. `1` = include
  population-specific index in likelihood; `0` = exclude. Default: all
  zeros.

- srv_idx_type:

  Character vector `[n_srv_fleets]` specifying the index type per fleet.
  One of `"biom"` (biomass), `"abd"` (abundance), or `"none"` (no index
  for that fleet). Converted to integer codes (`1`, `0`, `999`) before
  storage.

- ObsSrvAgeComps:

  Observed survey age compositions, array
  `[n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]`.
  Values may be counts or proportions on a comparable scale.

- UseSrvAgeComps:

  Binary indicator array
  `[n_regions × n_years × n_seas × n_srv_fleets]`. `1` = fit age
  compositions; `0` = exclude.

- ObsSrvLenComps:

  Observed survey length compositions, array
  `[n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]`.
  Only validated when `input_list$data$fit_lengths = 1` in `$data`.

- UseSrvLenComps:

  Binary indicator array
  `[n_regions × n_years × n_seas × n_srv_fleets]`. `1` = fit length
  compositions; `0` = exclude.

- ISS_SrvAgeComps:

  Input sample sizes for survey age compositions, array
  `[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`, or `NULL`
  to derive automatically by summing `ObsSrvAgeComps` across the age
  dimension each year, respecting `SrvAgeComps_Type`.

- ISS_SrvLenComps:

  Input sample sizes for survey length compositions, same structure as
  `ISS_SrvAgeComps`, or `NULL` for automatic derivation from
  `ObsSrvLenComps`.

- SrvAgeComps_LikeType:

  Character vector `[n_srv_fleets]` specifying the likelihood for survey
  age compositions. One of `"none"`, `"Multinomial"`,
  `"Dirichlet-Multinomial"`, `"iid-Logistic-Normal"`,
  `"1d-Logistic-Normal"`, `"2d-Logistic-Normal"`. Converted to integer
  codes (`999`, `0`–`4`) before storage.

- SrvLenComps_LikeType:

  Character vector `[n_srv_fleets]` specifying the likelihood for survey
  length compositions. Same options as `SrvAgeComps_LikeType`.

- SrvAgeComps_Type:

  Character vector defining the survey age composition structure per
  fleet and year range. Each element follows the format
  `"<type>_Year_<start>-<end>_Fleet_<fleet>"`. Use `"terminal"` in place
  of the end year to extend to the final model year. Valid types:

  `"agg"`

  :   Aggregated across regions and sexes. Not compatible with
      `"2d-Logistic-Normal"`.

  `"spltRspltS"`

  :   Split by region and sex.

  `"spltRjntS"`

  :   Split by region, joint across sexes.

  `"none"`

  :   No composition data used.

  Parsed into a `[n_years × n_srv_fleets]` integer matrix before
  storage. An error is raised if any cell remains `NA` after parsing,
  indicating an incomplete year range specification.

- SrvLenComps_Type:

  Character vector defining the survey length composition structure.
  Same format and options as `SrvAgeComps_Type`.

- ObsSrvAgeComps_pop:

  Observed population-specific survey age composition array
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]`.
  Required when any element of `UseSrvAgeComps_pop` is `1`.

- UseSrvAgeComps_pop:

  Binary indicator array
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`. `1` = fit
  population-specific age compositions; `0` = exclude. Default: all
  zeros.

- ISS_SrvAgeComps_pop:

  Input sample size array for population-specific survey age
  compositions
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`. If
  `NULL` (default), computed automatically by summing
  `ObsSrvAgeComps_pop` within each population–year–fleet–season–region
  cell according to `SrvAgeComps_pop_Type`.

- ObsSrvLenComps_pop:

  Observed population-specific survey length composition array
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]`.
  Required when `input_list$data$fit_lengths == 1` and any element of
  `UseSrvLenComps_pop` is `1`.

- UseSrvLenComps_pop:

  Binary indicator array
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`. `1` = fit
  population-specific length compositions; `0` = exclude. Default: all
  zeros.

- ISS_SrvLenComps_pop:

  Input sample size array for population-specific survey length
  compositions
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`. If
  `NULL` (default), derived automatically from `ObsSrvLenComps_pop`.

- SrvAgeComps_pop_LikeType:

  Character vector of length `n_srv_fleets` specifying the likelihood
  for population-specific survey age compositions. Same options as
  `SrvAgeComps_LikeType`. Default: `"none"` for all fleets.

- SrvLenComps_pop_LikeType:

  Character vector of length `n_srv_fleets` specifying the likelihood
  for population-specific survey length compositions. Same options as
  `SrvLenComps_LikeType`. Default: `"none"` for all fleets.

- SrvAgeComps_pop_Type:

  Character vector defining the composition structure for
  population-specific survey age compositions. Same format and options
  as `SrvAgeComps_Type`. Default: `"none"` for all fleets across all
  years.

- SrvLenComps_pop_Type:

  Character vector defining the composition structure for
  population-specific survey length compositions. Same format and
  options as `SrvLenComps_Type`. Default: `"none"` for all fleets across
  all years.

- ...:

  Optional named starting values for overdispersion and correlation
  parameters.

## Value

The input `input_list` with survey data stored in `$data` (`ObsSrvIdx`,
`ObsSrvIdx_SE`, `UseSrvIdx`, `ObsSrvIdx_pop`, `ObsSrvIdx_pop_SE`,
`UseSrvIdx_pop`, `ObsSrvAgeComps`, `UseSrvAgeComps`, `ISS_SrvAgeComps`,
`ObsSrvLenComps`, `UseSrvLenComps`, `ISS_SrvLenComps`,
`ObsSrvAgeComps_pop`, `UseSrvAgeComps_pop`, `ISS_SrvAgeComps_pop`,
`ObsSrvLenComps_pop`, `UseSrvLenComps_pop`, `ISS_SrvLenComps_pop`,
`SrvAgeComps_LikeType`, `SrvLenComps_LikeType`,
`SrvAgeComps_pop_LikeType`, `SrvLenComps_pop_LikeType`,
`SrvAgeComps_Type`, `SrvLenComps_Type`, `SrvAgeComps_pop_Type`,
`SrvLenComps_pop_Type`, `srv_idx_type`); overdispersion and correlation
starting values in `$par`; and factor maps in `$map` for all pooled and
population-specific overdispersion and correlation parameter arrays.

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
[`Setup_Mod_Srvsel_and_Q()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md),
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
