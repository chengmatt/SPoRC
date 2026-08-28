# Set up observed survey indices and composition data

Ingests observed survey index, age composition, and length composition
data (both pooled and population-specific) into `input_list$data`,
initializes overdispersion and correlation starting values in
`input_list$par`, and constructs parameter maps via
[`do_comp_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_comp_theta_mapping.md)
and
[`do_comp_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_comp_corr_pars_mapping.md)
(called with `comp_prefix = "SrvAge"`/`"SrvLen"` and
`fleet_field = "n_srv_fleets"`). When `ISS_SrvAgeComps`,
`ISS_SrvLenComps`, `ISS_SrvAgeComps_pop`, or `ISS_SrvLenComps_pop` is
`NULL`, input sample sizes are derived automatically by summing observed
composition counts across the appropriate dimensions each year. Must be
called after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md)
and before model compilation.

## Usage

``` r
Setup_Mod_SrvIdx_and_Comps(
  input_list,
  ObsSrvIdx,
  ObsSrvIdx_SE,
  UseSrvIdx,
  ObsSrvIdxAA = NULL,
  UseSrvIdxAA = NULL,
  ObsSrvIdxAA_SE = NULL,
  ObsSrvIdxAA_pop = NULL,
  UseSrvIdxAA_pop = NULL,
  ObsSrvIdxAA_pop_SE = NULL,
  sigmaSrvIdxAA_key = NULL,
  sigmaSrvIdxAA_spec = "est",
  sigmaSrvIdxAA_pop_key = NULL,
  sigmaSrvIdxAA_pop_spec = "est",
  SrvIdxAA_Type = "spltRaggS",
  SrvIdxAA_pop_Type = "spltRaggS",
  SrvIdxAA_LikeType = "lognormal",
  SrvIdxAA_pop_LikeType = "lognormal",
  SrvIdxAA_sigma_form = "none",
  SrvIdxAA_pop_sigma_form = "none",
  AgeObsCorr_srv_idx = "iid",
  AgeObsCorr_srv_idx_pop = "iid",
  rho_srv_idx_spec = NULL,
  rho_srv_idx_pop_spec = NULL,
  sigmaSrvIdx_spec = "fix",
  sigmaSrvIdx_map = NULL,
  sigmaSrvIdx_pop_spec = "fix",
  sigmaSrvIdx_pop_map = NULL,
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
  srv_idx_ages = NULL,
  SrvAgeComps_bins = NULL,
  SrvLenComps_bins = NULL,
  Srv_caal_bins = NULL,
  SrvAgeComps_pop_bins = NULL,
  SrvLenComps_pop_bins = NULL,
  SrvIdx_LikeType = rep("lognormal", input_list$data$n_srv_fleets),
  SrvLenComps_sel = rep("age", input_list$data$n_srv_fleets),
  srv_waa_selected = rep(0, input_list$data$n_srv_fleets),
  SrvIdx_Cov = NULL,
  ObsSrv_caal = NULL,
  UseSrv_caal = NULL,
  ISS_Srv_caal = NULL,
  Srv_caal_LikeType = rep("none", input_list$data$n_srv_fleets),
  Srv_caal_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_srv_fleets,
    sep = ""),
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

- ObsSrvIdxAA:

  Observed survey index at age, an array with dimensions
  `[n_regions, n_years, n_seas, n_ages, n_sexes, n_srv_fleets]`.
  Supplying this fits the index at age directly, every age its own
  observation with its own catchability. The sex margin is required
  whatever the fleet reports: a stream summed over sexes carries its
  observation in sex slot one. A fleet uses this or the aggregated
  index, never both.

- UseSrvIdxAA:

  Integer array shaped like `ObsSrvIdxAA`, `1` where an observation is
  fit.

- ObsSrvIdxAA_SE, ObsSrvIdxAA_pop_SE:

  Reported standard errors shaped like their observation array, read
  only when `SrvIdxAA_sigma_form` asks for them. This is the parity the
  aggregated index already has: an index disaggregated by age keeps its
  survey-design errors.

- ObsSrvIdxAA_pop, UseSrvIdxAA_pop:

  Population-specific counterparts, with a leading population dimension.

- sigmaSrvIdxAA_key, sigmaSrvIdxAA_pop_key:

  Integer arrays `[n_ages, n_sexes, n_srv_fleets]` coupling the index at
  age observation error, the key matrix convention ICES assessments use.
  Equal entries share a parameter and `NA` excludes one. The sex margin
  is required; a key coupling the sexes repeats its entries across them.
  The age shape of catchability is not set here: an index fit age by age
  puts it in selectivity through the `"nonparfree"` form, which carries
  the height of the curve as well as its shape. See
  [`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md).

- sigmaSrvIdxAA_spec, sigmaSrvIdxAA_pop_spec:

  `"est"` (default) or `"fix"`.

- SrvIdxAA_Type, SrvIdxAA_pop_Type:

  Which margins the fleet reports separately: `"agg"`, `"spltRaggS"`
  (default), `"aggRspltS"` or `"spltRspltS"`. See
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- SrvIdxAA_LikeType, SrvIdxAA_pop_LikeType:

  `"lognormal"` (default) or `"normal"`, one setting for every fleet or
  one per fleet.

- SrvIdxAA_sigma_form, SrvIdxAA_pop_sigma_form:

  Where the observation error comes from: `"none"` (default), `"data"`,
  `"est_additive"` or `"est_quadrature"`.

- AgeObsCorr_srv_idx, AgeObsCorr_srv_idx_pop:

  Correlation across ages for the survey index at age, `"iid"`
  (default), `"1dar1"`, `"us"` or `"2dar1"`, one setting for every fleet
  or one per fleet. See
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- rho_srv_idx_spec, rho_srv_idx_pop_spec:

  How the correlation parameters are shared, over region, sex and fleet,
  using the package's spec strings. `NULL` (the default) gives one per
  fleet. See
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- sigmaSrvIdx_spec:

  Character string controlling the estimated component of the aggregated
  survey index observation error, one value per fleet. One of:

  `"fix"`

  :   The reported standard errors are used as they are and
      `ln_sigmaSrvIdx` is not estimated. The default.

  `"est_additive"`

  :   Total standard deviation is the reported standard error plus an
      estimated component, the additive extra standard deviation
      convention.

  `"est_quadrature"`

  :   Total standard deviation is the reported standard error and the
      estimated component added in quadrature, treating them as
      independent variances.

  `"est_replace"`

  :   An estimated standard deviation replaces the reported standard
      errors entirely, as several ICES assessments do.

  An estimated component is confounded with a likelihood weight, since a
  weight on a normal likelihood is the same statement as dividing the
  variance by that weight. `Setup_Mod_Weighting` warns when both are
  used. Fleets with a multivariate normal index likelihood take their
  scale from the supplied covariance and cannot carry one, which is an
  error rather than a silently unidentified parameter.

- sigmaSrvIdx_map:

  Optional integer vector of length `n_srv_fleets` giving the estimation
  groups for `ln_sigmaSrvIdx`. Fleets sharing a value share a parameter
  and `NA` holds a fleet at its starting value. Defaults to one free
  parameter per fleet. Use it when a reference assessment estimated some
  fleets and pinned others at a bound.

- sigmaSrvIdx_pop_spec:

  Character string controlling the estimated component of the
  population-specific survey index observation error, one value per
  fleet. One of:

  `"fix"`

  :   The reported standard errors are used as they are and
      `ln_sigmaSrvIdx_pop` is not estimated. The default.

  `"est_additive"`

  :   Total standard deviation is the reported standard error plus an
      estimated component, the additive extra standard deviation
      convention.

  `"est_quadrature"`

  :   Total standard deviation is the reported standard error and the
      estimated component added in quadrature, treating them as
      independent variances.

  `"est_replace"`

  :   An estimated standard deviation replaces the reported standard
      errors entirely, as several ICES assessments do.

  An estimated component is confounded with a likelihood weight, since a
  weight on a normal likelihood is the same statement as dividing the
  variance by that weight. `Setup_Mod_Weighting` warns when both are
  used. Fleets with a multivariate normal index likelihood take their
  scale from the supplied covariance and cannot carry one, which is an
  error rather than a silently unidentified parameter.

- sigmaSrvIdx_pop_map:

  Optional integer vector of length `n_srv_fleets` giving the estimation
  groups for `ln_sigmaSrvIdx_pop`. Fleets sharing a value share a
  parameter and `NA` holds a fleet at its starting value. Defaults to
  one free parameter per fleet. Use it when a reference assessment
  estimated some fleets and pinned others at a bound.

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
  One of `"biom"` (biomass), `"abd"` (abundance), `"recdev"`
  (recruitment deviations), or `"none"` (no index for that fleet).
  Converted to integer codes (`1`, `0`, `2`, `999`) before storage.

  A `"recdev"` fleet observes year class strength directly rather than
  any part of the population. Its predicted value is
  `q * (ln_RecDevs - mu)`, with `mu` the center the recruitment penalty
  asserts for that year, so it measures the anomaly rather than the
  deviation as stored; under a bias ramp the two differ. Such a fleet
  reads no numbers at age, so its selectivity, survey timing and weight
  at age are unused and its compositions should be left off. It requires
  `SrvIdx_LikeType = "normal"`, since deviations are signed, and
  `RecDevs_pen_center = "fixed"` in
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

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
  codes (`999`, `0`-`4`) before storage.

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
  `ObsSrvAgeComps_pop` within each population-year-fleet-season-region
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

- srv_idx_ages:

  Per-fleet selection of which ages contribute to the index total.
  Either a list with one element per survey fleet, where each element is
  a vector of ages or `NULL` for all ages, or an array
  `[n_ages x n_srv_fleets]` of 0/1 weights. Default `NULL` uses every
  age for every fleet. Restricting a fleet to a single age turns it into
  an index of that age alone, which is how an age-1 acoustic index is
  specified; the fleet's compositions are unaffected because the
  restriction applies to the index sum rather than to selectivity.

- SrvAgeComps_bins:

  Which age bins each survey fleet's age composition is fitted over.
  Supply a list with one element per fleet, each a vector of bin indices
  or `NULL` for all bins, or an `[n_obs_ages x n_srv_fleets]` array of
  0/1 weights. Both observed and expected compositions are restricted to
  the named bins and renormalized within them, so excluded bins are left
  out of the likelihood rather than being forced to be explained; this
  is how a fleet that only ages part of its age range is fitted. Indices
  refer to observed bins, that is after any ageing error has mapped
  model ages onto observed ones. The restriction applies whatever the
  composition type: for sex-joint comps the named bins are dropped from
  each sex's block, so the sex ratio the joint comps carry becomes the
  ratio within the fitted bins. Every fleet must retain at least two
  bins, since the proportion in a lone bin is one whatever the model
  predicts. Default `NULL`, which fits all bins for all fleets.

- SrvLenComps_bins:

  Which length bins each survey fleet's length composition is fitted
  over, in the same format as `SrvAgeComps_bins`. Indices refer to
  observed length bins, that is after any `LenBinMap` has mapped model
  bins onto observed ones.

- Srv_caal_bins:

  Which age bins each survey fleet's conditional age-at-length data are
  fitted over, in the same format as `SrvAgeComps_bins`. Applied to
  every length bin's row of ages alike.

- SrvAgeComps_pop_bins:

  Which age bins each survey fleet's population-specific age composition
  is fitted over, in the same format as `SrvAgeComps_bins`.

- SrvLenComps_pop_bins:

  Which length bins each survey fleet's population-specific length
  composition is fitted over, in the same format as `SrvAgeComps_bins`.

- SrvIdx_LikeType:

  Character vector `[n_srv_fleets]` giving the error structure of each
  survey index. Options are `"lognormal"` (default, the observation
  standard errors are on the log scale), `"normal"` (arithmetic scale),
  and `"mvn"` (multivariate normal on the arithmetic scale using a fixed
  covariance supplied through `SrvIdx_Cov`). One-step-ahead residuals
  are available only for lognormal fleets. A fleet's population-specific
  index stream follows the same choice for `"lognormal"` and `"normal"`,
  but stays lognormal under `"mvn"`, whose covariance describes the
  regional series only.

- SrvLenComps_sel:

  Character vector `[n_srv_fleets]`, whether a length-based selectivity
  is applied before or after the fish are spread over lengths. `"age"`
  (default) selects the index at age and spreads it afterwards;
  `"length"` spreads the numbers at each age over the key first and
  selects them length by length, so the survey sees the long fish of an
  age more often. The key is the survey's own, at `t_srv`. Requires
  length-based survey selectivity. Use `"length"` when selectivity is
  length based and the length compositions are what inform it.

- srv_waa_selected:

  Integer vector `[n_srv_fleets]` (0/1). With weight at age derived from
  growth and length-based selectivity, `1` makes a biomass index use the
  mean weight of the fish the survey sees at each age, \\\sum_l P(l
  \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)\\, instead of the
  population mean weight at that age. The survey twin of
  `fish_waa_selected`. Only applies to an index in weight
  (`srv_idx_type = "biom"`).

- SrvIdx_Cov:

  List with one element per survey fleet holding the fixed covariance
  matrix for fleets using `"mvn"`, and `NULL` otherwise. Each matrix
  must be square with one row per observation the fleet fits, ordered as
  the observations appear when scanning that fleet's `UseSrvIdx` slice
  in array order.

- ObsSrv_caal:

  Observed conditional age-at-length array
  `[n_regions x n_years x n_seas x n_lens x n_ages x n_sexes x n_srv_fleets]`.
  A CAAL observation is the age composition of the fish aged from one
  length bin, so the age margin of each length row is what gets fit.
  `NULL` (default) for a model with no CAAL data.

- UseSrv_caal:

  Use flags `[n_regions x n_years x n_seas x n_lens x n_srv_fleets]`.
  Length bins with no aged fish carry a zero and are skipped.

- ISS_Srv_caal:

  Input sample sizes
  `[n_regions x n_years x n_seas x n_lens x n_sexes x n_srv_fleets]`.
  Summed from `ObsSrv_caal` when `NULL`.

- Srv_caal_LikeType:

  Character vector of length `n_srv_fleets`. One of `"none"`,
  `"Multinomial"` or `"Dirichlet-Multinomial"`. The logistic-normal
  families are not available for CAAL, since a single length bin's age
  sample is small and mostly zeros, which the additive log-ratio
  transform cannot handle.

- Srv_caal_Type:

  Composition type specification, using the same
  `"CompType_Year_x-y_Fleet_z"` convention as the marginal compositions.

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
