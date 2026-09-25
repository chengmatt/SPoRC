# Set up observed survey indices and composition data

Sets the observed survey indices and compositions, pooled and
population-specific, with the overdispersion and correlation starting
values and the maps from
[`do_comp_theta_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_comp_theta_mapping.md)
and
[`do_comp_corr_pars_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_comp_corr_pars_mapping.md).
A `NULL` `ISS_*` argument is summed from the matching observed array
each year, following that data source's composition type. Call after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md).

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
  SrvIdx_seas_Type = NULL,
  SrvIdx_pop_seas_Type = NULL,
  SrvIdxAA_seas_Type = NULL,
  SrvIdxAA_pop_seas_Type = NULL,
  SrvAgeComps_seas_Type = NULL,
  SrvAgeComps_pop_seas_Type = NULL,
  SrvLenComps_seas_Type = NULL,
  SrvLenComps_pop_seas_Type = NULL,
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

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- ObsSrvIdx:

  Observed survey index array
  `[n_regions × n_years × n_seas × n_srv_fleets]`.

- ObsSrvIdx_SE:

  Lognormal standard errors for `ObsSrvIdx`, same dims.

- UseSrvIdx:

  Binary array dimensioned like `ObsSrvIdx`, `1` to include the index in
  the likelihood.

- ObsSrvIdxAA:

  Observed survey index at age
  `[n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_srv_fleets]`, the
  ages being the columns of the fleet's ageing error matrix, through
  which the predicted index at each model age is read before it is
  compared. Supplying this fits the index at age directly, every age its
  own observation with its own catchability. The sex dim is required
  whatever the fleet reports: a data source summed over sexes has its
  observation in sex slot one. A fleet uses this or the aggregated
  index, never both.

- UseSrvIdxAA:

  Integer array shaped like `ObsSrvIdxAA`, `1` where an observation is
  fit.

- ObsSrvIdxAA_SE, ObsSrvIdxAA_pop_SE:

  Reported standard errors shaped like their observation array, read
  only when `SrvIdxAA_sigma_form` asks for them.

- ObsSrvIdxAA_pop, UseSrvIdxAA_pop:

  Population-specific counterparts, with a leading population dim.

- sigmaSrvIdxAA_key, sigmaSrvIdxAA_pop_key:

  Integer arrays `[n_obs_ages, n_sexes, n_srv_fleets]` coupling the
  index at age observation error, the key matrix ICES assessments use.
  Equal entries share a parameter and `NA` excludes one. The sex dim is
  required; a key coupling the sexes repeats its entries across them.
  The age shape of catchability is not set here: an index fit age by age
  puts it in selectivity through the `"nonparfree"` form. See
  [`Setup_Mod_Srvsel_and_Q`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Srvsel_and_Q.md).

- sigmaSrvIdxAA_spec, sigmaSrvIdxAA_pop_spec:

  `"est"` (default) or `"fix"`.

- SrvIdxAA_Type, SrvIdxAA_pop_Type:

  Which dims the fleet reports separately: `"agg"`, `"spltRaggS"`
  (default), `"aggRspltS"` or `"spltRspltS"`, or year and fleet
  specifications such as `"spltRaggS_Year_1-20_Fleet_1"`. See
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- SrvIdxAA_LikeType, SrvIdxAA_pop_LikeType:

  `"lognormal"` (default) or `"normal"`, one setting for every fleet or
  one per fleet.

- SrvIdxAA_sigma_form, SrvIdxAA_pop_sigma_form:

  Where the observation error comes from: `"none"` (default), `"data"`,
  `"est_additive"` or `"est_quadrature"`.

- AgeObsCorr_srv_idx, AgeObsCorr_srv_idx_pop:

  Correlation across ages for the survey index at age: `"iid"`
  (default), `"1dar1"`, `"us"` or `"2dar1"`, one setting for every fleet
  or one per fleet. See
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- rho_srv_idx_spec, rho_srv_idx_pop_spec:

  How the correlation parameters are shared over region, sex and fleet.
  `NULL` (default) gives one per fleet. See
  [`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

- sigmaSrvIdx_spec, sigmaSrvIdx_pop_spec:

  The estimated component of the aggregated and population-specific
  survey index observation error, one value per fleet. `"fix"` (default)
  uses the reported standard errors as they are. `"est_additive"` adds
  an estimated component to them, `"est_quadrature"` adds it in
  quadrature, and `"est_replace"` replaces them, as several ICES
  assessments do. An estimated component is confounded with a likelihood
  weight, since a weight on a normal likelihood is the same statement as
  dividing the variance, and `Setup_Mod_Weighting` warns when both are
  used. A fleet with a multivariate normal index takes its scale from
  the supplied covariance and cannot have one, which is an error.

- sigmaSrvIdx_map, sigmaSrvIdx_pop_map:

  Optional integer vectors `[n_srv_fleets]` of estimation groups for
  `ln_sigmaSrvIdx` and `ln_sigmaSrvIdx_pop`. Fleets sharing a value
  share a parameter and `NA` holds a fleet at its starting value.
  Defaults to one free parameter per fleet.

- ObsSrvIdx_pop:

  Observed population-specific index array
  `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`.

- ObsSrvIdx_pop_SE:

  Lognormal standard errors for `ObsSrvIdx_pop`, same dims.

- UseSrvIdx_pop:

  Binary array dimensioned like `ObsSrvIdx_pop`. Default all zeros.

- srv_idx_type:

  Character vector `[n_srv_fleets]`: `"biom"`, `"abd"`, `"recdev"` or
  `"none"`, stored as `1`, `0`, `2` and `999`. A `"recdev"` fleet
  observes year class strength directly rather than any part of the
  population: its predicted value is `q * (ln_RecDevs - mu)`, with `mu`
  the center the recruitment penalty asserts for that year, so it
  measures the anomaly rather than the deviation as stored. It reads no
  numbers at age, so its selectivity, timing and weight at age are
  unused and its compositions should be left off. It requires
  `SrvIdx_LikeType = "normal"` and `RecDevs_pen_center = "fixed"` in
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md).

- ObsSrvAgeComps:

  Observed survey age compositions
  `[n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]`,
  counts or proportions on a comparable scale.

- UseSrvAgeComps:

  Binary array `[n_regions × n_years × n_seas × n_srv_fleets]`, `1` to
  fit the age compositions.

- ObsSrvLenComps:

  Observed survey length compositions
  `[n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]`.
  Only validated when `fit_lengths = 1`.

- UseSrvLenComps:

  Binary array `[n_regions × n_years × n_seas × n_srv_fleets]`, `1` to
  fit the length compositions.

- ISS_SrvAgeComps:

  Input sample sizes
  `[n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`. `NULL` sums
  `ObsSrvAgeComps` over ages.

- ISS_SrvLenComps:

  Input sample sizes, structured as `ISS_SrvAgeComps`. `NULL` sums
  `ObsSrvLenComps`.

- SrvAgeComps_LikeType:

  Character vector `[n_srv_fleets]`: `"none"`, `"Multinomial"`,
  `"Dirichlet-Multinomial"`, the three logistic-normal forms or their
  three `-miss0` counterparts, stored as `999` and `0`-`7`.

  The `miss0` forms drop the empty bins and renormalize the expected
  proportions over the bins that remain, rather than adding `addtocomp`
  to the zeros. Their standard deviation is divided by the square root
  of the input sample size, so the parameter is a per fish quantity and
  a year sampled harder is fit more tightly, and the change of variables
  from the log ratio is taken off so the result is a density on the
  composition itself. Their correlations run through the logistic
  function and are therefore positive, matching the autoregression they
  mirror. The `2d` form needs a composition joint across sexes.
  One-step-ahead residuals are not available for any of the three.

- SrvLenComps_LikeType:

  As `SrvAgeComps_LikeType`, for length compositions.

- SrvAgeComps_Type:

  Character vector of the composition structure per fleet and year
  range, each `"<type>_Year_<start>-<end>_Fleet_<fleet>"` with
  `"terminal"` allowed as the end year. Types are `"agg"` (aggregated
  across regions and sexes, not valid with `"2d-Logistic-Normal"`),
  `"spltRspltS"`, `"spltRjntS"` and `"none"`. Parsed into an
  `[n_years × n_srv_fleets]` integer matrix; a cell left `NA` means an
  incomplete year range and is an error.

- SrvLenComps_Type:

  As `SrvAgeComps_Type`, for length compositions.

- ObsSrvAgeComps_pop:

  Observed population-specific age composition array
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]`.
  Required when any `UseSrvAgeComps_pop` is `1`.

- UseSrvAgeComps_pop:

  Binary array `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`.
  Default all zeros.

- ISS_SrvAgeComps_pop:

  Input sample size array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`.
  `NULL` (default) sums `ObsSrvAgeComps_pop`.

- ObsSrvLenComps_pop:

  Observed population-specific length composition array
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_srv_fleets]`.
  Required when `fit_lengths == 1` and any `UseSrvLenComps_pop` is `1`.

- UseSrvLenComps_pop:

  Binary array `[n_pop × n_regions × n_years × n_seas × n_srv_fleets]`.
  Default all zeros.

- ISS_SrvLenComps_pop:

  Input sample size array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_srv_fleets]`.
  `NULL` (default) sums `ObsSrvLenComps_pop`.

- SrvAgeComps_pop_LikeType, SrvLenComps_pop_LikeType:

  Character vectors `[n_srv_fleets]` for the population-specific
  compositions, with the same options as `SrvAgeComps_LikeType`. Default
  `"none"`.

- SrvAgeComps_pop_Type, SrvLenComps_pop_Type:

  Composition structure for the population-specific compositions, in the
  same format as `SrvAgeComps_Type`. Default `"none"` for every fleet
  and year.

- srv_idx_ages:

  Which ages contribute to each fleet's index total, either a list with
  one element per fleet (a vector of ages, or `NULL` for all) or an
  `[n_ages x n_srv_fleets]` array of 0/1 weights. `NULL` (default) uses
  every age. Restricting a fleet to one age makes it an index of that
  age alone, which is how an age-1 acoustic index is specified; the
  compositions are unaffected, since the restriction applies to the
  index sum.

- SrvAgeComps_bins:

  Which age bins each fleet's age composition is fitted over, either a
  list with one element per fleet (bin indices, or `NULL` for all) or an
  `[n_obs_ages x n_srv_fleets]` array of 0/1 weights. Observed and
  expected are both restricted to the named bins and renormalized within
  them, so excluded bins leave the likelihood rather than being forced
  to be explained. Indices are observed bins, after any ageing error.
  For sex-joint comps the named bins are dropped from each sex's block,
  so the sex ratio becomes the ratio within the fitted bins. Every fleet
  must keep at least two bins. Default `NULL`, all bins.

- SrvLenComps_bins:

  Which length bins each fleet's length composition is fitted over, as
  `SrvAgeComps_bins`. Indices are observed length bins, after any
  `LenBinMap`.

- Srv_caal_bins:

  Which age bins each fleet's CAAL data are fitted over, as
  `SrvAgeComps_bins`, applied to every length bin's row of ages alike.

- SrvAgeComps_pop_bins, SrvLenComps_pop_bins:

  Which bins each fleet's population-specific age and length
  compositions are fitted over, as `SrvAgeComps_bins`.

- SrvIdx_LikeType:

  Character vector `[n_srv_fleets]` of each index's error structure:
  `"lognormal"` (default, standard errors on the log scale), `"normal"`
  (arithmetic scale), or `"mvn"` (multivariate normal on the arithmetic
  scale with a fixed covariance from `SrvIdx_Cov`). One-step-ahead
  residuals are available for lognormal fleets only. A fleet's
  population-specific index follows the same choice for the first two
  and stays lognormal under `"mvn"`, whose covariance describes the
  regional series alone.

- SrvIdx_seas_Type, SrvIdx_pop_seas_Type, SrvIdxAA_seas_Type,
  SrvIdxAA_pop_seas_Type, SrvAgeComps_seas_Type,
  SrvAgeComps_pop_seas_Type, SrvLenComps_seas_Type,
  SrvLenComps_pop_seas_Type:

  Whether a seasonal model reports this data source once a season or
  once a year, one value for every fleet or one per fleet. `"spltSeas"`
  (default) fits the observation against the prediction for the season
  it sits in; `"aggSeas"` sums the prediction over the year's seasons
  and fits one observation. Under `"aggSeas"` the observation stays in
  the season it was placed in, exactly one season per region and year
  may be on in the matching `Use` array, and the likelihood lands in
  that season. A survey measured at a point in time belongs in its own
  season with its own timing.

- SrvLenComps_sel:

  Character vector `[n_srv_fleets]`, whether length-based selectivity
  applies before or after the fish are spread over lengths. `"age"`
  (default) selects the index at age and spreads it afterwards;
  `"length"` spreads first and selects length by length, so the survey
  sees the long fish of an age more often. The key is the survey's own
  at `t_srv`, and length-based survey selectivity is required.

- srv_waa_selected:

  Integer vector `[n_srv_fleets]` (0/1). With weight at age derived from
  growth and length-based selectivity, `1` makes a biomass index use the
  mean weight of the fish the survey sees at each age, \\\sum_l P(l
  \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)\\, rather than the
  population mean weight. The survey twin of `fish_waa_selected`, and
  only read for an index in weight.

- SrvIdx_Cov:

  List with one element per fleet holding the fixed covariance for
  `"mvn"` fleets and `NULL` otherwise. Each matrix is square with one
  row per observation the fleet fits, ordered as they appear when
  scanning that fleet's `UseSrvIdx` slice in array order.

- ObsSrv_caal:

  Observed conditional age-at-length array
  `[n_regions x n_years x n_seas x n_lens x n_ages x n_sexes x n_srv_fleets]`.
  An observation is the age composition of the fish aged from one length
  bin, so the age dim of each length row is what is fit. `NULL`
  (default) for no CAAL data.

- UseSrv_caal:

  Use flags `[n_regions x n_years x n_seas x n_lens x n_srv_fleets]`.
  Length bins with no aged fish take a zero and are skipped.

- ISS_Srv_caal:

  Input sample sizes
  `[n_regions x n_years x n_seas x n_lens x n_sexes x n_srv_fleets]`.
  Summed from `ObsSrv_caal` when `NULL`.

- Srv_caal_LikeType:

  Character vector `[n_srv_fleets]`: `"none"`, `"Multinomial"` or
  `"Dirichlet-Multinomial"`. The logistic-normal families are not
  available for CAAL, since a length bin's age sample is small and
  mostly zeros.

- Srv_caal_Type:

  Composition type, in the same `"CompType_Year_x-y_Fleet_z"` vocabulary
  as the marginal compositions.

- ...:

  Optional starting values for the overdispersion and correlation
  parameters.

## Value

`input_list` with the survey observations, use flags, input sample sizes
and integer-coded likelihood and composition types in `$data`, the
overdispersion and correlation starting values in `$par`, and their
factor maps, pooled and population-specific, in `$map`.

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
