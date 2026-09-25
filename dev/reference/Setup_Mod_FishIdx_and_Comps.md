# Set up discards, fishery index, age composition, and length composition inputs

Sets the observed fishery indices and compositions, pooled and
population-specific, with their use flags, likelihood types, composition
types, input sample sizes and the overdispersion and correlation
starting values and maps. A `NULL` `ISS_*` argument is summed from the
matching observed array within each year, fleet, season and region cell,
following that data source's composition type. Call after
[`Setup_Mod_Catch_and_F`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Catch_and_F.md).

## Usage

``` r
Setup_Mod_FishIdx_and_Comps(
  input_list,
  ObsFishIdx,
  ObsFishIdx_SE,
  sigmaFishIdx_spec = "fix",
  sigmaFishIdx_map = NULL,
  sigmaFishIdx_pop_spec = "fix",
  sigmaFishIdx_pop_map = NULL,
  ObsFishIdx_pop = NULL,
  ObsFishIdx_pop_SE = NULL,
  UseFishIdx_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  fish_idx_type,
  t_fish = array(0, dim = c(input_list$data$n_regions, input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  FishLenComps_sel = rep("age", input_list$data$n_fish_fleets),
  fish_waa_selected = rep(0, input_list$data$n_fish_fleets),
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
  FishLenComps_bins = NULL,
  Fish_caal_bins = NULL,
  FishAgeComps_pop_bins = NULL,
  FishLenComps_pop_bins = NULL,
  FishIdx_LikeType = rep("lognormal", input_list$data$n_fish_fleets),
  FishIdx_seas_Type = NULL,
  FishIdx_pop_seas_Type = NULL,
  FishAgeComps_seas_Type = NULL,
  FishAgeComps_pop_seas_Type = NULL,
  FishLenComps_seas_Type = NULL,
  FishLenComps_pop_seas_Type = NULL,
  FishAgeComps_discard_seas_Type = NULL,
  FishAgeComps_discard_pop_seas_Type = NULL,
  FishLenComps_discard_seas_Type = NULL,
  FishLenComps_discard_pop_seas_Type = NULL,
  FishIdx_Cov = NULL,
  ObsFish_caal = NULL,
  UseFish_caal = NULL,
  ISS_Fish_caal = NULL,
  Fish_caal_LikeType = rep("none", input_list$data$n_fish_fleets),
  Fish_caal_Type = paste("none_Year_1-terminal_Fleet_", 1:input_list$data$n_fish_fleets,
    sep = ""),
  ObsFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, length(input_list$data$ages),
    input_list$data$n_sexes, input_list$data$n_fish_fleets)),
  UseFishAgeComps_discard = array(0, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ISS_FishAgeComps_discard = NULL,
  ObsFishLenComps_discard = array(0, dim = c(input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas, obs_len_bins(input_list),
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

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- ObsFishIdx:

  Observed fishery CPUE or biomass index array
  `[n_regions × n_years × n_seas × n_fish_fleets]`.

- ObsFishIdx_SE:

  Log-scale standard errors of `ObsFishIdx`, same dims.

- sigmaFishIdx_spec, sigmaFishIdx_pop_spec:

  The estimated component of the aggregated and population-specific
  fishery index observation error, one value per fleet. `"fix"`
  (default) uses the reported standard errors as they are.
  `"est_additive"` adds an estimated component to them,
  `"est_quadrature"` adds it in quadrature, and `"est_replace"` replaces
  them. An estimated component is confounded with a likelihood weight,
  since a weight on a normal likelihood is the same statement as
  dividing the variance, and `Setup_Mod_Weighting` warns when both are
  used. A fleet with a multivariate normal index takes its scale from
  the supplied covariance and cannot have one, which is an error.

- sigmaFishIdx_map, sigmaFishIdx_pop_map:

  Optional integer vectors `[n_fish_fleets]` of estimation groups for
  `ln_sigmaFishIdx` and `ln_sigmaFishIdx_pop`. Fleets sharing a value
  share a parameter and `NA` holds a fleet at its starting value.
  Defaults to one free parameter per fleet.

- ObsFishIdx_pop:

  Observed population-specific index array
  `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.

- ObsFishIdx_pop_SE:

  Lognormal standard errors for `ObsFishIdx_pop`, same dims.

- UseFishIdx_pop:

  Binary array dimensioned like `ObsFishIdx_pop`. Default all zeros.

- fish_idx_type:

  Character vector `[n_fish_fleets]`: `"biom"`, `"abd"`, or `"none"`.

- t_fish:

  Array `[n_regions x n_seas x n_fish_fleets]` of the fraction of the
  season elapsed when each index is observed. Numbers at age are decayed
  by `exp(-t_fish * ZAA)` before the index is formed, as `t_srv` does
  for surveys. Default `0`; use `0.5` for a mid-season index.

- FishLenComps_sel:

  Character vector `[n_fish_fleets]`, whether length-based selectivity
  applies before or after the fish are spread over lengths. `"age"`
  (default) selects at age and spreads afterwards, so the length
  composition within an age is the key's. `"length"` spreads first and
  selects length by length, so the long fish of an age are taken more
  often. The key is the fleet's own at `t_fish`, and length-based
  fishery selectivity is required. The two give different expected
  compositions.

- fish_waa_selected:

  Integer vector `[n_fish_fleets]` (0/1). With weight at age derived
  from growth and length-based selectivity, `1` makes the fleet's catch
  biomass use the mean weight of the fish it takes at each age, \\\sum_l
  P(l \mid a) s(l) w(l) / \sum_l P(l \mid a) s(l)\\, rather than the
  population mean weight. With flat or age-based selectivity the two
  agree.

- UseFishIdx:

  Binary array dimensioned like `ObsFishIdx`, `1` to include the index
  in the likelihood.

- ObsFishAgeComps:

  Observed fishery age composition array
  `[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`,
  raw counts or proportions; supply `ISS_FishAgeComps` for proportions.

- UseFishAgeComps:

  Binary array `[n_regions × n_years × n_seas × n_fish_fleets]`, `1` to
  fit the age compositions.

- ISS_FishAgeComps:

  Input sample size array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. `NULL`
  (default) sums `ObsFishAgeComps`.

- ObsFishLenComps:

  Observed fishery length composition array
  `[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.
  Only needed when `input_list$data$fit_lengths == 1`.

- UseFishLenComps:

  Binary array `[n_regions × n_years × n_seas × n_fish_fleets]`, `1` to
  fit the length compositions.

- ISS_FishLenComps:

  Input sample size array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. `NULL`
  (default) sums `ObsFishLenComps`.

- FishAgeComps_LikeType:

  Character vector `[n_fish_fleets]`: `"Multinomial"`,
  `"Dirichlet-Multinomial"`, `"iid-Logistic-Normal"`,
  `"1d-Logistic-Normal"`, `"2d-Logistic-Normal"`, the three `-miss0`
  forms, which drop the empty bins and scale the standard deviation by
  the input sample size, or `"none"`.

- FishLenComps_LikeType:

  As `FishAgeComps_LikeType`, for length compositions.

- FishAgeComps_Type:

  Character vector of the composition structure per fleet and time
  period, each `"<type>_Year_<start>-<end>_Fleet_<f>"` with `"terminal"`
  allowed as the end year. Types are `"agg"` (aggregated across regions
  and sexes, not valid with `"2d-Logistic-Normal"`), `"spltRspltS"`
  (split by region and sex), `"spltRjntS"` (split by region, joint
  across sexes) and `"none"`. For example
  `c("spltRjntS_Year_1-10_Fleet_1", "agg_Year_11-terminal_Fleet_1")`.

- FishLenComps_Type:

  As `FishAgeComps_Type`, for length compositions.

- ObsFishAgeComps_pop:

  Observed population-specific age composition array
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`.
  Required when any `UseFishAgeComps_pop` is `1`.

- UseFishAgeComps_pop:

  Binary array `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.
  Default all zeros.

- ISS_FishAgeComps_pop:

  Input sample size array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  `NULL` (default) sums `ObsFishAgeComps_pop`.

- ObsFishLenComps_pop:

  Observed population-specific length composition array
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.
  Required when `fit_lengths == 1` and any `UseFishLenComps_pop` is `1`.

- UseFishLenComps_pop:

  Binary array `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.
  Default all zeros.

- ISS_FishLenComps_pop:

  Input sample size array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  `NULL` (default) sums `ObsFishLenComps_pop`.

- FishAgeComps_pop_LikeType, FishLenComps_pop_LikeType:

  Character vectors `[n_fish_fleets]` for the population-specific
  compositions, with the same options as `FishAgeComps_LikeType`.
  Default `"none"`.

- FishAgeComps_pop_Type, FishLenComps_pop_Type:

  Composition structure for the population-specific compositions, in the
  same format as `FishAgeComps_Type`. Default `"none"` for every fleet
  and year.

- fish_idx_ages:

  Which ages contribute to each fleet's index total, either a list with
  one element per fleet (a vector of ages, or `NULL` for all) or an
  `[n_ages x n_fish_fleets]` array of 0/1 weights. `NULL` (default) uses
  every age. The compositions are unaffected.

- FishAgeComps_bins:

  Which age bins each fleet's age composition is fitted over, either a
  list with one element per fleet (bin indices, or `NULL` for all) or an
  `[n_obs_ages x n_fish_fleets]` array of 0/1 weights. Observed and
  expected are both restricted to the named bins and renormalized within
  them, so excluded bins leave the likelihood rather than being forced
  to be explained. Indices are observed bins, after any ageing error.
  For sex-joint comps the named bins are dropped from each sex's block,
  so the sex ratio becomes the ratio within the fitted bins. Every fleet
  must keep at least two bins. Default `NULL`, all bins.

- FishLenComps_bins:

  Which length bins each fleet's length composition is fitted over, as
  `FishAgeComps_bins`. Indices are observed length bins, after any
  `LenBinMap`.

- Fish_caal_bins:

  Which age bins each fleet's CAAL data are fitted over, as
  `FishAgeComps_bins`, applied to every length bin's row of ages alike.

- FishAgeComps_pop_bins, FishLenComps_pop_bins:

  Which bins each fleet's population-specific age and length
  compositions are fitted over, as `FishAgeComps_bins`.

- FishIdx_LikeType:

  Character vector `[n_fish_fleets]` of each index's error structure:
  `"lognormal"` (default, standard errors on the log scale), `"normal"`
  (arithmetic scale), or `"mvn"` (multivariate normal on the arithmetic
  scale with a fixed covariance from `FishIdx_Cov`). One-step-ahead
  residuals are available for lognormal fleets only. A fleet's
  population-specific index follows the same choice for the first two
  and stays lognormal under `"mvn"`, whose covariance describes the
  regional series alone.

- FishIdx_seas_Type, FishIdx_pop_seas_Type, FishAgeComps_seas_Type,
  FishAgeComps_pop_seas_Type, FishLenComps_seas_Type,
  FishLenComps_pop_seas_Type, FishAgeComps_discard_seas_Type,
  FishAgeComps_discard_pop_seas_Type, FishLenComps_discard_seas_Type,
  FishLenComps_discard_pop_seas_Type:

  Whether a seasonal model reports this data source once a season or
  once a year, one value for every fleet or one per fleet. `"spltSeas"`
  (default) fits the observation against the prediction for the season
  it sits in; `"aggSeas"` sums the prediction over the year's seasons
  and fits one observation. Under `"aggSeas"` the observation stays in
  the season it was placed in, exactly one season per region and year
  may be on in the matching `Use` array, and the likelihood lands in
  that season. An index measured at a point in time belongs in its own
  season.

- FishIdx_Cov:

  List with one element per fleet holding the fixed covariance for
  `"mvn"` fleets and `NULL` otherwise. Each matrix is square with one
  row per observation the fleet fits, ordered as they appear when
  scanning that fleet's `UseFishIdx` slice in array order.

- ObsFish_caal:

  Observed conditional age-at-length array
  `[n_regions x n_years x n_seas x n_lens x n_ages x n_sexes x n_fish_fleets]`.
  An observation is the age composition of the fish aged from one length
  bin, so the age dim of each length row is what is fit. `NULL`
  (default) for no CAAL data.

- UseFish_caal:

  Use flags `[n_regions x n_years x n_seas x n_lens x n_fish_fleets]`.
  Length bins with no aged fish take a zero and are skipped.

- ISS_Fish_caal:

  Input sample sizes
  `[n_regions x n_years x n_seas x n_lens x n_sexes x n_fish_fleets]`,
  the number aged within each length bin rather than the number
  measured. Summed from `ObsFish_caal` when `NULL`.

- Fish_caal_LikeType:

  Character vector `[n_fish_fleets]`: `"none"`, `"Multinomial"` or
  `"Dirichlet-Multinomial"`. The logistic-normal families are not
  available for CAAL, since a length bin's age sample is small and
  mostly zeros.

- Fish_caal_Type:

  Composition type, in the same `"CompType_Year_x-y_Fleet_z"` vocabulary
  as the marginal compositions.

- ObsFishAgeComps_discard:

  Observed discard age composition
  `[n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`,
  structured as `ObsFishAgeComps`.

- UseFishAgeComps_discard:

  Binary array `[n_regions × n_years × n_seas × n_fish_fleets]`.

- ISS_FishAgeComps_discard:

  Input sample size array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. `NULL`
  sums `ObsFishAgeComps_discard`.

- ObsFishLenComps_discard:

  Observed discard length composition
  `[n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.
  Required when `fit_lengths == 1`.

- UseFishLenComps_discard:

  Binary array `[n_regions × n_years × n_seas × n_fish_fleets]`.

- ISS_FishLenComps_discard:

  Input sample size array
  `[n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`. `NULL`
  sums `ObsFishLenComps_discard`.

- FishAgeComps_discard_LikeType, FishLenComps_discard_LikeType:

  Character vectors `[n_fish_fleets]` for the discard compositions, with
  the same options as `FishAgeComps_LikeType`.

- FishAgeComps_discard_Type, FishLenComps_discard_Type:

  Composition structure for the discard compositions, in the same format
  as `FishAgeComps_Type`.

- ObsFishAgeComps_discard_pop:

  Observed population-specific discard age composition
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`.

- UseFishAgeComps_discard_pop:

  Binary array `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.

- ISS_FishAgeComps_discard_pop:

  Input sample size array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  `NULL` sums `ObsFishAgeComps_discard_pop`.

- ObsFishLenComps_discard_pop:

  Observed population-specific discard length composition
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_sexes × n_fish_fleets]`.

- UseFishLenComps_discard_pop:

  Binary array `[n_pop × n_regions × n_years × n_seas × n_fish_fleets]`.

- ISS_FishLenComps_discard_pop:

  Input sample size array
  `[n_pop × n_regions × n_years × n_seas × n_sexes × n_fish_fleets]`.
  `NULL` sums `ObsFishLenComps_discard_pop`.

- FishAgeComps_discard_pop_LikeType, FishLenComps_discard_pop_LikeType:

  Character vectors `[n_fish_fleets]` for the population-specific
  discard compositions, with the same options as
  `FishAgeComps_LikeType`.

- FishAgeComps_discard_pop_Type, FishLenComps_discard_pop_Type:

  Composition structure for the population-specific discard
  compositions, in the same format as `FishAgeComps_Type`.

- ...:

  Optional starting values for the overdispersion and correlation
  parameters.

## Value

`input_list` with `$data`, `$par` and `$map` updated with the fishery
index and composition fields: the pooled and population-specific
observed arrays, the computed or supplied input sample sizes, the
integer-coded likelihood and composition type matrices, the
overdispersion parameters and their factor maps.

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
