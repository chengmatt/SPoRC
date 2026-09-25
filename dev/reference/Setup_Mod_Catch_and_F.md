# Set up fishing mortality, discard mortality, and catch observation inputs

Sets the observed catch and discards with their use flags, the fishing
mortality parameters (`ln_F_mean`, `ln_F_devs`) and their observation
and process error, the catch and discard at age data sources, and the
discard mortality rate parameters. Call after
[`Setup_Mod_Biologicals`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md).

## Usage

``` r
Setup_Mod_Catch_and_F(
  input_list,
  ObsCatch,
  ObsCatchAA = NULL,
  UseCatchAA = NULL,
  ObsCatchAA_SE = NULL,
  sigmaCAA_key = NULL,
  sigmaCAA_spec = "est",
  ObsDiscardAA = NULL,
  UseDiscardAA = NULL,
  ObsDiscardAA_SE = NULL,
  ObsDiscardAA_pop = NULL,
  UseDiscardAA_pop = NULL,
  ObsDiscardAA_pop_SE = NULL,
  ObsCatchAA_pop = NULL,
  UseCatchAA_pop = NULL,
  ObsCatchAA_pop_SE = NULL,
  sigmaCAA_pop_key = NULL,
  sigmaCAA_pop_spec = "est",
  sigmaDAA_key = NULL,
  sigmaDAA_spec = "est",
  sigmaDAA_pop_key = NULL,
  sigmaDAA_pop_spec = "est",
  CatchAA_Type = "spltRaggS",
  CatchAA_pop_Type = "spltRaggS",
  DiscardAA_Type = "spltRaggS",
  DiscardAA_pop_Type = "spltRaggS",
  Catch_seas_Type = NULL,
  Catch_pop_seas_Type = NULL,
  Discard_seas_Type = NULL,
  Discard_pop_seas_Type = NULL,
  CatchAA_seas_Type = NULL,
  CatchAA_pop_seas_Type = NULL,
  DiscardAA_seas_Type = NULL,
  DiscardAA_pop_seas_Type = NULL,
  CatchAA_LikeType = "lognormal",
  CatchAA_pop_LikeType = "lognormal",
  DiscardAA_LikeType = "lognormal",
  DiscardAA_pop_LikeType = "lognormal",
  CatchAA_sigma_form = "none",
  CatchAA_pop_sigma_form = "none",
  DiscardAA_sigma_form = "none",
  DiscardAA_pop_sigma_form = "none",
  AgeObsCorr_catch = "iid",
  AgeObsCorr_catch_pop = "iid",
  AgeObsCorr_discard = "iid",
  AgeObsCorr_discard_pop = "iid",
  rho_catch_spec = NULL,
  rho_catch_pop_spec = NULL,
  rho_discard_spec = NULL,
  rho_discard_pop_spec = NULL,
  UseCatch,
  catch_units = array("biom", dim = c(input_list$data$n_fish_fleets)),
  UseCatch_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ObsCatch_pop = NULL,
  Use_F_pen = 1,
  sigmaC_spec = "fix",
  sigmaC_pop_spec = "fix",
  sigmaF_spec = "fix",
  Fdev_model = "iid",
  Fdev_pen_center = "fixed",
  Fdev_rho_spec = "fix",
  ObsDiscard = NULL,
  UseDiscard = array(0, dim = c(input_list$data$n_regions, length(input_list$data$years),
    input_list$data$n_seas, input_list$data$n_fish_fleets)),
  discard_units = array("biom_frac", dim = c(input_list$data$n_fish_fleets)),
  UseDiscard_pop = array(0, dim = c(input_list$data$n_pop, input_list$data$n_regions,
    length(input_list$data$years), input_list$data$n_seas,
    input_list$data$n_fish_fleets)),
  ObsDiscard_pop = NULL,
  Use_dmr_pen = 0,
  sigmaD_spec = "fix",
  sigmaD_pop_spec = "fix",
  sigma_dmr_spec = "fix",
  dmr_mean_spec = "fix",
  dmr_dev_spec = "fix",
  ...,
  ln_F_mean_spec = "est"
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`.

- ObsCatch:

  Observed aggregated catch array
  `[n_regions x n_years x n_seas x n_fish_fleets]` in the units
  `catch_units` names. Where `UseCatch == 0` and no population-specific
  catch is used, an `NA` here is a missing observation: fishing is
  assumed to have continued and `Fmort` and `ln_F_devs` are estimated as
  usual. A recorded value, typically `0`, is a real closure: `Fmort` is
  forced to zero and no deviation is estimated. See
  [`Get_Fdev_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Fdev_PE_loglik.md).

- ObsCatchAA:

  Observed catch at age
  `[n_regions, n_years, n_seas, n_obs_ages, n_sexes, n_fish_fleets]`,
  the ages being the columns of the fleet's ageing error matrix, through
  which the predicted catch at each model age is read before it is
  compared. The sex dim is required whatever the fleet reports: a data
  source summed over sexes has its observation in sex slot one.
  Supplying this fits the catch at age directly, every age its own
  lognormal observation, in place of an aggregated catch with
  compositions, which is the native form for ICES age-structured
  assessments. The exact factorization of an at-age observation into a
  total and a composition holds for Poisson and multinomial but not
  lognormal, so a fleet must use one or the other and supplying both is
  an error. `NULL` (default) keeps the fleet on aggregated catch.

- UseCatchAA:

  Integer array shaped like `ObsCatchAA`, `1` where an observation is
  fit. A cell that is not fit is also not fished, so this governs
  closures the way `UseCatch` does.

- ObsCatchAA_SE, ObsDiscardAA_SE, ObsCatchAA_pop_SE,
  ObsDiscardAA_pop_SE:

  Reported standard errors shaped like their observation array, read
  only when that data source's `sigma_form` asks for them.

- sigmaCAA_key:

  Integer array `[n_obs_ages, n_sexes, n_fish_fleets]` coupling the
  catch at age observation error, the key matrix ICES assessments use.
  Equal entries share a parameter and `NA` excludes one. The sex dim is
  required; a key coupling the sexes repeats its entries across them.
  Along ages, `1 2 3 4 5` gives one sd per age, `1 1 2 2 2` gives sds by
  age group, and `1 1 1 1 1` gives one for the fleet. Defaults to one
  parameter per fleet. A parameter informed by fewer than two
  observations is refused, since an sd with a single observation drives
  the likelihood to negative infinity rather than failing outright.

- sigmaCAA_spec:

  `"est"` (default) or `"fix"`. Starting values go through `...` as
  `ln_sigmaCAA`.

- ObsDiscardAA, UseDiscardAA:

  Observed discard at age and its use flags, shaped like `ObsCatchAA`
  and read through the same fishery ageing error.

- ObsDiscardAA_pop, UseDiscardAA_pop, ObsCatchAA_pop, UseCatchAA_pop:

  Population-specific counterparts, with a leading population dim.

- sigmaCAA_pop_key, sigmaDAA_key, sigmaDAA_pop_key:

  Integer arrays coupling the observation error for the
  population-specific catch, the discards and the population-specific
  discards, following `sigmaCAA_key`. `sigmaDAA_key` is
  `[n_obs_ages, n_sexes, n_fish_fleets]`; the two population-specific
  keys take a leading population dim.

- sigmaCAA_pop_spec, sigmaDAA_spec, sigmaDAA_pop_spec:

  `"est"` or `"fix"`.

- CatchAA_Type, DiscardAA_Type, CatchAA_pop_Type, DiscardAA_pop_Type:

  Which dims the fleet reports separately, in the composition
  vocabulary, as one setting for every fleet, one per fleet, or year and
  fleet specifications such as `"spltRaggS_Year_1-20_Fleet_1"`. `"agg"`
  sums over regions and sexes, `"spltRaggS"` (default) splits regions
  and sums over sexes, `"aggRspltS"` does the reverse, and
  `"spltRspltS"` splits both. An observation summed over a dim belongs
  in slot one of it.

- Catch_seas_Type, Catch_pop_seas_Type, Discard_seas_Type,
  Discard_pop_seas_Type, CatchAA_seas_Type, CatchAA_pop_seas_Type,
  DiscardAA_seas_Type, DiscardAA_pop_seas_Type:

  Whether a seasonal model reports this data source once a season or
  once a year, one value for every fleet or one per fleet. `"spltSeas"`
  (default) fits the observation against the prediction for the season
  it sits in; `"aggSeas"` sums the prediction over the year's seasons
  and fits one observation, which is how a fleet that lands catch all
  year but reports one annual total is usually recorded. Under
  `"aggSeas"` the observation stays in the season it was placed in and
  exactly one season per region and year may be on in the matching `Use`
  array, since more than one would be fit against the same year total.
  Fishing mortality is still estimated season by season, so a fleet with
  one annual observation and free seasonal deviations leaves the split
  between seasons unidentified: share the deviations or fix the seasonal
  pattern.

- CatchAA_LikeType, DiscardAA_LikeType, CatchAA_pop_LikeType,
  DiscardAA_pop_LikeType:

  `"lognormal"` (default) or `"normal"`, one setting for every fleet or
  one per fleet.

- CatchAA_sigma_form, DiscardAA_sigma_form, CatchAA_pop_sigma_form,
  DiscardAA_pop_sigma_form:

  Where the observation error comes from. `"none"` (default) uses the
  estimated parameter alone, `"data"` the reported standard errors
  alone, and `"est_additive"` or `"est_quadrature"` both. Naming
  `"data"` holds the parameter fixed, since nothing reads it.

- AgeObsCorr_catch, AgeObsCorr_discard, AgeObsCorr_catch_pop,
  AgeObsCorr_discard_pop:

  Correlation across ages within a cell, one setting for every fleet or
  one per fleet. `"iid"` (default) treats ages as independent, `"1dar1"`
  correlates them as an AR(1) in age distance, `"us"` estimates an
  unstructured correlation, and `"2dar1"` correlates over ages and years
  jointly through a separable AR(1), which needs the fleet's observed
  ages and years to form a complete grid. A cell with one observed age
  falls back to independent. The population-specific data sources have
  their own settings. The index data sources are set in
  [`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md)
  and
  [`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md).

- rho_catch_spec, rho_discard_spec, rho_catch_pop_spec,
  rho_discard_pop_spec:

  How each data source's correlation parameters are shared, in the spec
  strings `sigmaF_spec` uses. The correlations sit over region, sex and
  fleet, with a leading population dim for the population-specific data
  sources, so `"est_shared_r_s"` gives one per fleet, `"est_shared_s"`
  one per region and fleet, `"est_shared_r_s_f"` a single value,
  `"est_all"` one per cell, and `"fix"` holds them. `NULL` (default)
  takes `"est_shared_r_s"`, or `"est_shared_p_r_s"` for the population
  data sources. The spec governs the across-age correlation, the
  across-year correlation and the unstructured matrix together, so
  fleets sharing under `"us"` share a whole matrix. A region, sex or
  population a fleet never observes has no parameter.

- UseCatch:

  Binary array dimensioned like `ObsCatch` controlling which aggregated
  catch observations enter the likelihood and whether `ln_F_devs` is
  estimated in each cell. `0` excludes the observation, unless
  `ObsCatch` is `NA` there, in which case the deviation is still
  estimated.

- catch_units:

  Character array `[n_fish_fleets]`: `"biom"` (default) or `"abd"`,
  stored as `0`/`1`.

- UseCatch_pop:

  Binary array dimensioned like `ObsCatch_pop`.

- ObsCatch_pop:

  Observed population-specific catch array
  `[n_pop x n_regions x n_years x n_seas x n_fish_fleets]`, in
  `catch_units`.

- Use_F_pen:

  Integer flag for the fishing mortality penalty on `ln_F_devs`. `1`
  (default) applies it.

- sigmaC_spec:

  Sharing structure for `ln_sigmaC`, the aggregated catch observation
  error sd. `"fix"` (default) holds it at its starting value,
  `log(0.01)` unless supplied through `...`, and warns when no starting
  value was given. Estimated options are `"est_shared_<dims>"` over any
  of `"r"` (regions), `"y"` (years), `"seas"` and `"f"` (fleets), e.g.
  `"est_shared_r_y_seas_f"`, or `"est_all"` for one parameter per cell.

- sigmaC_pop_spec:

  Sharing structure for `ln_sigmaC_pop`, as `sigmaC_spec` with an added
  population dim, e.g. `"est_shared_pop_r"` or
  `"est_shared_pop_r_y_seas_f"`.

- sigmaF_spec:

  Sharing structure for `ln_sigmaF`, the fishing mortality process error
  sd, following `sigmaC_spec`. `"fix"` (default) holds it at `log(1)`
  unless supplied through `...`, and warns.

- Fdev_model:

  Process error on `ln_F_devs`: `"iid"` (default), `"rw"` (the first
  catch-active year per region, season and fleet takes a diffuse
  \\N(0,5)\\), or `"ar1"` (that year is drawn from the stationary
  marginal, with `Fdev_rho_spec` setting the correlation). Catch-active
  years need not be contiguous under `"rw"` or `"ar1"`: the transition
  across a gap of \\d\\ closed years is taken over the elapsed gap, the
  same marginal as estimating the closed years and integrating them out.
  See
  [`Get_Fdev_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Fdev_PE_loglik.md).
  Warns under `Use_F_pen = 0` (the penalty is never evaluated),
  `sigmaF_spec = "fix"`, or, for `"ar1"`, `Fdev_rho_spec = "fix"`.

- Fdev_pen_center:

  Where the fishing mortality deviation penalty is centered. `"fixed"`
  (default) centers on zero, constraining the level and the spread.
  `"own_mean"` centers on the deviations' own mean, penalizing only
  their spread; the level is then already set by `ln_F_mean`, so it is
  not penalized twice, but the two are mutually unidentified unless one
  is fixed, which `ln_F_mean_spec = "fix"` does.

- Fdev_rho_spec:

  Sharing structure for `Fdev_rho`, following `sigmaF_spec`. Only read
  under `Fdev_model = "ar1"` and mapped entirely to `NA` otherwise.

- ObsDiscard:

  Observed aggregated discard array
  `[n_regions x n_years x n_seas x n_fish_fleets]` in `discard_units`.
  Default `NULL`.

- UseDiscard:

  Binary array dimensioned like `ObsDiscard`. Default all zeros.

- discard_units:

  Character array `[n_fish_fleets]`: `"abd"` (`0`), `"biom"` (`1`),
  `"abd_frac"` (`2`) or `"biom_frac"` (`3`, default).

- UseDiscard_pop:

  Binary array `[n_pop x n_regions x n_years x n_seas x n_fish_fleets]`.
  Default all zeros.

- ObsDiscard_pop:

  Observed population-specific discard array, same dims, in
  `discard_units`. Default `NULL`.

- Use_dmr_pen:

  Integer flag for the penalty on `logit_dmr_devs`. Default `0`. Must be
  `1` under `dmr_dev_spec = "est_all"` and `0` under `"fix"`.

- sigmaD_spec, sigmaD_pop_spec:

  Sharing structures for `ln_sigmaD` and `ln_sigmaD_pop`, the discard
  observation error sds, following `sigmaC_spec` and `sigmaC_pop_spec`.
  `"fix"` (default) holds them at `log(0.01)` and warns when no starting
  value was given.

- sigma_dmr_spec:

  Sharing structure for `ln_sigma_dmr`, the discard mortality rate
  process error sd, following `sigmaF_spec`. `"fix"` (default) holds it
  at `log(1)` and warns.

- dmr_mean_spec:

  Sharing structure for `logit_dmr_mean`. `"fix"` (default) holds it at
  `0`, a rate of 0.5 on the natural scale. See
  [`do_dmr_mean_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_dmr_mean_mapping.md).

- dmr_dev_spec:

  Sharing structure for `logit_dmr_devs`. `"fix"` (default) holds the
  deviations at zero; `"est_all"` estimates one in every fished cell and
  requires `Use_dmr_pen = 1`. See
  [`do_dmr_dev_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_dmr_dev_mapping.md).

- ...:

  Optional starting values for the catch and discard parameters.

- ln_F_mean_spec:

  `"est"` (default) or `"fix"`, matched by exact name only because it
  sits after `...`. `"fix"` maps `ln_F_mean` off at its starting value,
  `0` unless supplied through `...`, so the deviations hold all of log
  fishing mortality, `F = exp(ln_F_devs)`. It must be paired with
  `Fdev_pen_center = "own_mean"`, `Fdev_model = "rw"` or
  `Use_F_pen = 0`: an `"iid"` or `"ar1"` penalty centered on a fixed
  zero would shrink the deviations toward `F = 1`, so that combination
  is rejected at setup.

## Value

`input_list` with `$data`, `$par` and `$map` updated. `$data` gains
`ObsCatch`, `ObsCatch_pop`, `UseCatch`, `UseCatch_pop`, `Use_F_pen`,
`catch_units`, `Fdev_model`, `ObsDiscard`, `ObsDiscard_pop`,
`UseDiscard`, `UseDiscard_pop`, `Use_dmr_pen` and `discard_units`.
`$par` and `$map` both gain `ln_sigmaC`, `ln_sigmaC_pop`, `ln_sigmaF`,
`Fdev_rho`, `ln_F_mean`, `ln_F_devs`, `ln_sigmaD`, `ln_sigmaD_pop`,
`ln_sigma_dmr`, `logit_dmr_mean` and `logit_dmr_devs`.

## See also

Other Model Setup:
[`Setup_Mod_Biologicals()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Biologicals.md),
[`Setup_Mod_Dim()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md),
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
