# Set up fishing mortality, discard mortality, and catch observation inputs

Populates `input_list` with observed catch, catch usage indicators,
fishing mortality parameters (`ln_F_mean`, `ln_F_devs`), and
observation/process error structures (`ln_sigmaC`, `ln_sigmaC_pop`,
`ln_sigmaF`). Also populates discard observations, discard mortality
rate parameters (`logit_dmr_mean`, `logit_dmr_devs`), and discard
observation/process error structures (`ln_sigmaD`, `ln_sigmaD_pop`,
`ln_sigma_dmr`). Must be called after
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

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by upstream setup functions.

- ObsCatch:

  Observed aggregated catch array
  `[n_regions x n_years x n_seas x n_fish_fleets]`. Values should be in
  the units specified by `catch_units`. For a cell with `UseCatch == 0`
  (and no population-specific catch used), an `NA` entry here is treated
  as a genuinely missing observation; fishing is assumed to have
  continued and `Fmort`/ `ln_F_devs` are estimated normally for that
  year, whereas a true recorded value (typically `0`) is treated as a
  real closure: `Fmort` is forced to zero and no deviation is estimated.
  See
  [`Get_Fdev_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Fdev_PE_loglik.md).

- ObsCatchAA:

  Observed catch at age, an array with dimensions
  `[n_regions, n_years, n_seas, n_ages, n_sexes, n_fish_fleets]`. The
  sex margin is required whatever the fleet reports: a stream summed
  over sexes carries its observation in sex slot one. Supplying this
  fits the catch at age directly, every age its own lognormal
  observation, in place of an aggregated catch with compositions. This
  is the native form for ICES age-structured assessments. The two
  statements are not interchangeable: the exact factorization of an
  at-age observation into a total and a composition holds for Poisson
  and multinomial, not for lognormal, so a fleet must use one or the
  other and supplying both for the same fleet is an error. `NULL`
  (default) leaves the fleet on aggregated catch.

- UseCatchAA:

  Integer array shaped like `ObsCatchAA`, `1` where an observation is
  fit and `0` otherwise. A cell that is not fit is also not fished, so
  this governs closures the way `UseCatch` does for the aggregated
  stream.

- ObsCatchAA_SE, ObsDiscardAA_SE, ObsCatchAA_pop_SE,
  ObsDiscardAA_pop_SE:

  Reported standard errors shaped like their observation array, read
  only when the stream's `sigma_form` asks for them.

- sigmaCAA_key:

  Integer matrix `[n_ages, n_fish_fleets]` coupling the catch at age
  observation error, an integer key matrix, the convention ICES
  assessments use. Equal entries share a parameter and `NA` excludes
  one. This single structure covers every sharing pattern: `1 2 3 4 5`
  gives one standard deviation per age, `1 1 2 2 2` gives standard
  deviations by age group as several ICES assessments do, and
  `1 1 1 1 1` gives one for the fleet. Defaults to one parameter per
  fleet, shared across ages. A parameter informed by fewer than two
  observations is refused, since an observation error standard deviation
  with a single observation drives the likelihood to negative infinity
  rather than failing outright.

- sigmaCAA_spec:

  Character string, `"est"` (default) to estimate the coupled standard
  deviations, or `"fix"` to hold them at their starting values. Starting
  values are supplied through `...` as `ln_sigmaCAA`.

- ObsDiscardAA, UseDiscardAA:

  Observed discard at age and its use flags, shaped like `ObsCatchAA`.
  The discard counterpart of catch at age.

- ObsDiscardAA_pop, UseDiscardAA_pop, ObsCatchAA_pop, UseCatchAA_pop:

  Population-specific counterparts, with a leading population dimension.

- sigmaCAA_pop_key, sigmaDAA_key, sigmaDAA_pop_key:

  Integer matrices `[n_ages, n_fish_fleets]` coupling the observation
  error for the population-specific catch, the discards, and the
  population-specific discards, following the same convention as
  `sigmaCAA_key`.

- sigmaCAA_pop_spec, sigmaDAA_spec, sigmaDAA_pop_spec:

  `"est"` or `"fix"`.

- CatchAA_Type, DiscardAA_Type, CatchAA_pop_Type, DiscardAA_pop_Type:

  Which margins the fleet reports separately, following the composition
  vocabulary. Give it as one setting for every fleet, one per fleet, or
  as year and fleet specifications such as
  `"spltRaggS_Year_1-20_Fleet_1"` when the setting changes part way
  through the series. `"agg"` sums over regions and sexes, `"spltRaggS"`
  (default) splits regions and sums over sexes, `"aggRspltS"` does the
  reverse, and `"spltRspltS"` splits both. An observation summed over a
  margin belongs in slot one of it.

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
  unstructured correlation across ages, and `"2dar1"` correlates over
  ages and years jointly through a separable AR(1), which requires the
  fleet's observed ages and years to form a complete grid. A cell with a
  single observed age falls back to independent. The population-specific
  streams carry their own settings rather than borrowing the aggregated
  ones. The fishery and survey index streams are set in
  [`Setup_Mod_FishIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_FishIdx_and_Comps.md)
  and
  [`Setup_Mod_SrvIdx_and_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_SrvIdx_and_Comps.md).

- rho_catch_spec, rho_discard_spec, rho_catch_pop_spec,
  rho_discard_pop_spec:

  How each stream's correlation parameters are shared, using the same
  spec strings as `sigmaF_spec` and `Fdev_rho_spec`. The correlations
  sit over region, sex and fleet, with a leading population margin for
  the population-specific streams, so `"est_shared_r_s"` gives one per
  fleet, `"est_shared_s"` one per region and fleet, `"est_shared_r_s_f"`
  a single value, `"est_all"` one per cell, and `"fix"` holds them.
  `NULL` (the default) takes `"est_shared_r_s"`, or `"est_shared_p_r_s"`
  for the population streams, both one per fleet. The spec governs the
  across-age correlation, the across-year correlation and the
  unstructured matrix together, so fleets sharing under `"us"` share a
  whole matrix. A region, sex or population a fleet never observes
  carries no parameter, which is what holds the unused slots of a summed
  margin out.

- UseCatch:

  Binary indicator array
  `[n_regions x n_years x n_seas x n_fish_fleets]` controlling which
  aggregated catch observations enter the likelihood and whether
  `ln_F_devs` are estimated for each cell. `1` = use; `0` = exclude,
  unless `ObsCatch` is `NA` at that cell (see `ObsCatch` above), in
  which case `ln_F_devs` is still estimated as an ordinary active year
  despite not being fit against an observation.

- catch_units:

  Character array `[n_fish_fleets]` specifying catch units per fleet.
  `"biom"` = biomass (default); `"abd"` = abundance. Converted
  internally to `0`/`1` integer codes.

- UseCatch_pop:

  Binary indicator array
  `[n_pop x n_regions x n_years x n_seas x n_fish_fleets]` controlling
  which population-specific catch observations enter the likelihood. `1`
  = use; `0` = exclude.

- ObsCatch_pop:

  Observed population-specific catch array
  `[n_pop x n_regions x n_years x n_seas x n_fish_fleets]`. Values
  should be in the units specified by `catch_units`.

- Use_F_pen:

  Integer flag for applying a fishing mortality penalty to penalize
  large deviations in `ln_F_devs`. `1` = apply (default); `0` = do not
  apply.

- sigmaC_spec:

  Character string specifying the sharing structure for `ln_sigmaC`
  (aggregated catch observation error SD). Default `"fix"` holds
  `ln_sigmaC` at its starting value (`log(0.01)` unless overridden via
  `...`). Sharing options follow the convention `"est_shared_<dims>"`
  where `<dims>` is an underscore-separated list of dimensions to
  collapse: `"r"` (regions), `"y"` (years), `"seas"` (seasons), `"f"`
  (fleets), or any combination (e.g., `"est_shared_r_y"`,
  `"est_shared_r_y_seas_f"`). Use `"est_all"` for a fully independent
  parameter per cell. A warning is issued if `"fix"` is selected without
  providing a starting value in `...`.

- sigmaC_pop_spec:

  Character string specifying the sharing structure for `ln_sigmaC_pop`
  (population-specific catch observation error SD). Default `"fix"`
  holds `ln_sigmaC_pop` at its starting value (`log(0.01)` unless
  overridden via `...`). Sharing options follow the same convention as
  `sigmaC_spec` but with an additional population dimension: e.g.,
  `"est_shared_pop"` shares across populations, `"est_shared_pop_r"`
  shares across populations and regions, and
  `"est_shared_pop_r_y_seas_f"` collapses all dimensions into a single
  parameter. A warning is issued if `"fix"` is selected without
  providing a starting value in `...`.

- sigmaF_spec:

  Character string specifying the sharing structure for `ln_sigmaF`
  (fishing mortality process error SD). Default `"fix"` holds
  `ln_sigmaF` at its starting value (`log(1)`, i.e., \\\sigma_F = 1\\,
  unless overridden via `...`). A warning is issued if `"fix"` is
  selected without providing a starting value in `...`.

- Fdev_model:

  Character string specifying the process error structure for
  `ln_F_devs`. One of `"iid"` (default; independent deviations), `"rw"`
  (random walk; the first catch-active year per region/season/fleet is
  initialized with a diffuse \\N(0,5)\\ prior), or `"ar1"` (first-order
  autoregressive; the first catch-active year is drawn from its
  stationary marginal distribution, and `Fdev_rho_spec` controls the AR1
  correlation parameter). Catch-active years do not need to be
  contiguous for `"rw"` or `"ar1"`: the transition between two active
  years spanning a gap of \\d\\ closed years is taken over the elapsed
  gap directly (the same marginal transition as estimating deviations
  for the closed years and integrating them out, without actually
  estimating them), see
  [`Get_Fdev_PE_loglik`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Fdev_PE_loglik.md).
  A warning is issued if `"rw"` or `"ar1"` is selected but
  `Use_F_pen = 0` (the penalty is never evaluated, so the process
  structure has no effect), `sigmaF_spec = "fix"` (the process error SD
  is not estimated), or (for `"ar1"`) `Fdev_rho_spec = "fix"` (the
  correlation is not estimated), any of these may be intentional, but
  are common oversights when switching away from `"iid"`.

- Fdev_pen_center:

  Where the fishing mortality deviation penalty is centered. `"fixed"`
  (default) centers on zero, constraining both the level and the spread
  of the deviations. `"own_mean"` centers on the mean of the estimated
  deviations, penalizing only their spread and leaving the level free,
  which is what a sum of squares about the series' own mean amounts to.
  Under a mean-plus-deviations parameterization the level is already
  carried by `ln_F_mean`, so `"own_mean"` avoids penalizing it twice;
  note that it also leaves `ln_F_mean` and the deviations' level
  mutually unidentified unless one of them is fixed, which
  `ln_F_mean_spec = "fix"` does.

- Fdev_rho_spec:

  Character string specifying the sharing structure for the AR1
  correlation parameter `Fdev_rho`, following the same convention as
  `sigmaF_spec`. Only used when `Fdev_model = "ar1"`; ignored (and
  mapped entirely to `NA`) otherwise.

- ObsDiscard:

  Observed aggregated discard array
  `[n_regions x n_years x n_seas x n_fish_fleets]`. Values should be in
  the units specified by `discard_units`. Default: `NULL` (no discard
  observations).

- UseDiscard:

  Binary indicator array
  `[n_regions x n_years x n_seas x n_fish_fleets]` controlling which
  aggregated discard observations enter the likelihood. `1` = use; `0` =
  exclude. Default: all zeros.

- discard_units:

  Character array `[n_fish_fleets]` specifying discard units per fleet.
  `"abd"` = abundance (`0`), `"biom"` = biomass (`1`), `"abd_frac"` =
  abundance fraction (`2`), `"biom_frac"` = biomass fraction (`3`,
  default). Converted internally to integer codes.

- UseDiscard_pop:

  Binary indicator array
  `[n_pop x n_regions x n_years x n_seas x n_fish_fleets]` controlling
  which population-specific discard observations enter the likelihood.
  `1` = use; `0` = exclude. Default: all zeros.

- ObsDiscard_pop:

  Observed population-specific discard array
  `[n_pop x n_regions x n_years x n_seas x n_fish_fleets]`. Values
  should be in the units specified by `discard_units`. Default: `NULL`
  (no population-specific discard observations).

- Use_dmr_pen:

  Integer flag for applying a discard mortality rate penalty to penalize
  large deviations in `logit_dmr_devs`. `1` = apply; `0` = do not apply
  (default). Must be `1` when `dmr_dev_spec = "est_all"` and `0` when
  `dmr_dev_spec = "fix"`.

- sigmaD_spec:

  Character string specifying the sharing structure for `ln_sigmaD`
  (aggregated discard observation error SD). Default `"fix"` holds
  `ln_sigmaD` at its starting value (`log(0.01)` unless overridden via
  `...`). Sharing options follow the same convention as `sigmaC_spec`. A
  warning is issued if `"fix"` is selected without providing a starting
  value in `...`.

- sigmaD_pop_spec:

  Character string specifying the sharing structure for `ln_sigmaD_pop`
  (population-specific discard observation error SD). Default `"fix"`
  holds `ln_sigmaD_pop` at its starting value (`log(0.01)` unless
  overridden via `...`). Sharing options follow the same convention as
  `sigmaC_pop_spec`. A warning is issued if `"fix"` is selected without
  providing a starting value in `...`.

- sigma_dmr_spec:

  Character string specifying the sharing structure for `ln_sigma_dmr`
  (discard mortality rate process error SD). Default `"fix"` holds
  `ln_sigma_dmr` at its starting value (`log(1)` unless overridden via
  `...`). Sharing options follow the same convention as `sigmaF_spec`. A
  warning is issued if `"fix"` is selected without providing a starting
  value in `...`.

- dmr_mean_spec:

  Character string specifying the sharing/estimation structure for
  `logit_dmr_mean` (logit-scale mean discard mortality rate). Default
  `"fix"` holds at its starting value (`0`, i.e., DMR = 0.5 on the
  natural scale, unless overridden via `...`). See
  [`do_dmr_mean_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_dmr_mean_mapping.md)
  for sharing options.

- dmr_dev_spec:

  Character string specifying the sharing/estimation structure for
  `logit_dmr_devs` (logit-scale annual discard mortality rate
  deviations). Default `"fix"` holds deviations at zero (unless
  overridden via `...`). Use `"est_all"` to estimate a deviation in
  every fished cell; requires `Use_dmr_pen = 1`. See
  [`do_dmr_dev_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_dmr_dev_mapping.md)
  for sharing options.

- ...:

  Optional starting value overrides for catch and discard related
  parameters.

- ln_F_mean_spec:

  Character string, matched by exact name only because it sits after
  `...`. `"est"` (default, the previous and only behavior) or `"fix"`.
  `"fix"` maps `ln_F_mean` off at its starting value, which defaults to
  `0` under this spec unless supplied through `...`, so the deviations
  carry all of log fishing mortality: `F = exp(ln_F_devs)`, where it
  follows a free annual log-F parameterization. It must be paired with
  `Fdev_pen_center = "own_mean"` (penalize only the spread about the
  deviations' own mean), `Fdev_model = "rw"`, or `Use_F_pen = 0`: an
  `"iid"` or `"ar1"` penalty centered on a fixed zero mean would shrink
  the deviations toward `F = 1`, so that combination is rejected at
  setup. `"est"` keeps the mean-plus-deviations form, where the `"iid"`
  penalty shrinks each year toward the estimated average F.

## Value

The input `input_list` with `$data`, `$par`, and `$map` updated. Key
additions:

- `$data`:

  `ObsCatch`, `ObsCatch_pop`, `UseCatch`, `UseCatch_pop`, `Use_F_pen`,
  `catch_units`, `Fdev_model`, `ObsDiscard`, `ObsDiscard_pop`,
  `UseDiscard`, `UseDiscard_pop`, `Use_dmr_pen`, `discard_units`.

- `$par`:

  `ln_sigmaC`, `ln_sigmaC_pop`, `ln_sigmaF`, `Fdev_rho`, `ln_F_mean`,
  `ln_F_devs`, `ln_sigmaD`, `ln_sigmaD_pop`, `ln_sigma_dmr`,
  `logit_dmr_mean`, `logit_dmr_devs`.

- `$map`:

  `ln_sigmaC`, `ln_sigmaC_pop`, `ln_sigmaF`, `Fdev_rho`, `ln_F_mean`,
  `ln_F_devs`, `ln_sigmaD`, `ln_sigmaD_pop`, `ln_sigma_dmr`,
  `logit_dmr_mean`, `logit_dmr_devs`.

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
