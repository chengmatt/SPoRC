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
  ...
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

  Integer flag for applying a fishing mortality penalty to penalise
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

  Integer flag for applying a discard mortality rate penalty to penalise
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
