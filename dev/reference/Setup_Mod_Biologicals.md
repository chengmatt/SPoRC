# Set up biological inputs for the estimation model

Populates `input_list` with biological arrays and parameter structures
needed by the TMB/RTMB objective function: weight-at-age (spawning,
fishery, and survey), maturity-at-age, ageing error, the size-age
transition matrix (optional), small constants for numerical stability,
and the natural mortality block structure and mapping. Called after
[`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md).

## Usage

``` r
Setup_Mod_Biologicals(
  input_list,
  WAA,
  WAA_fish = NULL,
  WAA_srv = NULL,
  MatAA,
  addtocomp = 0.001,
  addtofishidx = 1e-04,
  addtosrvidx = 1e-04,
  addtotag = 1e-10,
  AgeingError = NULL,
  Use_M_prior = 0,
  M_prior = NA,
  fit_lengths = 0,
  SizeAgeTrans = NA,
  M_spec = "est_ln_M",
  M_popblk_spec = "constant",
  M_ageblk_spec = "constant",
  M_regionblk_spec = "constant",
  M_yearblk_spec = "constant",
  M_sexblk_spec = "constant",
  Fixed_natmort = NULL,
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map`, and `$verbose` sublists, as
  returned by
  [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md).

- WAA:

  Numeric array of spawning weight-at-age with dimensions
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]`. Used to
  compute spawning stock biomass. Also serves as the fallback for
  `WAA_fish` and `WAA_srv` when those are `NULL`.

- WAA_fish:

  Numeric array of fishery weight-at-age with dimensions
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_fish_fleets]`.
  If `NULL` (default), `WAA` is broadcast across all fishery fleets.

- WAA_srv:

  Numeric array of survey weight-at-age with dimensions
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes × n_srv_fleets]`.
  If `NULL` (default), `WAA` is broadcast across all survey fleets.

- MatAA:

  Numeric array of maturity-at-age proportions (\\\in \[0,1\]\\) with
  dimensions
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]`.

- addtocomp:

  Small constant added to composition proportions before likelihood
  evaluation to avoid `log(0)`. Default `1e-3`. Ignored when a
  logistic-normal likelihood is specified, as that family handles zeros
  internally.

- addtofishidx:

  Small constant added to fishery indices. Default `1e-4`.

- addtosrvidx:

  Small constant added to survey indices. Default `1e-4`.

- addtotag:

  Small constant added to tag recovery observations. Default `1e-10`.

- AgeingError:

  Ageing error (age–age transition) array mapping true modelled ages to
  observed age bins. Accepted forms:

  2D matrix `[n_model_ages × n_obs_ages]`

  :   Time-invariant ageing error; replicated internally across all
      years.

  3D array `[n_years × n_model_ages × n_obs_ages]`

  :   Time-varying ageing error.

  `NULL` (default)

  :   An identity matrix is constructed, assuming modelled and observed
      age bins are identical. If observed bins are a subset of modelled
      ages (e.g., observed ages 2–10 vs. modelled ages 1–10), supply a
      shifted identity matrix such as
      `diag(1, n_model_ages)[, obs_age_index]` to avoid a dimensional
      mismatch.

- Use_M_prior:

  Integer flag to apply a lognormal prior on natural mortality. `0` = no
  prior (default); `1` = apply prior.

- M_prior:

  Data frame of prior hyperparameters for natural mortality, with one
  row per unique block combination. Required columns:

  `popblk`, `regionblk`, `yearblk`, `ageblk`, `sexblk`

  :   Block indices identifying which parameter the prior applies to.

  `mu`

  :   Prior mean in natural (untransformed) space.

  `sd`

  :   Prior standard deviation.

  Example for a single shared prior:

      M_prior <- data.frame(
        popblk = 1, regionblk = 1, yearblk = 1,
        ageblk = 1, sexblk = 1,
        mu = 0.085, sd = 0.05
      )

  Only used when `Use_M_prior = 1`.

- fit_lengths:

  Integer flag for fitting length compositions. `0` = no (default); `1`
  = yes. Requires a valid `SizeAgeTrans` array.

- SizeAgeTrans:

  Numeric array of size-at-age transition probabilities
  (column-stochastic; each age column sums to 1) with dimensions
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_ages × n_sexes]`.
  Required when `fit_lengths = 1`; ignored otherwise.

- M_spec:

  Character string controlling natural mortality estimation. One of:

  `"est_ln_M"` (default)

  :   Estimate `ln_M` across the defined blocks.

  `"fix"`

  :   Fix mortality to `Fixed_natmort`; `ln_M` parameters are mapped to
      `NA` and not passed to the optimizer.

- M_popblk_spec:

  Blocking structure for `ln_M` across populations. Either `"constant"`
  (default; single shared value) or a list of integer index vectors
  defining population groups, e.g., `list(1, 2)` for population-specific
  M.

- M_ageblk_spec:

  Blocking structure across ages. Either `"constant"` (default) or a
  list of integer index vectors, e.g., `list(1:5, 6:10)`.

- M_regionblk_spec:

  Blocking structure across regions. Either `"constant"` (default) or a
  list of integer index vectors, e.g., `list(1:3, 4:5)`.

- M_yearblk_spec:

  Blocking structure across years. Either `"constant"` (default) or a
  list of integer index vectors, e.g., `list(1:10, 11:30)`.

- M_sexblk_spec:

  Blocking structure across sexes. Either `"constant"` (default; shared
  across sexes) or a list of integer index vectors, e.g., `list(1, 2)`
  for sex-specific M.

- Fixed_natmort:

  Numeric array of fixed natural mortality rates with dimensions
  `[n_pop × n_regions × n_years × n_ages × n_sexes]`. Note the absence
  of an `n_seas` dimension. Required when `M_spec = "fix"`; ignored
  otherwise.

- ...:

  Optional starting value overrides passed by name. Currently
  recognised:

  `ln_M`

  :   Array of log-scale starting values for natural mortality,
      dimensioned
      `[n_popblks × n_regionblks × n_yearblks × n_ageblks × n_sexblks]`.
      Defaults to `log(0.5)` for all blocks if not supplied.

  All `...` arguments are silently ignored when `M_spec = "fix"`.

## Value

The input `input_list` with `$data`, `$par`, and `$map` sublists
updated. Key additions include `$data$WAA`, `$data$WAA_fish`,
`$data$WAA_srv`, `$data$MatAA`, `$data$AgeingError`, `$data$M_blocks`,
`$par$ln_M`, and `$map$ln_M`.

## See also

Other Model Setup:
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
[`Setup_Mod_Tagging()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Tagging.md),
[`Setup_Mod_Weighting()`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
