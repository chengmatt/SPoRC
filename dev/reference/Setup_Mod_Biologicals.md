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
  addtocomp = NULL,
  comp_const_obs = NULL,
  addtofishidx = NULL,
  addtosrvidx = NULL,
  addtotag = NULL,
  AgeingError = NULL,
  AgeingError_fish = NULL,
  AgeingError_srv = NULL,
  Use_M_prior = 0,
  M_prior = NA,
  fit_lengths = 0,
  SizeAgeTrans = NA,
  SizeAgeTrans_fish = NULL,
  SizeAgeTrans_srv = NULL,
  do_caal = 0,
  growth_model = "none",
  growth_spec = "est_all",
  growth_fix = NULL,
  growth_tv_model = NULL,
  growth_tv_years = NULL,
  growth_tv_link = "log",
  growth_par_bounds = NULL,
  growth_tv_sigma_spec = "fix",
  growth_tv_spec = "est_all",
  growth_tv_type = "curve",
  growth_rw_init_sigma = 5,
  growth_semipar = "none",
  growth_semipar_spec = "fix",
  growth_semipar_ages = NULL,
  growth_semipar_years = NULL,
  LenBinMap = NULL,
  growth_A1 = NULL,
  growth_A2 = NULL,
  growth_len_lower = NULL,
  growth_L0 = NULL,
  growth_cv_type = "len",
  growth_sd_type = "cv",
  growth_dist = "normal",
  growth_plus_group = "mixture",
  waa_model = "data",
  wt_len_pars = NULL,
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
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]`. When
  `rec_lag = 0` (age-0 recruitment, set via
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)),
  maturity at the recruit age (the first age class) must be exactly `0`
  for all populations, regions, years, seasons, and sexes, an error is
  raised otherwise. Requires
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
  to have been called first so `rec_lag` is already set.

- addtocomp:

  **Deprecated here**, pass it to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
  instead, which now owns this constant along with every other
  likelihood weight. Still accepted for backward compatibility: if
  supplied, it is forwarded to `Setup_Mod_Weighting` with a message
  rather than applied here directly. (Small constant added to
  composition proportions before likelihood evaluation to avoid
  `log(0)`; default `1e-3` in `Setup_Mod_Weighting`. Ignored when a
  logistic-normal likelihood is specified, as that family handles zeros
  internally.)

- comp_const_obs:

  **Deprecated here**, pass it to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
  instead. Still accepted for backward compatibility (forwarded with a
  message). Integer switch (`0` or `1`) controlling where `addtocomp` is
  applied in the multinomial likelihood, not a constant to be tuned. `1`
  (default in `Setup_Mod_Weighting`) adds it to the observed proportions
  that weight the multinomial as well as inside the logarithms, so the
  likelihood is stationary exactly at `pred = obs`. `0` weights by the
  raw observed proportions. The Dirichlet-multinomial sanity check that
  used to read it here (inside `Setup_Mod_FishIdx_and_Comps`/
  `Setup_Mod_SrvIdx_and_Comps`) now runs inside `Setup_Mod_Weighting`
  once the final value is known.

- addtofishidx:

  **Deprecated here**, pass it to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
  instead. Still accepted for backward compatibility (forwarded with a
  message). Small constant added to fishery indices; default `1e-4` in
  `Setup_Mod_Weighting`.

- addtosrvidx:

  **Deprecated here**, pass it to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
  instead. Still accepted for backward compatibility (forwarded with a
  message). Small constant added to survey indices; default `1e-4` in
  `Setup_Mod_Weighting`.

- addtotag:

  **Deprecated here**, pass it to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md)
  instead. Still accepted for backward compatibility (forwarded with a
  message). Small constant added to tag recovery observations; default
  `1e-10` in `Setup_Mod_Weighting`.

- AgeingError:

  Ageing error (age-age transition) array mapping true modeled ages to
  observed age bins. Each row is one model age's share across the
  observed bins and sums to one, or to zero to drop that model age from
  the observations. This is the age-axis twin of `LenBinMap`: the
  likelihood applies the two identically and validates them identically,
  so read either one for the other. It changes which bins the
  compositions are recorded on; to leave observed bins out of the
  likelihood without changing the bins themselves, use the `*_bins`
  arguments instead. Accepted forms:

  2D matrix `[n_model_ages × n_obs_ages]`

  :   Time-invariant ageing error; replicated internally across all
      years.

  3D array `[n_years × n_model_ages × n_obs_ages]`

  :   Time-varying ageing error.

  `NULL` (default)

  :   An identity matrix is constructed, assuming modeled and observed
      age bins are identical. If observed bins are a subset of modeled
      ages (e.g., observed ages 2-10 vs. modeled ages 1-10), supply a
      shifted identity matrix such as
      `diag(1, n_model_ages)[, obs_age_index]` to avoid a dimensional
      mismatch.

- AgeingError_fish:

  Optional fleet-specific ageing error for the fishery fleets, for when
  the fleets do not read ages the same way. Accepted forms: a 3D array
  `[n_model_ages × n_obs_ages × n_fish_fleets]` for a time-invariant
  matrix per fleet, a 4D array
  `[n_years × n_model_ages × n_obs_ages × n_fish_fleets]` for a
  time-varying one, or `NULL` (default), which gives every fishery fleet
  the shared `AgeingError`. Each fleet's slice is validated the same way
  `AgeingError` is, and every fleet must land on the same observed age
  bins, since the observed composition arrays carry one age dimension
  shared across fleets.

- AgeingError_srv:

  Optional fleet-specific ageing error for the survey fleets, in the
  same forms as `AgeingError_fish`, with `n_srv_fleets` in place of
  `n_fish_fleets`. `NULL` (default) gives every survey fleet the shared
  `AgeingError`.

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
  Required when `fit_lengths = 1`; ignored otherwise. The shared key
  every fleet reads unless `SizeAgeTrans_fish`/`SizeAgeTrans_srv`
  override it for that fleet type.

- SizeAgeTrans_fish, SizeAgeTrans_srv:

  Optional per-fleet size-at-age transition arrays, dimensioned like
  `SizeAgeTrans` with an added trailing fleet dimension
  (`n_fish_fleets`/`n_srv_fleets`). `NULL` (default) reads every fleet's
  key from the shared `SizeAgeTrans`. Only meaningful with
  `growth_model = "none"`; a growth model already derives one key per
  fleet, at that fleet's own timing, and rejects these to avoid mixing
  two sources for the same key. This is the fixed-data counterpart of
  `Setup_Sim_Biologicals`'s
  `SizeAgeTrans_fish_input`/`SizeAgeTrans_srv_input`, and of
  `WAA_fish`/`WAA_srv` overriding the shared `WAA`.

- do_caal:

  Integer flag for building the joint arrays at length and age. `0` = no
  (default); `1` = yes. Requires `fit_lengths = 1`. Turning this on adds
  `Fish_caal`, `Fish_caal_discard` and `Srv_caal` to the report, holding
  predicted retained catch, discards and survey index jointly by length
  and age.

- growth_model:

  Character. `"none"` (default) keeps `SizeAgeTrans` and the
  weight-at-age arrays as data. `"vb_schnute"` builds the size-age
  transition from estimable von Bertalanffy parameters in Schnute's
  form: length `L1` at reference age `growth_A1`, length `L2` at
  `growth_A2`, rate `K`, and CVs of length at age `CV1` and `CV2` at the
  two reference ages. Growth below `growth_A1` is linear from
  `growth_L0` at age zero, the CV interpolates between the two
  references, and the plus group carries an adjustment for fish older
  than the accumulator age. `"richards"` is the same curve with a sixth
  parameter, the Richards coefficient `rho`, applied to the lengths
  raised to that power (`rho = 1` recovers the von Bertalanffy form).
  Requires `fit_lengths = 1`; `SizeAgeTrans` is then ignored and may be
  `NA`.

- growth_spec:

  Character. How the growth parameters are estimated: `"est_all"`
  (default, one set per population, region and sex), `"est_shared_r"`
  (shared across regions), `"est_shared_s"` (shared across sexes),
  `"est_shared_r_s"` (one set per population), or `"fix"`.

- growth_fix:

  Logical vector, one entry per growth parameter, naming which of L1,
  L2, K, CV1, CV2 (and rho) stay at their starting values whatever
  `growth_spec` says.

- growth_tv_model:

  Time variation of the growth parameters. `NULL` (default) holds every
  parameter constant. Otherwise a character vector naming a structure
  per parameter, either of length `n_gpars` in the parameter order or
  named by parameter (`L1`, `L2`, `K`, `CV1`, `CV2`, `rho`) with the
  rest constant, each one of `"none"`, `"iid"` (independent annual
  deviations) or `"rw"` (a random walk). A varying parameter gets a
  deviation series `ln_growth_devs` and a log sigma in the first stream
  of `growth_pe_pars`.

- growth_tv_years:

  Years the deviations are active in, calendar years. `NULL` (default)
  for every model year, a vector applied to every varying parameter, or
  a list named by parameter. Deviations outside the range are held at
  zero.

- growth_tv_link:

  Character, the scale a deviation enters on. `"log"` (default)
  multiplies the parameter by \\e^{\delta}\\; `"logit"` keeps it inside
  `growth_par_bounds`, \\P_y = lo + (hi -
  lo)\\\mathrm{logit}^{-1}(\mathrm{logit}((P - lo)/(hi - lo)) +
  \delta_y)\\, so the parameter approaches a bound however large the
  deviation instead of crossing it.

- growth_par_bounds:

  Matrix `[n_gpars x 2]` of lower and upper bounds, natural scale,
  required under the logit link.

- growth_tv_sigma_spec:

  Character, `"fix"` (default) holds the process error standard
  deviations of the deviations at their starting values, `"est"`
  estimates them. Both read the first stream of `growth_pe_pars`, one
  slot per growth parameter.

- growth_tv_spec:

  Character, how the deviations are shared across strata, with the same
  vocabulary as `growth_spec`: `"est_all"` (default), `"est_shared_r"`,
  `"est_shared_s"` or `"est_shared_r_s"`.

- growth_tv_type:

  Character. `"curve"` (default) reads every year's size at age off that
  year's curve. `"cohort"` carries size at age forward cohort by cohort:
  each year every cohort grows by the increment the current year's
  parameters imply from the size it reached, ages still in the linear
  phase keep the length at `growth_A1` their birth year's parameters
  gave them, the first age past `growth_A1` is placed on the current
  year's curve, and the plus group's size blends the cohort entering it
  with the fish already there by their numbers at age. The CV at age is
  then held at the first year's sizes. The propagation starts in the
  first year any deviation is active; every earlier year sits on the
  first year's curve.

- growth_rw_init_sigma:

  Standard deviation given to the first year of a random walk on a
  growth parameter, as `srvsel_rw_init_sigma` for selectivity. Default
  `5`.

- growth_semipar:

  Character. Semi-parametric growth: a year-by-age surface of deviations
  on mean length at age, multiplying the parametric curve, so the curve
  stays the parametric part and the deviations hold departures from it.
  `"none"` (default) keeps growth purely parametric; otherwise one of
  `"iid"`, `"rw"` (a random walk over years within an age), `"3dmarg"`
  or `"3dcond"` (a three-dimensional Gaussian Markov random field over
  age, year and cohort, on the marginal or conditional variance), or
  `"2dar1"` (a separable first-order autoregression over ages and
  years). The same process error forms the selectivity deviations use,
  so a growth surface and a selectivity surface are scored the same way.
  The spread at age follows the deviated mean, which leaves the
  coefficient of variation at age to the parametric part.

- growth_semipar_spec:

  Character, whether the second stream of `growth_pe_pars` is estimated.
  Whether the process error hyperparameters are estimated (`"est"`) or
  held at their starting values (`"fix"`, the default). The deviations
  themselves are always estimated.

- growth_semipar_ages:

  Ages the deviations are estimated over, as ages (not indices). `NULL`
  (default) uses every age. Ages outside the set are held at zero, which
  is how a surface is restricted to the ages the length data actually
  inform.

- growth_semipar_years:

  Years the deviations are estimated over, calendar years. `NULL`
  (default) uses every year.

- LenBinMap:

  Optional matrix `[n_lens x n_obs_lens]` mapping the model's length
  bins onto the bins the length compositions are recorded on, for
  compositions on coarser bins than the model carries (a population of 1
  cm bins fit to 5 cm compositions, say). Observed length compositions
  are then dimensioned by `n_obs_lens` and the expected compositions are
  mapped through it inside the likelihood. This is the length-axis twin
  of `AgeingError`: the likelihood applies the two identically and
  validates them identically, so read either one for the other. Each row
  is one model bin's share across the observed bins and sums to one, or
  to zero to drop that model bin from the observations. It changes which
  bins the compositions are recorded on; to leave observed bins out of
  the likelihood without changing the bins themselves, use the
  `*LenComps_bins` arguments instead. `NULL` (default) fits the
  compositions on the model bins.

- growth_A1, growth_A2:

  Reference ages for `L1` and `L2`. `growth_A2 = "Linf"` instead makes
  `L2` the asymptotic length itself, with no second reference age to
  solve it from.

- growth_len_lower:

  Numeric vector of the lower edges of the length bins. `lens` in
  `Setup_Mod_Dim` are bin midpoints; the key is built on the edges.

- growth_L0:

  Length at age zero anchoring the linear phase. Defaults to
  `growth_len_lower[1]`.

- growth_cv_type:

  Character, `"len"` (default) interpolates the CV on mean length
  between `L1` and `L2`, `"age"` on age.

- growth_sd_type:

  Character, `"cv"` (default) scales the mean by the CV parameters,
  `"sd"` reads them as standard deviations.

- growth_dist:

  Character, `"normal"` (default) or `"lognormal"` distribution of
  length at age.

- growth_plus_group:

  Character. `"mixture"` (default) takes the plus group's mean length as
  the survivorship-weighted mixture of the ages it holds, their numbers
  declining at an assumed 0.2 per year and their length rising from the
  curve at the accumulator age to the asymptote; `"curve"` reads the
  curve at the accumulator age.

- waa_model:

  Character. Where weight at age comes from. `"data"` (default) reads
  `WAA`, `WAA_fish` and `WAA_srv` from the arguments of the same name.
  `"wt_len"` builds them from the size-age key and the weight-length
  relationship \\W = a L^b\\ applied at the bin midpoints, so weight at
  age carries the spread of length at age rather than being the weight
  of the mean length; the spawning weight uses the key at spawning time
  and each fleet's weight the key at that fleet's timing, `t_fish` or
  `t_srv`. Under `"wt_len"`, `WAA` may be `NULL`, and reference point
  and projection code still read `data$WAA`, so copy the reported arrays
  into the data list before calling them.

- wt_len_pars:

  Weight-length parameters \\a, b\\ in \\W = a L^b\\, a vector of two or
  an array `[n_pop x n_regions x n_sexes x 2]`. Required when
  `waa_model = "wt_len"`.

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
  recognized:

  `ln_M`

  :   Array of log-scale starting values for natural mortality,
      dimensioned
      `[n_popblks × n_regionblks × n_yearblks × n_ageblks × n_sexblks]`.
      Defaults to `log(0.5)` for all blocks if not supplied.

  `ln_growth_pars`

  :   Array of log-scale starting values for the growth parameters,
      dimensioned `[n_pop × n_regions × n_sexes × n_gpars]` in the order
      `L1, L2, K, CV1, CV2` and, under the Richards form, `rho`.
      Defaults to the ends of the length bins with a rate of `0.15` and
      CVs of `0.1`, so supply your own for any real model.

  `growth_pe_pars`

  :   Array of process error starting values for both growth deviation
      streams, dimensioned
      `[n_pop × n_regions × max(4, n_ages, n_gpars) × n_sexes × 2]`. The
      first stream holds one log sigma per growth parameter for the
      time-varying deviations; the second holds the semi-parametric
      surface's correlations by age, year and cohort in slots one to
      three and a log scale in slot four for the correlated forms, or
      one log sigma per age for `"iid"` and `"rw"`. Defaults to
      `log(0.1)` for the first stream and `log(0.05)` with correlations
      of `0.3` for the second. Slots a form does not read are mapped
      off.

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
