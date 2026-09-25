# Set up biological inputs for the estimation model

Sets weight-at-age, maturity-at-age, ageing error, the size-age
transition and any growth model, the numbers-at-age state, and the
natural mortality blocks and mapping. Call after
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
  M_seasblk_spec = "constant",
  M_sexblk_spec = "constant",
  Fixed_natmort = NULL,
  NAA_re = "none",
  NAA_re_ages = NULL,
  NAA_re_years = NULL,
  NAA_re_seasons = "annual",
  NAA_re_season = "iid",
  NAA_re_season_spec = "est_all",
  NAA_re_where = NULL,
  NAA_pe_spec = "est_all",
  NAA_sigma_spec = "est",
  NAA_re_region = "iid",
  NAA_re_region_spec = "est_all",
  NAA_re_pop = "iid",
  NAA_re_sex = "iid",
  NAA_sigma_popblk_spec = "constant",
  NAA_sigma_regionblk_spec = "constant",
  NAA_sigma_yearblk_spec = "constant",
  NAA_sigma_seasblk_spec = "constant",
  NAA_sigma_ageblk_spec = "constant",
  NAA_sigma_sexblk_spec = "constant",
  ...
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, `$map` and `$verbose`, as returned by
  [`Setup_Mod_Dim`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Dim.md).

- WAA:

  Spawning weight-at-age array
  `[n_pop × n_regions × n_years × n_seas × n_ages × n_sexes]`, also the
  fallback for `WAA_fish` and `WAA_srv`.

- WAA_fish:

  Fishery weight-at-age array, `WAA` with a trailing `n_fish_fleets`
  dim. `NULL` (default) reads `WAA` for every fleet.

- WAA_srv:

  Survey weight-at-age array, `WAA` with a trailing `n_srv_fleets` dim.
  `NULL` (default) reads `WAA` for every fleet.

- MatAA:

  Maturity-at-age array in \\\[0,1\]\\, dimensioned like `WAA`. Maturity
  at the first age must be exactly `0` under `rec_lag = 0`, so
  [`Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
  must have been called first.

- addtocomp:

  Deprecated here, pass it to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md).
  Still forwarded with a message. The constant added to composition
  proportions to avoid `log(0)`, ignored by the logistic normal.

- comp_const_obs:

  Deprecated here, pass it to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md).
  Still forwarded with a message. Integer switch for where `addtocomp`
  enters the multinomial: `1` adds it to the observed proportions that
  weight the likelihood as well as inside the logarithms, so the
  likelihood is stationary at `pred = obs`, `0` weights by the raw
  observed proportions.

- addtofishidx, addtosrvidx, addtotag:

  Deprecated here, pass them to
  [`Setup_Mod_Weighting`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Weighting.md).
  Still forwarded with a message. Constants added to the fishery
  indices, survey indices and tag recoveries.

- AgeingError:

  Ageing error array mapping true model ages onto observed age bins.
  Each row is one model age's share across the observed bins, summing to
  one, or to zero to drop that age. The age-axis twin of `LenBinMap`:
  the likelihood applies and validates the two identically. It sets
  which bins `ObsCatchAA`, `ObsDiscardAA` and `ObsSrvIdxAA` are
  dimensioned by; use the `*_bins` arguments to leave bins out of the
  likelihood instead. A `[n_model_ages × n_obs_ages]` matrix is
  time-invariant and expanded across years, a
  `[n_years × n_model_ages × n_obs_ages]` array is time-varying, and
  `NULL` (default) builds an identity matrix. For observed bins that are
  a subset of the model ages, supply a shifted identity such as
  `diag(1, n_model_ages)[, obs_age_index]`.

- AgeingError_fish:

  Optional per-fleet ageing error for the fishery fleets, either
  `[n_model_ages × n_obs_ages × n_fish_fleets]` or with a leading
  `n_years` dim, or `NULL` (default) to read the shared `AgeingError`.
  Each slice is validated as `AgeingError` is, and every fleet must land
  on the same observed bins. Read by a fleet's age compositions and its
  catch and discards at age.

- AgeingError_srv:

  As `AgeingError_fish` with `n_srv_fleets` in place of `n_fish_fleets`.
  Read by a fleet's age compositions and its index at age.

- Use_M_prior:

  Integer flag for a lognormal prior on natural mortality. `0` (default)
  or `1`.

- M_prior:

  Data frame of prior hyperparameters, one row per block combination,
  with columns `popblk`, `regionblk`, `yearblk`, `ageblk`, `sexblk`,
  `mu` on the natural scale, `sd`, and optionally `seasblk` (left out,
  it reads the block covering season one). Only used when
  `Use_M_prior = 1`.

- fit_lengths:

  Integer flag for fitting length compositions, `0` (default) or `1`.
  Requires a valid `SizeAgeTrans`.

- SizeAgeTrans:

  Size-at-age transition array
  `[n_pop × n_regions × n_years × n_seas × n_lens × n_ages × n_sexes]`,
  column-stochastic over ages. Required when `fit_lengths = 1`. Read by
  every fleet unless overridden.

- SizeAgeTrans_fish, SizeAgeTrans_srv:

  Optional per-fleet size-at-age arrays, dimensioned like `SizeAgeTrans`
  with a trailing fleet dim. `NULL` (default) reads the shared array.
  Only meaningful under `growth_model = "none"`; a growth model already
  derives one key per fleet at that fleet's timing and refuses these.

- do_caal:

  Integer flag for building the joint arrays at length and age, `0`
  (default) or `1`. Requires `fit_lengths = 1` and adds `Fish_caal`,
  `Fish_caal_discard` and `Srv_caal` to the report.

- growth_model:

  `"none"` (default) keeps `SizeAgeTrans` and the weight-at-age arrays
  as data. `"vb_schnute"` builds the size-age key from estimable von
  Bertalanffy parameters in Schnute's form: length `L1` at `growth_A1`,
  `L2` at `growth_A2`, rate `K`, and CVs `CV1` and `CV2` at the two
  reference ages. `"richards"` adds the coefficient `rho`, with
  `rho = 1` recovering von Bertalanffy. Both require `fit_lengths = 1`
  and ignore `SizeAgeTrans`.

- growth_spec:

  How the growth parameters are estimated: `"est_all"` (default, one set
  per population, region and sex), `"est_shared_r"`, `"est_shared_s"`,
  `"est_shared_r_s"`, or `"fix"`.

- growth_fix:

  Logical vector, one entry per growth parameter, naming which of L1,
  L2, K, CV1, CV2 and rho stay at their starting values whatever
  `growth_spec` says.

- growth_tv_model:

  Time variation of the growth parameters. `NULL` (default) holds every
  parameter constant. Otherwise a character vector of length `n_gpars`
  in parameter order, or named by parameter, each `"none"`, `"iid"`,
  `"rw"`, or `"dsem"`. A varying parameter gets a deviation series
  `ln_growth_devs` and a log sigma in the first data source of
  `growth_pe_pars`; under `"dsem"` the density comes from
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md)
  and that sigma stays at its start.

- growth_tv_years:

  Calendar years the deviations are active in. `NULL` (default) for
  every model year, a vector for every varying parameter, or a list
  named by parameter. Deviations outside the range are kept at zero.

- growth_tv_link:

  The scale a deviation enters on. `"log"` (default) multiplies the
  parameter by \\e^{\delta}\\; `"logit"` keeps it inside
  `growth_par_bounds`, so the parameter approaches a bound instead of
  crossing it.

- growth_par_bounds:

  Matrix `[n_gpars x 2]` of lower and upper bounds on the natural scale,
  required under the logit link.

- growth_tv_sigma_spec:

  `"fix"` (default) holds the process error sds of the deviations at
  their starting values, `"est"` estimates them. Both read the first
  data source of `growth_pe_pars`, one slot per growth parameter.

- growth_tv_spec:

  How the deviations are shared across strata, in the `growth_spec`
  vocabulary: `"est_all"` (default), `"est_shared_r"`, `"est_shared_s"`
  or `"est_shared_r_s"`.

- growth_tv_type:

  `"curve"` (default) reads every year's size at age off that year's
  curve. `"cohort"` advances size at age cohort by cohort: each cohort
  grows by the increment the current year's parameters imply from the
  size it reached, ages still in the linear phase keep their birth
  year's length at `growth_A1`, the first age past `growth_A1` is placed
  on the current year's curve, and the plus group blends the entering
  cohort with the fish already there by numbers at age. The CV at age
  stays at the first year's sizes, and propagation starts in the first
  year any deviation is active.

- growth_rw_init_sigma:

  Standard deviation given to the first year of a random walk on a
  growth parameter, as `srvsel_rw_init_sigma` for selectivity. Default
  `5`.

- growth_semipar:

  Semi-parametric growth: a year by age surface of deviations
  multiplying the parametric curve, so the deviations move mean length
  around it. `"none"` (default) keeps growth parametric; otherwise
  `"iid"`, `"rw"` over years within an age, `"3dmarg"` or `"3dcond"` (a
  Gaussian Markov random field over age, year and cohort on the marginal
  or conditional variance), `"2dar1"` (separable over ages and years),
  or `"dsem"` (density from
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md),
  one series per age, which refuses `growth_semipar_spec = "est"`). The
  spread at age follows the deviated mean, leaving the CV at age to the
  parametric part.

- growth_semipar_spec:

  Whether the second data source of `growth_pe_pars` is estimated
  (`"est"`) or kept at its starting values (`"fix"`, default). The
  deviations themselves are always estimated.

- growth_semipar_ages:

  Ages the deviations are estimated over, as ages rather than indices.
  `NULL` (default) uses every age; ages outside the set stay at zero.

- growth_semipar_years:

  Calendar years the deviations are estimated over. `NULL` (default)
  uses every year.

- LenBinMap:

  Optional matrix `[n_lens x n_obs_lens]` mapping the model's length
  bins onto the bins the compositions are recorded on, for compositions
  on coarser bins than the model has. Each row is one model bin's share
  across the observed bins, summing to one, or to zero to drop that bin.
  The length-axis twin of `AgeingError`, applied and validated
  identically. Use the `*LenComps_bins` arguments to leave bins out of
  the likelihood instead. `NULL` (default) fits on the model bins.

- growth_A1, growth_A2:

  Reference ages for `L1` and `L2`. `growth_A2 = "Linf"` makes `L2` the
  asymptotic length itself.

- growth_len_lower:

  Lower edges of the length bins. `lens` in `Setup_Mod_Dim` are
  midpoints; the key is built on the edges.

- growth_L0:

  Length at age zero anchoring the linear phase. Defaults to
  `growth_len_lower[1]`.

- growth_cv_type:

  `"len"` (default) interpolates the CV on mean length between `L1` and
  `L2`, `"age"` on age.

- growth_sd_type:

  `"cv"` (default) scales the mean by the CV parameters, `"sd"` reads
  them as standard deviations.

- growth_dist:

  `"normal"` (default) or `"lognormal"` length at age.

- growth_plus_group:

  `"mixture"` (default) takes the plus group's mean length as the
  survivorship-weighted mixture of the ages it holds, their numbers
  declining at an assumed 0.2 per year; `"curve"` reads the curve at the
  accumulator age.

- waa_model:

  Where weight at age comes from. `"data"` (default) reads `WAA`,
  `WAA_fish` and `WAA_srv`. `"wt_len"` builds them from the size-age key
  and \\W = a L^b\\ at the bin midpoints, so weight at age holds the
  spread of length at age; the spawning weight uses the key at spawning
  time and each fleet's weight the key at `t_fish` or `t_srv`. Under
  `"wt_len"` `WAA` may be `NULL`, and reference point and projection
  code still read `data$WAA`, so copy the reported arrays into the data
  list before calling them.

- wt_len_pars:

  Weight-length parameters \\a, b\\, a vector of two or an array
  `[n_pop x n_regions x n_sexes x 2]`. Required under
  `waa_model = "wt_len"`.

- M_spec:

  Natural mortality estimation. `"est_ln_M"` (default) estimates `ln_M`
  across the blocks; `"fix"` holds mortality at `Fixed_natmort` and maps
  `ln_M` off.

- M_popblk_spec:

  Blocking for `ln_M` across populations, either `"constant"` (default)
  or a list of integer index vectors, e.g. `list(1, 2)`.

- M_ageblk_spec:

  Blocking across ages, `"constant"` (default) or a list of integer
  index vectors, e.g. `list(1:5, 6:10)`.

- M_regionblk_spec:

  Blocking across regions, `"constant"` (default) or a list of integer
  index vectors, e.g. `list(1:3, 4:5)`.

- M_yearblk_spec:

  Blocking across years, `"constant"` (default) or a list of integer
  index vectors, e.g. `list(1:10, 11:30)`.

- M_seasblk_spec:

  Blocking across seasons, `"constant"` (default, one rate all year) or
  a list of integer index vectors, e.g. `list(1:2, 3:4)`. Blocks hold
  rates per year, so two half-year seasons at `0.2` and `0.4` accumulate
  an annual `0.3`. Only identifiable off within-year data (seasonal
  catch, seasonal comps, or surveys in more than one season), and even
  then the annual total comes back much better than the split, so prefer
  fixing the split and estimating the level. Warns for a single season
  model.

- M_sexblk_spec:

  Blocking across sexes, `"constant"` (default) or a list of integer
  index vectors, e.g. `list(1, 2)`.

- Fixed_natmort:

  Fixed natural mortality array, either
  `[n_pop × n_regions × n_years × n_ages × n_sexes]` or the same with
  `n_seas` between years and ages; the 5d form is expanded across
  seasons. Values are rates per year either way, so mortality in a
  season is the rate times `seasdur`. Required when `M_spec = "fix"`.

- NAA_re:

  State-space numbers at age: the log numbers become parameters for ages
  two and older including the plus group, and the deterministic
  mortality and ageing step becomes the prediction they are penalized
  against. `"none"` (default) keeps numbers deterministic; otherwise
  `"iid"`, `"1dar1_a"` over ages, `"1dar1_y"` over years, `"2dar1"`
  separable over both, `"3dcond"` or `"3dmarg"` over age, year and
  cohort, or `"dsem"`, which takes the density from
  [`Setup_Mod_DSEM`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_DSEM.md)
  one series per state age, leaves `ln_sigmaNAA` at its start and
  refuses `NAA_sigma_spec = "est"`. Each series is the log state with
  the log deterministic prediction as its mean. Age one belongs to
  `ln_RecDevs` and year one at ages two and older to `ln_InitDevs`, so
  the three partition the numbers at age. Which cells are estimated is
  set by `map$ln_NAA` and `data$n_est_naa_re`, never by `dim(ln_NAA)`.
  The state covers the assessment years only, so a forecast from
  [`Do_Population_Projection`](https://chengmatt.github.io/SPoRC/dev/reference/Do_Population_Projection.md)
  omits this process error while the closed loop operating model
  projects the state forward.

- NAA_re_ages:

  Ages the state is estimated over, matched against
  `input_list$data$ages` by value. A model on ages `0:4` takes
  `c(1, 2, 3, 4)` for the full state. `NULL` (default) uses `ages[-1]`.
  Must be a contiguous run.

- NAA_re_years:

  Calendar years the state is estimated over, matched against
  `input_list$data$years` by value. `NULL` (default) uses `years[-1]`.
  Must be a contiguous run.

- NAA_re_seasons:

  Seasons the state is estimated over. `"annual"` (default) puts a state
  at season one only, leaving the numbers within a year deterministic.
  `"all"` puts one at the start of every season, and an integer vector
  selects specific seasons, which need not be contiguous. Use it when
  only some seasons have observations, since a season with no data
  returns its prior as its posterior. The age, year and cohort
  correlations in `NAA_pe_pars` have no season dim, so every active
  season shares them within a population, region and sex; only the
  standard deviation varies by season, through `NAA_sigma_seasblk_spec`.

- NAA_re_season:

  Correlation across seasons within a year. `"iid"` (default) leaves the
  seasonal innovations independent; `"us"` estimates an unstructured
  correlation, \\n_k(n_k-1)/2\\ parameters over the \\n_k\\ active
  seasons. Needs more than one active season.

- NAA_re_season_spec:

  How the season correlations are shared, taking the same values as
  `NAA_re_region_spec`.

- NAA_re_where:

  Integer matrix `[population, region]`, `1` where the state runs and
  `0` where a population never occupies that region. `NULL` (default)
  gives every cell a state. A cell holding no fish has an undefined
  lognormal state, so `0` drops it from the map and the penalty; such
  cells need the region and population correlations off.

- NAA_pe_spec:

  How the age, year and cohort correlations in `NAA_pe_pars` are shared.
  `"est_all"` (default) gives a free set per population, region and sex.
  `"est_shared_p"`, `"est_shared_r"` and `"est_shared_s"` share one dim,
  `"est_shared_p_r"`, `"est_shared_p_s"` and `"est_shared_r_s"` two, and
  `"est_shared_p_r_s"` gives one set for the model. `"fix"` holds them
  at their starting values. Sharing a correlation is not correlating the
  innovations: regions that share \\\rho\\ still get independent shocks,
  whereas `NAA_re_region = "us"` makes the shocks covary.

- NAA_sigma_spec:

  Whether the process error standard deviations are estimated (`"est"`,
  default) or kept at their starting values (`"fix"`). The states
  themselves are always estimated.

- NAA_re_region:

  Correlation across regions, composed with the age and year grid.
  `"iid"` (default) leaves regions independent, `"us"` estimates an
  unstructured correlation of \\n_r(n_r-1)/2\\ parameters.

- NAA_re_region_spec:

  How the region correlations are shared. `"est_all"` (default) gives a
  free matrix per population and sex, `"est_shared_p"` and
  `"est_shared_s"` share one dim, `"est_shared_p_s"` gives a single
  matrix, and `"fix"` holds them.

- NAA_re_pop, NAA_re_sex:

  Correlation across populations and across sexes, composed with the
  region, age and year structures. `"iid"` (default) leaves them
  independent, `"us"` estimates an unstructured correlation. Both are
  global, so a two-sex model spends one parameter on
  `NAA_re_sex = "us"`.

- NAA_sigma_popblk_spec, NAA_sigma_regionblk_spec,
  NAA_sigma_yearblk_spec, NAA_sigma_seasblk_spec, NAA_sigma_ageblk_spec,
  NAA_sigma_sexblk_spec:

  Blocking for the process error standard deviation, each `"constant"`
  (default) or a list of integer vectors, exactly as the `M_*blk_spec`
  arguments. Blocking shares a standard deviation and never removes a
  cell from the state. Only `NAA_re = "iid"` admits one varying over
  years or ages; every other form is separable or Markov in a dim. The
  season dim is the exception, being whitened outside the age and year
  density, and is ruled out only by `NAA_re_season = "us"`.

- ...:

  Optional starting values by name. `ln_M` is dimensioned
  `[n_popblks × n_regionblks × n_yearblks × n_seasblks × n_ageblks × n_sexblks]`
  and defaults to `log(0.5)`; a 5d array from an older script works when
  there is one season block. `ln_growth_pars` is
  `[n_pop × n_regions × n_sexes × n_gpars]` in the order
  `L1, L2, K, CV1, CV2` and `rho`, defaulting to the ends of the length
  bins with a rate of `0.15` and CVs of `0.1`, so supply your own for
  any real model. `growth_pe_pars` is
  `[n_pop × n_regions × max(4, n_ages, n_gpars) × n_sexes × 2]`: the
  first data source holds one log sigma per growth parameter for the
  time-varying deviations, the second the semi-parametric surface's
  correlations by age, year and cohort in slots one to three with a log
  scale in slot four, or one log sigma per age under `"iid"` and `"rw"`.
  Slots a form does not read are mapped off. All `...` arguments are
  ignored when `M_spec = "fix"`.

## Value

`input_list` with `$data`, `$par` and `$map` updated, including
`$data$WAA`, `$data$WAA_fish`, `$data$WAA_srv`, `$data$MatAA`,
`$data$AgeingError`, `$data$M_blocks`, `$par$ln_M` and `$map$ln_M`.

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
