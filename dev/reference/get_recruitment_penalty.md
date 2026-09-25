# Recruitment and initial age deviation penalties

The two deviation penalties the recruitment section owes, gathered so
the objective reads them in one place: the initial age deviations from
[`get_init_devs_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_init_devs_penalty.md)
and the recruitment deviations from
[`get_rec_devs_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_rec_devs_penalty.md).
Called once from the "Recruitment (Penalty)" section of `SPoRC_rtmb.R`.

## Usage

``` r
get_recruitment_penalty(
  n_pop,
  n_regions,
  n_ages,
  n_est_rec_devs,
  rec_region_prop_spec,
  rec_region_prop,
  equil_init_age_strc,
  ln_InitDevs,
  init_age_devs_shared,
  ln_sigmaR,
  bias_ramp,
  sigmaR_switch,
  ln_RecDevs,
  sigmaR2_early,
  sigmaR2_late,
  do_rec_bias_ramp,
  map_ln_RecDevs = NULL,
  RecDevs_model = 1,
  RecDevs_rho = NULL,
  RecDevs_rw_init_sigma = 5,
  RecDevs_pen_center = 0,
  InitDevs_pen_center = 0,
  init_devs_pen_use = NULL,
  Use_init_sex_pen = 0,
  ln_sigma_init_sex = 0,
  init_bias_ramp = NULL,
  map_ln_InitDevs = NULL
)
```

## Arguments

- n_pop, n_regions, n_ages:

  Dimension sizes.

- rec_region_prop_spec:

  Integer switch; when `1`, populations and regions with a fixed zero
  recruitment proportion are skipped.

- rec_region_prop:

  Array `[pop, region]` of recruitment regional apportionment.

- equil_init_age_strc:

  Integer switch naming which initial age deviations are penalized (`0`:
  none, `1`: all but the plus group, `2`: all, `3`: the shared subset).

- ln_InitDevs:

  Array `[pop, region, age, sex]` of initial age deviations. A 3-D
  `[pop, region, age]` array, the layout before the sex dim existed, is
  accepted and treated as one shared curve.

- init_age_devs_shared:

  Integer vector of shared initial age deviation indices, read when
  `equil_init_age_strc == 3`.

- ln_sigmaR:

  Array `[early/late, pop, region]` of log sigma. Initial ages read the
  early one.

- bias_ramp:

  Numeric vector `[year]` of bias ramp adjustment factors.

- sigmaR_switch:

  Integer year index the deviations switch from the early to the late
  sigma at.

- ln_RecDevs:

  Array `[pop, region, year]` of recruitment deviations.

- sigmaR2_early, sigmaR2_late:

  Arrays `[pop, region]` of squared sigma, used for the bias-corrected
  mean.

- do_rec_bias_ramp:

  Integer switch enabling the bias ramp log sigma term.

- map_ln_RecDevs:

  Array `[pop, region, year]` mirroring `map$ln_RecDevs`. `NA` cells are
  fixed rather than estimated and go unpenalized, and cells sharing a
  level split one penalty. `NULL` penalizes every cell in full.

- RecDevs_model:

  Integer process error structure: `1` independent, `2` random walk, `3`
  AR1.

- RecDevs_rho:

  Array `[pop, region]` of unconstrained AR1 correlations, transformed
  to \\(-1, 1)\\ here. Read under `RecDevs_model = 3`.

- RecDevs_rw_init_sigma:

  Standard deviation given to year one of a random walk. Default `5`,
  which leaves the level of the series effectively free; `NA` starts the
  walk at zero under its own sigma. Read under `RecDevs_model = 2`.

- RecDevs_pen_center:

  Integer. `1` centers on the deviations' own weighted mean, `0` on the
  bias-corrected mean. Read under `RecDevs_model = 1` only.

- InitDevs_pen_center:

  Integer. `1` centers on the deviations' own weighted mean, `0` on the
  bias-corrected mean.

- init_devs_pen_use:

  Array of 0/1 matching `ln_InitDevs`, naming which cells are penalized.
  Sexes sharing one parameter keep only the first sex's copy flagged so
  the shared parameter is not penalized twice; sex-specific deviations
  flag every sex. `NULL` penalizes only the first sex's slice, which is
  the pre-sex-dim behavior.

- Use_init_sex_pen:

  Integer (0/1). Whether each later sex's deviations are tied to the
  first sex's through a Gaussian on their difference at every penalized
  age. Only meaningful when the sexes have their own curves.

- ln_sigma_init_sex:

  Log standard deviation of that tie.

- init_bias_ramp:

  Numeric vector of length `n_ages - 1`, the bias ramp read at the year
  each initial age was born (deviation index `1 - age`). `NULL` reads
  the first model year's ramp value at every age.

- map_ln_InitDevs:

  Numeric array matching `ln_InitDevs` of map levels (`NA` where fixed).
  Cells sharing a level hold one parameter and split one penalty between
  them. `NULL` penalizes every cell in full.

## Value

List with `Init_Rec_nLL` and `Init_Sex_nLL` (arrays
`[pop, region, age, sex]`) and `Rec_nLL` (array `[pop, region, year]`),
each holding negative log-likelihood penalties and zero where nothing is
penalized.
