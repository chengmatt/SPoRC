# Recruitment and initial age deviation penalties

Population/region-specific IID penalties on initial age deviations
(`ln_InitDevs`) and on early/late recruitment deviations (`ln_RecDevs`),
including the Methot & Taylor bias-ramp adjustment. Called once from the
"Recruitment (Penalty)" section of `SPoRC_rtmb.R`.

## Usage

``` r
get_recruitment_penalty(
  n_pop,
  n_regions,
  n_ages,
  n_est_rec_devs,
  rec_dd,
  natal_region,
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

- n_pop, n_regions, n_ages, n_est_rec_devs:

  Dimension sizes.

- rec_dd:

  Integer recruitment density-dependence switch (used only to pick
  `sigma_idx` when `n_pop == 1`).

- natal_region:

  Integer vector `[pop]` of natal region indices.

- rec_region_prop_spec:

  Integer switch; when `1`, populations/regions with a fixed zero
  recruitment proportion are skipped.

- rec_region_prop:

  Array `[pop, region]` of recruitment regional apportionment.

- equil_init_age_strc:

  Integer switch selecting which initial age deviations are penalized
  (`1`: all but plus group, `2`: all, `3`: shared subset).

- ln_InitDevs:

  Array `[pop, region, age, sex]` of initial age deviations. A 3-D
  `[pop, region, age]` array (the layout before the sex dimension
  existed) is accepted and treated as one shared curve.

- init_age_devs_shared:

  Integer vector of shared initial age deviation indices (used when
  `equil_init_age_strc == 3`).

- ln_sigmaR:

  Array `[early/late, pop, region]` of log-sigma for recruitment
  deviations.

- bias_ramp:

  Numeric vector `[year]` of bias-ramp adjustment factors.

- sigmaR_switch:

  Integer year index at which recruitment deviations switch from the
  early to the late sigma regime.

- ln_RecDevs:

  Array `[pop, region, year]` of recruitment deviations.

- sigmaR2_early, sigmaR2_late:

  Arrays `[pop, region]` of squared sigma used for the bias-ramp mean
  offset.

- do_rec_bias_ramp:

  Integer switch enabling the bias-ramp log-sigma adjustment.

- map_ln_RecDevs:

  Array `[pop, region, year]` mirroring `map$ln_RecDevs`; cells that are
  `NA` are fixed rather than estimated and are left unpenalized. `NULL`
  penalizes every cell.

- init_devs_pen_use:

  Array of 0/1 matching `ln_InitDevs`, naming which cells are penalized.
  Sexes sharing one parameter keep only the first sex's copy flagged so
  the shared parameter is not penalized twice; sex-specific deviations
  flag every sex. `NULL` penalizes only the first sex's slice, which is
  the pre-sex-dimension behaviour.

- Use_init_sex_pen:

  Integer (0/1). Whether each later sex's initial age deviations are
  tied to the first sex's through a Gaussian on their difference at
  every penalized age. Only meaningful when the sexes carry their own
  curves.

- ln_sigma_init_sex:

  Log standard deviation of that tie.

- init_bias_ramp:

  Numeric vector of length `n_ages - 1`, the bias ramp read at the year
  each initial age was born (deviation index `1 - age`). `NULL` uses the
  first model year's ramp value for every age, the previous behaviour.

- map_ln_InitDevs:

  Numeric array matching `ln_InitDevs`, the map levels (`NA` where
  fixed). Cells sharing a level hold one parameter and split one penalty
  between them. `NULL` penalizes every cell.

## Value

List with elements `Init_Rec_nLL` (array `[pop, region, age, sex]`),
`Init_Sex_nLL` (the same layout, the between-sex tie, zero for the first
sex and whenever the tie is off) and `Rec_nLL` (array
`[pop, region, year]`), each holding negative log-likelihood penalties
(0 where not penalized).
