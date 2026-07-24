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
  do_rec_bias_ramp
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

  Array `[pop, region, age]` of initial age deviations.

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

## Value

List with elements `Init_Rec_nLL` (array `[pop, region, age]`) and
`Rec_nLL` (array `[pop, region, year]`), each holding negative
log-likelihood penalties (0 where not penalized).
