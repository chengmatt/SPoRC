# Recruitment deviation penalties

Penalties on `ln_RecDevs` by population and region, independent or as a
process over time. Called from
[`get_recruitment_penalty`](https://chengmatt.github.io/SPoRC/dev/reference/get_recruitment_penalty.md).

## Usage

``` r
get_rec_devs_penalty(
  n_pop,
  n_regions,
  n_est_rec_devs,
  rec_region_prop_spec,
  rec_region_prop,
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
  RecDevs_pen_center = 0
)
```

## Arguments

- n_pop, n_regions, n_est_rec_devs:

  Dimension sizes.

- rec_region_prop_spec:

  Integer switch; `1` skips populations and regions with a fixed zero
  recruitment proportion.

- rec_region_prop:

  Array `[pop, region]` of the regional apportionment.

- ln_sigmaR:

  Array `[early/late, pop, region]` of log sigma.

- bias_ramp:

  Numeric vector `[year]` of bias ramp factors, between 0 and 1.

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

## Value

Array `[pop, region, year]` of negative log-likelihood penalties, zero
where nothing is penalized.

## Details

Under `RecDevs_model = 1` the deviations are independent, on an early
and a late sigma either side of `sigmaR_switch`. Both regimes are the
same normal penalty read at a different sigma over different years, with
an extra log sigma term when the bias ramp is on; that term is what
makes the sigma estimable, since without it a larger sigma always
reduces the penalty. Under `RecDevs_model = 2` or `3` each year is
centered on the one before it, so neither the bias ramp nor the own-mean
center applies. The equations are in the model equations vignette.
