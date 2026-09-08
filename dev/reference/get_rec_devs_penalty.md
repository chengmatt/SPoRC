# Recruitment deviation penalties

Population and region specific penalties on the recruitment deviations
(`ln_RecDevs`), independent or as a process over time. Called from
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

  Integer switch; when `1`, populations and regions with a fixed zero
  recruitment proportion are skipped.

- rec_region_prop:

  Array `[pop, region]` of recruitment regional apportionment.

- ln_sigmaR:

  Array `[early/late, pop, region]` of log sigma.

- bias_ramp:

  Numeric vector `[year]` of bias ramp adjustment factors.

- sigmaR_switch:

  Integer year index at which the deviations switch from the early to
  the late sigma regime.

- ln_RecDevs:

  Array `[pop, region, year]` of recruitment deviations.

- sigmaR2_early, sigmaR2_late:

  Arrays `[pop, region]` of squared sigma used for the bias-corrected
  mean.

- do_rec_bias_ramp:

  Integer switch enabling the bias ramp log sigma term.

- map_ln_RecDevs:

  Array `[pop, region, year]` mirroring `map$ln_RecDevs`. Cells that are
  `NA` are fixed rather than estimated and go unpenalized; cells sharing
  a level split one penalty. `NULL` penalizes every cell in full.

- RecDevs_model:

  Integer process error structure: `1` independent, `2` random walk, `3`
  AR1.

- RecDevs_rho:

  Array `[pop, region]` of unconstrained AR1 correlations, transformed
  to \\(-1, 1)\\ here. Read when `RecDevs_model = 3`.

- RecDevs_rw_init_sigma:

  Standard deviation given to year one of a random walk. Default `5`,
  which leaves the level of the series effectively free. `NA` starts the
  walk at zero under its own sigma. Read when `RecDevs_model = 2`.

- RecDevs_pen_center:

  Integer. `1` centers on the deviations' own weighted mean, `0` on the
  bias-corrected mean. Read under `RecDevs_model = 1` only.

## Value

Array `[pop, region, year]` of negative log-likelihood penalties, zero
where nothing is penalized.

## Details

Under `RecDevs_model = 1` the deviations are independent and split into
an early and a late sigma regime at `sigmaR_switch`. Both regimes are
the same penalty read at a different sigma over different years, \\-\log
\phi(\varepsilon_y \mid \mu_y, \sigma\_{R,k})\\ with an extra \\-(1 -
b_y / 2)\log \sigma\_{R,k}\\ when the bias ramp is on, where

- \\\varepsilon_y\\ is the deviation in year \\y\\, log scale,
  estimated.

- \\k\\ is 1 for years before `sigmaR_switch` and 2 from it on.

- \\\sigma\_{R,k} = \exp(\code{ln_sigmaR\[k,p,r\]})\\ is that regime's
  sigma, log scale.

- \\b_y\\ is the Methot and Taylor bias ramp in year \\y\\, between 0
  and 1, data.

- \\\mu_y = -\sigma\_{R,k}^2 b_y / 2\\ is the bias-corrected mean, or
  the deviations' own weighted mean over the regime's years.

The \\\log \sigma\\ term is what makes the sigma estimable: without it a
larger sigma always reduces the penalty. Under `RecDevs_model = 2` or
`3` the deviations are a walk or an AR1 process instead, each year
centered on the one before it, so neither the bias ramp nor the own-mean
center applies.
