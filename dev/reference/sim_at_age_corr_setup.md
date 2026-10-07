# The correlation an operating model draws one at-age data source's residuals under

Stores the form and the unconstrained correlations under the names the
estimation model holds them, so a fit's values pass straight through and
both sides read them through `rho_trans` and `build_us_corr`.

## Usage

``` r
sim_at_age_corr_setup(
  sim_list,
  tag,
  corr,
  rho,
  rho_year,
  us_pars,
  n_fleets,
  pop = FALSE
)
```

## Arguments

- sim_list:

  Simulation list with `n_regions`, `n_sexes` and `n_obs_ages`.

- tag:

  `"catch"`, `"discard"` or `"srv_idx"`.

- corr:

  Form per fleet, `"iid"`, `"1dar1"`, `"us"` or `"2dar1"`, or their
  codes 0 to 3. `NULL` is `"iid"`.

- rho, rho_year:

  Unconstrained correlations across ages and across years,
  `[n_regions, n_sexes, n_fleets]`, or `NULL` for zero.

- us_pars:

  Unconstrained unstructured correlation parameters, pair by region by
  sex by fleet, or `NULL` for zero.

- n_fleets:

  Number of fleets of this data source.

- pop:

  Logical. `TRUE` for a population-specific data source, whose
  parameters have a leading population dim.

## Value

`sim_list` with `AgeObsCorr_<tag>`, `trans_rho_<tag>`,
`trans_rho_<tag>_year` and `trans_rho_<tag>_us`.
