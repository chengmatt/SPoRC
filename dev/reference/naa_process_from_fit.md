# The fit's numbers-at-age process, as the operating model draws it

The sd by cell, read through the fit's sigma blocks; one correlation per
dim of the age-year field, a per-cell parameter averaged before it is
transformed; and the unstructured correlations across populations,
regions, sexes and seasons in the lower-triangle order
`draw_naa_innovations` reads. `simulation_self_test` draws the years
after `n_cond_yrs` from it and `condition_closed_loop_simulations` the
projection years.

## Usage

``` r
naa_process_from_fit(data, pars)
```

## Arguments

- data:

  Data list of the fit.

- pars:

  Parameter list at the fitted values.

## Value

List of `sigmaNAA` `[n_pop, n_regions, n_yrs, n_seas, n_ages, n_sexes]`
over the fitted years, `naa_rho`, the four margin switches `NAA_re_pop`,
`NAA_re_region`, `NAA_re_sex` and `NAA_re_season`, and their
correlations `naa_pop_corr`, `naa_region_corr`, `naa_sex_corr` and
`naa_season_corr`.
