# Observed and predicted conditional age-at-length proportions

A conditional age-at-length observation is an age composition within one
length bin, so each (year, season, length bin) row is normalized on its
own through the same
[`Restrc_Comps`](https://chengmatt.github.io/SPoRC/dev/reference/Restrc_Comps.md)
the marginal compositions use. The predicted rows come from the report's
joint arrays summed over populations, with the year's ageing error
applied, so the proportions here are the ones the likelihood fits.

## Usage

``` r
get_caal_prop(data, rep)
```

## Arguments

- data:

  Data list from the fitted model.

- rep:

  Report list from the fitted model.

## Value

A list with `Obs_Fish_caal`, `Pred_Fish_caal`, `Obs_Srv_caal` and
`Pred_Srv_caal`, each dimensioned by
`n_regions, n_years, n_seas, n_lens, n_obs_ages, n_sexes, n_fleets`, or
`NULL` for a fleet type carrying no conditional age-at-length data.
