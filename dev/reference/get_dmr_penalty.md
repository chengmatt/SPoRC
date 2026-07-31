# Discard mortality rate deviation penalty

IID penalty on every estimated discard mortality rate deviation, called
once from the "Discard Mortality Rate (Penalty)" section of
`SPoRC_rtmb.R`.

## Usage

``` r
get_dmr_penalty(
  logit_dmr_devs,
  ln_sigma_dmr,
  map_logit_dmr_devs,
  n_fish_fleets,
  n_yrs,
  n_regions,
  n_seas
)
```

## Arguments

- logit_dmr_devs:

  Array `[region, year, season, fish_fleet]` of discard mortality rate
  deviations on the logit scale.

- ln_sigma_dmr:

  Array `[region, season, fish_fleet]` of log-sigma for the deviation
  penalty.

- map_logit_dmr_devs:

  Array `[region, year, season, fish_fleet]` mirroring
  `$map$logit_dmr_devs`: an estimation index where a deviation is
  estimated, `NA` where it is fixed.

- n_fish_fleets, n_yrs, n_regions, n_seas:

  Dimension sizes.

## Value

Array `[region, year, season, fish_fleet]` of negative log-likelihood
penalties (0 where the deviation is not estimated).

## Details

The penalized set is read off `map_logit_dmr_devs` rather than
recomputed, so it always matches what is estimated, including deviations
mapped off by hand after setup. By default
[`do_dmr_dev_mapping`](https://chengmatt.github.io/SPoRC/dev/reference/do_dmr_dev_mapping.md)
estimates a deviation in every cell the objective does not treat as a
true closure. Discard observations are not the boundary: `dmr` is
identified through total mortality (`ZAA`) wherever a cell is fished and
retention is less than one, and it cancels out of `PredDiscard` itself.
