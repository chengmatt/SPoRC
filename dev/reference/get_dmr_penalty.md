# Discard mortality rate deviation penalty

IID penalty on discard mortality rate deviations for
region-year-season-fleet cells with active discard data, called once
from the "Discard Mortality Rate (Penalty)" section of `SPoRC_rtmb.R`.

## Usage

``` r
get_dmr_penalty(
  logit_dmr_devs,
  ln_sigma_dmr,
  UseDiscard,
  UseDiscard_pop,
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

- UseDiscard:

  Array `[region, year, season, fish_fleet]` flagging aggregated discard
  observations in use.

- UseDiscard_pop:

  Array `[pop, region, year, season, fish_fleet]` flagging
  population-specific discard observations in use.

- n_fish_fleets, n_yrs, n_regions, n_seas:

  Dimension sizes.

## Value

Array `[region, year, season, fish_fleet]` of negative log-likelihood
penalties (0 where discard data are not in use).
