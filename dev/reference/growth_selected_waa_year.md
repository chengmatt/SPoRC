# Selection-weighted weight at age for one year

A length-selective gear does not take fish evenly across an age, so the
mean weight of what it takes is the weight averaged over the fleet's key
re-weighted by selectivity, not the population mean. Applied only to the
fleets that ask for it. Used for fishery and survey fleets alike.

## Usage

``` r
growth_selected_waa_year(
  WAA_fleet,
  SizeAgeTrans_fleet,
  sel_l,
  wt_len_pars,
  len_mid,
  waa_selected,
  y,
  n_pop,
  n_regions,
  n_seas,
  n_sexes
)
```

## Arguments

- WAA_fleet:

  Array `[pop, region, year, season, age, sex, fleet]` of the fleet
  type's weight at age, returned with year `y` overwritten for the
  fleets named in `waa_selected`.

- SizeAgeTrans_fleet:

  Each fleet's size-age key.

- sel_l:

  Array `[region, year, len, sex, fleet]` of selectivity at length.

- wt_len_pars:

  Array `[pop, region, sex, 2]` of weight-length parameters.

- len_mid:

  Bin midpoints the weight-length relationship is read at.

- waa_selected:

  Integer vector `[fleet]` (0/1).

- y:

  Year index.

- n_pop, n_regions, n_seas, n_sexes:

  Dimensions.

## Value

`WAA_fleet` with year `y` updated.
