# Seasons a fit holds its year totals in

For each data source a fleet reports once a year, the season its `Use`
array turns on in each year, which is where the operating model draws
that year's total (`seas_agg_slot` in
[`Setup_Sim_Fishing`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md)
and
[`Setup_Sim_Survey`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md)).
A year with no observation takes season one, and years past the fit
repeat the last fitted year.

## Usage

``` r
seas_agg_slot_list(data, n_yrs, platform)
```

## Arguments

- data:

  Data list of the fitted model.

- n_yrs:

  Number of years the operating model runs.

- platform:

  `"fish"` or `"srv"`.

## Value

Named list of integer matrices `[n_yrs, n_fleets]`, one for each data
source with a fleet reporting once a year.
