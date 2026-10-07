# Store the season each fleet's year total is drawn into

Checks and stores the `seas_agg_slot` entries a
[`Setup_Sim_Fishing`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Fishing.md)
or
[`Setup_Sim_Survey`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Sim_Survey.md)
call was given, leaving any from the other call in place.

## Usage

``` r
store_seas_agg_slot(sim_list, seas_agg_slot, n_fleets)
```

## Arguments

- sim_list:

  Simulation list with `n_yrs` and `n_seas`.

- seas_agg_slot:

  Named list of integer matrices `[n_yrs, n_fleets]`, or `NULL`.

- n_fleets:

  Number of fleets of these data sources.

## Value

`sim_list` with `seas_agg_slot` updated.
