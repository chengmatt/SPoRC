# Season each fleet's year total of a data source is drawn into

The fit keeps a year total in whichever season its `Use` array names,
and
[`simulation_self_test`](https://chengmatt.github.io/SPoRC/dev/reference/simulation_self_test.md)
and
[`condition_closed_loop_simulations`](https://chengmatt.github.io/SPoRC/dev/reference/condition_closed_loop_simulations.md)
hand that season to the operating model as `seas_agg_slot`. Without it
the total goes in season one.

## Usage

``` r
seas_agg_season(sim_env, data_name, y, n_fleets)
```

## Arguments

- sim_env:

  Simulation environment.

- data_name:

  Data source, such as `"Catch"` or `"FishAgeComps_pop"`.

- y:

  Year.

- n_fleets:

  Number of fleets of the data source.

## Value

Integer vector, one season per fleet.
