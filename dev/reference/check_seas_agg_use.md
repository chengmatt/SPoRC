# Check that a seasonally aggregated data source has one observation per year

Under `"aggSeas"` the prediction is a year total, so more than one
season turned on in a region and year would fit that same total twice.

## Usage

``` r
check_seas_agg_use(use_arr, seas_agg, arg_name)
```

## Arguments

- use_arr:

  Use array with region, year, season and fleet in its last four dims. A
  population-specific array is allowed to have a leading dim.

- seas_agg:

  Integer vector from
  [`parse_seas_agg_spec`](https://chengmatt.github.io/SPoRC/dev/reference/parse_seas_agg_spec.md).

- arg_name:

  Name of the `Use` argument, used in error messages.

## Value

`NULL`, invisibly. Called for the error it raises.
