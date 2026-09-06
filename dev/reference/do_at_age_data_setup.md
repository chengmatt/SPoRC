# Store one at-age data source's observations, use flags and standard errors

Shared by the three data sources and their population-specific forms. An
absent data source is given a zeroed array so the objective can index it
unconditionally.

## Usage

``` r
do_at_age_data_setup(
  input_list,
  obs,
  use,
  se,
  data_source,
  fleet_field,
  pop = FALSE
)
```

## Arguments

- input_list:

  Named list with `$data`.

- obs, use, se:

  Observation, use and reported standard error arrays, any of them
  `NULL`.

- data_source:

  Data source tag naming the data elements, e.g. `"CatchAA"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- pop:

  Logical. `TRUE` for the population-specific data source.

## Value

`input_list` with the data source's three data elements set.
