# Store one at-age stream's observations, use flags and standard errors

Shared by the four streams and their population-specific forms. An
absent stream is given a zeroed array so the objective can index it
unconditionally.

## Usage

``` r
do_at_age_data_setup(
  input_list,
  obs,
  use,
  se,
  stream,
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

- stream:

  Stream tag naming the data elements, e.g. `"CatchAA"`.

- fleet_field:

  `"n_fish_fleets"` or `"n_srv_fleets"`.

- pop:

  Logical. `TRUE` for the population-specific stream.

## Value

`input_list` with the stream's three data elements set.
