# Truncate Model Inputs for Retrospective Diagnostics

Removes the last `j` years from the model inputs and updates the data,
parameter arrays, maps, block structures and anything else dimensioned
by the number of years, so the result can be handed straight to the
model as one retrospective peel. Called by
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md).

## Usage

``` r
truncate_yr(j, data, parameters, mapping)
```

## Arguments

- j:

  Integer terminal years to remove. `0` returns the full dataset.

- data:

  List of model data supplied to the RTMB model.

- parameters:

  List of model parameters.

- mapping:

  List of parameter mappings used during estimation.

## Value

A list of `retro_data`, `retro_parameters` and `retro_mapping`, each
truncated to the shortened series.
