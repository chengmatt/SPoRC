# Truncate Model Inputs for Retrospective Diagnostics

Internal helper used by
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/dev/reference/do_retrospective.md)
to truncate model inputs when conducting retrospective diagnostics. The
function removes the last `j` years from the terminal portion of the
time series and updates all associated data objects, parameter arrays,
and parameter mappings so that their dimensions remain internally
consistent.

## Usage

``` r
truncate_yr(j, data, parameters, mapping)
```

## Arguments

- j:

  Integer specifying the number of terminal years to remove from the
  dataset. A value of `0` returns the full dataset with no truncation.

- data:

  List containing model data supplied to the RTMB model.

- parameters:

  List containing model parameters supplied to the RTMB model.

- mapping:

  List defining parameter mappings used during estimation.

## Value

A list containing truncated versions of the RTMB inputs:

- `retro_data` – Modified data list with terminal years removed.

- `retro_parameters` – Parameter list truncated to match the shortened
  time series.

- `retro_mapping` – Mapping list updated to match truncated parameter
  dimensions.

## Details

Specifically, the function adjusts the model `data`, `parameters`, and
`mapping` lists used by the RTMB model by:

- Truncating the `years` vector.

- Removing terminal years from observations (catch, indices, and
  composition data).

- Truncating time-varying parameter arrays (e.g., recruitment
  deviations, fishing mortality deviations, selectivity deviations,
  movement parameters).

- Updating parameter mappings to match the truncated parameter
  dimensions.

- Adjusting block structures and auxiliary objects that depend on the
  number of modeled years.

The resulting objects can be passed directly to the model to fit a
retrospective peel.
