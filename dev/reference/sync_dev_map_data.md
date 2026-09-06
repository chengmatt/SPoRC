# Refresh the map mirrors kept in the data list

Several deviation penalties key on a copy of their parameter's factor
map kept in the data list under `map_<par>`, since the map itself is
applied by `RTMB::MakeADFun` and is invisible inside the objective.
Those copies are written at setup, so a map edited by hand afterwards
would otherwise leave the penalty evaluating deviations that are no
longer estimated. Rebuilding the mirrors from the map immediately before
the model is constructed keeps the two in step, with the map treated as
authoritative.

## Usage

``` r
sync_dev_map_data(data, mapping)
```

## Arguments

- data:

  Named list of model data, as passed to `RTMB::MakeADFun`.

- mapping:

  Named list of factor maps, as passed to `RTMB::MakeADFun`.

## Value

`data` with every `map_<par>` element refreshed from `mapping[[par]]`.

## Details

A mirror whose parameter has no entry in `mapping`, or whose length no
longer matches (as when a caller has truncated one but not the other),
is left untouched.
