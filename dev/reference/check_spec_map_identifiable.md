# Check that an observation-error spec leaves every parameter identified

An observation-error standard deviation needs more than one observation
to be estimable. Given one, the likelihood is unbounded: the standard
deviation collapses towards zero on whatever residual the model can fit
exactly and the `log(sigma)` term runs to negative infinity. The
optimiser reports convergence, so nothing about the fit announces the
problem.

## Usage

``` r
check_spec_map_identifiable(
  map,
  use,
  spec,
  dims,
  dim_abbrev,
  what = "parameter",
  min_obs = 2,
  warn_obs = 5
)
```

## Arguments

- map:

  Factor map returned by
  [`build_shared_spec_map`](https://chengmatt.github.io/SPoRC/dev/reference/build_shared_spec_map.md).

- use:

  Array of the same dimensions as the parameter, non-zero where an
  observation informs that cell.

- spec:

  The spec string, used in messages.

- dims:

  Named integer vector of array dimensions.

- dim_abbrev:

  Named character vector of dimension abbreviations.

- what:

  Character label for the parameter, used in messages.

- min_obs:

  Integer below which a group raises an error.

- warn_obs:

  Integer below which a group raises a warning.

## Value

`invisible(NULL)`. Called for its error and warning side effects.

## Details

This is not confined to `"est_all"`. Any spec that leaves a dimension
free when the corresponding observation array is one cell deep in the
other dimensions has the same failure, so which specs are safe depends
on the model's dimensions rather than on the spec string alone. The
check therefore counts the observations actually informing each
estimation group instead of rejecting particular spec names.
