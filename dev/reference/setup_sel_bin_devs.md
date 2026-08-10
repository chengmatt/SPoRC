# Set up bin-override selectivity deviations for one selectivity stream

Creates the bin-override deviation parameter array, its factor map, and
its process-error hyperparameters, and records which bins each fleet
overrides. Bins named here take a free annual selectivity value instead
of whatever the fleet's functional form produces, which lets an
otherwise parametric curve carry a handful of freely estimated bins.

## Usage

``` r
setup_sel_bin_devs(
  input_list,
  bin_dev_bins,
  pe_model,
  prefix,
  n_fleets,
  bins,
  starting_values = list()
)
```

## Arguments

- input_list:

  Named list with `$data`, `$par`, and `$map`.

- bin_dev_bins:

  List with one element per fleet, each a vector of bins to override or
  `NULL` for none, or `NULL` for no overrides anywhere.

- pe_model:

  Character vector `[n_fleets]` giving the process-error structure of
  the override deviations: `"none"`, `"iid"`, or `"rw"`.

- prefix:

  One of `"fish"`, `"ret"`, `"srv"`.

- n_fleets:

  Integer. Number of fleets in this stream.

- bins:

  Integer. Number of age or length bins.

- starting_values:

  Named list of user-supplied starting values.

## Value

`input_list` with the parameter, map, and data entries added.
