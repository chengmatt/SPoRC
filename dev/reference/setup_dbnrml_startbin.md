# Validate the double normal start bin per fleet

Validate the double normal start bin per fleet

## Usage

``` r
setup_dbnrml_startbin(x, n_fleets, n_bins, what)
```

## Arguments

- x:

  `NULL` or an integer vector of start bins, one per fleet.

- n_fleets:

  Number of fleets.

- n_bins:

  Number of selectivity bins.

- what:

  Argument name for messages.

## Value

Integer vector `[n_fleets]`, ones when `x` is `NULL`.
