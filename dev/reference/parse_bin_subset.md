# Build a per-fleet bin-selection array

Turns a per-fleet specification of which bins are used into the
`[bin, fleet]` array of 0/1 weights the objective function reads. Shared
by the ages that contribute to an index (`fish_idx_ages`,
`srv_idx_ages`) and by the observed bins a composition is fitted over
(the `*_bins` arguments), so both spellings behave identically. Accepts
a list with one element per fleet, where each element is a vector of bin
indices or `NULL` for all bins, or an array already in `[bin, fleet]`
form.

## Usage

``` r
parse_bin_subset(idx_bins, n_bins, n_fleets, what)
```

## Arguments

- idx_bins:

  List, array, or `NULL`. Per-fleet bin selection.

- n_bins:

  Integer. Number of bins the selection indexes into: model ages for the
  index arguments, observed composition bins for the `*_bins` arguments.

- n_fleets:

  Integer. Number of fleets.

- what:

  Character. Name used in error messages.

## Value

Array `[n_bins x n_fleets]` of 0/1 weights.
