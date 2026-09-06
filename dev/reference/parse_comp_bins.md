# Parse and report the observed bins a composition data source is fitted over

Wraps
[`parse_bin_subset`](https://chengmatt.github.io/SPoRC/dev/reference/parse_bin_subset.md)
and records which fleets ended up restricted, so the setup messages name
the bins a data source is fitted over rather than leaving it implicit.
Used for every `*_bins` argument, age and length, retained and
discarded, marginal and conditional, so a restriction reads the same way
whichever data source it was set on.

## Usage

``` r
parse_comp_bins(bins, n_bins, n_fleets, what)
```

## Arguments

- bins:

  List, array, or `NULL`. Per-fleet bin selection.

- n_bins:

  Integer. Number of observed bins the data source is recorded on, that
  is after any ageing error or length bin map.

- n_fleets:

  Integer. Number of fleets.

- what:

  Character. Argument name, used in messages and errors.

## Value

Array `[n_bins x n_fleets]` of 0/1 weights.
