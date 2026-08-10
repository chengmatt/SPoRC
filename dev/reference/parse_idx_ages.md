# Build an index age-selection array

Turns a per-fleet specification of which ages contribute to an index
into the `[age, fleet]` array of 0/1 weights the objective function
uses. Accepts a list with one element per fleet, where each element is a
vector of ages or `NULL` for all ages, or an array already in
`[age, fleet]` form.

## Usage

``` r
parse_idx_ages(idx_ages, n_ages, n_fleets, what)
```

## Arguments

- idx_ages:

  List, array, or `NULL`. Per-fleet age selection.

- n_ages:

  Integer. Number of model ages.

- n_fleets:

  Integer. Number of fleets.

- what:

  Character. Name used in error messages.

## Value

Array `[n_ages x n_fleets]` of 0/1 weights.
