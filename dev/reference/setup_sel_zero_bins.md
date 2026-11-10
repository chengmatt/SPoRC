# Flag the bins a fleet does not select

Flag the bins a fleet does not select

## Usage

``` r
setup_sel_zero_bins(x, n_fleets, n_bins, what)
```

## Arguments

- x:

  `NULL` or a list with one element per fleet, each the bins whose
  selectivity is zero, or `NULL` for none.

- n_fleets:

  Number of fleets.

- n_bins:

  Number of selectivity bins.

- what:

  Argument name for messages.

## Value

0/1 matrix `[n_bins, n_fleets]`, or `NULL` when `x` is `NULL`.
