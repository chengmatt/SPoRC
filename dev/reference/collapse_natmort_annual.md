# Collapse a seasonal natural mortality array to an annual rate

Mortality in a season is the rate times its duration, so the annual
total is the duration weighted sum over seasons. Read by anything
stepping a whole year at once: equilibrium seeds and the plus group
geometric series.

## Usage

``` r
collapse_natmort_annual(natmort, seasdur, seas_dim = 4)
```

## Arguments

- natmort:

  Array of rates with a season dim, either the full model grid or the
  smaller one the reference points use.

- seasdur:

  Season durations, summing to one.

- seas_dim:

  Position of the season dim.

## Value

`natmort` with seasons summed out. A constant rate comes back unchanged.
