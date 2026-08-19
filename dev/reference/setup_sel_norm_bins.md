# Build the selectivity standardization window

Records which bins the mean-one standardization averages over for each
fleet. Every bin is the default; a fleet whose catchability is defined
against only part of the bin range standardizes over that part instead.

## Usage

``` r
setup_sel_norm_bins(input_list, sel_norm_bins, prefix, n_fleets, bins)
```

## Arguments

- input_list:

  Input list to append to

- sel_norm_bins:

  List with one element per fleet, or NULL for every bin

- prefix:

  One of "fish", "ret" or "srv"

- n_fleets:

  Number of fleets

- bins:

  Number of bins
