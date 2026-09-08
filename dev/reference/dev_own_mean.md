# A deviation series' own weighted mean

The center a penalty takes when the level of the series is left to the
rest of the model rather than fixed at the bias-corrected mean.
Penalizing about it constrains only the spread, which is what a sum of
squares about the series' own mean amounts to. Fewer than two penalized
cells leaves no spread to measure, so the center falls back to zero.

## Usage

``` r
dev_own_mean(devs, wt)
```

## Arguments

- devs:

  Vector of deviations, on the log scale.

- wt:

  Numeric vector the same length, zero where the cell is not penalized.

## Value

The weighted mean, or `0` when fewer than two cells are penalized.
