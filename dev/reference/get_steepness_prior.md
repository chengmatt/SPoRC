# Scaled beta prior on steepness

Called once from the "Steepness (Prior)" section of `SPoRC_rtmb.R`.

## Usage

``` r
get_steepness_prior(h_prior, h_trans)
```

## Arguments

- h_prior:

  Data frame with columns `pop`, `region`, `mu` (prior mean steepness,
  natural scale), `sd` (prior SD, natural scale) — one row per penalized
  parameter.

- h_trans:

  Array `[pop, region]` of steepness on its transformed (0.2, 1) scale.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `h_prior`.
