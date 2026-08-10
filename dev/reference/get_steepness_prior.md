# Scaled beta prior on steepness

Called once from the "Steepness (Prior)" section of `SPoRC_rtmb.R`.

## Usage

``` r
get_steepness_prior(h_prior, h_trans)
```

## Arguments

- h_prior:

  Data frame with optional columns `lb` and `ub` giving the beta's
  support (defaulting to `0.2` and `1`), and columns `pop`, `region`,
  `mu` (prior mean steepness, natural scale), `sd` (prior SD, natural
  scale) with one row per penalized parameter.

- h_trans:

  Array `[pop, region]` of steepness on its transformed (0.2, 1) scale.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `h_prior`.
