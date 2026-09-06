# Normal prior on natural mortality

Called once from the "Natural Mortality (Prior)" section of
`SPoRC_rtmb.R`.

## Usage

``` r
get_natmort_prior(M_prior, ln_M, M_blocks)
```

## Arguments

- M_prior:

  Data frame with columns `popblk`, `regionblk`, `yearblk`, `ageblk`,
  `sexblk` (block indices into `M_blocks`), `mu` (prior mean, natural
  scale), `sd` (prior SD, log scale), one row per penalized parameter.
  An optional `seasblk` column names the season block the prior applies
  to; without it the prior reads the block covering the first season,
  which is every season when mortality does not vary within the year.

- ln_M:

  Vector of estimated log natural mortality values, indexed by
  `M_blocks`.

- M_blocks:

  Array `[pop, region, year, season, age, sex]` mapping each
  population/region/year/season/age/sex cell to an index into `ln_M`.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `M_prior`.
