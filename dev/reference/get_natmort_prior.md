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
  `sexblk` (block numbers, read as
  `ln_M[popblk, regionblk, yearblk, seasblk, ageblk, sexblk]`), `mu`
  (prior mean, natural scale), `sd` (prior SD, log scale), one row per
  penalized parameter. An optional `seasblk` column names the season
  block the prior applies to; without it the prior reads the block
  covering the first season, which is every season when mortality does
  not vary within the year.

- ln_M:

  Array of log natural mortality, one value per block combination and
  dimensioned by the blocks. A five-dim array without seasons is read as
  one season block.

- M_blocks:

  Not read. The prior reads `ln_M` by block number.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `M_prior`.
