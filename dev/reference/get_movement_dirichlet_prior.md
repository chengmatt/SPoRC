# Dirichlet prior on movement rates

Called once from the "Movement Rates (Prior)" section of `SPoRC_rtmb.R`.

## Usage

``` r
get_movement_dirichlet_prior(Movement_prior, Movement)
```

## Arguments

- Movement_prior:

  Data frame with columns `pop`, `region_from`, `year`, `seas`, `age`,
  `sex`, and `alpha` (list column of Dirichlet concentration vectors) —
  one row per penalized movement-from vector.

- Movement:

  Array `[pop, region_from, region_to, year, season, age, sex]` of
  movement rates.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `Movement_prior`.
