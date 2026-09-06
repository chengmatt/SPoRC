# Dirichlet prior on movement rates

Called once from the "Movement Rates (Prior)" section of `SPoRC_rtmb.R`.

## Usage

``` r
get_movement_dirichlet_prior(Movement_prior, Movement, Mrate = NULL)
```

## Arguments

- Movement_prior:

  Data frame with columns `pop`, `region_from`, `year`, `seas`, `age`,
  `sex`, and `alpha` (list column of Dirichlet concentration vectors),
  one row per penalized movement-from vector. For CTMC movement, `seas`
  selects which season's generator supplies the annual fractions.

- Movement:

  Array `[pop, region_from, region_to, year, season, age, sex]` of
  movement rates.

- Mrate:

  Array of instantaneous movement rates (generator \\\dot{Q}\\),
  dimensioned like `Movement` and stored unscaled by season duration.
  `NULL` for unstructured movement, in which case `Movement` is used
  as-is.

## Value

Numeric scalar negative log-likelihood contribution, summed across all
rows of `Movement_prior`.

## Details

For CTMC movement (`move_type = 1`) the prior is placed on the *annual*
movement fractions \\\exp(\dot{\mathbf{Q}})\\, not on the seasonal
fractions \\\exp(\dot{\mathbf{Q}}\\\Delta t)\\ stored in `Movement`.
Those differ once `ctmc_scale_by_seasdur = 1`, and the difference is not
benign: a season's movement matrix approaches the identity as the season
shortens, so a fixed `alpha` silently becomes a much stronger constraint
as `n_seas` grows. On a three-region test setup the same `alpha = 3`
prior cost 1.04 nLL units at `n_seas = 1` but 9.91 at `n_seas = 12`.
Evaluating on the annual matrix makes `alpha` mean the same thing
regardless of seasonal structure, and matches how such priors are
elicited ("what fraction of fish move per year"). Under
`ctmc_scale_by_seasdur = 0` the two coincide, so nothing changes there.

Unstructured movement (`move_type = 0`) has no generator; its `Movement`
entries are transition fractions in their own right and are used
directly, unchanged.
