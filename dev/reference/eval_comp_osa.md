# Evaluate OSA composition negative log-likelihood from a flat tracked vector

Walks the group order
[`pack_comp_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_osa.md)
used, keeping the pointer `k` in step with the packed slice lengths, and
evaluates the multinomial, Dirichlet-multinomial or logistic-normal
likelihood for each region, sex, fleet and season block.

## Usage

``` r
eval_comp_osa(
  nLL_arr,
  tracked,
  ExpArrFn,
  UseArr,
  TypeMat,
  LikeTypeVec,
  ISSArr,
  lnThetaArr,
  lnThetaAggVec,
  LNcorrArr,
  LNcorrAggVec,
  n_regions,
  n_yrs,
  n_seas,
  n_fleets,
  n_sexes,
  n_model_bins,
  n_obs_bins,
  age_or_len,
  AgeingErrorFn,
  addtocomp,
  BinsArr = NULL,
  family = "discrete",
  zero_init = TRUE,
  pop = FALSE,
  n_pop = 1
)
```

## Arguments

- nLL_arr:

  Array receiving the negative log-likelihood contributions.

- tracked:

  Flat tracked OBS vector.

- ExpArrFn:

  Function returning the expected proportions for `(p, y, seas, f)`.

- UseArr:

  Region-use flags.

- TypeMat:

  Composition type matrix.

- LikeTypeVec:

  Likelihood type per fleet.

- ISSArr:

  Input sample sizes.

- lnThetaArr, lnThetaAggVec:

  Log overdispersion and its aggregated counterpart.

- LNcorrArr, LNcorrAggVec:

  Logistic-normal correlation parameters and their aggregated
  counterpart.

- n_regions, n_yrs, n_seas, n_fleets, n_sexes, n_pop:

  Model dimensions.

- n_model_bins, n_obs_bins:

  Numbers of model and observed bins.

- age_or_len:

  Flag for an age-based or length-based composition.

- AgeingErrorFn:

  Function `(y, f)` returning that year and fleet's ageing error matrix,
  or the length bin map, which ignores both. It is fleet specific
  because a fishery and a survey need not read ages the same way.

- addtocomp:

  Small constant added to the proportions before normalization.

- BinsArr:

  Optional `[n_obs_bins x n_fleets]` 0/1 array naming the observed bins
  each fleet is fitted over, or `NULL` (default) for all bins. Must be
  the array handed to
  [`pack_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_osa.md),
  since the strides walked here are sized on it.

- family:

  `"discrete"` or `"continuous"`.

- zero_init:

  Logical; whether the nLL array is zeroed on entry.

- pop:

  Logical; `TRUE` accounts for the population layer.

## Value

`nLL_arr` with the evaluated values.

## Details

The slice lengths have to match the packer exactly. A discrete family
takes `n_fit_bins` under comp type 0 and `n_ru x n_fit_bins x n_sexes`
under types 1 and 2. A logistic-normal family takes one fewer bin, since
the tracked `Obs` vector arrives already transformed with its reference
bin dropped: `n_fit_bins - 1` under type 0,
`n_ru x (n_fit_bins - 1) x n_sexes` under type 1, and
`n_ru x (n_fit_bins x n_sexes - 1)` under type 2. `n_fit_bins` comes
from `BinsArr` and equals `n_obs_bins` when the fleet fits every bin.
