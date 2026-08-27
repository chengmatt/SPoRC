# Evaluate OSA composition negative log-likelihood from a flat tracked vector

Walks the same group order used by
[`pack_comp_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_osa.md),
keeping the pointer `k` synchronized with the packed slice lengths.
Evaluates the multinomial, Dirichlet-multinomial, or logistic-normal
likelihood for each region/sex/fleet/season block.

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

  Array receiving negative log-likelihood contributions.

- tracked:

  Flat tracked OBS vector.

- ExpArrFn:

  Function returning expected proportions for (p,y,seas,f).

- UseArr:

  Region-use flags.

- TypeMat:

  Composition type matrix.

- LikeTypeVec:

  Likelihood type per fleet.

- ISSArr:

  Input sample sizes.

- lnThetaArr:

  Log overdispersion.

- lnThetaAggVec:

  Aggregated log overdispersion.

- LNcorrArr:

  LN correlation parameters.

- LNcorrAggVec:

  Aggregated LN correlation parameters.

- n_regions:

  Total number of structural regions.

- n_yrs:

  Number of model years.

- n_seas:

  Number of seasons per year.

- n_fleets:

  Total number of fishing fleets.

- n_sexes:

  Number of biological sexes.

- n_model_bins:

  Number of internal model bins.

- n_obs_bins:

  Number of observational bins.

- age_or_len:

  Flag indicating age-based or length-based composition.

- AgeingErrorFn:

  Function `(y, f)` returning the ageing error matrix for a given year
  and fleet, or the length bin map, which ignores both. Fleet specific
  because a fishery and a survey need not read ages the same way.

- addtocomp:

  Small constant added to proportions before normalization.

- BinsArr:

  Optional `[n_obs_bins x n_fleets]` 0/1 array naming the observed bins
  each fleet is fitted over, or `NULL` (default) for all bins. Must be
  the same array handed to
  [`pack_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_osa.md),
  since the strides walked here are sized on it.

- family:

  Character string specifying the likelihood type, either "discrete" or
  "continuous".

- zero_init:

  Logical; whether to zero out the nLL array on entry.

- pop:

  Logical; if TRUE, evaluations account for the population structure
  layer.

- n_pop:

  Number of population structures or pools.

## Value

Updated `nLL_arr` containing the evaluated negative log-likelihood
values.

## Details

Slice lengths must match the packer exactly:

Discrete (LikeType 0,1):

- Comp_Type 0: `n_fit_bins`

- Comp_Type 1/2: `n_ru x n_fit_bins x n_sexes`

Logistic-normal (LikeType 2,3,4):

- Comp_Type 0: `n_fit_bins - 1`

- Comp_Type 1: `n_ru x (n_fit_bins - 1) x n_sexes`

- Comp_Type 2: `n_ru x (n_fit_bins x n_sexes - 1)`

`n_fit_bins` is the number of bins the fleet is fitted over, taken from
`BinsArr` and equal to `n_obs_bins` when the fleet fits every bin.

These reduced lengths reflect that the tracked `Obs` vector is already
ALR-transformed (the last reference bin is dropped).
