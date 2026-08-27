# Evaluate CAAL likelihoods from a tracked OBS vector (OSA)

Walks the flat vector built by
[`pack_caal_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_caal_osa.md)
in the same order and evaluates each length bin's age composition
through
[`Get_Comp_Likelihoods_OSA`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Comp_Likelihoods_OSA.md),
so the observations stay on the tape as the tracked quantity
[`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
needs.

## Usage

``` r
eval_caal_osa(
  nLL_arr,
  tracked,
  ExpArrFn,
  UseArr,
  TypeMat,
  LikeTypeVec,
  ISSArr,
  lnThetaArr,
  lnThetaAggVec,
  n_regions,
  n_yrs,
  n_seas,
  n_lens,
  n_fleets,
  n_sexes,
  n_model_bins,
  n_obs_bins,
  AgeingErrorFn,
  addtocomp,
  BinsArr = NULL
)
```

## Arguments

- nLL_arr:

  Array \\\[region \times year \times season \times len \times sex
  \times fleet\]\\ receiving negative log-likelihood contributions.

- tracked:

  Flat tracked OBS vector from
  [`pack_caal_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_caal_osa.md).

- ExpArrFn:

  Function of `(y, seas, l, f)` returning the expected joint numbers at
  length and age for that group, indexed \\\[region \times age \times
  sex\]\\.

- UseArr:

  Use flags \\\[region \times year \times season \times len \times
  fleet\]\\.

- TypeMat:

  Composition type matrix \\\[year \times fleet\]\\.

- LikeTypeVec:

  Likelihood type per fleet.

- ISSArr:

  Input sample sizes \\\[region \times year \times season \times len
  \times sex \times fleet\]\\.

- lnThetaArr:

  Log overdispersion \\\[region \times sex \times fleet\]\\.

- lnThetaAggVec:

  Aggregated log overdispersion, one per fleet.

- n_regions, n_yrs, n_seas, n_lens, n_fleets, n_sexes:

  Dimension sizes.

- n_model_bins:

  Number of model age bins.

- n_obs_bins:

  Number of observed age bins.

- AgeingErrorFn:

  Function `(y, f)` returning the ageing error matrix for a given year
  and fleet.

- addtocomp:

  Small constant added to proportions before normalization.

- BinsArr:

  Optional `[n_obs_bins x n_fleets]` 0/1 array naming the observed age
  bins each fleet is fitted over, or `NULL` (default) for all bins. Must
  be the same array handed to
  [`pack_caal_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_caal_osa.md),
  since the strides walked here are sized on it.

## Value

Updated `nLL_arr`.
