# Evaluate one composition data source through the direct likelihood

Walks every year, season, fleet and, for a population-specific data
source, every population in one composition data source, and calls
[`Get_Comp_Likelihoods`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Comp_Likelihoods.md)
on each cell that is fit. One call stands for what used to be written
out separately for retained fishery, discarded fishery and survey
compositions, for ages and for lengths, and again for the regional and
the population-specific data source of each.

## Usage

``` r
get_comp_source_nLL(
  nLL_arr,
  ObsArr,
  ExpArr,
  UseArr,
  ISSArr,
  WtArr,
  TypeMat,
  LikeTypeVec,
  lnThetaArr,
  lnThetaAggVec,
  LNcorrArr,
  LNcorrAggVec,
  age_or_len,
  n_model_bins,
  comp_bins_spec,
  AgeingErrorArr = NULL,
  LenBinMap_lik = NA,
  LenBinMap_fn = NULL,
  n_pop,
  n_regions,
  n_yrs,
  n_seas,
  n_fleets,
  n_sexes,
  pop = FALSE,
  addtocomp = 0,
  comp_const_obs = 1,
  do_internal_comp_osa = FALSE,
  tracked_discrete = NULL,
  tracked_continuous = NULL
)
```

## Arguments

- nLL_arr:

  Container for this data source's negative log likelihood, dimensioned
  region by year by season by sex by fleet, with a leading population
  dimension when `pop` is `TRUE`.

- ObsArr:

  Observed compositions.

- ExpArr:

  Predicted numbers the compositions are formed from: catch at age or
  length, discards at age or length, or survey available numbers.

- UseArr:

  Integer array flagging which cells are fit.

- ISSArr:

  Input sample size.

- WtArr:

  Likelihood weight applied to a multinomial.

- TypeMat:

  Integer codes naming how each year and fleet's compositions are split,
  over year by fleet.

- LikeTypeVec:

  Integer likelihood per fleet.

- lnThetaArr, lnThetaAggVec:

  Dirichlet-multinomial dispersion, for the split and the aggregated
  compositions.

- LNcorrArr, LNcorrAggVec:

  Logistic normal correlation parameters, for the split and the
  aggregated compositions.

- age_or_len:

  `0` for age compositions, `1` for length.

- n_model_bins:

  Number of model bins, `n_ages` or `n_lens`.

- comp_bins_spec:

  Array of 0/1 flags over bin by fleet naming the bins this data source
  is fit over. All ones means the whole range.

- AgeingErrorArr:

  Ageing error over year by model bin by observed bin by fleet, read
  only for age compositions.

- LenBinMap_lik:

  Model bin to observed bin map, read only for length compositions.

- LenBinMap_fn:

  Function of year and fleet returning the model bin to observed bin
  map, read only for length compositions on the OSA route.

- n_pop, n_regions, n_yrs, n_seas, n_fleets, n_sexes:

  Model dimensions.

- pop:

  Logical. `TRUE` for the population-specific data source, whose arrays
  have a leading population dimension and are never summed over
  populations.

- addtocomp:

  Small constant added to a composition.

- comp_const_obs:

  Constant the observations are scaled by.

- do_internal_comp_osa:

  Logical. `TRUE` hands the data source to
  [`eval_comp_source_osa`](https://chengmatt.github.io/SPoRC/dev/reference/eval_comp_source_osa.md),
  which reads the vectors
  [`pack_comp_source_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_source_osa.md)
  built and the call site registered.

- tracked_discrete, tracked_continuous:

  Registered observation vectors, read only on the OSA route and `NULL`
  where no fleet uses that family.

## Value

`nLL_arr`, filled in.
