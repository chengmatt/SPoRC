# Evaluate one composition data source from its registered OSA observations

Runs
[`eval_comp_osa`](https://chengmatt.github.io/SPoRC/dev/reference/eval_comp_osa.md)
over the discrete fleets and then the continuous ones, on the vectors
[`pack_comp_source_osa`](https://chengmatt.github.io/SPoRC/dev/reference/pack_comp_source_osa.md)
built and the call site registered. The discrete pass zeroes the
container it fills; the continuous pass only does so when there were no
discrete fleets, so the two families add rather than overwrite.

## Usage

``` r
eval_comp_source_osa(
  nLL_arr,
  tracked_discrete,
  tracked_continuous,
  ObsArr,
  ExpArr,
  UseArr,
  TypeMat,
  LikeTypeVec,
  ISSArr,
  lnThetaArr,
  lnThetaAggVec,
  LNcorrArr,
  LNcorrAggVec,
  age_or_len,
  n_model_bins,
  comp_bins_spec,
  AgeingErrorArr = NULL,
  LenBinMap_fn = function(y, f) NA,
  n_pop,
  n_regions,
  n_yrs,
  n_seas,
  n_fleets,
  n_sexes,
  pop = FALSE,
  addtocomp = 0,
  seas_agg = 0
)
```

## Arguments

- nLL_arr:

  Container for this data source's negative log likelihood, dimensioned
  region by year by season by sex by fleet, with a leading population
  dimension when `pop` is `TRUE`.

- tracked_discrete, tracked_continuous:

  Registered observation vectors, or `NULL` where no fleet uses that
  family.

- ObsArr:

  Observed compositions.

- ExpArr:

  Predicted numbers the compositions are formed from: catch at age or
  length, discards at age or length, or survey available numbers.

- UseArr:

  Integer array flagging which cells are fit.

- TypeMat:

  Integer codes naming how each year and fleet's compositions are split,
  over year by fleet.

- LikeTypeVec:

  Integer likelihood per fleet.

- ISSArr:

  Input sample size.

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

- LenBinMap_fn:

  Function `(y, f)` returning the model bin to observed bin map, read
  only for length compositions. The map varies by neither year nor fleet
  and takes both so it is read the same way `AgeingErrorArr` is. `NA`
  leaves the observed bins as the model bins.

- n_pop, n_regions, n_yrs, n_seas, n_fleets, n_sexes:

  Model dimensions.

- pop:

  Logical. `TRUE` for the population-specific data source, whose arrays
  have a leading population dimension and are never summed over
  populations.

- addtocomp:

  Small constant added to a composition.

- seas_agg:

  Integer vector, one per fleet. `1` builds the predicted composition
  from every season of the year summed together, `0` from the season the
  observation sits in.

## Value

`nLL_arr`, filled in.
