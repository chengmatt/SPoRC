# Pack observed composition data into a single flat OBS vector (OSA)

Produces the flat tracked OBS vector required by
[`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html).
Population is the outermost dimension, so the entire result is one
continuous vector with a single pointer.

## Usage

``` r
pack_comp_osa(
  ObsArr,
  ISSArr,
  WtArr,
  UseArr,
  TypeMat,
  LikeTypeVec,
  n_yrs,
  n_seas,
  n_fleets,
  n_sexes,
  addtocomp,
  family = "discrete",
  pop = FALSE,
  n_pop = 1,
  return_labels = FALSE,
  BinsArr = NULL
)
```

## Arguments

- ObsArr:

  Observed proportions or counts.

- ISSArr:

  Input sample sizes.

- WtArr:

  Optional weighting for multinomial.

- UseArr:

  Region-use flags.

- TypeMat:

  Composition type matrix (0,1,2).

- LikeTypeVec:

  Likelihood type per fleet.

- n_yrs:

  Number of model years.

- n_seas:

  Number of seasons per year.

- n_fleets:

  Total number of fishing fleets.

- n_sexes:

  Number of biological sexes.

- addtocomp:

  Small constant added to proportions before normalization.

- family:

  Character string specifying the likelihood type, either "discrete" or
  "continuous".

- pop:

  Logical; if TRUE, the population dimension is treated as the outermost
  layer.

- n_pop:

  Number of population structures or pools.

- return_labels:

  Logical; if TRUE, also builds a per-element label data.frame
  identifying the origin (pop, region, year, season, fleet, sex, bin,
  comp_type, likelihood_type, family, last_in_group) of every entry in
  the tracked vector, in the same order. Intended for post-hoc
  relabeling of
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  residuals (see \[get_osa()\]); left `FALSE` (default) inside the model
  itself to avoid the extra bookkeeping cost.

- BinsArr:

  Optional `[n_obs_bins x n_fleets]` 0/1 array naming the observed bins
  each fleet is fitted over, or `NULL` (default) for all bins.
  Restricted fleets pack a shorter block, and `eval_comp_osa` must be
  handed the same array so its strides stay in step with the packer.

## Value

If `return_labels = FALSE` (default): flat OBS vector, or `NULL` if no
fleet of this family is present (unchanged behavior). If
`return_labels = TRUE`: a list with elements `vec` (the flat OBS vector)
and `labels` (a data.frame with one row per element of `vec`), or `NULL`
if no fleet of this family is present.

## Details

Discrete families (LikeType 0,1):

- Multinomial (0): counts = round(prop x ISS x Wt)

- Dirichlet-multinomial (1): counts = round(prop x ISS)

Continuous families (LikeType 2,3,4): logistic-normal The ALR transform
is performed here, because the tracked OBS vector cannot be modified
later. Proportions receive `+addtocomp`, are renormalized, then
transformed to `log(p_k / p_K)` for k = 1..K-1. The last bin is the ALR
reference and is dropped:

- Comp_Type 0: length = `n_obs_bins - 1`

- Comp_Type 1: length = `n_ru x (n_obs_bins - 1) x n_sexes`

- Comp_Type 2: joint ALR of the full \[bin x sex\] stack -\> length =
  `n_obs_bins x n_sexes - 1` (Joint drops one reference for the whole
  stack -\> length `n_obs_bins * n_sexes - 1`)

The resulting vector is ordered region-fastest so that the likelihood
evaluator can use simple strided indexing.
