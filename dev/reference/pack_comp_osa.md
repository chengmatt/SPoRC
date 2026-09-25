# Pack observed composition data into a single flat OBS vector (OSA)

Builds the flat tracked OBS vector
[`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html)
needs, ordered region-fastest with population outermost so the evaluator
can use strided indexing over one continuous vector.

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

  Optional weighting for the multinomial.

- UseArr:

  Region-use flags.

- TypeMat:

  Composition type matrix (0, 1, 2).

- LikeTypeVec:

  Likelihood type per fleet.

- n_yrs, n_seas, n_fleets, n_sexes, n_pop:

  Model dimensions.

- addtocomp:

  Small constant added to the proportions before normalization.

- family:

  `"discrete"` or `"continuous"`.

- pop:

  Logical; `TRUE` treats the population dim as the outermost layer.

- return_labels:

  Logical; `TRUE` also builds a per-element label data frame giving the
  pop, region, year, season, fleet, sex, bin, comp_type,
  likelihood_type, family and last_in_group of every entry, in the same
  order, for relabeling
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  residuals afterwards (see \[get_osa()\]). Left `FALSE` (default)
  inside the model to avoid the extra tracking cost.

- BinsArr:

  Optional `[n_obs_bins x n_fleets]` 0/1 array naming the observed bins
  each fleet is fitted over, or `NULL` (default) for all bins. A
  restricted fleet packs a shorter block, and `eval_comp_osa` must be
  handed the same array so its strides stay in step.

## Value

The flat OBS vector, or, under `return_labels = TRUE`, a list of `vec`
and `labels`. `NULL` when no fleet of this family is present.

## Details

The discrete families pack counts: the multinomial as
`round(prop x ISS x Wt)` and the Dirichlet-multinomial as
`round(prop x ISS)`. The logistic-normal families pack the additive log
ratio of the observation, which has to happen here because a tracked OBS
vector cannot be changed later: proportions take `+addtocomp`, are
renormalized, and become `log(p_k / p_K)` for k = 1..K-1, with the last
bin the reference and dropped. A block is then `n_obs_bins - 1` long
under comp type 0, `n_ru x (n_obs_bins - 1) x n_sexes` under type 1, and
`n_obs_bins * n_sexes - 1` under type 2, which takes one joint reference
for the whole bin by sex stack.
