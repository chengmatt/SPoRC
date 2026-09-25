# Pack conventional-tag observations for OSA

Packs the conventional tag observations of one likelihood family into a
flat vector for one-step-ahead analysis.

## Usage

``` r
pack_tag_osa(
  family,
  like_type,
  obs_recap,
  pred_recap,
  tagged_fish,
  conv_tag_release_indicator,
  conv_tag_max_liberty,
  n_conv_tag_cohorts,
  n_yrs,
  n_seas,
  n_regions,
  n_fish_fleets,
  n_pop_pool,
  n_age_pool,
  n_sex_pool,
  pop_pool,
  age_pool,
  sex_pool,
  use_fish_tagging,
  conv_tag_mixing_period,
  addtotag,
  return_labels = FALSE
)
```

## Arguments

- family:

  Character, `"count"` or `"comp"`.

- like_type:

  Integer likelihood type code (0-5).

- obs_recap, pred_recap:

  Observed and predicted recapture arrays, the second used for scaling.

- tagged_fish:

  Numbers of tagged fish released.

- conv_tag_release_indicator:

  Matrix of the release region, year and season of each cohort.

- conv_tag_max_liberty:

  Maximum years at liberty to evaluate.

- n_conv_tag_cohorts, n_yrs, n_seas, n_regions, n_fish_fleets:

  Model dimensions.

- n_pop_pool, n_age_pool, n_sex_pool:

  Numbers of pooling groups.

- pop_pool, age_pool, sex_pool:

  Lists of the index pools themselves.

- use_fish_tagging:

  Vector flagging the fleets with tagging data.

- conv_tag_mixing_period:

  Minimum seasons at liberty before tags are modeled.

- addtotag:

  Small constant added to avoid zeros.

- return_labels:

  Logical; `TRUE` also builds a per-element label data frame giving the
  family, like_type, cohort release region, year and season, recovery
  year and season, fleet, region, pools, is_tail and last_in_group of
  every entry, in the same order, for relabeling
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  residuals afterwards (see \[get_osa()\]). Left `FALSE` (default)
  inside the model to avoid the extra tracking cost.

## Value

A list with `vec`, the flat vector of packed observations or `NULL` when
there are no events; `grp_end`, the end index of each composition group,
empty under `family == "count"`; `lengths`, the per-group lengths; and
`labels`, `NULL` unless `return_labels = TRUE`.

## Details

Under `family == "count"` each valid event contributes one integer per
`[fleet, pop_pool, region, age_pool, sex_pool]` cell with
`use_fish_tagging[f] == 1`, in `(f, p, r, a, s)` order within an event
and events in `tag_grid` order. That mirrors how
[`get_conv_tag_likelihoods()`](https://chengmatt.github.io/SPoRC/dev/reference/get_conv_tag_likelihoods.md)
accumulates, a separate Poisson or negative binomial term per cell
summed afterwards, so packing one count per region and fleet pre-summed
across pools would evaluate a different likelihood whenever more than
one pool is used.

Under `family == "comp"` each valid event contributes one composition
vector. The release-conditioned forms (`like_type` 2 and 4) take the
recapture cells in `(f, p, a, s, r)` order plus a non-recapture tail,
with counts `round(prop * n_tags_released)`; the recapture-conditioned
forms (3 and 5) take the recapture cells alone, with counts
`round(prop * n_tags_recap)`.
