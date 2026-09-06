# Pack conventional-tag observations for OSA

Packs conventional-tag observations into a flat vector suitable for
one-step-ahead (OSA) analysis for a single likelihood family.

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

  Character, either `"count"` or `"comp"`.

- like_type:

  Integer likelihood type code (0-5).

- obs_recap:

  Observed recapture array.

- pred_recap:

  Predicted recapture array (used for scaling).

- tagged_fish:

  Array of numbers of tagged fish released.

- conv_tag_release_indicator:

  Matrix giving release region, year, season for each cohort.

- conv_tag_max_liberty:

  Maximum years at liberty to evaluate.

- n_conv_tag_cohorts:

  Number of conventional tag cohorts.

- n_yrs:

  Total number of modeled years.

- n_seas:

  Number of seasons per year.

- n_regions:

  Number of spatial regions.

- n_fish_fleets:

  Number of fishing fleets.

- n_pop_pool:

  Number of population pooling groups.

- n_age_pool:

  Number of age pooling groups.

- n_sex_pool:

  Number of sex pooling groups.

- pop_pool:

  List of population index pools.

- age_pool:

  List of age index pools.

- sex_pool:

  List of sex index pools.

- use_fish_tagging:

  Vector indicating which fleets use conventional tagging.

- conv_tag_mixing_period:

  Minimum seasons at liberty before tags are modeled.

- addtotag:

  Small constant added to avoid zeros.

- return_labels:

  Logical; if TRUE, also builds a per-element label data.frame
  identifying the origin (family, like_type, tag cohort/release
  region-year-season, recovery year/season, fleet, region, pop/age/sex
  pool, is_tail, last_in_group) of every entry in `vec`, in the same
  order. Intended for post-hoc relabeling of
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  residuals (see \[get_osa()\]); left `FALSE` (default) inside the model
  itself to avoid the extra residual tracking cost.

## Value

A list with components:

- `vec`: flat numeric/AD vector of packed observations, or `NULL` if no
  events.

- `grp_end`: integer vector of end indices for each composition group
  (empty for `family == "count"`).

- `lengths`: integer vector of per-group lengths.

- `labels`: data.frame with one row per element of `vec` (`NULL` unless
  `return_labels = TRUE`).

## Details

For `family == "count"`, each valid event contributes one integer
observation per `[fleet, pop_pool, region, age_pool, sex_pool]` cell
with `use_fish_tagging[f] == 1`, in `(f, p, r, a, s)` loop order within
each event, and events in `tag_grid` order. This mirrors
[`get_conv_tag_likelihoods()`](https://chengmatt.github.io/SPoRC/dev/reference/get_conv_tag_likelihoods.md)'s
exact accumulation structure: a separate Poisson/NB term is fit per
pool/age/sex cell and their -log-likelihoods summed, so packing one
count per `[region, fleet]` pre-summed across pools would evaluate a
different (non-equivalent) likelihood whenever more than one pool is
used.

For `family == "comp"`, each valid event contributes one composition
vector:

- release-conditioned (`like_type` 2, 4): recap cells in
  `(f, p, a, s, r)` loop order plus a non-recapture tail; counts are
  `round(prop * n_tags_released)`.

- recapture-conditioned (`like_type` 3, 5): recap cells only,
  conditioned on total recaptures; counts are
  `round(prop * n_tags_recap)`.
