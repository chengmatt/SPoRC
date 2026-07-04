# Evaluate conventional-tag OSA likelihood from a tracked vector

Evaluates the one-step-ahead (OSA) negative log-likelihood for a single
conventional-tag likelihood family from a flat tracked observation
vector. Fills `nLL_arr` in the same grid order used by
[`pack_tag_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/pack_tag_osa.md).

## Usage

``` r
eval_tag_osa(
  nLL_arr,
  tracked,
  family,
  like_type,
  pred_recap,
  obs_recap,
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
  ln_theta = 0,
  zero_init = TRUE
)
```

## Arguments

- nLL_arr:

  Array of negative log-likelihood values to update.

- tracked:

  Flat tracked observation vector produced by
  [`pack_tag_osa()`](https://chengmatt.github.io/SPoRC/dev/reference/pack_tag_osa.md).

- family:

  Character, either `"count"` or `"comp"`.

- like_type:

  Integer likelihood type code (0-5).

- pred_recap:

  Predicted recapture array.

- obs_recap:

  Observed recapture array. Used to derive the recapture total \\N\\ for
  recapture-conditioned Dirichlet-multinomial (`like_type` 5).

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

- ln_theta:

  Log overdispersion parameter for the Dirichlet-multinomial composition
  family.

- zero_init:

  Logical; if `TRUE`, zero `nLL_arr` on entry.

## Value

Updated `nLL_arr` array with OSA negative log-likelihood contributions
filled in loop/grid order.

## Details

The Dirichlet-multinomial concentration matches the fitting likelihood:
the total \\N\\ is the number of tags released for release-conditioned
families (`like_type` 2, 4) and the total recaptures for
recapture-conditioned families (`like_type` 3, 5). Both totals are taken
from raw data (`tagged_fish`, `obs_recap`) rather than from the tracked
observation, since summing an OSA-tagged S4 object is not permitted.
