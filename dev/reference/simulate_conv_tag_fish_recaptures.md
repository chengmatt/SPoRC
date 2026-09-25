# Simulate conventional tag recaptures for fishery fleets

Draws observed recapture counts for one liberty, season and cohort cell
from the predicted recaptures, under the Poisson, negative binomial, or
the release- and recovery-conditioned multinomial and
Dirichlet-multinomial. Dims absent from `tag_recaptures_attr` are summed
over and the recaptures are written into index 1 of each.

## Usage

``` r
simulate_conv_tag_fish_recaptures(
  conv_fish_tag_like,
  tag_recaptures_attr,
  conv_tagged_fish,
  pred_conv_tag_fish_recap,
  obs_conv_tag_fish_recap,
  ln_conv_fish_tag_theta,
  ry,
  rseas,
  tc,
  sim,
  n_pop,
  n_regions,
  n_ages,
  n_sexes,
  n_fish_fleets
)
```

## Arguments

- conv_fish_tag_like:

  Integer likelihood: `0` Poisson, `1` negative binomial, `2` and `3`
  the release- and recovery-conditioned multinomial, `4` and `5` the two
  Dirichlet-multinomials.

- tag_recaptures_attr:

  Which dims are attended in the recapture likelihood, built from `"p"`,
  `"a"` and `"s"` joined by underscores. Region and fleet are always
  kept, and the rest are summed over into index 1.

- conv_tagged_fish:

  Released tagged fish
  `[n_conv_tag_cohorts × n_pop × n_ages × n_sexes × n_sims]`, the
  release sample size for the release-conditioned forms.

- pred_conv_tag_fish_recap:

  Predicted recaptures
  `[conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop × n_regions × n_ages × n_sexes × n_fish_fleets × n_sims]`.

- obs_conv_tag_fish_recap:

  Observed recaptures on the same dims, written in place at the
  `[ry, rseas, tc, ...]` slice.

- ln_conv_fish_tag_theta:

  Log overdispersion: the negative binomial size is
  `exp(ln_conv_fish_tag_theta)` and the Dirichlet-multinomial
  concentration `exp(ln_conv_fish_tag_theta) × N × p`. Ignored by the
  Poisson and the multinomial.

- ry, rseas, tc, sim:

  Years at liberty, recovery season, tag cohort and replicate indices.

- n_pop, n_regions, n_ages, n_sexes, n_fish_fleets:

  Dimension sizes.

## Value

`obs_conv_tag_fish_recap` with the draws filled in at
`[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim]`, the
summed dims fixed at index 1.

## Details

The release-conditioned forms express the predictions as proportions of
the tags released, appending a not-recaptured bin to complete the
probability vector before the draw and removing it afterwards. The
recovery-conditioned forms condition on the total predicted recaptures
and need no such bin.
