# Simulate conventional tag recaptures for fishery fleets

Draws observed tag recapture counts for a single liberty–season–cohort
cell from predicted recapture arrays, supporting six likelihood
structures: Poisson, negative binomial, and release- or
recovery-conditioned multinomial and Dirichlet-multinomial. Dimensions
absent from `tag_recaptures_attr` are marginalised by summing over them,
and all recaptures are placed into index 1 of the corresponding
dimension in the output array.

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

  Integer. Likelihood for tag recaptures: `0` = Poisson, `1` = negative
  binomial, `2` = multinomial (release-conditioned), `3` = multinomial
  (recovery-conditioned), `4` = Dirichlet-multinomial
  (release-conditioned), `5` = Dirichlet-multinomial
  (recovery-conditioned).

- tag_recaptures_attr:

  Character string specifying which biological dimensions are attended
  in the recapture likelihood. Built from any combination of `"p"`
  (population), `"a"` (age), and `"s"` (sex), joined by underscores.
  Region and fleet are always retained. Unattended dimensions are
  marginalised and output into index 1.

- conv_tagged_fish:

  Array of released tagged fish
  `[n_conv_tag_cohorts × n_pop × n_ages × n_sexes × n_sims]`. Used as
  the release sample size for release-conditioned likelihoods.

- pred_conv_tag_fish_recap:

  Array of predicted recaptures
  `[conv_tag_max_liberty × n_seas × n_conv_tag_cohorts × n_pop × n_regions × n_ages × n_sexes × n_fish_fleets × n_sims]`.

- obs_conv_tag_fish_recap:

  Array of observed recaptures with the same dimensions as
  `pred_conv_tag_fish_recap`. Simulated values are written in-place at
  the `[ry, rseas, tc, ...]` slice.

- ln_conv_fish_tag_theta:

  Numeric. Log overdispersion: negative binomial size =
  `exp(ln_conv_fish_tag_theta)`; Dirichlet-multinomial concentration =
  `exp(ln_conv_fish_tag_theta) × N × p`.

- ry:

  Integer. Years-at-liberty index (first dimension of recapture arrays).

- rseas:

  Integer. Recovery season index.

- tc:

  Integer. Tag cohort index.

- sim:

  Integer. Simulation replicate index.

- n_pop:

  Integer. Number of populations.

- n_regions:

  Integer. Number of regions.

- n_ages:

  Integer. Number of age classes.

- n_sexes:

  Integer. Number of sexes.

- n_fish_fleets:

  Integer. Number of fishery fleets.

## Value

The `obs_conv_tag_fish_recap` array with simulated recaptures filled in
at `[ry, rseas, tc, pop_idx, reg_idx, age_idx, sex_idx, flt_idx, sim]`.
Marginalised dimensions are fixed at index 1.

## Details

For release-conditioned likelihoods (`2`, `4`), predicted recaptures are
expressed as proportions of total tags released. A "not-recaptured" bin
is appended to complete the probability vector before drawing and
removed before assignment. For recovery-conditioned likelihoods (`3`,
`5`), the draw is conditioned on the total predicted recapture count
with no not-recaptured bin needed. The overdispersion parameter
`ln_conv_fish_tag_theta` governs the negative-binomial size parameter
and the Dirichlet-multinomial concentration scaling, and is ignored for
Poisson and multinomial likelihoods.
