# Conventional Tagging Likelihood Block

Computes the negative log-likelihood contribution from conventional fish
tagging data under Poisson, negative binomial, multinomial, or
Dirichlet-multinomial likelihoods.

## Usage

``` r
get_conv_tag_likelihoods(
  n_conv_tag_cohorts,
  conv_tag_release_indicator,
  conv_tag_max_liberty,
  n_yrs,
  n_seas,
  conv_tag_mixing_period,
  n_fish_fleets,
  use_conv_fish_tagging,
  n_conv_tag_pop_pool,
  n_regions,
  n_conv_tag_age_pool,
  n_conv_tag_sex_pool,
  conv_tag_pop_pool,
  conv_tag_age_pool,
  conv_tag_sex_pool,
  conv_fish_tag_like,
  conv_fish_tag_nLL,
  obs_conv_tag_fish_recap,
  pred_conv_tag_fish_recap,
  addtotag,
  ln_conv_fish_tag_theta,
  conv_tagged_fish,
  zero_init = TRUE
)
```

## Arguments

- n_conv_tag_cohorts:

  Number of conventional tag cohorts.

- conv_tag_release_indicator:

  Matrix giving release region, year, season for each cohort.

- conv_tag_max_liberty:

  Maximum years at liberty to evaluate.

- n_yrs:

  Total number of modeled years.

- n_seas:

  Number of seasons per year.

- conv_tag_mixing_period:

  Minimum seasons at liberty before tags are modeled.

- n_fish_fleets:

  Number of fishing fleets.

- use_conv_fish_tagging:

  Vector indicating which fleets use conventional tagging.

- n_conv_tag_pop_pool:

  Number of population pooling groups.

- n_regions:

  Number of spatial regions.

- n_conv_tag_age_pool:

  Number of age pooling groups.

- n_conv_tag_sex_pool:

  Number of sex pooling groups.

- conv_tag_pop_pool:

  List of population index pools.

- conv_tag_age_pool:

  List of age index pools.

- conv_tag_sex_pool:

  List of sex index pools.

- conv_fish_tag_like:

  Likelihood type indicator (0-5).

- conv_fish_tag_nLL:

  Array of negative log-likelihood values to update.

- obs_conv_tag_fish_recap:

  Observed recapture array.

- pred_conv_tag_fish_recap:

  Predicted recapture array.

- addtotag:

  Small constant added to avoid zeros.

- ln_conv_fish_tag_theta:

  Log overdispersion parameter for NB/Dirichlet-multinomial.

- conv_tagged_fish:

  Array of numbers of tagged fish released.

- zero_init:

  Whether to zero the nLL array on entry.

## Value

Updated \`conv_fish_tag_nLL\` array.

## Details

This function loops over tag cohorts, recapture years, seasons, fleets,
population pools, regions, age pools, and sex pools, accumulating
likelihood contributions into \`conv_fish_tag_nLL\`.
