# Conditional Age-at-Length Likelihood

Computes the negative log-likelihood contribution for conditional
age-at-length (CAAL) data for a single year, season and fleet. A CAAL
observation is the age composition of the fish sampled from one length
bin, so each length bin is treated as an independent age composition and
evaluated through
[`Get_Comp_Likelihoods`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Comp_Likelihoods.md).

## Usage

``` r
Get_CAAL_Likelihoods(
  Exp,
  Obs,
  ISS,
  Wt_Mltnml,
  ln_theta,
  ln_theta_agg,
  Comp_Type,
  Likelihood_Type,
  n_regions,
  n_lens,
  n_model_bins,
  n_obs_bins,
  n_sexes,
  AgeingError,
  use,
  addtocomp,
  comp_bins = NULL,
  comp_const_obs = 1
)
```

## Arguments

- Exp:

  Expected joint numbers at length and age (a `Fish_caal`,
  `Fish_caal_discard` or `Srv_caal` slice), indexed by \\\[region \times
  len \times model\\age \times sex\]\\.

- Obs:

  Observed CAAL counts or proportions indexed by \\\[region \times len
  \times observed\\age \times sex\]\\.

- ISS:

  Input sample size indexed by \\\[region \times len \times sex\]\\.
  This is the number aged within the length bin, not the number
  measured.

- Wt_Mltnml:

  Multinomial weighting indexed by \\\[region \times len \times sex\]\\.

- ln_theta:

  Log overdispersion for the Dirichlet-multinomial, indexed by
  \\\[region \times sex\]\\. One value per region and sex is shared
  across length bins, since the bins come from one length-stratified
  sample rather than from independent surveys.

- ln_theta_agg:

  Log overdispersion used when `Comp_Type = 0`.

- Comp_Type:

  Integer specifying the composition parameterization, as in
  [`Get_Comp_Likelihoods`](https://chengmatt.github.io/SPoRC/dev/reference/Get_Comp_Likelihoods.md)
  (0 aggregated, 1 split by region and sex, 2 joint across sexes and
  split by region).

- Likelihood_Type:

  Integer specifying the likelihood family. Only `0` (multinomial) and
  `1` (Dirichlet-multinomial) are supported.

- n_regions:

  Number of regions modeled.

- n_lens:

  Number of length bins.

- n_model_bins:

  Number of age bins used internally in the model.

- n_obs_bins:

  Number of observed age bins.

- n_sexes:

  Number of sexes modeled.

- AgeingError:

  Ageing error matrix mapping model ages to observed ages.

- use:

  Integer matrix \\\[region \times len\]\\ indicating which cells have
  observations (`1` = use, `0` = ignore). A length bin with no aged fish
  in any region is skipped entirely.

- addtocomp:

  Small constant added to compositions to avoid numerical issues when
  zeros are present.

- comp_bins:

  Integer vector of age bins the composition is fitted over, or `NULL`
  (default) for all of them.

- comp_const_obs:

  Integer (0 or 1). Whether `addtocomp` is added to the observed
  proportions that weight the multinomial.

## Value

Array \\\[region \times len \times sex\]\\ of negative log-likelihood
contributions. Length bins with no observations stay at zero.
