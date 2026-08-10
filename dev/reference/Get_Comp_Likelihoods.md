# Composition Data Likelihood

Computes the negative log-likelihood contribution for composition data
(age or length) for a single year and fleet. The function supports
multiple composition parameterizations and likelihood families commonly
used in stock assessment models, including multinomial,
Dirichlet-multinomial, and logistic-normal likelihoods.

## Usage

``` r
Get_Comp_Likelihoods(
  Exp,
  Obs,
  ISS,
  Wt_Mltnml,
  ln_theta_agg,
  ln_theta,
  LN_corr_pars = 0,
  LN_corr_pars_agg = 0,
  Comp_Type,
  Likelihood_Type,
  n_regions,
  n_model_bins,
  n_obs_bins,
  n_sexes,
  age_or_len,
  AgeingError,
  use,
  addtocomp,
  comp_bins = NULL,
  comp_const_obs = 1
)
```

## Arguments

- Exp:

  Expected composition values (e.g., predicted catch-at-age or survey
  age compositions), structured as an array indexed by \\\[region \times
  model\\bins \times sex\]\\.

- Obs:

  Observed composition counts indexed by \\\[region \times
  observed\\bins \times sex\]\\.

- ISS:

  Input sample size for the composition data, indexed by \\\[region
  \times sex\]\\.

- Wt_Mltnml:

  Multinomial weighting applied to the effective sample size, indexed by
  \\\[region \times sex\]\\.

- ln_theta_agg:

  Log overdispersion parameter used when compositions are aggregated
  (`Comp_Type = 0`).

- ln_theta:

  Log overdispersion parameters used for Dirichlet-multinomial or
  logistic-normal likelihoods, indexed by \\\[region \times sex\]\\.

- LN_corr_pars:

  Logistic-normal correlation parameters used for correlated
  logistic-normal likelihoods, dimensioned by \\\[region \times sex
  \times parameters\]\\.

- LN_corr_pars_agg:

  Logistic-normal correlation parameters used when compositions are
  aggregated.

- Comp_Type:

  Integer specifying the composition parameterization:

  - `0`: Aggregated compositions across sexes and regions.

  - `1`: Compositions split by sex and region (no implicit sex or region
    ratio information).

  - `2`: Joint compositions across sexes but split by region (implicit
    sex ratio information).

- Likelihood_Type:

  Integer specifying the likelihood family:

  - `0`: Multinomial.

  - `1`: Dirichlet-multinomial.

  - `2`: Logistic-normal with independent bins.

  - `3`: Logistic-normal with AR(1) correlation across bins.

  - `4`: Logistic-normal with AR(1) correlation across bins and constant
    correlation across sexes.

- n_regions:

  Number of regions modeled.

- n_model_bins:

  Number of composition bins used internally in the model.

- n_obs_bins:

  Number of observed composition bins.

- n_sexes:

  Number of sexes modeled.

- age_or_len:

  Indicator for composition type:

  - `0`: Age compositions.

  - `1`: Length compositions.

- AgeingError:

  Ageing error matrix used to map model age bins to observed age bins.

- use:

  Integer vector indicating which regions have observations (`1` = use
  data, `0` = ignore).

- addtocomp:

  Small constant added to compositions to avoid numerical issues when
  zeros are present.

- comp_bins:

  Integer vector of bins the composition is fitted over, or `NULL`
  (default) for all of them. Both the observed and expected compositions
  are restricted to these bins and renormalized within them, so bins
  outside the range are left out of the likelihood rather than being
  forced to be explained. Indices refer to observed bins, that is after
  any ageing error has mapped model bins onto observed ones.

- comp_const_obs:

  Integer (0 or 1). Whether `addtocomp` is added to the observed
  proportions that weight the multinomial, as well as inside the
  logarithms. `1` (default) is the unbiased choice: the stationary point
  of the kernel is exactly `p = obs`. `0` weights by the raw observed
  proportions, which is what the ADMB templates do, and sharpens the
  expected composition away from the data by `n * addtocomp`. The
  difference is not cosmetic once effective sample sizes are large: on
  the EBS pollock bridge, switching from 0 to 1 moves estimated SSB by a
  median of 8.8 percent. Use `0` only to reproduce an ADMB model.

## Details

Expected and observed compositions are provided as arrays indexed by
region, composition bin (age or length), and sex. The function
optionally applies ageing error, aggregates compositions depending on
the parameterization type, and evaluates the likelihood for each region
and/or sex.
