# Gives negative log liklelihood values for composition data for a given year and a given fleet (fishery or survey)

Gives negative log liklelihood values for composition data for a given
year and a given fleet (fishery or survey)

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
  comp_agg_type,
  addtocomp
)
```

## Arguments

- Exp:

  Expected values (catch at age or survey index at age) indexed for a
  given year and fleet (structured as a matrix by age and sex)

- Obs:

  Observed values (catch at age or survey index at age) indexed for a
  given year and fleet (structured as a matrix by age and sex)

- ISS:

  Input sample size indexed for a given year and fleet (structured as a
  vector w/ sexes)

- Wt_Mltnml:

  Mutlinomial weight (if any) for a given fleet (structured as a vector
  w/ sexes)

- ln_theta_agg:

  Log overdispersion parameter if comp_type == 0, but we want to
  estsimate either a dirichlet or multinomial

- ln_theta:

  Log theta overdispersion for Dirichlet mutlinomial (scalar or vector
  depending on if 'Split' or 'Joint')

- LN_corr_pars:

  Logistic normal correlation parameters (dimensioned by n_regions,
  n_sexes, and 3 parameters)

- LN_corr_pars_agg:

  Logistic normal correlation parameters if comps are aggregated (just
  dimensioned by length of 1 value)

- Comp_Type:

  Composition Parameterization Type (== 0, aggregated comps by sex, ==
  1, split comps by sex and region (no implicit sex and region ratio
  information), == 2, joint comps across sexes but split by region
  (implicit sex ratio information, but not region information))

- Likelihood_Type:

  Composition Likelihood Type (== 0, Multinomial, == 1 Dirichlet
  Multinomial)

- n_regions:

  number of regions modeled

- n_model_bins:

  Number of bins used in the model

- n_obs_bins:

  Number of observed composition bins

- n_sexes:

  Number of sexes modeled

- age_or_len:

  Age or length comps (== 0, Age, == 1, Length)

- AgeingError:

  Ageing Error matrix

- use:

  Vector of 0s and 1s corresponding to regions (==0, don't have obs and
  dont' use, ==1, have obs and use)

- comp_agg_type:

  How to aggregate data (if aggregating)

- addtocomp:

  Small constant to add to composition data

## Value

Returns negative log likelihood for composition data (age and/or length)
