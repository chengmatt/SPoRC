# Composition Data Likelihood (OSA variant)

Evaluates the multinomial (0), Dirichlet-multinomial (1) and
logistic-normal (2 iid, 3 AR1, 4 2D-AR1) composition likelihoods on one
flat tracked OBS vector for
[`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html),
respecting the reduced logistic-normal block lengths the packer used.

## Usage

``` r
Get_Comp_Likelihoods_OSA(
  Exp,
  Obs,
  ISS,
  ln_theta,
  ln_theta_agg,
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
  comp_bins = NULL
)
```

## Arguments

- Exp:

  Expected proportions `[n_regions × n_model_bins × n_sexes]`.

- Obs:

  Flat tracked observation vector, already transformed for the logistic
  normal.

- ISS:

  Input sample size `[n_regions × n_sexes]`.

- ln_theta:

  Log overdispersion `[n_regions × n_sexes]`.

- ln_theta_agg:

  Log overdispersion scalar for the aggregated comps.

- LN_corr_pars:

  Logistic-normal correlation parameters `[n_regions × n_sexes × 3]`.

- LN_corr_pars_agg:

  Its aggregated counterpart.

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
  For length compositions, either `NA` (the model and observed bins
  coincide) or a matrix `[n_model_bins x n_obs_bins]` mapping the
  model's length bins onto the observed ones, the same way, when the
  compositions are recorded on coarser bins than the model has.

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
  any ageing error or length bin map has mapped model bins onto observed
  ones. The restriction applies to every composition type: for the
  sex-joint comps (`Comp_Type = 2`) the named bins are dropped from each
  sex's block of the joint stack, so the sex ratio the joint comps have
  is the ratio within the fitted bins. Logistic-normal covariances are
  built over all observed bins and then cut down to the fitted ones, so
  a gap in `comp_bins` still counts towards the AR1 lag between the bins
  on either side of it.

## Details

The tracked `Obs` vector is never reshaped. Everything on the
expectation side (`Exp`, `ISS`, `ln_theta`, `LN_corr_pars`, the ageing
error) is reshaped and filtered by `use` exactly as in fitting.

Because [`RTMB::OBS()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)
cannot be altered after tracking, the additive log ratio transform of
the observation happens in the packer, so `Obs[idx]` arrives already
transformed with its last bin dropped. Here only the expectation is
transformed, `Sigma` is built with its last row and column dropped, and
the multivariate normal density is evaluated. A block is
`n_fit_bins - 1` long under comp type 0, the same per region and sex
under type 1, and `n_fit_bins * n_sexes - 1` per region under type 2,
which takes one joint reference. `n_fit_bins` is the number of bins
`comp_bins` names, equal to `n_obs_bins` when the fleet fits every bin,
and the packer applies the same restriction before transforming, so the
reference is the last fitted bin rather than the last observed one.
