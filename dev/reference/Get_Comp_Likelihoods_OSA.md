# Composition Data Likelihood (OSA variant)

Computes multinomial (0), Dirichlet-multinomial (1), and logistic-normal
(2 iid, 3 AR1, 4 2D‑AR1) composition likelihoods for one‑step‑ahead
(OSA) residuals using
[`RTMB::oneStepPredict`](https://rdrr.io/pkg/RTMB/man/OSA-residuals.html).
The function evaluates the likelihood for a single flat tracked OBS
vector, respecting the reduced logistic‑normal block lengths used during
packing.

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
  addtocomp
)
```

## Arguments

- Exp:

  Expected proportions \[n_regions × n_model_bins × n_sexes\].

- Obs:

  Flat tracked observation vector (already ALR‑transformed for LN).

- ISS:

  Input sample size \[n_regions × n_sexes\].

- ln_theta:

  Log overdispersion \[n_regions × n_sexes\].

- ln_theta_agg:

  Log overdispersion scalar for aggregated comps.

- LN_corr_pars:

  LN correlation parameters \[n_regions × n_sexes × 3\].

- LN_corr_pars_agg:

  LN aggregated correlation scalar(s).

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

## Details

The tracked `Obs` vector is \*\*never reshaped\*\*. All expectation‑side
quantities (`Exp`, `ISS`, `ln_theta`, `LN_corr_pars`, ageing error) are
reshaped and filtered by `use`, exactly as in fitting.

\*\*Logistic‑normal note:\*\* Because
[`RTMB::OBS()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html) cannot
be altered after tracking, the additive‑log‑ratio (ALR) transform of the
\*observation\* is performed in the packer. Thus, `Obs[idx]` is
\*\*already ALR‑transformed\*\* (last bin dropped). Here we only
ALR‑transform the expectation, construct the covariance matrix `Sigma`
(dropping its last row/column), and evaluate the multivariate normal
density.

Reduced LN block lengths: \* Comp_Type 0: `n_obs_bins - 1` \* Comp_Type
1: `n_obs_bins - 1` per region/sex \* Comp_Type 2:
`n_obs_bins * n_sexes - 1` per region (one joint reference)
