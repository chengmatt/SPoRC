# Evaluate an age-disaggregated observation likelihood

Every at-age observation stream in the model routes through here:
retained catch at age, discard catch at age, and the fishery and survey
indices at age, each in an aggregated and a population-specific form.
They differ only in which array supplies the prediction and which
parameter supplies the standard deviation, so the likelihood itself is
written once.

## Usage

``` r
get_at_age_nLL(obs_log, pred_log, sigma, corr_type = 0, rho = 0)
```

## Arguments

- obs_log:

  Numeric vector of log observations for the ages present in one cell,
  already offset by any constant.

- pred_log:

  Vector of log predictions, matching `obs_log`.

- sigma:

  Vector of standard deviations, matching `obs_log`.

- corr_type:

  Integer. `0` is `"iid"`, `1` is `"1dar1"`.

- rho:

  Correlation for `corr_type = 1`, on the natural scale.

## Value

Numeric vector the length of `obs_log`, the negative log likelihood
contribution of each age.

## Details

Observations are lognormal on the natural scale. Ages within one cell
may be independent, or correlated as an AR(1) across ages, as ICES
age-structured assessments allow. A correlated cell contributes its
whole density to the first age present, leaving the remaining ages at
zero, so the returned array still sums to the cell's contribution.
