# Evaluate an age-disaggregated observation likelihood

Every at-age data source in the model routes through here: retained
catch at age, discard catch at age, and the survey index at age, each in
an aggregated and a population-specific form. They differ only in which
array supplies the prediction and which parameter supplies the standard
deviation, so the likelihood itself is written once.

## Usage

``` r
get_at_age_nLL(
  obs_t,
  pred_t,
  sigma,
  corr_type = 0,
  rho = 0,
  ages = NULL,
  corr_mat = NULL
)
```

## Arguments

- obs_t:

  Numeric vector of observations for the ages present in one cell, on
  the scale the fleet's likelihood is written on.

- pred_t:

  Vector of predictions, matching `obs_t`.

- sigma:

  Vector of standard deviations, matching `obs_t`.

- corr_type:

  Integer. `0` is `"iid"`, `1` is `"1dar1"`, `2` is `"us"`.

- rho:

  Correlation for `corr_type = 1`, on the natural scale.

- ages:

  Integer vector of the ages present, used to space the AR(1). The
  position in the vector is assumed when this is `NULL`.

- corr_mat:

  Correlation matrix for `corr_type = 2`, already subset to the ages
  present.

## Value

Numeric vector the length of `obs_t`, the negative log likelihood
contribution of each age.

## Details

Ages within one cell may be independent, correlated as an AR(1) across
ages, or correlated through an unstructured matrix, the three structures
ICES age-structured assessments allow. A correlated cell contributes its
whole density to the first age present, leaving the remaining ages at
zero, so the returned vector still sums to the cell's contribution.

An AR(1) is a statement about age distance, not about position in the
observed vector. A fleet that observes ages 2, 3, 5 and 6 has a gap, and
lag one across that gap is not lag one, so the covariance is built from
the ages themselves whenever they are not consecutive. Consecutive ages
take the autoregressive recursion instead, which is the same density at
lower cost.
