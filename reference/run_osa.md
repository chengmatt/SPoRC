# Internal function to compute OSA residuals for a single composition slice

Backend function used by \[get_osa()\] to calculate one-step-ahead (OSA)
residuals from observed and expected composition data.

## Usage

``` r
run_osa(
  obs,
  exp,
  N = NULL,
  DM_theta = NULL,
  LN_Sigma = NULL,
  fleet,
  index,
  years,
  index_label,
  comp_like
)
```

## Arguments

- obs:

  Matrix of observed compositions (rows = years, columns = bins)

- exp:

  Matrix of expected compositions (same shape as obs)

- N:

  Sample size for multinomial/Dirichlet-multinomial

- DM_theta:

  Dirichlet-multinomial overdispersion parameter(s)

- LN_Sigma:

  Logistic-normal covariance matrix

- fleet:

  Fleet identifier

- index:

  Vector of composition bins (last entry dropped for logistic-normal)

- years:

  Vector of years corresponding to rows of obs/exp

- index_label:

  Character describing ages or lengths

- comp_like:

  Integer specifying likelihood:

  0

  :   multinomial

  1

  :   Dirichlet-multinomial

  2-4

  :   logistic-normal variants

## Value

A list with one element:

- res:

  Data frame of OSA residuals with columns fleet, index_label, year,
  index, resid

## Details

For multinomial and Dirichlet-multinomial, observed proportions are
scaled to counts. For logistic-normal, bins with zeros are removed,
compositions renormalized, residuals computed in log-ratio space, and
mapped back to original bins.
