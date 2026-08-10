# Convert an index covariance matrix to common-factor parameters

A multivariate normal index likelihood is supplied as a full covariance
over the observation vector, which cannot be drawn from one year at a
time and does not extend past the years it covers, so it is unusable for
closed loop as given. This decomposes it into a marginal scale and a
single common factor,
`obs_t = pred_t + d_t (lambda_t u + sqrt(1 - lambda_t^2) e_t)`, with `u`
shared across years and `e_t` independent. Both are then drawn per year,
and a projection year past the end of the covariance simply reuses the
mean loading and scale with the same `u`.

## Usage

``` r
cov_to_factor(S)
```

## Arguments

- S:

  Covariance matrix over a fleet's observation vector.

## Value

List with `d` (marginal sd by observation) and `lambda` (factor loading
by observation, in (-1, 1)).
