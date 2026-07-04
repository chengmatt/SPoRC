# Compute shape parameters for a scaled beta distribution

Converts a mean and standard deviation expressed on the `[low, high]`
scale to the \\\alpha\\ and \\\beta\\ shape parameters of a beta
distribution defined on `[0, 1]` after location-scale transformation.
Used in SPoRC to construct beta priors for steepness (\\h\\) bounded to
\\\[0.2, 1\]\\.

## Usage

``` r
get_beta_scaled_pars(low, high, mu, sigma)
```

## Arguments

- low:

  Numeric. Lower bound of the parameter support.

- high:

  Numeric. Upper bound of the parameter support.

- mu:

  Numeric. Prior mean on the `[low, high]` scale.

- sigma:

  Numeric. Prior standard deviation on the `[low, high]` scale. Must
  satisfy \\\sigma^2 \< \bar{\mu}(1-\bar{\mu})\\ where \\\bar{\mu} =
  (\mu - \text{low}) / (\text{high} - \text{low})\\.

## Value

Numeric vector `c(alpha, beta, low, scale)` where `scale = high - low`.
