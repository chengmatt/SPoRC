# Draw an observed index given its likelihood family

The estimation model supports lognormal, arithmetic-scale normal, and
multivariate normal index likelihoods. The simulator drew lognormal
error unconditionally, so a self-test of a normal or MVN index generated
data under the wrong error structure and then reported the mismatch as
estimation bias.

## Usage

``` r
draw_index_obs(true, se, like_type = 0, d = NULL, lambda = NULL, u = NULL)
```

## Arguments

- true:

  True (error-free) index value(s).

- se:

  Observation standard deviation, on the log scale for lognormal and the
  arithmetic scale for normal. Unused for MVN, which takes its scale
  from the covariance instead; note the two are not interchangeable, the
  pollock trawl covariance has a diagonal about twice the reported SEs.

- like_type:

  0 lognormal, 1 normal, 2 multivariate normal.

- d, lambda:

  Common-factor scale and loading for this observation, from
  [`cov_to_factor`](https://chengmatt.github.io/SPoRC/dev/reference/cov_to_factor.md).
  MVN only.

- u:

  Shared factor draw for this fleet and replicate, kept constant across
  years. MVN only.
