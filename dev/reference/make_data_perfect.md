# Make the data effectively exact

Every observation gets a standard deviation of `obs_sd` and every
composition a sample size of `iss`. Process error is left alone, being
part of the truth. That covers the population-specific at-age sources,
the estimated part of an index sd, and a multivariate normal index,
whose covariance keeps its correlations at a marginal sd of `obs_sd`.
Tag recaptures are counts and keep their error.

## Usage

``` r
make_data_perfect(sim_list, obs_sd = 0.001, iss = 1e+06)
```

## Arguments

- sim_list:

  The simulation list, once the setup routines have filled it.

- obs_sd:

  Observation standard deviation to impose. Default `1e-3`.

- iss:

  Input sample size to impose on every composition. Default `1e6`.

## Value

`sim_list` with its observation error replaced.
