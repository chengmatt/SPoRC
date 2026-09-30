# Mean of a covariate observation, the grid cell through its link

Identity, exp, inverse logit or inverse cloglog, by link code.

## Usage

``` r
get_dsem_cov_mean(x, link)
```

## Arguments

- x:

  Grid cells, on the link scale.

- link:

  Link code: 0 identity, 1 log, 2 logit, 3 cloglog.

## Value

The mean of the observation, on the scale it is observed.
