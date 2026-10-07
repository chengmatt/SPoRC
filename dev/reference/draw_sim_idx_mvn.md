# Draw each multivariate normal index fleet's errors over its covariance

The cells a fleet's covariance covers are drawn together, one vector per
replicate, `t(chol(Cov)) %*% z`, so the draws have the covariance the
estimation model's `dmvnorm` reads exactly. Cells outside it, a closed
loop's projection years, keep the common-factor approximation of
[`cov_to_factor`](https://chengmatt.github.io/SPoRC/dev/reference/cov_to_factor.md).

## Usage

``` r
draw_sim_idx_mvn(sim_env)
```

## Arguments

- sim_env:

  Simulation environment holding `fish_idx_mvn` or `srv_idx_mvn` from
  [`build_idx_factor`](https://chengmatt.github.io/SPoRC/dev/reference/build_idx_factor.md).

## Value

`invisible(NULL)`; `sim_env` gains `fish_idx_mvn_eps` and
`srv_idx_mvn_eps`, one `[covariance row, replicate]` matrix per mvn
fleet, `NULL` for the others.
