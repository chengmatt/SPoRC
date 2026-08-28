# Draw parameter deviations from the joint covariance

Draw parameter deviations from the joint covariance

## Usage

``` r
refpt_draw_devs(cov_obj, n_draw)
```

## Arguments

- cov_obj:

  Output of
  [`refpt_par_cov`](https://chengmatt.github.io/SPoRC/dev/reference/refpt_par_cov.md).

- n_draw:

  Integer. Number of draws.

## Value

Numeric matrix `[n_draw, n_par]` of mean zero deviations.
