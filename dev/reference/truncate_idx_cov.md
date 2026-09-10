# Cut Fixed Index Covariances and Index Weights to a Retrospective Peel

A multivariate normal index likelihood reads one covariance row per
fitted observation, ordered the way the fleet's use flags are scanned
(region fastest, then year, then season). A peel removes observations,
so the covariance is cut to the rows and columns those cells pick out,
the marginal over what the peel kept. Year-indexed index weights are cut
the same way.

## Usage

``` r
truncate_idx_cov(retro_data, data)
```

## Arguments

- retro_data:

  Truncated data list from
  [`truncate_yr`](https://chengmatt.github.io/SPoRC/dev/reference/truncate_yr.md),
  with any data lags already applied to its use flags.

- data:

  Full data list the peel was taken from.

## Value

`retro_data` with its index covariances and index weights matching the
observations the peel kept.

## Details

Call this once every use flag has settled, so a lagged data source is
peeled the same way a terminal year is.
