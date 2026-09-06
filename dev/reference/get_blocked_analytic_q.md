# Solve a catchability analytically within each time block

A catchability that is a pure scaling nuisance can be concentrated out
of the likelihood instead of estimated. The solve is done separately
WITHIN each catchability time block, so a blocked catchability gets one
solved value per block. A single block, which is the default, reduces to
one value for the whole series.

## Usage

``` r
get_blocked_analytic_q(q_type, obs_vec, pred_vec, yr_obs, blk_yr, n_blk)
```

## Arguments

- q_type:

  1 for arithmetic scaling, 2 for geometric.

- obs_vec:

  Numeric vector of observed index values.

- pred_vec:

  Vector of predicted index values before scaling, same length.

- yr_obs:

  Integer vector of the year each observation belongs to.

- blk_yr:

  Integer vector `[n_yrs]` of each year's block.

- n_blk:

  Number of blocks this fleet and region use.

## Value

Vector `[n_yrs]` of the solved catchability by year.

## Details

Two forms, matching `q_type`:

- arithmetic (1):

  \\\hat q_b = \sum\_{i \in b} o_i / \sum\_{i \in b} p_i\\, the ratio of
  the block's mean observed to its mean predicted. The observation
  counts cancel, so only the two sums are needed.

- geometric (2):

  \\\log \hat q_b = \sum\_{i \in b}(\log o_i - \log p_i) / n_b\\, the
  block's mean log ratio, which is the lognormal maximum likelihood
  value of \\q\\.

A block with no observations has nothing to solve from and is given the
pooled value over every observation instead.
