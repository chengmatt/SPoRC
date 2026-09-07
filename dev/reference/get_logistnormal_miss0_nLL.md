# Evaluate a logistic normal composition with the zeros dropped

The zeros are dropped rather than nudged with a constant, and the
expected proportions are renormalized over the bins that remain, so the
density is the one the positive part of the composition actually has.
The input sample size scales the variance, which is how a year sampled
harder comes to be fit more tightly than a year sampled lightly, and the
change of variables from the log ratio is taken off so the result is a
density on the composition rather than on its transform.

## Usage

``` r
get_logistnormal_miss0_nLL(
  obs,
  pred,
  ln_sigma,
  ISS,
  corr_type = 0,
  trans_rho = 0,
  lag_bins = NULL,
  corr_mat = NULL
)
```

## Arguments

- obs:

  Observed proportions in one cell, zeros included.

- pred:

  Predicted proportions, matching `obs`.

- ln_sigma:

  Log-scale standard deviation before the sample size scaling.

- ISS:

  Input sample size for this cell.

- corr_type:

  Integer. `0` independent bins, `1` first-order autoregressive across
  bins, `2` whatever correlation `corr_mat` states, which is how the
  separable bin by sex structure is passed in.

- trans_rho:

  Unconstrained correlation, mapped through the logistic function so the
  correlation is positive.

- lag_bins:

  Bin number of each element of `obs`, so the autoregression is spaced
  by bin rather than by position when the composition is fit over a
  restricted set of bins. Position is assumed when this is `NULL`. Read
  only by `corr_type = 1`.

- corr_mat:

  Correlation matrix over the whole of `obs`, before any bin is dropped.
  Required by `corr_type = 2` and ignored otherwise; it is cut down to
  the bins that were seen here rather than by the caller, so the
  correlation between two bins is the one the full structure gives them.

## Value

The negative log likelihood of the cell, a scalar.

## Details

A cell with fewer than two positive bins has no log ratio to take and
contributes nothing.
