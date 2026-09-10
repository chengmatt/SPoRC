# Index likelihoods for the fleets not fitted on the log scale

The regional index likelihood for every fleet whose observations stay
untransformed: arithmetic-scale normal, and multivariate normal.
Lognormal fleets are left to
[`eval_index_osa_nLL`](https://chengmatt.github.io/SPoRC/dev/reference/eval_index_osa_nLL.md).

## Usage

``` r
get_index_regional_nLL(
  nLL_arr,
  Use,
  Obs,
  Pred,
  SD,
  LikeType,
  Cov,
  seas_Type,
  const,
  n_fleets
)
```

## Arguments

- nLL_arr:

  Container `[region, year, season, fleet]`.

- Use:

  Indicator array of the cells that are fit.

- Obs:

  Observed index, untransformed.

- Pred:

  Predicted index `[pop, region, year, season, fleet]`.

- SD:

  Index standard deviation.

- LikeType:

  Integer vector, one likelihood code per fleet.

- Cov:

  List of covariance matrices, one per fleet, read only by a
  multivariate normal.

- seas_Type:

  Integer vector, whether each fleet is fit as a season total.

- const:

  Constant added inside the log.

- n_fleets:

  Number of fleets.

## Value

`nLL_arr` with the fitted cells filled.

## Details

Call this before that one. Registering an observation for OSA residuals
binds the flattened log-scale vector back to the observation's own name,
since [`RTMB::OBS`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)
files residuals under the name it is called on, so the untransformed
array is gone once the lognormal fleets have been done.

The fitted cells are flattened into vectors because a multivariate
normal reads the whole series at once. Its total lands in the first cell
and the rest stay zero, so the fleet's likelihood still sums correctly.
