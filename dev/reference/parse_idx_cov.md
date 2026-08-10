# Validate fixed index covariance matrices for a multivariate normal likelihood

Checks each supplied covariance once at setup and returns it ready for
[`dmvnorm`](https://rdrr.io/pkg/RTMB/man/MVgauss.html). Each matrix must
be square with one row per observation the fleet actually fits, ordered
the way the observations appear when scanning the fleet's use flags in
array order (region varies fastest, then year, then season).

## Usage

``` r
parse_idx_cov(cov_list, like_type_vals, use_arr, n_fleets, what)
```

## Arguments

- cov_list:

  List with one element per fleet, each either a covariance matrix or
  `NULL`, or `NULL` for no matrices at all.

- like_type_vals:

  Integer vector of index likelihood codes; only fleets coded `2`
  require a matrix.

- use_arr:

  Array `[region, year, season, fleet]` of use flags.

- n_fleets:

  Integer. Number of fleets.

- what:

  Character. Name used in error messages.

## Value

List with one element per fleet, holding the validated covariance matrix
for fleets using the multivariate normal and `NULL` otherwise.

## Details

The checks are not decorative.
[`RTMB::dmvnorm`](https://rdrr.io/pkg/RTMB/man/MVgauss.html) reads only
the lower triangle without verifying symmetry, and returns `NaN`
silently when the covariance is not positive definite, so either mistake
would otherwise surface as an unexplained `NaN` objective rather than a
setup error.
